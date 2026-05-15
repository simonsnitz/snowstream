"""Precompute Snowprint predictions for every representative of a TF family.

Usage:
  # Run all stages end-to-end (resumable: rerun is idempotent):
  python -m scripts.precompute_family tetr

  # Run a single stage (e.g. just the prediction loop):
  python -m scripts.precompute_family tetr --stage predict

  # Smoke-test on a small slice:
  python -m scripts.precompute_family tetr --max-clusters 5

  # Re-run a stage that previously completed:
  python -m scripts.precompute_family tetr --stage cluster --force

  # Skip PaperBLAST (e.g. while it's blocked behind Cloudflare):
  python -m scripts.precompute_family tetr --no-paperblast

Stages run in this order; each is idempotent (skips if its output exists,
unless --force is passed):

  1. members      — fetch UniProt IDs annotated with the family's InterPro entry
  2. sequences    — fetch FASTA sequences for those IDs
  3. cluster      — MMseqs2 cluster at 50 % identity (configurable)
  4. representatives — pick a rep per cluster (groovDB > PaperBLAST > centroid)
  5. dmnd         — build representatives.fasta + representatives.dmnd
  6. predict      — run the full Snowprint pipeline for each rep, append JSONL

The predict stage reads `SNOWPRINT_DIAMOND_DB` to pick which BLAST DB to query
(default: the small bHTH DB shipped in `databases/`). Set this to your local
NCBI nr Diamond DB before running the full TetR precompute.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.precompute import (  # noqa: E402
    cluster as cluster_mod,
    cluster_search as cluster_search_mod,
    dmnd as dmnd_mod,
    groovdb as groovdb_mod,
    interpro as interpro_mod,
    members_dmnd as members_dmnd_mod,
    members_predict as members_predict_mod,
    paperblast as paperblast_mod,
    predict as predict_mod,
    representatives as representatives_mod,
    sequences as sequences_mod,
)

ALL_STAGES = (
    "members",
    "sequences",
    "cluster",
    "representatives",
    "dmnd",
    "cluster_search",
    "predict",
    "members_predict",
    "members_dmnd",
)


def load_manifest(family_key: str) -> tuple[dict, Path]:
    fam_dir = PROJECT_ROOT / "families" / family_key
    manifest_path = fam_dir / "family.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"no manifest at {manifest_path}")
    return json.loads(manifest_path.read_text()), fam_dir


def stage_members(manifest: dict, fam_dir: Path, force: bool) -> None:
    interpro_mod.fetch_members(
        interpro_id=manifest["interpro_id"],
        output_path=fam_dir / "members.txt",
        force=force,
    )


def stage_sequences(_manifest: dict, fam_dir: Path, force: bool) -> None:
    sequences_mod.fetch_sequences(
        members_path=fam_dir / "members.txt",
        output_path=fam_dir / "members.fasta",
        force=force,
    )


def stage_cluster(manifest: dict, fam_dir: Path, force: bool) -> None:
    cluster_mod.cluster(
        fasta_path=fam_dir / "members.fasta",
        output_path=fam_dir / "clusters.tsv",
        identity=manifest.get("cluster_identity_threshold", 0.5),
        force=force,
    )


def stage_representatives(
    manifest: dict,
    fam_dir: Path,
    force: bool,
    use_paperblast: bool,
    max_clusters: int | None,
) -> None:
    representatives_path = fam_dir / "representatives.json"
    if representatives_path.exists() and not force:
        logging.info("representatives.json already exists; skipping")
        return

    clusters = cluster_mod.parse_clusters(fam_dir / "clusters.tsv")
    if max_clusters is not None:
        # Deterministic: take the first N centroids alphabetically.
        keep = sorted(clusters.keys())[:max_clusters]
        clusters = {k: clusters[k] for k in keep}
        logging.info("limited to %d clusters for smoke test", len(clusters))

    sequences = _read_sequences(fam_dir / "members.fasta")
    groov_dump = groovdb_mod.load(PROJECT_ROOT / "families" / "_resources" / "groovdb.json")

    pblast_cache = PROJECT_ROOT / "families" / "_resources" / "paperblast_cache" if use_paperblast else None

    def progress(i: int, n: int, centroid: str) -> None:
        if i % 100 == 0:
            logging.info("rep selection %d/%d (cluster %s)", i, n, centroid)

    records = representatives_mod.build_representatives(
        clusters=clusters,
        sequences=sequences,
        groovdb_dump=groov_dump,
        groovdb_entry_prefix=manifest.get("groovdb_entry_prefix", ""),
        paperblast_cache_dir=pblast_cache,
        use_paperblast=use_paperblast,
        on_progress=progress,
    )
    representatives_mod.write_representatives(records, representatives_path)
    logging.info("wrote %d representatives to %s", len(records), representatives_path)


def stage_dmnd(_manifest: dict, fam_dir: Path, force: bool) -> None:
    representatives = representatives_mod.read_representatives(fam_dir / "representatives.json")
    dmnd_mod.write_representatives_fasta(
        members_fasta=fam_dir / "members.fasta",
        representatives_json=representatives,
        output_path=fam_dir / "representatives.fasta",
    )
    dmnd_mod.build_dmnd(
        fasta_path=fam_dir / "representatives.fasta",
        dmnd_path=fam_dir / "representatives.dmnd",
        force=force,
    )


def stage_cluster_search(_manifest: dict, fam_dir: Path, force: bool, max_clusters: int | None) -> None:
    representatives = representatives_mod.read_representatives(fam_dir / "representatives.json")
    if max_clusters is not None:
        representatives = representatives[:max_clusters]
    cluster_search_mod.search_cluster_homologs(
        representatives_fasta=fam_dir / "representatives.fasta",
        members_fasta=fam_dir / "members.fasta",
        representatives=representatives,
        output_path=fam_dir / "cluster_homologs.json",
        force=force,
    )


def stage_members_predict(
    _manifest: dict,
    fam_dir: Path,
    max_members: int | None,
    workers: int,
    batch_size: int,
) -> None:
    out_path = fam_dir / "members_predictions.jsonl"
    cluster_predictions = fam_dir / "predictions.jsonl"

    # 1. Lift everything we can from the existing per-cluster predictions.jsonl
    if cluster_predictions.exists():
        n_lifted = members_predict_mod.lift_from_clusters(cluster_predictions, out_path)
        logging.info("lifted %d records from %s", n_lifted, cluster_predictions.name)
    else:
        logging.warning("no %s found; nothing to lift", cluster_predictions)

    # 2. Process the gap (members in members.fasta but not yet in the JSONL)
    member_ids = list(members_predict_mod.iter_member_ids(fam_dir / "members.fasta"))
    logging.info("members in fasta: %d", len(member_ids))

    def progress(i: int, total: int, _elapsed: int) -> None:
        if (i + 1) % 50 == 0 or i + 1 == total:
            logging.info("[members_predict batch %d/%d]", i + 1, total)

    members_predict_mod.predict_members(
        member_ids=member_ids,
        output_path=out_path,
        batch_size=batch_size,
        workers=workers,
        max_entries=max_members,
        on_progress=progress,
    )


def stage_members_dmnd(_manifest: dict, fam_dir: Path, force: bool) -> None:
    members_dmnd_mod.build_members_dmnd(
        members_fasta=fam_dir / "members.fasta",
        predictions_path=fam_dir / "members_predictions.jsonl",
        output_fasta=fam_dir / "members_with_promoters.fasta",
        output_dmnd=fam_dir / "members_with_promoters.dmnd",
        force=force,
    )


def stage_predict(_manifest: dict, fam_dir: Path, max_clusters: int | None, workers: int) -> None:
    representatives = representatives_mod.read_representatives(fam_dir / "representatives.json")
    sequences = _read_sequences(fam_dir / "members.fasta")
    cluster_homologs = json.loads((fam_dir / "cluster_homologs.json").read_text())

    def progress(i: int, n: int, uid: str, status: str) -> None:
        logging.info("[predict %d/%d] %s — %s", i + 1, n, uid, status)

    predict_mod.predict_from_clusters(
        representatives=representatives,
        sequences=sequences,
        cluster_homologs=cluster_homologs,
        jsonl_path=fam_dir / "predictions.jsonl",
        max_entries=max_clusters,
        workers=workers,
        on_progress=progress,
    )


def _read_sequences(fasta_path: Path) -> dict[str, str]:
    """Parse a FASTA into {uniprot_id: sequence}."""
    if not fasta_path.exists():
        return {}
    out: dict[str, str] = {}
    current_id: str | None = None
    current_chunks: list[str] = []
    for line in fasta_path.read_text().splitlines():
        if line.startswith(">"):
            if current_id is not None:
                out[current_id] = "".join(current_chunks)
            body = line[1:].strip()
            if "|" in body:
                parts = body.split("|")
                current_id = parts[1] if len(parts) >= 2 and parts[0] in ("sp", "tr") else body.split()[0]
            else:
                current_id = body.split()[0]
            current_chunks = []
        else:
            current_chunks.append(line.strip())
    if current_id is not None:
        out[current_id] = "".join(current_chunks)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("family", help="Family key (matches families/<key>/)")
    parser.add_argument(
        "--stage",
        choices=ALL_STAGES,
        help="Run only this stage (default: all stages in order)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run a stage even if its output already exists",
    )
    parser.add_argument(
        "--max-clusters",
        type=int,
        help="Limit rep-selection and predict stages to N clusters (smoke test)",
    )
    parser.add_argument(
        "--no-paperblast",
        dest="use_paperblast",
        action="store_false",
        help="Skip PaperBLAST queries during representative selection",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel workers for predict / members_predict (default: 1, sequential)",
    )
    parser.add_argument(
        "--max-members",
        type=int,
        help="Limit members_predict to N gap members (smoke test)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=50,
        help="Members per batch in members_predict (one NCBI IPG call per batch)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="DEBUG-level logging")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s",
    )
    logging.info("SNOWPRINT_DIAMOND_DB=%s", os.environ.get("SNOWPRINT_DIAMOND_DB", "<default bHTH>"))

    manifest, fam_dir = load_manifest(args.family)
    logging.info("family: %s (%s) → %s", manifest["name"], manifest["interpro_id"], fam_dir)

    stages_to_run = [args.stage] if args.stage else list(ALL_STAGES)
    for name in stages_to_run:
        logging.info("=== stage: %s ===", name)
        if name == "members":
            stage_members(manifest, fam_dir, args.force)
        elif name == "sequences":
            stage_sequences(manifest, fam_dir, args.force)
        elif name == "cluster":
            stage_cluster(manifest, fam_dir, args.force)
        elif name == "representatives":
            stage_representatives(manifest, fam_dir, args.force, args.use_paperblast, args.max_clusters)
        elif name == "dmnd":
            stage_dmnd(manifest, fam_dir, args.force)
        elif name == "cluster_search":
            stage_cluster_search(manifest, fam_dir, args.force, args.max_clusters)
        elif name == "predict":
            stage_predict(manifest, fam_dir, args.max_clusters, args.workers)
        elif name == "members_predict":
            stage_members_predict(manifest, fam_dir, args.max_members, args.workers, args.batch_size)
        elif name == "members_dmnd":
            stage_members_dmnd(manifest, fam_dir, args.force)


if __name__ == "__main__":
    main()
