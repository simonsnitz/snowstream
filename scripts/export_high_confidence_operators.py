"""Export a JSON file of high-confidence centroid operator predictions.

Filters the precomputed `families/<key>/members_predictions.jsonl` to
centroids meeting BOTH of:
  1. number of homologs (used in the V0 consensus calculation) ≥ 5
  2. operator consensus score > 75 for V0 OR V2.6 OR both

For each passing centroid, emits:
  * uniprot_id            — primary identifier
  * ncbi_accession        — RefSeq accession (looked up from the inverted
                            refseq_to_uniprot.json map; None if absent)
  * protein_sequence      — from members.fasta
  * genome                — NCBI nucleotide accession the operon was
                            extracted from
  * n_homologs            — count used in the consensus
  * operator_v0 / operator_v26 — each with:
      * motif                — consensus motif string (mixed-case: lowercase
                                 flanks + uppercase palindromic core)
      * consensus_score
      * rerank_score         — V2.6 only
      * native_for_centroid  — operator extracted from the centroid's own
                                 promoter
      * frequency_matrix     — PPM as [[pA, pC, pG, pT], ...] per position,
                                 driven directly by the aligned operators —
                                 ready to feed logojs-react's DNALogo
      * aligned_operators    — list of {uniprot_id, operator, align_score}
                                 for every homolog whose promoter aligned
                                 above the score cutoff (this is what built
                                 the consensus)

Since the JSONL only stores summary fields (motif + score + native), the
per-homolog aligned operators and the frequency matrix are regenerated
here by calling operator_fetch inline for each qualifying centroid. This
is fast (~1 minute with threading for the ~5k passing entries) because
the promoters are already in memory from the JSONL scan.

The output JSON is a list, sorted by uniprot_id for stable diffs across
re-runs. Default output path is
`algorithms/benchmark/datasets/tetr_high_confidence_operators.json`
so it sits alongside the other benchmark datasets.
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterator, Optional

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Re-use the exact per-version specs the precompute stage uses so the
# regenerated frequency_matrix + aligned_operators match what would be in
# the JSONL if we ever grew the block schema to include them.
from scripts.precompute.compute_operators import (  # noqa: E402
    V0_SPEC,
    V26_SPEC,
    _build_homolog_input,
)


def load_fasta(path: Path) -> dict[str, str]:
    """Parse members.fasta into {uniprot_id: sequence}. The FASTA headers
    use either `sp|UID|NAME` / `tr|UID|NAME` (UniProtKB) or bare UID, so
    handle both shapes."""
    out: dict[str, str] = {}
    current_id: str | None = None
    current_seq: list[str] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    out[current_id] = "".join(current_seq)
                body = line[1:].strip()
                if "|" in body:
                    parts = body.split("|")
                    if len(parts) >= 2 and parts[0] in ("sp", "tr"):
                        current_id = parts[1]
                    else:
                        current_id = body.split()[0]
                else:
                    current_id = body.split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
    if current_id is not None:
        out[current_id] = "".join(current_seq)
    return out


def invert_refseq_map(refseq_to_uniprot_path: Path) -> dict[str, str]:
    """Build a uniprot_id → first_refseq map from the inverse file. The
    on-disk map is refseq → uniprot; inverting collapses multi-refseq
    entries to whichever appears first, which is fine for display."""
    if not refseq_to_uniprot_path.exists():
        return {}
    forward = json.loads(refseq_to_uniprot_path.read_text())
    out: dict[str, str] = {}
    for refseq, uid in forward.items():
        out.setdefault(uid, refseq)
    return out


def load_latest_records(jsonl_path: Path) -> Iterator[dict]:
    """Yield the latest record per uniprot_id."""
    by_uid: dict[str, dict] = {}
    with jsonl_path.open() as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            uid = r.get("uniprot_id")
            if uid:
                by_uid[uid] = r
    yield from by_uid.values()


def passes_filters(rec: dict, min_homologs: int, score_cutoff: float) -> bool:
    v0 = rec.get("operator_v0")
    v26 = rec.get("operator_v26")
    if not (isinstance(v0, dict) and isinstance(v26, dict)):
        return False
    if (v0.get("n_homologs_used") or 0) < min_homologs:
        return False
    s0 = v0.get("consensus_score") or 0
    s26 = v26.get("consensus_score") or 0
    return s0 > score_cutoff or s26 > score_cutoff


def _render_motif(motif) -> Optional[str]:
    if not motif:
        return None
    if isinstance(motif, str):
        return motif if motif != "None" else None
    if isinstance(motif, list):
        try:
            return "".join(x["base"] for x in motif if isinstance(x, dict))
        except Exception:
            return None
    return None


def _version_block(spec, homolog_input: list[dict]) -> dict:
    """Run operator_fetch for one version and pack the parts the frontend
    needs into a JSON-serialisable dict."""
    try:
        result = spec.fetch(homolog_input, spec.fetch_params)
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}

    rerank_score = result.get("rerank_score")
    if rerank_score in (None, float("-inf")):
        rerank_score = None

    # aligned_seqs is what fed the consensus. Trim to the compact
    # {uniprot_id, operator, align_score} shape for the frontend.
    aligned = []
    for h in (result.get("aligned_seqs") or []):
        if not isinstance(h, dict):
            continue
        score = h.get("Align score")
        aligned.append({
            "uniprot_id": h.get("Uniprot Id"),
            "operator": h.get("Predicted operator"),
            "align_score": (round(float(score), 2)
                             if isinstance(score, (int, float)) else score),
        })

    # PPM: list of [A, C, G, T] per position, values in [0, 1]. Round to 3
    # decimals — plenty for a sequence logo and cuts file size roughly in
    # half vs full-precision floats.
    freq = result.get("frequency_matrix") or []
    freq_rounded = [
        [round(float(x), 3) for x in row] if isinstance(row, list) else row
        for row in freq
    ]

    cons = result.get("consensus_score")
    if isinstance(cons, (int, float)):
        cons = round(float(cons), 3)
    if rerank_score is not None:
        rerank_score = round(float(rerank_score), 4)

    return {
        "motif": _render_motif(result.get("motif")),
        "consensus_score": cons,
        "rerank_score": rerank_score,
        "native_for_centroid": result.get("native_operator"),
        "frequency_matrix": freq_rounded,
        "aligned_operators": aligned,
    }


def build_entry(rec: dict,
                 sequences: dict[str, str],
                 uid_to_refseq: dict[str, str],
                 cluster_members: list[dict],
                 records: dict[str, dict]) -> dict:
    """Assemble the per-centroid output entry — runs V0 and V2.6
    operator_fetch inline so the frequency matrix + per-homolog operators
    can be captured (the compact JSONL block only stores summary fields)."""
    uid = rec["uniprot_id"]
    homolog_input = _build_homolog_input(uid, cluster_members, records)
    return {
        "uniprot_id": uid,
        "ncbi_accession": uid_to_refseq.get(uid),
        "protein_sequence": sequences.get(uid),
        "genome": rec.get("genome"),
        "n_homologs": rec.get("operator_v0", {}).get("n_homologs_used"),
        "operator_v0": _version_block(V0_SPEC, homolog_input),
        "operator_v26": _version_block(V26_SPEC, homolog_input),
    }


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--family", default="tetr",
                    help="Family key (matches families/<key>/) — default tetr")
    p.add_argument("--min-homologs", type=int, default=5,
                    help="Minimum n_homologs (default 5)")
    p.add_argument("--score-cutoff", type=float, default=75.0,
                    help="Consensus-score cutoff; an entry passes if V0 OR "
                         "V2.6 exceeds it (default 75)")
    p.add_argument("--output", type=Path,
                    help="Output JSON path (default: "
                         "algorithms/benchmark/datasets/"
                         "<family>_high_confidence_operators.json)")
    p.add_argument("--workers", type=int, default=8,
                    help="Threads to run V0 + V2.6 operator_fetch inline "
                         "over qualifying entries (default 8)")
    args = p.parse_args()

    root = Path(__file__).resolve().parent.parent
    fam_dir = root / "families" / args.family
    jsonl_path = fam_dir / "members_predictions.jsonl"
    fasta_path = fam_dir / "members.fasta"
    refseq_map_path = fam_dir / "refseq_to_uniprot.json"
    cluster_homologs_path = fam_dir / "cluster_homologs.json"

    out_path = args.output or (
        root / "algorithms" / "benchmark" / "datasets"
        / f"{args.family}_high_confidence_operators.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading {fasta_path} …", flush=True)
    sequences = load_fasta(fasta_path)
    print(f"  {len(sequences):,} sequences indexed", flush=True)

    print(f"Loading {refseq_map_path} …", flush=True)
    uid_to_refseq = invert_refseq_map(refseq_map_path)
    print(f"  {len(uid_to_refseq):,} UniProt → RefSeq mappings", flush=True)

    print(f"Loading {cluster_homologs_path} …", flush=True)
    cluster_homologs = json.loads(cluster_homologs_path.read_text())
    print(f"  {len(cluster_homologs):,} clusters", flush=True)

    print(f"Scanning {jsonl_path} …", flush=True)
    records: dict[str, dict] = {}
    qualifying: list[dict] = []
    for rec in load_latest_records(jsonl_path):
        uid = rec.get("uniprot_id")
        if uid:
            records[uid] = rec
        if passes_filters(rec, args.min_homologs, args.score_cutoff):
            qualifying.append(rec)
    print(f"  {len(qualifying):,} entries pass filters "
           f"(n_homologs ≥ {args.min_homologs}, "
           f"consensus > {args.score_cutoff})")

    print(f"Running V0 + V2.6 operator_fetch inline with {args.workers} workers "
           f"(reconstructs frequency_matrix + aligned_operators) …", flush=True)
    entries: list[dict] = []
    lock = threading.Lock()
    completed = 0
    total = len(qualifying)

    def _one(rec: dict) -> dict:
        cm = cluster_homologs.get(rec["uniprot_id"]) or []
        return build_entry(rec, sequences, uid_to_refseq, cm, records)

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(_one, r): r for r in qualifying}
        for fut in as_completed(futures):
            try:
                entry = fut.result()
            except Exception as exc:
                print(f"  worker raised: {exc}", flush=True)
                continue
            with lock:
                entries.append(entry)
                completed += 1
                if completed % 500 == 0 or completed == total:
                    print(f"  {completed:,}/{total:,}", flush=True)

    missing_sequence = sum(1 for e in entries if e["protein_sequence"] is None)
    if missing_sequence:
        print(f"  ⚠  {missing_sequence:,} entries had no sequence in members.fasta")

    entries.sort(key=lambda e: e["uniprot_id"])

    # Compact JSON — this file is a data payload consumed by the frontend
    # Vite build; pretty-printing was doubling its size.
    out_path.write_text(json.dumps(entries, separators=(",", ":")))
    size_kb = out_path.stat().st_size / 1024
    print(f"Wrote {len(entries):,} entries to {out_path} "
           f"({size_kb:,.0f} KB)")


if __name__ == "__main__":
    main()
