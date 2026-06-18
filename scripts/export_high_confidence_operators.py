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
  * operator_v0           — { motif, consensus_score, native_for_centroid }
  * operator_v26          — { motif, consensus_score, rerank_score,
                              native_for_centroid }

The output JSON is a list, sorted by uniprot_id for stable diffs across
re-runs. Default output path is
`algorithms/benchmark/datasets/tetr_high_confidence_operators.json`
so it sits alongside the other benchmark datasets.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterator


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


def build_entry(rec: dict,
                 sequences: dict[str, str],
                 uid_to_refseq: dict[str, str]) -> dict:
    uid = rec["uniprot_id"]
    v0 = rec["operator_v0"]
    v26 = rec["operator_v26"]
    return {
        "uniprot_id": uid,
        "ncbi_accession": uid_to_refseq.get(uid),
        "protein_sequence": sequences.get(uid),
        "genome": rec.get("genome"),
        "n_homologs": v0.get("n_homologs_used"),
        "operator_v0": {
            "motif": v0.get("motif"),
            "consensus_score": v0.get("consensus_score"),
            "native_for_centroid": v0.get("native_for_centroid"),
        },
        "operator_v26": {
            "motif": v26.get("motif"),
            "consensus_score": v26.get("consensus_score"),
            "rerank_score": v26.get("rerank_score"),
            "native_for_centroid": v26.get("native_for_centroid"),
        },
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
    args = p.parse_args()

    root = Path(__file__).resolve().parent.parent
    fam_dir = root / "families" / args.family
    jsonl_path = fam_dir / "members_predictions.jsonl"
    fasta_path = fam_dir / "members.fasta"
    refseq_map_path = fam_dir / "refseq_to_uniprot.json"

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

    print(f"Scanning {jsonl_path} …", flush=True)
    entries: list[dict] = []
    missing_sequence = 0
    for rec in load_latest_records(jsonl_path):
        if not passes_filters(rec, args.min_homologs, args.score_cutoff):
            continue
        entry = build_entry(rec, sequences, uid_to_refseq)
        if entry["protein_sequence"] is None:
            missing_sequence += 1
        entries.append(entry)

    entries.sort(key=lambda e: e["uniprot_id"])
    print(f"  {len(entries):,} entries pass filters "
           f"(n_homologs ≥ {args.min_homologs}, "
           f"consensus > {args.score_cutoff})")
    if missing_sequence:
        print(f"  ⚠  {missing_sequence:,} entries had no sequence in members.fasta")

    out_path.write_text(json.dumps(entries, indent=2))
    size_kb = out_path.stat().st_size / 1024
    print(f"Wrote {len(entries):,} entries to {out_path} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
