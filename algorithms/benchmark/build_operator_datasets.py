"""Build two operator datasets for evaluating candidate-scoring functions.

POSITIVES — every literature-validated operator from the benchmark dataset.
  Source: algorithms/benchmark/datasets/tetr_137.json
  Output: algorithms/benchmark/datasets/operator_positives.json

NEGATIVES — predicted-motif sequences from a previous benchmark run that
  diverged from the known operator (algorithm output identity < 50%).
  These are the "spurious" palindromes our current operator_fetch picks
  when it gets the wrong answer — i.e. they look palindromic enough to
  win on the current consensus_score metric but they aren't real operators.
  Source: any benchmark run's per_protein.json
  Output: algorithms/benchmark/datasets/operator_negatives.json

A good scoring function will score positives higher than negatives. The
score_classifier.py harness measures that.

Each record is a dict with at minimum:
  {"sequence": str, "ncbi_accession": str, "alias": str, "source": str}

Negatives also carry:
  {"v0_identity_to_known": float, "v0_consensus_score": float}

Usage:
  python -m algorithms.benchmark.build_operator_datasets
  python -m algorithms.benchmark.build_operator_datasets \
      --benchmark-run algorithms/benchmark/results/2026-05-22_14-55-10/
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATASETS_DIR = Path(__file__).resolve().parent / "datasets"
DEFAULT_BENCHMARK = (DATASETS_DIR / "tetr_137.json")
DEFAULT_RESULTS_DIR = (Path(__file__).resolve().parent / "results"
                       / "2026-05-22_14-55-10")
POSITIVES_OUT = DATASETS_DIR / "operator_positives.json"
NEGATIVES_OUT = DATASETS_DIR / "operator_negatives.json"

# Predicted motifs whose alignment to the known operator is below this
# percent identity count as "negative controls" — the algorithm picked
# something that's clearly not the real operator.
NEGATIVE_IDENTITY_THRESHOLD = 50.0


def _clean(seq: str | None) -> str | None:
    """Keep only ACGT (case preserved). Returns None if no valid chars."""
    if not seq or not isinstance(seq, str):
        return None
    out = re.sub(r"[^ACGTacgt]", "", seq)
    return out or None


def build_positives(benchmark_path: Path) -> list[dict]:
    """One row per known operator (proteins with multiple operators contribute
    multiple rows). Same sequence appearing on multiple proteins is kept
    once per protein for accurate provenance."""
    raw = json.loads(benchmark_path.read_text())
    out: list[dict] = []
    for entry in raw:
        acc = entry.get("ncbi_accession")
        alias = entry.get("alias")
        for op in entry.get("operators") or []:
            seq = _clean(op.get("sequence"))
            if not seq:
                continue
            out.append({
                "sequence": seq,
                "ncbi_accession": acc,
                "alias": alias,
                "doi": op.get("doi"),
                "source": op.get("source"),
            })
    return out


def build_negatives(results_dir: Path,
                     identity_threshold: float = NEGATIVE_IDENTITY_THRESHOLD,
                     version: str = "v0") -> list[dict]:
    """Predicted motifs from a benchmark run where the algorithm's output
    is clearly wrong (low identity to the known operator). These are the
    spurious palindromes our current operator_fetch picks."""
    pp = json.loads((results_dir / "per_protein.json").read_text())
    out: list[dict] = []
    for r in pp:
        if not r.get("homologs_fetched"):
            continue
        vblock = r.get("versions", {}).get(version, {})
        motif = vblock.get("predicted_motif")
        seq = _clean(motif)
        if not seq:
            continue
        identity = vblock.get("metrics", {}).get(
            "known_operator_in_predicted_motif", 0)
        if identity >= identity_threshold:
            continue
        out.append({
            "sequence": seq,
            "ncbi_accession": r["ncbi_accession"],
            "alias": r.get("alias"),
            "source": f"{version}_predicted_motif",
            "v0_identity_to_known": identity,
            "v0_consensus_score": vblock.get("consensus_score"),
        })
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--benchmark", type=Path, default=DEFAULT_BENCHMARK,
                    help="Benchmark dataset JSON (positives source)")
    p.add_argument("--benchmark-run", type=Path, default=DEFAULT_RESULTS_DIR,
                    help="Per-protein benchmark results dir (negatives source)")
    p.add_argument("--version", default="v0",
                    help="Which version's predicted_motifs to source negatives from")
    p.add_argument("--identity-threshold", type=float,
                    default=NEGATIVE_IDENTITY_THRESHOLD,
                    help="Predicted motifs below this identity to the known "
                         "operator count as negatives")
    args = p.parse_args()

    DATASETS_DIR.mkdir(parents=True, exist_ok=True)

    pos = build_positives(args.benchmark)
    neg = build_negatives(args.benchmark_run, args.identity_threshold,
                          args.version)

    POSITIVES_OUT.write_text(json.dumps(pos, indent=2))
    NEGATIVES_OUT.write_text(json.dumps(neg, indent=2))

    pos_lens = [len(p["sequence"]) for p in pos]
    neg_lens = [len(n["sequence"]) for n in neg]
    print(f"Positives: {len(pos)} sequences from {args.benchmark.name}")
    print(f"  length range: {min(pos_lens)}-{max(pos_lens)} bp (median "
          f"{sorted(pos_lens)[len(pos_lens)//2]})")
    print(f"Negatives: {len(neg)} sequences from {args.benchmark_run.name} "
          f"({args.version}, identity < {args.identity_threshold}%)")
    print(f"  length range: {min(neg_lens)}-{max(neg_lens)} bp (median "
          f"{sorted(neg_lens)[len(neg_lens)//2]})")
    print(f"\nSaved positives → {POSITIVES_OUT}")
    print(f"Saved negatives → {NEGATIVES_OUT}")


if __name__ == "__main__":
    main()
