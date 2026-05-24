"""Evaluate candidate-scoring functions on operator positive/negative sets.

For every scorer in operator_scorers.SCORERS:
  * Compute a score for each positive (known operator from tetr_137.json)
    and each negative (spurious predicted motif from a benchmark run).
  * Discrimination statistics — ROC AUC, Mann-Whitney U, Cohen's d, plus
    mean/median by class.
  * Rank scorers by AUC.

Outputs two files:
  * results/<timestamp>/scorers_report.md       — full detail per scorer
  * algorithms/benchmark/LATEST_SCORERS.md      — committed compact summary

Usage:
  python -m algorithms.benchmark.score_classifier
  python -m algorithms.benchmark.score_classifier \
      --positives algorithms/benchmark/datasets/operator_positives.json \
      --negatives algorithms/benchmark/datasets/operator_negatives.json
"""

from __future__ import annotations

import argparse
import datetime
import json
import math
import statistics
import sys
from pathlib import Path
from typing import Callable

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from algorithms.benchmark.operator_scorers import SCORERS  # noqa: E402

DATASETS_DIR = Path(__file__).resolve().parent / "datasets"
RESULTS_ROOT = Path(__file__).resolve().parent / "results"
LATEST_PATH  = Path(__file__).resolve().parent / "LATEST_SCORERS.md"
DEFAULT_POSITIVES = DATASETS_DIR / "operator_positives.json"
DEFAULT_NEGATIVES = DATASETS_DIR / "operator_negatives.json"


# --- statistics ----------------------------------------------------------

def rank_scores(scores: list[float]) -> list[float]:
    """Return average-rank ordinal positions (1-indexed). Ties get the
    average of their tied positions."""
    indexed = sorted(enumerate(scores), key=lambda x: x[1])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(indexed):
        j = i
        while j + 1 < len(indexed) and indexed[j + 1][1] == indexed[i][1]:
            j += 1
        avg_rank = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[indexed[k][0]] = avg_rank
        i = j + 1
    return ranks


def auc_from_ranks(pos_scores: list[float], neg_scores: list[float]) -> float:
    """ROC AUC via the Mann-Whitney U identity:
        AUC = (U_pos) / (n_pos × n_neg)
    where U_pos = R_pos - n_pos(n_pos+1)/2 and R_pos is the rank sum of
    the positive class."""
    n_pos = len(pos_scores)
    n_neg = len(neg_scores)
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    all_scores = pos_scores + neg_scores
    ranks = rank_scores(all_scores)
    r_pos = sum(ranks[:n_pos])
    u_pos = r_pos - n_pos * (n_pos + 1) / 2
    return u_pos / (n_pos * n_neg)


def cohens_d(pos_scores: list[float], neg_scores: list[float]) -> float:
    """Standardised mean difference using pooled SD. Higher = larger
    effect size in the positives > negatives direction."""
    if len(pos_scores) < 2 or len(neg_scores) < 2:
        return float("nan")
    mp = statistics.mean(pos_scores)
    mn = statistics.mean(neg_scores)
    sp = statistics.stdev(pos_scores)
    sn = statistics.stdev(neg_scores)
    pooled = math.sqrt(((len(pos_scores) - 1) * sp ** 2
                          + (len(neg_scores) - 1) * sn ** 2)
                         / (len(pos_scores) + len(neg_scores) - 2))
    if pooled == 0:
        return float("nan")
    return (mp - mn) / pooled


def threshold_separation(pos_scores: list[float],
                          neg_scores: list[float]) -> dict:
    """Find the score threshold that maximises Youden's J = TPR - FPR.
    Returns the threshold, TPR/FPR/precision at that threshold, and the
    fraction of positives correctly classified at that threshold."""
    candidates = sorted(set(pos_scores + neg_scores))
    best = None
    n_pos = len(pos_scores)
    n_neg = len(neg_scores)
    for t in candidates:
        tp = sum(1 for s in pos_scores if s >= t)
        fp = sum(1 for s in neg_scores if s >= t)
        if tp == 0 and fp == 0:
            continue
        tpr = tp / n_pos
        fpr = fp / n_neg
        j = tpr - fpr
        precision = tp / max(1, tp + fp)
        if best is None or j > best["youden_j"]:
            best = {"threshold": t, "tpr": tpr, "fpr": fpr,
                     "precision": precision, "youden_j": j,
                     "tp": tp, "fp": fp}
    return best or {"threshold": 0, "tpr": 0, "fpr": 0,
                     "precision": 0, "youden_j": 0, "tp": 0, "fp": 0}


# --- driver --------------------------------------------------------------

def evaluate_scorer(name: str, fn: Callable[[str], float],
                     positives: list[dict], negatives: list[dict]) -> dict:
    pos_scores = [fn(p["sequence"]) for p in positives]
    neg_scores = [fn(n["sequence"]) for n in negatives]
    return {
        "name": name,
        "n_pos": len(pos_scores),
        "n_neg": len(neg_scores),
        "auc": auc_from_ranks(pos_scores, neg_scores),
        "cohens_d": cohens_d(pos_scores, neg_scores),
        "pos_mean": statistics.mean(pos_scores) if pos_scores else 0.0,
        "pos_median": statistics.median(pos_scores) if pos_scores else 0.0,
        "neg_mean": statistics.mean(neg_scores) if neg_scores else 0.0,
        "neg_median": statistics.median(neg_scores) if neg_scores else 0.0,
        "threshold": threshold_separation(pos_scores, neg_scores),
        "pos_scores": pos_scores,
        "neg_scores": neg_scores,
    }


def render_markdown(results: list[dict], pos_meta: dict, neg_meta: dict,
                     run_info: dict) -> str:
    lines: list[str] = []
    lines.append("# Operator-scorer classifier benchmark")
    lines.append("")
    lines.append(f"- **Run**: `{run_info['timestamp']}`")
    lines.append(f"- **Positives**: {pos_meta['count']} known operators from "
                 f"`{pos_meta['source']}`")
    lines.append(f"- **Negatives**: {neg_meta['count']} spurious predicted "
                 f"motifs from `{neg_meta['source']}` "
                 f"(identity to known < {neg_meta['threshold']}%)")
    lines.append("")
    lines.append("Each scorer returns a single float per sequence. Higher "
                 "should mean more operator-like. We rank by ROC AUC — "
                 "AUC = 0.5 is random, 1.0 is perfect separation.")
    lines.append("")

    # Ranked summary table
    lines.append("## Ranking (by ROC AUC)")
    lines.append("")
    lines.append("| Rank | Scorer | AUC | Cohen's d | Best threshold | TP@thr | FP@thr | Precision |")
    lines.append("|---:|---|---:|---:|---:|---:|---:|---:|")
    ranked = sorted(results, key=lambda r: r["auc"], reverse=True)
    for i, r in enumerate(ranked, start=1):
        thr = r["threshold"]
        lines.append(
            f"| {i} | `{r['name']}` | "
            f"{r['auc']:.3f} | {r['cohens_d']:.2f} | "
            f"{thr['threshold']:.3f} | "
            f"{thr['tp']}/{r['n_pos']} ({100*thr['tpr']:.0f}%) | "
            f"{thr['fp']}/{r['n_neg']} ({100*thr['fpr']:.0f}%) | "
            f"{100*thr['precision']:.0f}% |"
        )
    lines.append("")

    # Per-scorer detail
    lines.append("## Per-scorer detail")
    lines.append("")
    lines.append("| Scorer | Positives mean (median) | Negatives mean (median) | Δ mean |")
    lines.append("|---|---:|---:|---:|")
    for r in ranked:
        lines.append(
            f"| `{r['name']}` | "
            f"{r['pos_mean']:.3f} ({r['pos_median']:.3f}) | "
            f"{r['neg_mean']:.3f} ({r['neg_median']:.3f}) | "
            f"{r['pos_mean'] - r['neg_mean']:+.3f} |"
        )
    lines.append("")

    # Top scorer textual note
    if ranked:
        top = ranked[0]
        lines.append(f"## Top scorer: `{top['name']}`")
        lines.append("")
        lines.append(
            f"AUC = **{top['auc']:.3f}**, Cohen's d = **{top['cohens_d']:.2f}**.  "
            f"At the optimal threshold (`{top['threshold']['threshold']:.3f}`) "
            f"it captures **{top['threshold']['tp']}/{top['n_pos']} "
            f"({100*top['threshold']['tpr']:.0f}%)** of positives while only "
            f"letting through **{top['threshold']['fp']}/{top['n_neg']} "
            f"({100*top['threshold']['fpr']:.0f}%)** of negatives. "
            f"Precision at this threshold: **{100*top['threshold']['precision']:.0f}%**."
        )
        lines.append("")

    return "\n".join(lines)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--positives", type=Path, default=DEFAULT_POSITIVES)
    p.add_argument("--negatives", type=Path, default=DEFAULT_NEGATIVES)
    p.add_argument("--negatives-threshold", type=float, default=50.0,
                    help="Identity threshold used to build negatives (for the report)")
    p.add_argument("--output", type=Path, default=None,
                    help="Output directory (default: results/<timestamp>)")
    args = p.parse_args()

    positives = json.loads(args.positives.read_text())
    negatives = json.loads(args.negatives.read_text())

    results = [evaluate_scorer(name, fn, positives, negatives)
                for name, fn, _ in SCORERS]

    ts = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    out_dir = args.output or (RESULTS_ROOT / f"scorers_{ts}")
    out_dir.mkdir(parents=True, exist_ok=True)

    run_info = {"timestamp": ts,
                 "n_scorers": len(SCORERS),
                 "positives_path": str(args.positives),
                 "negatives_path": str(args.negatives)}
    pos_meta = {"count": len(positives), "source": args.positives.name}
    neg_meta = {"count": len(negatives), "source": args.negatives.name,
                 "threshold": args.negatives_threshold}

    # Save raw per-scorer scores for ad-hoc analysis
    raw = {"run_info": run_info, "positives_meta": pos_meta,
            "negatives_meta": neg_meta,
            "scorers": [{k: v for k, v in r.items() if k != "raw"}
                         for r in results]}
    (out_dir / "scorers.json").write_text(json.dumps(raw, indent=2))

    md = render_markdown(results, pos_meta, neg_meta, run_info)
    (out_dir / "scorers_report.md").write_text(md)
    LATEST_PATH.write_text(md)

    print(f"Wrote scorers.json + scorers_report.md to {out_dir}")
    print(f"Wrote LATEST_SCORERS.md to {LATEST_PATH.parent}")

    # Brief stdout ranking
    ranked = sorted(results, key=lambda r: r["auc"], reverse=True)
    print(f"\nRanking by AUC ({len(positives)} positives vs {len(negatives)} negatives):")
    for i, r in enumerate(ranked, start=1):
        print(f"  {i:>2}. {r['name']:<35} AUC={r['auc']:.3f}  d={r['cohens_d']:+.2f}")


if __name__ == "__main__":
    main()
