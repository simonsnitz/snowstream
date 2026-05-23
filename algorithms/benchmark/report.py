"""Render a benchmark results directory into a markdown report + PNG charts.

Reads `per_protein.json`, `summary.json`, `run_info.json` from a results
directory and writes:
  * `report.md` — detailed comparison tables + chart references (lives in
    the gitignored results/ dir alongside the PNGs)
  * `chart_<metric>.png` — distribution plots:
      - identity metrics: cumulative count above each threshold, per version
      - count metrics: histogram of values per version
  * `algorithms/benchmark/LATEST.md` — committed, text-only summary that
    overwrites on every run so the most recent result is always visible
    in VSCode/GitHub without spelunking into gitignored dirs

Usage:
  python -m algorithms.benchmark.report path/to/results/dir
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .metrics import IDENTITY_THRESHOLDS


def _fmt_pct(n: int, total: int) -> str:
    if total == 0:
        return "-"
    return f"{n} ({100*n/total:.1f}%)"


def _emit_summary_markdown(summary: dict, run_info: dict,
                            per_protein: list[dict]) -> str:
    lines: list[str] = []
    versions = list(summary["versions"].keys())
    n_total = run_info["dataset_size"]
    n_resolved = sum(1 for r in per_protein if r.get("homologs_fetched"))

    lines.append(f"# Benchmark report — {run_info['timestamp']}")
    lines.append("")
    lines.append(f"- **Dataset**: `{run_info['dataset_path']}`  "
                 f"(N = {n_total}; {n_resolved} resolved via snowstream)")
    lines.append(f"- **Versions compared**: " + ", ".join(f"`{v}`" for v in versions))
    if run_info.get("max_proteins"):
        lines.append(f"- **Smoke-test limit**: first {run_info['max_proteins']} proteins")
    lines.append("")

    lines.append("## Identity-metric counts at each threshold")
    lines.append("")
    lines.append("For each identity metric, the cell shows the count of proteins "
                 "whose best alignment score met or exceeded the threshold.")
    lines.append("")

    # One table per identity metric — versions as columns, thresholds as rows
    for mdef in summary["metric_definitions"]:
        if mdef["kind"] != "identity":
            continue
        name = mdef["name"]
        label = mdef["label"]
        any_block = next(iter(summary["versions"].values()))["metrics"][name]
        denom = any_block["n"]
        lines.append(f"### {label}  (n = {denom})")
        lines.append("")
        header = "| Threshold | " + " | ".join(versions) + " |"
        sep    = "|---|" + "|".join(["---:"] * len(versions)) + "|"
        lines.append(header)
        lines.append(sep)
        for t in IDENTITY_THRESHOLDS:
            row = [f"≥ {t}%"]
            for v in versions:
                c = summary["versions"][v]["metrics"][name][f"count_ge_{t}"]
                row.append(_fmt_pct(c, denom))
            lines.append("| " + " | ".join(row) + " |")
        # Mean line
        mean_row = ["mean"]
        for v in versions:
            mean_row.append(f"{summary['versions'][v]['metrics'][name]['mean']:.1f}%")
        lines.append("| " + " | ".join(mean_row) + " |")
        lines.append("")

    lines.append("## Count metrics")
    lines.append("")
    for mdef in summary["metric_definitions"]:
        if mdef["kind"] != "count":
            continue
        name = mdef["name"]
        label = mdef["label"]
        lines.append(f"### {label}")
        lines.append("")
        header = "| Statistic | " + " | ".join(versions) + " |"
        sep    = "|---|" + "|".join(["---:"] * len(versions)) + "|"
        lines.append(header)
        lines.append(sep)
        for stat in ("mean", "min", "max"):
            row = [stat]
            for v in versions:
                row.append(str(summary["versions"][v]["metrics"][name][stat]))
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")

    # Per-protein head-to-head: only show rows where versions differ
    if len(versions) >= 2:
        lines.append("## Per-protein head-to-head (only proteins where versions disagree)")
        lines.append("")
        algo_key = "known_operator_in_predicted_motif"
        rows: list[tuple] = []
        for r in per_protein:
            if not r.get("homologs_fetched"):
                continue
            vmetrics = r["versions"]
            scores = []
            for v in versions:
                vm = vmetrics.get(v, {}).get("metrics", {})
                scores.append(vm.get(algo_key, 0))
            if max(scores) - min(scores) > 1:
                rows.append((r["ncbi_accession"], r.get("alias") or "", scores))

        if rows:
            # Δ = latest version − earliest version (positive = improvement)
            last_v_idx = len(versions) - 1
            header = ("| NCBI | Alias | " + " | ".join(versions)
                       + f" | Δ ({versions[last_v_idx]}−{versions[0]}) |")
            sep    = "|---|---|" + "|".join(["---:"] * len(versions)) + "|---:|"
            lines.append(header)
            lines.append(sep)
            # Sort by Δ (largest gain first; regressions at the bottom)
            def _delta(row):
                _, _, s = row
                return s[last_v_idx] - s[0]
            rows.sort(key=_delta, reverse=True)
            for acc, alias, scores in rows:
                cells = ["{:.1f}%".format(s) for s in scores]
                delta = scores[last_v_idx] - scores[0]
                lines.append("| " + " | ".join(
                    [acc, alias] + cells + [f"{delta:+.1f}%"]) + " |")
            lines.append("")
            # Aggregate counts of wins/losses/ties between first and last version
            wins = sum(1 for _, _, s in rows if s[last_v_idx] - s[0] > 1)
            losses = sum(1 for _, _, s in rows if s[last_v_idx] - s[0] < -1)
            lines.append(f"**Net direction**: {versions[last_v_idx]} wins on "
                         f"{wins}, loses on {losses} (out of {len(rows)} "
                         f"proteins with disagreement > 1%).")
            lines.append("")
        else:
            lines.append("_(no per-protein differences > 1%; all versions agree)_")
            lines.append("")

    lines.append("## Charts")
    lines.append("")
    for mdef in summary["metric_definitions"]:
        png = f"chart_{mdef['name']}.png"
        lines.append(f"- `{mdef['label']}` — see `{png}`")
    lines.append("")

    return "\n".join(lines)


def _emit_charts(summary: dict, out_dir: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"  (matplotlib unavailable: {e}) — skipping PNG charts")
        return

    versions = list(summary["versions"].keys())
    colors = plt.cm.viridis([i / max(1, len(versions) - 1)
                               for i in range(len(versions))])

    for mdef in summary["metric_definitions"]:
        name = mdef["name"]
        label = mdef["label"]
        if mdef["kind"] == "identity":
            fig, ax = plt.subplots(figsize=(7, 4))
            xs = list(IDENTITY_THRESHOLDS)
            width = 0.8 / len(versions)
            for i, v in enumerate(versions):
                counts = [summary["versions"][v]["metrics"][name][f"count_ge_{t}"]
                          for t in xs]
                positions = [x + (i - (len(versions) - 1) / 2) * width
                             for x in xs]
                ax.bar(positions, counts, width=width, label=v, color=colors[i])
            ax.set_xticks(xs)
            ax.set_xticklabels([f"≥ {t}%" for t in xs])
            ax.set_ylabel("# proteins")
            ax.set_title(label)
            ax.legend(loc="upper right")
            ax.grid(True, axis="y", alpha=0.3)
            fig.tight_layout()
            fig.savefig(out_dir / f"chart_{name}.png", dpi=120)
            plt.close(fig)
        elif mdef["kind"] == "count":
            fig, ax = plt.subplots(figsize=(7, 4))
            values_per_v = [summary["versions"][v]["metrics"][name]["values"]
                             for v in versions]
            all_vals = [x for vs in values_per_v for x in vs]
            if all_vals:
                lo, hi = min(all_vals), max(all_vals)
                bins = max(10, min(30, hi - lo + 1))
            else:
                bins = 10
            for i, (v, vals) in enumerate(zip(versions, values_per_v)):
                ax.hist(vals, bins=bins, alpha=0.5, label=v, color=colors[i])
            ax.set_xlabel(label)
            ax.set_ylabel("# proteins")
            ax.set_title(label)
            ax.legend()
            ax.grid(True, axis="y", alpha=0.3)
            fig.tight_layout()
            fig.savefig(out_dir / f"chart_{name}.png", dpi=120)
            plt.close(fig)


def _emit_latest_markdown(summary: dict, run_info: dict,
                           per_protein: list[dict]) -> str:
    """Compact, text-only summary suitable for committing as LATEST.md.

    Same threshold-count tables as the full report (these are already small),
    drops the "## Charts" section, and shows the top 10 wins + top 10 losses
    from the head-to-head instead of the full sorted list. Designed to stay
    under ~5 KB regardless of dataset size.
    """
    lines: list[str] = []
    versions = list(summary["versions"].keys())
    n_total = run_info["dataset_size"]
    n_resolved = sum(1 for r in per_protein if r.get("homologs_fetched"))

    dataset_name = Path(run_info["dataset_path"]).name
    lines.append("# Latest benchmark result")
    lines.append("")
    lines.append("> Auto-generated on every `algorithms.benchmark.run`. "
                 "For per-run history (with PNG charts) see "
                 "`algorithms/benchmark/results/<timestamp>/`.")
    lines.append("")
    lines.append(f"- **Run**: `{run_info['timestamp']}`")
    lines.append(f"- **Dataset**: `{dataset_name}`  "
                 f"(N = {n_total}; {n_resolved} resolved)")
    lines.append(f"- **Versions**: " + ", ".join(f"`{v}`" for v in versions))
    if run_info.get("max_proteins"):
        lines.append(f"- **Smoke-test limit**: first {run_info['max_proteins']} proteins")
    lines.append("")

    # Identity metrics — same threshold tables as the full report
    lines.append("## Identity-metric threshold counts")
    lines.append("")
    for mdef in summary["metric_definitions"]:
        if mdef["kind"] != "identity":
            continue
        name = mdef["name"]
        label = mdef["label"]
        any_block = next(iter(summary["versions"].values()))["metrics"][name]
        denom = any_block["n"]
        lines.append(f"### {label}  (n = {denom})")
        lines.append("")
        header = "| Threshold | " + " | ".join(versions) + " |"
        sep    = "|---|" + "|".join(["---:"] * len(versions)) + "|"
        lines.append(header)
        lines.append(sep)
        for t in IDENTITY_THRESHOLDS:
            row = [f"≥ {t}%"]
            for v in versions:
                c = summary["versions"][v]["metrics"][name][f"count_ge_{t}"]
                row.append(_fmt_pct(c, denom))
            lines.append("| " + " | ".join(row) + " |")
        mean_row = ["mean"]
        for v in versions:
            mean_row.append(f"{summary['versions'][v]['metrics'][name]['mean']:.1f}%")
        lines.append("| " + " | ".join(mean_row) + " |")
        lines.append("")

    # Count metrics — compact one-liner per metric
    lines.append("## Count metrics")
    lines.append("")
    header = "| Metric | " + " | ".join(versions) + " |"
    sep    = "|---|" + "|".join(["---:"] * len(versions)) + "|"
    lines.append(header)
    lines.append(sep)
    for mdef in summary["metric_definitions"]:
        if mdef["kind"] != "count":
            continue
        name = mdef["name"]
        label = mdef["label"]
        row = [f"{label} (mean)"]
        for v in versions:
            row.append(f"{summary['versions'][v]['metrics'][name]['mean']:.1f}")
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")

    # Head-to-head (capped at top 10 wins + top 10 losses)
    if len(versions) >= 2:
        last_v = versions[-1]
        first_v = versions[0]
        algo_key = "known_operator_in_predicted_motif"
        rows: list[tuple] = []
        for r in per_protein:
            if not r.get("homologs_fetched"):
                continue
            vm = r["versions"]
            scores = [vm.get(v, {}).get("metrics", {}).get(algo_key, 0)
                      for v in versions]
            if max(scores) - min(scores) > 1:
                rows.append((r["ncbi_accession"], r.get("alias") or "",
                              scores))

        if rows:
            lines.append(f"## Head-to-head: {last_v} − {first_v} "
                         f"(predicted-motif identity to known operator)")
            lines.append("")
            rows.sort(key=lambda x: x[2][-1] - x[2][0], reverse=True)
            wins = [r for r in rows if r[2][-1] - r[2][0] > 1]
            losses = [r for r in rows if r[2][-1] - r[2][0] < -1]
            cap = 10

            def _table(subset, title):
                if not subset:
                    return
                lines.append(f"### {title} (showing {min(len(subset), cap)}"
                             f" of {len(subset)})")
                lines.append("")
                header = "| NCBI | Alias | " + " | ".join(versions) + f" | Δ |"
                sep    = "|---|---|" + "|".join(["---:"] * len(versions)) + "|---:|"
                lines.append(header)
                lines.append(sep)
                for acc, alias, scores in subset[:cap]:
                    cells = [f"{s:.1f}%" for s in scores]
                    delta = scores[-1] - scores[0]
                    lines.append("| " + " | ".join(
                        [acc, alias] + cells + [f"{delta:+.1f}%"]) + " |")
                lines.append("")

            _table(wins, f"{last_v} wins")
            _table(list(reversed(losses)), f"{last_v} losses")

            lines.append(f"**Net direction**: {last_v} wins on {len(wins)}, "
                         f"loses on {len(losses)} (out of "
                         f"{len(rows)} proteins where versions disagree by "
                         f">1%).")
            lines.append("")
        else:
            lines.append("_(all versions agree within 1% on every protein)_")
            lines.append("")

    return "\n".join(lines)


# Path of the committed-to-repo summary, written on every run.
LATEST_PATH = Path(__file__).resolve().parent / "LATEST.md"


def render(results_dir: Path) -> None:
    run_info = json.loads((results_dir / "run_info.json").read_text())
    summary = json.loads((results_dir / "summary.json").read_text())
    per_protein = json.loads((results_dir / "per_protein.json").read_text())

    # Per-run detailed report (with chart references) lives in the
    # gitignored results/ dir alongside the PNGs.
    md = _emit_summary_markdown(summary, run_info, per_protein)
    (results_dir / "report.md").write_text(md)

    _emit_charts(summary, results_dir)

    # Committed text-only summary — overwritten on every run.
    latest_md = _emit_latest_markdown(summary, run_info, per_protein)
    LATEST_PATH.write_text(latest_md)

    print(f"Wrote report.md and chart_*.png to {results_dir}")
    print(f"Wrote {LATEST_PATH.name} to {LATEST_PATH.parent}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("results_dir", type=Path,
                    help="Directory containing per_protein.json, summary.json, run_info.json")
    args = p.parse_args()
    render(args.results_dir)


if __name__ == "__main__":
    main()
