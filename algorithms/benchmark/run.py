"""Benchmark runner — entrypoint for the algorithm comparison harness.

Usage:
  python -m algorithms.benchmark.run                 # all versions, default dataset
  python -m algorithms.benchmark.run --versions v0 v1
  python -m algorithms.benchmark.run --dataset path/to/file.json
  python -m algorithms.benchmark.run --max-proteins 10

What it does:
  1. Loads the dataset (default: benchmark/datasets/tetr_137.json).
  2. For each protein, queries snowstream's /api/predict once to grab the
     smart-lookup homolog set (each homolog already has its operon cached).
  3. For each requested version, runs the version-aware pipeline on that
     fixed homolog set, computes every metric in metrics.METRICS, and
     records the per-protein result.
  4. Writes results to algorithms/benchmark/results/<timestamp>/per_protein.json
     plus a summary.json. The report.py script then renders markdown +
     PNG visuals.

Why we pull homologs from snowstream instead of running operon_fetch fresh:
the live members DB already has ~484k operons cached. Re-running operon_fetch
for every protein × every homolog would be tens of thousands of NCBI calls
per benchmark run. By default we reuse the cached operons and benchmark
only the parts of the pipeline that differ between versions. A future
operon_fetch version would need its own dedicated benchmark mode that
re-runs operon_fetch (a flag here can flip that).
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import sys
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from algorithms import pipeline  # noqa: E402
from algorithms.benchmark import metrics as metrics_mod  # noqa: E402
from algorithms.versions import VERSIONS, get_version  # noqa: E402

API_BASE = os.environ.get("SNOWSTREAM_API_BASE", "http://127.0.0.1:8000")
DEFAULT_DATASET = (Path(__file__).parent / "datasets" / "tetr_137.json")
RESULTS_ROOT   = Path(__file__).parent / "results"


# Same defaults the Streamlit app and deployed Snowprint use.
DEFAULT_PROMOTER_PARAMS = {
    "min_length": 80,
    "max_length": 800,
    "min_internal_length": 60,  # V1-only; ignored by V0
}
DEFAULT_OPERATOR_PARAMS = {
    "extension_length": 5,
    "win_score": 2,
    "loss_score": -2,
    "spacer_penalty": {str(i): v for i, v in zip(range(21),
        [4, 4, 4, 4, 4, 2, 2, 0, 0, -2, -2, -4, -4, -6, -6, -8, -8, -10, -10, -12, -12])},
    "gap_open": -100,
    "gap_extend": 0,
    "align_match": 2,
    "align_mismatch": -0.5,
    "min_operator_length": 5,
    "max_operator_length": 15,
    "seq_to_align": None,
    "search_method": "Look for inverted repeats",
}


# --- snowstream client (for getting homologs+operons) --------------------

def _query_snowstream(input_method: str, value: str,
                       timeout: float = 60) -> list[dict]:
    """POST to /api/predict and return the smart_lookup_hit homologs (each
    has uniprot_id, promoter, operon, protein_index, genome)."""
    try:
        payload = json.dumps({
            "input_method": input_method, "input_value": value,
            "database": "local_diamond",
        }).encode()
        req = urllib.request.Request(
            f"{API_BASE}/api/predict", data=payload,
            headers={"Content-Type": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = resp.read().decode("utf-8", errors="replace")
    except (urllib.error.HTTPError, urllib.error.URLError, OSError):
        return []
    except Exception as e:
        # IncompleteRead may carry a partial body
        partial = getattr(e, "partial", b"")
        body = partial.decode("utf-8", errors="replace") if partial else None
        if not body:
            return []
    for line in body.splitlines():
        if not line.startswith("data: "):
            continue
        try:
            ev = json.loads(line[6:])
        except json.JSONDecodeError:
            continue
        if ev.get("type") in ("smart_lookup_hit", "cached"):
            return ev.get("result", {}).get("homologs", []) or []
    return []


def _fetch_homologs(entry: dict) -> list[dict]:
    """Try RefSeq accession first, fall back to UniProt id."""
    acc = entry.get("ncbi_accession")
    uid = entry.get("uniprot_id")
    if acc:
        h = _query_snowstream("RefSeq", acc)
        if h:
            return h
    if uid:
        h = _query_snowstream("Uniprot", uid)
        if h:
            return h
    return []


# --- main loop -----------------------------------------------------------

def _process_one(entry: dict, version_names: list[str],
                  params: dict) -> dict:
    """Run the snowstream fetch + every requested version's pipeline for a
    single protein. Returns the per-protein record (same shape as before).
    Designed to be safe to call from a thread pool — no shared mutable
    state inside."""
    acc = entry.get("ncbi_accession")
    alias = entry.get("alias") or ""

    t0 = time.time()
    homologs = _fetch_homologs(entry)
    fetch_secs = round(time.time() - t0, 2)
    if not homologs:
        return {
            "ncbi_accession": acc, "alias": alias,
            "sources": entry.get("sources"),
            "uniprot_id": entry.get("uniprot_id"),
            "homologs_fetched": False,
            "versions": {},
            "_fetch_secs": fetch_secs,
        }

    record: dict = {
        "ncbi_accession": acc, "alias": alias,
        "sources": entry.get("sources"),
        "uniprot_id": entry.get("uniprot_id"),
        "homologs_fetched": True,
        "n_homologs_returned": len(homologs),
        "homolog_fetch_seconds": fetch_secs,
        "versions": {},
    }

    for vname in version_names:
        v = get_version(vname)
        t1 = time.time()
        try:
            out = pipeline.run(v, homologs, params)
        except Exception as e:
            out = {"error": f"{type(e).__name__}: {e}",
                    "operator_result": None,
                    "n_homologs": len(homologs),
                    "n_homologs_with_promoter": 0,
                    "query_promoter": None,
                    "all_homolog_promoters": [],
                    "all_query_candidates": []}
        run_secs = round(time.time() - t1, 2)
        m = metrics_mod.compute_all(entry, out)
        record["versions"][vname] = {
            "metrics": m,
            "run_seconds": run_secs,
            "n_query_candidates": len(out.get("all_query_candidates") or []),
            "n_homologs_with_promoter": out.get("n_homologs_with_promoter", 0),
            "predicted_motif": _render_motif(
                (out.get("operator_result") or {}).get("motif")),
            "consensus_score":
                (out.get("operator_result") or {}).get("consensus_score"),
            "error": out.get("error"),
        }
    return record


def run(dataset_path: Path, version_names: list[str],
        max_proteins: int | None, out_dir: Path, workers: int = 1) -> None:
    dataset = json.loads(dataset_path.read_text())
    if max_proteins:
        dataset = dataset[:max_proteins]

    promoter_params = DEFAULT_PROMOTER_PARAMS
    operator_params = DEFAULT_OPERATOR_PARAMS
    full_params = {**promoter_params, **operator_params}

    print(f"Dataset:   {dataset_path}  ({len(dataset)} proteins)", flush=True)
    print(f"Versions:  {', '.join(version_names)}", flush=True)
    print(f"Workers:   {workers}", flush=True)
    print(f"Output:    {out_dir}", flush=True)
    print()

    per_protein: list[dict] = [None] * len(dataset)  # type: ignore
    write_lock = threading.Lock()
    completed = 0
    total = len(dataset)

    def _emit(i: int, entry: dict, record: dict) -> None:
        nonlocal completed
        acc = entry.get("ncbi_accession") or ""
        alias = entry.get("alias") or ""
        with write_lock:
            per_protein[i] = record
            completed += 1
            if record.get("homologs_fetched"):
                scores = " ".join(
                    f"{v}={record['versions'][v]['metrics']['known_operator_in_predicted_motif']:.0f}%"
                    for v in version_names)
                print(f"[{completed:>3}/{total}] {acc:<18} {alias:<14} → algo: {scores}",
                      flush=True)
            else:
                print(f"[{completed:>3}/{total}] {acc:<18} {alias:<14} → no homologs",
                      flush=True)
            # Snapshot to disk so a crash doesn't lose progress
            (out_dir / "per_protein.json").write_text(
                json.dumps([r for r in per_protein if r is not None], indent=2))

    if workers <= 1:
        for i, entry in enumerate(dataset):
            record = _process_one(entry, version_names, full_params)
            _emit(i, entry, record)
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_process_one, entry, version_names, full_params): i
                        for i, entry in enumerate(dataset)}
            for fut in as_completed(futures):
                i = futures[fut]
                entry = dataset[i]
                try:
                    record = fut.result()
                except Exception as e:
                    record = {
                        "ncbi_accession": entry.get("ncbi_accession"),
                        "alias": entry.get("alias"),
                        "homologs_fetched": False,
                        "versions": {},
                        "_error": f"{type(e).__name__}: {e}",
                    }
                _emit(i, entry, record)

    per_protein_clean = [r for r in per_protein if r is not None]
    summary = _aggregate(per_protein_clean, version_names)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nWrote per_protein.json and summary.json to {out_dir}")


def _render_motif(motif) -> str | None:
    if not motif:
        return None
    if isinstance(motif, str):
        return motif
    try:
        return "".join(x["base"] for x in motif if isinstance(x, dict))
    except Exception:
        return None


def _aggregate(per_protein: list[dict], version_names: list[str]) -> dict:
    """Build the high-level summary the report renders. For each version,
    for each metric, capture (a) the raw values, (b) per-threshold counts
    for identity metrics, (c) means for count metrics."""
    out: dict = {"versions": {}, "metric_definitions": []}
    for name, _fn, kind, label in metrics_mod.METRICS:
        out["metric_definitions"].append(
            {"name": name, "kind": kind, "label": label})

    for v in version_names:
        version_block: dict = {"metrics": {}}
        for name, _fn, kind, label in metrics_mod.METRICS:
            values = [r["versions"][v]["metrics"][name]
                      for r in per_protein
                      if r.get("homologs_fetched") and v in r.get("versions", {})]
            entry: dict = {"label": label, "kind": kind, "n": len(values),
                            "values": values}
            if kind == "identity":
                for t in metrics_mod.IDENTITY_THRESHOLDS:
                    entry[f"count_ge_{t}"] = sum(
                        1 for x in values if x >= t)
                entry["mean"] = round(
                    sum(values) / len(values), 2) if values else 0
            elif kind == "count":
                entry["mean"] = round(
                    sum(values) / len(values), 2) if values else 0
                entry["min"] = min(values) if values else 0
                entry["max"] = max(values) if values else 0
            version_block["metrics"][name] = entry
        out["versions"][v] = version_block
    return out


# --- CLI -----------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--versions", nargs="+", default=list(VERSIONS.keys()),
                    help=f"Versions to benchmark (default: all). Known: "
                         f"{', '.join(VERSIONS.keys())}")
    p.add_argument("--dataset", type=Path, default=DEFAULT_DATASET,
                    help=f"Path to benchmarking dataset JSON (default: "
                         f"{DEFAULT_DATASET.name})")
    p.add_argument("--max-proteins", type=int, default=None,
                    help="Limit to first N proteins (smoke test)")
    p.add_argument("--output", type=Path, default=None,
                    help="Output directory (default: results/<timestamp>)")
    p.add_argument("--render-report", action="store_true",
                    help="Also render the markdown + PNG report after the run")
    p.add_argument("--workers", type=int, default=4,
                    help="Process N proteins concurrently (default: 4). Set to 1 "
                         "for fully sequential. Each worker still calls NCBI "
                         "serially; we just parallelise across proteins.")
    args = p.parse_args()

    for v in args.versions:
        if v not in VERSIONS:
            raise SystemExit(f"unknown version: {v}; known: {list(VERSIONS)}")
    if not args.dataset.exists():
        raise SystemExit(f"dataset not found: {args.dataset}")

    ts = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    out_dir = args.output or (RESULTS_ROOT / ts)
    out_dir.mkdir(parents=True, exist_ok=True)
    # Record what was run
    (out_dir / "run_info.json").write_text(json.dumps({
        "timestamp": ts,
        "dataset_path": str(args.dataset),
        "dataset_size": len(json.loads(args.dataset.read_text())),
        "versions": args.versions,
        "max_proteins": args.max_proteins,
    }, indent=2))

    run(args.dataset, args.versions, args.max_proteins, out_dir,
         workers=args.workers)

    if args.render_report:
        from algorithms.benchmark import report as report_mod
        report_mod.render(out_dir)


if __name__ == "__main__":
    main()
