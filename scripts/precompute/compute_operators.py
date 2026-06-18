"""Stage: compute predicted operators per cluster centroid using V2.6.

For each cluster centroid that has a populated promoter, runs the V2.6
operator_fetch (widened candidate pool + AT-based selector at gc_medium_weak
strength) using the cluster's members as the homolog set. Appends a new
JSONL record with all original fields preserved plus the operator data,
following the operon_retry / promoter_retry pattern. The byte-offset
index in backend/families.py automatically prefers the new record.

Why per-centroid (not per-member):
  Smart-lookup returns the matched protein's record at query time. For
  TetR the family DB has ~36k cluster centroids with promoters; precomputing
  operators on those covers every cluster. A future stage can propagate
  the centroid's operator out to its non-centroid cluster members so that
  smart-lookup hits on non-centroids also return a precomputed operator.

What gets stored, per centroid:
  * operator_fetch_version       — "v2.6"
  * operator_motif               — consensus motif sequence string (e.g.
                                   "ataatAAACGGAGAGTTATCCGTTTgtcaa")
  * operator_native_for_centroid — the operator extracted from the centroid's
                                   own promoter (lowercase flanks + uppercase
                                   core, same as v0 output)
  * operator_consensus_score     — V0's selection metric, kept for reference
  * operator_rerank_score        — V2.6's selection metric (gc_medium_weak)
  * operator_rerank_scorer       — "gc_medium_weak"
  * operator_n_homologs_used     — count of cluster members that contributed
                                   a valid alignment (subset of the
                                   homolog list — some align below cutoff)
  * operator_n_candidates_eval'd — count of palindrome candidates that
                                   produced a non-empty consensus
  * operator_computed_at         — ISO timestamp
  * source                       — "compute_operators_v2.6"

Skip behaviour: centroids whose latest JSONL record already has
operator_fetch_version=="v2.6" are skipped, so re-running is idempotent.
"""

from __future__ import annotations

import json
import logging
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from algorithms.operator_fetch.v1 import fetch as operator_v1_fetch  # noqa: E402

log = logging.getLogger(__name__)

# V2.6 = widened candidate pool + gc_medium_weak selector.
# Other operator params match DEFAULT_OPERATOR_PARAMS in
# algorithms/benchmark/run.py.
_V26_PARAMS = {
    "candidate_strategy": "widened",
    "rerank_scorer": "gc_medium_weak",
    "extension_length": 5,
    "win_score": 2,
    "loss_score": -2,
    "spacer_penalty": {str(i): v for i, v in zip(range(21),
        [4, 4, 4, 4, 4, 2, 2, 0, 0, -2, -2,
         -4, -4, -6, -6, -8, -8, -10, -10, -12, -12])},
    "gap_open": -100,
    "gap_extend": 0,
    "align_match": 2,
    "align_mismatch": -0.5,
    "min_operator_length": 5,
    "max_operator_length": 15,
    "seq_to_align": None,
    "search_method": "Look for inverted repeats",
}

_VERSION = "v2.6"
_SOURCE = "compute_operators_v2.6"

# Cap homologs per centroid so very-large clusters don't blow up runtime.
# Smart-lookup's default at query time is 100, which is a reasonable upper
# bound here too.
_MAX_HOMOLOGS_PER_CLUSTER = 100


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


def _load_latest_records(jsonl_path: Path) -> dict[str, dict]:
    """Latest-record-per-uid scan. Skip JSON-broken lines."""
    out: dict[str, dict] = {}
    if not jsonl_path.exists():
        return out
    with jsonl_path.open() as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            uid = r.get("uniprot_id")
            if uid:
                out[uid] = r
    return out


def _build_homolog_input(centroid: str, cluster_members: list[dict],
                          records: dict[str, dict]) -> list[dict]:
    """Centroid first (the V2.6 anchor); then other cluster members that
    have a populated promoter in the precomputed DB. Each dict shape
    matches what operator_fetch expects: {"Uniprot Id", "promoter"}."""
    out: list[dict] = []
    c_rec = records.get(centroid) or {}
    c_promoter = c_rec.get("promoter")
    if c_promoter:
        out.append({"Uniprot Id": centroid, "promoter": c_promoter})

    for m in cluster_members:
        uid = m.get("uniprot_id")
        if not uid or uid == centroid:
            continue
        rec = records.get(uid)
        if not rec:
            continue
        prom = rec.get("promoter")
        if not prom:
            continue
        out.append({"Uniprot Id": uid, "promoter": prom})
        if len(out) >= _MAX_HOMOLOGS_PER_CLUSTER:
            break
    return out


def _operator_record(centroid_rec: dict, result: dict) -> dict:
    """Build the augmented JSONL record: all original fields + operator
    fields + new source/timestamp. `centroid_rec` is the latest existing
    record for this centroid in the JSONL (so operon/promoter/etc. carry
    through unchanged)."""
    motif_str = _render_motif(result.get("motif"))
    rerank_score = result.get("rerank_score")
    if rerank_score in (None, float("-inf")):
        rerank_score = None

    out = {
        # Preserve everything from the original record (operon, promoter,
        # protein_index, genome, etc.). Override only `source` and
        # `computed_at`, and add the new operator_* fields.
        **{k: v for k, v in centroid_rec.items()
            if k not in ("source", "computed_at")},
        "source": _SOURCE,
        "computed_at": datetime.now(timezone.utc).isoformat(),
        "operator_fetch_version": _VERSION,
        "operator_motif": motif_str,
        "operator_native_for_centroid": result.get("native_operator"),
        "operator_consensus_score": result.get("consensus_score"),
        "operator_rerank_score": rerank_score,
        "operator_rerank_scorer": result.get("rerank_scorer"),
        "operator_n_homologs_used": result.get("num_seqs"),
        "operator_n_candidates_evaluated": result.get("n_candidates_evaluated"),
        "operator_computed_at": datetime.now(timezone.utc).isoformat(),
    }
    return out


def find_candidates(
    jsonl_path: Path,
    cluster_homologs_path: Path,
) -> tuple[list[str], dict[str, dict], dict[str, list[dict]]]:
    """Return:
      * list of cluster centroid uniprot_ids that need operator computation
        (centroid has a promoter, and its latest record doesn't already
        carry operator_fetch_version == v2.6)
      * mapping uniprot_id → latest_record  (for building homolog inputs)
      * mapping centroid_uid → cluster_member_list
    """
    if not cluster_homologs_path.exists():
        log.warning("cluster_homologs.json missing: %s", cluster_homologs_path)
        return [], {}, {}
    cluster_homologs: dict[str, list[dict]] = json.loads(
        cluster_homologs_path.read_text())
    records = _load_latest_records(jsonl_path)

    candidates: list[str] = []
    for centroid in cluster_homologs.keys():
        rec = records.get(centroid)
        if not rec or not rec.get("promoter"):
            continue
        if rec.get("operator_fetch_version") == _VERSION:
            continue  # already done — resume support
        candidates.append(centroid)
    return candidates, records, cluster_homologs


def compute_operators(
    jsonl_path: Path,
    cluster_homologs_path: Path,
    workers: int = 4,
    max_records: Optional[int] = None,
    on_progress: Optional[Callable[[int, int, int], None]] = None,
) -> int:
    """For each pending centroid, run V2.6 and append a new JSONL record.
    Returns the count of records appended. Idempotent on re-run."""
    log.info("scanning %s for centroids needing operator computation",
              jsonl_path)
    candidates, records, cluster_homologs = find_candidates(
        jsonl_path, cluster_homologs_path)
    log.info("  %d cluster centroids need V2.6 operators", len(candidates))
    if max_records is not None:
        candidates = candidates[:max_records]
        log.info("  limited to %d for this run", len(candidates))
    if not candidates:
        return 0

    write_lock = threading.Lock()
    appended = 0
    successes = 0
    completed = 0
    total = len(candidates)

    def one(centroid: str):
        try:
            homologs = _build_homolog_input(
                centroid, cluster_homologs.get(centroid, []) or [], records)
        except Exception as exc:
            return centroid, None, f"build_homologs: {exc}"
        if not homologs:
            return centroid, None, "no_homologs_with_promoter"
        try:
            result = operator_v1_fetch(homologs, _V26_PARAMS)
        except Exception as exc:
            return centroid, None, f"operator_fetch: {exc}"
        return centroid, result, None

    with jsonl_path.open("a") as fh, ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(one, c): c for c in candidates}
        for fut in as_completed(futures):
            try:
                centroid, result, err = fut.result()
            except Exception as exc:
                log.warning("worker raised: %s", exc)
                with write_lock:
                    completed += 1
                continue

            with write_lock:
                completed += 1
                if result and not err:
                    centroid_rec = records.get(centroid) or {}
                    rec = _operator_record(centroid_rec, result)
                    fh.write(json.dumps(rec) + "\n")
                    fh.flush()
                    appended += 1
                    if result.get("num_seqs", 0):
                        successes += 1
                if on_progress and (completed % 500 == 0 or completed == total):
                    on_progress(completed, total, successes)

    log.info("appended %d / %d centroids (%d with non-empty consensus)",
             appended, total, successes)
    return appended
