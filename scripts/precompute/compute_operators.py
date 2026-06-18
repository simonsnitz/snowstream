"""Stage: compute predicted operators per cluster centroid.

For each cluster centroid with a populated promoter, runs an operator-finding
algorithm using the cluster's members as the homolog set, and appends a new
JSONL record with all original fields preserved plus a version-keyed nested
block (e.g. `operator_v26`, `operator_v0`). The algorithm version is in the
*key name itself* — so the same centroid can carry several blocks side-by-
side (`operator_v0` + `operator_v26` + ... future ones) without renaming
anything, and consumers can tell at a glance which algorithm produced
which fields.

Currently registered versions:

  * `v0`   — legacy `src/fetch_operator.fetch_operator`, narrow candidate
             pool (max-IR-tied palindromes), winner picked by V0's
             consensus_score.
  * `v2.6` — `algorithms/operator_fetch/v1.fetch` with
             `candidate_strategy="widened"` and
             `rerank_scorer="gc_medium_weak"` — the PR #16 winner.

What gets stored per centroid (using v2.6 as an example):

    "operator_v26": {
        "version":              "v2.6",
        "algorithm":            "operator_fetch.v1 (widened + gc_medium_weak)",
        "candidate_strategy":   "widened",
        "rerank_scorer":        "gc_medium_weak",
        "motif":                "...consensus motif string...",
        "native_for_centroid":  "...operator in centroid's own promoter...",
        "consensus_score":      <float>,
        "rerank_score":         <float>,   # v2.6 only
        "n_homologs_used":      <int>,
        "n_candidates_evaluated": <int>,   # v2.6 only
        "computed_at":          "<ISO timestamp>"
    }

V0's block has the same shape minus `candidate_strategy`, `rerank_scorer`,
`rerank_score`, and `n_candidates_evaluated` (V0 doesn't track those).

Top-level `source` records which version produced the record, e.g.
`"compute_operators_v2.6"` or `"compute_operators_v0"`.

Skip behaviour: centroids whose latest JSONL record already has a non-null
block for the requested version are skipped, so re-running each version
is idempotent independently. Re-running v2.6 doesn't touch v0 blocks
that may already be present on the same record, and vice versa.
"""

from __future__ import annotations

import json
import logging
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Optional

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from algorithms.operator_fetch.v0 import fetch as operator_v0_fetch  # noqa: E402
from algorithms.operator_fetch.v1 import fetch as operator_v1_fetch  # noqa: E402

log = logging.getLogger(__name__)

# Shared operator-finding params. Both V0 and V2.6 use the same numeric
# settings; V2.6 layers on its own `candidate_strategy` / `rerank_scorer`
# overrides in V26_PARAMS.
_BASE_PARAMS = {
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

V0_PARAMS = dict(_BASE_PARAMS)
V26_PARAMS = {
    **_BASE_PARAMS,
    "candidate_strategy": "widened",
    "rerank_scorer": "gc_medium_weak",
}

# Cap homologs per centroid so very-large clusters don't blow up runtime.
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


# --- Per-version spec ---------------------------------------------------

@dataclass(frozen=True)
class VersionSpec:
    """Everything that varies between operator-finding versions."""
    version: str                       # "v0", "v2.6", ...
    algorithm_label: str               # human-readable description
    operator_key: str                  # JSONL block key, e.g. "operator_v26"
    source: str                        # top-level `source` value
    fetch: Callable[[list[dict], dict], dict]
    fetch_params: dict
    block_builder: Callable[[dict], dict]
    # Legacy operator-* flat-shape keys produced by an earlier cut of this
    # stage (V2.6 only). Stripped from new records to keep things clean.
    legacy_flat_keys: tuple[str, ...] = field(default_factory=tuple)


def _v0_block(result: dict) -> dict:
    """Build the `operator_v0` payload from a legacy fetch_operator result."""
    return {
        "version": V0_SPEC.version,
        "algorithm": V0_SPEC.algorithm_label,
        "motif": _render_motif(result.get("motif")),
        "native_for_centroid": result.get("native_operator"),
        "consensus_score": result.get("consensus_score"),
        "n_homologs_used": result.get("num_seqs"),
        "computed_at": datetime.now(timezone.utc).isoformat(),
    }


def _v26_block(result: dict) -> dict:
    """Build the `operator_v26` payload from a V2.6 fetch result."""
    rerank_score = result.get("rerank_score")
    if rerank_score in (None, float("-inf")):
        rerank_score = None
    return {
        "version": V26_SPEC.version,
        "algorithm": V26_SPEC.algorithm_label,
        "candidate_strategy": V26_PARAMS["candidate_strategy"],
        "rerank_scorer": V26_PARAMS["rerank_scorer"],
        "motif": _render_motif(result.get("motif")),
        "native_for_centroid": result.get("native_operator"),
        "consensus_score": result.get("consensus_score"),
        "rerank_score": rerank_score,
        "n_homologs_used": result.get("num_seqs"),
        "n_candidates_evaluated": result.get("n_candidates_evaluated"),
        "computed_at": datetime.now(timezone.utc).isoformat(),
    }


V0_SPEC = VersionSpec(
    version="v0",
    algorithm_label="operator_fetch.v0 (legacy)",
    operator_key="operator_v0",
    source="compute_operators_v0",
    fetch=operator_v0_fetch,
    fetch_params=V0_PARAMS,
    block_builder=_v0_block,
)

V26_SPEC = VersionSpec(
    version="v2.6",
    algorithm_label="operator_fetch.v1 (widened + gc_medium_weak)",
    operator_key="operator_v26",
    source="compute_operators_v2.6",
    fetch=operator_v1_fetch,
    fetch_params=V26_PARAMS,
    block_builder=_v26_block,
    legacy_flat_keys=(
        "operator_fetch_version",
        "operator_motif",
        "operator_native_for_centroid",
        "operator_consensus_score",
        "operator_rerank_score",
        "operator_rerank_scorer",
        "operator_n_homologs_used",
        "operator_n_candidates_evaluated",
        "operator_computed_at",
    ),
)

SPECS: dict[str, VersionSpec] = {V0_SPEC.version: V0_SPEC,
                                  V26_SPEC.version: V26_SPEC}


# --- I/O helpers --------------------------------------------------------

def _has_block(rec: dict, spec: VersionSpec) -> bool:
    """True iff the record already carries operator data for this version,
    either in the new nested shape or — for v2.6 only — the legacy flat
    shape from the very first cut of this stage. Both count as 'already
    done' for resume / idempotency."""
    block = rec.get(spec.operator_key)
    if isinstance(block, dict) and block.get("version") == spec.version:
        return True
    if (spec.version == V26_SPEC.version
            and rec.get("operator_fetch_version") == spec.version):
        return True
    return False


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
    """Centroid first (the algorithm's anchor); then other cluster members
    that have a populated promoter in the precomputed DB."""
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


def _operator_record(centroid_rec: dict, result: dict,
                      spec: VersionSpec) -> dict:
    """Build the augmented JSONL record: all original fields preserved
    (operon, promoter, protein_index, genome, and any pre-existing
    `operator_v*` blocks from other versions) plus this version's nested
    block + new source/computed_at."""
    strip = {"source", "computed_at", spec.operator_key, *spec.legacy_flat_keys}
    return {
        **{k: v for k, v in centroid_rec.items() if k not in strip},
        "source": spec.source,
        "computed_at": datetime.now(timezone.utc).isoformat(),
        spec.operator_key: spec.block_builder(result),
    }


def find_candidates(
    jsonl_path: Path,
    cluster_homologs_path: Path,
    spec: VersionSpec,
) -> tuple[list[str], dict[str, dict], dict[str, list[dict]]]:
    """Return:
      * list of centroid uniprot_ids that need operator computation for
        this version (centroid has a promoter; its latest record doesn't
        already carry a block for this version)
      * latest_record by uniprot_id (for building homolog inputs)
      * cluster_member list by centroid uid
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
        if _has_block(rec, spec):
            continue
        candidates.append(centroid)
    return candidates, records, cluster_homologs


# --- Driver -------------------------------------------------------------

def compute_operators(
    jsonl_path: Path,
    cluster_homologs_path: Path,
    version: str = V26_SPEC.version,
    workers: int = 4,
    max_records: Optional[int] = None,
    on_progress: Optional[Callable[[int, int, int], None]] = None,
) -> int:
    """For each pending centroid, run the requested version's operator
    finder and append a new JSONL record with the version's nested block.
    Returns the count of records appended. Idempotent on re-run."""
    if version not in SPECS:
        raise KeyError(f"unknown operator version {version!r}; "
                        f"choose from {list(SPECS)}")
    spec = SPECS[version]

    log.info("scanning %s for centroids needing %s operators",
              jsonl_path, spec.version)
    candidates, records, cluster_homologs = find_candidates(
        jsonl_path, cluster_homologs_path, spec)
    log.info("  %d cluster centroids need %s operators",
              len(candidates), spec.version)
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
            result = spec.fetch(homologs, spec.fetch_params)
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
                    rec = _operator_record(centroid_rec, result, spec)
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


# --- Migration ----------------------------------------------------------

def migrate_legacy_records(jsonl_path: Path) -> int:
    """One-shot migration for records written by the very first cut of
    this stage (when fields were flat `operator_fetch_version` /
    `operator_motif` / etc. rather than the nested `operator_v26` block).

    For each centroid whose latest record has `operator_fetch_version ==
    "v2.6"` but no `operator_v26` block, append a new record that
    repackages the existing operator data into the nested shape. The
    byte-offset index then prefers the latest (nested) record. Returns
    the count of records appended.
    """
    spec = V26_SPEC
    records = _load_latest_records(jsonl_path)
    to_migrate: list[tuple[str, dict]] = []
    for uid, rec in records.items():
        if isinstance(rec.get(spec.operator_key), dict):
            continue
        if rec.get("operator_fetch_version") != spec.version:
            continue
        to_migrate.append((uid, rec))
    log.info("migration: %d records to migrate to nested %s",
             len(to_migrate), spec.operator_key)
    if not to_migrate:
        return 0
    appended = 0
    with jsonl_path.open("a") as fh:
        for uid, rec in to_migrate:
            block = {
                "version": spec.version,
                "algorithm": spec.algorithm_label,
                "candidate_strategy": spec.fetch_params["candidate_strategy"],
                "rerank_scorer": (rec.get("operator_rerank_scorer")
                                   or spec.fetch_params["rerank_scorer"]),
                "motif": rec.get("operator_motif"),
                "native_for_centroid": rec.get("operator_native_for_centroid"),
                "consensus_score": rec.get("operator_consensus_score"),
                "rerank_score": rec.get("operator_rerank_score"),
                "n_homologs_used": rec.get("operator_n_homologs_used"),
                "n_candidates_evaluated":
                    rec.get("operator_n_candidates_evaluated"),
                "computed_at": rec.get("operator_computed_at"),
            }
            strip = {"source", "computed_at", *spec.legacy_flat_keys}
            new_rec = {
                **{k: v for k, v in rec.items() if k not in strip},
                "source": spec.source,
                "computed_at": datetime.now(timezone.utc).isoformat(),
                spec.operator_key: block,
            }
            fh.write(json.dumps(new_rec) + "\n")
            appended += 1
    log.info("migration: appended %d nested records", appended)
    return appended
