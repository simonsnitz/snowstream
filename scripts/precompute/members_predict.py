"""Stage: per-member operon+promoter precompute for an entire InterPro family.

Reframes the precompute around individual UniProt members rather than cluster
representatives. The output (`members_predictions.jsonl`) has one record per
processed member with its operon + promoter (or nulls if NCBI couldn't resolve
them). The downstream `members_dmnd` stage filters to members with a non-null
promoter and builds a Diamond DB from their sequences, which becomes the new
smart-lookup index.

Resume: re-running skips members already present in the JSONL (by uniprot_id).

Reuse: existing per-cluster `predictions.jsonl` already contains operon+promoter
for ~380k members (as cluster homologs). We "lift" those into the member-keyed
JSONL first to avoid re-doing work. Lifted records are tagged `source: "lifted"`;
freshly computed ones are tagged `source: "computed"`.
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
from typing import Callable, Iterable, Optional

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from backend.pipeline import run_pipeline_from_homologs  # noqa: E402
from backend.schemas import PromoterParams  # noqa: E402

log = logging.getLogger(__name__)

# Per-batch wall-clock cap: a stuck NCBI call shouldn't be allowed to hold up
# a worker indefinitely. Same mechanism as scripts/precompute/predict.py's
# SNOWPRINT_PER_REP_TIMEOUT — runs the batch in a daemon thread and abandons
# it on timeout; in-flight HTTP calls eventually exit via the 30s request
# timeout in src/http_utils.py.
PER_BATCH_TIMEOUT_SEC = float(os.environ.get("SNOWPRINT_PER_BATCH_TIMEOUT", "600"))


# Per-member record carries enough to (a) rebuild a homolog dict for the
# operator finder, and (b) drive the Diamond DB index. Identity/coverage
# are intentionally omitted — they'll be filled in at query-time relative
# to the user's actual query, not relative to some cluster rep. genome /
# start / stop / strand are also omitted: the regulator's coords live
# inside operon.operon[protein_index] and the genome accession is on
# operon.genome — top-level copies would just duplicate that data.
def _empty_record(uid: str, source: str) -> dict:
    return {
        "uniprot_id": uid,
        "operon": None,
        "promoter": None,
        "source": source,
        "computed_at": datetime.now(timezone.utc).isoformat(),
    }


def _record_from_homolog(h: dict, source: str, computed_at: str | None) -> dict:
    return {
        "uniprot_id": h["uniprot_id"],
        "operon": h.get("operon"),
        "promoter": h.get("promoter"),
        "source": source,
        "computed_at": computed_at or datetime.now(timezone.utc).isoformat(),
    }


def _existing_ids(jsonl_path: Path) -> set[str]:
    """Members already in the JSONL — but timed-out batches are excluded so a
    subsequent run will re-try them (NCBI may have been transiently slow)."""
    if not jsonl_path.exists():
        return set()
    out: set[str] = set()
    with jsonl_path.open() as fh:
        for line in fh:
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if rec.get("source") == "computed_timeout":
                continue
            uid = rec.get("uniprot_id")
            if uid:
                out.add(uid)
    return out


def lift_from_clusters(
    cluster_predictions: Path,
    output_path: Path,
) -> int:
    """One pass over the cluster-based predictions.jsonl; for each homolog
    encountered, emit a per-member record IF that member isn't already in
    `output_path`. Returns the number of newly-lifted records."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    already = _existing_ids(output_path)
    appended = 0
    seen_in_this_pass: set[str] = set()
    with cluster_predictions.open() as fin, output_path.open("a") as fout:
        for line in fin:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if "error" in r:
                continue
            for h in r.get("homologs", []):
                uid = h.get("uniprot_id")
                if not uid or uid in already or uid in seen_in_this_pass:
                    continue
                seen_in_this_pass.add(uid)
                rec = _record_from_homolog(h, source="lifted", computed_at=r.get("computed_at"))
                fout.write(json.dumps(rec) + "\n")
                appended += 1
    log.info("lifted %d new per-member records to %s", appended, output_path)
    return appended


def _process_batch(batch_ids: list[str], promoter_params: dict) -> list[dict]:
    """Run one batch through the post-BLAST pipeline. Returns one record per
    input ID — populated if NCBI resolved the member, empty otherwise.

    The batch runs in a daemon thread bounded by PER_BATCH_TIMEOUT_SEC; if it
    doesn't return in time the records are marked computed_timeout and we
    abandon the in-flight thread (it will exit on its own via the 30s
    HTTP-level timeouts)."""
    homologs = [
        {"Uniprot Id": uid, "identity": None, "coverage": None}
        for uid in batch_ids
    ]
    box: dict = {}

    def _runner() -> None:
        try:
            box["result"] = run_pipeline_from_homologs(
                homologs, get_coordinates_method="batch", promoter_params=promoter_params
            )
        except BaseException as e:
            box["exc"] = e

    t = threading.Thread(target=_runner, daemon=True)
    t.start()
    t.join(timeout=PER_BATCH_TIMEOUT_SEC)

    if t.is_alive():
        log.warning("batch of %d exceeded %.0fs cap; abandoning", len(batch_ids), PER_BATCH_TIMEOUT_SEC)
        return [_empty_record(uid, source="computed_timeout") for uid in batch_ids]
    if "exc" in box:
        log.exception("batch crashed: %s", box["exc"])
        return [_empty_record(uid, source="computed_failed") for uid in batch_ids]

    result = box["result"]
    if "error" in result:
        return [_empty_record(uid, source="computed_failed") for uid in batch_ids]

    by_id = {h["uniprot_id"]: h for h in result["homologs"]}
    out: list[dict] = []
    for uid in batch_ids:
        h = by_id.get(uid)
        if h is None:
            out.append(_empty_record(uid, source="computed_missing"))
        else:
            out.append(_record_from_homolog(h, source="computed", computed_at=None))
    return out


def predict_members(
    member_ids: Iterable[str],
    output_path: Path,
    batch_size: int = 50,
    workers: int = 1,
    max_entries: Optional[int] = None,
    on_progress: Optional[Callable[[int, int, int], None]] = None,
) -> int:
    """For every member in `member_ids` not already in `output_path`, run the
    post-BLAST pipeline in batches and append per-member records. Returns
    the number of records appended."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    already = _existing_ids(output_path)
    todo = [uid for uid in member_ids if uid not in already]
    if max_entries is not None:
        todo = todo[:max_entries]
    log.info("members to process: %d (workers=%d, batch=%d)", len(todo), workers, batch_size)

    if not todo:
        return 0

    batches = [todo[i : i + batch_size] for i in range(0, len(todo), batch_size)]
    promoter_params = PromoterParams().model_dump()
    write_lock = threading.Lock()
    appended = 0
    completed_batches = 0

    def emit(records: list[dict]) -> None:
        nonlocal appended
        with write_lock:
            for rec in records:
                fh.write(json.dumps(rec) + "\n")
            fh.flush()
            appended += len(records)

    with output_path.open("a") as fh:
        if workers <= 1:
            for i, batch in enumerate(batches):
                t0 = time.time()
                records = _process_batch(batch, promoter_params)
                emit(records)
                if on_progress:
                    on_progress(i, len(batches), int(time.time() - t0))
        else:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                future_to_batch = {pool.submit(_process_batch, b, promoter_params): b for b in batches}
                for fut in as_completed(future_to_batch):
                    try:
                        records = fut.result()
                    except Exception as exc:
                        batch = future_to_batch[fut]
                        log.exception("worker raised for batch of %d: %s", len(batch), exc)
                        records = [_empty_record(uid, source="computed_failed") for uid in batch]
                    emit(records)
                    completed_batches += 1
                    if on_progress:
                        on_progress(completed_batches - 1, len(batches), 0)

    log.info("appended %d per-member records to %s", appended, output_path)
    return appended


def iter_member_ids(members_fasta: Path) -> Iterable[str]:
    """Yield UniProt accessions from a FASTA file's headers."""
    with members_fasta.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            body = line[1:].strip()
            if "|" in body:
                parts = body.split("|")
                if len(parts) >= 2 and parts[0] in ("sp", "tr"):
                    yield parts[1]
                    continue
            yield body.split()[0]
