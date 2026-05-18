"""Stage: retry fetch_promoter for records whose first attempt failed
in the "should have worked" bucket — where the algorithm picked a valid
intergenic gap (80-799 bp) but NCBI's eFetch returned nothing usable.

About half of the 200k operon-but-no-promoter records fall here. The chosen
positions are biologically sensible; the failure was almost certainly a
transient NCBI issue (the same "silently drops responses under concurrency"
behaviour we saw on the IPG side).

For each such record we re-run `fetch_promoter` directly. The streamlit
in-memory cache key includes the homolog dict so distinct records don't
collide, and we retry the NCBI call with exponential backoff inside a
small wrapper. Successfully-recovered records are appended to the JSONL
with `source: "promoter_retry"` so the byte-offset index picks them up
as the latest write.

Pre-filtering by chosen gap size avoids spending NCBI calls on records
that would fail the length filter regardless (gap < 80 or > 800 bp).
"""

from __future__ import annotations

import json
import logging
import os
import random
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.fetch_promoter import fetch_promoter  # noqa: E402
from src.http_utils import ncbi_get  # noqa: E402
from backend.schemas import PromoterParams  # noqa: E402

log = logging.getLogger(__name__)

PER_RETRY_TIMEOUT = float(os.environ.get("SNOWPRINT_PROMOTER_RETRY_TIMEOUT", "60"))


def _chosen_gap(operon: list, protein_index: int, min_length: int = 80) -> Optional[int]:
    """Simulate fetch_promoter's choice (the live, un-fixed logic) and
    return the size of the intergenic gap it would request from NCBI,
    or None if the algorithm wouldn't pick a position at all."""
    if not operon or protein_index is None:
        return None
    n = len(operon)
    if protein_index < 0 or protein_index >= n:
        return None
    reg = operon[protein_index]
    direction = reg.get("direction")
    if direction == "+":
        if protein_index == 0:
            return None
        upstream = list(reversed(operon[0:protein_index]))
        try:
            immediate_gap = int(operon[protein_index]["start"]) - int(operon[protein_index - 1]["stop"])
        except (KeyError, ValueError, TypeError):
            return None
        index = protein_index
        for i in upstream:
            if i.get("direction") == "-":
                try:
                    return int(operon[index]["start"]) - int(i["stop"])
                except (KeyError, ValueError, TypeError):
                    return None
            if immediate_gap > min_length:
                return immediate_gap
            if index == 1:
                return None
            index -= 1
    elif direction == "-":
        if protein_index >= n - 1:
            return None
        downstream = operon[protein_index + 1:]
        try:
            immediate_gap = int(operon[protein_index + 1]["start"]) - int(operon[protein_index]["stop"])
        except (KeyError, ValueError, TypeError):
            return None
        index = protein_index
        for i in downstream:
            if i.get("direction") == "+":
                try:
                    return int(i["start"]) - int(operon[index]["stop"])
                except (KeyError, ValueError, TypeError):
                    return None
            if immediate_gap > 100:  # original hardcoded threshold for '-' regulators
                return immediate_gap
            if index == n - 2:
                return None
            index += 1
    return None


def _new_record(orig: dict, promoter: str) -> dict:
    """Build the per-member JSONL record with the recovered promoter,
    preserving the operon / protein_index / genome from the original."""
    return {
        "uniprot_id": orig["uniprot_id"],
        "genome": orig.get("genome"),
        "promoter": promoter,
        "source": "promoter_retry",
        "computed_at": datetime.now(timezone.utc).isoformat(),
        "protein_index": orig.get("protein_index"),
        "operon": orig.get("operon"),
    }


def _try_promoter(record: dict, params: dict, max_attempts: int = 3) -> Optional[str]:
    """Call fetch_promoter with simple retry/backoff. Returns the
    promoter string on success, None on permanent failure."""
    homolog = {
        "operon": record["operon"],
        "protein_index": record["protein_index"],
        "genome": record["genome"],
    }
    last_exc: Optional[BaseException] = None
    for attempt in range(max_attempts):
        try:
            prom = fetch_promoter(homolog, params)
            if prom:
                return prom
        except BaseException as exc:
            last_exc = exc
        if attempt < max_attempts - 1:
            time.sleep((1.5 ** attempt) + random.uniform(0, 0.5))
    if last_exc is not None:
        log.debug("retry failed for %s: %s", record.get("uniprot_id"), last_exc)
    return None


def find_candidates(
    jsonl_path: Path,
    min_gap: int = 80,
    max_gap: int = 800,
) -> list[dict]:
    """Walk the JSONL (latest-record-per-uid), keep records where the
    chosen intergenic gap falls in [min_gap, max_gap] but the promoter
    is null. These are the "should-have-worked" candidates."""
    latest: dict[str, dict] = {}
    if not jsonl_path.exists():
        return []
    with jsonl_path.open() as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            uid = r.get("uniprot_id")
            if uid:
                latest[uid] = r
    out: list[dict] = []
    for r in latest.values():
        if r.get("promoter") or not r.get("operon"):
            continue
        if not r.get("genome") or r.get("protein_index") is None:
            continue
        gap = _chosen_gap(r["operon"], r["protein_index"], min_length=min_gap)
        if gap is None:
            continue
        if min_gap <= gap <= max_gap:
            out.append(r)
    return out


def retry_promoters(
    jsonl_path: Path,
    workers: int = 4,
    max_records: Optional[int] = None,
    on_progress: Optional[Callable[[int, int, int], None]] = None,
) -> int:
    """Run the retry pass. Returns the number of records appended."""
    log.info("scanning %s for retry candidates", jsonl_path)
    candidates = find_candidates(jsonl_path)
    log.info("  %d candidates (gap in [80, 800])", len(candidates))
    if max_records is not None:
        candidates = candidates[:max_records]
        log.info("  limited to %d for this run", len(candidates))

    if not candidates:
        return 0

    params = PromoterParams().model_dump()
    write_lock = threading.Lock()
    appended = 0
    recovered = 0
    completed = 0

    def one(r):
        prom = _try_promoter(r, params)
        return r, prom

    with jsonl_path.open("a") as fh, ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(one, r): r for r in candidates}
        for fut in as_completed(futures):
            try:
                r, prom = fut.result()
            except Exception as exc:
                log.warning("worker raised: %s", exc)
                continue
            with write_lock:
                completed += 1
                if prom:
                    fh.write(json.dumps(_new_record(r, prom)) + "\n")
                    fh.flush()
                    appended += 1
                    recovered += 1
                if on_progress and (completed % 500 == 0 or completed == len(candidates)):
                    on_progress(completed, len(candidates), recovered)

    log.info("recovered %d / %d candidates → %d records appended", recovered, len(candidates), appended)
    return appended
