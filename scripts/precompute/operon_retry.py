"""Stage: retry operon resolution for records that ended up with no operon.

Symmetric counterpart to `promoter_retry`. Most of the ~250k no-operon records
in members_predictions.jsonl failed because of the same NCBI-silently-drops-
responses-under-concurrency behaviour we saw on the promoter side: a hand-
sampled batch from the "had IPG coords but no operon" bucket recovered 20/20
on a quiet retry.

For each empty record we:
  1. Find a RefSeq proxy.
       * `source == "tsv_proxy"` records already carry `proxy_accession` — use it.
       * Other sources (`computed`, `lifted`) need a lookup via the inverted
         RefSeq → UniProt map produced by the `xref_index` stage.
  2. Hit NCBI IPG for the proxy → CDS coords (with retry).
  3. Run `acc2operon` for those coords (with retry — acc2operon already
     tries 4 flanking-window sizes internally, so retrying the whole call
     handles the "all four windows came back empty under load" case).
  4. If the operon resolves, run `fetch_promoter` while we're here.
  5. Append a new record (`source == "operon_retry"`) so the byte-offset
     index picks up the populated copy.

The retry uses light per-record parallelism (workers default 4) rather than
batching, so transient drops don't poison a whole batch of accessions.

Records whose UniProt accession isn't a value in the RefSeq → UniProt map
(i.e. UniProt had no RefSeq xref for them at TSV-export time) are skipped
here — they need an EMBL-fallback proxy source, or a fresh live UniProt
lookup, which is left to a separate stage.
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

import xmltodict

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.accID2operon import acc2operon  # noqa: E402
from src.fetch_promoter import fetch_promoter  # noqa: E402
from src.http_utils import ncbi_get  # noqa: E402
from backend.schemas import PromoterParams  # noqa: E402

log = logging.getLogger(__name__)

PER_CALL_TIMEOUT = float(os.environ.get("SNOWPRINT_OPERON_RETRY_TIMEOUT", "60"))

IPG_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=protein&id={acc}&rettype=ipg"
)


def _invert_refseq_map(refseq_to_uniprot: dict[str, str]) -> dict[str, str]:
    """Build uniprot → one-of-its-refseqs from the existing refseq → uniprot map.

    A UniProt accession typically has several RefSeqs; we only need one to feed
    into IPG. The xref_index stage already collapsed RefSeq → UniProt as 1:1
    within a family, so iterating and taking the first hit per UniProt is fine.
    """
    out: dict[str, str] = {}
    for ref, uid in refseq_to_uniprot.items():
        out.setdefault(uid, ref)
    return out


def _ipg_coords(acc: str, max_attempts: int = 3) -> Optional[dict]:
    """Fetch IPG for a single accession with retry. Returns
    `{"Genome", "Start", "Stop", "Strand"}` on success, else None."""
    last_exc: Optional[BaseException] = None
    for attempt in range(max_attempts):
        try:
            resp = ncbi_get(IPG_URL.format(acc=acc), timeout=PER_CALL_TIMEOUT)
            if not resp.ok:
                raise RuntimeError(f"http {resp.status_code}")
            parsed = xmltodict.parse(resp.text)
            report = parsed.get("IPGReportSet", {}).get("IPGReport", {})
            proteins = report.get("ProteinList", {}).get("Protein") if report else None
            if not proteins:
                return None
            if isinstance(proteins, list):
                proteins = proteins[0]
            cds = proteins.get("CDSList", {}).get("CDS")
            if isinstance(cds, list):
                cds = cds[0]
            if not cds:
                return None
            return {
                "Genome": cds["@accver"],
                "Start": cds["@start"],
                "Stop": cds["@stop"],
                "Strand": cds["@strand"],
            }
        except BaseException as exc:
            last_exc = exc
            if attempt < max_attempts - 1:
                time.sleep((1.5 ** attempt) + random.uniform(0, 0.5))
    log.debug("ipg failed for %s: %s", acc, last_exc)
    return None


def _try_operon(coords: dict, uniprot_id: str, max_attempts: int = 3) -> Optional[dict]:
    """Call acc2operon with retry. acc2operon already tries 4 flanking-window
    sizes internally, so per-call retry here only covers the "every window
    timed out" case. Returns the nested operon dict or None."""
    h = {"Uniprot Id": uniprot_id, **coords}
    last_exc: Optional[BaseException] = None
    for attempt in range(max_attempts):
        try:
            op = acc2operon(h)
            if op and op.get("operon"):
                return op
        except BaseException as exc:
            last_exc = exc
        if attempt < max_attempts - 1:
            time.sleep((1.5 ** attempt) + random.uniform(0, 0.5))
    if last_exc is not None:
        log.debug("acc2operon failed for %s: %s", uniprot_id, last_exc)
    return None


def _new_record(uid: str, proxy: str, op: dict, promoter: Optional[str]) -> dict:
    return {
        "uniprot_id": uid,
        "genome": op.get("genome"),
        "promoter": promoter,
        "source": "operon_retry",
        "computed_at": datetime.now(timezone.utc).isoformat(),
        "proxy_accession": proxy,
        "protein_index": op.get("protein_index"),
        "operon": op.get("operon"),
    }


def find_candidates(
    jsonl_path: Path,
    uid_to_refseq: dict[str, str],
) -> list[tuple[str, str]]:
    """Walk the JSONL (latest-record-per-uid), keep records with no operon
    that we can attach a RefSeq proxy to. Returns `[(uniprot_id, proxy_accession), ...]`.
    """
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

    out: list[tuple[str, str]] = []
    no_proxy = 0
    for uid, r in latest.items():
        if r.get("operon"):
            continue
        proxy = r.get("proxy_accession") or uid_to_refseq.get(uid)
        if not proxy:
            no_proxy += 1
            continue
        out.append((uid, proxy))
    if no_proxy:
        log.info("  %d empties have no RefSeq proxy available — skipping", no_proxy)
    return out


def retry_operons(
    jsonl_path: Path,
    refseq_map_path: Optional[Path] = None,
    workers: int = 4,
    max_records: Optional[int] = None,
    on_progress: Optional[Callable[[int, int, int], None]] = None,
) -> int:
    """Run the retry pass. Returns the number of records appended."""
    if refseq_map_path and refseq_map_path.exists():
        refseq_to_uniprot = json.loads(refseq_map_path.read_text())
        uid_to_refseq = _invert_refseq_map(refseq_to_uniprot)
        log.info("loaded %d UniProt → RefSeq entries from %s",
                 len(uid_to_refseq), refseq_map_path.name)
    else:
        uid_to_refseq = {}
        log.warning("no refseq_to_uniprot.json; only records with stored "
                    "proxy_accession will be retried")

    log.info("scanning %s for operon-retry candidates", jsonl_path)
    candidates = find_candidates(jsonl_path, uid_to_refseq)
    log.info("  %d candidates", len(candidates))
    if max_records is not None:
        candidates = candidates[:max_records]
        log.info("  limited to %d for this run", len(candidates))

    if not candidates:
        return 0

    promoter_params = PromoterParams().model_dump()
    write_lock = threading.Lock()
    appended = 0
    recovered_op = 0
    recovered_prom = 0
    completed = 0
    total = len(candidates)

    def one(item):
        uid, proxy = item
        coords = _ipg_coords(proxy)
        if not coords:
            return uid, proxy, None, None
        op = _try_operon(coords, uid)
        if not op:
            return uid, proxy, None, None
        try:
            prom = fetch_promoter(op, promoter_params)
        except BaseException:
            prom = None
        return uid, proxy, op, prom

    with jsonl_path.open("a") as fh, ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(one, c): c for c in candidates}
        for fut in as_completed(futures):
            try:
                uid, proxy, op, prom = fut.result()
            except Exception as exc:
                log.warning("worker raised: %s", exc)
                continue
            with write_lock:
                completed += 1
                if op:
                    fh.write(json.dumps(_new_record(uid, proxy, op, prom)) + "\n")
                    fh.flush()
                    appended += 1
                    recovered_op += 1
                    if prom:
                        recovered_prom += 1
                if on_progress and (completed % 500 == 0 or completed == total):
                    on_progress(completed, total, recovered_op)

    log.info(
        "recovered %d operons (+ %d promoters) from %d candidates → %d records appended",
        recovered_op, recovered_prom, total, appended,
    )
    return appended
