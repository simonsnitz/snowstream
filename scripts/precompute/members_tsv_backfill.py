"""Stage: backfill empty per-member records by proxying through UniProt's
RefSeq cross-references from an InterPro TSV export.

Many of the empty records in members_predictions.jsonl are empty because the
in-memory pipeline's `uniprot2EMBL` only inspects the first
`uniProtKBCrossReferences` entry and silently drops anything that doesn't
expose a `ProteinId` property. UniProt's downloadable TSV (Entry / RefSeq /
EMBL columns) reliably contains the cross-references we need; piping each
empty record's RefSeq accession into NCBI's IPG endpoint recovers most of
them.

For each empty record we:
  1. Look up its RefSeq (preferred) or protein-shaped EMBL accession in the TSV.
  2. Batch the accessions through NCBI IPG (with per-batch retry: NCBI silently
     drops some accessions from batched responses under concurrency, and a
     second pass typically recovers them).
  3. For each report that comes back with usable CDS coordinates, run
     acc2operon → fetch_promoter — same path as the regular per-member pipeline.
  4. Append a new per-member record (source="tsv_proxy", with the proxy
     accession tracked) to the JSONL. The byte-offset index in
     backend/families.py automatically prefers the later record, so reads
     pick up the populated copy.

Existing populated records and any record with `source == "tsv_proxy"` are
left alone — only empties are re-processed.
"""

from __future__ import annotations

import json
import logging
import os
import re
import subprocess  # noqa: F401  (kept for ergonomics; not used directly)
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

# Per-record operon+promoter calls run in their own thread pool inside each
# batch. Matches the SNOWPRINT_OPERON_WORKERS knob used by the regular
# per-member pipeline so we get the same intra-batch parallelism (otherwise
# acc2operon and fetch_promoter run serially for each of the up-to-10 records
# in a batch — observed to limit throughput to ~3-4 records/sec).
_INNER_OPERON_WORKERS = int(os.environ.get("SNOWPRINT_OPERON_WORKERS", "6"))

# An EMBL accession that points at a single protein looks like "CAA12345" or
# "AAB67890.1" — three letters, 5–7 digits, optional version. Longer
# alphanumeric blobs are WGS contigs/scaffolds and IPG can't resolve them.
_EMBL_PROTEIN_RE = re.compile(r"^[A-Z]{3}\d{5,7}(\.\d+)?$")

IPG_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=protein&id={ids}&rettype=ipg"
)


# --- I/O helpers ---------------------------------------------------------

def load_tsv(tsv_path: Path) -> dict[str, tuple[list[str], list[str]]]:
    """Parse the InterPro UniProt TSV → {uniprot_id: (refseqs, protein_embls)}.

    EMBL accessions are filtered to ones that look like protein IDs; WGS
    contigs are dropped. The Sequence column (if present) is ignored here
    — separate stage handles sequence indexing.
    """
    mapping: dict[str, tuple[list[str], list[str]]] = {}
    with tsv_path.open() as fh:
        fh.readline()  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            uid, refseq_col, embl_col = parts[0], parts[1].strip(";"), parts[2].strip(";")
            refs = [x for x in refseq_col.split(";") if x] if refseq_col else []
            embls_all = [x for x in embl_col.split(";") if x] if embl_col else []
            embls_protein = [x for x in embls_all if _EMBL_PROTEIN_RE.match(x)]
            mapping[uid] = (refs, embls_protein)
    return mapping


def find_empty_ids(jsonl_path: Path) -> list[str]:
    """Return the uniprot_ids of records with no operon and no promoter."""
    empty: list[str] = []
    if not jsonl_path.exists():
        return empty
    with jsonl_path.open() as fh:
        # Latest-wins: track the most recent record per uid, then filter.
        latest: dict[str, dict] = {}
        for line in fh:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            uid = r.get("uniprot_id")
            if uid:
                latest[uid] = r
    for uid, r in latest.items():
        if not r.get("operon") and not r.get("promoter"):
            empty.append(uid)
    return empty


def _pick_proxy(refs: list[str], embls: list[str]) -> tuple[Optional[str], str]:
    if refs:
        return refs[0], "refseq"
    if embls:
        return embls[0], "embl"
    return None, "none"


# --- IPG batching --------------------------------------------------------

def _ipg_call(accs: list[str], timeout: float) -> tuple[dict[str, dict], str]:
    if not accs:
        return {}, "ok"
    url = IPG_URL.format(ids=",".join(accs))
    try:
        resp = ncbi_get(url, timeout=timeout)
    except Exception as exc:
        return {}, f"exc:{type(exc).__name__}"
    if not resp.ok:
        return {}, f"http_{resp.status_code}"
    try:
        parsed = xmltodict.parse(resp.text)
    except Exception as exc:
        return {}, f"parse:{type(exc).__name__}"
    reports = parsed.get("IPGReportSet", {}).get("IPGReport", [])
    if isinstance(reports, dict):
        reports = [reports]
    out: dict[str, dict] = {}
    for r in reports:
        if isinstance(r, dict):
            acc = r.get("@product_acc")
            if acc:
                out[acc] = r
    return out, "ok"


def ipg_batch_with_retry(
    accs: list[str], max_attempts: int = 3, timeout: float = 120.0
) -> dict[str, dict]:
    """Fetch IPG reports for `accs`. Re-issue calls for any accessions the
    server silently omitted (common under concurrency)."""
    if not accs:
        return {}
    pending = set(accs)
    result: dict[str, dict] = {}
    for attempt in range(max_attempts):
        if not pending:
            break
        got, _status = _ipg_call(list(pending), timeout=timeout)
        result.update(got)
        pending -= set(got.keys())
        if pending and attempt < max_attempts - 1:
            time.sleep(1.0 + attempt)  # gentle backoff
    return result


def _coords_from_report(rep: dict) -> Optional[dict]:
    if "ProteinList" not in rep:
        return None
    try:
        protein = rep["ProteinList"]["Protein"]
        if isinstance(protein, list):
            protein = protein[0]
        cds = protein["CDSList"]["CDS"]
        if isinstance(cds, list):
            cds = cds[0]
        return {
            "Genome": cds["@accver"],
            "Start": cds["@start"],
            "Stop": cds["@stop"],
            "Strand": cds["@strand"],
        }
    except (KeyError, TypeError, IndexError):
        return None


# --- Per-record processing ----------------------------------------------

def _build_record(
    uid: str,
    proxy_acc: Optional[str],
    proxy_kind: str,
    rep: Optional[dict],
    promoter_params: dict,
) -> dict:
    """Run the post-IPG pipeline for a single record. Returns the per-member
    JSONL record shape (matches the flat schema in members_predict.py)."""
    base = {
        "uniprot_id": uid,
        "promoter": None,
        "source": "tsv_proxy",
        "computed_at": datetime.now(timezone.utc).isoformat(),
        "proxy_accession": proxy_acc,
        "proxy_kind": proxy_kind,
    }
    if rep is None or proxy_acc is None:
        return base
    coords = _coords_from_report(rep)
    if not coords:
        return base
    h = {"Uniprot Id": uid, **coords}
    try:
        op = acc2operon(h)
    except Exception:
        op = None
    if not op:
        return {**base, "genome": coords["Genome"]}
    try:
        prom = fetch_promoter(op, promoter_params)
    except Exception:
        prom = None
    return {
        **base,
        "genome": op.get("genome"),
        "promoter": prom,
        "protein_index": op.get("protein_index"),
        "operon": op.get("operon"),
    }


def _process_batch(
    items: list[tuple[str, str, str]],
    promoter_params: dict,
    max_attempts: int,
    timeout: float,
) -> list[dict]:
    """`items`: list of (uniprot_id, proxy_acc, proxy_kind). Returns one
    record per input. The per-record acc2operon + fetch_promoter calls run
    in a small thread pool inside the batch so a slow operon lookup doesn't
    serialise the whole batch."""
    proxies = [it[1] for it in items if it[1]]
    rep_by_acc = ipg_batch_with_retry(proxies, max_attempts=max_attempts, timeout=timeout)

    def one(item):
        uid, proxy, kind = item
        rep = rep_by_acc.get(proxy) if proxy else None
        return _build_record(uid, proxy, kind, rep, promoter_params)

    if _INNER_OPERON_WORKERS > 1 and len(items) > 1:
        with ThreadPoolExecutor(max_workers=_INNER_OPERON_WORKERS) as pool:
            return list(pool.map(one, items))
    return [one(it) for it in items]


# --- Main entrypoint ----------------------------------------------------

def backfill_via_tsv(
    jsonl_path: Path,
    tsv_path: Path,
    workers: int = 4,
    batch_size: int = 10,
    max_attempts: int = 3,
    timeout: float = 120.0,
    max_empties: Optional[int] = None,
    on_progress: Optional[Callable[[int, int], None]] = None,
) -> int:
    """End-to-end backfill. Appends new records (source="tsv_proxy") to
    `jsonl_path`. Returns the number of records appended."""
    log.info("loading TSV mapping from %s", tsv_path)
    mapping = load_tsv(tsv_path)
    log.info("  %d entries", len(mapping))

    log.info("scanning %s for empty records", jsonl_path)
    empties = find_empty_ids(jsonl_path)
    log.info("  %d empties found", len(empties))

    items: list[tuple[str, str, str]] = []
    kinds = {"refseq": 0, "embl": 0, "none": 0}
    for uid in empties:
        refs, embls = mapping.get(uid, ([], []))
        proxy, kind = _pick_proxy(refs, embls)
        kinds[kind] += 1
        if proxy:
            items.append((uid, proxy, kind))
    log.info("proxy types: %s — will process %d", kinds, len(items))

    if max_empties is not None:
        items = items[:max_empties]
        log.info("limited to %d for this run", len(items))

    if not items:
        return 0

    promoter_params = PromoterParams().model_dump()
    batches = [items[i : i + batch_size] for i in range(0, len(items), batch_size)]
    log.info("running %d batches across %d workers", len(batches), workers)

    write_lock = threading.Lock()
    appended = 0
    completed = 0

    with jsonl_path.open("a") as fh, ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(_process_batch, b, promoter_params, max_attempts, timeout): b
            for b in batches
        }
        for fut in as_completed(futures):
            try:
                records = fut.result()
            except Exception as exc:
                batch = futures[fut]
                log.exception("batch crashed: %s", exc)
                records = [
                    _build_record(uid, proxy, kind, None, promoter_params)
                    for uid, proxy, kind in batch
                ]
            with write_lock:
                for rec in records:
                    fh.write(json.dumps(rec) + "\n")
                fh.flush()
                appended += len(records)
                completed += 1
                if on_progress and (completed % 50 == 0 or completed == len(batches)):
                    on_progress(completed, len(batches))

    log.info("appended %d records to %s", appended, jsonl_path)
    return appended
