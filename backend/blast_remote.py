"""Remote BLAST against NCBI's nr database via the BLAST URL API.

Three steps: PUT (submit) → poll status → GET XML results. The job typically
takes 5–30 minutes so the caller is expected to provide a progress callback
that the polling loop calls every iteration.

NCBI nr hits land as RefSeq / GenBank / EMBL accessions, not UniProt
accessions. The downstream Snowprint pipeline (uniprot2EMBL, fetch_promoter,
etc.) expects UniProt IDs, so we translate via UniProt's xref query and drop
hits that don't map.

References:
- https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html
- https://www.uniprot.org/help/api_queries (xref:RefSeq:xxx)
"""

from __future__ import annotations

import logging
import re
import time
from io import StringIO
from typing import Callable, Optional

import pandas as pd
import requests
from Bio.Blast import NCBIXML

log = logging.getLogger(__name__)

NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

# How often to poll NCBI for job status, and how long to wait overall.
DEFAULT_POLL_INTERVAL = 30.0
DEFAULT_MAX_WAIT = 60 * 60.0  # 1 hour


ProgressCallback = Callable[[str, dict], None]


def _emit(on_progress: Optional[ProgressCallback], event: str, **payload) -> None:
    if on_progress:
        try:
            on_progress(event, payload)
        except Exception:  # never let a callback take the BLAST down
            log.exception("progress callback raised; continuing")


def submit(query_seq: str, http=requests) -> tuple[str, int]:
    """Submit a blastp-vs-nr job. Returns (RID, RTOE_seconds).

    RTOE is NCBI's estimated time-of-execution; we use it for the first sleep
    so we don't poll before the job has had a chance to start.
    """
    response = http.post(
        NCBI_BLAST_URL,
        data={
            "CMD": "Put",
            "PROGRAM": "blastp",
            "DATABASE": "nr",
            "QUERY": query_seq,
        },
        timeout=60,
    )
    response.raise_for_status()
    rid_match = re.search(r"RID\s*=\s*(\S+)", response.text)
    rtoe_match = re.search(r"RTOE\s*=\s*(\d+)", response.text)
    if not rid_match:
        raise RuntimeError(f"NCBI BLAST submit returned no RID; body starts: {response.text[:300]!r}")
    rid = rid_match.group(1)
    rtoe = int(rtoe_match.group(1)) if rtoe_match else 0
    return rid, rtoe


def poll(rid: str, http=requests) -> str:
    """Returns one of 'WAITING', 'READY', 'FAILED', 'UNKNOWN'."""
    response = http.get(
        NCBI_BLAST_URL,
        params={"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"},
        timeout=60,
    )
    response.raise_for_status()
    text = response.text
    if "Status=READY" in text:
        return "READY"
    if "Status=WAITING" in text:
        return "WAITING"
    if "Status=FAILED" in text:
        return "FAILED"
    if "Status=UNKNOWN" in text:
        return "UNKNOWN"
    return "UNKNOWN"


def fetch_xml(rid: str, http=requests) -> str:
    response = http.get(
        NCBI_BLAST_URL,
        params={"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"},
        timeout=180,
    )
    response.raise_for_status()
    return response.text


def parse_xml(xml: str, ident_cutoff: float, cov_cutoff: float) -> pd.DataFrame:
    """Parse a BLAST XML result into a (RefSeq-accession-keyed) DataFrame.

    Columns mirror what `src.blast.blast` produces from diamond:
        Uniprot Id, Identity, Coverage
    Despite the column name, the values here are the *NCBI accessions*; the
    caller translates to real UniProt IDs separately.
    """
    record = NCBIXML.read(StringIO(xml))
    rows: list[tuple[str, float, float]] = []
    query_len = record.query_length or 1
    for alignment in record.alignments:
        if not alignment.hsps:
            continue
        hsp = alignment.hsps[0]
        identity_pct = (hsp.identities / hsp.align_length) * 100 if hsp.align_length else 0.0
        coverage_pct = (hsp.align_length / query_len) * 100 if query_len else 0.0
        if identity_pct < ident_cutoff or coverage_pct < cov_cutoff:
            continue
        accession = alignment.accession or _extract_accession_from_id(alignment.hit_id)
        if accession:
            rows.append((accession, round(identity_pct, 2), round(coverage_pct, 2)))
    return pd.DataFrame(rows, columns=["Uniprot Id", "Identity", "Coverage"])


def _extract_accession_from_id(hit_id: str) -> str:
    """`hit_id` looks like 'gi|123|ref|WP_X.1|' — pull out the bare accession."""
    parts = hit_id.split("|")
    for db, acc in zip(parts[0::2], parts[1::2]):
        if db in ("ref", "sp", "tr", "gb", "emb"):
            return acc.rstrip(".")
    return parts[-1] if parts else hit_id


def translate_to_uniprot(accessions: list[str], http=requests) -> dict[str, str]:
    """Map RefSeq/GenBank accessions to UniProt accessions via UniProt xref query.

    Batched (~50 IDs per query) to keep URLs sane. Hits that can't be mapped
    are simply absent from the returned dict — the caller drops them.
    """
    if not accessions:
        return {}
    out: dict[str, str] = {}
    BATCH = 50
    seen = set()
    for i in range(0, len(accessions), BATCH):
        batch = [a for a in accessions[i : i + BATCH] if a not in seen]
        if not batch:
            continue
        seen.update(batch)
        # Strip version suffixes (".1") so the xref match isn't picky.
        stripped = [a.split(".")[0] for a in batch]
        query = " OR ".join(f"xref:{a}" for a in stripped)
        try:
            response = http.get(
                UNIPROT_SEARCH,
                params={
                    "query": query,
                    "fields": "accession,xref_refseq",
                    "format": "json",
                    "size": str(min(500, BATCH * 2)),
                },
                timeout=60,
            )
            if not response.ok:
                continue
            data = response.json()
        except (requests.RequestException, ValueError):
            continue
        for entry in data.get("results", []):
            uniprot_acc = entry.get("primaryAccession")
            if not uniprot_acc:
                continue
            for xref in entry.get("uniProtKBCrossReferences", []):
                if xref.get("database") != "RefSeq":
                    continue
                xref_id = (xref.get("id") or "").split(".")[0]
                # Match against any input accession with the same stem.
                for original in batch:
                    if original.split(".")[0] == xref_id and original not in out:
                        out[original] = uniprot_acc
                        break
    return out


def blast_against_nr(
    query_seq: str,
    params: dict,
    on_progress: Optional[ProgressCallback] = None,
    poll_interval: float = DEFAULT_POLL_INTERVAL,
    max_wait: float = DEFAULT_MAX_WAIT,
    http=requests,
) -> pd.DataFrame:
    """End-to-end remote BLAST → parsed DataFrame with UniProt-accession IDs.

    Returns the same columns/shape as `src.blast.blast` so the rest of the
    pipeline can consume it without branching.
    """
    rid, rtoe = submit(query_seq, http=http)
    _emit(on_progress, "submitted", rid=rid, rtoe=rtoe)
    log.info("BLAST submitted: rid=%s rtoe=%s", rid, rtoe)

    elapsed = 0.0
    initial = max(min(rtoe, poll_interval), 5.0)
    time.sleep(initial)
    elapsed += initial

    while elapsed < max_wait:
        status = poll(rid, http=http)
        _emit(on_progress, "polling", rid=rid, status=status, elapsed=int(elapsed))
        if status == "READY":
            break
        if status == "FAILED":
            raise RuntimeError(f"NCBI BLAST job {rid} failed")
        time.sleep(poll_interval)
        elapsed += poll_interval
    else:
        raise TimeoutError(f"NCBI BLAST job {rid} did not complete within {int(max_wait)}s")

    _emit(on_progress, "fetching", rid=rid)
    xml = fetch_xml(rid, http=http)
    df = parse_xml(xml, params["ident_cutoff"], params["cov_cutoff"])
    if df.empty:
        return df

    _emit(on_progress, "translating", count=len(df))
    accessions = df["Uniprot Id"].tolist()
    translation = translate_to_uniprot(accessions, http=http)
    df["Uniprot Id"] = df["Uniprot Id"].map(translation)
    df = df.dropna(subset=["Uniprot Id"]).reset_index(drop=True)
    _emit(on_progress, "translated", mapped=len(df), unmapped=len(accessions) - len(df))
    return df
