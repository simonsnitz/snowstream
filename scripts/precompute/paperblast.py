"""PaperBLAST remote query + on-disk response cache.

PaperBLAST (papers.genomics.lbl.gov) takes a protein sequence and returns
matches against proteins mentioned in published papers. We use it to find
literature-backed cluster representatives.

Practical notes:

- The endpoint sits behind Cloudflare's browser challenge. Plain `requests`
  with browser-like headers usually gets through, but if Cloudflare flips
  to a JS challenge the request will return HTML for the challenge page
  instead of results. Callers should treat parse failure as a soft miss
  (no literature evidence) rather than a hard error.
- We cache responses on disk by SHA-1 of the query sequence under
  `families/_resources/paperblast_cache/` so reruns of the precompute
  don't re-query.
- Hits are parsed loosely from the HTML using regex; we extract the
  subject UniProt ID, percent identity, percent coverage, and any
  associated DOIs/years. The exact HTML format is fragile; if it changes
  upstream the parser returns no hits and the cluster falls back to its
  MMseqs2 centroid.
"""

from __future__ import annotations

import hashlib
import json
import logging
import re
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

import requests

log = logging.getLogger(__name__)

PAPERBLAST_URL = "https://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi"
USER_AGENT = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"
)


@dataclass
class Paper:
    doi: Optional[str]
    year: Optional[int]


@dataclass
class PaperBlastHit:
    subject_uniprot_id: str
    identity_pct: float
    coverage_pct: float
    papers: list[Paper]


def _query_hash(sequence: str) -> str:
    return hashlib.sha1(sequence.strip().upper().encode()).hexdigest()


def query(
    sequence: str,
    cache_dir: Path,
    http=requests,
    timeout: float = 60.0,
) -> list[PaperBlastHit]:
    """Query PaperBLAST with a protein sequence; return parsed hits.

    Empty list means either no hits, an unreachable server, or a Cloudflare
    challenge page — callers shouldn't distinguish.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{_query_hash(sequence)}.html"

    if cache_path.exists():
        html = cache_path.read_text()
    else:
        try:
            response = http.get(
                PAPERBLAST_URL,
                params={"query": sequence, "more": ""},
                headers={"User-Agent": USER_AGENT, "Accept": "text/html"},
                timeout=timeout,
            )
            if response.status_code != 200:
                log.warning("PaperBLAST returned HTTP %d", response.status_code)
                return []
            html = response.text
        except requests.RequestException as exc:
            log.warning("PaperBLAST request failed: %s", exc)
            return []
        cache_path.write_text(html)
        # Be polite — don't hammer the endpoint.
        time.sleep(1.0)

    if "Just a moment" in html or "Enable JavaScript" in html:
        log.warning("PaperBLAST returned a Cloudflare challenge; treating as no hits")
        return []
    return _parse(html)


# Crude HTML parser — PaperBLAST pages list each hit in a section that
# starts with the UniProt accession and percent identity, followed by
# linked papers. We extract by regex; if the format shifts upstream this
# returns [] and the cluster falls back to its default representative.

_HIT_BLOCK_RE = re.compile(
    r"\b(?P<acc>[A-Z][A-Z0-9]{5,9})\b[^\n]{0,200}?"
    r"(?P<identity>\d+(?:\.\d+)?)\s*%\s*identity[^\n]*?"
    r"(?P<coverage>\d+(?:\.\d+)?)\s*%\s*coverage",
    re.IGNORECASE,
)
_DOI_RE = re.compile(r"\b(10\.\d{4,9}/[^\s\"<>]+)", re.IGNORECASE)
_YEAR_RE = re.compile(r"\b(19|20)\d{2}\b")


def _parse(html: str) -> list[PaperBlastHit]:
    hits: list[PaperBlastHit] = []
    for match in _HIT_BLOCK_RE.finditer(html):
        # Look at a window of HTML following the hit for paper metadata.
        end = match.end()
        window = html[end : end + 4000]
        # End the window at the next hit start, if any.
        next_match = _HIT_BLOCK_RE.search(window)
        if next_match:
            window = window[: next_match.start()]

        papers = []
        seen_dois = set()
        for doi_match in _DOI_RE.finditer(window):
            doi = doi_match.group(1).rstrip(".,)\"")
            if doi in seen_dois:
                continue
            seen_dois.add(doi)
            year_match = _YEAR_RE.search(window[max(0, doi_match.start() - 200) : doi_match.start()])
            year = int(year_match.group(0)) if year_match else None
            papers.append(Paper(doi=doi, year=year))

        hits.append(
            PaperBlastHit(
                subject_uniprot_id=match.group("acc"),
                identity_pct=float(match.group("identity")),
                coverage_pct=float(match.group("coverage")),
                papers=papers,
            )
        )
    return hits


def to_jsonable(hit: PaperBlastHit) -> dict:
    return {
        "subject_uniprot_id": hit.subject_uniprot_id,
        "identity": hit.identity_pct,
        "coverage": hit.coverage_pct,
        "papers": [asdict(p) for p in hit.papers],
    }


def to_jsonable_list(hits: list[PaperBlastHit]) -> list[dict]:
    return [to_jsonable(h) for h in hits]
