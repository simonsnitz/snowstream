"""Stage: fetch all UniProt IDs annotated with a given InterPro entry.

Output: families/<key>/members.txt — one UniProt accession per line.
Resume: skip if the output file already exists (use --force to override).

InterPro's REST API paginates with a `next` URL; we follow it until exhausted.
The API can be slow/flaky on large families, so we retry transient errors.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Iterable

import requests

log = logging.getLogger(__name__)

INTERPRO_API_BASE = "https://www.ebi.ac.uk/interpro/api"


def _iter_members(interpro_id: str, page_size: int = 200, http=requests) -> Iterable[str]:
    """Yield UniProt accessions for every protein annotated with `interpro_id`."""
    url = f"{INTERPRO_API_BASE}/protein/UniProt/entry/InterPro/{interpro_id}/?page_size={page_size}"
    while url:
        response = _get_with_retries(url, http=http)
        data = response.json()
        for entry in data.get("results", []):
            acc = entry.get("metadata", {}).get("accession")
            if acc:
                yield acc
        url = data.get("next")


def _get_with_retries(url: str, http=requests, max_attempts: int = 5, base_delay: float = 1.0):
    """GET with exponential backoff for 5xx and connection errors."""
    last_exc: Exception | None = None
    for attempt in range(max_attempts):
        try:
            response = http.get(url, timeout=60)
            if response.status_code < 500:
                response.raise_for_status()
                return response
            last_exc = requests.HTTPError(f"{response.status_code} for {url}")
        except (requests.ConnectionError, requests.Timeout, requests.HTTPError) as exc:
            last_exc = exc
        delay = base_delay * (2**attempt)
        log.warning("retry %d after %.1fs (last error: %s)", attempt + 1, delay, last_exc)
        time.sleep(delay)
    raise RuntimeError(f"InterPro request failed after {max_attempts} attempts: {last_exc}")


def fetch_members(interpro_id: str, output_path: Path, force: bool = False, http=requests) -> int:
    """Materialise members.txt for the family. Returns the count written."""
    if output_path.exists() and not force:
        existing = output_path.read_text().strip().splitlines()
        log.info("members.txt already exists (%d entries); skipping", len(existing))
        return len(existing)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = output_path.with_suffix(output_path.suffix + ".tmp")
    count = 0
    with tmp.open("w") as fh:
        for acc in _iter_members(interpro_id, http=http):
            fh.write(acc + "\n")
            count += 1
            if count % 1000 == 0:
                log.info("fetched %d members so far", count)
    tmp.replace(output_path)
    log.info("wrote %d members to %s", count, output_path)
    return count
