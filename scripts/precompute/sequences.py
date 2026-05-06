"""Stage: fetch FASTA sequences for the UniProt IDs in members.txt.

Output: families/<key>/members.fasta
Resume: skip if output exists. Within a single run, batches are streamed to a
        .tmp file and atomically renamed at the end.

UniProt's stream endpoint accepts an OR-joined accession query of arbitrary
length, but the URL has practical limits, so we batch ~250 IDs per request.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Iterable

import requests

log = logging.getLogger(__name__)

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
BATCH_SIZE = 250


def _batches(items: list[str], n: int) -> Iterable[list[str]]:
    for i in range(0, len(items), n):
        yield items[i : i + n]


def _fetch_batch(accessions: list[str], http=requests, max_attempts: int = 5) -> str:
    query = " OR ".join(f"accession:{a}" for a in accessions)
    last_exc: Exception | None = None
    for attempt in range(max_attempts):
        try:
            response = http.get(
                UNIPROT_STREAM,
                params={"query": query, "format": "fasta"},
                timeout=120,
            )
            if response.status_code < 500:
                response.raise_for_status()
                return response.text
            last_exc = requests.HTTPError(f"{response.status_code}")
        except (requests.ConnectionError, requests.Timeout, requests.HTTPError) as exc:
            last_exc = exc
        delay = 2.0 * (2**attempt)
        log.warning("uniprot batch retry %d after %.1fs (%s)", attempt + 1, delay, last_exc)
        time.sleep(delay)
    raise RuntimeError(f"UniProt stream failed after {max_attempts} attempts: {last_exc}")


def fetch_sequences(
    members_path: Path, output_path: Path, force: bool = False, http=requests
) -> int:
    """Fetch sequences for every UniProt ID in members_path → write a single FASTA.

    Returns the number of sequences written.
    """
    if output_path.exists() and not force:
        # Cheap-and-cheerful count of '>' lines.
        n = sum(1 for line in output_path.open() if line.startswith(">"))
        log.info("members.fasta already exists (%d sequences); skipping", n)
        return n

    if not members_path.exists():
        raise FileNotFoundError(f"missing prerequisite: {members_path}")
    accessions = [a.strip() for a in members_path.read_text().splitlines() if a.strip()]
    log.info("fetching %d sequences in batches of %d", len(accessions), BATCH_SIZE)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = output_path.with_suffix(output_path.suffix + ".tmp")
    written = 0
    with tmp.open("w") as fh:
        for i, batch in enumerate(_batches(accessions, BATCH_SIZE)):
            fasta = _fetch_batch(batch, http=http)
            fh.write(fasta)
            if not fasta.endswith("\n"):
                fh.write("\n")
            written += fasta.count(">")
            if (i + 1) % 10 == 0:
                log.info("batch %d: %d sequences so far", i + 1, written)
    tmp.replace(output_path)
    log.info("wrote %d sequences to %s", written, output_path)
    return written
