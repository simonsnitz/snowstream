"""groovDB lookup helper.

The groovDB dump is downloaded once (manually, via the button on
https://www.groov.bio/about/programmatic-access) and saved to
`families/_resources/groovdb.json`. This module loads it into a
`{uniprot_id: entry}` dict and exposes two helpers:

- `load(path)` — read the dump from disk
- `lookup(dump, uniprot_id)` — return the entry or None

The exact JSON shape from groovDB isn't fully specified upstream, so we try
a few common shapes (top-level list of records, or a top-level dict whose
values are records) and look for any of `uniprot_id`, `uniprotID`,
`uniprot`, or `accession` in each record.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Optional

log = logging.getLogger(__name__)


_UNIPROT_KEYS = ("uniprot_id", "uniprotID", "uniprot", "accession", "Uniprot Id")


def load(path: Path) -> dict[str, dict[str, Any]]:
    """Parse the groovDB dump into {uniprot_id: full_entry}."""
    if not path.exists():
        log.info("groovDB dump not found at %s; lookups will return None", path)
        return {}

    raw = json.loads(path.read_text())
    records: list[dict[str, Any]]
    if isinstance(raw, list):
        records = raw
    elif isinstance(raw, dict):
        # Either {uniprot: entry} or {category: [entries]}
        first_value = next(iter(raw.values()), None)
        if isinstance(first_value, dict):
            records = list(raw.values())
        elif isinstance(first_value, list):
            records = [item for sublist in raw.values() for item in sublist if isinstance(item, dict)]
        else:
            records = []
    else:
        records = []

    by_uniprot: dict[str, dict[str, Any]] = {}
    for rec in records:
        uid = _extract_uniprot(rec)
        if uid:
            by_uniprot[uid] = rec
    log.info("loaded %d groovDB entries from %s", len(by_uniprot), path)
    return by_uniprot


def _extract_uniprot(rec: dict[str, Any]) -> Optional[str]:
    for key in _UNIPROT_KEYS:
        val = rec.get(key)
        if isinstance(val, str) and val.strip():
            return val.strip()
    return None


def lookup(dump: dict[str, dict[str, Any]], uniprot_id: str) -> Optional[dict[str, Any]]:
    return dump.get(uniprot_id)


def make_url(entry_prefix: str, uniprot_id: str) -> str:
    """Construct a groovDB entry URL given the family-specific prefix."""
    return entry_prefix.rstrip("/") + "/" + uniprot_id
