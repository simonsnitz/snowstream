"""Family-level precomputed cache.

Each family lives at `families/<key>/` and is described by `family.json`. Its
`predictions.jsonl` is a streaming-friendly per-line record of cached
predictions, one per representative protein, keyed by UniProt ID.

Lookups are O(1) after a one-pass index build: we scan the JSONL once at
load time, parsing only enough of each line to extract `uniprot_id`, and
remember the byte offset. A direct lookup then seeks to that offset and
parses the full record. Re-loading is required when the JSONL grows, so
in dev we re-stat the file's mtime on each lookup; in prod the index can
be rebuilt explicitly.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

FAMILIES_DIR = Path(__file__).resolve().parent.parent / "families"


@dataclass
class Family:
    key: str
    name: str
    interpro_id: str
    groovdb_entry_prefix: str
    cluster_identity_threshold: float
    smart_lookup: dict
    dir: Path
    predictions_path: Path
    representatives_fasta: Path
    representatives_dmnd: Path
    # Per-member dataset (newer path: smart-lookup BLASTs against this).
    members_predictions_path: Path
    members_with_promoters_dmnd: Path

    def manifest_dict(self) -> dict:
        """Public manifest view for `GET /api/families` (no filesystem paths)."""
        return {
            "key": self.key,
            "name": self.name,
            "interpro_id": self.interpro_id,
            "cluster_identity_threshold": self.cluster_identity_threshold,
            "smart_lookup": self.smart_lookup,
            "predictions_count": _count_predictions(self.predictions_path),
        }


@dataclass
class _Index:
    mtime: float
    offsets: dict[str, int]  # uniprot_id -> byte offset of line


def _count_predictions(path: Path) -> int:
    if not path.exists():
        return 0
    n = 0
    with path.open("rb") as fh:
        for _ in fh:
            n += 1
    return n


def _build_index(path: Path) -> _Index:
    """One pass through the JSONL, recording byte offsets per uniprot_id.

    Only parses each line far enough to extract uniprot_id; the rest of
    the record is parsed lazily on lookup.
    """
    offsets: dict[str, int] = {}
    if not path.exists():
        return _Index(mtime=0.0, offsets=offsets)
    mtime = path.stat().st_mtime
    with path.open("rb") as fh:
        offset = fh.tell()
        line = fh.readline()
        while line:
            try:
                # The records are small enough that parsing the whole
                # line is cheaper than a regex / partial parser, and we
                # avoid mis-matching any field-named 'uniprot_id' inside
                # nested data.
                rec = json.loads(line)
                uid = rec.get("uniprot_id")
                if uid:
                    offsets[uid] = offset
            except (json.JSONDecodeError, ValueError):
                pass
            offset = fh.tell()
            line = fh.readline()
    return _Index(mtime=mtime, offsets=offsets)


class FamilyRegistry:
    """In-memory registry of families and their JSONL byte-offset indexes.

    Two parallel indexes per family — one for the legacy per-rep
    `predictions.jsonl`, one for the newer per-member
    `members_predictions.jsonl`. Both are rebuilt lazily on mtime change.
    """

    def __init__(self, root: Path = FAMILIES_DIR) -> None:
        self.root = root
        self._families: dict[str, Family] = {}
        self._indexes: dict[str, _Index] = {}
        self._member_indexes: dict[str, _Index] = {}
        self._reload()

    def _reload(self) -> None:
        self._families.clear()
        self._indexes.clear()
        self._member_indexes.clear()
        if not self.root.exists():
            return
        for child in sorted(self.root.iterdir()):
            if not child.is_dir() or child.name.startswith("_"):
                continue
            manifest_path = child / "family.json"
            if not manifest_path.exists():
                continue
            try:
                manifest = json.loads(manifest_path.read_text())
            except (json.JSONDecodeError, OSError):
                continue
            fam = Family(
                key=manifest["key"],
                name=manifest["name"],
                interpro_id=manifest["interpro_id"],
                groovdb_entry_prefix=manifest.get("groovdb_entry_prefix", ""),
                cluster_identity_threshold=float(manifest.get("cluster_identity_threshold", 0.5)),
                smart_lookup=manifest.get("smart_lookup", {}),
                dir=child,
                predictions_path=child / manifest.get("predictions_path", "predictions.jsonl"),
                representatives_fasta=child / manifest.get("representatives_fasta", "representatives.fasta"),
                representatives_dmnd=child / manifest.get("representatives_dmnd", "representatives.dmnd"),
                members_predictions_path=child / manifest.get("members_predictions_path", "members_predictions.jsonl"),
                members_with_promoters_dmnd=child / manifest.get("members_with_promoters_dmnd", "members_with_promoters.dmnd"),
            )
            self._families[fam.key] = fam

    def list(self) -> list[Family]:
        return list(self._families.values())

    def get(self, key: str) -> Optional[Family]:
        return self._families.get(key)

    def lookup(self, family_key: str, uniprot_id: str) -> Optional[dict[str, Any]]:
        fam = self.get(family_key)
        if not fam:
            return None
        index = self._fresh_index(fam)
        offset = index.offsets.get(uniprot_id)
        if offset is None:
            return None
        with fam.predictions_path.open("rb") as fh:
            fh.seek(offset)
            line = fh.readline()
        try:
            return json.loads(line)
        except (json.JSONDecodeError, ValueError):
            return None

    def _fresh_index(self, fam: Family) -> _Index:
        """Return an index, rebuilding it if the JSONL has changed on disk."""
        cur_mtime = fam.predictions_path.stat().st_mtime if fam.predictions_path.exists() else 0.0
        index = self._indexes.get(fam.key)
        if index is None or index.mtime != cur_mtime:
            index = _build_index(fam.predictions_path)
            self._indexes[fam.key] = index
        return index

    def lookup_member(self, family_key: str, uniprot_id: str) -> Optional[dict[str, Any]]:
        """Per-member lookup against `members_predictions.jsonl`. Same
        byte-offset / mtime-invalidation pattern as the per-rep `lookup`."""
        fam = self.get(family_key)
        if not fam:
            return None
        index = self._fresh_member_index(fam)
        offset = index.offsets.get(uniprot_id)
        if offset is None:
            return None
        with fam.members_predictions_path.open("rb") as fh:
            fh.seek(offset)
            line = fh.readline()
        try:
            return json.loads(line)
        except (json.JSONDecodeError, ValueError):
            return None

    def _fresh_member_index(self, fam: Family) -> _Index:
        cur_mtime = (
            fam.members_predictions_path.stat().st_mtime
            if fam.members_predictions_path.exists()
            else 0.0
        )
        index = self._member_indexes.get(fam.key)
        if index is None or index.mtime != cur_mtime:
            index = _build_index(fam.members_predictions_path)
            self._member_indexes[fam.key] = index
        return index

    def reload(self) -> None:
        """Re-scan the families/ directory (e.g. after a new family is added)."""
        self._reload()


# Module-level singleton, lazily initialised by FastAPI on first import.
registry = FamilyRegistry()
