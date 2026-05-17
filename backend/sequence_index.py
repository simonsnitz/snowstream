"""Local sequence resolution for query inputs.

Replaces the NCBI/UniProt round trip on every `/api/predict` submission with
a byte-offset lookup against the family's `members.fasta` (UniProt input) or
its `refseq_to_uniprot.json` cross-reference (RefSeq input). NCBI fetch is
retained as a fallback for IDs that aren't in any family's index — e.g. when
a user queries a protein outside our precompute scope.

Indexes are loaded lazily on first use and cached at module level. The FASTA
byte-offset scan is a single-pass over the file at first access and produces
a `{uniprot_id: offset}` dict (~5 MB RAM per 650k-member family).
"""

from __future__ import annotations

import json
import logging
import threading
from pathlib import Path
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from .families import Family, FamilyRegistry
    from .schemas import PredictRequest

log = logging.getLogger(__name__)


def _normalise_uniprot_id(header: str) -> Optional[str]:
    """Pull the bare accession out of a FASTA header line body
    (`sp|P12345|NAME ...` → `P12345`)."""
    body = header.strip()
    if "|" in body:
        parts = body.split("|")
        if len(parts) >= 2 and parts[0] in ("sp", "tr"):
            return parts[1]
    return body.split()[0] if body else None


class SequenceIndex:
    """Per-family sequence lookup. Thread-safe — uses fresh file handles
    on each `get_*` call, never shares an open file between calls."""

    def __init__(self, fasta_path: Path, refseq_map_path: Optional[Path] = None) -> None:
        self._fasta_path = fasta_path
        self._fasta_offsets: dict[str, int] = {}
        self._refseq_to_uid: dict[str, str] = {}
        self._lock = threading.Lock()
        self._loaded = False
        self._refseq_map_path = refseq_map_path

    def _ensure_loaded(self) -> None:
        if self._loaded:
            return
        with self._lock:
            if self._loaded:
                return
            self._fasta_offsets = self._scan_fasta(self._fasta_path)
            if self._refseq_map_path and self._refseq_map_path.exists():
                try:
                    self._refseq_to_uid = json.loads(self._refseq_map_path.read_text())
                except (json.JSONDecodeError, OSError) as exc:
                    log.warning("could not read %s: %s", self._refseq_map_path, exc)
            self._loaded = True
            log.info(
                "sequence index loaded: %d uniprot, %d refseq mappings",
                len(self._fasta_offsets),
                len(self._refseq_to_uid),
            )

    @staticmethod
    def _scan_fasta(fasta_path: Path) -> dict[str, int]:
        offsets: dict[str, int] = {}
        if not fasta_path.exists():
            return offsets
        with fasta_path.open("rb") as fh:
            offset = fh.tell()
            line = fh.readline()
            while line:
                if line.startswith(b">"):
                    header = line[1:].decode("ascii", errors="ignore")
                    uid = _normalise_uniprot_id(header)
                    if uid:
                        offsets[uid] = offset
                offset = fh.tell()
                line = fh.readline()
        return offsets

    def get_by_uniprot(self, uid: str) -> Optional[str]:
        self._ensure_loaded()
        offset = self._fasta_offsets.get(uid)
        if offset is None:
            return None
        with self._fasta_path.open("rb") as fh:
            fh.seek(offset)
            fh.readline()  # consume header
            chunks: list[str] = []
            for line in fh:
                if line.startswith(b">"):
                    break
                chunks.append(line.strip().decode("ascii", errors="ignore"))
        return "".join(chunks) or None

    def get_by_refseq(self, acc: str) -> Optional[str]:
        self._ensure_loaded()
        uid = self._refseq_to_uid.get(acc)
        if uid is None:
            return None
        return self.get_by_uniprot(uid)


# Module-level cache, keyed by family.key.
_indexes: dict[str, SequenceIndex] = {}
_indexes_lock = threading.Lock()


def get_index(family: "Family") -> SequenceIndex:
    with _indexes_lock:
        idx = _indexes.get(family.key)
        if idx is None:
            idx = SequenceIndex(
                fasta_path=family.dir / "members.fasta",
                refseq_map_path=family.dir / "refseq_to_uniprot.json",
            )
            _indexes[family.key] = idx
        return idx


def resolve_query_sequence(
    req: "PredictRequest",
    registry: "FamilyRegistry",
    fallback_refseq: Optional[callable] = None,
    fallback_uniprot: Optional[callable] = None,
) -> Optional[str]:
    """Resolve `req.input_value` to a protein sequence using local indexes
    first, falling back to the provided NCBI/UniProt fetchers if no family
    has the ID.

    `fallback_refseq` and `fallback_uniprot` are the existing
    `accID2sequence` and `uniprotID2sequence` callables — passed in to keep
    this module from importing src.* (which drags in streamlit at import
    time)."""
    if req.input_method == "Protein sequence":
        return req.input_value

    for family in registry.list():
        idx = get_index(family)
        if req.input_method == "Uniprot":
            seq = idx.get_by_uniprot(req.input_value)
        elif req.input_method == "RefSeq":
            seq = idx.get_by_refseq(req.input_value)
        else:
            seq = None
        if seq:
            return seq

    if req.input_method == "RefSeq" and fallback_refseq is not None:
        return fallback_refseq(req.input_value)
    if req.input_method == "Uniprot" and fallback_uniprot is not None:
        return fallback_uniprot(req.input_value)
    return req.input_value
