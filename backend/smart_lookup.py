"""Smart-lookup pre-flight for /api/predict.

Before running the full Snowprint pipeline, BLAST the input sequence against
each family's `representatives.dmnd` (built by PR 2's precompute). If the top
hit clears that family's identity / coverage thresholds, return the cached
prediction for the matched representative wrapped with match metadata. The
client treats this just like a `cached` event, except the banner also shows
the matched UniProt ID and the family's characterisation evidence.

Skipping this is controlled by `req.force`: when `force=True` the user wants a
fresh full pipeline run, bypassing both the smart-lookup and the
parameter-hash cache.
"""

from __future__ import annotations

import logging
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.blast import accID2sequence, uniprotID2sequence  # noqa: E402

from .families import Family, FamilyRegistry, registry as default_registry
from .schemas import PredictRequest

log = logging.getLogger(__name__)


@dataclass
class SmartLookupHit:
    """A successful smart-lookup, wrapping the cached prediction record."""

    family_key: str
    uniprot_id: str
    identity_pct: float
    coverage_pct: float
    record: dict  # the JSONL line for the matched representative


def _resolve_query_sequence(req: PredictRequest) -> Optional[str]:
    """Return the protein sequence implied by `req`, or None on lookup failure."""
    if req.input_method == "RefSeq":
        return accID2sequence(req.input_value)
    if req.input_method == "Uniprot":
        return uniprotID2sequence(req.input_value)
    return req.input_value


def _diamond_top_hit(
    query_seq: str,
    dmnd_path: Path,
    diamond_bin: str = "diamond",
) -> Optional[dict]:
    """Run `diamond blastp` for a single query sequence; return the best hit
    as `{uniprot_id, identity_pct, coverage_pct}` or None if there are no hits.
    """
    if not dmnd_path.exists():
        return None

    query_path = None
    output_path = None
    try:
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
            fh.write(f">query\n{query_seq}\n")
            query_path = Path(fh.name)
        output_path = Path(tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name)

        cmd = [
            diamond_bin,
            "blastp",
            "-d",
            str(dmnd_path.with_suffix("")),
            "-q",
            str(query_path),
            "-o",
            str(output_path),
            "--outfmt",
            "6",
            "sseqid",
            "pident",
            "qcovhsp",
            "--max-target-seqs",
            "1",
            "--quiet",
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, timeout=60)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            log.warning("smart-lookup diamond failed for %s: %s", dmnd_path, exc)
            return None

        line = output_path.read_text().splitlines()
        if not line:
            return None
        first = line[0].strip()
        if not first:
            return None
        sseqid, pident, qcovhsp = first.split("\t")[:3]
        return {
            "uniprot_id": sseqid,
            "identity_pct": float(pident),
            "coverage_pct": float(qcovhsp),
        }
    finally:
        for p in (query_path, output_path):
            if p is not None and p.exists():
                try:
                    os.unlink(p)
                except OSError:
                    pass


def _meets_thresholds(hit: dict, family: Family) -> bool:
    thresholds = family.smart_lookup or {}
    min_identity = thresholds.get("min_identity_pct", 50.0)
    min_coverage = thresholds.get("min_coverage_pct", 95.0)
    return hit["identity_pct"] >= min_identity and hit["coverage_pct"] >= min_coverage


def lookup(
    req: PredictRequest,
    registry: FamilyRegistry = default_registry,
) -> Optional[SmartLookupHit]:
    """Try every configured family in registration order; return the first hit.

    Returns None when:
      * `req.force` is true (caller asked for a fresh run),
      * the input sequence can't be resolved,
      * no family's `representatives.dmnd` exists yet,
      * the top hit doesn't clear that family's thresholds,
      * the matched UniProt ID isn't in the family's predictions.jsonl.
    """
    if req.force:
        return None

    seq = _resolve_query_sequence(req)
    if not seq:
        return None

    for family in registry.list():
        if not family.representatives_dmnd.exists():
            continue

        # diamond runs from the project root by convention (matches the rest
        # of the pipeline so its database paths resolve consistently).
        prev_cwd = os.getcwd()
        os.chdir(PROJECT_ROOT)
        try:
            hit = _diamond_top_hit(seq, family.representatives_dmnd)
        finally:
            os.chdir(prev_cwd)
        if hit is None:
            continue
        if not _meets_thresholds(hit, family):
            continue

        record = registry.lookup(family.key, hit["uniprot_id"])
        if record is None:
            log.warning(
                "smart-lookup hit %s in family %s but no record in predictions.jsonl",
                hit["uniprot_id"],
                family.key,
            )
            continue

        return SmartLookupHit(
            family_key=family.key,
            uniprot_id=hit["uniprot_id"],
            identity_pct=hit["identity_pct"],
            coverage_pct=hit["coverage_pct"],
            record=record,
        )
    return None


def to_event_payload(hit: SmartLookupHit) -> dict:
    """Format a SmartLookupHit as the `result` block emitted in the SSE event."""
    return {
        "matched_via": {
            "family": hit.family_key,
            "uniprot_id": hit.uniprot_id,
            "identity": round(hit.identity_pct / 100.0, 4),
            "coverage": round(hit.coverage_pct / 100.0, 4),
        },
        "result": hit.record,
    }
