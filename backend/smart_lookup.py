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


def _resolve_query_sequence(req: PredictRequest, registry: FamilyRegistry) -> Optional[str]:
    """Return the protein sequence implied by `req`. Tries each family's local
    sequence index first (FASTA byte-offset for UniProt input, RefSeq → UniProt
    cross-ref then FASTA for RefSeq input); falls back to NCBI/UniProt fetch
    only when the ID isn't in any family's precomputed members."""
    from . import sequence_index  # local import to avoid src.* import cycle
    return sequence_index.resolve_query_sequence(
        req,
        registry,
        fallback_refseq=accID2sequence,
        fallback_uniprot=uniprotID2sequence,
    )


def _normalise_uniprot_id(raw: str) -> str:
    """Strip a UniProt 'sp|P12345|NAME' / 'tr|P12345|NAME' style header to its
    bare accession. Diamond reports sseqid as whatever appears after `>` in
    the database's FASTA — for UniProt-sourced DBs that's the full pipe-
    delimited form, which doesn't match the bare accessions used as keys in
    our JSONL records."""
    if "|" in raw:
        parts = raw.split("|")
        if len(parts) >= 2 and parts[0] in ("sp", "tr"):
            return parts[1]
    return raw


def _diamond_hits(
    query_seq: str,
    dmnd_path: Path,
    diamond_bin: str = "diamond",
    max_target_seqs: int = 25,
) -> list[dict]:
    """Run `diamond blastp` for a single query sequence; return up to
    `max_target_seqs` hits in diamond's bit-score order as a list of
    `{uniprot_id, identity_pct, coverage_pct}` dicts. Empty list on
    failure or no hits — the caller filters by family thresholds.
    """
    if not dmnd_path.exists():
        return []

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
            str(max_target_seqs),
            "--quiet",
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, timeout=60)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            log.warning("smart-lookup diamond failed for %s: %s", dmnd_path, exc)
            return []

        hits: list[dict] = []
        for line in output_path.read_text().splitlines():
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                hits.append(
                    {
                        "uniprot_id": _normalise_uniprot_id(parts[0]),
                        "identity_pct": float(parts[1]),
                        "coverage_pct": float(parts[2]),
                    }
                )
            except ValueError:
                continue
        return hits
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


def _synthesize_record(
    family_key: str,
    qualifying_hits: list[dict],
    member_records: list[dict],
) -> dict:
    """Build a result record from the per-member predictions of the qualifying
    diamond hits. Shape matches the per-rep `predictions.jsonl` records that
    the existing frontend already knows how to render:

        {uniprot_id, centroid_uniprot_id, cluster, evidence, computed_at,
         input, protein_info, homologs: [{uniprot_id, identity, coverage,
                                          genome, start, stop, strand,
                                          operon, promoter}, ...]}

    `homologs` come from member_records (operon + promoter), enriched with
    the per-query identity/coverage from the diamond hits.
    """
    top_id = qualifying_hits[0]["uniprot_id"]
    homologs: list[dict] = []
    for hit, rec in zip(qualifying_hits, member_records):
        homologs.append(
            {
                "uniprot_id": hit["uniprot_id"],
                "identity": hit["identity_pct"],
                "coverage": hit["coverage_pct"],
                "genome": rec.get("genome"),
                "promoter": rec.get("promoter"),
                "protein_index": rec.get("protein_index"),
                "operon": rec.get("operon"),
            }
        )
    return {
        "uniprot_id": top_id,
        "centroid_uniprot_id": top_id,
        "cluster": None,
        "evidence": {"source": "members_db"},
        "computed_at": None,
        "input": None,
        "protein_info": None,
        "homologs": homologs,
    }


def lookup(
    req: PredictRequest,
    registry: FamilyRegistry = default_registry,
) -> Optional[SmartLookupHit]:
    """Try every configured family in registration order; return the first hit.

    Returns None when:
      * `req.force` is true (caller asked for a fresh run),
      * the input sequence can't be resolved,
      * no family's `members_with_promoters.dmnd` exists yet,
      * no diamond hit clears that family's thresholds AND has a record in
        `members_predictions.jsonl`.

    The smart-lookup index is now per-member (built from
    `members_predictions.jsonl`), not per-cluster-rep. We collect every
    diamond hit clearing the family's thresholds, pull each one's cached
    operon+promoter, and synthesize a result record whose `homologs` list
    feeds the operator finder directly. The top hit (by diamond bit score)
    becomes the `matched_via` anchor for the UI banner.
    """
    if req.force:
        return None

    seq = _resolve_query_sequence(req, registry)
    if not seq:
        return None

    for family in registry.list():
        if not family.members_with_promoters_dmnd.exists():
            continue

        max_homologs = int(family.smart_lookup.get("max_homologs", 30)) if family.smart_lookup else 30

        # diamond runs from the project root by convention (matches the rest
        # of the pipeline so its database paths resolve consistently).
        prev_cwd = os.getcwd()
        os.chdir(PROJECT_ROOT)
        try:
            hits = _diamond_hits(
                seq,
                family.members_with_promoters_dmnd,
                max_target_seqs=max(max_homologs * 3, 50),
            )
        finally:
            os.chdir(prev_cwd)

        # Collect every hit that clears the family's thresholds AND has a
        # cached record. Stop once we have max_homologs of them.
        qualifying: list[dict] = []
        records: list[dict] = []
        for hit in hits:
            if not _meets_thresholds(hit, family):
                continue
            rec = registry.lookup_member(family.key, hit["uniprot_id"])
            if rec is None:
                log.warning(
                    "smart-lookup hit %s in family %s but no record in members_predictions.jsonl",
                    hit["uniprot_id"],
                    family.key,
                )
                continue
            qualifying.append(hit)
            records.append(rec)
            if len(qualifying) >= max_homologs:
                break

        if not qualifying:
            continue

        anchor = qualifying[0]
        record = _synthesize_record(family.key, qualifying, records)
        return SmartLookupHit(
            family_key=family.key,
            uniprot_id=anchor["uniprot_id"],
            identity_pct=anchor["identity_pct"],
            coverage_pct=anchor["coverage_pct"],
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
