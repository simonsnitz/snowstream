"""Backend pipeline: BLAST -> genome coordinates -> operons/promoters.

Yields SSE-shaped events as the pipeline runs. Wraps the existing
synchronous functions in `src/` with `asyncio.to_thread` so the event
loop isn't blocked.
"""
import asyncio
import json
import os
import sys
from pathlib import Path
from typing import AsyncIterator

import requests

# blast.py expects to find ../databases/ relative to its own location, but
# the diamond binary is invoked from the current working directory. So we
# chdir to the project root for the duration of the pipeline.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.blast import blast, accID2sequence, uniprotID2sequence  # noqa: E402
from src.get_genome_coordinates import (  # noqa: E402
    get_genome_coordinates,
    get_genome_coordinates_batch,
)
from src.accID2operon import acc2operon  # noqa: E402
from src.fetch_promoter import fetch_promoter  # noqa: E402
from src.troubleshoot import troubleshoot  # noqa: E402

from . import cache
from . import smart_lookup as smart_lookup_mod
from .blast_remote import blast_against_nr
from .schemas import PredictRequest


def _sse(event_type: str, **payload) -> str:
    return f"data: {json.dumps({'type': event_type, **payload})}\n\n"


def _filter_redundant(blast_df) -> list[dict]:
    homologs = []
    seen = set()
    for _, row in blast_df.iterrows():
        key = (row["Identity"], row["Coverage"])
        if key in seen:
            continue
        seen.add(key)
        homologs.append(
            {
                "Uniprot Id": row["Uniprot Id"],
                "identity": row["Identity"],
                "coverage": row["Coverage"],
            }
        )
    return homologs


def _all_homologs(blast_df) -> list[dict]:
    return [
        {
            "Uniprot Id": row["Uniprot Id"],
            "identity": row["Identity"],
            "coverage": row["Coverage"],
        }
        for _, row in blast_df.iterrows()
    ]


def _uniprot_protein_info(uniprot_id: str) -> dict | None:
    url = (
        f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
        "?format=json&fields=sequence,organism_name,protein_name"
    )
    try:
        response = requests.get(url, timeout=10)
        if not response.ok:
            return None
        data = response.json()
        return {
            "annotation": data["proteinDescription"]["recommendedName"]["fullName"]["value"],
            "organism": data["organism"]["scientificName"],
            "lineage": data["organism"]["lineage"],
        }
    except (requests.RequestException, KeyError, ValueError):
        return None


def _normalize_homolog(h: dict) -> dict:
    """Convert internal homolog dict (mixed casing) to API-friendly snake_case."""
    return {
        "uniprot_id": h.get("Uniprot Id"),
        "identity": h.get("identity"),
        "coverage": h.get("coverage"),
        "genome": h.get("Genome") or h.get("genome"),
        "start": h.get("start"),
        "stop": h.get("stop"),
        "strand": h.get("strand"),
        "operon": h.get("operon"),
        "promoter": h.get("promoter"),
    }


# --- Per-stage helpers (pure-sync; reused by both run_pipeline and
# run_pipeline_core). They each touch the network or call external binaries,
# so the SSE wrapper invokes them via asyncio.to_thread. ---


def _resolve_query_sequence(req: PredictRequest) -> str | None:
    """Fetch the protein sequence implied by `req`. Returns None on lookup failure."""
    if req.input_method == "RefSeq":
        return accID2sequence(req.input_value)
    if req.input_method == "Uniprot":
        return uniprotID2sequence(req.input_value)
    return req.input_value


def _blast_stage(req: PredictRequest, on_progress=None):
    """Run BLAST + UniProt info lookup. Returns (homologs, protein_info, error_or_none).

    Branches on `req.database`: the local Diamond DB (default) vs remote NCBI nr.
    """
    if req.database == "nr_remote":
        seq = _resolve_query_sequence(req)
        if not seq:
            return [], None, {"kind": "input_lookup_failed", "message": f"could not resolve sequence for {req.input_value}"}
        try:
            blast_df = blast_against_nr(
                seq,
                req.blast_params.model_dump(),
                on_progress=on_progress,
            )
        except (RuntimeError, TimeoutError) as exc:
            return [], None, {"kind": "remote_blast_failed", "message": str(exc)}
    else:
        blast_df = blast(
            req.input_value,
            req.input_method,
            req.blast_params.model_dump(),
            500,
        )
    if blast_df.empty:
        mode = None
        if req.input_method != "Protein sequence":
            mode, _ = troubleshoot(req.input_method, req.input_value)
        return (
            [],
            None,
            {
                "kind": mode or "no_results",
                "message": (
                    "Protein is not from Bacteria. Snowprint only works for bacterial proteins."
                    if mode == "not bacteria"
                    else "BLAST returned no results."
                ),
            },
        )
    homologs = (
        _filter_redundant(blast_df)
        if req.blast_params.filter_redundant
        else _all_homologs(blast_df)
    )
    homologs = homologs[: req.blast_params.max_homologs]
    protein_info = _uniprot_protein_info(blast_df.iloc[0]["Uniprot Id"])
    return homologs, protein_info, None


def _coords_one(homolog: dict) -> dict | None:
    return get_genome_coordinates(homolog)


def _coords_batch(homologs: list[dict]) -> list[dict] | None:
    return get_genome_coordinates_batch(homologs)


def _operon_and_promoter(homolog: dict, promoter_params: dict) -> dict:
    homolog["operon"] = acc2operon(homolog)
    try:
        homolog["promoter"] = fetch_promoter(homolog["operon"], promoter_params)
    except Exception:
        homolog["promoter"] = None
    return homolog


def run_pipeline_core(req: PredictRequest) -> dict:
    """Synchronous full-pipeline run. No SSE, no caching.

    Returns either {input, protein_info, homologs} on success, or
    {error: {stage, kind, message}} on failure. Used by the precompute
    script so it can drive the pipeline without parsing SSE or polluting
    the per-query parameter-hash cache.
    """
    prev_cwd = os.getcwd()
    os.chdir(PROJECT_ROOT)
    try:
        homologs, protein_info, err = _blast_stage(req)
        if err is not None:
            return {"error": {"stage": "blast", **err}}

        if req.get_coordinates_method == "batch":
            result = _coords_batch(homologs)
            if result is None:
                return {
                    "error": {
                        "stage": "coordinates",
                        "kind": "batch_failed",
                        "message": "Failed fetching genome coordinates in batch mode.",
                    }
                }
            homologs = [h for h in result if h is not None]
        else:
            updated = [_coords_one(h) for h in homologs]
            homologs = [h for h in updated if h is not None]
        homologs = [h for h in homologs if "Genome" in h]

        promoter_params = req.promoter_params.model_dump()
        homologs = [_operon_and_promoter(h, promoter_params) for h in homologs]

        return {
            "input": req.model_dump(),
            "protein_info": protein_info,
            "homologs": [_normalize_homolog(h) for h in homologs],
        }
    finally:
        os.chdir(prev_cwd)


async def run_pipeline(req: PredictRequest) -> AsyncIterator[str]:
    """Run the full pipeline, yielding SSE-formatted strings."""

    # 1. Smart-lookup pre-flight: BLAST input against precomputed family
    #    representatives. If we hit a known protein above-threshold, return its
    #    cached prediction with match metadata. `req.force` skips this.
    if not req.force:
        smart_hit = await asyncio.to_thread(smart_lookup_mod.lookup, req)
        if smart_hit is not None:
            payload = smart_lookup_mod.to_event_payload(smart_hit)
            yield _sse("smart_lookup_hit", **payload)
            return

    # 2. Parameter-hash cache check (per-query exact-match cache).
    key = cache.compute_key(req)
    if not req.force:
        hit = cache.read(key)
        if hit is not None:
            yield _sse("cached", cache_key=key, result=hit)
            return

    # Need to chdir for diamond to find its database (`../databases/...`).
    prev_cwd = os.getcwd()
    os.chdir(PROJECT_ROOT)
    try:
        # === Stage 1: BLAST ===
        yield _sse("stage_started", stage="blast", database=req.database)

        # For remote nr we'll be polling NCBI for many minutes; use a heartbeat
        # loop so the SSE connection stays alive and the user sees progress.
        # For local diamond the inner call returns in seconds and we never tick.
        blast_task = asyncio.create_task(asyncio.to_thread(_blast_stage, req))
        elapsed = 0.0
        heartbeat_interval = 15.0
        while not blast_task.done():
            try:
                await asyncio.wait_for(asyncio.shield(blast_task), timeout=heartbeat_interval)
            except asyncio.TimeoutError:
                elapsed += heartbeat_interval
                yield _sse(
                    "progress",
                    stage="blast",
                    database=req.database,
                    elapsed_seconds=int(elapsed),
                    message="BLAST in progress",
                )
        homologs, protein_info, err = await blast_task
        if err is not None:
            yield _sse("error", stage="blast", **err)
            return
        yield _sse(
            "stage_done",
            stage="blast",
            homologs=[_normalize_homolog(h) for h in homologs],
            protein_info=protein_info,
        )

        # === Stage 2: Genome coordinates ===
        yield _sse("stage_started", stage="coordinates")
        if req.get_coordinates_method == "batch":
            result = await asyncio.to_thread(_coords_batch, homologs)
            if result is None:
                yield _sse(
                    "error",
                    stage="coordinates",
                    kind="batch_failed",
                    message="Failed fetching genome coordinates in batch mode. Try individually.",
                )
                return
            homologs = [h for h in result if h is not None]
        else:
            updated = []
            for i, h in enumerate(homologs):
                yield _sse(
                    "progress",
                    stage="coordinates",
                    current=i,
                    total=len(homologs),
                    uniprot_id=h["Uniprot Id"],
                )
                updated.append(await asyncio.to_thread(_coords_one, h))
            homologs = [h for h in updated if h is not None]

        homologs = [h for h in homologs if "Genome" in h]
        yield _sse(
            "stage_done",
            stage="coordinates",
            homologs=[_normalize_homolog(h) for h in homologs],
        )

        # === Stage 3: Operons + promoters ===
        yield _sse("stage_started", stage="operons")
        promoter_params = req.promoter_params.model_dump()
        for i, h in enumerate(homologs):
            yield _sse(
                "progress",
                stage="operons",
                current=i,
                total=len(homologs),
                uniprot_id=h["Uniprot Id"],
            )
            await asyncio.to_thread(_operon_and_promoter, h, promoter_params)

        result_payload = {
            "input": req.model_dump(),
            "protein_info": protein_info,
            "homologs": [_normalize_homolog(h) for h in homologs],
        }
        saved = cache.write(key, result_payload)
        yield _sse("complete", cache_key=key, result=saved)
    finally:
        os.chdir(prev_cwd)
