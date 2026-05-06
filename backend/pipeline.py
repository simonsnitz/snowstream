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

from src.blast import blast  # noqa: E402
from src.get_genome_coordinates import (  # noqa: E402
    get_genome_coordinates,
    get_genome_coordinates_batch,
)
from src.accID2operon import acc2operon  # noqa: E402
from src.fetch_promoter import fetch_promoter  # noqa: E402
from src.troubleshoot import troubleshoot  # noqa: E402

from . import cache
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


async def run_pipeline(req: PredictRequest) -> AsyncIterator[str]:
    """Run the full pipeline, yielding SSE-formatted strings."""

    # Cache check
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
        yield _sse("stage_started", stage="blast")
        blast_df = await asyncio.to_thread(
            blast,
            req.input_value,
            req.input_method,
            req.blast_params.model_dump(),
            500,
        )

        if blast_df.empty:
            mode = None
            if req.input_method != "Protein sequence":
                mode, _ = await asyncio.to_thread(
                    troubleshoot, req.input_method, req.input_value
                )
            yield _sse(
                "error",
                stage="blast",
                kind=mode or "no_results",
                message=(
                    "Protein is not from Bacteria. Snowprint only works for bacterial proteins."
                    if mode == "not bacteria"
                    else "BLAST returned no results."
                ),
            )
            return

        homologs = (
            _filter_redundant(blast_df)
            if req.blast_params.filter_redundant
            else _all_homologs(blast_df)
        )
        homologs = homologs[: req.blast_params.max_homologs]

        protein_info = await asyncio.to_thread(
            _uniprot_protein_info, blast_df.iloc[0]["Uniprot Id"]
        )

        yield _sse(
            "stage_done",
            stage="blast",
            homologs=[_normalize_homolog(h) for h in homologs],
            protein_info=protein_info,
        )

        # === Stage 2: Genome coordinates ===
        yield _sse("stage_started", stage="coordinates")

        if req.get_coordinates_method == "batch":
            result = await asyncio.to_thread(get_genome_coordinates_batch, homologs)
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
                updated.append(await asyncio.to_thread(get_genome_coordinates, h))
            homologs = [h for h in updated if h is not None]

        homologs = [h for h in homologs if "Genome" in h]
        yield _sse(
            "stage_done",
            stage="coordinates",
            homologs=[_normalize_homolog(h) for h in homologs],
        )

        # === Stage 3: Operons + promoters ===
        yield _sse("stage_started", stage="operons")
        for i, h in enumerate(homologs):
            yield _sse(
                "progress",
                stage="operons",
                current=i,
                total=len(homologs),
                uniprot_id=h["Uniprot Id"],
            )
            h["operon"] = await asyncio.to_thread(acc2operon, h)
            try:
                h["promoter"] = await asyncio.to_thread(
                    fetch_promoter, h["operon"], req.promoter_params.model_dump()
                )
            except Exception:
                h["promoter"] = None

        result_payload = {
            "input": req.model_dump(),
            "protein_info": protein_info,
            "homologs": [_normalize_homolog(h) for h in homologs],
        }
        saved = cache.write(key, result_payload)
        yield _sse("complete", cache_key=key, result=saved)
    finally:
        os.chdir(prev_cwd)
