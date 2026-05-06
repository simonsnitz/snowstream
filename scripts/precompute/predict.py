"""Stage: run the full Snowprint pipeline for each cluster representative,
appending one record per representative to families/<key>/predictions.jsonl.

Resumable & idempotent: at startup we scan the existing JSONL once to record
which UniProt IDs are already present, and skip them. Each new record is
flushed to disk immediately so a Ctrl-C doesn't lose progress beyond the
in-flight entry.

The pipeline is invoked via `backend.pipeline.run_pipeline_core`, which is the
shared synchronous entrypoint. The caller is expected to set
`SNOWPRINT_DIAMOND_DB` to the local NCBI nr Diamond DB before running this
stage; otherwise the default small bHTH DB is used (fine for testing).
"""

from __future__ import annotations

import json
import logging
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional

# Import the shared pipeline core from backend/.
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from backend.pipeline import run_pipeline_core  # noqa: E402
from backend.schemas import (  # noqa: E402
    BlastParams,
    PredictRequest,
    PromoterParams,
)

log = logging.getLogger(__name__)


def _existing_uniprot_ids(jsonl_path: Path) -> set[str]:
    if not jsonl_path.exists():
        return set()
    out: set[str] = set()
    with jsonl_path.open() as fh:
        for line in fh:
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            uid = rec.get("uniprot_id")
            if uid:
                out.add(uid)
    return out


def _build_request(uniprot_id: str, sequence: str) -> PredictRequest:
    """Construct a PredictRequest from a representative's sequence.

    We feed the protein sequence directly (rather than the UniProt ID) so the
    BLAST input is identical regardless of whether the representative is in
    UniProt's flat files at the time of precompute. Default BLAST + promoter
    params; batch coordinate fetch.
    """
    return PredictRequest(
        input_method="Protein sequence",
        input_value=sequence,
        blast_params=BlastParams(),
        promoter_params=PromoterParams(),
        get_coordinates_method="batch",
        force=False,
    )


def predict_all(
    representatives: list[dict],
    sequences: dict[str, str],
    jsonl_path: Path,
    max_entries: Optional[int] = None,
    on_progress: Optional[Callable[[int, int, str, str], None]] = None,
) -> int:
    """Run the pipeline for every representative not already in the JSONL.

    `representatives` is the list produced by the rep-selection stage; each
    entry's `uniprot_id` and `evidence` are preserved verbatim in the
    resulting JSONL record.

    Returns the number of *new* records appended (not the total).
    """
    jsonl_path.parent.mkdir(parents=True, exist_ok=True)
    already = _existing_uniprot_ids(jsonl_path)
    log.info("predictions.jsonl already has %d entries", len(already))

    todo = [r for r in representatives if r["uniprot_id"] not in already]
    if max_entries is not None:
        todo = todo[:max_entries]
    log.info("will run pipeline for %d new representatives", len(todo))

    appended = 0
    with jsonl_path.open("a") as fh:
        for i, rep in enumerate(todo):
            uid = rep["uniprot_id"]
            seq = sequences.get(uid)
            if not seq:
                log.warning("missing sequence for %s; skipping", uid)
                continue
            if on_progress:
                on_progress(i, len(todo), uid, "starting")

            t0 = time.time()
            try:
                result = run_pipeline_core(_build_request(uid, seq))
            except Exception as exc:
                log.exception("pipeline crashed for %s: %s", uid, exc)
                if on_progress:
                    on_progress(i, len(todo), uid, f"crashed: {exc}")
                continue

            if "error" in result:
                log.warning("pipeline reported error for %s: %s", uid, result["error"])
                if on_progress:
                    on_progress(i, len(todo), uid, f"error: {result['error'].get('kind')}")
                # Still record the failure so we don't retry endlessly.
                record = {
                    "uniprot_id": uid,
                    "centroid_uniprot_id": rep.get("centroid_uniprot_id", uid),
                    "cluster": rep.get("cluster", {}),
                    "evidence": rep.get("evidence", {}),
                    "computed_at": datetime.now(timezone.utc).isoformat(),
                    "error": result["error"],
                }
            else:
                record = {
                    "uniprot_id": uid,
                    "centroid_uniprot_id": rep.get("centroid_uniprot_id", uid),
                    "cluster": rep.get("cluster", {}),
                    "evidence": rep.get("evidence", {}),
                    "computed_at": datetime.now(timezone.utc).isoformat(),
                    "input": result["input"],
                    "protein_info": result.get("protein_info"),
                    "homologs": result["homologs"],
                }

            fh.write(json.dumps(record) + "\n")
            fh.flush()
            appended += 1
            elapsed = time.time() - t0
            if on_progress:
                on_progress(i, len(todo), uid, f"done in {elapsed:.1f}s")
    log.info("appended %d records to %s", appended, jsonl_path)
    return appended
