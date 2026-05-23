"""Per-protein metrics for the benchmark harness.

A metric is a function `(dataset_entry, pipeline_output, ctx) -> value` —
where `ctx` is a small helper dict carrying things like the list of
candidate promoters that the pipeline considered.

Metrics fall into two natural shapes:
  * **counts** — integers (e.g. `n_homologs`). Reported as raw values
    plus distribution histograms.
  * **identities** — percentages 0-100. Reported as raw values plus
    threshold-cross counts at [60, 70, 80, 90].

The report generator picks up the metric type from `kind` so adding a new
one just requires appending to the METRICS list.

Identity thresholds are centralised here too — the report uses the same
buckets across versions for consistency.
"""

from __future__ import annotations

from typing import Callable

from .._shared import best_window_match


IDENTITY_THRESHOLDS = (60, 70, 80, 90)


# --- individual metric functions -----------------------------------------

def _all_known_operators(entry: dict) -> list[str]:
    return [o["sequence"] for o in entry.get("operators", []) if o.get("sequence")]


def n_homologs(entry, out, ctx) -> int:
    return int(out.get("n_homologs") or 0)


def n_homologs_with_promoter(entry, out, ctx) -> int:
    return int(out.get("n_homologs_with_promoter") or 0)


def known_operator_in_query_promoter(entry, out, ctx) -> float:
    """Best identity (across all known operators × both orientations) of any
    known operator against the query's chosen promoter."""
    promoter = out.get("query_promoter") or ""
    if not promoter:
        return 0.0
    best = 0.0
    for op in _all_known_operators(entry):
        m = best_window_match(op, promoter)
        if m["identity_pct"] > best:
            best = m["identity_pct"]
    return best


def known_operator_in_any_homolog_promoter(entry, out, ctx) -> float:
    """Best identity across *all* homolog promoters (query + others)."""
    promoters = []
    qp = out.get("query_promoter")
    if qp:
        promoters.append(qp)
    promoters.extend(out.get("all_homolog_promoters") or [])
    if not promoters:
        return 0.0
    best = 0.0
    for op in _all_known_operators(entry):
        for p in promoters:
            m = best_window_match(op, p)
            if m["identity_pct"] > best:
                best = m["identity_pct"]
    return best


def known_operator_in_predicted_motif(entry, out, ctx) -> float:
    """Best identity of any known operator against the predicted consensus
    motif from operator_fetch (the actual algorithm output)."""
    result = out.get("operator_result") or {}
    motif = result.get("motif")
    if not motif:
        return 0.0
    if isinstance(motif, list):
        motif_str = "".join(x["base"] for x in motif if isinstance(x, dict))
    elif isinstance(motif, str):
        motif_str = motif
    else:
        return 0.0
    if not motif_str:
        return 0.0
    best = 0.0
    for op in _all_known_operators(entry):
        m = best_window_match(op, motif_str)
        if m["identity_pct"] > best:
            best = m["identity_pct"]
    return best


# --- metric registry -----------------------------------------------------

# Each entry: (name, function, kind, label)
# `kind` ∈ {"count", "identity"} — drives how the report aggregates the
# metric. Add new metrics by appending here; the report picks them up
# automatically.
METRICS: list[tuple[str, Callable, str, str]] = [
    ("n_homologs",
     n_homologs, "count",
     "# homologs"),
    ("n_homologs_with_promoter",
     n_homologs_with_promoter, "count",
     "# homologs with promoter"),
    ("known_operator_in_query_promoter",
     known_operator_in_query_promoter, "identity",
     "Known operator in query promoter"),
    ("known_operator_in_any_homolog_promoter",
     known_operator_in_any_homolog_promoter, "identity",
     "Known operator in any homolog promoter"),
    ("known_operator_in_predicted_motif",
     known_operator_in_predicted_motif, "identity",
     "Known operator in predicted motif"),
]


def compute_all(entry: dict, pipeline_out: dict, ctx: dict | None = None) -> dict:
    ctx = ctx or {}
    return {name: fn(entry, pipeline_out, ctx) for name, fn, _, _ in METRICS}
