"""promoter_fetch v1 — multi-candidate intergenic enumeration.

Same input contract as v0, but returns a *list* of candidate intergenic
sequences rather than a single string. The list always includes:
  1. The V0 primary (whatever legacy fetch returned, if anything) — so
     downstream consumers that pick the best candidate by score can never
     regress relative to V0.
  2. Any intra-operon same-direction gaps that exceed `min_internal_length`
     (default 60 bp; configurable via params["min_internal_length"]).
  3. The gap at the operon boundary (the gap to a divergent gene if one
     is present upstream, otherwise the last same-direction gap walked).

Each candidate is fetched from NCBI separately and length-filtered against
`params["min_length"]` / `params["max_length"]`. Candidates that fail the
filter (too short or too long) are dropped silently.

The downstream pipeline runner picks the candidate that produces the
highest-scoring operator from operator_fetch.
"""

from __future__ import annotations

from typing import Optional

from .._shared import ncbi_efetch_nuccore
from . import v0 as v0_mod


def _enumerate_candidate_coords(operon: list[dict], reg_idx: int,
                                 min_internal: int = 60) -> list[dict]:
    """Walk upstream from the regulator and emit every plausible intergenic
    gap. Returns a list of {start, stop, length, kind, before_idx, after_idx}.

    `kind` is one of:
      * "primary" — first same-direction gap > min_internal encountered
                    walking upstream (matches the legacy `regType=2` slot)
      * "internal" — any subsequent same-direction gap > min_internal
                    encountered before hitting the operon boundary
      * "boundary" — the gap to the divergent gene where the operon ends,
                    OR the last same-direction gap if we walk off the end
                    of the operon without hitting a divergent gene
    Overlapping genes (stop <= start) are skipped.
    """
    out: list[dict] = []
    if not operon or reg_idx is None:
        return out
    n = len(operon)
    if reg_idx < 0 or reg_idx >= n:
        return out

    def _add(start: int, stop: int, kind: str, before: int, after: int) -> bool:
        if stop <= start:
            return False
        out.append({"start": int(start), "stop": int(stop),
                    "length": int(stop - start), "kind": kind,
                    "before_idx": before, "after_idx": after})
        return True

    direction = operon[reg_idx].get("direction")
    primary_recorded = False

    if direction == "+":
        for j in range(reg_idx - 1, -1, -1):
            if operon[j].get("direction") == "-":
                try:
                    _add(int(operon[j]["stop"]),
                         int(operon[j + 1]["start"]),
                         "boundary", j, j + 1)
                except Exception:
                    pass
                break
            if j > 0:
                try:
                    start = int(operon[j - 1]["stop"])
                    stop = int(operon[j]["start"])
                    length = stop - start
                except Exception:
                    continue
                if length > min_internal:
                    kind = "primary" if not primary_recorded else "internal"
                    if _add(start, stop, kind, j - 1, j):
                        primary_recorded = True

    elif direction == "-":
        for j in range(reg_idx + 1, n):
            if operon[j].get("direction") == "+":
                try:
                    _add(int(operon[j - 1]["stop"]),
                         int(operon[j]["start"]),
                         "boundary", j - 1, j)
                except Exception:
                    pass
                break
            if j < n - 1:
                try:
                    start = int(operon[j]["stop"])
                    stop = int(operon[j + 1]["start"])
                    length = stop - start
                except Exception:
                    continue
                if length > min_internal:
                    kind = "primary" if not primary_recorded else "internal"
                    if _add(start, stop, kind, j, j + 1):
                        primary_recorded = True

    return out


def fetch(operon_data: dict, params: dict) -> list[str]:
    """Return a list of candidate intergenic sequences. May be empty."""
    operon = operon_data["operon"]
    reg_idx = operon_data["protein_index"]
    genome_id = operon_data["genome"]

    min_len = params.get("min_length", 80)
    max_len = params.get("max_length", 800)
    min_internal = params.get("min_internal_length", 60)

    out: list[str] = []
    seen: set[str] = set()

    def _maybe_add(seq: Optional[str]) -> None:
        if not seq:
            return
        if len(seq) < min_len or len(seq) > max_len:
            return
        key = seq.upper()
        if key in seen:
            return
        seen.add(key)
        out.append(seq)

    # 1. Always include the V0 primary first (preserves legacy behaviour)
    _maybe_add(v0_mod.fetch(operon_data, params))

    # 2. Walk and fetch each enumerated candidate
    for c in _enumerate_candidate_coords(operon, reg_idx, min_internal):
        seq = ncbi_efetch_nuccore(genome_id, c["start"], c["stop"],
                                   strand=1, rettype="fasta")
        _maybe_add(seq)

    return out
