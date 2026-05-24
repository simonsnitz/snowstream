"""Version-aware end-to-end pipeline runner.

Given a version bundle (from versions.py), a list of homologs (each with an
operon already computed), and the standard parameter dicts, this runs the
operator-prediction pipeline:

    for each homolog: promoter_fetch(operon, params) -> str or list[str]
    if all promoters are strings:
        operator_fetch(homologs_with_promoters, params)
    else (multi-candidate mode):
        for each query-promoter candidate:
            operator_fetch(...) with that candidate as the query promoter
            and the V0 primary for each other homolog
        keep the scenario with the highest consensus_score

We deliberately don't re-run operon_fetch here — the benchmark harness
short-circuits by reusing snowstream's precomputed members DB (the operons
are already correct and cached). When a new operon_fetch version needs
testing, swap how `homologs_with_operons` are produced upstream.
"""

from __future__ import annotations

from typing import Optional

from .promoter_fetch import v0 as promoter_v0  # legacy primary fallback


def _normalise_promoter(out) -> list[str]:
    """promoter_fetch may return None, a single string, or a list of strings.
    Normalise to a list of (zero or more) non-empty strings."""
    if out is None:
        return []
    if isinstance(out, str):
        return [out] if out else []
    if isinstance(out, list):
        return [s for s in out if s]
    raise TypeError(f"promoter_fetch returned unexpected type: {type(out)}")


def _operon_data(homolog: dict) -> dict:
    """Build the operon_data dict promoter_fetch expects. Passes the
    snowstream-cached promoter (if any) through as `cached_promoter` so
    V0/V1 can short-circuit the NCBI eFetch when the precomputed value
    is already available — saves ~100 HTTP round-trips per protein."""
    return {
        "operon": homolog.get("operon") or [],
        "protein_index": homolog.get("protein_index"),
        "genome": homolog.get("genome"),
        "cached_promoter": homolog.get("promoter"),
    }


def _promoter_for(homolog: dict, promoter_fetch_fn, params: dict) -> Optional[str]:
    """Resolve a single promoter for a non-query homolog. If the version's
    promoter_fetch returns multiple candidates, pick the V0 primary (the
    most legacy-faithful choice) so we don't combinatorially explode the
    search space."""
    operon_data = _operon_data(homolog)
    if not operon_data["operon"] or operon_data["protein_index"] is None:
        # Even without an operon, if snowstream has a cached promoter
        # we can use it (the smart-lookup record may be valid even when
        # this benchmark-side reconstruction lacks fields).
        cached = operon_data.get("cached_promoter")
        return cached if cached else None
    primaries = _normalise_promoter(promoter_fetch_fn(operon_data, params))
    if primaries:
        return primaries[0]
    legacy = promoter_v0.fetch(operon_data, params)
    return legacy if isinstance(legacy, str) and legacy else None


def run(version: dict, homologs: list[dict], params: dict) -> dict:
    """Run the version's pipeline on the given homologs. The first homolog
    is treated as the query — its promoter is the reference for the
    palindrome search.

    Returns:
        {
          "version_description": str,
          "n_homologs": int,
          "n_homologs_with_promoter": int,
          "query_promoter": str | None,    # the promoter that won (best score)
          "all_query_candidates": list[str],
          "operator_result": dict,         # raw operator_fetch output
        }
    """
    promoter_fetch_fn = version["promoter_fetch"]
    operator_fetch_fn = version["operator_fetch"]

    if not homologs:
        return {
            "version_description": version.get("description"),
            "n_homologs": 0, "n_homologs_with_promoter": 0,
            "query_promoter": None, "all_query_candidates": [],
            "operator_result": None,
        }

    query = homologs[0]
    others = homologs[1:]

    # 1. Resolve query's candidate promoters (uses cached_promoter shortcut)
    query_operon = _operon_data(query)
    query_candidates = _normalise_promoter(promoter_fetch_fn(query_operon, params))

    # 2. Resolve a single promoter for each other homolog (no combinatorial
    #    explosion: use the version's first candidate, falling back to V0)
    other_with_promoters: list[dict] = []
    for h in others:
        p = _promoter_for(h, promoter_fetch_fn, params)
        if p:
            other_with_promoters.append({
                "Uniprot Id": h.get("uniprot_id") or h.get("Uniprot Id"),
                "promoter": p,
            })

    # 3. For each query candidate, run operator_fetch and keep the best
    best_result = None
    best_score = -1.0
    best_query_promoter = None
    for candidate in query_candidates:
        full_set = [{
            "Uniprot Id": query.get("uniprot_id") or query.get("Uniprot Id"),
            "promoter": candidate,
        }] + other_with_promoters
        try:
            result = operator_fetch_fn(full_set, params)
        except Exception:
            continue
        score = float(result.get("consensus_score") or 0)
        if score > best_score:
            best_score = score
            best_result = result
            best_query_promoter = candidate

    n_with_promoter = (1 if best_query_promoter else 0) + len(other_with_promoters)
    return {
        "version_description": version.get("description"),
        "n_homologs": len(homologs),
        "n_homologs_with_promoter": n_with_promoter,
        "query_promoter": best_query_promoter,
        "all_query_candidates": query_candidates,
        "all_homolog_promoters": [h["promoter"] for h in other_with_promoters],
        "operator_result": best_result,
    }
