"""promoter_fetch v0 — original single-candidate algorithm.

Port of src/fetch_promoter.py with two cleanups:
  * Removed the streamlit @cache_data decorator (irrelevant outside the
    legacy UI; the live pipeline doesn't need it).
  * Removed the global module-level streamlit import.

Logic is otherwise byte-identical to the legacy version, including the
known `regType=2` loop quirk where the gap test uses the regulator index
instead of the walking `index` (we intentionally preserve this so V0 is
faithful to the legacy behaviour we benchmarked against).

Returns a single intergenic-region string, or None if the algorithm
can't find a region that meets `params["min_length"]` / `params["max_length"]`.
"""

from __future__ import annotations

from typing import Optional

from .._shared import ncbi_efetch_nuccore


def fetch(operon_data: dict, params: dict) -> Optional[str]:
    # Fast path: snowstream's smart-lookup already cached the V0 promoter
    # for every member in the precomputed members DB. Reuse it when present
    # and length-valid — saves an NCBI eFetch round-trip per homolog.
    cached = operon_data.get("cached_promoter")
    if cached and isinstance(cached, str):
        if params["min_length"] <= len(cached) <= params["max_length"]:
            return cached

    operon = operon_data["operon"]
    reg_idx = operon_data["protein_index"]
    genome_id = operon_data["genome"]

    start_pos = None
    stop_pos = None

    if operon[reg_idx]["direction"] == "+":
        query_genes = list(reversed(operon[0:reg_idx]))
        index = reg_idx
        if len(query_genes) == 0:
            return None
        for i in query_genes:
            if i["direction"] == "-":
                start_pos = i["stop"]
                stop_pos = operon[index]["start"]
                break
            else:
                start = operon[reg_idx - 1]["stop"]
                stop = operon[reg_idx]["start"]
                test_length = int(stop) - int(start)
                if test_length > params["min_length"]:
                    start_pos = start
                    stop_pos = stop
                    break
                else:
                    if index == 1:
                        return None
                    index -= 1

    elif operon[reg_idx]["direction"] == "-":
        query_genes = operon[reg_idx + 1:]
        index = reg_idx
        if len(query_genes) == 0:
            return None
        for i in query_genes:
            if i["direction"] == "+":
                stop_pos = i["start"]
                start_pos = operon[index]["stop"]
                break
            else:
                start = operon[reg_idx]["stop"]
                stop = operon[reg_idx + 1]["start"]
                test_length = int(stop) - int(start)
                if test_length > 100:  # hardcoded in legacy code
                    start_pos = start
                    stop_pos = stop
                    break
                else:
                    if index == len(operon) - 2:
                        return None
                    index += 1

    if start_pos is None or stop_pos is None:
        return None
    intergenic = ncbi_efetch_nuccore(
        genome_id, int(start_pos), int(stop_pos), strand=1, rettype="fasta")
    if not intergenic:
        return None
    if (len(intergenic) > params["max_length"]
            or len(intergenic) < params["min_length"]):
        return None
    return intergenic
