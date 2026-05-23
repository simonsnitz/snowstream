"""Named algorithm bundles. Each version pins one function per algorithm
class. The benchmark harness in `benchmark/` iterates over these.

Add new versions here. Keep the keys in chronological order — they're
used to derive `v0 → v1 → v2 …` comparison axes in the report.
"""

from __future__ import annotations

from . import operon_fetch, promoter_fetch, operator_fetch


VERSIONS: dict[str, dict] = {
    "v0": {
        "operon_fetch":   operon_fetch.v0.fetch,
        "promoter_fetch": promoter_fetch.v0.fetch,
        "operator_fetch": operator_fetch.v0.fetch,
        "description": ("Original algorithms — legacy acc2operon + "
                         "single-candidate promoter fetch + "
                         "inverted-repeat operator finder. Reference point "
                         "for all future comparisons."),
    },
    "v1": {
        "operon_fetch":   operon_fetch.v0.fetch,
        "promoter_fetch": promoter_fetch.v1.fetch,
        "operator_fetch": operator_fetch.v0.fetch,
        "description": ("V0 + multi-candidate promoter enumeration. The "
                         "promoter_fetch returns the legacy primary "
                         "candidate plus any intra-operon same-direction "
                         "gaps that exceed min_internal_length (default "
                         "60 bp). The pipeline runner tries each candidate "
                         "as the query promoter and keeps the result with "
                         "the highest consensus_score. By construction "
                         "this can never regress V0."),
    },
}


def list_versions() -> list[str]:
    """Stable-ordered list of version names."""
    return list(VERSIONS.keys())


def get_version(name: str) -> dict:
    if name not in VERSIONS:
        available = ", ".join(VERSIONS)
        raise KeyError(f"unknown version {name!r}; available: {available}")
    return VERSIONS[name]
