"""Named algorithm bundles. Each version pins one function per algorithm
class (plus optional pre-bound params via a wrapper). The benchmark
harness in `benchmark/` iterates over these.

Version naming:
  * Major versions (v0, v1, v2, v3) encode which algorithm classes change.
  * Decimal sub-versions (v2.1, v2.2, ...) are ablations within an
    algorithm change — e.g. the V2.x matrix factors operator_fetch v1
    into its two component changes (widened candidate pool vs new
    selection metric) and tests three strengths of the AT/GC penalty.

Add new versions here. Keep the keys in roughly chronological order.
"""

from __future__ import annotations

from typing import Callable

from . import operon_fetch, promoter_fetch, operator_fetch


def _bind_operator_v1(candidate_strategy: str, rerank_scorer: str) -> Callable:
    """Pre-bind operator_fetch.v1's tunable params (candidate_strategy,
    rerank_scorer) so each V2.x bundle is a self-contained function. The
    rest of the operator_fetch params come from DEFAULT_OPERATOR_PARAMS
    in benchmark/run.py — we only override the two ablation knobs.
    """
    base_fn = operator_fetch.v1.fetch

    def fetch(homologs, params):
        return base_fn(homologs, {
            **params,
            "candidate_strategy": candidate_strategy,
            "rerank_scorer": rerank_scorer,
        })
    fetch.__name__ = f"operator_v1[{candidate_strategy}/{rerank_scorer}]"
    return fetch


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
                         "60 bp). Operator_fetch is unchanged."),
    },

    # --- V2.x ablation matrix --------------------------------------------
    # Factor operator_fetch v1 into its two independent changes:
    #   - candidate_strategy: "legacy" (max-IR-tied) vs "widened" (top-25)
    #   - rerank_scorer:      "consensus_score" (V0 metric) vs gc_weak /
    #                          gc_medium_weak / gc_strong (new AT-based
    #                          gradient)
    # The 7 combinations span every meaningful permutation.

    "v2.1": {
        "operon_fetch":   operon_fetch.v0.fetch,
        "promoter_fetch": promoter_fetch.v0.fetch,
        "operator_fetch": _bind_operator_v1("widened", "consensus_score"),
        "description": ("Widen candidate pool only (top-25 unique by IR "
                         "score); keep V0's consensus_score selector. "
                         "Measures the effect of widening alone."),
    },
    "v2.2": {
        "operon_fetch":   operon_fetch.v0.fetch,
        "promoter_fetch": promoter_fetch.v0.fetch,
        "operator_fetch": _bind_operator_v1("legacy", "gc_weak"),
        "description": ("Legacy narrow candidate pool + AT-based selector "
                         "with WEAK GC penalty (0.2·AT + 0.5·length). "
                         "Measures the effect of the new selector alone "
                         "at the lightest AT preference."),
    },
    "v2.3": {
        "operon_fetch":   operon_fetch.v0.fetch,
        "promoter_fetch": promoter_fetch.v0.fetch,
        "operator_fetch": _bind_operator_v1("legacy", "gc_medium_weak"),
        "description": ("Legacy narrow candidate pool + AT-based selector "
                         "with MEDIUM-WEAK GC penalty (0.5·AT + 0.3·length)."),
    },
    "v2.4": {
        "operon_fetch":   operon_fetch.v0.fetch,
        "promoter_fetch": promoter_fetch.v0.fetch,
        "operator_fetch": _bind_operator_v1("legacy", "gc_strong"),
        "description": ("Legacy narrow candidate pool + AT-based selector "
                         "with STRONG GC penalty (1.0·AT + 0.3·length, "
                         "the original `combined_v2`). Tests whether the "
                         "new selector helps even without widening."),
    },
    "v2.5": {
        "operon_fetch":   operon_fetch.v0.fetch,
        "promoter_fetch": promoter_fetch.v0.fetch,
        "operator_fetch": _bind_operator_v1("widened", "gc_weak"),
        "description": ("Widened pool + AT-based selector with WEAK GC "
                         "penalty. Combined ablation, lightest AT bias."),
    },
    "v2.6": {
        "operon_fetch":   operon_fetch.v0.fetch,
        "promoter_fetch": promoter_fetch.v0.fetch,
        "operator_fetch": _bind_operator_v1("widened", "gc_medium_weak"),
        "description": ("Widened pool + AT-based selector with MEDIUM-WEAK "
                         "GC penalty. Likely the best blend — strong "
                         "enough to avoid spurious GC palindromes, weak "
                         "enough to not over-bias against GC-rich-genome "
                         "operators."),
    },
    "v2.7": {
        "operon_fetch":   operon_fetch.v0.fetch,
        "promoter_fetch": promoter_fetch.v0.fetch,
        "operator_fetch": _bind_operator_v1("widened", "gc_strong"),
        "description": ("Widened pool + AT-based selector with STRONG GC "
                         "penalty. Equivalent to the original V2 from "
                         "PR #15 (closed); reproduced here as the strong "
                         "end of the V2.x gradient. May regress on "
                         "GC-rich hosts (Streptomyces, Mycobacterium)."),
    },

    # --- V2.8-V2.11: intermediate GC-penalty scorers ----------------------
    # Densify the GC-penalty gradient on the widened-pool side, where the
    # action is. V2.5 (weak) → V2.6 (medium-weak) → V2.7 (strong) showed
    # V2.6 was best; these four fill in points between the coarse trio.

    "v2.8": {
        "operon_fetch":   operon_fetch.v0.fetch,
        "promoter_fetch": promoter_fetch.v0.fetch,
        "operator_fetch": _bind_operator_v1("widened", "gc_weak_plus"),
        "description": ("Widened pool + 0.3·AT + 0.4·length. Between V2.5 "
                         "(weak) and V2.6 (medium-weak), closer to weak."),
    },
    "v2.9": {
        "operon_fetch":   operon_fetch.v0.fetch,
        "promoter_fetch": promoter_fetch.v0.fetch,
        "operator_fetch": _bind_operator_v1("widened", "gc_weak_strong"),
        "description": ("Widened pool + 0.4·AT + 0.35·length. Between V2.5 "
                         "and V2.6, closer to medium-weak."),
    },
    "v2.10": {
        "operon_fetch":   operon_fetch.v0.fetch,
        "promoter_fetch": promoter_fetch.v0.fetch,
        "operator_fetch": _bind_operator_v1("widened", "gc_medium"),
        "description": ("Widened pool + 0.65·AT + 0.3·length. Between V2.6 "
                         "and V2.7, closer to medium-weak."),
    },
    "v2.11": {
        "operon_fetch":   operon_fetch.v0.fetch,
        "promoter_fetch": promoter_fetch.v0.fetch,
        "operator_fetch": _bind_operator_v1("widened", "gc_medium_strong"),
        "description": ("Widened pool + 0.8·AT + 0.3·length. Between V2.6 "
                         "and V2.7, closer to strong."),
    },
}

# Combined V3.x bundles (V1 multi-candidate promoter + V2.x operator) were
# benchmarked and removed. V3.6 was a strict-but-tiny improvement over V2.6
# (3 wins / 0 losses / 132 within ±1pp; mean +0.1pp); the per-query NCBI
# overhead for alternative promoter candidates wasn't worth the marginal
# accuracy gain. V2.6 chosen as the production-ready peak. See the
# operator-fetch-ablation PR thread for the supporting benchmark.


def list_versions() -> list[str]:
    """Stable-ordered list of version names."""
    return list(VERSIONS.keys())


def get_version(name: str) -> dict:
    if name not in VERSIONS:
        available = ", ".join(VERSIONS)
        raise KeyError(f"unknown version {name!r}; available: {available}")
    return VERSIONS[name]
