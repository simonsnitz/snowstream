"""operator_fetch v1 — widened candidate enumeration + re-rank by AT/length.

Drop-in replacement for v0 with two changes that share a single
motivation: the legacy inverted-repeat score actively misleads candidate
selection (AUC 0.188 against the held-out positive/negative set in
algorithms/benchmark/LATEST_SCORERS.md). Real TetR operators are
AT-rich and operator-sized; spurious palindromes that fool the v0 picker
are GC-rich and often the wrong length.

What changes vs v0:

1. **Candidate enumeration is widened.** The legacy `findBestPalindrome`
   only returns palindromes tied at the *global maximum* inverted-repeat
   score, which can collapse the candidate pool to a single (and often
   wrong) sequence. v1 enumerates *all* imperfect palindromes across the
   configured size range, dedupes by sequence, and keeps the top-K by IR
   score (default K = 25).

2. **Final winner is picked by `combined_v2` (AT% + length-Gaussian),
   not consensus_score.** For each candidate the existing align-to-
   homologs / consensus-motif machinery still runs (so the returned
   motif, num_seqs, frequency_matrix etc. all have v0-identical
   semantics). But the *selection* between candidates uses the AT-rich,
   operator-sized score from the scorer-classifier benchmark, which had
   AUC 0.916 on the same dataset.

Non-inverted-repeats search modes ("Align an input sequence", "Scan
entire promoter region") fall through to v0 unchanged.

The returned dict has the same shape as v0, plus one new field:
  * `rerank_score` — the combined_v2 score of the winning candidate.

Tunable params (all optional, sensible defaults):
  * `top_k_palindromes` (default 25) — candidate pool size
  * `rerank_scorer`     (default "combined_v2") — name from
                         algorithms.benchmark.operator_scorers.SCORERS
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Callable, Optional

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Reuse the legacy primitives — only the selection layer changes.
from src.fetch_operator import (  # noqa: E402
    complement,
    findImperfectPalindromes,
    findOperatorInIntergenic,
    getConsensus,
    get_consensus_score,
    generate_frequency_matrix,
)

# Default rerank scorer comes from the benchmark registry so the choice
# is centrally documented and easy to swap.
from algorithms.benchmark.operator_scorers import SCORERS  # noqa: E402

from . import v0 as v0_mod


_DEFAULT_TOP_K = 25
_DEFAULT_RERANK = "combined_v2"


# --- Helpers -------------------------------------------------------------

def _scorer_by_name(name: str) -> Callable[[str], float]:
    for n, fn, _ in SCORERS:
        if n == name:
            return fn
    raise KeyError(f"unknown rerank scorer {name!r}; choose from "
                    f"{[n for n, _, _ in SCORERS]}")


def _enumerate_palindromes(intergenic: str, shortest: int, longest: int,
                            win_score: int, loss_score: int,
                            spacer_penalty: dict) -> list[dict]:
    """Return *every* imperfect palindrome across the size range as
    `{seq, score}` dicts — no per-size or global max filter applied.

    Equivalent to the inner loop of `findBestPalindrome` but without the
    `if score == max_score` filter that collapses the candidate pool.
    """
    out: list[dict] = []
    intergenic = intergenic.upper()
    for size in range(shortest, longest):
        ops = findImperfectPalindromes(
            intergenic, size, win_score, loss_score, spacer_penalty)
        if ops:
            out.extend(ops)
    return out


def _top_k_unique_palindromes(intergenic: str, shortest: int, longest: int,
                               win_score: int, loss_score: int,
                               spacer_penalty: dict, k: int) -> list[dict]:
    """Take the union of `findImperfectPalindromes` outputs across all
    sizes, dedupe by uppercase seq, and return the top-K by score.

    `findImperfectPalindromes` itself filters to max-score-per-size, so
    the input pool is already curated — we just stop applying the
    per-size restriction and the global-max restriction that
    `findBestPalindrome` adds on top.
    """
    pool = _enumerate_palindromes(intergenic, shortest, longest,
                                    win_score, loss_score, spacer_penalty)
    seen: set[str] = set()
    unique: list[dict] = []
    for p in sorted(pool, key=lambda x: x["score"], reverse=True):
        key = (p.get("seq") or "").upper()
        if not key or key == "NONE" or key in seen:
            continue
        seen.add(key)
        unique.append(p)
        if len(unique) >= k:
            break
    return unique


def _consensus_seq(consensus: dict) -> str:
    md = consensus.get("motif_data") or []
    return "".join(x["base"] for x in md if isinstance(x, dict))


# --- Public interface ----------------------------------------------------

def fetch(homolog_metadata: list[dict], params: dict) -> dict:
    """Same I/O contract as v0. See module docstring for what changes."""

    # Non-inverted-repeats modes are unchanged from v0.
    if params.get("search_method") != "Look for inverted repeats":
        return v0_mod.fetch(homolog_metadata, params)

    ext_length = params["extension_length"]
    acc = homolog_metadata[0]["Uniprot Id"]
    regulated_seqs = [h["promoter"] for h in homolog_metadata]
    reference_candidates = [s for s in regulated_seqs if s]
    if not reference_candidates:
        # No usable promoter — defer to v0 which handles the failure path.
        return v0_mod.fetch(homolog_metadata, params)
    reference_seq = reference_candidates[0]

    top_k = int(params.get("top_k_palindromes", _DEFAULT_TOP_K))
    rerank_name = params.get("rerank_scorer", _DEFAULT_RERANK)
    rerank_fn = _scorer_by_name(rerank_name)

    candidates = _top_k_unique_palindromes(
        intergenic=reference_seq,
        shortest=params["min_operator_length"],
        longest=params["max_operator_length"],
        win_score=params["win_score"],
        loss_score=params["loss_score"],
        spacer_penalty=params["spacer_penalty"],
        k=top_k,
    )

    operator_data = {
        "Uniprot Id": str(acc),
        "aligned_seq": "None",
        "num_seqs": "None",
        "consensus_score": 0,
        "rerank_score": -math.inf,
        "rerank_scorer": rerank_name,
        "n_candidates_evaluated": 0,
        "motif": "None",
        "aligned_seqs": "None",
        "intergenic": reference_seq,
    }

    for candidate in candidates:
        metrics: list[dict] = []
        for h in homolog_metadata:
            promoter = h.get("promoter")
            if not promoter:
                continue
            op = findOperatorInIntergenic(promoter, candidate["seq"], params)
            if op is None:
                continue
            metrics.append({
                "Uniprot Id": h["Uniprot Id"],
                "Predicted operator": op["operator"],
                "Align score": op["score"],
            })
        if not metrics:
            continue

        operator_data["n_candidates_evaluated"] += 1

        consensus = getConsensus(metrics)
        consensus_score = get_consensus_score(
            candidate["seq"], consensus, ext_length)

        # Render the consensus to a sequence and score it with the new metric.
        motif_str = _consensus_seq(consensus)
        rerank_score = rerank_fn(motif_str) if motif_str else float("-inf")

        # Also extract the per-query operator (same as v0)
        op_seq = findOperatorInIntergenic(
            reference_seq, candidate["seq"], params)
        native_operator = (op_seq["operator"] if op_seq is not None
                            else candidate["seq"])

        # Selection metric is the rerank_score — NOT consensus_score.
        if rerank_score > operator_data["rerank_score"]:
            operator_data["rerank_score"] = rerank_score
            operator_data["consensus_score"] = consensus_score
            operator_data["native_operator"] = native_operator
            operator_data["consensus_seq"] = motif_str
            operator_data["num_seqs"] = consensus["num_seqs"]
            operator_data["motif"] = consensus["motif_data"]
            operator_data["aligned_seqs"] = metrics
            operator_data["frequency_matrix"] = generate_frequency_matrix(metrics)

    return operator_data
