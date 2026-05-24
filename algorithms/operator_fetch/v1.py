"""operator_fetch v1 — parameterised candidate enumeration + selection.

Two independent knobs:

  * **candidate_strategy** — "legacy" (use `findBestPalindrome`, V0
    behaviour: palindromes tied at the global max IR score, often just
    1-5 candidates) vs "widened" (top-K unique palindromes by IR score,
    default K=25). Wider pool gives the selector access to lower-IR-score
    real operators that V0 filters out.

  * **rerank_scorer** — name from `algorithms.benchmark.operator_scorers.SCORERS`,
    *or* the special string "consensus_score" to fall back to V0's legacy
    selection metric (cross-homolog conservation of the palindrome
    alignment). Recommended scorers from the gradient:
      - "gc_weak"        : 0.2 · AT + 0.5 · length
      - "gc_medium_weak" : 0.5 · AT + 0.3 · length
      - "gc_strong"      : 1.0 · AT + 0.3 · length
    `gc_strong` was originally named `combined_v2` in PR #14.

The ablation matrix in `algorithms/versions.py` (V2.1 through V2.7)
exercises every meaningful combination of these two knobs against the
TetR-137 benchmark. The legacy combination (candidate_strategy="legacy",
rerank_scorer="consensus_score") reproduces V0 exactly.

The returned dict has the same shape as v0, plus three new fields:
  * `rerank_score`           — selector score of the winning candidate
  * `rerank_scorer`          — which scorer was used (for provenance)
  * `n_candidates_evaluated` — how many candidates produced a valid
                               cross-homolog consensus

Non-inverted-repeats search modes ("Align an input sequence", "Scan
entire promoter region") fall through to v0 unchanged.

Tunable params (optional, sensible defaults):
  * `candidate_strategy` (default "widened"): "legacy" | "widened"
  * `top_k_palindromes`  (default 25): pool size when widened
  * `rerank_scorer`      (default "gc_strong"): scorer name OR
                          the literal "consensus_score" to use V0's
                          selection metric
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
    findBestPalindrome,
    findImperfectPalindromes,
    findOperatorInIntergenic,
    getConsensus,
    get_consensus_score,
    generate_frequency_matrix,
)

from algorithms.benchmark.operator_scorers import SCORERS  # noqa: E402

from . import v0 as v0_mod


_DEFAULT_TOP_K = 25
_DEFAULT_CANDIDATE_STRATEGY = "widened"
_DEFAULT_RERANK = "gc_strong"
_CONSENSUS_SCORE_SENTINEL = "consensus_score"


# --- Helpers -------------------------------------------------------------

def _scorer_by_name(name: str) -> Callable[[str], float]:
    for n, fn, _ in SCORERS:
        if n == name:
            return fn
    raise KeyError(f"unknown rerank scorer {name!r}; choose from "
                    f"{[n for n, _, _ in SCORERS]} (or "
                    f"{_CONSENSUS_SCORE_SENTINEL!r} to use V0's metric)")


def _enumerate_palindromes(intergenic: str, shortest: int, longest: int,
                            win_score: int, loss_score: int,
                            spacer_penalty: dict) -> list[dict]:
    """Return *every* imperfect palindrome across the size range as
    `{seq, score}` dicts — no per-size or global max filter applied."""
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
    """Union across all sizes, dedupe by uppercase seq, top-K by score."""
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


def _legacy_candidates(intergenic: str, shortest: int, longest: int,
                        win_score: int, loss_score: int,
                        spacer_penalty: dict) -> list[dict]:
    """Reproduce V0's narrow candidate pool: palindromes tied at the global
    max IR score across the size range. Often just 1-5 candidates."""
    try:
        ops = findBestPalindrome(
            intergenic=intergenic.upper(), shortest=shortest, longest=longest,
            winScore=win_score, lossScore=loss_score, sPenalty=spacer_penalty)
    except Exception:
        return []
    if not ops:
        return []
    # Dedupe by seq even within legacy output — sometimes ties produce
    # duplicates in different orientations / spacer offsets.
    seen: set[str] = set()
    out: list[dict] = []
    for p in ops:
        key = (p.get("seq") or "").upper()
        if not key or key == "NONE" or key in seen:
            continue
        seen.add(key)
        out.append(p)
    return out


def _consensus_seq(consensus: dict) -> str:
    md = consensus.get("motif_data") or []
    return "".join(x["base"] for x in md if isinstance(x, dict))


# --- Public interface ----------------------------------------------------

def fetch(homolog_metadata: list[dict], params: dict) -> dict:
    """Same I/O contract as v0. See module docstring for tunable params."""

    # Non-inverted-repeats modes are unchanged from v0.
    if params.get("search_method") != "Look for inverted repeats":
        return v0_mod.fetch(homolog_metadata, params)

    ext_length = params["extension_length"]
    acc = homolog_metadata[0]["Uniprot Id"]
    regulated_seqs = [h["promoter"] for h in homolog_metadata]
    reference_candidates = [s for s in regulated_seqs if s]
    if not reference_candidates:
        return v0_mod.fetch(homolog_metadata, params)
    reference_seq = reference_candidates[0]

    strategy = params.get("candidate_strategy", _DEFAULT_CANDIDATE_STRATEGY)
    top_k = int(params.get("top_k_palindromes", _DEFAULT_TOP_K))
    rerank_name = params.get("rerank_scorer", _DEFAULT_RERANK)

    use_consensus_score_selector = (rerank_name == _CONSENSUS_SCORE_SENTINEL)
    rerank_fn: Optional[Callable[[str], float]] = (
        None if use_consensus_score_selector
        else _scorer_by_name(rerank_name))

    if strategy == "legacy":
        candidates = _legacy_candidates(
            intergenic=reference_seq,
            shortest=params["min_operator_length"],
            longest=params["max_operator_length"],
            win_score=params["win_score"],
            loss_score=params["loss_score"],
            spacer_penalty=params["spacer_penalty"],
        )
    elif strategy == "widened":
        candidates = _top_k_unique_palindromes(
            intergenic=reference_seq,
            shortest=params["min_operator_length"],
            longest=params["max_operator_length"],
            win_score=params["win_score"],
            loss_score=params["loss_score"],
            spacer_penalty=params["spacer_penalty"],
            k=top_k,
        )
    else:
        raise ValueError(f"unknown candidate_strategy {strategy!r}; "
                          f"expected 'legacy' or 'widened'")

    operator_data = {
        "Uniprot Id": str(acc),
        "aligned_seq": "None",
        "num_seqs": "None",
        "consensus_score": 0,
        "rerank_score": -math.inf,
        "rerank_scorer": rerank_name,
        "candidate_strategy": strategy,
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

        motif_str = _consensus_seq(consensus)
        if use_consensus_score_selector:
            # Reproduce V0's selection metric exactly.
            selector_value = float(consensus_score)
        else:
            selector_value = (rerank_fn(motif_str) if motif_str
                               else float("-inf"))

        op_seq = findOperatorInIntergenic(
            reference_seq, candidate["seq"], params)
        native_operator = (op_seq["operator"] if op_seq is not None
                            else candidate["seq"])

        if selector_value > operator_data["rerank_score"]:
            operator_data["rerank_score"] = selector_value
            operator_data["consensus_score"] = consensus_score
            operator_data["native_operator"] = native_operator
            operator_data["consensus_seq"] = motif_str
            operator_data["num_seqs"] = consensus["num_seqs"]
            operator_data["motif"] = consensus["motif_data"]
            operator_data["aligned_seqs"] = metrics
            operator_data["frequency_matrix"] = generate_frequency_matrix(metrics)

    return operator_data
