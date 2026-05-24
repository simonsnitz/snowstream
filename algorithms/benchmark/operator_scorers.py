"""Candidate-scoring functions for the operator-fetch classifier benchmark.

Each scorer is a function `f(seq: str) -> float`. Higher score = "more
operator-like" — we want positives to score higher than negatives.

Two registries:
  * SCORERS — individual feature scorers (one biological signal each).
  * COMBINED_SCORERS — weighted combinations of the individual scorers
    that we want to evaluate as composite scoring functions.

Adding a new scorer:
  1. Define a function taking `seq: str` and returning a float.
  2. Append `(name, fn, description)` to SCORERS.
  3. Run `python -m algorithms.benchmark.score_classifier`.

Sequence-intrinsic only — these scorers operate on a single candidate
string without any homolog alignment context. Cross-homolog conservation
features will need a separate scorer family that takes the full set of
homolog promoters as input (deferred to V2 work).
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path
from typing import Callable

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Reuse the legacy palindrome finder for the "palindrome strength" scorer.
from src.fetch_operator import findBestPalindrome  # noqa: E402

# Default spacer penalty matrix (same as deployed Snowprint / streamlit UI).
_DEFAULT_SPACER_PENALTY = {
    str(i): v for i, v in zip(range(21),
        [4, 4, 4, 4, 4, 2, 2, 0, 0, -2, -2,
         -4, -4, -6, -6, -8, -8, -10, -10, -12, -12])}


def _normalize(seq: str) -> str:
    """Uppercase + strip to ACGT only."""
    return re.sub(r"[^ACGT]", "", seq.upper())


# --- Individual scorers --------------------------------------------------

def at_content(seq: str) -> float:
    """Fraction A+T. TetR-family operators tend to be AT-rich; spurious
    palindromes that fool the current scorer are often GC-rich."""
    s = _normalize(seq)
    if not s:
        return 0.0
    return (s.count("A") + s.count("T")) / len(s)


def gc_penalty(seq: str) -> float:
    """Negative of GC content — directly equivalent to at_content but
    sign-flipped. Included for clarity in the report; AUC will be
    identical to at_content."""
    return -1.0 * (1.0 - at_content(seq))


def length_preference(seq: str, target: int = 22, sigma: float = 8.0) -> float:
    """Gaussian centered at `target` bp. Typical TetR operators are
    ~18-25 bp; the predicted motifs from operator_fetch (with the default
    5-bp extension) are usually a bit longer. Returns a value in (0, 1]."""
    s = _normalize(seq)
    if not s:
        return 0.0
    return math.exp(-((len(s) - target) ** 2) / (2 * sigma ** 2))


def palindrome_strength(seq: str) -> float:
    """Score of the best inverted-repeat palindrome found in the sequence,
    using the legacy `findBestPalindrome` machinery. This is essentially
    the score the current operator_fetch picks by — it's our baseline
    for comparison."""
    s = _normalize(seq)
    if len(s) < 12:
        return 0.0
    try:
        ops = findBestPalindrome(
            intergenic=s, shortest=5, longest=15,
            winScore=2, lossScore=-2, sPenalty=_DEFAULT_SPACER_PENALTY)
    except Exception:
        return 0.0
    if not ops:
        return 0.0
    return float(max(o["score"] for o in ops))


def palindrome_arm_at(seq: str) -> float:
    """AT content of the palindrome arms only (excluding the spacer in the
    middle). Real operators have AT-rich arms; GC-rich spurious palindromes
    score this low."""
    s = _normalize(seq)
    if len(s) < 12:
        return 0.0
    try:
        ops = findBestPalindrome(
            intergenic=s, shortest=5, longest=15,
            winScore=2, lossScore=-2, sPenalty=_DEFAULT_SPACER_PENALTY)
    except Exception:
        return 0.0
    if not ops:
        return 0.0
    # The legacy returned `seq` looks like ARM + spacer.lower() + ARM
    arm_seqs = []
    for op in ops:
        op_seq = op.get("seq", "")
        # Arms are the uppercase runs at start and end
        upper_prefix = ""
        for c in op_seq:
            if c.isupper():
                upper_prefix += c
            else:
                break
        upper_suffix = ""
        for c in reversed(op_seq):
            if c.isupper():
                upper_suffix = c + upper_suffix
            else:
                break
        arm_seqs.append(upper_prefix + upper_suffix)
    arms = max(arm_seqs, key=len) if arm_seqs else ""
    if not arms:
        return 0.0
    return (arms.count("A") + arms.count("T")) / len(arms)


def core_vs_flank_at(seq: str, core_frac: float = 0.5) -> float:
    """AT content of the central `core_frac` of the sequence minus AT of
    the flanks. Positive = AT-richer middle than edges (operator-like).

    Helps when the sequence is the predicted motif from operator_fetch —
    the algorithm extends ±5 bp around the palindrome, so the flanks are
    intergenic context. Real operators have AT-rich centres and somewhat
    GC-richer flanking promoter sequence."""
    s = _normalize(seq)
    if len(s) < 16:
        return 0.0
    core_len = max(8, int(len(s) * core_frac))
    start = (len(s) - core_len) // 2
    core = s[start:start + core_len]
    flank = s[:start] + s[start + core_len:]
    if not flank:
        return 0.0
    core_at = (core.count("A") + core.count("T")) / len(core)
    flank_at = (flank.count("A") + flank.count("T")) / len(flank)
    return core_at - flank_at


def cg_dinucleotide_penalty(seq: str) -> float:
    """Negative count of CG dinucleotides, normalised by sequence length.
    CG dinucleotides are under-represented in TetR-family operators (and
    in most prokaryotic regulatory regions) but are abundant in the
    GC-rich palindromes that fool the current scorer."""
    s = _normalize(seq)
    if len(s) < 2:
        return 0.0
    cg_count = sum(1 for i in range(len(s) - 1) if s[i:i + 2] == "CG")
    return -cg_count / max(1, len(s) - 1)


def tandem_repeat_penalty(seq: str) -> float:
    """Penalty for long mononucleotide runs (CCCC, AAAA, etc.). Pure-base
    runs in spurious palindromes inflate inverted-repeat scores without
    being biologically meaningful."""
    s = _normalize(seq)
    if not s:
        return 0.0
    longest_run = 0
    cur_base = None
    cur_len = 0
    for b in s:
        if b == cur_base:
            cur_len += 1
        else:
            cur_base = b
            cur_len = 1
        longest_run = max(longest_run, cur_len)
    return -max(0, longest_run - 3)


def shannon_entropy(seq: str) -> float:
    """Per-base Shannon entropy in bits (max = 2 for ACGT uniform).
    Operators have moderate complexity; tandem-repeat palindromes have
    artificially low entropy."""
    s = _normalize(seq)
    if not s:
        return 0.0
    n = len(s)
    counts = {b: s.count(b) for b in "ACGT"}
    h = 0.0
    for c in counts.values():
        if c == 0:
            continue
        p = c / n
        h -= p * math.log2(p)
    return h


def best_palindrome_spacer_score(seq: str) -> float:
    """The legacy spacer-length term from the existing scorer, applied to
    the best palindrome. Higher (less negative) for spacers of 0-7 bp —
    the typical operator spacer length."""
    s = _normalize(seq)
    if len(s) < 12:
        return 0.0
    try:
        ops = findBestPalindrome(
            intergenic=s, shortest=5, longest=15,
            winScore=2, lossScore=-2, sPenalty=_DEFAULT_SPACER_PENALTY)
    except Exception:
        return 0.0
    if not ops:
        return 0.0
    # Read the spacer length out of the returned palindrome string
    best = max(ops, key=lambda o: o["score"])
    op_seq = best.get("seq", "")
    # spacer = lowercase region in middle
    lower = "".join(c for c in op_seq if c.islower())
    return _DEFAULT_SPACER_PENALTY.get(str(len(lower)), -12)


# --- Combined scorers ----------------------------------------------------

def combined_v1(seq: str) -> float:
    """Linear combination of the strongest individual features.

    Weights chosen a priori (not tuned on the test set) based on the
    biological story:
      + 1.0  · at_content         — primary signal
      + 0.5  · core_vs_flank_at   — operator-shape signal
      + 0.3  · length_preference  — typical-size sanity check
      - 0.05 · |palindrome_strength|  — discourage runaway palindrome scores
      + 0.5  · cg_dinucleotide_penalty
    """
    return (
        1.0  * at_content(seq)
        + 0.5  * core_vs_flank_at(seq)
        + 0.3  * length_preference(seq)
        - 0.05 * abs(palindrome_strength(seq))
        + 0.5  * cg_dinucleotide_penalty(seq)
    )


def combined_v3(seq: str) -> float:
    """AT content + core-vs-flank + CG penalty. Drops length and palindrome
    strength to see whether shape-only features carry the signal."""
    return (
        1.0 * at_content(seq)
        + 0.7 * core_vs_flank_at(seq)
        + 0.5 * cg_dinucleotide_penalty(seq)
    )


# --- GC-penalty gradient (used by operator_fetch v1.x ablation matrix) ----
#
# Three increasingly aggressive variants of the AT-content / GC-penalty
# preference, ordered weak → medium_weak → strong. `gc_strong` is the
# scorer the first cut of operator_fetch v1 used (originally named
# `combined_v2` in PR #14). The weaker variants test whether less-
# aggressive AT preference avoids regressions on GC-rich-genome operators
# (Streptomyces, Mycobacterium, etc.) where the real operator isn't AT-rich.
#
# Both at_content and length_preference return values in [0, 1], so the
# coefficients describe the relative weighting directly.

def gc_weak(seq: str) -> float:
    """Light AT bias — length dominates; AT acts mostly as tiebreaker."""
    return 0.2 * at_content(seq) + 0.5 * length_preference(seq)


def gc_medium_weak(seq: str) -> float:
    """Moderate AT bias — AT and length contribute roughly equally."""
    return 0.5 * at_content(seq) + 0.3 * length_preference(seq)


def gc_strong(seq: str) -> float:
    """Strong AT bias — AT dominates. Equivalent to the old `combined_v2`
    from PR #14 that drove the first cut of operator_fetch v1. May over-
    apply on GC-rich-genome hosts."""
    return 1.0 * at_content(seq) + 0.3 * length_preference(seq)


# Backward-compatibility alias — `combined_v2` was the name PR #14's
# score-classifier benchmark + the first operator_fetch v1 used. Keep
# resolvable so existing benchmark JSON / docs still load, but prefer
# `gc_strong` going forward.
combined_v2 = gc_strong


# --- Registries ----------------------------------------------------------

# Each entry: (name, function, description)
SCORERS: list[tuple[str, Callable[[str], float], str]] = [
    ("at_content",                       at_content,
     "Fraction of A+T bases (higher = more operator-like)"),
    ("gc_penalty",                       gc_penalty,
     "Negative GC content (sign-flipped at_content; AUC identical)"),
    ("length_preference",                length_preference,
     "Gaussian centred at 22 bp; penalises very short or very long"),
    ("palindrome_strength",              palindrome_strength,
     "Best inverted-repeat score (the current algorithm's pick metric)"),
    ("palindrome_arm_at",                palindrome_arm_at,
     "AT content of the palindrome arms only"),
    ("core_vs_flank_at",                 core_vs_flank_at,
     "AT% of central 50% minus AT% of flanks"),
    ("cg_dinucleotide_penalty",          cg_dinucleotide_penalty,
     "Negative count of CG dinucleotides per base"),
    ("tandem_repeat_penalty",            tandem_repeat_penalty,
     "Penalty for long mononucleotide runs (AAAA, CCCC, ...)"),
    ("shannon_entropy",                  shannon_entropy,
     "Per-position Shannon entropy (max 2 bits)"),
    ("best_palindrome_spacer_score",     best_palindrome_spacer_score,
     "Spacer-length term from legacy scoring (favours 0-7 bp spacers)"),
    ("combined_v1",                      combined_v1,
     "Weighted combo: AT, core-vs-flank, length, -palindrome, -CG"),
    ("combined_v3",                      combined_v3,
     "AT + core-vs-flank + CG penalty (shape-only)"),
    ("gc_weak",                          gc_weak,
     "Light AT bias (0.2 · AT + 0.5 · length) — GC-penalty gradient: weak"),
    ("gc_medium_weak",                   gc_medium_weak,
     "Moderate AT bias (0.5 · AT + 0.3 · length) — GC-penalty gradient: medium-weak"),
    ("gc_strong",                        gc_strong,
     "Strong AT bias (1.0 · AT + 0.3 · length) — GC-penalty gradient: strong "
     "(was `combined_v2` in PR #14)"),
]
