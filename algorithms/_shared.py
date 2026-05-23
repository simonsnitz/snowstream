"""Helpers shared across algorithm versions and the benchmark harness.

This module centralises NCBI HTTP access and the alignment scorer so that
each algorithm version can focus on its own logic. None of the per-version
files should call NCBI directly — go through `ncbi_efetch_nuccore` here.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.http_utils import ncbi_get  # noqa: E402  (re-exported for versions)

__all__ = [
    "ncbi_efetch_nuccore",
    "rc",
    "best_window_match",
]


def ncbi_efetch_nuccore(genome_id: str, start: int, stop: int,
                         strand: int = 1, rettype: str = "fasta") -> str | None:
    """Thin wrapper around NCBI eFetch on the nuccore database.

    `strand` is 1 (forward) or 2 (reverse-complement); matches the live
    fetch_promoter convention of always passing strand=1 unless the caller
    explicitly wants the reverse complement.
    """
    if stop <= start:
        return None
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
           f"db=nuccore&id={genome_id}&seq_start={start}&seq_stop={stop}"
           f"&strand={strand}&rettype={rettype}")
    try:
        r = ncbi_get(url, timeout=30)
    except Exception:
        return None
    if not r.ok:
        return None
    lines = r.text.split("\n")
    return "".join(lines[1:]).strip()


def rc(s: str) -> str:
    """Reverse-complement a DNA string. Preserves case."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G",
            "a": "t", "t": "a", "g": "c", "c": "g",
            "N": "N", "n": "n"}
    return "".join(comp.get(b, b) for b in reversed(s))


def best_window_match(needle: str, haystack: str, min_overlap: int = 8) -> dict:
    """Slide `needle` across `haystack` (forward and RC of needle). Return
    the offset with the most matching bases.

    Used by the benchmark harness to score alignments between known operators
    and predicted/extracted regions. Ungapped because the inputs are short
    (20-50 bp operators against 80-800 bp promoters) and gaps tend to inflate
    spurious matches.
    """
    if not needle or not haystack:
        return _empty_match()
    haystack_u = haystack.upper()
    best = None
    for orient, q in (("fwd", needle.upper()), ("rc", rc(needle).upper())):
        nq = len(q)
        for offset in range(-nq + 1, len(haystack_u)):
            p_start = max(0, offset)
            p_end = min(len(haystack_u), offset + nq)
            q_start = p_start - offset
            q_end = q_start + (p_end - p_start)
            aligned_len = p_end - p_start
            if aligned_len < min_overlap:
                continue
            ha = haystack_u[p_start:p_end]
            na = q[q_start:q_end]
            matches = sum(1 for x, y in zip(ha, na) if x == y)
            cand = {
                "orientation": orient,
                "offset": offset,
                "matches": matches,
                "aligned_len": aligned_len,
                "identity_pct": round(100 * matches / aligned_len, 1),
                "coverage_pct": round(100 * aligned_len / nq, 1),
                "haystack_excerpt": ha,
                "needle_excerpt": na,
                "_rank": (matches, matches - (aligned_len - matches)),
            }
            if best is None or cand["_rank"] > best["_rank"]:
                best = cand
    if best is None:
        return _empty_match()
    best.pop("_rank", None)
    return best


def _empty_match() -> dict:
    return {"identity_pct": 0, "matches": 0, "aligned_len": 0,
            "coverage_pct": 0, "orientation": "fwd", "offset": 0,
            "haystack_excerpt": "", "needle_excerpt": ""}
