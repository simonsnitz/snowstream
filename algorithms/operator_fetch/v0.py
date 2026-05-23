"""operator_fetch v0 — original inverted-repeat palindrome finder.

Port of src/fetch_operator.py with the streamlit decorator and unused
imports removed. Algorithm logic is byte-identical.

Wraps the legacy `fetch_operator` so each algorithm version can be tested
in isolation. We re-export the legacy implementation rather than copying
hundreds of lines — the algorithm itself is unchanged in V0.
"""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.fetch_operator import fetch_operator as _legacy_fetch_operator  # noqa: E402


def fetch(homologs: list[dict], params: dict) -> dict:
    """Original Snowprint operator finder. See module docstring for I/O."""
    return _legacy_fetch_operator(homologs, params)
