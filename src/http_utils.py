"""Shared HTTP helpers for pipeline calls.

Centralises three concerns the rest of `src/` shouldn't repeat:
  * NCBI rate limiting via a global semaphore (3/sec anonymous, 10/sec with key).
  * Optional NCBI API key (env: `NCBI_API_KEY`) appended to eutils URLs.
  * A safe default `timeout` so a hung server can't stall a worker indefinitely.

Concurrent NCBI request cap defaults to 8 (with key) or 2 (without). Override
with `SNOWPRINT_NCBI_CONCURRENT`.
"""

from __future__ import annotations

import os
import threading

import requests

NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "").strip()
DEFAULT_TIMEOUT = float(os.environ.get("SNOWPRINT_HTTP_TIMEOUT", "30"))

_default_ncbi_concurrent = "8" if NCBI_API_KEY else "2"
_NCBI_SEM = threading.BoundedSemaphore(
    int(os.environ.get("SNOWPRINT_NCBI_CONCURRENT", _default_ncbi_concurrent))
)


def _with_api_key(url: str) -> str:
    if not NCBI_API_KEY or "api_key=" in url:
        return url
    sep = "&" if "?" in url else "?"
    return f"{url}{sep}api_key={NCBI_API_KEY}"


def ncbi_get(url: str, timeout: float | None = None) -> requests.Response:
    """GET an NCBI eutils URL — adds api_key, enforces global concurrency cap,
    applies a bounded timeout. Blocks (briefly) when the semaphore is full."""
    final_url = _with_api_key(url)
    with _NCBI_SEM:
        return requests.get(final_url, timeout=timeout or DEFAULT_TIMEOUT)


def http_get(url: str, timeout: float | None = None) -> requests.Response:
    """GET a non-NCBI URL with a bounded timeout."""
    return requests.get(url, timeout=timeout or DEFAULT_TIMEOUT)
