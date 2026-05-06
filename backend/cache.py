import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from .schemas import PredictRequest

CACHE_DIR = Path(__file__).resolve().parent.parent / "cache"


def _canonical_input(req: PredictRequest) -> dict:
    # Only fields that affect the output. get_coordinates_method does not.
    # `database` is intentionally only included when it differs from the
    # default — so caches written before the database choice was introduced
    # remain reachable for `database="local_diamond"` requests.
    payload: dict = {
        "input_method": req.input_method,
        "input_value": req.input_value,
        "blast_params": req.blast_params.model_dump(),
        "promoter_params": req.promoter_params.model_dump(),
    }
    if req.database != "local_diamond":
        payload["database"] = req.database
    return payload


def compute_key(req: PredictRequest) -> str:
    payload = json.dumps(_canonical_input(req), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def _path(key: str) -> Path:
    return CACHE_DIR / f"{key}.json"


def read(key: str) -> Optional[dict]:
    p = _path(key)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except (json.JSONDecodeError, OSError):
        return None


def write(key: str, payload: dict) -> dict:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    payload = {
        **payload,
        "cache_key": key,
        "cached_at": datetime.now(timezone.utc).isoformat(),
    }
    _path(key).write_text(json.dumps(payload, indent=2))
    return payload


def delete(key: str) -> bool:
    p = _path(key)
    if p.exists():
        p.unlink()
        return True
    return False
