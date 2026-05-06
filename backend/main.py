from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse

from . import cache
from .families import registry as family_registry
from .pipeline import run_pipeline
from .schemas import PredictRequest

app = FastAPI(title="Snowprint API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/api/health")
def health():
    return {"ok": True}


@app.post("/api/predict")
async def predict(req: PredictRequest):
    return StreamingResponse(
        run_pipeline(req),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


@app.get("/api/cache/{key}")
def cache_get(key: str):
    hit = cache.read(key)
    if hit is None:
        raise HTTPException(status_code=404, detail="not found")
    return hit


@app.delete("/api/cache/{key}")
def cache_delete(key: str):
    deleted = cache.delete(key)
    if not deleted:
        raise HTTPException(status_code=404, detail="not found")
    return {"deleted": True}


@app.get("/api/families")
def families_list():
    """Public manifest of every configured precomputed family."""
    return [f.manifest_dict() for f in family_registry.list()]


@app.get("/api/cached-by-uniprot/{family_key}/{uniprot_id}")
def cached_by_uniprot(family_key: str, uniprot_id: str):
    """Direct lookup of a precomputed prediction by representative UniProt ID."""
    if family_registry.get(family_key) is None:
        raise HTTPException(status_code=404, detail="unknown family")
    record = family_registry.lookup(family_key, uniprot_id)
    if record is None:
        raise HTTPException(status_code=404, detail="not found")
    return record
