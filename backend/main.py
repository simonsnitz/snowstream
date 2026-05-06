from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse

from . import cache
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
