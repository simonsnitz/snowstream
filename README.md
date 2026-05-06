# snowstream

Predicts the DNA sequence that a prokaryotic transcription factor binds to.

## Architecture

- **`src/`** — original Python pipeline (BLAST → genome coordinates → operons → promoters). Untouched from the Streamlit version.
- **`backend/`** — FastAPI app that wraps `src/` with an SSE endpoint and a per-query JSON cache. The expensive part of the pipeline lives here.
- **`frontend/`** — Vite + React + MUI app. The operator-extraction step (Smith-Waterman alignment, consensus motif, inverted-repeat search) runs entirely in the browser, so tweaking search parameters re-renders instantly without re-hitting the backend.
- **`cache/`** — gitignored, holds `<key>.json` per cached query.
- **`streamlit_app.py`** — original Streamlit UI; still runs.

The pipeline cache key hashes the inputs that affect the predicted promoters (input value, BLAST params, promoter params). Search-method and alignment params do *not* affect the cache because operator extraction happens client-side.

## Prereqs

- Python 3.10 (the pinned dependencies don't build on 3.12).
- `diamond` on `$PATH` (`brew install diamond`).
- A diamond database at `../databases/bHTH_RefSeq.dmnd` (relative to this repo). The `src/blast.py` code looks for it there.
- Node 18+ for the frontend.

## Develop

Two processes — backend on `:8000`, frontend on `:5173`. Vite proxies `/api/*` to the backend.

```bash
# one-time setup
python3.10 -m venv .venv
.venv/bin/pip install -r requirements.txt
(cd frontend && npm install)

# terminal 1 — backend
.venv/bin/uvicorn backend.main:app --port 8000 --reload

# terminal 2 — frontend
(cd frontend && npm run dev)
```

Open <http://localhost:5173>.

## Operator-extraction methods

Methods live in `frontend/src/operator/methods/` and are registered in `frontend/src/operator/index.js`. Each method exports `{ label, params, run(promoters, params) → candidates }`. The shared `pipeline.js` aligns each candidate against every homolog promoter (Smith-Waterman with affine gaps), builds a consensus motif, and picks the best-scoring candidate.

To add a method (e.g. BioMSA): drop a file in `methods/`, register it in `index.js`. The UI exposes it automatically via the search-method dropdown in advanced options.

## Verifying the JS port against Python

The JS Smith-Waterman + consensus implementation is regression-tested against fixtures generated from `Bio.pairwise2` and `src/fetch_operator.py`. Run:

```bash
node frontend/src/operator/__tests__/verify.mjs
```

To regenerate fixtures after a Python-side change, recreate `backend/_dump_reference.py` from git history and run it (the script is a one-off and has been removed from the repo).

## Production

`npm run build` in `frontend/` produces `frontend/dist/`. Easiest deployment is to mount that directory as static files from FastAPI (not configured yet — uses two processes in dev).

## Original Streamlit app

```bash
streamlit run streamlit_app.py
```

Still functional; will be removed once the React/FastAPI version is at parity.
