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

Backend tests:

```bash
.venv/bin/python -m unittest discover -s backend -p 'test_*.py'
```

## Family precomputed caches

Each family lives at `families/<key>/` with a tracked `family.json` manifest
and a gitignored `predictions.jsonl` of cached predictions (one per
representative protein, keyed by UniProt ID). Lookups are O(1) via an
in-memory byte-offset index that rebuilds when the JSONL changes on disk.

Endpoints:

- `GET /api/families` — list configured families with their manifests.
- `GET /api/cached-by-uniprot/{family}/{uniprot_id}` — direct lookup of a
  cached prediction.

### One-time: build a local NCBI nr Diamond DB

The precompute BLASTs each representative against NCBI's full nr database. You
need a local Diamond DB built from nr — this is a one-time, hours-long
operation that produces a ~50 GB index.

```bash
# Pick somewhere with ~250 GB free for the download + intermediate + final DB.
cd /path/to/big/disk

# Download nr (gzipped FASTA, ~150 GB).
curl -O https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

# Build the Diamond DB (~hours; ~50 GB output).
diamond makedb --in nr.gz -d nr
```

Then export the path before running the precompute:

```bash
export SNOWPRINT_DIAMOND_DB=/path/to/big/disk/nr.dmnd
```

If `SNOWPRINT_DIAMOND_DB` is unset, the live `/api/predict` and the
precompute both fall back to the small `bHTH_RefSeq.dmnd` shipped under
`databases/`.

### Smart-lookup against precomputed family caches

Before running the full pipeline, `/api/predict` BLASTs the query sequence
against every family's `representatives.dmnd`. If the top hit clears the
family's `smart_lookup` thresholds (`min_identity_pct` and `min_coverage_pct`
in `family.json`; default 50/95), we return that representative's cached
prediction immediately, wrapped with a `matched_via` block containing the
family key, the matched UniProt ID, and the identity/coverage percentages.

The frontend renders this in the cache banner: *"Matched cached
representative P43506 (87 % identity, 99 % coverage) — family tetr.
Characterized in groovDB ↗"* (or a list of papers for PaperBLAST-backed
representatives, or no characterisation line for default-backed ones).

The "Re-run with full BLAST" button on that banner sets `force=true` in the
request, which bypasses both the smart-lookup *and* the parameter-hash
cache, forcing a fresh pipeline run against whichever database is
configured.

The smart-lookup never runs when the family has no `representatives.dmnd`
yet (i.e. before PR 2's precompute is finished). It also never runs when
`force=true` — so existing caches and the parameter-hash flow still work
exactly as before.

### Database choice in the live UI

`/api/predict` accepts a `database` field (`"local_diamond"` by default,
`"nr_remote"` to BLAST against NCBI's nr remotely). The Advanced options
drawer surfaces this as a radio in the BLAST section.

The remote path submits a job to NCBI's BLAST URL API, polls every 30 s
(emitting heartbeat SSE events so the connection stays alive), fetches
the BLAST XML, and translates RefSeq accessions to UniProt accessions
via UniProt's xref query. Hits without a UniProt mapping are dropped.
Wall-clock per query is typically 5–30 minutes; if NCBI is slow it can
exceed an hour, in which case the request times out and the user can
retry. Closing the browser tab kills the SSE stream but the BLAST job
continues to completion at NCBI — re-submitting the same query means
re-running, since the cache key includes the `database` choice.

### One-time: download the groovDB dump

groovDB-known regulators are preferred as cluster representatives. Download
the full dump (a single JSON file, button on the
[programmatic-access page](https://www.groov.bio/about/programmatic-access))
and save it to `families/_resources/groovdb.json`. The precompute reads it
into an in-memory `{uniprot_id: entry}` table.

If the dump is missing the precompute still runs; representatives without a
groovDB or PaperBLAST hit fall back to the MMseqs2 centroid with no
characterisation evidence attached.

### Running the precompute

```bash
# All stages, end-to-end (resumable: rerun is idempotent):
python -m scripts.precompute_family tetr

# A single stage:
python -m scripts.precompute_family tetr --stage representatives

# Smoke-test on a small slice (5 clusters):
python -m scripts.precompute_family tetr --max-clusters 5

# Re-run a stage that previously completed:
python -m scripts.precompute_family tetr --stage cluster --force

# Skip PaperBLAST (it's behind Cloudflare's JS challenge — you may need to
# rerun representative selection later when that's resolved):
python -m scripts.precompute_family tetr --no-paperblast
```

Stages, in order:

1. **members** — fetch UniProt accessions annotated with the family's InterPro entry
2. **sequences** — batch-fetch FASTA sequences for those accessions
3. **cluster** — MMseqs2 cluster at 50 % identity (configurable via `family.json`)
4. **representatives** — pick a rep per cluster: groovDB > PaperBLAST > MMseqs2 centroid
5. **dmnd** — extract `representatives.fasta` and build `representatives.dmnd`
6. **predict** — run the full Snowprint pipeline per rep, append to `predictions.jsonl`

Each stage's output is a file in `families/<key>/`; the next stage picks it
up. Crash-resume: rerunning the script skips any stage whose output already
exists, and the `predict` stage resumes by scanning `predictions.jsonl` for
already-written UniProt IDs.

Run the tests:

```bash
.venv/bin/python -m unittest \
    backend.test_families \
    scripts.precompute.test_representatives \
    scripts.precompute.test_other_stages
```

#### Known caveats

- **PaperBLAST is behind Cloudflare.** Plain HTTP requests with a browser
  User-Agent often get through, but the endpoint can flip to a JS challenge.
  When that happens the precompute logs a warning and treats the cluster as
  having no PaperBLAST evidence — the representative falls back to the
  MMseqs2 centroid (or to a groovDB hit if present). Re-running the
  `representatives` stage with `--force` once PaperBLAST is reachable will
  pick up the literature.
- **PaperBLAST HTML parsing is regex-based.** If their output format changes,
  hits will silently stop being recognised; the cluster falls back to its
  centroid. Tests cover the empty and challenge-page paths, not the format.

## Production

`npm run build` in `frontend/` produces `frontend/dist/`. Easiest deployment is to mount that directory as static files from FastAPI (not configured yet — uses two processes in dev).

## Original Streamlit app

```bash
streamlit run streamlit_app.py
```

Still functional; will be removed once the React/FastAPI version is at parity.
