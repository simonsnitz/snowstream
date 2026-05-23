# algorithms/ — versioned algorithms + benchmark harness

This directory is the canonical home for the three core search algorithms that
make up the Snowstream operator-prediction pipeline. Each algorithm has its
own subdirectory containing one file per version, and a top-level
`versions.py` wires them into named bundles (V0, V1, …) that the benchmark
harness runs end-to-end against a standardised dataset.

The goal: every meaningful change to one of the algorithms ships as a new
version, with its own benchmark run, so we can track performance over time,
compare changes head-to-head, and roll back if a "fix" turns out to regress
something.

## Layout

```
algorithms/
├── README.md                       — this file
├── versions.py                     — VERSIONS dict mapping name → bundle
├── _shared.py                      — alignment scoring + NCBI helpers
├── operon_fetch/
│   ├── __init__.py
│   └── v0.py                       — original acc2operon (genome fetch + walk)
├── promoter_fetch/
│   ├── __init__.py
│   ├── v0.py                       — original single-candidate (legacy)
│   └── v1.py                       — multi-candidate enumeration
├── operator_fetch/
│   ├── __init__.py
│   └── v0.py                       — original palindrome-finder
├── pipeline.py                     — version-aware end-to-end runner
└── benchmark/
    ├── __init__.py
    ├── run.py                      — entrypoint: python -m algorithms.benchmark.run
    ├── metrics.py                  — per-protein metric computation
    ├── report.py                   — markdown + PNG report generation
    ├── datasets/
    │   └── tetr_137.json           — default 137-protein dataset
    └── results/                    — timestamped output folders (gitignored)
```

## How a version is defined

`versions.py` declares each version as a bundle of (operon_fetch, promoter_fetch,
operator_fetch) function references plus metadata describing what changed:

```python
VERSIONS = {
    "v0": {
        "operon_fetch":   operon_fetch.v0.fetch,
        "promoter_fetch": promoter_fetch.v0.fetch,
        "operator_fetch": operator_fetch.v0.fetch,
        "description":    "Original algorithms (legacy behaviour).",
    },
    "v1": {
        "operon_fetch":   operon_fetch.v0.fetch,
        "promoter_fetch": promoter_fetch.v1.fetch,
        "operator_fetch": operator_fetch.v0.fetch,
        "description":    "V0 + multi-candidate promoter enumeration.",
    },
}
```

The pipeline.py runner is version-aware — it inspects what shape the
`promoter_fetch` returns (single string vs list of candidates) and dispatches
accordingly. This lets V1's multi-candidate promoter feed V0's operator_fetch
without changing the operator_fetch interface.

## Adding a new version

1. **Implement the change** as a new file in the relevant `*_fetch/` subdir
   (e.g. `promoter_fetch/v2.py`).
2. **Register it** in `versions.py` as `"v2"` with a description.
3. **Run the benchmark**:
   ```bash
   python -m algorithms.benchmark.run --versions v0 v1 v2
   ```
4. **Review the report** in `algorithms/benchmark/results/<timestamp>/`.
5. **Commit** the new version + the benchmark report.

## Running the benchmark

```bash
# Default: V0 + V1 on the TetR 137-protein dataset
python -m algorithms.benchmark.run

# Specific versions
python -m algorithms.benchmark.run --versions v0 v1

# Different dataset
python -m algorithms.benchmark.run --dataset path/to/my_dataset.json

# Limit for smoke testing
python -m algorithms.benchmark.run --max-proteins 10
```

Output goes to `algorithms/benchmark/results/<YYYY-MM-DD_HH-MM-SS>/`:
- `per_protein.json` — every protein's metrics for every version
- `summary.json` — aggregated counts at each identity threshold
- `report.md` — human-readable comparison
- `chart_*.png` — distribution histograms per metric

## Datasets

Default: `benchmark/datasets/tetr_137.json` — 137 TetR-family proteins with
literature-backed operators, sourced from groovDB + the Snowprint paper.

Custom datasets must be JSON with this shape:

```json
[
  {
    "ncbi_accession": "WP_001224188.1",
    "uniprot_id": "P0A9R7",          // optional
    "alias": "GbaA",                 // optional
    "operators": [
      {"sequence": "ATAAACGGAGAGTTATCCGTTTGT", "doi": "10.1016/..."}
    ],
    "sources": ["snowprint"]          // optional provenance tag
  },
  ...
]
```

Pass `--dataset path/to/file.json` to use a custom one.

## Metrics

For each protein × version the benchmark records:

| Metric | What it measures |
|---|---|
| `n_homologs` | Number of homologs returned by snowstream smart-lookup |
| `n_homologs_with_promoter` | Subset where promoter_fetch returned non-empty |
| `query_promoter_match` | Best alignment (id%) of any known operator vs the query's promoter |
| `any_homolog_match` | Best alignment across all homolog promoters |
| `algorithm_match` | Best alignment of the predicted motif vs known operator |

Aggregated at thresholds `[60%, 70%, 80%, 90%]` per metric per version.

Extending the metric set: add a new function to `benchmark/metrics.py`
(takes `(per_protein_record, version_output) -> float | int`) and register
it in the `METRICS` list. The report generator picks up new metrics
automatically — they appear as new rows in the comparison table.
