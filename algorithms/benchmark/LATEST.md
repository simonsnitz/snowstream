# Latest benchmark result

> Auto-generated on every `algorithms.benchmark.run`. For per-run history (with PNG charts) see `algorithms/benchmark/results/<timestamp>/`.

- **Run**: `2026-05-22_14-55-10`
- **Dataset**: `tetr_137.json`  (N = 137; 135 resolved)
- **Versions**: `v0`, `v1`

## Identity-metric threshold counts

### Known operator in query promoter  (n = 135)

| Threshold | v0 | v1 |
|---|---:|---:|
| ≥ 60% | 103 (76.3%) | 104 (77.0%) |
| ≥ 70% | 92 (68.1%) | 93 (68.9%) |
| ≥ 80% | 85 (63.0%) | 87 (64.4%) |
| ≥ 90% | 79 (58.5%) | 81 (60.0%) |
| mean | 82.7% | 83.4% |

### Known operator in any homolog promoter  (n = 135)

| Threshold | v0 | v1 |
|---|---:|---:|
| ≥ 60% | 125 (92.6%) | 125 (92.6%) |
| ≥ 70% | 109 (80.7%) | 111 (82.2%) |
| ≥ 80% | 98 (72.6%) | 99 (73.3%) |
| ≥ 90% | 84 (62.2%) | 85 (63.0%) |
| mean | 88.5% | 88.8% |

### Known operator in predicted motif  (n = 135)

| Threshold | v0 | v1 |
|---|---:|---:|
| ≥ 60% | 69 (51.1%) | 73 (54.1%) |
| ≥ 70% | 58 (43.0%) | 58 (43.0%) |
| ≥ 80% | 46 (34.1%) | 46 (34.1%) |
| ≥ 90% | 31 (23.0%) | 31 (23.0%) |
| mean | 65.1% | 65.4% |

## Count metrics

| Metric | v0 | v1 |
|---|---:|---:|
| # homologs (mean) | 76.1 | 76.1 |
| # homologs with promoter (mean) | 75.1 | 75.8 |

## Head-to-head: v1 − v0 (predicted-motif identity to known operator)

### v1 wins (showing 10 of 10)

| NCBI | Alias | v0 | v1 | Δ |
|---|---|---:|---:|---:|
| WP_012480398.1 | MfsR | 41.7% | 72.7% | +31.0% |
| AAC01726.1 | RifQ | 30.8% | 45.5% | +14.7% |
| WP_012392470.1 | FdmR | 50.0% | 61.5% | +11.5% |
| WP_003107570.1 | DesT | 93.8% | 100.0% | +6.2% |
| WP_011060270.1 | psrA | 55.0% | 60.0% | +5.0% |
| WP_011015509.1 | McbR | 56.2% | 61.1% | +4.9% |
| ABK70852.1 | DarR | 33.3% | 37.5% | +4.2% |
| WP_050417930.1 | MbdR | 41.2% | 44.8% | +3.6% |
| AMR55379.1 | NdpR | 39.4% | 41.7% | +2.3% |
| AAB62296.1 | CymR | 42.9% | 45.0% | +2.1% |

### v1 losses (showing 8 of 8)

| NCBI | Alias | v0 | v1 | Δ |
|---|---|---:|---:|---:|
| WP_013222564.1 | RifQ | 52.9% | 38.5% | -14.4% |
| AAD13556.1 | LanK | 72.7% | 63.6% | -9.1% |
| AFA75826.1 | LcpRVH2 | 54.5% | 45.5% | -9.0% |
| WP_003856101.1 | AcnR | 50.0% | 44.4% | -5.6% |
| WP_009949582.1 | PccD | 57.1% | 52.4% | -4.7% |
| WP_011030795.1 | ScbR2 | 41.4% | 38.5% | -2.9% |
| WP_033220308.1 | VarR | 30.0% | 27.3% | -2.7% |
| WP_067435216.1 | SLCG_2919 | 34.8% | 33.3% | -1.5% |

**Net direction**: v1 wins on 10, loses on 8 (out of 18 proteins where versions disagree by >1%).
