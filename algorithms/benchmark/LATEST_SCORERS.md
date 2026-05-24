# Operator-scorer classifier benchmark

- **Run**: `2026-05-24_11-52-57`
- **Positives**: 203 known operators from `operator_positives.json`
- **Negatives**: 42 spurious predicted motifs from `operator_negatives.json` (identity to known < 50.0%)

Each scorer returns a single float per sequence. Higher should mean more operator-like. We rank by ROC AUC — AUC = 0.5 is random, 1.0 is perfect separation.

## Ranking (by ROC AUC)

| Rank | Scorer | AUC | Cohen's d | Best threshold | TP@thr | FP@thr | Precision |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | `combined_v2` | 0.916 | 2.15 | 0.514 | 191/203 (94%) | 9/42 (21%) | 96% |
| 2 | `combined_v1` | 0.893 | 2.06 | -0.436 | 179/203 (88%) | 10/42 (24%) | 95% |
| 3 | `at_content` | 0.868 | 1.96 | 0.414 | 182/203 (90%) | 10/42 (24%) | 95% |
| 4 | `gc_penalty` | 0.868 | 1.96 | -0.586 | 182/203 (90%) | 10/42 (24%) | 95% |
| 5 | `palindrome_arm_at` | 0.842 | 1.74 | 0.389 | 182/203 (90%) | 11/42 (26%) | 94% |
| 6 | `cg_dinucleotide_penalty` | 0.797 | 1.36 | -0.081 | 136/203 (67%) | 8/42 (19%) | 94% |
| 7 | `combined_v3` | 0.780 | 1.23 | 0.226 | 188/203 (93%) | 18/42 (43%) | 91% |
| 8 | `length_preference` | 0.734 | 0.96 | 0.458 | 167/203 (82%) | 15/42 (36%) | 92% |
| 9 | `shannon_entropy` | 0.727 | 1.30 | 1.802 | 176/203 (87%) | 17/42 (40%) | 91% |
| 10 | `tandem_repeat_penalty` | 0.688 | 0.70 | 0.000 | 171/203 (84%) | 19/42 (45%) | 90% |
| 11 | `best_palindrome_spacer_score` | 0.535 | 0.20 | 4.000 | 174/203 (86%) | 33/42 (79%) | 84% |
| 12 | `core_vs_flank_at` | 0.325 | -0.58 | 0.667 | 1/203 (0%) | 0/42 (0%) | 100% |
| 13 | `palindrome_strength` | 0.188 | -1.27 | 2.000 | 203/203 (100%) | 42/42 (100%) | 83% |

## Per-scorer detail

| Scorer | Positives mean (median) | Negatives mean (median) | Δ mean |
|---|---:|---:|---:|
| `combined_v2` | 0.784 (0.791) | 0.412 (0.416) | +0.372 |
| `combined_v1` | -0.067 (-0.009) | -0.716 (-0.799) | +0.649 |
| `at_content` | 0.578 (0.591) | 0.292 (0.223) | +0.285 |
| `gc_penalty` | -0.422 (-0.409) | -0.708 (-0.777) | +0.285 |
| `palindrome_arm_at` | 0.586 (0.600) | 0.287 (0.203) | +0.299 |
| `cg_dinucleotide_penalty` | -0.064 (-0.059) | -0.153 (-0.146) | +0.089 |
| `combined_v3` | 0.490 (0.500) | 0.237 (0.206) | +0.253 |
| `length_preference` | 0.688 (0.755) | 0.399 (0.325) | +0.289 |
| `shannon_entropy` | 1.895 (1.926) | 1.719 (1.738) | +0.175 |
| `tandem_repeat_penalty` | -0.236 (0.000) | -0.690 (-1.000) | +0.454 |
| `best_palindrome_spacer_score` | 3.635 (4.000) | 3.429 (4.000) | +0.207 |
| `core_vs_flank_at` | -0.079 (-0.067) | 0.030 (0.054) | -0.109 |
| `palindrome_strength` | 15.586 (16.000) | 21.333 (21.000) | -5.747 |

## Top scorer: `combined_v2`

AUC = **0.916**, Cohen's d = **2.15**.  At the optimal threshold (`0.514`) it captures **191/203 (94%)** of positives while only letting through **9/42 (21%)** of negatives. Precision at this threshold: **96%**.
