# Herald-conditioned decoding of the erasure-converted rotated surface code

<!-- Summary: 3 sentences (session 2). -->

## Results

### Thresholds

<!-- BEGIN:thresholds -->
| `r_e` | decoder | per-shot (7,11) crossing, % [95% CI] | bracket, % | sign changes (unsaturated) | per-shot collapse `d ≥ 7`, % [95% CI] | gap collapse − pairwise (pp) | per-round collapse `d ≥ 7` (secondary), % [95% CI] |
|---|---|---|---|---|---|---|---|
| 0 | blind_mwpm | 0.98 [0.95, 1.02] | [0.91, 1.00] | 2 | 0.98 [0.96, 1.00] | -0.004 | 1.44 [1.36, 1.58] |
| 0 | herald_mwpm | 0.98 [0.93, 1.05] | [0.91, 1.00] | 1 | 0.96 [0.94, 0.98] | -0.018 | 1.37 [1.33, 1.46] |
| 0.5 | blind_mwpm | 1.10 [1.03, 1.14] | [1.10, 1.18] | 1 | 1.09 [1.07, 1.11] | -0.017 | 1.48 [1.44, 1.80] |
| 0.5 | herald_mwpm | 1.60 [1.50, 1.69] | [1.57, 1.68] | 1 | 1.57 [1.54, 1.60] | -0.025 | 2.33 [2.21, 2.73] |
| 0.98 | blind_mwpm | 1.19 [1.15, 1.23] | [1.14, 1.23] | 1 | 1.18 [1.16, 1.21] | -0.010 | not resolved (unresolved_ci) |
| 0.98 | herald_mwpm | 3.84 [3.60, 4.00] | [3.60, 3.87] | 1 | 3.83 [3.74, 3.94] | -0.016 | not resolved (bound_pin) |

CI endpoints can move in the second decimal with the NumPy version.
<!-- END:thresholds -->

![Per-shot threshold panels](figures/threshold_panels.png)

### Herald minus blind

<!-- BEGIN:herald_vs_blind -->
| `r_e` | herald (7,11), % | blind (7,11), % | herald − blind (pp) [95% CI] | replicates missing or saturated | note |
|---|---|---|---|---|---|
| 0 | 0.98 [0.93, 1.05] | 0.98 [0.95, 1.02] | 0.00 [-0.06, 0.07] | 0 / 1000 | control: no heralds |
| 0.5 | 1.60 [1.50, 1.69] | 1.10 [1.03, 1.14] | 0.49 [0.40, 0.60] | 0 / 1000 |  |
| 0.98 | 3.84 [3.60, 4.00] | 1.19 [1.15, 1.23] | 2.65 [2.43, 2.80] | 0 / 1000 |  |

CI endpoints can move in the second decimal with the NumPy version.
<!-- END:herald_vs_blind -->

### Ablation

<!-- BEGIN:ablation -->
| `p`, % | `r_e` | `d` | p_L(blind) / p_L(herald) | herald Λ_d | blind Λ_d |
|---|---|---|---|---|---|
| 0.50 | 0.5 | 3 | 1.3× | 5.11 | 3.29 |
| 0.50 | 0.5 | 5 | 2.1× | 5.63 | 3.02 |
| 0.50 | 0.5 | 7 | 3.9× | — | 3.00 |
| 0.50 | 0.5 | 9 | > 3× † | — | 2.42 |
| 0.50 | 0.5 | 11 | > 4× † | — | — |
| 0.50 | 0.98 | 3 | 2.1× | 21.16 | 3.62 |
| 0.50 | 0.98 | 5 | 12.5× | — | 3.53 |
| 0.50 | 0.98 | 7 | > 16× † | — | 3.31 |
| 0.50 | 0.98 | 9 | > 17× † | — | 2.75 |
| 0.50 | 0.98 | 11 | > 10× † | — | — |
| 1.00 | 0.5 | 3 | 1.3× | 2.48 | 1.72 |
| 1.00 | 0.5 | 5 | 1.9× | 2.26 | 1.54 |
| 1.00 | 0.5 | 7 | 2.8× | 1.99 | 1.35 |
| 1.00 | 0.5 | 9 | 4.1× | 2.00 | 1.30 |
| 1.00 | 0.5 | 11 | 6.4× | — | — |
| 1.00 | 0.98 | 3 | 2.1× | 8.06 | 1.93 |
| 1.00 | 0.98 | 5 | 8.8× | 7.76 | 1.66 |
| 1.00 | 0.98 | 7 | 41.2× | — | 1.52 |
| 1.00 | 0.98 | 9 | > 103× † | — | 1.44 |
| 1.00 | 0.98 | 11 | > 339× † | — | — |

100,000 shots per cell, seed 0. † fewer than 50 herald errors: >= 95% one-sided lower confidence bound, not a point ratio. Λ_d = per-round p_L(d) / p_L(d+2); — below 50 errors at d or d+2.

CI endpoints can move in the second decimal with the NumPy version.
<!-- END:ablation -->

![Ablation](figures/ablation.png)

## Validation against Wu et al.

<!-- BEGIN:wu_validation -->
| config | `R_e` | `p_m` | herald (9,13), % [95% CI] | sign changes (unsaturated) | Wu et al. | ours − Wu (pp) | Wu inside CI |
|---|---|---|---|---|---|---|---|
| wu_fig3_re000 | 0 | 0.0 | 1.05 [1.03, 1.12] | 1 | 0.937(4)% | 0.11 | False |
| wu_fig3_re098 | 0.98 | 0.0 | 4.72 [4.49, 4.79] | 1 | 4.15(2)% | 0.57 | False |
| wu_fig3b_re050 | 0.5 | 0.0 | 1.72 [1.66, 1.87] | 1 | not tabulated | — | — |
| wu_fig3b_re090 | 0.9 | 0.0 | 3.48 [3.27, 3.55] | 1 | not tabulated | — | — |
| wu_fig3b_re100 | 1 | 0.0 | 5.08 [4.99, 5.24] | 1 | 5.13% | -0.05 | True |
| wu_si_pm005 | 0.98 | 0.005 | 3.97 [3.89, 4.19] | 3 | 3.80(2)% | 0.17 | False |
| wu_si_pmp | 0.98 | p | 2.66 [2.62, 2.71] | 1 | 2.85(1)% | -0.19 | False |

CI endpoints can move in the second decimal with the NumPy version.
<!-- END:wu_validation -->

![Wu validation](figures/wu_validation.png)

Full comparison: [docs/VALIDATION.md](docs/VALIDATION.md).

## Model

### Noise conventions

### Metric

## Limitations

## Reproduce

## Layout

## References

## Citation

## License

MIT — see [LICENSE](LICENSE).
