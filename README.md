# Herald-conditioned decoding of the erasure-converted rotated surface code

This repository simulates a distance-`d` rotated surface-code Z memory in which a fraction `r_e` of each two-qubit gate's error probability `p` becomes heralded erasure.
A herald-conditioned matching decoder sets the edges of each heralded erasure to weight 0, and a blind decoder ignores the heralds.
At `r_e = 0.98` the per-shot `d = 7` and `d = 11` curves cross at 3.84% [3.60, 4.00] for herald and 1.19% [1.15, 1.23] for blind (bootstrap 95% CI, `scripts/headline_results.py`).

## Results

Every rate in this section is a logical failure probability; [Metric](#metric) defines `P_L` (per shot) and `p_L` (per round).
Each table below is rendered from `results/*.json` by `scripts/render_readme_tables.py`.
Methods are in [docs/METHODS.md](docs/METHODS.md).

### Thresholds

The headline threshold estimate is the first crossing of the `d = 7` and `d = 11` curves of `P_L` against `p`.
It is located by linear interpolation of `ln P_L(11) − ln P_L(7)` in `ln p` between two grid points, listed as the bracket.
Its 95% CI is a percentile interval of an independent-binomial parametric bootstrap with `n_boot = 1000` and `seed = 0`.
The CI covers binomial noise only and excludes interpolation error inside the bracket.
The values come from `scripts/headline_results.py`, which writes `results/headline.json`.

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

The per-shot collapse column is a finite-size scaling fit over `d ∈ {7, 9, 11}` with the ansatz of Wang et al. [[9]](#references).
It lies below the pairwise crossing in all six conditions, by 0.004 to 0.025 percentage points (pp; gap column).
We report the gap as a systematic difference between two estimators, not as an error of either one.
The per-round column is a secondary metric (see [Metric](#metric)).

At `r_e = 0.5` the blind (7,11) estimate, 1.10% [1.03, 1.14], sits at the lower edge of its bracket [1.10, 1.18].
Its CI extends below the bracket.
The blind (9,11) crossing is 1.03% [1.02, 1.19] (same estimator, `results/headline.json`).
The `r_e = 0.5` blind value should be quoted only with its CI.
Its per-round CI upper endpoint, 1.80%, coincides with the sweep grid point `p = 0.018`.

At `r_e = 0.98` the per-shot collapse fit resolves both thresholds: herald 3.83% [3.74, 3.94] and blind 1.18% [1.16, 1.21].
The per-round collapse fit at `r_e = 0.98` resolves neither threshold.
Its herald fit pins the scaling exponent `ν` at the optimiser bound 6 (`bound_pin`).
Its blind fit has a 95% CI 2.37 times wider than its own fitted threshold `p_th` (`unresolved_ci`, `results/headline.json`).
This is an observation about two fits; it makes no claim about which metric is correct.

![Per-shot threshold panels](figures/threshold_panels.png)

*Per-shot `P_L` against `p` for `d ∈ {3, 5, 7, 9, 11}`, one panel per `r_e`, herald solid and blind dashed.*
*Vertical lines mark each decoder's (7,11) crossing, and shaded bands its 95% bootstrap CI.*
*A decoder whose crossing is missing or saturated carries no line.*

### Herald minus blind

The table gives herald − blind of the per-shot (7,11) crossings, with a 95% CI from an unpaired bootstrap (`n_boot = 1000`, `seed = 0`, `scripts/headline_results.py`).

<!-- BEGIN:herald_vs_blind -->
| `r_e` | herald (7,11), % | blind (7,11), % | herald − blind (pp) [95% CI] | replicates missing or saturated | note |
|---|---|---|---|---|---|
| 0 | 0.98 [0.93, 1.05] | 0.98 [0.95, 1.02] | 0.00 [-0.06, 0.07] | 0 / 1000 | control: no heralds |
| 0.5 | 1.60 [1.50, 1.69] | 1.10 [1.03, 1.14] | 0.49 [0.40, 0.60] | 0 / 1000 |  |
| 0.98 | 3.84 [3.60, 4.00] | 1.19 [1.15, 1.23] | 2.65 [2.43, 2.80] | 0 / 1000 |  |

CI endpoints can move in the second decimal with the NumPy version.
<!-- END:herald_vs_blind -->

At `r_e = 0` there are no heralds, and the difference, 0.00 pp [−0.06, 0.07], is consistent with zero as a control requires.
At `r_e = 0.5` the difference is 0.49 pp [0.40, 0.60].
At `r_e = 0.98` the difference is 2.65 pp [2.43, 2.80].
No bootstrap replicate lacked a crossing or had a saturated one, at any `r_e`.

The sweeps sample each decoder on its own shots, so the bootstrap resamples the two decoders independently.
A paired CI would need both decoders on shared shots.
With positively correlated failures on shared shots a paired CI would be narrower, so this CI is conservative relative to that paired CI.
We make no causal claim about why the blind crossing moves less with `r_e` than the herald crossing.

### Ablation

`scripts/ablation_table.py` samples each cell's shots once and decodes the same detection events with both decoders.
The ablation is therefore paired by construction and needs no threshold fit.
It does not depend on how the noise budget is normalised across `r_e`.
Each cell has 100,000 shots with seed 0, at `p = 0.5%` and `p = 1.0%` (`results/ablation.json`).
The ratio column is the per-round `p_L(blind) / p_L(herald)`.
`Λ_d` is one decoder's per-round `p_L(d) / p_L(d+2)`, the per-round suppression factor (see [docs/METHODS.md](docs/METHODS.md#lambda)).

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

The ratio is quoted as a point value only where herald has at least 50 errors.
Below that gate (†) the cell shows a lower bound, not a point ratio.
The bound is blind's Wilson lower limit divided by herald's Wilson upper limit, both mapped through the per-round conversion.
Each limit is an endpoint of a two-sided 95% Wilson interval, which is a one-sided 97.5% limit.
By Bonferroni the quoted value is then a ≥ 95% one-sided lower bound on that cell's ratio.
The bound is per cell and not simultaneous across the table.
It ignores the positive correlation from decoding identical shots, which only makes it looser.
At single-digit counts the Wilson "95%" is a nominal approximation with no guaranteed coverage.

The gated cells rest on few herald errors (`h_err` in `results/ablation.json`).
At `p = 1.0%` and `r_e = 0.98` the `d = 9` cell has 23 herald errors and the `d = 11` cell has 3.
**The `d = 11` bound > 339× rests on 3 herald errors and is weak evidence.**
It reads as "herald errors were too rare to measure at `d = 11`", not as "suppression is 339×".
At `p = 0.5%` and `r_e = 0.98` the `d = 7`, 9 and 11 cells have 9, 1 and 0 herald errors.
At `p = 0.5%` and `r_e = 0.5` the `d = 9` and 11 cells have 42 and 9 herald errors.

For each pair of `p` and `r_e`, the ungated ratios increase with `d`.
In every cell that reports both, herald `Λ_d` exceeds blind `Λ_d`.
How far the increase continues in the gated cells is bounded from below by this run but not measured by it.
Measuring it needs more shots per cell, not larger distances.

![Per-round Lambda against p](figures/lambda_vs_p.png)

*Per-round `Λ = p_L(d) / p_L(d+2)` against `p` from the threshold sweeps, with bootstrap CIs.*
*A point is shown only where both `d` and `d + 2` have at least 50 observed errors.*

## Validation against Wu et al.

The `wu2022` convention reproduces the noise model of Wu et al. [[6]](#references), defined in their Methods, Sec. VII.A ("Error correcting code simulations").
`R_e` is the erasure fraction of each two-qubit gate's error probability `p`.
`p_m` is the ancilla preparation and measurement flip probability.
The table gives the herald crossing of the `d = 9` and `d = 13` per-shot curves (`scripts/wu_validation.py`, `results/wu_validation.json`).

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

Wu's value lies inside our 95% CI for `R_e = 1` only, and outside it for the four other tabulated configs.
The `d` pair (9,13) was set by a pre-registered memory rule ([WU_DPAIR_DECISION.md](docs/remediation/WU_DPAIR_DECISION.md)).
Wu et al. used the (11,15) pair.
Wu simulate the planar XZZX code [[10]](#references); we simulate the rotated CSS code.
Wu state that the XZZX code "performs identically to the standard surface code for the case of unbiased noise".
Both erasure models here are unbiased, and Wu's `d × d` data-qubit array has the same count as our rotated patch.
Wu decode with weighted Union-Find, which they note "is optimal for pure erasure errors" [[11]](#references).
We decode with MWPM and set heralded edges to weight 0.
The full comparison and the other differences are in [docs/VALIDATION.md](docs/VALIDATION.md).

The repo-model thresholds under [Results](#results) are not on a comparable footing with Wu's.
The repo model adds measurement, reset and idle noise that Wu's model does not have.
So no numerical agreement with Wu is claimed for them at any `r_e`, including the `r_e = 0` baseline.
The two models agree only structurally: in both, the herald threshold rises with the converted fraction.

## Model

The circuit is a rotated surface-code memory [[1,2,3]](#references) in the Z basis, built in Stim [[4]](#references).
It runs `T = d` rounds of syndrome extraction with circuit-level noise at gate, preparation, measurement and idle locations.
Decoding uses PyMatching [[5]](#references) through `sinter`.
Circuit, schedule and detector details are in [docs/METHODS.md](docs/METHODS.md#circuit-and-schedule).

### Noise conventions

The repo convention (`data/*.csv`) holds the non-identity error budget of each two-qubit gate constant across `r_e`.
A fraction `r_e` of that budget becomes a heralded erasure: the hardware flags which qubit was disturbed and when.
`HERALDED_ERASE` replaces the qubit with `I/2`, so all four Paulis are equally likely and the erasure is unbiased.
Per two-qubit gate on `(a, b)`, and at the other noise locations:

| Channel | Rate | Where |
|---|---|---|
| `DEPOLARIZE2` | `p·(1 − r_e)` | on `(a, b)`, residual Pauli |
| `HERALDED_ERASE` | `(2/3)·p·r_e` | on `a`, and independently on `b`; appends a herald bit |
| `X_ERROR` / `X_ERROR` / `DEPOLARIZE1` | `p_meas` / `p_reset` / `p_idle` (each `= p`) | before `M`/`MR`, after `R`, on idle qubits |

An erasure is a non-identity error only ¾ of the time, so `2·q·¾ = p·r_e` gives `q = (2/3)·p·r_e`.
To first order in `p` this holds the per-gate non-identity probability at `p` for every `r_e` (`tests/test_noise_model.py::test_error_budget_invariance`).
So `r_e` converts the budget rather than shrinking it, and `p` is an iso-noise axis.

`r_e` is the two-qubit-gate fraction, not the circuit-wide heralded fraction.
Measurement and reset errors are never converted, and idle errors only with `convert_idle`.
The committed sweeps convert two-qubit-gate errors only.
At `r_e = 0.98` the heralded share of the error mass is 0.45, not 0.98 (exact DEM sum, `d = 5`, `p = 0.02`, `scripts/heralded_fraction.py`).

The `wu2022` convention (`data/wu/*.csv`) follows Wu et al. [[6]](#references).
Each two-qubit gate has exactly one of three outcomes.
With probability `p_p = p·(1 − R_e)` it takes a Pauli uniform over the 15 non-identity two-qubit Paulis.
With probability `p_e = p·R_e` it erases both qubits, leaving each in `I/2` with one herald.
Otherwise it has no error.
Ancilla preparation and measurement flip with probability `p_m`, and there is no idle noise.
Wu's Figs. 3–4 set `p_m = 0` and "neglect idle errors", so two-qubit gates are their only error locations.
Their `R_e` is therefore also their circuit-wide erasure fraction, whereas the repo convention's `r_e` is not.
The repo convention thus has noise locations that Wu's model lacks, so its `r_e = 0` baseline is not compared with Wu's 0.937%.

### Metric

The headline metric is `P_L`, the logical failure probability of one shot after `T = d` rounds.
The per-round rate `p_L = ½·[1 − (1 − 2·P_L)^(1/T)]` is a secondary metric.
At `T = d` the per-round crossing of two distances never lies below their per-shot crossing (`tests/test_crossing.py::test_round_negative_wherever_shot_nonpositive`).
Which of the two is closer to the asymptotic threshold at finite `d` is an empirical question this repository does not settle.
Thresholds are per-shot, while `Λ` is per-round; each appearance names its metric.

## Limitations

This is an idealized erasure model.
Each assumption below differs from hardware.
The size of each effect on the results was not measured.

- Erasure is injected after the CX, independently per qubit.
  There is no mid-gate erasure that propagates through the rest of the gate.
  There is no correlated two-qubit erasure in the repo convention: the two qubits of a gate erase independently.
  Real leakage happens during the gate and can corrupt both qubits together.
- The erased qubit keeps participating in all later rounds.
  `HERALDED_ERASE` replaces the qubit with `I/2` for that instant only, and the qubit works normally from the next tick.
  A real ¹⁷¹Yb erasure [[7]](#references) is an atom leaked out of the qubit subspace, which stays broken until it is re-prepared.
  A single erasure could then degrade many later rounds, which is not modelled.
- Heralds are perfect: instantaneous, 100% reliable, with no false positives and no latency.
  Real fluorescence heralding has finite fidelity, false positives and negatives, and a latency the decoder must tolerate.
- The erasure is unbiased: `HERALDED_ERASE` is `I/2`, with all four Paulis equally likely.
  The Z-biased erasure of Sahay et al. [[8]](#references) would give a further threshold advantage that is not modelled.
  The advantage reported here comes from herald-conditioned decoding alone.
- The memory is Z-basis only, with a single logical observable.
  The observable is one `OBSERVABLE_INCLUDE` over `logical_z_support(d)`, and only Z-type detectors close at the time boundaries.
  Every threshold here is a Z-basis memory threshold under that observable.
  It is not a full X and Z circuit-level threshold, and not a logical-gate or computation threshold.
  Every number in this repository carries the "Z-memory" qualifier.
- The decoder is MWPM with heralded edges at weight 0, which is not maximum-likelihood for erasure.
  The gap to a maximum-likelihood erasure decoder, and whether it grows with `r_e`, was not measured.
  A peeling or Union-Find decoder [[11,12]](#references) would measure it.

Future directions are in [docs/FUTURE_WORK.md](docs/FUTURE_WORK.md).
They include imperfect and delayed heralds, other hardware, logical algorithms and alternative decoders.

## Reproduce

All results and figures are computed from the committed sweeps `data/baseline_pauli.csv`, `data/erasure_r50.csv` and `data/erasure_r98.csv`, and from `data/wu/*.csv`.
`sinter.collect` takes no seed, so the sweeps cannot be regenerated shot for shot ([docs/METHODS.md](docs/METHODS.md#seeds)).
Wall times are given only where a run was recorded; the machine is an Apple M2 with 8 cores.

| step | command | wall time |
|---|---|---|
| install | `uv sync --group dev` | not timed |
| fast tests | `uv run pytest -m "not slow" -q` | not timed |
| headline thresholds | `uv run python scripts/headline_results.py` | not timed |
| pairwise crossings and fits, both metrics | `uv run python scripts/threshold_table.py --out results/threshold_table_dense.json` | not timed |
| ablation | `uv run python scripts/ablation_table.py` | not timed |
| heralded fraction | `uv run python scripts/heralded_fraction.py` | not timed |
| Wu validation table, JSON, figure | `uv run python scripts/wu_validation.py` | 29 s ([docs/VALIDATION.md](docs/VALIDATION.md#reproduce)) |
| figures | `uv run python -m erasure_qec.analysis.plotting --data-dir data --figures-dir figures --ablation-json results/ablation.json` | not timed |
| README tables | `uv run python scripts/render_readme_tables.py` | not timed |

`scripts/ablation_table.py` is the slow step: it decodes 100,000 shots per cell from circuits, with no CSV shortcut.
Its `d = 9` and `d = 11` cells dominate the runtime, and `--shots` trades runtime against the 50-error gate.
It refuses to run on a dirty work tree.

Collection appends to the CSVs and is resumable.
The repo-model sweeps use `uv run python experiments/collect_threshold_sweep.py experiments/configs/<name>.yaml --workers 8` (not timed).
The Wu sweeps use `experiments/collect_wu.py`; their recorded wall times are in [docs/VALIDATION.md](docs/VALIDATION.md#reproduce).

Two legacy scripts are not part of the current results.
`scripts/paired_separation.py` reproduces its recorded output (`docs/remediation/baseline_outputs/paired_separation.txt`) only at commit `da47e3a`.
`scripts/audit_checks.py` reproduces the archived audit table in [docs/archive/AUDIT_HISTORY.md](docs/archive/AUDIT_HISTORY.md), and its expected column is the pre-fix baseline.
Neither reproduces a number published in this README.

## Layout

```
src/erasure_qec/
├── config.py        # NoiseParams, ExperimentConfig, YAML loading
├── circuits/        # layout.py (geometry), scheduling.py (CX schedules), builder.py (emits stim)
├── noise/           # model.py (channel rates), injector.py (repo convention), wu.py (wu2022)
├── decoding/        # dem_partition.py, tagged_partition.py, herald_matching.py, sinter_adapter.py
└── analysis/        # statistics.py, crossing.py, threshold_fit.py, dem_stats.py, synthetic.py, plotting.py
experiments/         # collection drivers and configs/ (repo convention) and configs/wu/ (wu2022)
scripts/             # analyses that write results/*.json, the README table renderer, dev probes
results/             # JSON outputs with git describe, software versions and data sha256
data/                # committed sweeps: *.csv (repo convention), wu/*.csv (wu2022)
figures/             # committed figures
docs/                # METHODS.md, VALIDATION.md, FUTURE_WORK.md, dem_worksheet.md, remediation/, archive/
tests/               # pytest suite; fixtures/ holds pinned CSVs
```

Pinned regression snapshots of real sweeps live at `tests/fixtures/real_*.csv`.
The committed synthetic fixtures `tests/fixtures/synthetic_*.csv` drive the analysis and plotting tests and are not measurements.
They are generated from the same collapse variable the fitter inverts, and the fitter must reject the adversarial one.
The DEM hand-verification worksheet is [docs/dem_worksheet.md](docs/dem_worksheet.md).
History up to `da47e3a` is in [docs/archive/](docs/archive/README.md), and changes are in [CHANGELOG.md](CHANGELOG.md).

## References

1. Dennis, Kitaev, Landahl & Preskill, *Topological quantum memory*, J. Math. Phys. 43, 4452 (2002), [arXiv:quant-ph/0110143](https://arxiv.org/abs/quant-ph/0110143). The surface-code memory.
2. Bombin & Martin-Delgado, *Optimal resources for topological two-dimensional stabilizer codes: comparative study*, Phys. Rev. A 76, 012305 (2007), [arXiv:quant-ph/0703272](https://arxiv.org/abs/quant-ph/0703272). The rotated layout.
3. Tomita & Svore, *Low-distance surface codes under realistic quantum noise*, Phys. Rev. A 90, 062320 (2014), [arXiv:1404.3747](https://arxiv.org/abs/1404.3747). The low-distance rotated layout.
4. Gidney, *Stim: a fast stabilizer circuit simulator*, Quantum 5, 497 (2021), [arXiv:2103.02202](https://arxiv.org/abs/2103.02202).
5. Higgott & Gidney, *Sparse Blossom: correcting a million errors per core second with minimum-weight matching*, Quantum 9, 1600 (2025), [arXiv:2303.15933](https://arxiv.org/abs/2303.15933).
6. Wu, Kolkowitz, Puri & Thompson, *Erasure conversion for fault-tolerant quantum computing in alkaline earth Rydberg atom arrays*, Nat. Commun. 13, 4657 (2022), [arXiv:2201.03540](https://arxiv.org/abs/2201.03540). Noise model: Sec. III and Methods, Sec. VII.A.
7. Ma et al., *High-fidelity gates and mid-circuit erasure conversion in an atomic qubit*, Nature 622, 279 (2023).
8. Sahay, Jin, Claes, Thompson & Puri, *High-threshold codes for neutral-atom qubits with biased erasure errors*, Phys. Rev. X 13, 041013 (2023), [arXiv:2302.03063](https://arxiv.org/abs/2302.03063).
9. Wang, Harrington & Preskill, *Confinement-Higgs transition in a disordered gauge theory and the accuracy threshold for quantum memory*, Ann. Phys. 303, 31 (2003), [arXiv:quant-ph/0207088](https://arxiv.org/abs/quant-ph/0207088).
10. Bonilla Ataides, Tuckett, Bartlett, Flammia & Brown, *The XZZX surface code*, Nat. Commun. 12, 2172 (2021), [arXiv:2009.07851](https://arxiv.org/abs/2009.07851).
11. Delfosse & Zémor, *Linear-time maximum likelihood decoding of surface codes over the quantum erasure channel*, Phys. Rev. Research 2, 033042 (2020), [arXiv:1703.01517](https://arxiv.org/abs/1703.01517).
12. Delfosse & Nickerson, *Almost-linear time decoding algorithm for topological codes*, Quantum 5, 595 (2021), [arXiv:1709.06218](https://arxiv.org/abs/1709.06218).

## Citation

Cite the repository and the commit you used:

```bibtex
@software{erasure_qec,
  author = {Naman},
  title  = {Herald-conditioned decoding of the erasure-converted rotated surface code},
  url    = {https://github.com/albnewt97/erasure-qec},
  note   = {version 0.2.0, commit <hash>},
  year   = {2026}
}
```

## License

MIT — see [LICENSE](LICENSE).
