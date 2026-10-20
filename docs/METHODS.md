# Methods

This document states how every result in [README.md](../README.md) and [VALIDATION.md](VALIDATION.md) is computed.
Symbols: `d` is the code distance, `T` the number of syndrome rounds, and `p` the two-qubit-gate error probability.
`r_e` (repo convention) and `R_e` (`wu2022`) are the erasure fractions of `p`.
`p_e` and `p_p` are the `wu2022` erasure and Pauli probabilities, and `p_m` its ancilla flip probability.
`P_L` is the per-shot logical failure probability, `p_L` the per-round rate, and `Λ` a per-round suppression factor.
References in brackets are listed at the [end](#references).

## Circuit and schedule

The circuit is a rotated surface-code memory [1, 2] in the Z basis, with a `d × d` data-qubit array and `d² − 1` ancillas.
The low-distance rotated layout follows Tomita & Svore [3].
`src/erasure_qec/circuits/builder.py` emits every circuit, and every run uses `T = d` rounds.
Each round has 8 ticks: ancilla reset, H on X-ancillas, four CX layers, H on X-ancillas, and `MR` on all ancillas.
After the last round the data qubits are measured, the time-closing Z detectors are formed, and one logical observable is declared.
The observable is one `OBSERVABLE_INCLUDE` over `logical_z_support(d)`.

The CX order is hook-safe: X-ancillas visit their neighbours in the order NE, NW, SE, SW, and Z-ancillas in NE, SE, NW, SW.
`shortest_graphlike_error()` checks this from the builder, with no Monte Carlo (`tests/test_distance_invariants.py`, `figures/hook_regression.png`).
With the hook-safe schedule the graphlike distance equals `d`.
With a deliberately broken schedule it drops to `⌈(d+1)/2⌉`.

## Noise conventions

### Constant non-identity budget

`src/erasure_qec/noise/model.py` maps `NoiseParams` to channel rates, and `src/erasure_qec/noise/injector.py` places them.
Per two-qubit gate the residual Pauli is `DEPOLARIZE2(p·(1 − r_e))`.
Each of the two qubits independently receives `HERALDED_ERASE((2/3)·p·r_e)`, which appends one herald bit.
`HERALDED_ERASE` replaces the qubit with `I/2`, so it is a non-identity error with probability ¾.
To first order in `p` the per-gate non-identity probability is then `p` for every `r_e`.
The exact probability is `p` at `r_e = 0` and `p − p²/4` at `r_e = 1` (`tests/test_noise_model.py::test_error_budget_invariance`).
`X_ERROR(p)` precedes every `M` and `MR`, follows every `R`, and `DEPOLARIZE1(p)` acts on qubits idle in a CX layer.
With `convert_idle` the idle budget is split the same way, as `DEPOLARIZE1(p·(1 − r_e))` plus `HERALDED_ERASE((4/3)·p·r_e)`.
The committed sweeps use `convert_idle = False`.

### wu2022

`src/erasure_qec/noise/wu.py` implements the model of Wu et al. [4], Sec. III and Methods, Sec. VII.A ("Error correcting code simulations").
Each two-qubit gate has exactly one outcome.
With probability `p_p = p·(1 − R_e)` it takes a Pauli uniform over the 15 non-identity two-qubit Paulis.
With probability `p_e = p·R_e` both qubits are erased: one herald fires and each qubit is left in `I/2`.
Otherwise nothing happens.
Ancilla preparation and measurement flip with probability `p_m`, which acts on ancillas only.
There is no idle noise and no single-qubit gate noise.
Data-qubit preparation is noiseless and the final data readout is perfect; Wu et al. do not fix either.

The sampling circuit makes Pauli and erasure exclusive with one `CORRELATED_ERROR` / `ELSE_CORRELATED_ERROR` chain per gate.
Its herald is a flag qubit measured once per gate pair.
The decoding circuit has the same detector layout.
In it each erasure location is four independent tagged generators `X_ERROR` and `Z_ERROR` on each qubit, each at `p_e / 2`.
The decoding circuit treats Pauli and erasure as independent channels, which is an approximation of order `p_p·p_e`.

## DEM partition

### Single-qubit components

`src/erasure_qec/decoding/dem_partition.py` splits the detector error model (DEM) into a Pauli sub-DEM and a herald table.
It consumes `detector_error_model(decompose_errors=True, flatten_loops=True, approximate_disjoint_errors=True)`, flattened.
Every herald detector carries a fourth sentinel coordinate, `(x, y, t, 1)`, and the partition classifies detectors by it.
A mechanism with no herald goes unchanged into the Pauli sub-DEM.
A mechanism with one herald is split at its `^` separators into graphlike components.
Each non-herald component becomes one conditioned edge under that herald.
A mechanism with two or more heralds raises an error; the noise model cannot produce one.
The partition was verified by hand on the `d = 3`, `T = 2` probes and reconciled against Stim ([dem_worksheet.md](dem_worksheet.md)).
The reconciliation includes a case where two erased data qubits fire only one detector.

### Tagged partition (wu2022)

`src/erasure_qec/decoding/tagged_partition.py` partitions the `wu2022` decoding DEM.
Each erasure generator carries a tag `wu_e{k}:{Xa|Za|Xb|Zb}`, where `k` is the herald index.
A tagged component becomes a conditioned edge under herald `k`, so the herald detector never appears in the mechanism.
An untagged mechanism that touches a herald raises an error.

## Decoders

`herald_mwpm` and `wu_herald_mwpm` run MWPM (PyMatching [5]) on the Pauli sub-DEM plus the conditioned edges of the heralds that fired in that shot.
Conditional on its herald, a conditioned edge has probability ½, so its weight `ln((1 − ½)/½)` is 0.
Herald-free shots decode in a batch on one static matcher.
Heralded shots are grouped by their fired-herald signature, with one matcher built per signature.
`blind_mwpm` and `wu_blind_mwpm` zero the herald columns and decode on the static full-DEM matcher.
The blind decoders therefore price each erasure at its unconditional marginal probability.

`tests/test_herald_decoder.py::test_forced_erasure_pattern_blind_fails_herald_corrects` checks a forced-erasure pattern.
Shots are sampled with probability-1 probe erasures on `(1,3)` and `(1,5)` (`probe_q = 1.0`).
Both decoders are compiled from a circuit with `probe_q = 0.02`, a deliberate mis-specification.
It makes the blind decoder price those two locations at a small unconditional marginal, as a decoder unaware of the erasures would.
It costs the herald decoder nothing, because the heralds set those edges to weight 0 in every shot.
The test asserts that herald fails 0 of 2048 shots and blind fails between 15% and 35% of them.
No exact blind count is quoted, because stim's seeded sampling reproduces only on the same machine and stim version.

## Metrics

### Per-shot and per-round failure

`P_L = errors / shots` is the logical failure probability of one shot after `T = d` rounds.
It is the headline metric.
The per-round rate is `p_L = ½·[1 − (1 − 2·P_L)^(1/T)]`, which inverts `P_L = ½·[1 − (1 − 2·p_L)^T]` (`src/erasure_qec/analysis/statistics.py`).
It is a secondary metric.

### Ordering of per-shot and per-round crossings

At `T = d` the per-round crossing of two distances never lies below their per-shot crossing.
`tests/test_crossing.py::test_round_negative_wherever_shot_nonpositive` proves this from the conversion and checks it on synthetic curves.
Which crossing is closer to the asymptotic threshold at finite `d` is an empirical question.
Neither metric is called biased or wrong here.

### Lambda

`Λ_d = p_L(d) / p_L(d + 2)` at fixed `p` is a per-round quantity, and every appearance labels it so.
`Λ > 1` means distance suppresses the per-round rate, the below-threshold behaviour reviewed by Fowler et al. [6].
The `Λ` notation is not from that paper; the citation supports the suppression statement only.
`lambda_factor` bootstraps `Λ` with binomial resampling at both distances.
It excludes replicates whose `d + 2` draw has zero errors rather than propagating `+∞`, and counts them as `n_boot_dropped`.
The exclusion can only lower the upper CI endpoint.
The `lambda_p` values in `experiments/configs/erasure_r50.yaml` (1.5%) and `erasure_r98.yaml` (2.5%) lie above the per-shot blind (7,11) thresholds in `results/headline.json` (1.10% [1.03, 1.14] at `r_e = 0.5`, 1.19% [1.15, 1.23] at `r_e = 0.98`), and no committed result uses `experiments/collect_lambda_scan.py` output, since no `data/*_lambda.csv` is committed and every `Λ` here comes from the threshold sweeps or `results/ablation.json`.

## Crossing estimator

### Pairwise crossing and bracket

`src/erasure_qec/analysis/crossing.py` defines `g(p) = ln P_L(p, d2) − ln P_L(p, d1)` for `d2 > d1`.
The crossing is the first sign change of `g` from negative to non-negative on the grid.
It is located by linear interpolation of `g` in `ln p` between the two bracketing grid points.
A grid point with zero errors at either distance is skipped.
The crossing depends on the distance pair; the (5,7), (7,9) and (9,11) values are in `results/headline.json` (`drift`).
The headline pair is (7,11), and the collapse fit uses `d ≥ 7`, so `d ∈ {3, 5}` enter neither.

### Bootstrap

Each CI is a 95% percentile interval of a parametric bootstrap with `n_boot = 1000` and `seed = 0`.
Each replicate draws `errors ~ Binomial(shots, P_L)` independently at every point and recomputes the estimate.
The interval covers binomial noise only; it excludes interpolation error inside the bracket.
CI endpoints print to two decimals in percent.
Overlap of two marginal 95% CIs is not a 95% test of their difference.
The herald − blind difference is therefore bootstrapped directly (`scripts/headline_results.py`).
Each of its replicates resamples herald, then blind, from one `default_rng(seed)` stream.
A replicate in which either crossing is missing or saturated is counted, not used.
`tests/test_threshold_fit.py::test_bootstrap_ci_covers_known_p_th_at_nominal_rate` checks the collapse-fit CI calibration.
Over 40 seeded synthetic realizations it asserts coverage of the known threshold between 0.80 and 1.0 (per-round metric, `n_boot = 100`).

### Saturated-crossing rule

A crossing whose bracket has `P_L ≥ 0.4` at either distance at either bracket point is reported as "saturated, not a threshold", with no value.
A bootstrap replicate whose crossing is saturated counts as having no crossing.
The cap 0.4 is `SATURATION_CAP` of `threshold_fit.py`.
`estimate_crossing`, the spread-minimising crossing of the collapse fit, is unreliable on saturated data.
It only seeds the collapse-fit window; no quoted threshold is its raw value.

### Sign-change counts

`n_sign_changes` counts sign changes of `g` between consecutive retained grid points.
A value above 1 means the curves cross more than once on the grid; only the first crossing is reported.
`n_sign_changes_unsaturated` counts only the steps with `P_L < 0.4` at both distances at both points.
It flags repeated crossings in the saturated tail, where `g` flips sign on sampling noise.

## Collapse fit

`fit_threshold` (`src/erasure_qec/analysis/threshold_fit.py`) fits `m(p, d) = A + B·x + C·x²` with `x = (p − p_th)·d^(1/ν)`.
Here `m` is `P_L` (metric "shot") or `p_L` (metric "round"), `p_th` is the threshold and `ν` the scaling exponent.
The ansatz follows Wang, Harrington & Preskill [7].
The fit is weighted least squares over the five parameters, with weights from the Wilson half-width of `m`.
It uses only points with `p_center/1.5 ≤ p ≤ 1.5·p_center` and `m < 0.4`, where `p_center` comes from `estimate_crossing`.
The 0.4 cap acts on the fitted metric.
For the per-round fit a near-crossing point at `d = 11`, `T = 11` can have `P_L ≈ 0.45` and still enter the window.
A fit needs at least 3 degrees of freedom, so at least 8 points, and at least 2 distances.
`ν` is bounded to [0.5, 6]; a fit landing on a bound is reported as not converged (`bound_pin`).

A fit is quoted as a threshold only when the code marks it `resolved`.
`resolved` means converged and a 95% CI narrower than `p_th` itself (relative width below 1.0).
In `results/headline.json` every resolved fit has relative width at most 0.243, from the per-round `r_e = 0.5` blind fit.
The per-round `r_e = 0.98` blind fit has relative width 2.37 and is `unresolved_ci`.
The per-round `r_e = 0.98` herald fit pins `ν` at 6 and is `bound_pin`.
At `r_e = 0.98` the per-shot fit resolves both decoders, at 3.83% [3.74, 3.94] herald and 1.18% [1.16, 1.21] blind.
This is an observation; it does not say which metric is correct.

The fit's bootstrap re-runs the whole pipeline per replicate: resampling, `estimate_crossing`, window selection and the weighted fit.
It keeps replicates whose fit pins at a bound, while `bootstrap_threshold_difference` discards them.
The per-round `r_e = 0.98` blind verdict `unresolved_ci` rests on that retention.
`ν` is weakly constrained by three distances, so `p_th` is better determined than `ν`.
Every fit reports `χ²/dof` in `results/headline.json`, so a poor fit is visible.
For caps in {0.2, 0.25, 0.3, 0.4} the per-round fit result was insensitive to the cap, as measured at `da47e3a`.

## Densification and edge rules

The repo-convention grid is 13 log-spaced points in [1e-3, 1e-1] (`experiments/configs/*.yaml`).
`extra_p` adds linearly spaced points around the per-round crossings.
`extra_p_shot` adds points around the per-shot crossings by a rule fixed before the dense data were collected.
For each decoder, with `c` its per-shot (7,11) crossing in `results/threshold_table_committed.json`, the candidates are `round(c·f·10⁴)/10⁴`.
The factors are `f ∈ {0.80, 0.87, 0.93, 1.00, 1.07, 1.15, 1.25}`.
A candidate within 2·10⁻⁴ of an existing point, or of the previous kept candidate, is dropped.
All bands are merged into the grid (union, sorted, de-duplicated).

The `wu2022` grids are either fixed steps or `round(c·f, 4)` with `c` the (5,7) crossing of `scripts/wu_coarse_scan.py`.
If a (9,13) herald bracket touches a grid end, 3 points are added beyond that end, at most once per config.
That edge rule was fixed before any Wu data ([WU_DPAIR_DECISION.md](remediation/WU_DPAIR_DECISION.md#edge-rule-for-scaled-grids-decided-2026-10-16-before-any-wu-data)).

## Statistics gates

### 50-error gate

The ablation and the `Λ` figure use the same constant, 50 observed errors, but not the same rule.
`figure_lambda` requires at least 50 errors at both `d` and `d + 2`.
`scripts/ablation_table.py` gates its ratio on the herald count alone, because herald is the small count there.
Each rule applies the constant to the count that can become small in its statistic.
In `results/ablation.json` blind counts run from 52 to 6290, so both rules keep the same ratio cells.
The two rules are different predicates and would diverge on data where blind is the small count.

### Wilson bounds

A gated ablation cell reports blind's Wilson lower limit over herald's Wilson upper limit, both mapped to per-round rates.
Each limit is an endpoint of a two-sided 95% Wilson interval, so each is a one-sided 97.5% limit.
By Bonferroni the ratio is a ≥ 95% one-sided lower bound on that cell's ratio.
The bound is per cell, unpaired, and nominal at small counts, where Wilson coverage is not guaranteed.

## Heralded fraction

`scripts/heralded_fraction.py` reports the share of DEM error mass that fires at least one herald.
It is a deterministic sum over the DEM, at `d = 5` and `p = 0.02` (`docs/remediation/baseline_outputs/heralded_fraction.txt`).
With gate-only conversion it is 0.23 at `r_e = 0.5` and 0.45 at `r_e = 0.98`.
With `convert_idle` it is 0.32 and 0.63.
The sum excludes no-op herald branches from numerator and denominator.
`HERALDED_ERASE` also heralds on its identity branch, which flips a herald detector and nothing else.
Such a branch causes no syndrome and no logical error, so it is not part of the error budget.
Counting it (`count_no_op_heralds=True`) gives 0.30 and 0.55 gate-only, and 0.42 and 0.72 with `convert_idle`.
That definition overstates conversion.
No threshold, ablation or difference result depends on this quantity.

## Reproducibility

`tests/test_plotting.py::test_figures_regenerate_byte_stable` checks run-to-run determinism within one environment.
The same CSVs rendered twice give byte-identical PNGs.
That is not a cross-environment guarantee.
PNG bytes depend on the matplotlib and FreeType build, and `lambda_vs_p.png` draws bootstrap CIs that move with the NumPy RNG.
On a different matplotlib or NumPy the figures should be visually identical but not byte-identical.

### Commands and wall times

Wall times are given only where a run was recorded; the machine is an Apple M2 with 8 cores.

| step | command | wall time |
|---|---|---|
| install | `uv sync --group dev` | not timed |
| fast tests | `uv run pytest -m "not slow" -q` | not timed |
| headline thresholds | `uv run python scripts/headline_results.py` | not timed |
| pairwise crossings and fits, both metrics | `uv run python scripts/threshold_table.py --out results/threshold_table_dense.json` | not timed |
| ablation | `uv run python scripts/ablation_table.py` | not timed |
| heralded fraction | `uv run python scripts/heralded_fraction.py` | not timed |
| Wu validation table, JSON, figure | `uv run python scripts/wu_validation.py` | 29 s ([VALIDATION.md](VALIDATION.md#reproduce)) |
| figures | `uv run python -m erasure_qec.analysis.plotting --data-dir data --figures-dir figures --ablation-json results/ablation.json` | not timed |
| README tables | `uv run python scripts/render_readme_tables.py` | not timed |

`scripts/ablation_table.py` is the slow step: it decodes 100,000 shots per cell from circuits, with no CSV shortcut.
Its `d = 9` and `d = 11` cells dominate the runtime, and `--shots` trades runtime against the 50-error gate.
It refuses to run on a dirty work tree.

Collection appends to the CSVs and is resumable.
The repo-model sweeps use `uv run python experiments/collect_threshold_sweep.py experiments/configs/<name>.yaml --workers 8` (not timed).
The Wu sweeps use `experiments/collect_wu.py`; their recorded wall times are in [VALIDATION.md](VALIDATION.md#reproduce).

Two legacy scripts are not part of the current results.
`scripts/paired_separation.py` reproduces its recorded output (`docs/remediation/baseline_outputs/paired_separation.txt`) only at commit `da47e3a`.
`scripts/audit_checks.py` reproduces the archived audit table in [archive/AUDIT_HISTORY.md](archive/AUDIT_HISTORY.md), and its expected column is the pre-fix baseline.
Neither reproduces a number published in README.md.

### Seeds

Every analysis bootstrap uses `seed = 0` through `numpy.random.default_rng`.
The ablation samples with seed 0 through stim, whose seeded sampling reproduces only on the same machine and stim version.
`sinter.collect` takes no seed, so the sweeps in `data/` cannot be regenerated shot for shot.

### NumPy-dependent CI endpoints

The analysis is deterministic given the CSV, the seed and the NumPy build.
Point estimates, `ν`, `χ²/dof` and every converged or resolved verdict carry no RNG.
Bootstrap CI endpoints depend on the NumPy RNG and can move in the second decimal with the NumPy version.
`results/*.json` record the software versions they were produced with.

### sinter text pickling

sinter sends circuits and DEMs to its workers by pickle, which stim implements through text with about 6 significant digits.
Sampled probabilities can therefore differ from their nominal values by about 1e-6 relative.
No result depends on this.

## Model assumptions

The physical assumptions (instantaneous perfect heralds, independent per-qubit erasure, unbiased erasure, Z memory) are listed in [README.md](../README.md#limitations).
The analysis assumes that each (decoder, `d`, `p`) task's errors are binomial in its shots.
The `wu2022` decoder assumes independent Pauli and erasure channels, while the sampler makes them exclusive.

## Frozen structural quirks

### Redundant R after MR

Each round starts with `R` on all ancillas, although the previous round's `MR` has already reset them.
Reset noise `X_ERROR(p)` follows that `R`, so each ancilla takes one reset flip and one measurement flip per round.

### Idle noise only in CX layers

`DEPOLARIZE1` idle noise is applied only to qubits idle during one of the four CX layers.
Qubits idle during the reset, H and `MR` ticks receive no idle noise.

### Golden hashes

`tests/test_circuit_golden.py` pins `sha256(str(circuit))` for six repo-convention circuits under stim 1.16.0.
The committed sweeps were collected from these circuits, so the two quirks above are frozen with them.
A failing digest means the circuits changed and the committed data no longer match them.

## References

1. Dennis, Kitaev, Landahl & Preskill, Topological quantum memory, J. Math. Phys. 43, 4452 (2002).
2. Bombin & Martin-Delgado, Optimal resources for topological two-dimensional stabilizer codes: comparative study, Phys. Rev. A 76, 012305 (2007).
3. Tomita & Svore, Low-distance surface codes under realistic quantum noise, Phys. Rev. A 90, 062320 (2014).
4. Wu, Kolkowitz, Puri & Thompson, Erasure conversion for fault-tolerant quantum computing in alkaline earth Rydberg atom arrays, Nat. Commun. 13, 4657 (2022).
5. Higgott & Gidney, Sparse Blossom: correcting a million errors per core second with minimum-weight matching, Quantum 9, 1600 (2025).
6. Fowler, Mariantoni, Martinis & Cleland, Surface codes: Towards practical large-scale quantum computation, Phys. Rev. A 86, 032324 (2012).
7. Wang, Harrington & Preskill, Confinement-Higgs transition in a disordered gauge theory and the accuracy threshold for quantum memory, Ann. Phys. 303, 31 (2003).
