# Validation against Wu et al. (2022) — draft

Reference: Wu, Kolkowitz, Puri & Thompson, *Erasure conversion for fault-tolerant quantum
computing in alkaline earth Rydberg atom arrays*, Nat. Commun. 13, 4657 (2022). Every number below is printed by `scripts/wu_validation.py` and stored in
`results/wu_validation.json`, unless another source is named.

## What was simulated

The noise model follows Wu et al. (convention `wu2022`, `src/erasure_qec/noise/wu.py`).
- Each two-qubit gate has exactly one of three outcomes: a Pauli uniform over the 15 non-identity two-qubit Paulis with probability p(1 − R_e); an erasure of both qubits (each left in I/2, heralded) with probability p·R_e; or no error.
- Ancilla initialisation and measurement each flip with probability p_m; p_m acts on ancillas only.
- There is no idle noise and no single-qubit gate noise.
- Data-qubit initialisation is noiseless and the final data readout is perfect.

The circuit is a rotated CSS surface-code Z memory with T = d rounds and d ∈ {7, 9, 13}. Each (d, p) task ran to 100 000 shots or 1000 logical failures, whichever came first. The seven configs (`experiments/configs/wu/*.yaml`) are re000 (R_e = 0), re050 (0.5), re090 (0.9), re098 (0.98) and re100 (1.0), all with p_m = 0, plus pm005 (R_e = 0.98, p_m = 0.005) and pmp (R_e = 0.98, p_m = p).

The re050 grid was extended once by the pre-registered edge rule ([WU_DPAIR_DECISION.md](remediation/WU_DPAIR_DECISION.md#edge-rule-for-scaled-grids-decided-2026-10-16-before-any-wu-data)). On its 11-point base grid the (9,13) herald curves did not cross, and g = ln P_L(13) − ln P_L(9) was +0.044 at the lowest point, p = 1.76% (the `g` list in the JSON). That is the rule's lower-edge case, so three points below the grid were added: p = 1.46%, 1.56% and 1.66% (commits 0ee7559 and 2668f5e). On the extended grid the bracket is interior, and no other config's herald bracket touches an edge.

## How it was decoded

`wu_herald_mwpm` runs MWPM (PyMatching 2.4) on the decoding circuit's error model, with the edges of every heralded erasure set to weight 0 in that shot. `wu_blind_mwpm` runs MWPM on the same model with the heralds ignored. The comparison with Wu et al. uses `wu_herald_mwpm`.

The headline estimate is the first crossing of the per-shot P_L curves of d = 9 and d = 13. It is the first g < 0 → g ≥ 0 sign change, interpolated linearly in (ln p, g).
Its 95% CI is a percentile interval from a parametric bootstrap with independent binomial resampling (n_boot = 1000, seed = 0). The CI covers binomial noise only; the interpolation error inside the bracket is not included. The (7,9) crossing is a cross-check. `fit_threshold(metric="shot", d_min=7)` gives a secondary estimate. When a curve pair has no crossing on the grid, the table gives a status in place of a value. `n_sign_changes_unsaturated` counts the sign changes among the p at which P_L_shot < 0.4 at both distances. "Wu inside CI" tests Wu's central value against our 95% CI and does not use Wu's stated uncertainty.

## Results

Per-shot pairwise crossings, % (95% bootstrap percentile CI, n_boot = 1000, seed = 0); sc = n_sign_changes, sc_u = n_sign_changes_unsaturated, nx = replicates without a crossing.

| config | R_e | p_m | decoder | (9,13) | bracket | sc / sc_u / nx | (7,9) | bracket | sc / sc_u / nx |
|---|---|---|---|---|---|---|---|---|---|
| wu_fig3_re000 | 0.0 | 0.0 | wu_herald_mwpm | 1.049 [1.034, 1.122] | [1.000, 1.050] | 1 / 1 / 0 | 1.095 [1.036, 1.166] | [1.050, 1.100] | 1 / 1 / 0 |
| wu_fig3b_re050 | 0.5 | 0.0 | wu_herald_mwpm | 1.720 [1.655, 1.866] | [1.660, 1.760] | 1 / 1 / 0 | 1.712 [1.556, 1.779] | [1.660, 1.760] | 1 / 1 / 0 |
| wu_fig3b_re050 | 0.5 | 0.0 | wu_blind_mwpm | no crossing in grid (below) | — | 0 / 0 / 1000 | no crossing in grid (below) | — | 0 / 0 / 1000 |
| wu_fig3b_re090 | 0.9 | 0.0 | wu_herald_mwpm | 3.481 [3.267, 3.550] | [3.440, 3.580] | 1 / 1 / 0 | 3.372 [3.303, 3.489] | [3.290, 3.440] | 1 / 1 / 0 |
| wu_fig3b_re090 | 0.9 | 0.0 | wu_blind_mwpm | no crossing in grid (below) | — | 0 / 0 / 705 | no crossing in grid (below) | — | 0 / 0 / 999 |
| wu_fig3_re098 | 0.98 | 0.0 | wu_herald_mwpm | 4.717 [4.486, 4.794] | [4.700, 4.800] | 1 / 1 / 0 | 4.380 [4.169, 4.642] | [4.300, 4.400] | 3 / 3 / 0 |
| wu_fig3_re098 | 0.98 | 0.0 | wu_blind_mwpm | 3.505 [3.419, 4.600] | [3.500, 3.600] | 6 / 0 / 9 | no crossing in grid (below) | — | 0 / 0 / 143 |
| wu_fig3b_re100 | 1.0 | 0.0 | wu_herald_mwpm | 5.076 [4.993, 5.239] | [5.000, 5.100] | 1 / 1 / 0 | 4.936 [4.750, 5.032] | [4.900, 5.000] | 1 / 1 / 0 |
| wu_fig3b_re100 | 1.0 | 0.0 | wu_blind_mwpm | 4.810 [4.309, 5.413] | [4.800, 4.900] | 5 / 0 / 0 | 5.839 [4.211, 5.600] | [5.800, 5.900] | 2 / 0 / 0 |
| wu_si_pm005 | 0.98 | 0.005 | wu_herald_mwpm | 3.971 [3.891, 4.193] | [3.900, 4.000] | 3 / 3 / 0 | 4.050 [3.588, 4.109] | [4.000, 4.100] | 1 / 1 / 0 |
| wu_si_pm005 | 0.98 | 0.005 | wu_blind_mwpm | 3.900 [3.313, 4.368] | [3.800, 3.900] | 4 / 0 / 20 | 4.435 [3.509, 4.570] | [4.400, 4.500] | 2 / 0 / 225 |
| wu_si_pmp | 0.98 | p | wu_herald_mwpm | 2.656 [2.617, 2.707] | [2.600, 2.700] | 1 / 1 / 0 | 2.496 [2.440, 2.766] | [2.400, 2.500] | 3 / 3 / 0 |
| wu_si_pmp | 0.98 | p | wu_blind_mwpm | no crossing in grid (below) | — | 0 / 0 / 312 | no crossing in grid (below) | — | 0 / 0 / 606 |

fit_threshold(metric = shot, d_min = 7), % (95% bootstrap percentile CI, n_boot = 1000); gap = collapse − (9,13) pairwise, pp.

| config | decoder | p_th | ν | χ²/dof | converged | resolved | reason | gap (pp) |
|---|---|---|---|---|---|---|---|---|
| wu_fig3_re000 | wu_herald_mwpm | 1.091 [1.077, 1.107] | 1.189 | 1.580 | True | True | ok | 0.042 |
| wu_fig3b_re050 | wu_herald_mwpm | 1.720 [1.701, 1.738] | 1.079 | 1.241 | True | True | ok | 0.000 |
| wu_fig3b_re050 | wu_blind_mwpm | 1.460 [1.460, 1.460] | 0.553 | 9.027 | False | False | bound_pin | — |
| wu_fig3b_re090 | wu_herald_mwpm | 3.393 [3.362, 3.427] | 1.175 | 2.296 | True | True | ok | -0.088 |
| wu_fig3b_re090 | wu_blind_mwpm | — | — | — | False | False | insufficient_data | — |
| wu_fig3_re098 | wu_herald_mwpm | 4.608 [4.568, 4.649] | 1.251 | 4.304 | True | True | ok | -0.109 |
| wu_fig3_re098 | wu_blind_mwpm | — | — | — | False | False | insufficient_data | — |
| wu_fig3b_re100 | wu_herald_mwpm | 5.002 [4.979, 5.026] | 0.969 | 1.387 | True | True | ok | -0.073 |
| wu_fig3b_re100 | wu_blind_mwpm | — | — | — | False | False | insufficient_data | — |
| wu_si_pm005 | wu_herald_mwpm | 3.966 [3.934, 4.000] | 1.262 | 3.627 | True | True | ok | -0.005 |
| wu_si_pm005 | wu_blind_mwpm | — | — | — | False | False | insufficient_data | — |
| wu_si_pmp | wu_herald_mwpm | 2.625 [2.611, 2.638] | 0.945 | 8.539 | True | True | ok | -0.031 |
| wu_si_pmp | wu_blind_mwpm | 2.200 [2.200, 2.900] | 4.099 | 24.886 | False | False | bound_pin | — |

Comparison with Wu et al. (wu_herald_mwpm, (9,13) per-shot crossing).

| config | ours, % [95% CI] | Wu et al. | source | ours − Wu (pp) | Wu inside CI | edge rule |
|---|---|---|---|---|---|---|
| wu_fig3_re000 | 1.049 [1.034, 1.122] | 0.937(4)% | Fig. 3a | 0.112 | False | none |
| wu_fig3b_re050 | 1.720 [1.655, 1.866] | not tabulated | not tabulated | — | — | none |
| wu_fig3b_re090 | 3.481 [3.267, 3.550] | not tabulated | not tabulated | — | — | none |
| wu_fig3_re098 | 4.717 [4.486, 4.794] | 4.15(2)% | Fig. 3a | 0.567 | False | none |
| wu_fig3b_re100 | 5.076 [4.993, 5.239] | 5.13% | Sec. III text | -0.054 | True | none |
| wu_si_pm005 | 3.971 [3.891, 4.193] | 3.80(2)% | SI S6 | 0.171 | False | none |
| wu_si_pmp | 2.656 [2.617, 2.707] | 2.85(1)% | SI S6 | -0.194 | False | none |

Gap to Wu per config, from the (9,13) herald crossing and its 95% bootstrap CI:
- re000: ours − Wu = +0.112 pp, and Wu's 0.937% lies outside our CI [1.034, 1.122]%. Not explained.
- re098: ours − Wu = +0.567 pp, and Wu's 4.15% lies outside our CI [4.486, 4.794]%. Not explained.
- re100: ours − Wu = −0.054 pp, and Wu's 5.13% lies inside our CI [4.993, 5.239]%. Our value agrees with Wu's to within our 95% CI. Wu give no uncertainty for 5.13%.
- pm005: ours − Wu = +0.171 pp, and Wu's 3.80% lies outside our CI [3.891, 4.193]%. Not explained.
- pmp: ours − Wu = −0.194 pp, and Wu's 2.85% lies outside our CI [2.617, 2.707]%. Not explained.
- re050 and re090: Wu do not tabulate a value, so there is no comparison.

Observations that qualify the table:
- Where sc > 1, only the first crossing is reported. This applies to herald (9,13) in pm005 and to herald (7,9) in re098 and pmp.
- For `wu_blind_mwpm` in re098, re100 and pm005, a (9,13) crossing lies in the grid, but sc_u = 0: no sign change occurs among the p with P_L_shot < 0.4 at both distances.
- The blind (7,9) point estimate in re100 (5.839%) lies outside its own bootstrap CI [4.211, 5.600]%.

![Wu validation figure](../figures/wu_validation.png)

*Figure (`figures/wu_validation.png`).* *(a)* The herald decoder's per-shot logical failure probability after T = d rounds vs p, for d = 9 (dashed) and d = 13 (solid), at R_e = 0 and 0.98, with Wilson 95% intervals. Diamonds mark our (9,13) crossings, with horizontal 95% bootstrap CI bars. Dashed vertical lines mark Wu's 0.937% and 4.15%. *(b)* The (9,13) herald crossing vs R_e at p_m = 0, with 95% CI (filled circles), and Wu's tabulated values at R_e = 0, 0.98 and 1.0 (open squares). Wu's figure gives no tabulated values at R_e = 0.5 and 0.9, and none are read off it.

## Differences from Wu et al.

The effect on the threshold of each item below was not measured here.
1. Finite-size window: we use the (9,13) crossing, where Wu used (11,15). The pair was set by a pre-registered memory rule ([WU_DPAIR_DECISION.md](remediation/WU_DPAIR_DECISION.md)). Effect not measured.
2. Decoder: we use MWPM with heralded edges at weight 0, where Wu used weighted Union-Find, which they describe as optimal for pure erasure. Effect not measured.
3. Code and circuit: we use the rotated CSS surface code with noiseless H gates, where Wu used the planar XZZX code with native CZ/CNOT. Wu state that the two perform identically under unbiased noise. Effect not measured.
4. Decoder noise model: our sampler is exclusive per gate, as in Wu. Our decoder's model treats Pauli and erasure as independent channels. Effect not measured.
5. Edge weights: Wu truncate ln(p') to integers. We pass unrounded weights to PyMatching, which maps them to its own internal integer scale. Effect not measured.

## Reproduce

Machine: Apple M2 (4 performance + 4 efficiency cores), 8 GiB RAM, macOS. Software (from the JSON): Python 3.14.6, stim 1.16.0, PyMatching 2.4.0, sinter 1.16.0, NumPy 2.5.1, SciPy 1.18.0, matplotlib 3.11.0.

| step | command | wall time |
|---|---|---|
| collect one config | `uv run python experiments/collect_wu.py experiments/configs/wu/<name>.yaml --workers 8` | re000 35 s; re050 4 min 49 s; re090 20 min 31 s; re100 18 min 6 s; pmp 30 min 27 s; re098 36 min 19 s; pm005 42 min 23 s; re050 extension 7 min 25 s |
| G5 check | `uv run python scripts/wu_g5_check.py` | 0.4 s |
| table, JSON, figure | `uv run python scripts/wu_validation.py` | 29 s |

The six batch wall times (re050 through pm005) are the start and end timestamps in the message of commit 1e71598. The re000 time and the re050 extension time come from the creation and last-write times of their log files, which are not committed.
