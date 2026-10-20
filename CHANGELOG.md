# Changelog

## v0.2.0 (unreleased)

Claim and caveat IDs (C, Y, K, A) refer to [CLAIM_LEDGER.md](docs/remediation/CLAIM_LEDGER.md) and [CAVEAT_LEDGER.md](docs/remediation/CAVEAT_LEDGER.md).
Numbered items refer to [P7_CARRY_FORWARD.md](docs/remediation/P7_CARRY_FORWARD.md).
v0.1 means the README at commit `da47e3a`.

### Changed claims

Every CI in this changelog is a bootstrap 95% percentile interval that covers binomial noise only; a pairwise-crossing CI excludes interpolation error inside the bracket.
CI endpoints can move in the second decimal with the NumPy version ([METHODS.md](docs/METHODS.md#numpy-dependent-ci-endpoints)).

- The headline threshold is now the per-shot (7,11) pairwise crossing, `P_L` after `T = d` rounds (`scripts/headline_results.py`).
  v0.1 quoted the per-round `d ≥ 7` collapse fit, which is now a secondary column (item 1).
- The thresholds and the herald − blind difference are computed on the densified sweeps (`extra_p_shot`, commits `6337619` and `090c341`).
- The v0.1 per-round values (C006, C018, C019, C020) are replaced by the rendered `thresholds` table.
  On the dense data the per-round `r_e = 0.5` herald fit is 2.33% [2.21, 2.73] (bootstrap 95% CI).
- At `r_e = 0.98` the per-shot collapse fit resolves both thresholds: herald 3.83% [3.74, 3.94], blind 1.18% [1.16, 1.21] (item 15; C010).
  The per-round fit still resolves neither (herald `bound_pin`, blind `unresolved_ci`).
- The per-shot collapse fit now sits next to each pairwise crossing, with the gap as a systematic column (item 5).
- `Λ` is labelled per-round at every appearance, and thresholds per-shot (item 16).
- The ablation now covers `p = 0.5%` as well as `p = 1.0%` (`results/ablation.json`, commit `05f59d1`).
  Its blind error counts now run from 52 to 6290, not 3016 to 6290 (C075).
- CI endpoints print to two decimals, with a footnote that they can move in the second decimal with the NumPy version (C017, K025).
- The comparison with Wu et al. now uses the `wu2022` convention and the (9,13) per-shot herald crossing (`docs/VALIDATION.md`).
  v0.1 compared a repo-model `d ≥ 7` collapse fit with Wu's reported thresholds, of which the Fig. 3a values come from Wu's (11,15) crossing (C036, item 13).

### Withdrawn claims

- "blind barely moves 1.44% -> 1.49%" (C002, C005, C044, C045, K007, K076, K077).
  These were per-round collapse fits on the pre-dense sweeps.
  The per-shot (7,11) blind crossing is 0.98% [0.95, 1.02] at `r_e = 0` and 1.10% [1.03, 1.14] at `r_e = 0.5`.
- "Delta = -0.83%" (C007, C027–C032, K029, K041, K044–K047, A007, A029–A036).
  This was a per-round collapse-fit difference with its tipping-point analysis, from `scripts/paired_separation.py` on the pre-dense sweeps.
  It is replaced by herald − blind of the per-shot (7,11) crossings: 0.49 pp [0.40, 0.60] at `r_e = 0.5`.
- "the r_e = 0.5 dagger footnote" (C003, C004, C021–C026, K008, K032–K039, A005, A006, A008; item 3).
  The footnote explained the per-round `r_e = 0.5` blind CI upper endpoint, 2.154%, with a 30-seed probe that was never committed.
  On the dense data that endpoint is 1.80%, which coincides with the grid point `p = 0.018` (`results/headline.json`).
- "the r_e = 0.98 threshold is not resolved" (K015, K031, K050, K136, K138, K139; item 15).
  It holds for the per-round fit only; the per-shot fit resolves both decoders at `r_e = 0.98`.
- "an unexplained sign problem" (C039–C042, K066–K070, K075).
  It compared the repo model, which has extra measurement, reset and idle noise, with Wu's model.
  That comparison is no longer made, and the proposed mechanism (the observable) was never tested.
  The comparison with Wu now uses the `wu2022` convention.
- "R_e defined exactly as our r_e" (K051, K052, K087).
  Wu's erasure removes both qubits in one event, exclusive with the Pauli, at probability `p·R_e`.
  The repo convention applies independent per-qubit erasures at `(2/3)·p·r_e` alongside `DEPOLARIZE2(p·(1 − r_e))`.
- "they never restrict to one observable or one basis" (K060, K061).
  v0.1 gave no quotation from Wu et al. supporting it.
- "follows Wu et al." (K080).
  The repo convention differs from Wu's model (see the previous entry); the `wu2022` convention follows it.
- "each assumption is systematically optimistic" (K120).
  Sahay et al. report a higher threshold under Z-biased erasure than under unbiased erasure for the XZZX code; this was not measured for the model here.
  The direction of the effect of each assumption on the results here was not measured.
- "models dual-rail superconducting cavities" (K003).
  The committed sweeps convert only two-qubit-gate errors, and no dual-rail-specific noise is modelled.
- "group: per-round paired-separation diagnostics" (C011–C014, C077, C078).
  These were seed spans and grid pins of the per-round `r_e = 0.98` blind CI from `scripts/paired_separation.py --full-stability` on the pre-dense sweeps.
  On the dense data that fit is published only as "not resolved (unresolved_ci)".
- "group: old literature-comparison numbers" (C037, C038).
  The noise-location counts of a `d = 5`, `T = 5` circuit had no committed producer.
- "group: history notes" (C046, C065).
  The 1.6× / 2.7× decomposition and the pre-2026-09-05 ablation bound are in [docs/archive/AUDIT_HISTORY.md](docs/archive/AUDIT_HISTORY.md).
- "group: ablation detail not in results/ablation.json" (C067, C068).
  The renderer prints only the lower bound for a gated cell, so the `d = 11` herald CI and its point ratio are not quoted.
- "group: throughput" (C081–C088).
  The decoder throughput numbers came from single runs on one machine, with no committed record of the medians.
- "group: timings and verification records at da47e3a" (C015, C056 timing, C089–C096, C094).
  These runtimes and the 2026-09-06 byte-identity check have no recorded run at the current commit; the README writes "not timed".
- "group: stale YAML comments" (Y002–Y005, Y007–Y016).
  The threshold values in `experiments/configs/*.yaml` comments are removed and point to `results/headline.json` (commit `2c2677e`).

### Fixed

- Per-shot crossings computed on the pre-dense grid had wide brackets, for example [1.00, 1.47]% for the `r_e = 0.5` blind (7,11) crossing (item 2).
  Their CIs, like the current ones, excluded interpolation error inside the bracket.
  Densification moved the `r_e = 0.5` blind (7,11) crossing from 1.181% [1.130, 1.233] to 1.104% [1.028, 1.136] (`results/threshold_table_committed.json` → `results/headline.json`).
  The new value lies outside the earlier 95% CI.
  Every pre-dense per-shot number is superseded, including "about 3.9%" for the `r_e = 0.98` herald crossing.
- At `r_e = 1` with no unheralded noise, `partition_flattened_dem` declared 0 observables and `HeraldMatchingDecoder` crashed (item 9).
  Commit `9484099` declares all observables; no committed sweep reached that case.
- Docstrings and `experiments/configs/erasure_r98.yaml` called the model "biased-erasure"; it is unbiased (commit `bf13a56`).
- `test_ablation_herald_beats_blind_on_biased_erasure_circuit` is renamed `test_ablation_herald_beats_blind_on_erasure_circuit` (commit `bf13a56`).
- ruff now targets Python 3.12, matching `requires-python` (commit `be1e9e0`).

### Archived

- `PLAN.md`, `docs/AUDIT.md` and `docs/DEVLOG.md` moved to [docs/archive/](docs/archive/README.md) as `PLAN_v1.md`, `AUDIT_HISTORY.md` and `DEVLOG.md` (commit `a82af5d`).
  Methods now live in [docs/METHODS.md](docs/METHODS.md) (K006).
- `scripts/paired_separation.py` reproduces its recorded output only at `da47e3a`, and `scripts/audit_checks.py` reproduces the archived audit table (item 14).
