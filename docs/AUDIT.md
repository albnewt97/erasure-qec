# Audit reproduction (Phase 0)

An external audit reported that the circuit builder, DEM partition, and the
herald-vs-blind ablation are correct, but that the noise-model semantics, the
threshold-fitting pipeline, and several README claims are wrong or overstated.
Before changing anything I reproduced every reported number independently.
`scripts/audit_checks.py` regenerates this table; run it with
`uv run python scripts/audit_checks.py`.

## Results

| Check | Expected (audit) | Measured (here) | Verdict |
|---|---|---|---|
| baseline `real_baseline_pauli.csv`, herald, `d_min=7` | p_th = 1.376%, ν = 1.45 | p_th = 1.376%, ν = 1.45 | confirmed (exact) |
| baseline herald, `d_min=None` | p_th = 1.708%, ν = 2.22 | p_th = 1.708%, ν = 2.22 | confirmed (exact) |
| `real_erasure_r50.csv`, herald, `d_min=7` | 2.185% ± 0.278% (README says 2.30% ± 0.15%) | 2.185% ± 0.278% | confirmed — README overstated |
| r50 herald, `d_min=None` | 3.162% = p_arr.max, χ²/dof ≈ 2.89, `converged=True` | 3.162%, pinned at p_arr.max, χ²/dof = 2.89, `converged=True` | confirmed — bound-pin reported as converged |
| r50 blind, `d_min=None` | ν = 3.75 (unphysical) | ν = 3.75 | confirmed |
| χ²/dof of `d_min=7` fits | 0.11–0.32 | 0.26, 0.30, 0.11, 0.32 | confirmed — over-parameterised (5 params, 4–19 dof) |
| max P_L_shot inside any fit window | 0.42–0.50 | 0.485 | confirmed — coin-flip data is being fitted |
| DEM heralded error mass, d=5, p=0.02, r_e=0.98 | 0.475 | 0.475 | confirmed — not 0.98 |
| syndrome density p=0.02, r_e=0/0.5/0.98 | 0.2149 / 0.1993 / 0.1837 (−15%) | 0.2144 / 0.1985 / 0.1833 | confirmed (seed differs; −15% holds) |
| non-identity Pauli prob per 2q gate | p → 0.875p → 0.755p | 1.000p → 0.871p → 0.752p | confirmed — budget shrinks with r_e |
| herald-free shot fraction, d=5, p=0.02, r_e=0.98 | ~13 / 40 000 | 20 / 40 000 (0.05%) | confirmed — fast path is ~0.03–0.05% of shots |
| ablation ratio p=1%, r_e=0.98 | d3: 1.9×, d5: 6.4×, d7: ~30× (n≈22) | d3: 1.9× (h=993,b=1836), d5: 6.7× (h=200,b=1308), d7: 29.3× (h=28,b=807) | confirmed — reproduces |

Small numeric differences (syndrome density, herald-free fraction, ablation
d5/d7) are sampling/seed noise; every finding holds. The non-identity Pauli
column here uses the exact per-gate expression
`1 − (1 − p(1−r_e))·(1 − ¾·p·r_e/2)²`, so it reads 0.871/0.752 rather than the
audit's linearised 0.875/0.755 — same effect.

## Root causes (to fix in Phase 1)

1. **Noise budget is not held constant across r_e.** The per-2q-gate
   non-identity error probability falls from `p` to `0.752p` as `r_e → 0.98`
   (`noise/model.py`: `DEPOLARIZE2(p(1−r_e))` + per-qubit `HERALDED_ERASE(p·r_e/2)`
   with no compensation for the ¾ non-identity fraction of an erasure). A
   threshold that "rises" partly because there is simply less noise at high
   r_e is not a clean measurement of erasure conversion.

2. **"r_e = 0.98" does not mean 98% of errors are heralded.** Only 47.5% of the
   DEM error mass is heralded at d=5, p=0.02, r_e=0.98, because measurement /
   reset / idle noise stays unheralded at rate `p` while the 2q-gate budget
   shrinks. The per-2q-gate heralded fraction is ~0.97, but the circuit-wide
   fraction is ~0.475, so "98% erasure conversion" is not what the circuit
   implements.

3. **The fit reports `converged=True` when the optimiser pins `p_th` at a bound.**
   `fit_threshold` treats any non-raising `curve_fit` call as converged, so a
   result sitting exactly on `p_arr.max()` with χ²/dof ≈ 2.89 is reported as a
   clean fit.

4. **The 5-parameter ansatz is over-parameterised for this data.** The `d_min=7`
   windows have 4–19 degrees of freedom against 5 free parameters, giving
   χ²/dof of 0.11–0.32 (fitting noise, not signal). ν is essentially
   unconstrained (e.g. r50 blind `d_min=None` gives ν = 3.75).

5. **The fit window admits coin-flip data.** The saturation cap is on the
   *per-round* p_L (0.4), which for T rounds corresponds to a *shot-level*
   P_L_shot near 0.5; points with P_L_shot up to 0.485 enter the fit.

## Reproducibility problem (beyond the audit list)

`data/` is git-ignored. The committed CSVs are only
`tests/fixtures/{baseline_pauli, erasure_r50, erasure_r98, real_baseline_pauli,
real_erasure_r50}.csv`. **There is no committed real r_e=0.98 sweep** — the
README's headline (threshold rising to 6.27% at r_e=0.98, the 4.5× / 2.7×
decomposition) is computed from the git-ignored `data/erasure_r98.csv`, which a
cloner cannot reproduce. Every real erasure sweep on disk was also generated
with the un-fixed noise model above, so its threshold is measured under the
shrinking-budget artifact. These claims cannot stand as written.

## Resolution (Phase 1)

The findings above describe the code *before* the fixes. They are left intact
as the historical record; `scripts/audit_checks.py`'s "expected" column is that
pre-fix baseline, so several of its rows now differ on `HEAD` **by design** (the
script confirms the fixes took effect). Post-fix state:

**Root cause 1 & 2 — noise budget (fixed, commit `6b6b87d`).** `noise/model.py`
now sizes the heralded-erase rate to `(2/3)·p·r_e`, so `2·q·¾ = p·r_e` and the
per-2q-gate non-identity budget is `p` for every `r_e`. Re-measured: non-identity
Pauli per gate is now `1.000p / 0.994p / 0.995p` at `r_e = 0 / 0.5 / 0.98` (was
`1.000p / 0.871p / 0.752p`), and syndrome density is flat at `0.2144 / 0.2088 /
0.2028` (was falling −15%). The circuit-wide heralded DEM mass at d=5, p=0.02,
r_e=0.98 rose to 0.547 (was 0.475), since the erasure rate is no longer
shrunk; it is still well below `r_e` because meas/reset/idle stay unheralded —
so "r_e" remains the *2q-gate* heralded fraction, not a circuit-wide one, and
the README must say so.

**Root causes 3, 4, 5 — fit pipeline (fixed, commits `0ce6ff5`, `df70ddc`).**
- Bound-pinning is now reported as `converged=False`. Both r50 and r98 herald
  `d_min=None` fits, which used to report `converged=True` pinned at a bound,
  now return `converged=False` with a message naming the pinned parameter.
- The fit requires ≥ 3 degrees of freedom (≥ 8 points for 5 params), so it
  cannot "converge" on a dof≈1 window.
- σ is now the ~1σ Wilson standard error (95% half-width ÷ Z_95). On the valid
  committed baseline this brings χ²/dof from ~0.2 to ~1.0 (the model now fits at
  the noise level rather than "too well"). `FitResult` carries a `chi2_dof`
  field.
- Saturated points are excluded on the *per-round* rate (the variable the
  collapse is fit in), not the shot-level rate. An initial shot-level cut
  (`0ce6ff5`) was reverted in `df70ddc`: at large T a near-crossing point still
  has `P_L_shot → ½`, so a shot cut deletes near-crossing large-d data and
  breaks high-threshold fits; imprecise near-saturation points are instead
  downweighted by their (now correct) 1σ Wilson errors. The result is
  cap-insensitive for caps in {0.2, 0.25, 0.3, 0.4}.
- The p_th initial guess is clamped into the windowed p-range so `curve_fit`
  cannot raise on an out-of-bounds guess (it now proceeds and reports
  bound-pinning instead).

**Bootstrap CI (fixed, commit `8777283`).** The parametric bootstrap measured a
*narrower* estimator than the point fit: it refit an unweighted model from the
point-estimate `popt` on a frozen window and dropped failed replicates. It now
re-runs the whole pipeline per replicate (resample all points → estimate_crossing
→ window → weighted fit), counts failures (`n_boot_failed`), and reports a
**percentile CI** (the p_th distribution is skewed). This exposed real
overconfidence — e.g. the r_e=0.98 blind fit's reported `± 0.083%` became a 95%
CI spanning ~`[1.5, 8]%` that shifts with the seed. Calibration is checked by
`test_bootstrap_ci_covers_known_p_th_at_nominal_rate` (90% coverage of a known
p_th over 40 seeded realizations).

**Re-measured thresholds (committed data, honest `d ≥ 7` pipeline; 95% bootstrap
percentile CI, `seed = 0`):**
| `r_e` | herald `p_th` [95% CI] | blind `p_th` [95% CI] | notes |
|---|---|---|---|
| 0.0 | 1.376% `[1.33, 1.47]` (χ²/dof 1.02) | 1.440% `[1.35, 1.57]` (1.15) | agree within CIs, as required with no heralds |
| 0.5 | 2.321% `[2.19, 2.64]` (0.62) | 1.491% `[1.42, 2.15]` (0.47) | herald point above blind; blind CI skewed up |
| 0.98 | **not resolved** (ν rails to bound) | **not resolved** (`[1.5, ~8]`, bistable) | at near-full conversion the fit resolves neither |

Both r_e=0.98 thresholds are **non-results** and reported as such — not forced.
Where the fit resolves, with the budget held constant the blind threshold barely
moves (1.44% → 1.49% from r_e=0 to 0.5), so the erasure gain is almost entirely
herald-conditioning; the old README's large "erasure conversion alone" benefit
was partly the shrinking-budget artifact.

**Reproducibility.** The stale (pre-fix) sweeps were moved to
`data/stale_old_model/` (never deleted). The r_e=0.5 and r_e=0.98 sweeps were
re-collected under the fixed model and committed (`data/*.csv`, un-ignored). The
strongest result — model-normalisation-independent and needing no threshold fit
— is the herald-vs-blind ablation, computed directly from circuits with a fixed
seed (`scripts/ablation_table.py`).

## Paired-decoder separation (Phase 4)

The README concluded the `r_e = 0.5` herald-vs-blind threshold separation was
"marginal at 95%" because the two *marginal* bootstrap CIs nearly touch. That is
not a valid significance test — overlapping (or nearly-touching) CIs do not
imply a non-significant difference, because `var(Δ)` is not the sum of the
marginal variances. The correct statistic is a bootstrap of
`Δ = p_th(blind) − p_th(herald)` itself (`bootstrap_threshold_difference`;
reproduce with `scripts/paired_separation.py`).

**1. Are the sweeps paired? No.** `sinter.collect` samples each decoder
independently: at a given `(p, d)` the herald and blind rows have *different*
shot counts (each decoder ran until it hit `max_errors = 1000`, and herald makes
fewer errors so runs more shots). Differing-shot-count `(p, d)` pairs: **12/90**
(r_e=0), **57/90** (r_e=0.5), **97/145** (r_e=0.98) — if the shots were shared
the counts would be identical. Moreover the CSV schema records only marginal
`(shots, errors)` per decoder, not the per-shot joint (herald-wrong × blind-wrong)
a paired bootstrap needs; so even shared shots could not be paired from this
data. **Real data therefore uses the unpaired difference bootstrap**, whose CI is
conservative (it exploits no correlation). `bootstrap_threshold_difference`
supports a paired mode (2×2 shared-shot joint counts) used by the synthetic
tests; the reported herald/blind bootstrap correlation on real data is ~0.02–0.04,
confirming pairing would have bought essentially nothing here.

**2. Marginal p_th 95% CIs, `n_boot` 200 (before) → 1000 (after), d ≥ 7:**

| `r_e` | decoder | p_th | CI @200 | CI @1000 |
|---|---|---|---|---|
| 0.0 | herald | 1.376% | [1.33, 1.47] | [1.33, 1.47] |
| 0.0 | blind | 1.440% | [1.35, 1.57] | [1.36, 1.58] |
| 0.5 | herald | 2.321% | [2.19, 2.64] | [2.18, 2.76] |
| 0.5 | blind | 1.491% | [1.42, 2.15] | [1.42, 2.15] |
| 0.98 | herald | not resolved | — | — |
| 0.98 | blind | 1.636% | [1.52, 6.99] | [1.52, 6.35] |

Raising `n_boot` mostly stabilised the noisy tails (the r_e=0.5 herald upper
endpoint 2.64% → 2.76%); it did not change any convergence verdict.

**3 + 4. Δ = p_th(blind) − p_th(herald), 95% CI, `n_boot = 1000`, seed 0:**

| `r_e` | Δ | 95% CI | excludes 0? | corr | n_paired_failed |
|---|---|---|---|---|---|
| 0.0 (control) | +0.064% | [−0.041, +0.203] | **no** | +0.022 | 51/1000 |
| 0.5 | −0.830% | [−1.240, −0.661] | **yes** | +0.038 | 102/1000 |
| 0.98 | non-result (herald fit does not converge) | — | — | — | — |

The r_e=0 **control passes**: Δ is consistent with zero, as it must be (herald
and blind are the same decoder with no heralds). The `n_paired_failed` counts
mean each Δ CI is conditioned on both decoders converging on the replicate.

**5. Seed stability** (r_e=0.5 Δ CI, 5 seeds): lower-endpoint spread 0.036%,
upper-endpoint spread 0.018%; **all five seeds exclude zero**. No unresolved
instability.

**Summary.** Under the correct statistic the `r_e = 0.5` separation **is
significant**: Δ = −0.83% with a 95% CI of [−1.24, −0.66] that excludes zero and
is stable across seeds. This is the first Phase to *restore* a claim rather than
weaken one — but the restoration comes entirely from using the difference
statistic, not from pairing (the sweeps are unpaired; correlation ≈ 0.04) and
not from raising `n_boot` (the Δ CI excluded zero at 200 too). The
deterministic ablation remains the primary evidence, as it needs no fit at all.

### Are the discarded replicates benign? (Phase-4 follow-up)

Each Δ CI above is conditioned on **both** decoders converging on the replicate,
and the discards are not negligible (102/1000 at r_e=0.5, 51/1000 at r_e=0), so
we checked whether they are missing-at-random with respect to Δ. Reproduce with
`scripts/paired_separation.py` (section 6).

**1. Failure breakdown by guard** (which decoder, which guard):

| `r_e` | discards | breakdown |
|---|---|---|
| 0.0 | 51/1000 | herald bound-pin 26, herald insufficient-window 8, blind bound-pin 11, blind insufficient-window 6 |
| 0.5 | 102/1000 | blind bound-pin 73, herald bound-pin 29, herald insufficient-window 1 |

At r_e=0.5 **~99% of discards are bound-pinning** — the resampled crossing sits
at/above the fit window, i.e. the fit *did* run and reported the crossing is out
of range. That is informative, not a sparse-window artefact.

**2. Missing-at-random diagnostic.** Compare the blind `p_th` of discarded vs
kept replicates (the concern's direction is blind pinning high):

| `r_e` | kept blind `p_th` (median / q90) | failed blind `p_th` (median / q10) | KS stat, p | verdict |
|---|---|---|---|---|
| 0.0 | 1.445% / 1.525% | 1.444% / 1.310% | 0.168, p = 0.16 | **missing-at-random** (no clustering) |
| 0.5 | 1.496% / 1.561% | 1.800% / 1.468% | 0.586, p = 3.6×10⁻³⁰ | **NOT missing-at-random** — failed-blind `p_th` is higher |

So at r_e=0.5 the discards **are** biased: they are disproportionately the
replicates where blind lands high (bound-pinning near the top of the grid),
which are the replicates with Δ nearest zero. Excluding them nudges the CI away
from zero. This is a genuine caveat and is reported wherever the claim appears.

**3. Failure rate before/after a mechanical fix.** We looked for a mechanical
cause per Task 1d. There is none to fix: the dominant guard is bound-pinning of
a genuinely *bistable* blind crossing (its marginal CI already spans [1.42, 2.15]
and its p_th shifts by grid steps under resampling), not a degenerate window
(one insufficient-window discard) or bounds mis-placed off the data. Widening the
p_th bounds beyond the data range would only let p_th extrapolate past the
observed grid — worse, not better. We therefore changed nothing (changing the
estimator to reduce discards would be chasing a better Δ, which we do not do).
Failure rate is **10.2% before and after** at r_e=0.5.

**4. Sensitivity to the discards.** The first version of this section leant on an
*imputation* (`directional_sensitivity_ci`); a later self-audit found it too weak
to trust, and it is replaced below by a tipping-point bound. Both are recorded
for reproducibility (`scripts/paired_separation.py`).

*Audit of the imputation (Task 1).* Imputing the 102 discards by resampling the
observed Δ draws in the top (least-negative) decile gave [−1.23, −0.64], which
**barely moved** from the observed [−1.24, −0.66]. That looked suspicious. The
imputed values are **not** a single point (63 distinct), so the resampling works,
but their quantiles are min/25/50/75/max = −0.73/−0.71/−0.69/−0.66/+0.38%: **75%
sit below −0.66%**, only 27 of 102 above it. The empirical top-decile is
density-weighted toward its *lower* edge (−0.66% is the 97.5th percentile, near
the TOP of the decile), so sampling it is only mildly pessimistic — not the
"plausible-pessimistic" case the earlier docstring claimed. The reported CI was
arithmetically correct; its *characterisation* was overstated. Because the
imputation's answer depends entirely on the chosen decile, it is replaced.

*Tipping-point bound (Task 2a).* Only **2** of the 898 successful draws have
Δ ≥ 0. The 97.5th percentile of 1000 draws (numpy 'linear' convention) reaches
zero once **26** are ≥ 0, so **the claim survives unless ≥ 24 of the 102 discards
would have given Δ ≥ 0** (verified against a direct construction in the tests).

*Partial-information implied Δ (Task 2b).* We do not have to guess: **101 of 102**
discards recorded *both* decoders' `p_th` (one converged, one bound-pinned), and
1 recorded blind only — **0 are unbounded** (Task 2c). Forming Δ = blind − herald
from the recorded values (filling the 1 missing herald from the successful
median) gives an implied-Δ distribution of min/25/50/75/max =
−1.73/−0.80/−0.64/−0.49/+2.42%, with **only 2 of 102 at Δ ≥ 0**. For a discard to
reach Δ ≥ 0 its blind `p_th` must reach herald's centre **2.37%**; the recorded
failed-blind median is **1.80%** and only **2** values reach 2.37%. Because a
blind bound-pin `p_th` is a *lower* bound on the true crossing, 2 is itself a
lower bound on the count — but reaching 24 would require blind's true crossing to
exceed herald's in 24 resamples where the recorded value sits ~0.57% below it,
which blind's marginal distribution (median 1.50%, 97.5th percentile 2.15%) does
not support.

**Verdict.** The conditioning is **not** benign at r_e=0.5 — the discards are not
missing-at-random and cluster at high blind `p_th` (that finding, item 2, stands
regardless). **The `r_e = 0.5` claim survives, unchanged in value, with the
caveat attached:** significant *conditional on convergence* (Δ = −0.83%, CI
[−1.24, −0.66]); 10.2% of replicates discarded and **not** missing-at-random; but
the tipping point is **24** discards at Δ ≥ 0 while the recorded partial
information implies **2**, and none are unbounded, so the margin is large. The
r_e=0 control has tipping point 0 (its Δ already includes zero, as a control
must) and missing-at-random discards. The unpaired difference bootstrap is
**conservative** (pairing could only have narrowed the interval), so clearing
zero under it is if anything a stronger result — and the deterministic ablation,
which needs no fit and no convergence filtering, remains the primary evidence.

## Threshold-fitting methodology (reference)

The full narrative behind the one-line summaries in the README.

**Per-round conversion.** Sinter reports a shot-level `P_L_shot` over `T = d`
rounds; the comparable per-round rate is `p_L = ½(1 − (1 − 2 P_L_shot)^{1/T})`,
with fixed points `P_L_shot = ½ → p_L = ½` and `T = 1 → p_L = P_L_shot`. The
collapse is fit in this per-round variable, which matters at high thresholds: at
`d = 11` (`T = 11`) a legitimate near-crossing point still has `P_L_shot ≈ 0.45`,
so the saturation cut acts on `p_L`, not `P_L_shot`.

**Finite-size ansatz.** Near the crossing the collapsed data fit the quadratic
finite-size-scaling ansatz of Wang, Harrington & Preskill: `p_L(p, d) = A + B x +
C x²`, `x = (p − p_th) d^{1/ν}`, for `(p_th, ν, A, B, C)` by weighted least
squares, weighting each point by its 1σ Wilson standard error (95% Wilson
half-width ÷ z₀.₉₇₅) in per-round space. Points at or above per-round `p_L = 0.4`
(saturating toward ½, outside the local ansatz) are excluded, and the fit reports
its χ²/dof.

**`estimate_crossing`.** The fit window centres on a data-driven crossing
estimate. That estimate is not reliable on saturated data: above threshold the
curves re-converge toward ½, producing a spurious second minimum of the
cross-distance spread, and a ragged high-`p` tail can pull it past the real
crossing. An ordering-inversion guard rejects candidates where `p_L` already
increases with `d`, but it is only a starting point for the window; every quoted
threshold comes from the `d ≥ 7` collapse fit, never the raw crossing estimate.

**Bootstrap.** Uncertainty is a 95% bootstrap percentile CI from a parametric
bootstrap that re-runs the *entire* pipeline per replicate: resample
`errors ~ Binomial(shots, P_L_shot)` over all points, then re-run the crossing
estimate, the window selection, and the same weighted fit. Failed replicates are
counted (`n_boot_failed`), not dropped; the interval is a percentile CI (the
`p_th` distribution is skewed), not `± σ`. `n_boot` defaults to 1000 (at 200 the
percentile tails were noisy). History: an earlier bootstrap refit an *unweighted*
model from the point estimate on a *frozen* window and dropped failures, all of
which reported the CI too narrow.

**Convergence vs resolution guards.** `fit_threshold` reports two booleans.
`converged=False` when the optimiser did not return a usable fit — insufficient
window points, a `curve_fit` raise, or a parameter *pinned* at a bound (the data
did not constrain it; e.g. `r_e = 0.98` herald rails ν to its upper bound).
`resolved=False` is the separate question "does the fit resolve a threshold": a
converged fit still fails it if the relative bootstrap CI width
`(ci_hi − ci_lo) / p_th ≥ 1.0` (the CI is at least as wide as the threshold
itself). The constant is chosen from the committed fits, not by taste: every
reported result has relative width ≤ 0.49 (`r_e = 0.5` blind), while the
`r_e = 0.98` blind non-result is ≈ 3.0 — a ~6× gap, so 1.0 rejects only that fit.
Every `FitResult` carries a structured `reason` code (`ok`, `insufficient_data`,
`curve_fit_error`, `bound_pin`, `unresolved_ci`). Only a `resolved` fit is quoted
as a `p_th`, in code (plotting draws no line for an unresolved fit) and in docs.
