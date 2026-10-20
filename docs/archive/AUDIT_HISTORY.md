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
CI whose upper endpoint runs 5.5–8.5% depending on the seed (`scripts/paired_separation.py
--full-stability`). Calibration is checked by
`test_bootstrap_ci_covers_known_p_th_at_nominal_rate` (90% coverage of a known
p_th over 40 seeded realizations).

**Re-measured thresholds (committed data, honest `d ≥ 7` pipeline; 95% bootstrap
percentile CI, `seed = 0`, **`n_boot = 200`** — the default at the time of Phase 1;
endpoints to one decimal — see "Reproducibility and environment sensitivity" for why
two is not justified):**

> **Superseded endpoint (noted 2026-09-06, table left as measured).** `n_boot` was
> later raised to 1000. Only one cell moves at one-decimal precision: the `r_e = 0.5`
> **herald** CI reads `[2.2, 2.6]` here at `n_boot = 200` and `[2.2, 2.8]` at
> `n_boot = 1000`. `README.md` and `docs/DEVLOG.md` quote the current `[2.2, 2.8]`;
> the two are the same measurement at different `n_boot`, not a disagreement. The
> before/after pair is tabulated in "Paired-decoder separation (Phase 4)" below.

| `r_e` | herald `p_th` [95% CI] | blind `p_th` [95% CI] | notes |
|---|---|---|---|
| 0.0 | 1.376% `[1.3, 1.5]` (χ²/dof 1.02) | 1.440% `[1.4, 1.6]` (1.15) | agree within CIs, as required with no heralds |
| 0.5 | 2.321% `[2.2, 2.6]` (0.62) | 1.491% `[1.4, 2.2]` (0.47) | herald point above blind; blind CI skewed up |
| 0.98 | **not resolved** (ν rails to bound) | **not resolved** (lower ≈ 1.5, upper 5.5–8.5 by seed — `paired_separation.py --full-stability`) | at near-full conversion the fit resolves neither |

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
| 0.0 | herald | 1.376% | [1.3, 1.5] | [1.3, 1.5] |
| 0.0 | blind | 1.440% | [1.4, 1.6] | [1.4, 1.6] |
| 0.5 | herald | 2.321% | [2.2, 2.6] | [2.2, 2.8] |
| 0.5 | blind | 1.491% | [1.4, 2.2] | [1.4, 2.2] |
| 0.98 | herald | not resolved | — | — |
| 0.98 | blind | 1.636% | [1.5, 7.0] | [1.5, 6.3] |

Raising `n_boot` mostly stabilised the noisy tails (the r_e=0.5 herald upper
endpoint 2.6% → 2.8%); it did not change any convergence verdict. Both r_e=0.98
blind endpoints are seed-dependent to the first digit (see "Reproducibility and
environment sensitivity") and are shown only to compare the two `n_boot`.

**3 + 4. Δ = p_th(blind) − p_th(herald), 95% CI, `n_boot = 1000`, seed 0:**

| `r_e` | Δ | 95% CI | excludes 0? | corr | n_paired_failed |
|---|---|---|---|---|---|
| 0.0 (control) | +0.064% | [−0.04, +0.20] | **no** | +0.022 | 51/1000 |
| 0.5 | −0.830% | [−1.2, −0.7] | **yes** | +0.038 | 102/1000 |
| 0.98 | non-result (herald fit does not converge) | — | — | — | — |

The r_e=0 **control passes**: Δ is consistent with zero, as it must be (herald
and blind are the same decoder with no heralds). The `n_paired_failed` counts
mean each Δ CI is conditioned on both decoders converging on the replicate.

**5. Seed stability** (r_e=0.5 Δ CI, 5 seeds — the default
`uv run python scripts/paired_separation.py`): lower-endpoint spread 0.036%,
upper-endpoint spread 0.018%; **all five seeds exclude zero**. No unresolved
instability. The ten-seed extension (`--full-stability`) is tabulated in
"Reproducibility and environment sensitivity" below; all 10 seeds exclude zero.

**Summary.** Under the correct statistic the `r_e = 0.5` separation **is
significant**: Δ = −0.83% with a 95% CI of [−1.2, −0.7] that excludes zero and
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
a genuinely *bistable* blind crossing (its marginal CI already spans [1.4, 2.2]
and its p_th shifts by grid steps under resampling), not a degenerate window
(one insufficient-window discard) or bounds mis-placed off the data. Widening the
p_th bounds beyond the data range would only let p_th extrapolate past the
observed grid — worse, not better. We therefore changed nothing (changing the
estimator to reduce discards would be chasing a better Δ, which we do not do).
Failure rate is **10.2% before and after** at r_e=0.5.
(Follow-up 2026-09-06: the "shifts by grid steps" observation here is the same
mechanism that sets the *marginal* CI upper endpoint, where -- unlike in this
difference bootstrap -- the pinned replicates are retained rather than discarded.
See "The marginal CI upper endpoints are bound pins (2026-09-06)".)

**4. Sensitivity to the discards.** The first version of this section leant on an
*imputation* (`directional_sensitivity_ci`); a later self-audit found it too weak
to trust, and it is replaced below by a tipping-point bound. Both are recorded
for reproducibility (`scripts/paired_separation.py`).

*Audit of the imputation (Task 1).* Imputing the 102 discards by resampling the
observed Δ draws in the top (least-negative) decile gave [−1.229, −0.636], which
**barely moved** from the observed [−1.240, −0.661]. That looked suspicious.
(Three decimals here only because both intervals come from the *same* seed-0
draws, so their difference is exact even though each endpoint drifts between
environments — see "Reproducibility and environment sensitivity".) The
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
[−1.2, −0.7]); 10.2% of replicates discarded and **not** missing-at-random; but
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

---

## Metric correction: no-op herald branches (2026-09-05)

**What changed.** `analysis/dem_stats.py::heralded_fraction` counted DEM error
mechanisms whose *only* target is a herald detector. `HERALDED_ERASE` replaces
the qubit with `I/2` and heralds on all four Pauli branches, identity included,
so the flattened decomposed DEM carries one such mechanism per herald detector —
at `d = 5, p = 0.02, r_e = 0.98` that is 800 mechanisms, each at ¼ of the herald
rate (`0.02·0.98/1.5 / 4 = 0.00326667`), confirming they are exactly the `I`
branch. These fire a herald, produce no syndrome detector and no logical
observable, and so flip nothing. Counting them in *both* the heralded mass and
the total inflated the reported fraction. They are now excluded from both sides.
The old definition is retained as `heralded_fraction(..., count_no_op_heralds=True)`
so the numbers below and elsewhere in this document stay reproducible.

**Old vs new** (`scripts/heralded_fraction.py`, `d = 5`, `p = 0.02`):

| `r_e` | `convert_idle` | no-op mass / total | old (counted) | new (excluded) |
|-------|----------------|--------------------|---------------|----------------|
| 0.5   | False          | 1.655 / 17.517     | 0.3045        | **0.2319**     |
| 0.98  | False          | 3.244 / 19.118     | 0.5468        | **0.4542**     |
| 0.5   | True           | 2.582 / 18.524     | 0.4175        | **0.3231**     |
| 0.98  | True           | 5.060 / 21.059     | 0.7198        | **0.6311**     |

**Direction.** The correction lowers every value. It does not weaken any claim:
it strengthens the existing point that `r_e` overstates circuit-wide conversion.
At `r_e = 0.98` the gate-only circuit heralds 45% of its error budget, not 55%
and not 98%. No threshold, ablation or separation result depends on this metric —
it is descriptive only, so no fitted number moves.

**Supersedes.** The `0.547` quoted in "Resolution (Phase 1)" above (and the
`0.475` in "Root causes" item 2, which additionally predates the noise-budget
fix) are pre-correction values under the old definition; they are left in place
as the historical record. The current value at that config is `0.4542`.

---

## Bound correction: the suppression table's `> N×` cells (2026-09-05)

**Defect.** `scripts/ablation_table.py::ablation_cell` computed the sub-threshold
suppression bound as

    ratio_lb = b_pl / h_pl_hi

— herald's Wilson 95% *upper* limit in the denominator, but blind's *point
estimate* in the numerator. The script docstring was candid about this ("blind
held at its well-measured point estimate"), but README.md presented the result as
"a Wilson lower bound `> N×`". **It was not a lower confidence bound on the
ratio.** Nothing constrained the numerator from below, so whenever blind's rate
was over-estimated by sampling the quoted value could exceed the true ratio. The
error is one-sided in the flattering direction: it inflates the headline.

**Fix.** Both ends are now bounded, each mapped through the monotone per-round
conversion:

    ratio_lb = per_round_p_l(wilson_lo(b_err, shots), d)
               / per_round_p_l(wilson_hi(h_err, shots), d)

**Coverage actually claimed.** `wilson_interval` is called at its default
`z = Z_95`, so it returns a *two-sided 95%* interval and each endpoint used is a
*one-sided 97.5%* limit. Bonferroni over the two gives joint coverage
`>= 1 - 0.025 - 0.025 = 95%`, and on that event `ratio_lb <= true ratio`. So each
cell carries a **>= 95% one-sided lower confidence bound on that cell's ratio** —
per-cell, **not** simultaneous over the ten cells, and not paired (the bound
discards the positive correlation from decoding identical shots, which is valid
under any dependence but loose). Wilson coverage is nominal, not exact; at
single-digit counts the "95%" is an approximation with no guaranteed coverage.

**Measured effect** (`uv run python scripts/ablation_table.py`, 100 000 shots,
seed 0, ~46 min single-core, 2026-09-05). Raw counts, herald/blind:

| `d` | `r_e = 0.5` h/b | `r_e = 0.98` h/b |
|---|---|---|
| 3 | 4828 / 6290 | 2765 / 5716 |
| 5 | 3261 / 6065 | 580 / 4924 |
| 7 | 2035 / 5521 | 105 / 4169 |
| 9 | 1322 / 5277 | 23 / 3543 |
| 11 | 810 / 4971 | 3 / 3016 |

Only the two sub-gate cells are quoted as bounds, and both move down slightly:

| cell | old (invalid) | new (valid) | point ratio |
|---|---|---|---|
| `d = 9`, `r_e = 0.98` | `> 106×` | **`> 103×`** | 159.1× |
| `d = 11`, `r_e = 0.98` | `> 352×` | **`> 339×`** | 1034.0× |

As predicted, the move is ~3%: `b_err` is in the thousands, so its Wilson
relative half-width is small. The point ratios and the eight above-gate cells are
unchanged (1.3/1.9/2.8/4.1/6.4 at `r_e = 0.5`; 2.1/8.8/41.2 at `r_e = 0.98`).

**Does this weaken a documented claim?** The numbers barely move, but the
*status* of the claim changes and the README now says so: what was published as a
bound was not one. Separately, and more importantly, the `d = 11` cell rests on
**3 herald errors in 100 000 shots** (herald shot-rate CI `[1.0e-5, 8.8e-5]`,
nearly a decade wide). `> 339×` is weak evidence, not a measurement, and the
README now carries that qualifier immediately adjacent to the number. The
README's previously stated counts (`d = 9` → 23, `d = 11` → 3) were checked
against this run and are **correct as stated**.

**Regression cover.** `tests/test_ablation_table.py` pins `ratio_lb <=` the
superseded expression `<=` the point ratio across a grid of count shapes
(including the zero-error edge and the real sub-gate counts), plus a strictness
test that fails if the numerator is reverted to the point estimate.

---

## Reproducibility and environment sensitivity (2026-09-06)

**Why this section exists.** The docs pin `n_boot = 1000, seed = 0` and quoted
bootstrap CI endpoints to two decimals, which reads as a promise of
bit-reproducibility. It is not one. `sinter.collect` accepts no seed, so the
committed sweeps in `data/` are not regenerable shot-for-shot; and the bootstrap
percentile endpoints move with the NumPy RNG and `np.percentile`
implementation, so even re-analysing the *same* CSVs on a different NumPy does
not return the same endpoints. This section states which quantities are exact,
which drift, and by how much.

**Environments compared.** Both re-analyse the committed `data/*.csv` with
`d_min = 7`, `n_boot = 1000`, `seed = 0`.

| | env A | env B |
|---|---|---|
| Python | not reported | 3.14.6 |
| NumPy | 2.4.4 | 2.5.1 |
| SciPy | 1.17.1 | 1.18.0 |
| stim | 1.16.0 | 1.16.0 |
| pymatching / sinter | not reported | 2.4.0 / 1.16.0 |

Env B is the run recorded below (`scripts/paired_separation.py`,
`scripts/audit_checks.py`, 2026-09-06); env A is an independent reviewer's
report of the same two scripts. Only the analysis path matters here — no
sampling is re-run — so NumPy and SciPy are the operative versions.

### Exactly reproducible across both environments

Everything downstream of the *point* fit carries no RNG at all — `curve_fit` on
fixed data — and agreed to every digit quoted:

| quantity | env B value | confirmed identical in env A? |
|---|---|---|
| `p_th` (r_e = 0 herald / blind) | 1.3764% / 1.4399% | yes (1.376 / 1.440) |
| `p_th` (r_e = 0.5 herald / blind) | 2.3209% / 1.4913% | yes (2.321 / 1.491) |
| `p_th` (r_e = 0.98 blind) | 1.6363% | yes (1.636) |
| Δ = p_th(blind) − p_th(herald), r_e = 0.5 | −0.8295% | yes (−0.830) |
| Δ, r_e = 0 (control) | +0.0635% | yes (+0.064) |
| ν (r_e = 0 / 0.5 herald) | 1.4492 / 1.9019 | not reported by env A |
| χ²/dof (r_e = 0 h/b, 0.5 h/b, 0.98 b) | 1.02 / 1.15, 0.62 / 0.47, 0.20 | not reported by env A |

The last two rows are RNG-free by construction (`curve_fit` and a residual sum
on fixed data), so they are expected to be exact, but env A did not report them
and that is not claimed as measured here.

**No verdict changed between the environments** — `converged`, `resolved`,
`excludes_zero`, which guard fired on a bound-pin, and the tipping-point
survival at r_e = 0.5 are identical in both. The r_e = 0.98 herald fit fails to
converge (ν railed to its bound) in both; the r_e = 0.98 blind fit converges but
is `resolved = False` in both.

### What drifts, and by how much

Two independent sources, measured separately. "env drift" is env A vs env B at
the *same* seed 0. "seed spread" is the range over seeds 0–9 within env B — the
Monte-Carlo error of the bootstrap itself, which is the larger effect
everywhere except the r_e = 0.5 Δ lower endpoint.

Reproduce the env-B and seed-spread columns with
`uv run python scripts/paired_separation.py --full-stability` (~4.7 min; the
default 5-seed run does not print them). The env-A column is an independent
reviewer's report and is not reproducible from this tree.

| quantity (%) | seed-0 env A | seed-0 env B | env drift | seed spread (10 seeds, env B) |
|---|---|---|---|---|
| r_e = 0 herald CI lo / hi | 1.33 / 1.47 | 1.331 / 1.468 | 0.00 / 0.00 | 0.005 / 0.011 |
| r_e = 0 blind CI lo / hi | 1.36 / 1.58 | 1.361 / 1.581 | 0.00 / 0.00 | 0.010 / 0.035 |
| r_e = 0.5 herald CI lo / hi | 2.19 / 2.77 | 2.183 / 2.758 | 0.01 / 0.01 | 0.012 / **0.090** |
| r_e = 0.5 blind CI lo / hi | 1.42 / 2.15 | 1.418 / 2.154 | 0.00 / 0.00 | 0.014 / **0.354** |
| r_e = 0.98 blind CI lo / hi | 1.52 / 6.61 | 1.517 / 6.347 | 0.00 / **0.26** | 0.015 / **3.000** |
| Δ CI (r_e = 0.5) lo / hi | −1.236 / −0.657 | −1.240 / −0.661 | 0.004 / 0.004 | 0.062 / 0.022 |
| Δ CI (r_e = 0) lo / hi | −0.041 / +0.208 | −0.041 / +0.203 | 0.000 / 0.005 | 0.008 / 0.038 |
| `n_paired_failed` (r_e = 0.5) | 101 | 102 | 1 | 77–115 |
| herald/blind bootstrap corr (r_e = 0.5) | +0.031 | +0.038 | 0.007 | −0.051 … +0.038 |
| tipping-point discards (r_e = 0.5) | 24 | 24 | 0 | 24–26 |

**The r_e = 0.98 blind upper endpoint is the outlier and is not a measurement.**
Over seeds 0–9 it takes 5.500, 5.653, 6.347, 6.734, 7.054, 7.067, 7.600, 7.606,
7.609, 8.500 — a 3-percentage-point span, and **three of the ten land exactly on
sweep grid points**: 5.500 (`p = 0.055`, seed 7), 7.600 (`p = 0.076`, seed 5) and
8.500 (`p = 0.085`, seed 4). The 97.5th percentile is there being set by
bound-pinned replicates rather than by a fitted crossing.

**Correction 2026-09-06 (later same day):** "land exactly on sweep grid points"
above was checked at 1e-6 tolerance (`round(u, 6) in grid`), not at float
precision. At full precision the endpoints agree with the grid values to ~1 ULP
(|diff| 7e-18 to 1.4e-17), **not** bit-identically -- the bounded least-squares
solver returns a pinned parameter to within rounding. The substance is unchanged
and in fact strengthened: see "The marginal CI upper endpoints are bound pins"
below, which identifies the pin as the window's **lower** edge, not the upper.

**Sharpened 2026-09-06 (`--full-stability`):** the two *extremes* of the quoted
5.5–8.5% span are precisely the grid-pinned seeds. So the span is not "the
spread of the fitted upper endpoint" — its width is set at both ends by where
the sweep grid happens to stop, and a denser or sparser grid would move it. Read
5.5–8.5% as "unbounded above the resolved region", not as an interval estimate.
This is the same fact the `resolved = False` guard already reports (`rel_ci_width`
2.43–4.26 across those seeds); no digit of that endpoint should be quoted as
though it were measured.

### Consequence for quoted precision

Two decimals on a marginal `p_th` CI endpoint is false precision: the env drift
alone moves the last digit by ±1 (2.76 vs 2.77; 6.35 vs 6.61), and the seed
spread moves it by 0.5–35 units of that digit (300 for the r_e = 0.98 blind
upper endpoint). **The marginal `p_th` CI endpoints and the r_e = 0.5 Δ CI are
therefore now quoted to one decimal** in this document, in `README.md`,
`docs/DEVLOG.md` and `PLAN.md`. The r_e = 0 control Δ CI keeps two decimals:
its endpoints are 0.04–0.20%, an order of magnitude smaller, so its drift
(0.005–0.008) still sits below the last digit. The Δ *point* estimate (−0.83%),
the `p_th` point estimates and χ²/dof keep their digits — they carry no RNG.
Seed-0 values remain in the tables of this section at full precision, labelled
as such, because they *are* the drift evidence rather than a claim about the
estimator.

**Follow-up (2026-09-06): the figure was missed.** The precision reduction above
was applied to the prose and tables but **not** to `analysis/plotting.py`, which
rendered each panel's CI with `{:.2f}` — so `figures/threshold_panels.png`, the
first figure in `README.md`, showed `95% CI [1.33, 1.47]` and `[2.18, 2.76]`
directly above README text stating that only one decimal is justified. The
format string is now `{:.1f}` and the panels read `[1.3, 1.5]` and `[2.2, 2.8]`,
matching the README table. `p_th` (`{:.2f}`) and ν are unchanged — they carry no
RNG. Only `threshold_panels.png` changed (md5 `b6296e8f…` → `bc9b0444…`); the
other three figures, `lambda_vs_p.png` (md5 `442e7f13…`) included, regenerate
byte-identical. This corrects a *display* inconsistency; no fitted value moves.

### What "reproducible" may and may not be claimed

- **May:** the analysis is deterministic given `(CSV, seed, NumPy build)`, and
  every point estimate, χ²/dof and verdict reproduces exactly across the two
  environments tested.
- **May not:** bit-identical CI endpoints. `sinter.collect` takes no seed, so
  the `data/` sweeps cannot be regenerated shot-for-shot at all; and re-running
  the *analysis* on a different NumPy shifts the percentile endpoints as
  tabulated above. Any statement that "every quoted threshold is reproducible"
  must be read as reproducible-to-the-stated-precision, not bit-for-bit.

**Defect found while measuring this (fixed).** `scripts/audit_checks.py::_fit_diagnostics`
built its χ² weights as the raw 95% Wilson half-width, while the fitter
(`_weighted_fit_arrays`) divides that half-width by `Z_95` for a true 1-σ
weight. The script's χ²/dof was therefore low by `Z_95² ≈ 3.84`, printing
`0.26 / 0.30 / 0.16 / 0.12` where the fitter's own weights give
`1.02 / 1.15 / 0.62 / 0.47`. The docs were right and the script was wrong; the
script now shares the fitter's convention. This does not touch the "Results"
row `χ²/dof of d_min=7 fits | 0.11–0.32`, which is the pre-fix historical
record measured with the pre-fix estimator.

## Conditioning in the Lambda bootstrap (2026-09-06)

The threshold fit counts and reports its discarded bootstrap replicates
(`n_boot_failed`, and the missing-at-random analysis above).
`statistics.py::lambda_factor` did the same class of conditioning **silently**:
its bootstrap loop kept a replicate only `if pb > 0.0`, dropping every draw
whose `d+2` denominator was zero, with no count and no comment. This section
records what that conditioning does and what was changed.

**Direction of the bias.** The discards are not missing-at-random, and unlike
the paired-decoder case the direction is knowable *a priori* rather than by
diagnostic: a dropped replicate is exactly one where the `d+2` code logged zero
errors, i.e. one whose unconditional `Λ` is `+∞`. The exclusion therefore
removes precisely the largest values, so the reported upper endpoint is biased
**low**, never high. Beyond a threshold this stops being a bias and becomes a
category error: once the dropped fraction exceeds the interval's own tail mass
`(1 − ci)/2`, the unconditional upper quantile *is* `+∞`, and the finite `high`
being returned is the endpoint of a different distribution — `Λ` conditioned on
the `d+2` code failing at least once.

**Scale on the committed data.** Measured over `data/{erasure_r50,erasure_r98,
baseline_pauli}.csv` at seed 0, 70 of the 520 `(decoder, r_e, d→d+2, p)` pairs
drop at least one replicate, and 43 of them drop **all 2000** (the `d+2` cell
logged zero errors, so `Λ` and its CI are already `NaN`). The count is itself
mildly seed-dependent at the margin — 71 / 43 at seed 12345. The instructive cases are the
partial drops. Worst example, `herald_mwpm`, `r_e = 0.98`, `d = 3→5`,
`p = 0.001` — 25 errors in 100k shots at `d = 3` against 1 error in 100k at
`d = 5`:

| quantity | value |
|---|---|
| replicates dropped | 745 / 2000 (37.2%) |
| reported CI (conditioned) | `[10.0, 56.1]`, point `Λ = 41.7` |
| honest unconditional 97.5th pct | `+∞` (37.2% ≫ the 2.5% upper tail) |

**Effect on the committed figures: none.** The `≥ 50`-observed-error gate in
`figure_lambda` (and the same gate in `scripts/ablation_table.py`) admits 425
pairs, and **every one of them drops zero replicates** — verified at two
different seeds. `figures/lambda_vs_p.png` regenerates byte-identical after
this change (md5 `442e7f13…` before and after). So this is a rigor defect in a
public API, not a correction to a quoted result; no headline number moves.
`scripts/ablation_table.py` needed no wiring — it derives its ratio bounds from
Wilson endpoints, not from this bootstrap.

**What changed.** `lambda_factor` now returns a `LambdaEstimate` (a subclass of
`Estimate`, so `.value/.low/.high` and existing callers are untouched) carrying
`n_boot`, `n_boot_dropped`, and the `ci` level, plus a `ci_conditioned`
property that is true exactly when the discards outweigh the tail the interval
claims to measure — the `Λ` analogue of `FitResult.resolved`. `figure_lambda`
counts conditioned points and, if any are ever plotted, says so in the figure's
corner note instead of drawing the bar unremarked; on the current data that
count is 0 and the note is unchanged.

**What deliberately did not change.** The estimator still *excludes* the
zero-denominator replicates rather than reweighting them or propagating `+∞`.
Reporting `high = +∞` for those points is arguably the more honest estimator,
but changing it is a change to a published quantity, not a bookkeeping fix, and
the gate means it would alter nothing currently plotted. It is recorded in
`docs/FUTURE_WORK.md` rather than done here.

## The marginal CI upper endpoints are bound pins (2026-09-06)

Follow-up to "Reproducibility and environment sensitivity". That section reported
the r_e = 0.98 blind upper endpoint as grid-pinned and treated the r_e = 0.5 rows
as ordinary bootstrap spread. **That was wrong for the r_e = 0.5 blind row, which
is published**, and the mechanism proposed for r_e = 0.98 was also wrong.

### Observation

30-seed sweep, `n_boot = 1000`, `d_min = 7`, env B (Python 3.14.6, NumPy 2.5.1,
SciPy 1.18.0):

| row | upper endpoint min | max | on grid (to ~1 ULP)? |
|---|---|---|---|
| r_e = 0 herald | 1.456984% | 1.467799% | max only |
| r_e = 0 blind | 1.562420% | 1.601264% | **no** (nearest 6e-4 away) |
| r_e = 0.5 herald | 2.704709% | 2.833190% | **no** (nearest 1e-3 away) |
| r_e = 0.5 blind | 1.800000% | 2.154435% | **both ends** |
| r_e = 0.98 blind | 5.500000% | 8.500000% | **both ends** |

The r_e = 0.5 blind endpoint takes only **four distinct values over 30 seeds**:
1.800000% (25 seeds), 1.817722% (1), 1.826583% (1), 2.154435% (3). `seed = 0` --
the value published everywhere in this repo -- is one of the three that give
2.154435%. Individual bootstrap draws *do* exceed it (max draw 3.14%), so this is
not a hard cap on the draws; the **percentile** is quantised onto discrete atoms.

### Mechanism (two wrong guesses first)

The first hypothesis was that `_weighted_fit_arrays` pins `p_th` at the window
**maximum** (`p_arr.max()`). Falsified: at seed 0 exactly **1 of 998** replicates
pins at the window max.

The mechanism is the **lower** edge. `estimate_crossing` on a resampled blind
replicate sometimes runs away into the saturated tail; `_select_window` then
builds a window whose lower edge is a high grid point; `curve_fit` drives `p_th`
down onto that lower bound and pins there. At seed 0: **74 lower-edge pins vs 1
upper-edge pin**, 7.5% of replicates.

**Falsification test.** If the endpoint is set by these pins, excluding
edge-pinned replicates must move it off the atom. It does:

| row | 97.5th pct, all replicates | excluding edge pins | shift |
|---|---|---|---|
| r_e = 0.5 blind (seed 0) | 2.154435% | 1.601044% | +0.553 pp |
| r_e = 0.5 blind (seed 1) | 1.800000% | 1.591298% | +0.209 pp |
| r_e = 0.98 blind (seed 0) | 6.346955%† | 3.974615%† | +2.372 pp |
| r_e = 0.5 **herald** (seed 0) | 2.757513% | 2.761549% | −0.004 pp |

† **Not reproducible from this tree (2026-09-07).** Both r_e = 0.98 values are
kept above as published; neither regenerates here. Re-running the probe on
`data/erasure_r98.csv`, blind, `d_min=7`, `window_factor=1.5`, `n_boot=1000`,
`seed=0` (Python 3.14.6 / NumPy 2.5.1 / SciPy 1.18.0) gives:

| column | published | this tree |
|---|---|---|
| 97.5th pct, all replicates | 6.346955% | **6.346751%** |
| excluding edge pins | 3.974615% | **3.457814%** |

The discrepancy is **unexplained**. One bounded reproduction pass was spent and
found nothing: sweeping `n_boot` ∈ {200, 500, 1000, 2000}, `window_factor` ∈
{1.3, 1.5, 2.0} and `d_min` ∈ {unset, 5, 7} reproduces neither value at any
combination, and neither does varying the *definition* of "edge pin" (excluding
lower-edge pins only → 6.719513%; upper-edge only → 4.048016%; edge **and** `nu`
bound pins → 3.180475%). The closest approach to 3.974615% anywhere in that scan
is 4.048016% (upper-edge pins only), which is not a match. The search was
stopped there rather than continued.

Note the shape, which is *not* "the whole row drifted". The **all-replicates**
column reproduces exactly for all three r_e = 0.5 rows (2.154435%, 1.800000%,
2.757513%) and misses only in the 4th decimal for r_e = 0.98. The
**excluding-edge-pins** column reproduces exactly *nowhere*: the r_e = 0.5 rows
land at 1.601160% (vs 1.601044%), 1.592775% (vs 1.591298%) and 2.761548% (vs
2.761549%) — agreeing to 3–4 significant figures, close enough to be the same
computation under a small tolerance or environment difference — while the
r_e = 0.98 row is off by 0.52 pp, far outside that band. So the r_e = 0.98
entries look like two numbers carried from different runs, not one drifted run.
**The conclusions of this section are unaffected**: the sign and rough size of
the shift (+2.4 pp, and +2.9 pp on this tree) are what the falsification test
turns on, and both are reproduced.

The herald row is the control: 1.0% edge pins and no movement. The mechanism is
specific to the blind fits, whose crossing is bistable.

### This is the same phenomenon as the Phase-4 discards

"Are the discarded replicates benign?" reports that at r_e = 0.5 the dominant
discard guard is `blind:bound_pin` (73 of 102 at seed 0) and that the discards
are not missing-at-random (KS p = 3.6e-30, failed-blind `p_th` clusters high). It
already noted the crossing "shifts by grid steps under resampling" but did not
connect that to the marginal CI.

It is one mechanism handled **oppositely in two code paths**:

- `bootstrap_threshold_difference` calls `fit_threshold(..., n_boot=0)` per
  replicate and keeps only `converged` ones; `converged = not pinned`, so
  bound-pinned replicates are **discarded** (they become `n_paired_failed`).
- `_bootstrap_pipeline`, which builds the *marginal* CI, applies **no pin check**
  -- it appends `core.popt[0]` for every replicate whose `curve_fit` returns. So
  the same pins are **retained**, and they set the upper endpoint.

Measured rates agree (73/1000 discarded as `blind:bound_pin` in the Δ bootstrap;
74/998 lower-edge pins in the marginal bootstrap) though the RNG streams differ,
so these are the same mechanism at the same rate, not literally the same
replicates.

### Consequence

No verdict moves: `converged` and `resolved` are True for all 30 seeds
(`rel_ci_width` 0.25–0.49, well inside the 1.0 guard), Δ is a separate difference
bootstrap that discards these replicates, `excludes_zero` is 10/10, and every
point estimate is RNG-free. The published `[1.4, 2.2]` is the **conservative**
(wider) end of the quantisation -- the error is not that the interval is too
narrow. But its upper endpoint must not be read as "the fit supports a blind
threshold as high as 2.2%": the modal 30-seed interval is `[1.4, 1.8]`. README,
`docs/DEVLOG.md` and this file now say so at each point of quotation.

Reproduce: `uv run python scripts/paired_separation.py --full-stability` for the
10-seed table; the 30-seed enumeration and the falsification test are the probes
described above.

## The bootstrap pin asymmetry is load-bearing (2026-09-07)

Follow-up to "The marginal CI upper endpoints are bound pins". That section
established that the two bootstrap code paths handle bound-pinned replicates
oppositely, and concluded "no verdict moves". That conclusion was measured for
r_e = 0.5, where it holds. **It does not hold for r_e = 0.98 blind**, whose
`resolved = False` verdict rests on the asymmetry.

### Observation

`data/erasure_r98.csv`, `blind_mwpm`, `d_min = 7`, `window_factor = 1.5`,
`n_boot = 1000`, seeds 0-29, env B (Python 3.14.6, NumPy 2.5.1, SciPy 1.18.0).
Point estimate `p_th = 1.636309%` (RNG-free; it is the `rel_ci` denominator and
is identical in both columns). "Pins excluded" applies the same test
`bootstrap_threshold_difference` applies — `converged = not pinned`, i.e. `p_th`
at a window edge or `nu` at a `_NU_BOUNDS` bound.

| seed | replicates | pinned | pin frac | rel_ci (pins kept) | resolved | rel_ci (pins excluded) | resolved |
|---|---|---|---|---|---|---|---|
| 0 | 1000 | 199 | 0.199 | 2.9515 | no | **1.0186** | no |
| 1 | 1000 | 207 | 0.207 | 3.3771 | no | 1.0255 | no |
| 2 | 1000 | 201 | 0.201 | 3.7229 | no | 0.9462 | **yes** |
| 3 | 1000 | 203 | 0.203 | 2.5235 | no | 1.2164 | no |
| 4 | 1000 | 207 | 0.207 | 4.2583 | no | 1.0907 | no |
| 5 | 1000 | 211 | 0.211 | 3.7162 | no | 0.9831 | **yes** |
| 6 | 1000 | 209 | 0.209 | 3.1830 | no | 1.0217 | no |
| 7 | 1000 | 214 | 0.214 | 2.4324 | no | 1.1578 | no |
| 8 | 1000 | 209 | 0.209 | 3.3872 | no | 1.0740 | no |
| 9 | 1000 | 187 | 0.187 | 3.7203 | no | 0.6255 | **yes** |
| 10 | 1000 | 211 | 0.211 | 3.7151 | no | 1.1639 | no |
| 11 | 1000 | 199 | 0.199 | 2.4344 | no | 0.3881 | **yes** |
| 12 | 1000 | 182 | 0.182 | 2.4263 | no | 0.9182 | **yes** |
| 13 | 1000 | 218 | 0.218 | 3.8465 | no | 1.0551 | no |
| 14 | 1000 | 192 | 0.192 | 3.5302 | no | 1.1527 | no |
| 15 | 1000 | 208 | 0.208 | 3.7148 | no | 1.0426 | no |
| 16 | 1000 | 230 | 0.230 | 3.7128 | no | 0.9671 | **yes** |
| 17 | 1000 | 215 | 0.215 | 4.2629 | no | 1.0413 | no |
| 18 | 1000 | 219 | 0.219 | 3.0807 | no | 1.3581 | no |
| 19 | 1000 | 221 | 0.221 | 2.4322 | no | 1.1758 | no |
| 20 | 1000 | 206 | 0.206 | 3.1173 | no | 0.8204 | **yes** |
| 21 | 1000 | 206 | 0.206 | 2.4324 | no | 1.0604 | no |
| 22 | 1000 | 209 | 0.209 | 2.6092 | no | 1.0534 | no |
| 23 | 1000 | 218 | 0.218 | 2.5661 | no | 1.3237 | no |
| 24 | 1000 | 201 | 0.201 | 3.7165 | no | 0.9716 | **yes** |
| 25 | 1000 | 207 | 0.207 | 3.1663 | no | 1.0609 | no |
| 26 | 1000 | 189 | 0.189 | 3.1654 | no | 1.0393 | no |
| 27 | 1000 | 182 | 0.182 | 2.4330 | no | 1.2696 | no |
| 28 | 1000 | 200 | 0.200 | 2.4312 | no | 1.0819 | no |
| 29 | 1000 | 232 | 0.232 | 3.1598 | no | 0.9792 | **yes** |

Summary over the 30 seeds:

| quantity | pins retained (today) | pins excluded |
|---|---|---|
| `rel_ci_width` range | 2.4263 – 4.2629 | 0.3881 – 1.3581 |
| `rel_ci_width` median | 3.1658 | 1.0480 |
| `resolved` (guard = 1.0) | **0 of 30** | **9 of 30** |
| bound-pin fraction | 18.2% – 23.2%, mean 20.6% | — |

No replicate failed at any seed (`n_boot_failed = 0`, 1000/1000 fit), so the
20.6% is a pin rate, not a fit-failure rate.

### ⚠ The failure is invisible at the seed a maintainer would test

At `seed = 0` — the seed this repo declares everywhere — the pins-excluded
`rel_ci_width` is **1.0186**, marginally **above** the 1.0 `_MAX_REL_CI_WIDTH`
guard. The verdict does **not** flip there. Anyone who makes the two bootstraps
consistent as an obvious cleanup, checks seed 0, and sees `resolved = False`
unchanged will conclude the change was inert. It is not: the verdict flips at
**9 of 30** seeds (2, 5, 9, 11, 12, 16, 20, 24, 29). Seed 0 clears the guard by
0.0186 -- the **second-smallest margin of all 30 seeds** (only seed 5, at 0.0169
on the resolving side, is closer). Seed 0 is close to the worst possible single
seed to test this change on.

### Reading

The honest reading survives under **either** estimator, but for a different
reason than the current one states. With pins retained, `rel_ci` is 2.4–4.3 and
the fit misses the guard by a factor of 2.4–4.3 — which reads as a settled,
comfortable non-result. With pins excluded, `rel_ci` clusters on the guard
itself (median 1.05, range 0.39–1.36) and the verdict is a coin flip in the
seed: 9 of 30 resolving is **instability, not resolution**. Under both
estimators, the r_e = 0.98 threshold is unresolvable from this data. What the
current wide margin misrepresents is *how* unresolvable: it makes a marginal,
seed-dependent call look like a decisive one.

### Not changed

The estimator is untouched. This section records the trap; it does not fix it.
Making the two paths consistent is a change to a published verdict's mechanism
and is listed in `docs/FUTURE_WORK.md`, not done here. The same warning is in
`src/erasure_qec/analysis/threshold_fit.py` at all three code sites a maintainer
would edit from: the "DELIBERATE ASYMMETRY" block in `_bootstrap_pipeline`, the
`converged = not pinned` line in `fit_threshold`, and the
`if hf.converged and bf.converged` discard in `bootstrap_threshold_difference`.

Reproduce: re-run `_bootstrap_pipeline`'s replicate loop over seeds 0-29 with
each replicate classified by `fit_threshold`'s pin test, and take the
`_CI_PERCENTILES` interval over each subset, divided by the point-fit `p_th`.
The pins-retained column is `fit_threshold(blind, d_min=7, n_boot=1000,
seed=s).rel_ci_width` verbatim (spot-checked equal at seeds 0, 2, 11).
