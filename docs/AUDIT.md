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

**4. Directional sensitivity** (Task 1e; imputation stated). Impute the 102
discards by resampling from the **top (least-negative) decile of the observed Δ
draws** (Δ ≥ −0.73%, which reaches up to +0.38%) — a plausible-pessimistic case,
not worst-case imputation at Δ = 0 (which, above a 2.5% failure rate, forces the
97.5th percentile into the imputed block and can never exclude zero, proving
nothing):

| imputation | r_e=0.5 imputed 95% CI | excludes 0? |
|---|---|---|
| none (observed) | [−1.24, −0.66] | yes |
| top decile of Δ (plausible pessimistic) | [−1.23, −0.64] | **yes** |
| all discards at the single max observed Δ (+0.38%) | [−1.23, +0.38] | no (degenerate worst case) |

**Verdict.** The conditioning is **not** benign at r_e=0.5 — the discards are not
missing-at-random and cluster at high blind `p_th`. **The `r_e = 0.5` claim
nonetheless survives, unchanged in value but with an added caveat:** the
separation is significant *conditional on convergence* (Δ = −0.83%, CI
[−1.24, −0.66]), 10.2% of replicates are discarded and are **not** missing-at-
random, and the significance is **robust to imputing the discards from the
least-negative decile of Δ** ([−1.23, −0.64] still excludes zero). It fails only
under the degenerate worst case (all discards at the most-positive observed Δ),
which is not evidence of non-significance. The r_e=0 control is clean: its
discards are missing-at-random and its imputed CI still includes zero. The
unpaired difference bootstrap is **conservative** (pairing could only have
narrowed the interval, not widened it), so clearing zero under it is if anything
a stronger result — and the deterministic ablation, which needs no fit and no
convergence filtering, remains the primary evidence.
