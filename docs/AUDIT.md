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

**Root causes 3, 4, 5 — fit pipeline (fixed, commit `0ce6ff5`).**
- Bound-pinning is now reported as `converged=False`. r50 herald `d_min=None`,
  which used to report `converged=True` pinned at `p_arr.max`, now returns
  `converged=False`.
- The fit requires ≥ 3 degrees of freedom (≥ 8 points for 5 params). r50 herald
  `d_min=7` no longer produces a (meaningless, dof≈1) fit — it returns
  `converged=False, "insufficient data"`. The committed r50 sweep's `d ≥ 7`
  grid is too sparse to fit the 5-parameter ansatz honestly.
- σ is now the ~1σ Wilson standard error (95% half-width ÷ Z_95). On the valid
  committed baseline this brings χ²/dof from ~0.2 to ~0.9–1.0 (i.e. the model
  now fits at the noise level rather than "too well"). `FitResult` carries a
  `chi2_dof` field.
- The fit window excludes coin-flip data on the *shot-level* rate; the max
  P_L_shot admitted to any baseline window dropped below 0.4.

**Re-measured baseline threshold (committed `real_baseline_pauli.csv`, fixed
pipeline), `d ≥ 7`:** herald `p_th = 1.380% ± 0.039%`, ν = 1.45 ± 0.41,
χ²/dof = 0.89; blind `p_th = 1.417% ± 0.045%`, ν = 1.93 ± 0.56, χ²/dof = 1.01.
The two agree within combined error bars, as required at `r_e = 0`.

**Reproducibility.** The stale (pre-fix) sweeps were moved to
`data/stale_old_model/` (never deleted). The erasure sweeps are being
re-collected under the fixed model; any erasure threshold that survives the
honest `d ≥ 7` pipeline is reported from a committed snapshot, and any that does
not converge is reported as a non-result. The reproducible, model-normalisation-
independent result is the herald-vs-blind ablation, computed directly from
circuits (deterministic seed, no CSV).
