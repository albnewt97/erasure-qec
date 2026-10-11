# Claim ledger

One row per numeric claim in `README.md` at `da47e3a` (C rows), then numeric claims in
`experiments/configs/*.yaml` comments (Y rows). Reference years, volumes and pages are skipped.

Conventions:

- **Line** is the README (or YAML) line range the claim spans.
- **Claim** is verbatim, at most 200 characters, with line breaks, runs of spaces, blockquote
  markers (`> `) and YAML comment markers (`# `) collapsed to a single space. A `|` inside a quoted table row is escaped as `\|`. A claim is cut at
  a clause boundary, never paraphrased; a long sentence is split across consecutive rows.
- **Producer** is the script or command that `README.md` or `docs/AUDIT.md` names for that number.
  `NONE` means neither document names one. Where a row mixes numbers with different producers, each
  part is listed. "external" marks a number quoted from a cited paper.
- **Status** is `open` for every row; later phases update it.

## README.md

| ID | Line | Claim (verbatim) | Producer | Status |
|---|---|---|---|---|
| C001 | 12-13 | a heralded site has conditional error probability ½, so `ln((1−p)/p) = 0` | NONE | open |
| C002 | 23-24 | **Blind** decoding (discard the herald bits) barely moves the threshold from baseline to `r_e = 0.5`: **1.44% → 1.49%** (95% CI `[1.4,1.6]` → `[1.4,2.2]`). | `uv run python scripts/paired_separation.py` (section 2) | open |
| C003 | 24-25 | That `2.2` upper endpoint is **not an interval estimate**: it is set by bound-pinned bootstrap replicates and lands on a discrete atom | NONE (AUDIT: 30-seed probe "described above", not committed) | open |
| C004 | 25-26 | `seed = 0` is one of only 3 of 30 seeds giving it, and 25 of 30 give `1.8` | NONE (AUDIT: 30-seed probe, not committed) | open |
| C005 | 27-28 | The point estimates (`1.44%`, `1.49%`) carry no RNG and are exact. | `uv run python scripts/paired_separation.py` (section 2) | open |
| C006 | 29 | **Herald-conditioned** decoding lifts it: **1.38% → 2.32%** (`[1.3,1.5]` → `[2.2,2.8]`) | `uv run python scripts/paired_separation.py` (section 2) | open |
| C007 | 29-31 | a **significant** separation (`Δ = blind − herald = −0.83%`, 95% CI `[−1.2,−0.7]`, excludes zero; see below) | `uv run python scripts/paired_separation.py` (section 3+4) | open |
| C008 | 32-33 | a measured **41×** at `d = 7` and a **> 103×** lower confidence bound at `d = 9`, both at `r_e = 0.98` | `uv run python scripts/ablation_table.py` | open |
| C009 | 33-34 | The larger `d = 11` bound rests on 3 herald errors and is weak evidence | `uv run python scripts/ablation_table.py` | open |
| C010 | 36 | At **`r_e = 0.98`** the collapse fit resolves **neither** threshold and the code says so (`resolved=False`) | `uv run python scripts/paired_separation.py --full-stability` (README L45-46) | open |
| C011 | 37-39 | the blind fit converges but its 95% CI runs from ~1.5% to an upper endpoint that lands anywhere in **5.5–8.5%** depending on the bootstrap seed (6.3% at `seed = 0`) | `uv run python scripts/paired_separation.py --full-stability` | open |
| C012 | 39 | 2.4–4.3× wider than the threshold itself, so it is not a resolved value | `uv run python scripts/paired_separation.py --full-stability` | open |
| C013 | 40-41 | **three of the ten seeds put it exactly on a sweep grid point** (`p = 0.055, 0.076, 0.085`) | `uv run python scripts/paired_separation.py --full-stability` | open |
| C014 | 42 | *both ends of the 5.5–8.5% span are themselves grid pins* | `uv run python scripts/paired_separation.py --full-stability` | open |
| C015 | 46 | `uv run python scripts/paired_separation.py --full-stability` (~4.7 min) | NONE (timing) | open |
| C016 | 50 | one curve per `d ∈ {3,5,7,9,11}`, one panel per `r_e` | `uv run python -m erasure_qec.analysis.plotting --data-dir data --figures-dir figures` (README L389) | open |
| C017 | 69 | 95% bootstrap percentile CI, `n_boot = 1000`, `seed = 0`, endpoints to **one decimal** | `uv run python scripts/paired_separation.py` (section 2, CI@1000) | open |
| C018 | 77 | \| 0.0 \| **1.38%** `[1.3,1.5]` (ν=1.45) \| **1.44%** `[1.4,1.6]` (ν=2.12) \| 1.02 / 1.15 \| control: no heralds → decoders agree within CIs \| | `scripts/paired_separation.py` (p_th, CI); `scripts/audit_checks.py` (herald ν=1.45; χ²/dof 1.02 / 1.15); blind ν=2.12: NONE | open |
| C019 | 78 | \| 0.5 \| **2.32%** `[2.2,2.8]` (ν=1.90) \| **1.49%** `[1.4,2.2]`† (ν=1.58) \| 0.62 / 0.47 \| herald above blind — separation **significant**, see Δ \| | `scripts/paired_separation.py` (p_th, CI); `scripts/audit_checks.py` (χ²/dof 0.62 / 0.47); ν=1.90, ν=1.58: NONE | open |
| C020 | 79 | \| 0.98 \| **not resolved** (ν-pinned) \| **not resolved** (CI `[1.5, 5.5–8.5]`, upper endpoint seed-dependent and grid-pinned — see above) \| — / 0.20 \| | `scripts/paired_separation.py` (verdicts, seed-0 CI); `--full-stability` (5.5–8.5 span); χ²/dof 0.20: NONE | open |
| C021 | 81-82 | † **The `r_e = 0.5` blind upper endpoint is not a fitted crossing, and `2.2` is an unlucky seed.** | NONE (AUDIT: 30-seed probe, not committed) | open |
| C022 | 82-84 | on ~7.5% of bootstrap replicates the resampled blind crossing estimate runs into the saturated tail | NONE | open |
| C023 | 85-87 | Over 30 seeds the endpoint takes **four distinct values**, and **25 of 30 seeds give 1.80%** — `seed = 0`, which this table quotes, is one of only **3 of 30** that give 2.154%. | NONE (AUDIT: 30-seed probe, not committed) | open |
| C024 | 87 | Excluding edge-pinned replicates moves it from 2.154% to 1.601%. | NONE (AUDIT: falsification-test probe, not committed) | open |
| C025 | 88 | So the modal interval at this precision is `[1.4, 1.8]`, not `[1.4, 2.2]` | NONE | open |
| C026 | 90-91 | The `p_th` point estimate (1.49%), ν and χ²/dof are unaffected — they carry no RNG — and no verdict changes (`converged` and `resolved` are True for all 30 seeds). | NONE (30 seeds); `--full-stability` covers 10 | open |
| C027 | 97 | **Δ = −0.83%, 95% CI `[−1.2, −0.7]`, excludes zero** | `uv run python scripts/paired_separation.py` (section 3+4) | open |
| C028 | 97-98 | stable across seeds (all of seeds 0–9 exclude zero — | `uv run python scripts/paired_separation.py --full-stability` | open |
| C029 | 99 | (control `r_e = 0` gives Δ consistent with zero) | `uv run python scripts/paired_separation.py` (section 3+4) | open |
| C030 | 102-103 | the ~10% discards are **not missing-at-random** (they cluster at high blind `p_th`) | `uv run python scripts/paired_separation.py` (section 6) | open |
| C031 | 104-105 | the upper endpoint reaches zero only if ≥ 24 of the 102 discards had `Δ ≥ 0` (seed 0; 24–26 of 77–115 across seeds 0–9, from | `scripts/paired_separation.py` (seed 0); `--full-stability` (seeds 0–9) | open |
| C032 | 106 | while their recorded `p_th` imply `Δ ≥ 0` for only 2 | `uv run python scripts/paired_separation.py` (section 6, [T2b]) | open |
| C033 | 112 | Wu et al. [[3]](#references) report thresholds rising **0.937% → 4.15%** at 98% conversion. | NONE (external) | open |
| C034 | 138 | Their Figs. 3–4 set `p_m = 0` and "neglect idle errors" | NONE (external) | open |
| C035 | 140-141 | at `r_e = 0.98` only ~45% of our DEM error mass is heralded (below) | `uv run python scripts/heralded_fraction.py` | open |
| C036 | 143-144 | Their thresholds come "from the crossing of `d = 11` and `d = 15`"; ours from a `d ≥ 7` collapse fit over `d ∈ {7, 9, 11}`. | NONE (external; ours: NONE) | open |
| C037 | 147-148 | a `d = 5, T = 5, p = 0.01` circuit emits 400 `DEPOLARIZE2(p)` two-qubit-gate locations (the same budget as theirs) *plus* 470 locations they do not have | NONE | open |
| C038 | 148-149 | 180 `DEPOLARIZE1(p)` idle and 290 `X_ERROR(p)` measure/reset. | NONE | open |
| C039 | 149-150 | we report **1.38% against their 0.937%** | `scripts/paired_separation.py` (1.38%); 0.937%: NONE (external) | open |
| C040 | 157-158 | **1.708%** (all `d`) → **1.509%** (`d ≥ 5`) → **1.376%** (`d ≥ 7`) → **1.300%** (`d ≥ 9`, unconverged, two distances) | `scripts/audit_checks.py` (1.708%, 1.376%); 1.509%, 1.300%: NONE | open |
| C041 | 158-159 | and it has not flattened at the top of our range, so a `d = 11`/`d = 15` crossing would sit lower still | NONE | open |
| C042 | 159 | That plausibly accounts for part of the 0.44-point gap, not obviously all of it. | NONE | open |
| C043 | 160-161 | their SI reports `p_m = p` only at `R_e = 0.98` (threshold falls 4.15% → 2.85(1)%) | NONE (external) | open |
| C044 | 167-168 | **Erasure conversion alone** (blind) moves the threshold `1.44% → 1.49%` from `r_e = 0` to `0.5` — close to a wash. | `uv run python scripts/paired_separation.py` (section 2) | open |
| C045 | 168 | **Herald-conditioned decoding** is where the benefit lives: `1.38% → 2.32%` | `uv run python scripts/paired_separation.py` (section 2) | open |
| C046 | 170-172 | the old `1.6×`-blind / `2.7×`-herald decomposition was an artifact of the pre-fix budget-shrinking model and does not survive | NONE | open |
| C047 | 183 | \| `DEPOLARIZE2` \| `p·(1 − r_e)` \| on `(a, b)` — residual Pauli \| | `tests/test_noise_model.py::test_error_budget_invariance` (named at README L188) | open |
| C048 | 184 | \| `HERALDED_ERASE` \| `(2/3)·p·r_e` \| on `a`, and independently on `b` — appends a herald bit \| | `tests/test_noise_model.py::test_error_budget_invariance` (named at README L188) | open |
| C049 | 185 | \| `X_ERROR` / `X_ERROR` / `DEPOLARIZE1` \| `p_meas` / `p_reset` / `p_idle` (`= p`) \| before `M`/`MR`, after `R`, on idle qubits \| | NONE | open |
| C050 | 187-188 | The `(2/3)·p·r_e` rate (an erasure is non-identity only ¾ of the time, so `2·q·¾ = p·r_e`) holds the per-gate non-identity budget at `p` for **every** `r_e` | `tests/test_noise_model.py::test_error_budget_invariance` (named at README L188) | open |
| C051 | 199-201 | Measured by `scripts/heralded_fraction.py` at `d = 5, p = 0.02`: gate-only **0.23 / 0.45** at `r_e = 0.5 / 0.98` (0.32 / 0.63 with `convert_idle`). | `uv run python scripts/heralded_fraction.py` | open |
| C052 | 201-202 | The committed sweeps use gate-only, so `r_e = 0.98` heralds **45%** of the mass, not 98%. | `uv run python scripts/heralded_fraction.py` | open |
| C053 | 204-206 | Counting them (the pre-2026-09-05 definition, still available as `count_no_op_heralds=True`) gave 0.30 / 0.55 and 0.42 / 0.72 — i.e. it overstated conversion. | `uv run python scripts/heralded_fraction.py` (legacy columns) | open |
| C054 | 221-223 | Shots are sampled with probability-1 probe erasures on `(1,3)`/`(1,5)` (`probe_q = 1.0`), but **both decoders are compiled from a circuit with `probe_q = 0.02`** | "M6 test" (README L221; file not named) | open |
| C055 | 229-230 | Blind fails on the ~1/4 of shots where both probes take an X-like Pauli — the test asserts 15–35% of 2048 shots — while herald fails 0. | "M6 test" (README L221; file not named) | open |
| C056 | 236-237 | `p_L(blind) / p_L(herald)` at fixed sub-threshold `p = 1.0%`, by distance (`scripts/ablation_table.py`, 100 000 shots, seed 0; ~46 min on one core). | `uv run python scripts/ablation_table.py`; timing: NONE | open |
| C057 | 237-239 | cells below the **≥ 50 observed-error gate** the Λ figure uses are reported as a lower confidence bound `> N×`, not a point ratio. | `uv run python scripts/ablation_table.py` | open |
| C058 | 243 | \| 3 \| 1.3× \| 2.1× \| | `uv run python scripts/ablation_table.py` | open |
| C059 | 244 | \| 5 \| 1.9× \| 8.8× \| | `uv run python scripts/ablation_table.py` | open |
| C060 | 245 | \| 7 \| 2.8× \| 41.2× \| | `uv run python scripts/ablation_table.py` | open |
| C061 | 246 | \| 9 \| 4.1× \| **> 103×** † \| | `uv run python scripts/ablation_table.py` | open |
| C062 | 247 | \| 11 \| 6.4× \| **> 339×** †‡ \| | `uv run python scripts/ablation_table.py` | open |
| C063 | 249 | (herald errors: `d = 9` → 23, `d = 11` → 3 of 100 000 shots at `r_e = 0.98`) | `uv run python scripts/ablation_table.py` | open |
| C064 | 251-254 | Each is an endpoint of a two-sided 95% Wilson interval — a one-sided 97.5% limit — so Bonferroni makes the quoted value a **≥ 95% one-sided lower bound on that cell's ratio**. | NONE (derivation in AUDIT "Bound correction") | open |
| C065 | 255-256 | (Before 2026-09-05 this column divided blind's *point estimate* by herald's Wilson upper limit, which bounded nothing; the corrected values are ~3% lower.) | `scripts/ablation_table.py` (corrected values); old values: NONE (pre-fix code) | open |
| C066 | 258-259 | **The `d = 11` cell rests on 3 herald errors in 100 000 shots — treat `> 339×` as weak evidence, not a measurement.** | `uv run python scripts/ablation_table.py` | open |
| C067 | 260-261 | the herald shot-rate CI there spans `[1.0e-5, 8.8e-5]`, nearly an order of magnitude. | `uv run python scripts/ablation_table.py` | open |
| C068 | 262 | The point ratio at that cell is 1034×, which is precisely why no point ratio is quoted. | `uv run python scripts/ablation_table.py` | open |
| C069 | 262-263 | The `d = 9` cell (23 errors) is also below the gate but far less fragile. | `uv run python scripts/ablation_table.py` | open |
| C070 | 265-266 | every `r_e = 0.5` cell and the `r_e = 0.98` cells up to `d = 7` clear the 50-error gate | `uv run python scripts/ablation_table.py` | open |
| C071 | 276 | `Λ > 1` below threshold means distance helps | `uv run python -m erasure_qec.analysis.plotting --data-dir data --figures-dir figures` (README L389) | open |
| C072 | 277-278 | larger `r_e` gives larger `Λ` | `uv run python -m erasure_qec.analysis.plotting --data-dir data --figures-dir figures` (README L389) | open |
| C073 | 279-280 | Points are shown only where **both** distances have ≥ 50 observed errors; the low-statistics tail is dropped. | `uv run python -m erasure_qec.analysis.plotting --data-dir data --figures-dir figures` (README L389) | open |
| C074 | 282-284 | `figure_lambda` requires ≥ 50 errors at *both* `d` and `d+2`; `ablation_table._fmt_ratio` gates on the **herald** count alone | NONE | open |
| C075 | 285-286 | On the committed sweeps the two coincide — blind counts run 3016–6290, far above 50 | `scripts/ablation_table.py` (AUDIT raw-count table, 2026-09-05; README L288-289) | open |
| C076 | 306 | Instantaneous, 100% reliable, zero false positives, zero latency. | NONE (model assumption) | open |
| C077 | 320-321 | The large codes saturate to ½ within one grid step above the ~6% crossing | NONE | open |
| C078 | 324 | **The effective crossing drifts with distance** (+0.33% at `r_e = 0`, larger at higher `r_e`) | NONE | open |
| C079 | 326 | **ν is only weakly constrained at `d ≤ 11`** (three large distances, coarse grid) | NONE | open |
| C080 | 337-338 | the hook-safe CX schedule keeps graphlike distance tracking `d`, a deliberately broken one halves it to `⌈(d+1)/2⌉` | `uv run python -m erasure_qec.analysis.plotting --data-dir data --figures-dir figures` (README L389; `figures/hook_regression.png`) | open |
| C081 | 343-344 | Measured on an **Apple M2 (8 cores, macOS 14.5), single-threaded, otherwise idle**, `d=5, p=3e-2`, 10⁴ shots | `uv run pytest tests/test_end_to_end_montecarlo.py -k throughput -s` | open |
| C082 | 345-346 | a **median 79.8k shots/s** herald-free vs a **median 1,578 shots/s** at `r_e = 0.98` — a **median 50.8× gap**. | `uv run pytest tests/test_end_to_end_montecarlo.py -k throughput -s` (single run per call; median over runs: NONE) | open |
| C083 | 346-347 | Over 10 consecutive runs (2026-09-06) the *observed* spreads were 77.3k–85.7k, 1,529–1,585 and 48.8–54.4× respectively (sd 2.5k, 17, 1.6). | `uv run pytest tests/test_end_to_end_montecarlo.py -k throughput -s` (single run per call; 10-run aggregation: NONE) | open |
| C084 | 348-349 | the previously published min–max (77.4k–84.2k, 1,542–1,593) was breached at both ends by this very sample | NONE | open |
| C085 | 352-353 | on this same machine under full CPU oversubscription the gap ranged 55.8–65.3× | NONE | open |
| C086 | 353-354 | so treat the ~51× as one machine's number, not a property of the decoder | `uv run pytest tests/test_end_to_end_montecarlo.py -k throughput -s` | open |
| C087 | 354-356 | Replacing the per-shot matcher rebuild with per-signature vectorized construction narrowed the gap from ~69× to ~51×; it did **not** bring it under the 50× target that motivated the work. | NONE (~69×); `tests/test_end_to_end_montecarlo.py -k throughput` (~51×) | open |
| C088 | 356-357 | (An earlier revision of this bullet claimed 48×; that was 80,000/1,650, a ratio of two independently rounded numbers rather than a measured one.) | NONE | open |
| C089 | 392 | uv run python scripts/heralded_fraction.py # circuit-wide heralded fraction vs r_e (~2 s) | NONE (timing) | open |
| C090 | 393 | uv run python scripts/paired_separation.py # Delta significance + discard diagnostics (~72 s) | NONE (timing) | open |
| C091 | 394 | uv run python scripts/paired_separation.py --full-stability # + the 10-seed stability table (~4.7 min) | NONE (timing) | open |
| C092 | 395 | uv run python scripts/audit_checks.py # docs/AUDIT.md table (pre-fix expectations vs HEAD) (~2 min) | NONE (timing) | open |
| C093 | 396 | uv run python scripts/ablation_table.py # sub-threshold suppression, direct from circuits (~46 min) | NONE (timing) | open |
| C094 | 403-404 | Verified 2026-09-06 to regenerate byte-identical to the committed PNGs (all four, incl. `lambda_vs_p.png` md5 `442e7f13…`) | `uv run python -m erasure_qec.analysis.plotting --data-dir data --figures-dir figures` (README L389) | open |
| C095 | 413 | it decodes 100k shots per cell over `d = 3..11` x `r_e = 0.5, 0.98` | `uv run python scripts/ablation_table.py` | open |
| C096 | 416-418 | rough wall-clock on 8 cores is ~10–20 min (`baseline_pauli`), ~1–2 h (`erasure_r50`), ~4–8 h (`erasure_r98`, slow-path dominated). | NONE (timing) | open |

## experiments/configs/*.yaml (comments)

| ID | Line | Claim (verbatim) | Producer | Status |
|---|---|---|---|---|
| Y001 | baseline_pauli.yaml:2-3 | p grid: 13 log-spaced points in [1e-3, 1e-1]. | NONE (config values) | open |
| Y002 | baseline_pauli.yaml:11-12 | Additive fit-resolution band: linearly spaced p-points straddling the ~1.5% Pauli crossing. | NONE | open |
| Y003 | baseline_pauli.yaml:12-13 | The 13-point log grid only lands ~3 points in the fit window, leaving nu under-constrained; | NONE | open |
| Y004 | baseline_pauli.yaml:13-14 | these densify the crossing so the finite-size scaling fit resolves nu ~ 1.5. | NONE | open |
| Y005 | baseline_pauli.yaml:20 | Fixed sub-threshold p for the Lambda scan (~0.5%, below the ~1% Pauli threshold). | NONE | open |
| Y006 | erasure_r50.yaml:1-2 | Threshold sweep, half-converted budget: r_e = 0.5 of the two-qubit error budget becomes heralded erasure, the rest stays unheralded depolarizing. | NONE (config values) | open |
| Y007 | erasure_r50.yaml:11-12 | Additive fit-resolution band bracketing the MEASURED r_e=0.5 crossing (~2.2%, from the d>=7 collapse fit). | NONE | open |
| Y008 | erasure_r50.yaml:12-13 | The 13-point log grid lands only ~1-3 points in the fit window, leaving nu under-constrained; | NONE | open |
| Y009 | erasure_r50.yaml:21 | Fixed sub-threshold p for the Lambda scan (below the r_e=0.5 threshold). | NONE (refers to `lambda_p: 0.015`, L22) | open |
| Y010 | erasure_r98.yaml:1-2 | Threshold sweep, erasure-dominated regime: r_e = 0.98 (biased-erasure hardware limit; threshold expected to slide toward ~4-5%). | NONE | open |
| Y011 | erasure_r98.yaml:12-13 | First band (0.030-0.050) was a prediction extrapolated from the measured r_e=0 (~1.4%) and r_e=0.5 (~2.2%) crossings | NONE | open |
| Y012 | erasure_r98.yaml:13-14 | it UNDERSHOT under the old (budget-shrinking) noise model, where the crossing sat in [0.050, 0.068]. | NONE | open |
| Y013 | erasure_r98.yaml:15-16 | The corrected constant-budget model (2/3*p*r_e herald rate) RAISES the erasure rate ~33% at r_e=0.98, shifting the herald crossing higher again, | NONE | open |
| Y014 | erasure_r98.yaml:16-17 | so a third band (0.067-0.085) brackets the fixed-model crossing. | NONE | open |
| Y015 | erasure_r98.yaml:18-19 | blind_mwpm crosses much lower, which the log grid + lower bands already resolve. | NONE | open |
| Y016 | erasure_r98.yaml:27 | Fixed sub-threshold p for the Lambda scan (below the ~4-5% r_e=0.98 threshold). | NONE (refers to `lambda_p: 0.025`, L28) | open |
