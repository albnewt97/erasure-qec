# Wu d-pair decision (pre-registered, Task 6)

Written and committed before any Wu data exist (`data/wu/` does not exist at this commit).
Every number below is printed by `scripts/wu_projection.py`; the full output is
`docs/remediation/baseline_outputs/wu_projection.txt`.
Machine: Apple M2 (4 performance + 4 efficiency cores), 8.0 GiB RAM; stim 1.16.0,
pymatching 2.4.0, sinter 1.16.0.

## Decision

**d-pair (9,13) for all seven configs, plus d = 7 as a cross-check. No config is deferred.**

Rule clause applied: 2.

- Clause 1 fails for (11,15) on memory.
  - The projected total time is 2.45 h, within the 12 h limit.
  - Three configs exceed the memory limit of 4.00 GiB (50% of 8.0 GiB):
    `wu_fig3_re098` 4.13 GiB, `wu_fig3b_re100` 4.32 GiB, `wu_si_pm005` 4.13 GiB.
- Clause 2 holds for (9,13).
  - The projected total is 1.28 h.
  - The largest per-config peak is 3.12 GiB (`wu_fig3b_re100`).

## Pre-registered rule (as given)

1. Use (11,15) for all configs if the total over all seven configs is ≤ 12 h and the peak memory
   of every config is ≤ 50% of this machine's RAM.
2. Otherwise use (9,13) for all configs if the same two conditions hold.
3. Otherwise use (9,13) and defer configs, in this order, until both hold: `wu_si_pmp`, then
   `wu_si_pm005`, then `wu_fig3b_re090`.

One d-pair applies to every config that runs.

## Projection

Total wall time per config is the sum of three terms:

- **Task builds:** the sum over tasks of (sampling build + decoding build + DEM), divided by 8.
  Tasks are built in an 8-process pool, one config at a time.
- **Sampling and decoding:** the sum over tasks and decoders of shots × (sample + decode
  seconds per shot), divided by 8 workers. sinter samples separately for each decoder.
- **Compile:** the sum over tasks and decoders of (decoder compile + sampler compile). All 8
  workers compile each task concurrently, so each compile adds its own duration once.

Peak memory per config is the larger of two phases:

- **Pool phase:** base RSS + held tasks + 8 × peak RSS of building one task.
- **Collect phase:** base RSS + held tasks + 8 × worker peak RSS.

In both phases, "held tasks" counts every task of the config: sinter holds them all at once.

### d-pair (11,15), with d = 7

| config | tasks | shots | build h | shots h | compile h | total h | peak RAM GiB | held tasks GiB | workers GiB |
|---|---|---|---|---|---|---|---|---|---|
| wu_fig3_re000 | 33 | 2150031 | 0.01 | 0.00 | 0.00 | 0.02 | 3.57 | 1.03 | 2.46 |
| wu_fig3_re098 | 51 | 6676446 | 0.02 | 0.61 | 0.01 | 0.64 | 4.13 | 1.60 | 2.46 |
| wu_fig3b_re100 | 57 | 3587178 | 0.02 | 0.35 | 0.01 | 0.38 | 4.32 | 1.79 | 2.46 |
| wu_si_pm005 | 51 | 5683698 | 0.02 | 0.52 | 0.01 | 0.55 | 4.13 | 1.60 | 2.46 |
| wu_si_pmp | 45 | 5679600 | 0.02 | 0.52 | 0.01 | 0.55 | 3.94 | 1.41 | 2.46 |
| wu_fig3b_re050 | 33 | 1079562 | 0.01 | 0.09 | 0.01 | 0.11 | 3.57 | 1.03 | 2.46 |
| wu_fig3b_re090 | 33 | 1975260 | 0.01 | 0.17 | 0.01 | 0.19 | 3.57 | 1.03 | 2.46 |
| **total** | | | | | | **2.45** | max 4.32 | | |

### d-pair (9,13), with d = 7

| config | tasks | shots | build h | shots h | compile h | total h | peak RAM GiB | held tasks GiB | workers GiB |
|---|---|---|---|---|---|---|---|---|---|
| wu_fig3_re000 | 33 | 2247039 | 0.01 | 0.00 | 0.00 | 0.01 | 2.60 | 0.72 | 1.79 |
| wu_fig3_re098 | 51 | 5226618 | 0.01 | 0.29 | 0.01 | 0.31 | 2.99 | 1.12 | 1.79 |
| wu_fig3b_re100 | 57 | 2386602 | 0.01 | 0.13 | 0.00 | 0.15 | 3.12 | 1.25 | 1.79 |
| wu_si_pm005 | 51 | 5713020 | 0.02 | 0.32 | 0.01 | 0.34 | 2.99 | 1.12 | 1.79 |
| wu_si_pmp | 45 | 4460646 | 0.01 | 0.25 | 0.01 | 0.27 | 2.86 | 0.99 | 1.79 |
| wu_fig3b_re050 | 33 | 1113642 | 0.01 | 0.06 | 0.00 | 0.07 | 2.60 | 0.72 | 1.79 |
| wu_fig3b_re090 | 33 | 2218584 | 0.01 | 0.12 | 0.00 | 0.13 | 2.60 | 0.72 | 1.79 |
| **total** | | | | | | **1.28** | max 3.12 | | |

"shots" counts both decoders where a config runs both.

### Assumptions

- **Grids:** as specified for P5a. The `wu_fig3b_re050` and `wu_fig3b_re090` grids use the
  planning values c(0.5) = 0.020 and c(0.9) = 0.036. These were chosen high, which makes the
  projection conservative.
- **Decoders:** `wu_herald_mwpm` and `wu_blind_mwpm` everywhere except `wu_fig3_re000`, which
  runs `wu_herald_mwpm` only.
- **Shots:** max_shots = 100000 and max_errors = 1000.
  - Shots per (d, p) = min(100000, ceil(1000 / P_L(p))).
  - P_L is interpolated log-linearly in p between 1000-shot herald pilots at the lowest and
    highest p of each grid, at d_hi.
  - The d_hi value of P_L is used for every d in the config.
  - One pilot end had zero failures (`wu_si_pmp`, d = 15, p = 0.022). It is set to
    P_L = 0.5 / 1000. That end is capped at 100000 shots either way.
- **Per-shot cost, build time and compile time:** each is interpolated linearly in p between the
  two pilot ends, at each d separately.
- **Held memory per task:** at each d, the larger of two measures of one `wu_fig3_re098` task at
  p = 0.042:
  - the slope of peak RSS over three tasks held in one process;
  - `len(str(sampling)) + len(str(DEM))`.

  The R_e = 0.98 circuit has the most instructions (31 branches per pair). It is used for all
  configs, which overstates the R_e = 0, 0.5, 0.9 and 1.0 configs.
- **Worker peak:** the largest measured worker peak at d_hi. This was 315 MB at d = 15
  (`wu_fig3_re098`, p = 0.05) and 229 MB at d = 13. Each measurement covers a fresh process that
  unpickles a task, compiles the sampler and both decoders, and decodes 1000 shots.

### Sensitivity (does not enter the rule)

Eight concurrent processes on this machine reached a parallel efficiency of 0.61 for task builds
(d = 11) and 0.53 for decoding (d = 13). Efficiency here is the time one job takes alone divided
by the mean time of 8 concurrent jobs.

With those efficiencies, the totals become 4.54 h for (11,15) and 2.37 h for (9,13). Both stay
within 12 h, so the time condition does not decide the outcome.

The outcome is decided by memory. The (11,15) excess is 0.13–0.32 GiB above the 4.00 GiB limit.
The memory model has two terms:

- **Workers:** 8 × the measured worker peak, giving 2.46 GiB at d_hi = 15.
- **Held tasks:** 1.60–1.79 GiB in the three failing configs.

## 6a profile (one `wu_fig3_re098` task at p = 0.042, fresh process per d)

The build columns are times after the 6b-ii change.

| d | build sampling s | build decoding s | DEM s | str() s | strong_id s | len(sampling) | str sampling MB | str DEM MB | peak RSS one task MB | held MB/task |
|---|---|---|---|---|---|---|---|---|---|---|
| 7 | 1.07 | 1.39 | 0.02 | 0.02 | 0.02 | 40533 | 1.5 | 0.5 | 91 | 8 |
| 9 | 2.27 | 2.98 | 0.03 | 0.04 | 0.05 | 89157 | 3.5 | 1.1 | 89 | 13 |
| 11 | 4.22 | 5.47 | 0.07 | 0.09 | 0.10 | 166301 | 6.5 | 2.1 | 90 | 22 |
| 13 | 6.91 | 9.04 | 0.12 | 0.13 | 0.17 | 278541 | 11.0 | 3.5 | 112 | 46 |
| 15 | 11.89 | 15.72 | 0.23 | 0.22 | 0.29 | 432453 | 17.1 | 5.7 | 96 | 66 |

## Pilot P_L (herald decoder, 1000 shots; planning-only, not a result)

These pilot data were used only for the shots model above. They are not in `data/` and are not
used in any result.

| config | p range | d = 13 failures | d = 15 failures |
|---|---|---|---|
| wu_fig3_re000 | 0.0075 / 0.0125 | 2 / 70 | 3 / 68 |
| wu_fig3_re098 | 0.034 / 0.050 | 5 / 116 | 1 / 130 |
| wu_fig3b_re100 | 0.042 / 0.060 | 16 / 298 | 9 / 307 |
| wu_si_pm005 | 0.030 / 0.046 | 4 / 102 | 3 / 137 |
| wu_si_pmp | 0.022 / 0.036 | 4 / 172 | 0 / 243 |
| wu_fig3b_re050 | 0.016 / 0.025 | 24 / 244 | 24 / 272 |
| wu_fig3b_re090 | 0.0288 / 0.045 | 10 / 205 | 11 / 246 |

## Speed-ups adopted (no emitted circuit changed)

- **6b-ii:** both Wu injectors emit stim text, parsed once per CX layer.
  - `tests/test_wu_text_build.py` keeps the per-instruction `append` construction as the
    reference. It asserts `str()` equality and exact `==` for d ∈ {3, 5}, R_e ∈ {0.5, 0.98, 1.0}
    and p_m ∈ {0, 0.005}, with and without `location_overrides`.
  - The d = 15 sampling build went from 157 s (`baseline_outputs/bench_wu.txt`) to 11.9 s
    (profile above).
- **6b-i:** tasks are built in an 8-process pool, one config at a time
  (`sinter_adapter.wu_tasks_parallel`).
  - `tests/test_wu_decoder.py::test_pool_built_tasks_match_serial_tasks` asserts three things:
    the `strong_id` of a pool-built task equals that of the serially built one, `str()` is
    equal, and the circuits and DEMs a sinter worker receives by pickle are equal.

## Notes for P5a (recorded now)

- **Distances:** d ∈ {7, 9, 13}, i.e. the cross-check plus d_lo and d_hi.
- **Stopping rule:** max_shots = 100000 and max_errors = 1000 for every Wu config.
- **Decoders:** pass `custom_decoders=WU_DECODERS` to `sinter.collect`.
- **`wu_fig3_re000` config comment:** it runs `wu_herald_mwpm` only. At R_e = 0 the herald
  decoder's predictions equal plain PyMatching on the same DEM (T3, 5000 shots, d = 5).
- **Edge rule (applied in P5b):** if a d-pair crossing's bracket touches either end of a grid,
  append 3 points beyond that end at the same step and re-collect that config once.
