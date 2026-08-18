# Herald-Conditioned Decoding of the Erasure-Converted Rotated Surface Code

This repo studies **erasure conversion** in a distance-`d` rotated surface-code memory [[1]](#references)
(Z basis, in [Stim](https://github.com/quantumlib/Stim) [[2]](#references)): a fraction `r_e` of each
two-qubit gate's error budget is turned into *heralded erasures* — the hardware flags *which* qubit was
disturbed and *when*, rather than letting the error pass silently. This models **¹⁷¹Yb neutral-atom
Rydberg gates** [[3,4]](#references) and **dual-rail superconducting cavities** [[5,6]](#references),
where leakage and photon-loss events are natively detected. A heralded error is far easier to correct
than a blind one: knowing its location lets a **herald-conditioned matching decoder**
(PyMatching [[7]](#references), wrapped as a `sinter.Decoder` [[2]](#references)) zero that edge's
weight — a heralded site has conditional error probability ½, so `ln((1−p)/p) = 0` and a correction
routes through it for free.

**Headline result (measured).** Holding the *per-gate error budget constant* across `r_e` (so `r_e`
converts errors rather than removing them — see [Noise channels](#noise-channels-5)), the improvement
from erasure comes almost entirely from **reading the herald bits**, not from the conversion itself:

- The **blind** decoder (which discards the herald bits and folds erasures into extra depolarizing
  noise) barely moves the threshold as `r_e` grows: **1.44% → 1.49% → 1.64%** at `r_e = 0 / 0.5 / 0.98`.
- The **herald-conditioned** decoder lifts it substantially: **1.38% → 2.32%** at `r_e = 0.5`, and the
  herald-vs-blind advantage compounds with code distance and erasure fraction (the
  [ablation](#the-herald-aware-vs-blind-ablation)) — reaching a large factor at `d = 11, r_e = 0.98`.

At **`r_e = 0.98`** the herald crossing sits near **~6%**, but the largest codes (`d = 9, 11`) saturate
to the ½ coin-flip limit within one grid step above it, so the `d ≥ 7` collapse fit **does not converge**
(ν rails to its bound). We report that as a **non-result**, not a number, and lead the near-full-conversion
story with the deterministic ablation instead. This corrects an earlier version of this README, whose
"6.27% ± 0.54%" headline was an artifact of a bug that reported a bound-pinned fit as converged, run on a
noise model whose per-gate budget shrank with `r_e` (see [`docs/AUDIT.md`](docs/AUDIT.md)).

![Threshold panels](figures/threshold_panels.png)

*Per-round `p_L` vs physical error rate `p`, one curve per distance `d ∈ {3,5,7,9,11}`, one panel per
erasure fraction; the dashed line is the fitted threshold (`d ≥ 7` collapse fit) with its bootstrap CI.
The `r_e = 0.98` herald panel carries no line — the fit does not converge (annotated "no fit").*

> **Data status.** These figures are generated from the **real Monte-Carlo sweeps committed in `data/`**
> (`baseline_pauli.csv`, `erasure_r50.csv`, `erasure_r98.csv`), collected under the fixed noise model.
> Every threshold quoted here is reproducible from that committed data by the code in this repo. The
> stale pre-fix sweeps are preserved (un-tracked) under `data/stale_old_model/`. Separate committed
> **synthetic fixtures** (`tests/fixtures/*.csv`) drive the deterministic analysis/plotting tests.

---

## Threshold scaling

The panels above are per-round logical error rate `p_L` vs physical error rate `p`. Each converged
panel's dashed line (with bootstrap CI band) is the fitted threshold `p_th` from the **asymptotic
`d ≥ 7` collapse fit** (see [Caveats](#caveats-and-limitations)); the inset is the collapse onto the
single scaling variable `x = (p − p_th)·d^{1/ν}`.

### Measured thresholds (`d ≥ 7` collapse fit)

| `r_e` | `herald_mwpm` `p_th` | `blind_mwpm` `p_th` | χ²/dof (h / b) | notes |
|---|---|---|---|---|
| 0.0 | **1.38% ± 0.11%** (ν = 1.45) | **1.44% ± 0.12%** (ν = 2.12) | 1.02 / 1.15 | control: no heralds, so the decoders agree within error bars (as required) |
| 0.5 | **2.32% ± 0.30%** (ν = 1.90) | **1.49% ± 0.11%** (ν = 1.58) | 0.62 / 0.47 | herald clearly above blind |
| 0.98 | **does not converge** | **1.64% ± 0.08%** (ν = 1.72) | — / 0.20 † | herald crossing ~6%, but `d = 9,11` saturate immediately above it → ν unconstrained |

All quoted `p_th ± σ` come from a seeded parametric bootstrap over the sinter counts. † The `r_e = 0.98`
blind fit has only 9 in-window points (dof = 4); its χ²/dof = 0.20 is a small-sample over-fit flag, so
read it as ~1.6% with a wide systematic, not a precise value. The `r_e = 0.98` **herald** threshold is a
genuine non-result under this pipeline: `fit_threshold(..., d_min=7)` returns `converged=False` with
`nu at upper bound`. The all-`d` fit "reaches" ~6–8.5% depending on the window, but mixes the
finite-size-corrected small distances and rests on the ragged, `max_errors`-limited data just above the
crossing, so we do not quote it as a threshold.

### Comparison to the literature

Wu et al. [[3]](#references) — the paper that introduced this erasure-conversion scheme for ¹⁷¹Yb —
report circuit-level surface-code thresholds rising from **0.937%** (no conversion) to **4.15%** at 98%
erasure conversion. Our measurements are **consistent in structure** — a several-fold gain from erasure
conversion, concentrated in herald-aware decoding — and comparable in magnitude at low `r_e` (our
baseline 1.38% vs their 0.937%, same order). We deliberately **do not** claim a matching 98% number: our
`d ≥ 7` herald fit does not converge at `r_e = 0.98`, so we have no clean threshold to compare there.
Where we can compare, the absolute values differ for honest reasons:

- **Noise model.** Our per-gate channel — `DEPOLARIZE2(p(1−r_e))` plus an independent
  `HERALDED_ERASE((2/3)·p·r_e)` on *each* qubit, sized to hold the per-gate non-identity budget at `p`
  for every `r_e` — need not match Wu et al.'s channel exactly; how the erasure and residual-Pauli
  weight are apportioned shifts the threshold.
- **Fitting choices.** We quote the asymptotic `d ≥ 7` collapse fit; our measured finite-size drift is
  **+0.33%** at `r_e = 0` (1.38% → 1.71% all-`d`) and larger at higher `r_e` (see
  [Caveats](#caveats-and-limitations)).
- **Measurement / reset / idle errors.** The uniform-sweep assumption `p_meas = p_reset = p_idle = p`
  is a modeling choice that moves the crossing.

### Where the gain comes from

With the per-gate budget held constant, the ablation splits cleanly and the split is the *opposite* of
what an "erasure conversion is intrinsically cheaper" intuition suggests:

- **Erasure conversion alone** — the `blind_mwpm` decoder, which discards the herald bits — moves the
  threshold hardly at all: **1.44% → 1.49% → 1.64%** across `r_e = 0 / 0.5 / 0.98`. Once the total
  per-gate error budget is fixed, replacing residual Pauli with (blind) erasure noise is close to a wash.
- **Herald-conditioned decoding** — reading each shot's herald bits and zeroing the conditioned edges —
  is where the benefit lives: **1.38% → 2.32%** at `r_e = 0.5` alone, and a growing sub-threshold
  suppression at higher `r_e` and distance (below). Knowing *where* the erasure occurred is the whole
  advantage.

(An earlier version of this README decomposed a `4.5×` gain into `1.6×` blind and `2.7×` herald. That
decomposition came from the pre-fix pipeline on the budget-shrinking noise model and does not survive
either correction — the blind factor in particular was largely the shrinking-budget artifact.)

---

## Exact definitions

### Noise channels (§5)

The biased-erasure model follows Wu et al. [[3]](#references); "biased" refers to the leakage favoring
one computational state, in the sense of Sahay et al. [[8]](#references). Per two-qubit gate on qubits
`(a, b)`, with physical rate `p` and erasure fraction `r_e`:

| Channel | Rate | Where |
|---|---|---|
| `DEPOLARIZE2` | `p·(1 − r_e)` | on `(a, b)` — residual Pauli component |
| `HERALDED_ERASE` | `(2/3)·p·r_e` | on `a`, and independently on `b` — appends a herald bit |
| `X_ERROR` | `p_meas` | immediately before every `M`/`MR` |
| `X_ERROR` | `p_reset` | immediately after every `R` |
| `DEPOLARIZE1` | `p_idle` | on qubits idle during a TICK |

`p_meas = p_reset = p_idle = p` for the uniform sweep. A `HERALDED_ERASE(q)` replaces the qubit with
`I/2` (each of `{I,X,Y,Z}` w.p. `q/4`) **and** records a herald bit; conditioned on that bit, each
non-trivial Pauli has probability ½. An erasure is a *non-identity* error only ¾ of the time (the `I`
outcome still heralds but causes no syndrome), so the two qubits contribute `2·q·¾` of non-identity
mass; setting that equal to the erasure share `p·r_e` gives the rate `q = (2/3)·p·r_e`. This holds the
per-gate non-identity error probability at `p` for **every** `r_e`, so `r_e` genuinely *converts* the
2q-gate error budget rather than shrinking it. (Measurement/reset/idle errors are not erasure-converted,
so the *circuit-wide* heralded fraction is below `r_e` — at `d = 5, p = 0.02, r_e = 0.98` about 55% of
the DEM error mass is heralded, not 98%.)

### Detector algebra (§3.4)

The surface-code detector/syndrome construction is standard, following Dennis et al. [[9]](#references)
and Fowler et al. [[10]](#references). Let `s_a(t)` be ancilla `a`'s outcome in round `t`, `m_q` the final
data measurement.

| Where | Detector | Declared for |
|---|---|---|
| Round 0 | `D_{a,0} = s_a(0)` | **Z-ancillas only** (data starts in `|0…0⟩`; X outcomes random) |
| Rounds 1…T−1 | `D_{a,t} = s_a(t) ⊕ s_a(t−1)` | all ancillas |
| Time close | `D_{a,T} = (⊕_{q∈∂a} m_q) ⊕ s_a(T−1)` | **Z-ancillas only** |

Syndrome detectors carry 3 coordinates `(x, y, t)`; **herald detectors carry a 4th sentinel
coordinate `(x, y, t, 1)`** — the key the DEM partition splits on.

### Per-round conversion (§10)

Sinter reports a shot-level `P_L_shot` over a `T`-round experiment; the comparable per-round rate is

```
p_L = ½·(1 − (1 − 2·P_L_shot)^(1/T))
```

with fixed points `P_L_shot = ½ → p_L = ½` and `T = 1 → p_L = P_L_shot`. The finite-size collapse is fit
in this per-round variable, which matters at high thresholds: at `d = 11` (`T = 11`) a legitimate
near-crossing point still has `P_L_shot ≈ 0.45`, so any saturation cut must act on `p_L`, not `P_L_shot`.

### Finite-size scaling ansatz (§10)

Near the crossing the collapsed data fit the quadratic finite-size-scaling ansatz of Wang, Harrington &
Preskill [[11]](#references):

```
p_L(p, d) = A + B·x + C·x²,     x = (p − p_th)·d^{1/ν}
```

for `(p_th, ν, A, B, C)` by weighted least squares, weighting each point by its **1σ Wilson standard
error** (the 95% Wilson half-width ÷ `z₀.₉₇₅`) in per-round space. Points at or above per-round
`p_L = 0.4` (saturating toward the ½ limit, outside the local ansatz) are excluded, and the fit reports
its **χ²/dof**; a fit whose `p_th` or `ν` pins at an optimiser bound is reported `converged=False` rather
than as a spurious number. `p_th ± σ` comes from a seeded parametric bootstrap over the sinter counts.

---

## The herald-aware vs. blind ablation

![Ablation](figures/ablation.png)

Same shots, two decoders: `herald_mwpm` (solid) reads the herald bits and zeroes the conditioned edges;
`blind_mwpm` (dashed) strips the herald columns and folds the erasures in as extra static depolarizing
noise. This comparison is **decoder-vs-decoder on identical shots**, so it is independent of how the
noise budget is normalised across `r_e` and needs no threshold fit — the most robust result here.
Herald-awareness wins at every distance below threshold.

**The forced-erasure story (from the M6 correctness test).** Erase the two off-row data qubits of the
`x = 1` vertical column (`(1,3)` and `(1,5)`) with probability-1 probe erasures. When both suffer an
X-type error (¼ of shots), the two X's meet at the `(0,4)` Z-check (which fires twice and cancels) and
terminate on the top boundary, leaving a syndrome of a **single** detector `D(2,2,1)`. That is
*exactly* the signature of a single X on `(1,1)` — a qubit that lies **on** the logical support and
whose correction flips `Z̄`. The two situations are indistinguishable from the syndrome alone:

- **blind** picks the cheaper one-edge explanation through `(1,1)` → flips the logical → **wrong**;
- **herald-aware** sees both erased edges at weight 0, so the two-edge no-flip route costs 0 → **right**.

On 2048 shots the blind decoder fails ~475 times; the herald-aware decoder fails 0.

### Sub-threshold suppression

How much herald-conditioning suppresses the logical rate, measured as the ratio
`p_L(blind) / p_L(herald)` at a fixed sub-threshold `p = 1.0%` (below every converged threshold above),
by distance. Computed directly from circuits with a fixed seed (`scripts/ablation_table.py`, 100 000
shots); regenerate with `uv run python scripts/ablation_table.py`.

| `d` | `r_e = 0.5` | `r_e = 0.98` |
|---|---|---|
| 3 | 1.3× | 2.1× |
| 5 | 1.9× | 8.8× |
| 7 | 2.8× | 41.2× |
| 9 | 4.1× | 159.1× |
| 11 | 6.4× | 1034× † |

*Raw herald/blind failure counts in 100 000 shots (for the low-statistics flag): at `r_e = 0.98`,
`d = 7` is 105/4169, `d = 9` is 23/3543, `d = 11` is 3/3016. † At `d = 11, r_e = 0.98` the herald
decoder produces only 3 errors (95% Wilson shot-rate CI `[1.0, 8.8]×10⁻⁵`), so 1034× is a low-statistics
lower bound, not a precise value.*

The advantage compounds with **both** distance and erasure fraction: at `r_e = 0.98` a large code makes
dramatically fewer logical errors with herald-conditioning than without, on identical shots. This is the
sub-threshold shadow of the growing herald advantage — the larger the code, the more the two exponential
suppression rates diverge. Cells marked † are low-statistics (few herald errors in 100 000 shots), so
the ratio there is a lower bound, not a precise value.

---

## Distance suppression (Λ factor)

![Lambda vs p](figures/lambda_vs_p.png)

`Λ = p_L(d) / p_L(d+2)` at fixed `p`, with bootstrap CIs. `Λ > 1` below threshold means adding two
rows of distance suppresses the logical rate; larger `r_e` gives larger `Λ` (heralds convert would-be
undetected errors into known-location ones). Points are shown only where **both** distances have ≥ 50
observed logical errors — the low-statistics tail (large `d`, low `p`, high `r_e`, where the ratio
divides a near-zero by a near-zero) is dropped, since its bootstrap CIs otherwise swamp the signal.

---

## Caveats and limitations

A few things worth being straight about.

**The `r_e = 0.98` herald threshold is not measured here.** At near-full conversion the herald crossing
sits near ~6%, where distances 9 and 11 saturate to the ½ coin-flip limit within one `p`-grid step, so
the `d ≥ 7` collapse window holds too few above-crossing large-`d` points to constrain `ν`; the fit rails
`ν` to its bound and is reported `converged=False`. Resolving it would need a denser, higher-statistics
`p`-grid straddling the crossing and possibly larger distances — future work. The ablation (above) is the
robust high-`r_e` result.

**The effective crossing drifts with distance.** Curves for small `d` carry the largest finite-size
corrections, so an all-distance collapse fit places the effective crossing above the asymptotic value —
+0.33% at `r_e = 0` (1.38% → 1.71%), and larger at higher `r_e`. Every threshold quoted here uses the
`d ≥ 7` fit (`fit_threshold(..., d_min=7)`), which drops `d ∈ {3,5}`.

**ν is only weakly constrained at `d ≤ 11`.** With three large distances and a coarse `p`-grid the
scaling exponent has wide bootstrap bars throughout; the thresholds `p_th` are far better determined than
`ν`. The fit reports `χ²/dof` so this is visible rather than hidden, and a fit that cannot constrain a
parameter pins it at a bound and is flagged non-converged rather than quoted.

**`estimate_crossing` is not reliable on saturated data.** Above threshold the curves re-converge toward
the ½ at-chance limit, producing a spurious second minimum of the cross-distance spread, and a
ragged/saturated high-`p` tail can pull the automatic estimate past the real crossing. An
ordering-inversion guard rejects candidates where `p_L` already increases with `d`, but it is only a
starting point for the fit window; every quoted threshold comes from the `d ≥ 7` collapse fit, never from
the raw crossing estimate alone.

---

## Hook regression: the schedule is load-bearing

![Hook regression](figures/hook_regression.png)

`shortest_graphlike_error()` length vs `d` for the hook-safe CX schedule (§3.2) and a deliberately
**broken** one — computed directly from the circuit builder, no Monte Carlo. Hook errors and the
CX-ordering that controls them are the surface-code "measurement gadget" analysis of
Dennis et al. [[9]](#references) and Fowler et al. [[10]](#references). The hook-safe schedule orients
every mid-window ancilla-fault hook *perpendicular* to the logical operator, so the graphlike distance
tracks `d`. The broken schedule (both bases sharing the ᴎ-order) runs the X-ancilla hooks
*parallel* to the vertical logical, halving the distance to `⌈(d+1)/2⌉`. This figure is the regression
test that proves the scheduling logic matters.

---

## DEM hand-verification

The one genuinely error-prone component is the DEM partition: splitting Stim's decomposed detector
error model into a herald-free Pauli sub-DEM plus a per-herald table of conditioned matching edges.
Rather than trust it, every mechanism was verified **by hand on paper** on the smallest non-trivial
instance (d=3, T=2, one `HERALDED_ERASE(0.25)`) and reconciled against Stim's actual output. The full
worksheet — prediction tables, reconciliation records, and the propagation diagrams — is in
[`docs/dem_worksheet.md`](docs/dem_worksheet.md).

Three probes pin down the three edge geometries a heralded erasure can produce:

- **Probe A — center data qubit `(3,3)`:** interior qubit; each Pauli component fires a *pair* of
  round-1 detectors. Confirms the herald splits into its own `^` component and the Y mechanism
  re-contributes the X and Z component edges (deduped, not a 4-detector hyperedge).
- **Probe B — corner data qubit `(1,1)`:** touched by only two checks, so each component is a
  **single-detector boundary edge**; the X component correctly carries the `L0` flip since `(1,1)` is
  on the logical support.
- **Probe C — mid-round ancilla erasure `(2,2)`:** the surprising one. The X component corrupts the
  ancilla's *own* measurement across consecutive rounds → a **time-like** edge `{D(2,2,1), D(2,2,2)}`;
  the Z component propagates to a single data neighbour. Tracing it by hand revealed a genuine
  cancellation — **two erased data qubits, but only one detector fires** — because one neighbour's
  check already completed earlier in the round and its propagated flip cancels against an identical
  flip a round later. That non-obvious result is exactly what the worksheet exists to catch.

---

## Decoder throughput

Two-tier dispatch (M6): herald-free shots are decoded in one batched call on the base matcher (the
fast path); heralded shots are grouped by their fired-herald *signature*, one reweighted matcher built
per distinct signature and decoded batched (the slow path).

| Path | Regime | Throughput (d=5, p=3e-2) |
|---|---|---|
| Fast | `r_e = 0`, no heralds | ~80,000 shots/s |
| Slow | `r_e = 0.98`, ~every shot heralded | ~1,650 shots/s |

**The 69× → 48× optimization.** The first slow-path implementation rebuilt a `pymatching.Matching`
per shot via Python `add_edge` calls: ~1,235 shots/s, a **69×** gap over the fast path — over the 50×
budget. Replacing it with herald-signature grouping (one matcher per distinct fired-herald set, built
vectorized from a precompiled sparse check matrix) brought it to ~1,650 shots/s, a **48×** gap. At this
operating point `E[fired heralds] ≈ 12` of ~800, so nearly every signature is unique and per-signature
matcher construction is the physical floor.

---

## Architecture

```
src/erasure_qec/
├── config.py              # frozen NoiseParams / ExperimentConfig (+ YAML I/O)
├── circuits/
│   ├── layout.py          # pure geometry (no stim imports)
│   ├── scheduling.py      # hook-safe + broken CX schedules (data)
│   └── builder.py         # the ONLY module that emits stim.Circuit
├── noise/
│   ├── model.py           # NoiseParams -> concrete channel rates
│   └── injector.py        # Null / PauliOnly / BiasedErasure injectors
├── decoding/
│   ├── dem_partition.py   # DEM -> (Pauli sub-DEM, herald->edges table)  [§6]
│   ├── herald_matching.py # per-shot reweighted matcher, fast/slow dispatch [§8]
│   └── sinter_adapter.py  # herald_mwpm / blind_mwpm sinter.Decoders  [§9]
└── analysis/
    ├── statistics.py      # per-round conversion, Wilson/bootstrap CIs, Λ  [§10]
    ├── threshold_fit.py    # finite-size scaling collapse fit  [§10]
    ├── synthetic.py        # ansatz-exact fixtures for the deterministic tests
    └── plotting.py         # all figures -> figures/  (deterministic)

experiments/
├── collect_threshold_sweep.py   # config -> sinter.Tasks -> data/<name>.csv
├── collect_lambda_scan.py       # fixed sub-threshold p, sweep d
├── make_synthetic_fixtures.py   # regenerate committed fixtures
└── configs/{baseline_pauli,erasure_r50,erasure_r98}.yaml

scripts/
├── audit_checks.py              # re-runnable reproduction of the external audit (docs/AUDIT.md)
└── ablation_table.py            # herald-vs-blind suppression table, direct from circuits
```

Build order was strictly milestone-by-milestone (M0–M9): the noiseless circuit and its determinism
tests came first; noise, DEM partitioning, and decoding only began once the distance/hook invariants
were green. See [PLAN.md](PLAN.md) for the full specification, and [`docs/AUDIT.md`](docs/AUDIT.md) for
the external audit that drove the noise-model and threshold-fit corrections.

---

## Reproducing

Setup (Python ≥ 3.12, [uv](https://github.com/astral-sh/uv)):

```bash
uv sync --group dev
uv run pytest -q                  # full suite (a few slow Monte-Carlo tests included)
uv run pytest -m "not slow" -q    # fast lane
```

Regenerate every figure from the committed real sweeps — byte-stable, deterministic:

```bash
uv run python -m erasure_qec.analysis.plotting --data-dir data --figures-dir figures
```

Reproduce the audit checks and the herald-vs-blind suppression table:

```bash
uv run python scripts/audit_checks.py       # docs/AUDIT.md table (pre-fix expectations vs HEAD)
uv run python scripts/ablation_table.py      # sub-threshold suppression, direct from circuits
```

Re-collect the real Monte-Carlo sweeps (resumable; re-run to accumulate). To collect under the fixed
model from scratch, remove the committed `data/<name>.csv` first:

```bash
uv run python experiments/collect_threshold_sweep.py experiments/configs/erasure_r50.yaml
uv run python experiments/collect_lambda_scan.py   experiments/configs/erasure_r50.yaml
```

**Runtime estimates** (per config, `max_shots = 1e5`, `max_errors = 1e3`, all cores). Cost is dominated
by the largest distances and by `r_e` (heralded shots take the slow path):

| Config | `r_e` | Rough wall-clock (8 cores) |
|---|---|---|
| `baseline_pauli` | 0.0 | ~10–20 min (all fast path) |
| `erasure_r50` | 0.5 | ~1–2 h |
| `erasure_r98` | 0.98 | ~4–8 h (slow path dominant at large `d`) |

Regenerate the synthetic fixtures (deterministic, no sampling):

```bash
uv run python experiments/make_synthetic_fixtures.py
```

## Related work and future directions

Directions this project does not yet build on — imperfect/delayed heralds, other hardware profiles
(trapped ions), logical algorithms, alternative decoders, and a denser high-`r_e` sweep to resolve the
`r_e = 0.98` herald threshold — with a concrete next step for each, are collected in
[`docs/FUTURE_WORK.md`](docs/FUTURE_WORK.md).

## References

1. Tomita, Y. & Svore, K. M. *Low-distance surface codes under realistic quantum noise.* Phys. Rev. A
   **90**, 062320 (2014). [arXiv:1404.3747](https://arxiv.org/abs/1404.3747)
2. Gidney, C. *Stim: a fast stabilizer circuit simulator.* Quantum **5**, 497 (2021).
   [arXiv:2103.02202](https://arxiv.org/abs/2103.02202)
3. Wu, Y., Kolkowitz, S., Puri, S. & Thompson, J. D. *Erasure conversion for fault-tolerant quantum
   computing in alkaline earth Rydberg atom arrays.* Nat. Commun. **13**, 4657 (2022).
   [arXiv:2201.03540](https://arxiv.org/abs/2201.03540)
4. Ma, S. et al. *High-fidelity gates and mid-circuit erasure conversion in an atomic qubit.* Nature
   **622**, 279 (2023).
5. Kubica, A. et al. *Erasure Qubits: Overcoming the T₁ Limit in Superconducting Circuits.* Phys. Rev. X
   **13**, 041022 (2023). [arXiv:2208.05461](https://arxiv.org/abs/2208.05461)
6. Teoh, J. D. et al. *Dual-rail encoding with superconducting cavities.* PNAS **120**, e2221736120
   (2023). [arXiv:2212.12077](https://arxiv.org/abs/2212.12077)
7. Higgott, O. & Gidney, C. *Sparse Blossom: correcting a million errors per core second with
   minimum-weight matching.* Quantum **9**, 1600 (2025).
   [arXiv:2303.15933](https://arxiv.org/abs/2303.15933)
8. Sahay, K., Jin, J., Claes, J., Thompson, J. D. & Puri, S. *High-Threshold Codes for Neutral-Atom
   Qubits with Biased Erasure Errors.* Phys. Rev. X **13**, 041013 (2023).
   [arXiv:2302.03063](https://arxiv.org/abs/2302.03063)
9. Dennis, E., Kitaev, A., Landahl, A. & Preskill, J. *Topological quantum memory.* J. Math. Phys.
   **43**, 4452 (2002). [arXiv:quant-ph/0110143](https://arxiv.org/abs/quant-ph/0110143)
10. Fowler, A. G., Mariantoni, M., Martinis, J. M. & Cleland, A. N. *Surface codes: Towards practical
    large-scale quantum computation.* Phys. Rev. A **86**, 032324 (2012).
11. Wang, C., Harrington, J. & Preskill, J. *Confinement-Higgs transition in a disordered gauge theory
    and the accuracy threshold for quantum memory.* Ann. Phys. **303**, 31 (2003).
    [arXiv:quant-ph/0207088](https://arxiv.org/abs/quant-ph/0207088)
12. Barrett, S. D. & Stace, T. M. *Fault Tolerant Quantum Computation with Very High Threshold for Loss
    Errors.* Phys. Rev. Lett. **105**, 200502 (2010).
    [arXiv:1005.2456](https://arxiv.org/abs/1005.2456)

## License

MIT — see [LICENSE](LICENSE).
