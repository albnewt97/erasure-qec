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

**Headline result (measured).** The circuit-level threshold rises from **1.4%** (pure Pauli) to **6.3%**
at `r_e = 0.98` — a **≈4.5× improvement** — of which erasure conversion *alone* (the blind decoder,
which ignores the herald bits) contributes **1.6×**; the remaining **2.7×** comes from the
herald-conditioned decoder.

![Threshold panels](figures/threshold_panels.png)

*Per-round `p_L` vs physical error rate `p`, one curve per distance `d ∈ {3,5,7,9,11}`, one panel per
erasure fraction; the dashed line is the fitted threshold (`d ≥ 7` collapse fit) with its bootstrap CI.*

> **Data status.** These figures are generated from the **real Monte-Carlo sweeps** in `data/`
> (gitignored; regenerate with the commands under [Reproducing](#reproducing)). All quoted thresholds
> are the asymptotic `d ≥ 7` finite-size-scaling collapse fits. The committed **synthetic fixtures**
> (`tests/fixtures/*.csv`) remain — they drive the deterministic, byte-stable analysis/plotting tests
> and stand in when the raw sweeps are absent.

---

## Threshold scaling

The panels above are per-round logical error rate `p_L` vs physical error rate `p`. Each panel's dashed
line (with bootstrap CI band) is the fitted threshold `p_th` from the **asymptotic `d ≥ 7` collapse fit**
(see [Caveats](#caveats-and-limitations)); the inset is the collapse onto the single scaling variable
`x = (p − p_th)·d^{1/ν}`.

### Measured thresholds (`d ≥ 7` collapse fit)

| `r_e` | `herald_mwpm` `p_th` | `blind_mwpm` `p_th` | notes |
|---|---|---|---|
| 0.0 | **1.38%** | **1.44%** | control: no heralds, so the decoders agree (as required) |
| 0.5 | **2.30% ± 0.15%** | **1.73% ± 0.10%** | herald and blind first diverge measurably |
| 0.98 | **6.27% ± 0.54%** (ν = 1.56) | **2.29% ± 0.18%** | near-full conversion (¹⁷¹Yb target); only sweep whose dense `p`-band pins ν |

### Comparison to the literature

Wu et al. [[3]](#references) — the paper that introduced this erasure-conversion scheme for ¹⁷¹Yb —
report circuit-level surface-code thresholds rising from **0.937%** (no conversion) to **4.15%** at 98%
erasure conversion. This work measures **1.38% → 6.27%** over the same range. The two are consistent in
structure (a several-fold threshold gain from near-total erasure conversion) and comparable in
magnitude, with our absolute values somewhat higher. This is stated as a **consistency check that
reproduces a known effect, not an improvement claim** — the offset is expected, and has honest
explanations:

- **Noise model.** Our per-gate channel — `DEPOLARIZE2(p(1−r_e))` plus an independent
  `HERALDED_ERASE(p·r_e/2)` on *each* of the two qubits — need not match Wu et al.'s channel exactly;
  small differences in how the erasure and residual-Pauli weight are apportioned shift the threshold.
- **Fitting choices.** We quote the asymptotic `d ≥ 7` collapse fit, and our own measured finite-size
  drift is **+0.3–0.4%** when the small distances are included (see
  [Caveats](#caveats-and-limitations)) — that alone accounts for a meaningful part of the gap.
- **Measurement / reset / idle errors.** The uniform-sweep assumption `p_meas = p_reset = p_idle = p`
  is a modeling choice; a different apportionment there moves the crossing.

### Where the gain comes from

The threshold climbs **1.4% → 6.3%** (≈4.5×) as `r_e → 0.98`, and the ablation splits that cleanly:

- **Erasure conversion alone** — the `blind_mwpm` decoder, which discards the herald bits and folds the
  erasures into extra static depolarizing noise — moves the threshold only **1.4% → 2.3%** (≈1.6×). This
  is the passive benefit of biased erasure: an erased qubit is depolarized, and depolarizing noise is
  cheaper to correct than the full Pauli budget it replaced.
- **Herald-conditioned decoding** — reading each shot's herald bits and zeroing the conditioned edges —
  lifts the `r_e = 0.98` threshold the rest of the way, **2.3% → 6.3%** (≈2.7×). Knowing *where* the
  erasure occurred is worth more than the conversion itself.

The full factor is `1.6× · 2.7× ≈ 4.5×`.

---

## Exact definitions

### Noise channels (§5)

The biased-erasure model follows Wu et al. [[3]](#references); "biased" refers to the leakage favoring
one computational state, in the sense of Sahay et al. [[8]](#references). Per two-qubit gate on qubits
`(a, b)`, with physical rate `p` and erasure fraction `r_e`:

| Channel | Rate | Where |
|---|---|---|
| `DEPOLARIZE2` | `p·(1 − r_e)` | on `(a, b)` — residual Pauli component |
| `HERALDED_ERASE` | `p·r_e/2` | on `a`, and independently on `b` — appends a herald bit |
| `X_ERROR` | `p_meas` | immediately before every `M`/`MR` |
| `X_ERROR` | `p_reset` | immediately after every `R` |
| `DEPOLARIZE1` | `p_idle` | on qubits idle during a TICK |

`p_meas = p_reset = p_idle = p` for the uniform sweep. A `HERALDED_ERASE(q)` replaces the qubit with
`I/2` (each of `{I,X,Y,Z}` w.p. `q/4`) **and** records a herald bit; conditioned on that bit, each
non-trivial Pauli has probability ½.

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

with fixed points `P_L_shot = ½ → p_L = ½` and `T = 1 → p_L = P_L_shot`.

### Finite-size scaling ansatz (§10)

Near the crossing the collapsed data fit the quadratic finite-size-scaling ansatz of Wang, Harrington &
Preskill [[11]](#references):

```
p_L(p, d) = A + B·x + C·x²,     x = (p − p_th)·d^{1/ν}
```

for `(p_th, ν, A, B, C)` by weighted least squares (Wilson half-widths as weights), fitting only points
within a configurable multiplicative window of a data-driven crossing estimate. `p_th ± σ` comes from a
seeded parametric bootstrap over the sinter counts.

---

## The herald-aware vs. blind ablation

![Ablation](figures/ablation.png)

Same shots, two decoders: `herald_mwpm` (solid) reads the herald bits and zeroes the conditioned edges;
`blind_mwpm` (dashed) strips the herald columns and folds the erasures in as extra static depolarizing
noise. Herald-awareness wins at every distance below threshold.

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
`p_L(blind) / p_L(herald)` at a fixed sub-threshold `p = 1.0%` (below every threshold above), by
distance:

| `d` | `r_e = 0.5` | `r_e = 0.98` |
|---|---|---|
| 3 | 1.2× | 1.8× |
| 5 | 1.8× | 6.4× |
| 7 | 2.5× | 22.6× |
| 9 | 3.2× | 144× |
| 11 | 4.5× | 469× † |

The advantage compounds with **both** distance and erasure fraction: at `r_e = 0.98` a distance-11 code
makes ~470× fewer logical errors with herald-conditioning than without, on identical shots. This is the
sub-threshold shadow of the threshold gap — the larger the code, the more the two exponential
suppression rates diverge. († At `d = 11, r_e = 0.98` the herald decoder produces only 2 errors in 10⁵
shots, so this ratio is a low-statistics lower bound, not a precise value.)

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

The effective crossing drifts with distance. Curves for small `d` carry the largest finite-size
corrections, so an all-distance collapse fit places the effective crossing above the asymptotic value —
+0.3–0.4% at `r_e = 0` (1.38% → 1.71%), and comparable or larger at higher `r_e`. Every threshold quoted
here uses the `d ≥ 7` fit (`fit_threshold(..., d_min=7)`), which drops `d ∈ {3,5}`; the panel
annotations do the same whenever ≥3 large distances are available.

ν is only weakly constrained at `d ≤ 11`. With three large distances and a coarse `p`-grid the scaling
exponent has wide bootstrap bars — `ν = 1.45 ± 1.2` (`r_e = 0`), `2.20 ± 0.9` (`r_e = 0.5`). The
exception is `r_e = 0.98`, where the dense `p`-band through the crossing pins `ν = 1.56 ± 0.48`; that
told me the fix for the earlier sweeps was denser sampling near the crossing, not more shots. The
thresholds `p_th` are far better determined than `ν` throughout.

`estimate_crossing` is not reliable on saturated data. Above threshold the curves re-converge toward the
½ at-chance limit, producing a spurious second minimum of the cross-distance spread, and a
ragged/saturated high-`p` tail can pull the automatic estimate past the real crossing — which is what
happened on the `r_e = 0.5` sweep. An ordering-inversion guard rejects candidates where `p_L` already
increases with `d`, but it is only a starting point for the fit window; every quoted threshold comes
from the `d ≥ 7` collapse fit, never from the raw crossing estimate alone.

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
```

Build order was strictly milestone-by-milestone (M0–M9): the noiseless circuit and its determinism
tests came first; noise, DEM partitioning, and decoding only began once the distance/hook invariants
were green. See [PLAN.md](PLAN.md) for the full specification.

---

## Reproducing

Setup (Python ≥ 3.12, [uv](https://github.com/astral-sh/uv)):

```bash
uv sync --group dev
uv run pytest -q                  # full suite (a few slow Monte-Carlo tests included)
uv run pytest -m "not slow" -q    # fast lane
```

Regenerate every figure — byte-stable, deterministic:

```bash
# from the real sweeps in data/ (what the figures above show):
uv run python -m erasure_qec.analysis.plotting --data-dir data --figures-dir figures
# from the committed synthetic fixtures (used by the deterministic plotting tests):
uv run python -m erasure_qec.analysis.plotting --data-dir tests/fixtures --figures-dir figures
```

Run the real Monte-Carlo sweeps (resumable; re-run to accumulate):

```bash
uv run python experiments/collect_threshold_sweep.py experiments/configs/erasure_r50.yaml
uv run python experiments/collect_lambda_scan.py   experiments/configs/erasure_r50.yaml
```

**Runtime estimates** (per config, `max_shots = 1e5`, `max_errors = 1e3`, all cores). The threshold
sweep is 5 distances × (13 log-spaced + a dense linear band through the crossing) `p` × 2 decoders
(~180–230 tasks). Cost is dominated by the largest distances and by `r_e` (heralded shots take the slow
path):

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
(trapped ions), logical algorithms, and alternative decoders — with a concrete next step for each, are
collected in [`docs/FUTURE_WORK.md`](docs/FUTURE_WORK.md).

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
