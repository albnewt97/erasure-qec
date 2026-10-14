# Caveat ledger

Checklist of every qualifying sentence in `README.md` at `da47e3a` (K rows), plus caveats stated only
in `docs/AUDIT.md` that qualify a number `README.md` publishes (A rows). A caveat may be moved or
reworded in later phases only if this ledger records its destination; it may not be dropped while
the number it qualifies is still published.

Conventions:

- **Line**: for K rows, the README line range of the sentence. For A rows, the README line(s) of the
  number the caveat qualifies, then the AUDIT line range the text is quoted from.
- **Text** is verbatim, at most 300 characters, with line breaks, runs of spaces and blockquote
  markers (`> `) collapsed to a single space. A `|` inside the text is escaped as `\|`. A sentence
  longer than 300 characters is cut at a clause boundary and continued in the next row; text is never
  paraphrased.
- **Category**: model | statistics | comparison | reproducibility | performance | scope.
- **Status** is `open`; **Destination** is blank until a later phase moves the caveat.

## README.md

| ID | Line | Text (verbatim) | Category | Status | Destination |
|---|---|---|---|---|---|
| K001 | 5-6 | with circuit-level noise [[1]](#references) (Z basis, in [Stim](https://github.com/quantumlib/Stim) [[2]](#references)) | scope | open | |
| K002 | 7-9 | a fraction `r_e` of each two-qubit gate's error budget becomes a *heralded erasure* — the hardware flags *which* qubit was disturbed and *when*. | model | open | |
| K003 | 9-10 | This models **¹⁷¹Yb neutral-atom Rydberg gates** [[3,4]](#references) and **dual-rail superconducting cavities** [[5,6]](#references). | model | open | |
| K004 | 12-13 | a heralded site has conditional error probability ½, so `ln((1−p)/p) = 0` and a correction routes through it for free. | model | open | |
| K005 | 15-16 | Everything here is a **Z-basis memory** threshold under a single observable — not a full X+Z circuit-level threshold (see [Modelling assumptions](#modelling-assumptions)). | scope | open | |
| K006 | 16-18 | The detailed statistical methodology, audit history, and full diagnostics live in [`docs/AUDIT.md`](docs/AUDIT.md); this README states the results and their caveats in short form. | scope | open | |
| K007 | 20-21 | **Headline (measured).** Holding the *per-gate error budget constant* across `r_e`, the gain from erasure is almost entirely from **reading the herald bits**, not from the conversion itself: | model | open | |
| K008 | 24-27 | That `2.2` upper endpoint is **not an interval estimate**: it is set by bound-pinned bootstrap replicates and lands on a discrete atom, `seed = 0` is one of only 3 of 30 seeds giving it, and 25 of 30 give `1.8` — see the [†footnote](#measured-thresholds-d--7-collapse-fit) on the results table. | statistics | open | |
| K009 | 27-28 | The point estimates (`1.44%`, `1.49%`) carry no RNG and are exact. | reproducibility | open | |
| K010 | 33-34 | The larger `d = 11` bound rests on 3 herald errors and is weak evidence — see the [table](#sub-threshold-suppression). | statistics | open | |
| K011 | 36-37 | At **`r_e = 0.98`** the collapse fit resolves **neither** threshold and the code says so (`resolved=False`): the herald fit rails ν to its bound (not converged), | statistics | open | |
| K012 | 37-39 | and the blind fit converges but its 95% CI runs from ~1.5% to an upper endpoint that lands anywhere in **5.5–8.5%** depending on the bootstrap seed (6.3% at `seed = 0`) — 2.4–4.3× wider than the threshold itself, so it is not a resolved value. | statistics | open | |
| K013 | 39-42 | The upper endpoint is not a measurement: **three of the ten seeds put it exactly on a sweep grid point** (`p = 0.055, 0.076, 0.085`), i.e. it is set by bound-pinned replicates rather than by a fitted crossing | statistics | open | |
| K014 | 42-43 | and *both ends of the 5.5–8.5% span are themselves grid pins*, so the width of that span is a property of the sweep grid, not a spread of fitted values. | statistics | open | |
| K015 | 43-44 | The near-full-conversion story therefore leads with the deterministic ablation. | scope | open | |
| K016 | 45-46 | Reproduce every number in this paragraph with `uv run python scripts/paired_separation.py --full-stability` (~4.7 min). | reproducibility | open | |
| K017 | 50-51 | dashed line = fitted `d ≥ 7` threshold with its bootstrap CI. | statistics | open | |
| K018 | 51 | Unresolved fits carry no line ("no fit"). | scope | open | |
| K019 | 53-54 | **Data status.** Figures are generated from the **committed real Monte-Carlo sweeps** in `data/` (`baseline_pauli.csv`, `erasure_r50.csv`, `erasure_r98.csv`). | reproducibility | open | |
| K020 | 54-57 | Point estimates, χ²/dof and every converged/resolved verdict reproduce **exactly** across NumPy/SciPy versions; bootstrap CI endpoints do **not** — they drift with the NumPy RNG, and `sinter.collect` takes no seed so the sweeps themselves cannot be regenerated shot-for-shot. | reproducibility | open | |
| K021 | 57-58 | See [`docs/AUDIT.md`](docs/AUDIT.md#reproducibility-and-environment-sensitivity-2026-09-06). | reproducibility | open | |
| K022 | 59 | Pinned regression snapshots live at `tests/fixtures/real_*.csv`. | reproducibility | open | |
| K023 | 59-62 | Committed **synthetic fixtures** (`tests/fixtures/synthetic_*.csv`) drive the byte-stable plotting/analysis tests and are **not measurements** — they are generated from the same collapse variable the fitter inverts (one is deliberately adversarial and the fitter must *reject* it). | scope | open | |
| K024 | 62-63 | Stale pre-fix sweeps are preserved untracked under `data/stale_old_model/`. | reproducibility | open | |
| K025 | 69-71 | 95% bootstrap percentile CI, `n_boot = 1000`, `seed = 0`, endpoints to **one decimal** — the second decimal moves with the seed and with the NumPy build (see [`docs/AUDIT.md`](docs/AUDIT.md#reproducibility-and-environment-sensitivity-2026-09-06)). | reproducibility | open | |
| K026 | 71-72 | `p_th`, ν and χ²/dof carry no RNG and are exact. | reproducibility | open | |
| K027 | 72-73 | A fit is only quoted as a threshold when the code marks it `resolved` (converged **and** CI narrower than `p_th`). | statistics | open | |
| K028 | 77 | control: no heralds → decoders agree within CIs | statistics | open | |
| K029 | 78 | herald above blind — separation **significant**, see Δ | statistics | open | |
| K030 | 79 | **not resolved** (CI `[1.5, 5.5–8.5]`, upper endpoint seed-dependent and grid-pinned — see above) | statistics | open | |
| K031 | 79 | fit resolves neither; use the ablation | scope | open | |
| K032 | 81-82 | † **The `r_e = 0.5` blind upper endpoint is not a fitted crossing, and `2.2` is an unlucky seed.** | statistics | open | |
| K033 | 82-85 | It is a *bound pin*: on ~7.5% of bootstrap replicates the resampled blind crossing estimate runs into the saturated tail, the fit window's **lower** edge lands on a high grid point, and `curve_fit` pins `p_th` there. | statistics | open | |
| K034 | 85 | Those pins form discrete atoms, and the 97.5th percentile falls inside one. | statistics | open | |
| K035 | 85-87 | Over 30 seeds the endpoint takes **four distinct values**, and **25 of 30 seeds give 1.80%** — `seed = 0`, which this table quotes, is one of only **3 of 30** that give 2.154%. | statistics | open | |
| K036 | 87 | Excluding edge-pinned replicates moves it from 2.154% to 1.601%. | statistics | open | |
| K037 | 88-89 | So the modal interval at this precision is `[1.4, 1.8]`, not `[1.4, 2.2]`; the published cell is the **wider, conservative** end of a quantised endpoint, not a threshold the fit supports at 2.2%. | statistics | open | |
| K038 | 90-91 | The `p_th` point estimate (1.49%), ν and χ²/dof are unaffected — they carry no RNG — and no verdict changes (`converged` and `resolved` are True for all 30 seeds). | statistics | open | |
| K039 | 92-93 | Reproduce with `uv run python scripts/paired_separation.py --full-stability`; the 30-seed enumeration is in [`docs/AUDIT.md`](docs/AUDIT.md#the-marginal-ci-upper-endpoints-are-bound-pins-2026-09-06). | reproducibility | open | |
| K040 | 95-96 | **Is the `r_e = 0.5` separation significant? Yes — and comparing the marginal CIs above is not how you test it** (`var(Δ)` is not the sum of the marginal variances). | statistics | open | |
| K041 | 96-98 | The valid test is a bootstrap of the difference `Δ = p_th(blind) − p_th(herald)`: **Δ = −0.83%, 95% CI `[−1.2, −0.7]`, excludes zero**, stable across seeds (all of seeds 0–9 exclude zero — | statistics | open | |
| K042 | 99 | (control `r_e = 0` gives Δ consistent with zero) | statistics | open | |
| K043 | 99-102 | Two caveats, both worked out in full in [`docs/AUDIT.md`](docs/AUDIT.md#paired-decoder-separation-phase-4): (i) the sweeps are sampled *independently* per decoder, so this is an *unpaired* difference bootstrap — which can only be wider than a paired one, making the result **conservative**; | statistics | open | |
| K044 | 102-105 | (ii) the CI conditions on both decoders converging, and the ~10% discards are **not missing-at-random** (they cluster at high blind `p_th`), but a tipping-point bound shows the claim survives — the upper endpoint reaches zero only if ≥ 24 of the 102 discards had `Δ ≥ 0` | statistics | open | |
| K045 | 105-106 | (seed 0; 24–26 of 77–115 across seeds 0–9, from `uv run python scripts/paired_separation.py --full-stability`), while their recorded `p_th` imply `Δ ≥ 0` for only 2. | statistics | open | |
| K046 | 106-107 | So: significant *conditional on convergence*. | statistics | open | |
| K047 | 107-108 | The [ablation](#the-herald-aware-vs-blind-ablation), which needs no fit, is the primary evidence. | scope | open | |
| K048 | 112-113 | **Our numbers are not on a comparable footing with theirs.** | comparison | open | |
| K049 | 113-115 | The comparison below is *structural only* — both show the threshold rising with the converted fraction — and no numerical agreement is claimed at any `r_e`, including `r_e = 0`. | comparison | open | |
| K050 | 115-116 | We also do not claim a matching 98% number: our `d ≥ 7` herald fit does not resolve at `r_e = 0.98`. | comparison | open | |
| K051 | 118 | The one axis that *does* line up is the conversion fraction. | comparison | open | |
| K052 | 118-120 | Their `R_e` is defined per two-qubit gate exactly as our `r_e` is: "each two-qubit gate experiences either a Pauli error with probability `p_p = p(1 − R_e)`, or an erasure with probability `p_e = pR_e`" (§III). | comparison | open | |
| K053 | 120 | Everything else differs: | comparison | open | |
| K054 | 122-124 | **Code.** They simulate the **planar XZZX** surface code (§III, after Bonilla Ataides et al. [[12]](#references)); we simulate the rotated surface code. | comparison | open | |
| K055 | 124-126 | They state it "performs identically to the standard surface code for the case of unbiased noise", and our `HERALDED_ERASE` is unbiased — so this is a caveat rather than a fatal mismatch, and their patch is a `d × d` data-qubit array, the same count as our rotated patch. | comparison | open | |
| K056 | 127-128 | **Decoder.** They use a **weighted Union-Find** decoder and note that "the UF decoder is optimal for pure erasure errors" (Delfosse & Zémor, PRR **2**, 033042). | comparison | open | |
| K057 | 128-129 | We use **MWPM with the heralded edges' weight set to 0**, which is *not* ML-optimal for erasure. | model | open | |
| K058 | 129-132 | This gap **grows with `r_e`**: it is smallest at our baseline and largest exactly where our high-`r_e` story lives, so our `r_e = 0.98` ablation should be read as a lower bound on what an ML erasure decoder would show. | model | open | |
| K059 | 133-134 | **Observable.** Ours is a **Z-basis memory scored on a single logical observable** (one `OBSERVABLE_INCLUDE` over `logical_z_support(d)`; only Z-type detectors close at the time boundaries). | scope | open | |
| K060 | 135-136 | Theirs is "the logical failure rate" from decoding the full syndrome — they never restrict to one observable or one basis. | comparison | open | |
| K061 | 136-137 | A single-observable Z-memory cannot fail in ways a full-code failure metric counts, so **ours is the more permissive metric**. | comparison | open | |
| K062 | 138-139 | **What else is noisy.** Their Figs. 3–4 set `p_m = 0` and "neglect idle errors", so two-qubit gates are their only error locations and `R_e` is *also* their circuit-wide erasure fraction. | comparison | open | |
| K063 | 139-141 | We add unheralded measurement, reset and idle noise at rate `p`, none of which is ever converted — so our `r_e` is **not** our circuit-wide fraction: at `r_e = 0.98` only ~45% of our DEM error mass is heralded (below). | model | open | |
| K064 | 141-142 | The two `0.98` columns therefore describe different circuits, and we read the shortfall in neither direction. | comparison | open | |
| K065 | 143-144 | **Finite-size window.** Their thresholds come "from the crossing of `d = 11` and `d = 15`"; ours from a `d ≥ 7` collapse fit over `d ∈ {7, 9, 11}`. | comparison | open | |
| K066 | 146-147 | **An unexplained sign problem.** At `r_e = 0` we carry **strictly more** noise than Wu et al. | comparison | open | |
| K067 | 149-150 | More noise should push the threshold **down**, yet we report **1.38% against their 0.937%**. | comparison | open | |
| K068 | 150 | We do not know why, and we do not read this gap as agreement. | comparison | open | |
| K069 | 152-153 | The leading candidate is the **observable** difference above: a single-observable Z-memory is a strictly more permissive failure criterion than a full-code logical failure rate. | comparison | open | |
| K070 | 153-155 | **We have not demonstrated this** — the repo has no X-memory or two-observable circuit to test it against — so it stays an **open item, not a conclusion**. | comparison | open | |
| K071 | 155-156 | The **finite-size window** is a second contributor, and it is measurable: refitting the committed baseline sweep at rising `d_min` moves the herald threshold monotonically *down* | comparison | open | |
| K072 | 158 | (`d ≥ 9`, unconverged, two distances) | statistics | open | |
| K073 | 158-159 | and it has not flattened at the top of our range, so a `d = 11`/`d = 15` crossing would sit lower still. | comparison | open | |
| K074 | 159 | That plausibly accounts for part of the 0.44-point gap, not obviously all of it. | comparison | open | |
| K075 | 159-161 | Neither paper has the noise-matched run that would settle it: their SI reports `p_m = p` only at `R_e = 0.98` (threshold falls 4.15% → 2.85(1)%), with no `R_e = 0` counterpart. | comparison | open | |
| K076 | 166-167 | With the budget held constant, the split is the *opposite* of an "erasure is intrinsically cheaper" intuition. | model | open | |
| K077 | 167-168 | **Erasure conversion alone** (blind) moves the threshold `1.44% → 1.49%` from `r_e = 0` to `0.5` — close to a wash. | scope | open | |
| K078 | 169-170 | We make **no causal claim** about *why* blind barely helps (both channels are Pauli; "depolarizing is cheaper to correct" would be hand-waving) | scope | open | |
| K079 | 170-172 | the old `1.6×`-blind / `2.7×`-herald decomposition was an artifact of the pre-fix budget-shrinking model and does not survive (re-derivation in [`docs/AUDIT.md`](docs/AUDIT.md)). | model | open | |
| K080 | 178 | The erasure model follows Wu et al. [[3]](#references). | model | open | |
| K081 | 178-179 | `HERALDED_ERASE` is **unbiased** — `I/2`, all four Paulis equal (the advantage is from decoding, not Pauli bias). | model | open | |
| K082 | 185 | \| `X_ERROR` / `X_ERROR` / `DEPOLARIZE1` \| `p_meas` / `p_reset` / `p_idle` (`= p`) \| before `M`/`MR`, after `R`, on idle qubits \| | model | open | |
| K083 | 187-189 | (an erasure is non-identity only ¾ of the time, so `2·q·¾ = p·r_e`) holds the per-gate non-identity budget at `p` for **every** `r_e` (tested by `test_error_budget_invariance`), so `r_e` *converts* the budget rather than shrinking it and `p` is itself the iso-noise axis. | model | open | |
| K084 | 189-191 | Detector algebra is standard; the one non-standard convention here is that **herald detectors carry a 4th sentinel coordinate `(x, y, t, 1)`** the DEM partition splits on ([PLAN.md](PLAN.md) §3.4). | model | open | |
| K085 | 195 | **`r_e` is the 2q-gate fraction, not the circuit-wide heralded fraction.** | scope | open | |
| K086 | 195-197 | Meas/reset are never converted (idle only if `convert_idle`), so the fraction of *total DEM error mass* heralded is well below `r_e`. | model | open | |
| K087 | 197-198 | Wu et al.'s `R_e` is on the *same* per-2q-gate axis as `r_e`, but because they set `p_m = p_idle = 0` their `R_e` is simultaneously their circuit-wide fraction, whereas ours is not. | comparison | open | |
| K088 | 199-201 | Measured by `scripts/heralded_fraction.py` at `d = 5, p = 0.02`: gate-only **0.23 / 0.45** at `r_e = 0.5 / 0.98` (0.32 / 0.63 with `convert_idle`). | scope | open | |
| K089 | 201-202 | The committed sweeps use gate-only, so `r_e = 0.98` heralds **45%** of the mass, not 98%. | scope | open | |
| K090 | 202-204 | These exclude the no-op herald branches: `HERALDED_ERASE` heralds on its identity branch too, so the DEM carries mechanisms that flip a herald detector and nothing else, and a branch that causes no syndrome and no logical error is not part of the error budget. | model | open | |
| K091 | 204-206 | Counting them (the pre-2026-09-05 definition, still available as `count_no_op_heralds=True`) gave 0.30 / 0.55 and 0.42 / 0.72 — i.e. it overstated conversion. | statistics | open | |
| K092 | 214-215 | `blind_mwpm` strips the herald columns and folds erasures into static depolarizing noise. | model | open | |
| K093 | 215-217 | This is decoder-vs-decoder **on identical shots** — independent of noise normalisation, needs no fit — the most robust result here. | scope | open | |
| K094 | 217-219 | **The ablation is paired by construction** (`scripts/ablation_table.py` samples once and decodes the same `dets` with both decoders); the sinter *threshold sweeps* are not (sampled per decoder), which is why the Δ test above uses an unpaired bootstrap while the ablation needs none. | statistics | open | |
| K095 | 221-223 | Shots are sampled with probability-1 probe erasures on `(1,3)`/`(1,5)` (`probe_q = 1.0`), but **both decoders are compiled from a circuit with `probe_q = 0.02`** — a deliberate mis-specification, and load-bearing: | model | open | |
| K096 | 223-225 | it is what makes blind price those two locations at a small unconditional marginal, as a decoder that does not know the erasures occurred would. | model | open | |
| K097 | 225-226 | (It costs the herald-aware decoder nothing, since the herald bits zero those edges per shot regardless of the compiled weight.) | model | open | |
| K098 | 230-232 | No exact count is quoted: stim's seeded sampling is reproducible only on the same machine and stim version, so the draw moves between environments. | reproducibility | open | |
| K099 | 237-239 | One statistical standard across the repo: cells below the **≥ 50 observed-error gate** the Λ figure uses are reported as a lower confidence bound `> N×`, not a point ratio. | statistics | open | |
| K100 | 249-250 | † Below the 50-error gate (herald errors: `d = 9` → 23, `d = 11` → 3 of 100 000 shots at `r_e = 0.98`), so what is quoted is a **lower confidence bound**, not a point value. | statistics | open | |
| K101 | 250-251 | It is blind's Wilson **lower** limit divided by herald's Wilson **upper** limit, both mapped through the per-round conversion. | statistics | open | |
| K102 | 251-254 | Each is an endpoint of a two-sided 95% Wilson interval — a one-sided 97.5% limit — so Bonferroni makes the quoted value a **≥ 95% one-sided lower bound on that cell's ratio**. | statistics | open | |
| K103 | 254-255 | It is per-cell, **not** simultaneous across the table, and not paired: it ignores the positive correlation from decoding identical shots, which only makes it looser. | statistics | open | |
| K104 | 255-256 | (Before 2026-09-05 this column divided blind's *point estimate* by herald's Wilson upper limit, which bounded nothing; the corrected values are ~3% lower.) | statistics | open | |
| K105 | 258-259 | ‡ **The `d = 11` cell rests on 3 herald errors in 100 000 shots — treat `> 339×` as weak evidence, not a measurement.** | statistics | open | |
| K106 | 259-261 | At single-digit counts Wilson's "95%" is a nominal approximation with no guaranteed actual coverage, and the herald shot-rate CI there spans `[1.0e-5, 8.8e-5]`, nearly an order of magnitude. | statistics | open | |
| K107 | 261-262 | The honest reading is "herald errors were too rare to measure at `d = 11`", not "suppression is 339×". | statistics | open | |
| K108 | 262 | The point ratio at that cell is 1034×, which is precisely why no point ratio is quoted. | statistics | open | |
| K109 | 262-263 | The `d = 9` cell (23 errors) is also below the gate but far less fragile. | statistics | open | |
| K110 | 265-266 | The advantage clearly compounds with both distance and erasure fraction over the range where the counts support it — every `r_e = 0.5` cell and the `r_e = 0.98` cells up to `d = 7` clear the 50-error gate. | scope | open | |
| K111 | 266-268 | How far that compounding continues at `d = 9` and `d = 11` is bounded from below by this run but not measured by it; settling it needs more shots, not more distance. | scope | open | |
| K112 | 277-278 | The `Λ` *notation* is not from that paper; the anchor is for the suppression statement only. | scope | open | |
| K113 | 279-280 | Points are shown only where **both** distances have ≥ 50 observed errors; the low-statistics tail is dropped. | statistics | open | |
| K114 | 282 | The suppression table uses the **same constant, not the same rule**. | statistics | open | |
| K115 | 282-284 | `figure_lambda` requires ≥ 50 errors at *both* `d` and `d+2`; `ablation_table._fmt_ratio` gates on the **herald** count alone, because herald is the fragile denominator there. | statistics | open | |
| K116 | 284-285 | Each applies the ≥ 50 threshold to whichever count can go small in that particular statistic. | statistics | open | |
| K117 | 285-288 | On the committed sweeps the two coincide — blind counts run 3016–6290, far above 50, so gating on herald alone excludes nothing the joint rule would have kept — but they are not the same predicate and would diverge on data where blind is the sparse side. | statistics | open | |
| K118 | 288-289 | (Blind-count range from the raw-count table recorded in [`docs/AUDIT.md`](docs/AUDIT.md) on 2026-09-05, not re-measured here — the sweep is ~46 min.) | reproducibility | open | |
| K119 | 295 | This is an idealized erasure model. | model | open | |
| K120 | 295-297 | Each assumption below is **systematically optimistic** — a real device is worse in each respect — so the thresholds and suppression ratios here are upper bounds on what this scheme delivers in hardware, not predictions of it. | model | open | |
| K121 | 299 | **Erasure is injected *after* the CX, independently per qubit.** | model | open | |
| K122 | 299-301 | There is no mid-gate erasure that propagates through the remainder of the gate, and no correlated two-qubit erasure — the two qubits of a gate erase independently. | model | open | |
| K123 | 301 | Real leakage happens *during* the gate and can corrupt both qubits together. | model | open | |
| K124 | 302 | **The erased qubit keeps participating perfectly in all later rounds.** | model | open | |
| K125 | 302-303 | `HERALDED_ERASE` replaces the qubit with `I/2` for that instant only; it is a fully-functional qubit again on the next tick. | model | open | |
| K126 | 303-305 | A real ¹⁷¹Yb erasure is an atom *leaked out of the qubit subspace* — it stays broken until physically re-prepared, so a single erasure should degrade many subsequent rounds, which we do not model. | model | open | |
| K127 | 306 | **Heralds are perfect.** | model | open | |
| K128 | 306 | Instantaneous, 100% reliable, zero false positives, zero latency. | model | open | |
| K129 | 306-308 | Real fluorescence heralding has finite fidelity, false positives/negatives, and a latency the decoder would have to tolerate. | model | open | |
| K130 | 309 | **The erasure is unbiased.** | model | open | |
| K131 | 309 | `HERALDED_ERASE` is `I/2` (all four Paulis equal). | model | open | |
| K132 | 309-311 | The Z-biased leakage of Sahay et al. [[8]](#references) would give a *further* threshold advantage that we deliberately do not model; our advantage is from herald-conditioned decoding alone. | model | open | |
| K133 | 312 | **Z-memory only, single logical observable.** | scope | open | |
| K134 | 312-314 | Every "threshold" here is the **Z-basis memory** threshold under a single observable — not a full X+Z circuit-level threshold, and not a logical-gate or computation threshold. | scope | open | |
| K135 | 314 | Read every number in this repo with the "Z-memory" qualifier attached. | scope | open | |
| K136 | 320 | **Neither `r_e = 0.98` threshold is measured.** | statistics | open | |
| K137 | 320-322 | The large codes saturate to ½ within one grid step above the ~6% crossing, so the herald fit rails ν to its bound and the blind fit's CI is wider than its own `p_th`; both are `resolved=False` in code. | statistics | open | |
| K138 | 322-323 | Resolving them needs a denser, higher-statistics `p`-grid (future work). | scope | open | |
| K139 | 323 | The ablation is the robust high-`r_e` result. | scope | open | |
| K140 | 324-325 | **The effective crossing drifts with distance** (+0.33% at `r_e = 0`, larger at higher `r_e`), so every threshold uses the asymptotic `d ≥ 7` fit, which drops `d ∈ {3,5}`. | statistics | open | |
| K141 | 326-327 | **ν is only weakly constrained at `d ≤ 11`** (three large distances, coarse grid); `p_th` is far better determined than `ν`, and the fit reports χ²/dof so this is visible. | statistics | open | |
| K142 | 328-330 | **`estimate_crossing` is unreliable on saturated data**; it only seeds the fit window, and every quoted threshold comes from the `d ≥ 7` collapse fit, never the raw crossing estimate (mechanism in [`docs/AUDIT.md`](docs/AUDIT.md#threshold-fitting-methodology-reference)). | statistics | open | |
| K143 | 336-337 | `shortest_graphlike_error()` vs `d` from the builder (no Monte Carlo) | scope | open | |
| K144 | 339-341 | The DEM partition was verified by hand on the `d=3, T=2` instance and reconciled against Stim (incl. a non-obvious two-erasure/one-detector cancellation): | scope | open | |
| K145 | 343-344 | Measured on an **Apple M2 (8 cores, macOS 14.5), single-threaded, otherwise idle**, `d=5, p=3e-2`, 10⁴ shots: | performance | open | |
| K146 | 346-347 | Over 10 consecutive runs (2026-09-06) the *observed* spreads were 77.3k–85.7k, 1,529–1,585 and 48.8–54.4× respectively (sd 2.5k, 17, 1.6). | performance | open | |
| K147 | 347-349 | **Those spreads are what 10 runs happened to produce, not a predictive interval** — an 11th run may fall outside them, and the previously published min–max (77.4k–84.2k, 1,542–1,593) was breached at both ends by this very sample. | performance | open | |
| K148 | 349-350 | Quote the median; treat the range as a scatter estimate. | performance | open | |
| K149 | 352-354 | Both the rates *and the gap* are hardware- and load-dependent — on this same machine under full CPU oversubscription the gap ranged 55.8–65.3× — so treat the ~51× as one machine's number, not a property of the decoder. | performance | open | |
| K150 | 354-356 | Replacing the per-shot matcher rebuild with per-signature vectorized construction narrowed the gap from ~69× to ~51×; it did **not** bring it under the 50× target that motivated the work. | performance | open | |
| K151 | 356-357 | (An earlier revision of this bullet claimed 48×; that was 80,000/1,650, a ratio of two independently rounded numbers rather than a measured one.) | performance | open | |
| K152 | 388 | # (byte-identical to the committed PNGs *within one environment* -- see note below): | reproducibility | open | |
| K153 | 391 | # reproduce the analyses (wall-clock measured single-core on an Apple M2, 2026-09-06): | performance | open | |
| K154 | 399-401 | **On "byte-stable".** `tests/test_plotting.py::test_figures_regenerate_byte_stable` pins *run-to-run* determinism: the same CSVs rendered twice in the same environment give byte-identical PNGs. | reproducibility | open | |
| K155 | 401 | That is **not** a cross-environment guarantee. | reproducibility | open | |
| K156 | 401-403 | PNG bytes depend on the matplotlib/FreeType build, and `lambda_vs_p.png` draws bootstrap CI bars, which drift with the NumPy RNG exactly as the threshold CIs do. | reproducibility | open | |
| K157 | 403-406 | Verified 2026-09-06 to regenerate byte-identical to the committed PNGs (all four, incl. `lambda_vs_p.png` md5 `442e7f13…`) on **macOS 14.5 arm64, Python 3.14.6, NumPy 2.5.1, SciPy 1.18.0, stim 1.16.0, pymatching 2.4.0, sinter 1.16.0** | reproducibility | open | |
| K158 | 406-410 | matching "env B" of the [`docs/AUDIT.md`](docs/AUDIT.md#reproducibility-and-environment-sensitivity-2026-09-06) reproducibility table on every version that table records, plus **matplotlib 3.11.0**, which that table does not record and which is the dependency PNG bytes are most sensitive to. | reproducibility | open | |
| K159 | 410-411 | On a different matplotlib or NumPy, expect the figures to be visually identical but not byte-identical. | reproducibility | open | |
| K160 | 413-414 | `ablation_table.py` is the slow one: it decodes 100k shots per cell over `d = 3..11` x `r_e = 0.5, 0.98` with no CSV shortcut, and `d = 9`/`d = 11` dominate. | performance | open | |
| K161 | 414 | `--shots` trades runtime for the >= 50-error gate. | performance | open | |
| K162 | 416-418 | Re-collecting the real sweeps is resumable (remove `data/<name>.csv` first to start clean); rough wall-clock on 8 cores is ~10–20 min (`baseline_pauli`), ~1–2 h (`erasure_r50`), ~4–8 h (`erasure_r98`, slow-path dominated). | performance | open | |
| K163 | 423-425 | Future directions (imperfect/delayed heralds, other hardware, logical algorithms, alternative decoders, and a denser high-`r_e` sweep to resolve the `r_e = 0.98` herald threshold) are in [`docs/FUTURE_WORK.md`](docs/FUTURE_WORK.md). | scope | open | |

## docs/AUDIT.md caveats not repeated in README.md

| ID | Line | Text (verbatim) | Category | Status | Destination |
|---|---|---|---|---|---|
| A001 | README 40-41; AUDIT 562-564 | "land exactly on sweep grid points" above was checked at 1e-6 tolerance (`round(u, 6) in grid`), not at float precision. | statistics | open | |
| A002 | README 40-41; AUDIT 564-566 | At full precision the endpoints agree with the grid values to ~1 ULP (\|diff\| 7e-18 to 1.4e-17), **not** bit-identically -- the bounded least-squares solver returns a pinned parameter to within rounding. | statistics | open | |
| A003 | README 40-43; AUDIT 733-734 | † **Not reproducible from this tree (2026-09-07).** Both r_e = 0.98 values are kept above as published; neither regenerates here. | reproducibility | open | |
| A004 | README 40-43; AUDIT 743 | The discrepancy is **unexplained**. | reproducibility | open | |
| A005 | README 87; AUDIT 755-758 | The **excluding-edge-pins** column reproduces exactly *nowhere*: the r_e = 0.5 rows land at 1.601160% (vs 1.601044%), 1.592775% (vs 1.591298%) and 2.761548% (vs 2.761549%) | reproducibility | open | |
| A006 | README 85-87, 91-93; AUDIT 801-803 | Reproduce: `uv run python scripts/paired_separation.py --full-stability` for the 10-seed table; the 30-seed enumeration and the falsification test are the probes described above. | reproducibility | open | |
| A007 | README 83; AUDIT 785-788 | Measured rates agree (73/1000 discarded as `blind:bound_pin` in the Δ bootstrap; 74/998 lower-edge pins in the marginal bootstrap) though the RNG streams differ, so these are the same mechanism at the same rate, not literally the same replicates. | statistics | open | |
| A008 | README 85; AUDIT 707-709 | Individual bootstrap draws *do* exceed it (max draw 3.14%), so this is not a hard cap on the draws; the **percentile** is quantised onto discrete atoms. | statistics | open | |
| A009 | README 36-39, 79, 320-322; AUDIT 809-811 | That conclusion was measured for r_e = 0.5, where it holds. | statistics | open | |
| A010 | README 36-39, 79, 320-322; AUDIT 810-811 | **It does not hold for r_e = 0.98 blind**, whose `resolved = False` verdict rests on the asymmetry. | statistics | open | |
| A011 | README 36, 79; AUDIT 861 | \| `resolved` (guard = 1.0) \| **0 of 30** \| **9 of 30** \| | statistics | open | |
| A012 | README 36, 79; AUDIT 869-871 | At `seed = 0` — the seed this repo declares everywhere — the pins-excluded `rel_ci_width` is **1.0186**, marginally **above** the 1.0 `_MAX_REL_CI_WIDTH` guard. | statistics | open | |
| A013 | README 39; AUDIT 882-884 | With pins retained, `rel_ci` is 2.4–4.3 and the fit misses the guard by a factor of 2.4–4.3 — which reads as a settled, comfortable non-result. | statistics | open | |
| A014 | README 39; AUDIT 884-886 | With pins excluded, `rel_ci` clusters on the guard itself (median 1.05, range 0.39–1.36) and the verdict is a coin flip in the seed: 9 of 30 resolving is **instability, not resolution**. | statistics | open | |
| A015 | README 39; AUDIT 888-889 | What the current wide margin misrepresents is *how* unresolvable: it makes a marginal, seed-dependent call look like a decisive one. | statistics | open | |
| A016 | README 72-73; AUDIT 363-366 | The constant is chosen from the committed fits, not by taste: every reported result has relative width ≤ 0.49 (`r_e = 0.5` blind), while the `r_e = 0.98` blind non-result is ≈ 3.0 — a ~6× gap, so 1.0 rejects only that fit. | statistics | open | |
| A017 | README 72-73; AUDIT 96-97 | The fit requires ≥ 3 degrees of freedom (≥ 8 points for 5 params), so it cannot "converge" on a dof≈1 window. | statistics | open | |
| A018 | README 75-79; AUDIT 334-336 | Points at or above per-round `p_L = 0.4` (saturating toward ½, outside the local ansatz) are excluded, and the fit reports its χ²/dof. | statistics | open | |
| A019 | README 75-79; AUDIT 326-328 | at `d = 11` (`T = 11`) a legitimate near-crossing point still has `P_L_shot ≈ 0.45`, so the saturation cut acts on `p_L`, not `P_L_shot`. | statistics | open | |
| A020 | README 75-79; AUDIT 107-108 | The result is cap-insensitive for caps in {0.2, 0.25, 0.3, 0.4}. | statistics | open | |
| A021 | README 69; AUDIT 121-123 | Calibration is checked by `test_bootstrap_ci_covers_known_p_th_at_nominal_rate` (90% coverage of a known p_th over 40 seeded realizations). | statistics | open | |
| A022 | README 29, 78; AUDIT 131-134 | the `r_e = 0.5` **herald** CI reads `[2.2, 2.6]` here at `n_boot = 200` and `[2.2, 2.8]` at `n_boot = 1000`. | statistics | open | |
| A023 | README 29, 78; AUDIT 701 | \| r_e = 0.5 herald \| 2.704709% \| 2.833190% \| **no** (nearest 1e-3 away) \| | statistics | open | |
| A024 | README 69-71; AUDIT 531-534 | "seed spread" is the range over seeds 0–9 within env B — the Monte-Carlo error of the bootstrap itself, which is the larger effect everywhere except the r_e = 0.5 Δ lower endpoint. | reproducibility | open | |
| A025 | README 54-55, 71-72; AUDIT 519-521 | The last two rows are RNG-free by construction (`curve_fit` and a residual sum on fixed data), so they are expected to be exact, but env A did not report them and that is not claimed as measured here. | reproducibility | open | |
| A026 | README 54-55; AUDIT 606-608 | **May:** the analysis is deterministic given `(CSV, seed, NumPy build)`, and every point estimate, χ²/dof and verdict reproduces exactly across the two environments tested. | reproducibility | open | |
| A027 | README 54-55; AUDIT 500-502 | Only the analysis path matters here — no sampling is re-run — so NumPy and SciPy are the operative versions. | reproducibility | open | |
| A028 | README 54-55; AUDIT 538-539 | The env-A column is an independent reviewer's report and is not reproducible from this tree. | reproducibility | open | |
| A029 | README 101-102; AUDIT 176-178 | the reported herald/blind bootstrap correlation on real data is ~0.02–0.04, confirming pairing would have bought essentially nothing here. | statistics | open | |
| A030 | README 97-98; AUDIT 209-211 | lower-endpoint spread 0.036%, upper-endpoint spread 0.018%; **all five seeds exclude zero**. | statistics | open | |
| A031 | README 102-103; AUDIT 248-250 | So at r_e=0.5 the discards **are** biased: they are disproportionately the replicates where blind lands high (bound-pinning near the top of the grid), which are the replicates with Δ nearest zero. | statistics | open | |
| A032 | README 102-103; AUDIT 250-251 | Excluding them nudges the CI away from zero. | statistics | open | |
| A033 | README 104-105; AUDIT 288-290 | The 97.5th percentile of 1000 draws (numpy 'linear' convention) reaches zero once **26** are ≥ 0, so **the claim survives unless ≥ 24 of the 102 discards would have given Δ ≥ 0** (verified against a direct construction in the tests). | statistics | open | |
| A034 | README 106; AUDIT 295-296 | Forming Δ = blind − herald from the recorded values (filling the 1 missing herald from the successful median) | statistics | open | |
| A035 | README 106; AUDIT 300-301 | Because a blind bound-pin `p_th` is a *lower* bound on the true crossing, 2 is itself a lower bound on the count | statistics | open | |
| A036 | README 99; AUDIT 313-314 | The r_e=0 control has tipping point 0 (its Δ already includes zero, as a control must) and missing-at-random discards. | statistics | open | |
| A037 | README 251-254; AUDIT 436-437 | Wilson coverage is nominal, not exact; at single-digit counts the "95%" is an approximation with no guaranteed coverage. | statistics | open | |
| A038 | README 276; AUDIT 678-680 | The estimator still *excludes* the zero-denominator replicates rather than reweighting them or propagating `+∞`. | statistics | open | |
| A039 | README 276, 279; AUDIT 660-663 | The `≥ 50`-observed-error gate in `figure_lambda` (and the same gate in `scripts/ablation_table.py`) admits 425 pairs, and **every one of them drops zero replicates** — verified at two different seeds. | statistics | open | |
| A040 | README 276; AUDIT 648-649 | The count is itself mildly seed-dependent at the margin — 71 / 43 at seed 12345. | statistics | open | |
| A041 | README 200-202; AUDIT 397-399 | No threshold, ablation or separation result depends on this metric — it is descriptive only, so no fitted number moves. | scope | open | |
