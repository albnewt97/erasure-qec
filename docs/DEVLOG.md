# Development log

A narrative of the decisions and the bugs, for a reader who won't open the source. Every claim points
at something committed — a test, a CSV, a figure, or the hand-verification worksheet.

## Why this problem

I wanted to answer one question: if a fraction of the gate errors are *heralded* — the hardware flags
which qubit was disturbed and when — how much does that buy a surface-code memory, and how much of the
benefit is the physics versus the decoder? The hardware context is current: ¹⁷¹Yb neutral-atom Rydberg
gates and dual-rail superconducting cavities both natively flag leakage and photon loss. Erasure
conversion turns those flags into `HERALDED_ERASE` events. The interesting part is the
herald-conditioned decoder that reads the herald bits per shot and zeroes the corresponding
matching-graph edges. That part is not available off the shelf: PyMatching decodes a fixed graph, so the
per-shot reweighting and the DEM bookkeeping behind it were the real work. (I originally wrote that the
"passive benefit — a depolarized qubit is cheaper to correct than the Pauli budget it replaced — is
free"; that was hand-waving. Both the residual DEPOLARIZE2 and the erasure's I/2 are Pauli channels, and
once the per-gate budget is held constant the blind decoder barely improves at all — see "The result".)

## Noiseless first

I refused to add any noise until the noiseless circuit was provably correct. The M2 gate
(`tests/test_detectors_deterministic.py::test_all_detectors_deterministically_zero`) samples 1000 shots
of the noiseless circuit and asserts every detector stays zero. This sounds trivial; it isn't. Every
detector-algebra mistake — a wrong measurement-record offset, a miscounted round, a detector wired to
the wrong ancilla — surfaces here as a spurious detection event, for free, with no statistics to wade
through. By the time noise went in, the skeleton was known-good.

## The schedule is load-bearing

The order in which each ancilla touches its data neighbours decides whether a mid-round fault becomes a
harmless weight-2 error or a "hook" that halves the code distance. The hook-safe order (X-checks in
Z-shape, Z-checks in ᴎ-shape) orients hooks perpendicular to the logical operator. I proved this matters
with a deliberately-broken fixture: `tests/test_distance_invariants.py` builds the code with both bases
sharing one order and asserts `shortest_graphlike_error()` collapses from `d` to `⌈d/2⌉` — pinned
exactly as `{3:2, 5:3, 7:4}`. `figures/hook_regression.png` shows the two curves diverging. Without that
regression I'd have no evidence the schedule does anything at all.

## Hand-verifying the DEM partition

The one component I did not trust to code review is the split of Stim's decomposed detector error model
into a Pauli sub-DEM plus a per-herald table of conditioned edges. I verified it on paper first, on the
smallest nontrivial circuit (d=3, T=2, one erasure), across three probes — a centre data qubit, a corner
(boundary) qubit, and a mid-round ancilla — predicting each DEM line by hand and reconciling against
Stim line by line (`docs/dem_worksheet.md`). The payoff was Probe C: a Z error copied onto two data
qubits fired only one detector, purely because of CX layer ordering — one neighbour's check reads the
error in a later layer (fires), the other had already read it before injection (never sees it). That
is the concrete demonstration that hook analysis is schedule-dependent, not a slogan. The visual
evidence is committed at `docs/figures/probe_mid_round_ancilla_2_2_detslice.svg`.

## The Design-2 partition decision

A heralded Y error decomposes in Stim as `X-component ^ Z-component ^ herald`. I chose to split at the
`^` separators and emit one conditioned edge per component, rather than one hyperedge per whole
mechanism (`src/erasure_qec/decoding/dem_partition.py`). The reason is that PyMatching is graphlike: a
4-detector Y hyperedge can never enter the matching graph, but once its two constituent component edges
both go to weight zero, the Y case is covered automatically. Whole-mechanism hyperedges would have
needed special handling the graphlike matcher can't use anyway.

## The throughput optimization, honestly

At r_e=0.98 nearly every shot is heralded, so the per-shot slow path dominates. I measured before
optimizing: ~80,000 shots/s on the herald-free fast path versus ~1,235 shots/s rebuilding a matcher per
shot with Python `add_edge` calls — a 69× gap (`README.md`, "Decoder throughput"). My first idea was to
group shots by herald *signature* and build one matcher per distinct signature. It barely helped: at
d=5 the expected number of fired heralds is ~12 of ~800, so signatures are near-unique and there's
almost nothing to amortize. The real win came from making each construction cheap — building the
reweighted matcher with a vectorized `pymatching.Matching.from_check_matrix` over a precompiled sparse
check matrix instead of per-edge Python loops (`src/erasure_qec/decoding/herald_matching.py`). That got
it to ~1,650 shots/s (48×). Grouping stayed because it's free, but it wasn't the lever.

## Three bugs in my own analysis code

`estimate_crossing` — which locates the threshold to seed the fit — failed three times, each documented
in `src/erasure_qec/analysis/threshold_fit.py`. First it returned the low edge of the p-grid: deep
sub-threshold every curve vanishes, so the *absolute* spread across distances is trivially smallest
there. Fix: minimize *relative* spread, normalized by the mean. Then it returned an *above*-threshold
value on the r_e=0.5 sweep, where the saturated high-p tail (all curves near ½) faked a second flat
region. Fix: an ordering-inversion guard that rejects any p where p_L already increases with d — still
imperfect, and the docstring says so. Third, the fit window wasn't excluding saturated points at all; it
does now. Both crossing bugs are pinned by regression tests:
`tests/test_threshold_fit.py::test_estimate_crossing_on_real_baseline_is_near_threshold` and
`::test_estimate_crossing_r50_not_pulled_by_saturated_tail`.

## Finite-size drift, diagnosed not hidden

Including the small distances (d=3,5) pulls the effective crossing up by +0.3–0.4% and inflates ν toward
~2.5, because small codes carry the largest finite-size corrections. Rather than quietly pick whichever
fit looked best, every threshold I quote uses the d≥7 collapse fit (`fit_threshold(..., d_min=7)`), and
the drift is written down as a systematic in `README.md` under "Caveats and limitations."

## An external audit, and what it cost me

An external reviewer checked the repo and found the circuit builder, DEM partition, and the
herald-vs-blind ablation correct — but two things wrong that the tests hadn't caught, both of which had
inflated my headline. `docs/AUDIT.md` records the reproduction and the fixes; the short version:

**The noise budget shrank with r_e.** I sized the heralded-erase rate at `p·r_e/2`, forgetting that an
erasure is a non-identity error only ¾ of the time. So the total per-gate non-identity probability fell
from `p` to `~0.75p` as `r_e → 0.98`. Part of my "rising threshold" was just *less noise* at high r_e —
not a measurement of erasure conversion at all. Fix: rate `(2/3)·p·r_e`, so `2·q·¾ = p·r_e` and the
per-gate budget is held at `p` for every r_e (commit `6b6b87d`). This invalidated every erasure sweep on
disk, so I re-collected r_e=0.5 and r_e=0.98 under the fixed model and committed them.

**The fit reported bound-pinned results as converged.** `curve_fit` not raising was treated as success,
so an r_e=0.98 fit sitting exactly on the top of its p-range came back "converged" with a clean-looking
number. And I was using the 95% Wilson half-width as the 1σ weight, deflating χ²/dof by ~3.8× so nothing
ever looked over-parameterised. Fixes (commits `0ce6ff5`, `df70ddc`): report bound-pinning as
`converged=False`; use the real 1σ error; require ≥3 dof; add a `chi2_dof` field. I also tried a
shot-level saturation cut and had to revert it — at large T a legitimate near-crossing point still has
`P_L_shot → ½`, so the cut deleted the near-crossing large-d data and broke exactly the high-threshold
fit I cared about. The cap belongs on the per-round rate.

## A fourth bug: the bootstrap measured the wrong estimator

One more that the audit didn't list but the same scrutiny caught: my CI was too narrow because the
bootstrap wasn't the same estimator as the fit. It refit an *unweighted* model, starting from the
point-estimate parameters, on a *frozen* window, and silently dropped replicates that failed — all four
shrink the reported σ. I rewrote it to run the whole pipeline per replicate (resample every point,
re-run the crossing estimate, the window, and the *weighted* fit), count failures, and report a
percentile CI since the p_th distribution is skewed. A coverage test now checks the 95% CI actually
covers a known p_th ~95% of the time. This mattered: the r_e=0.98 *blind* fit that looked like
`1.64% ± 0.08%` has an honest 95% CI spanning roughly `[1.5, 8]%` and shifting with the seed — its
crossing is bistable, so it is not a resolved threshold either.

## The result, honestly

From the re-collected fixed-model sweeps in `data/`, `figures/threshold_panels.png` (d≥7 collapse fits;
95% bootstrap percentile CI, seed 0):

| r_e | herald | blind | χ²/dof |
|---|---|---|---|
| 0 | 1.38% `[1.33, 1.47]` | 1.44% `[1.35, 1.57]` | ~1.0 |
| 0.5 | 2.32% `[2.19, 2.64]` | 1.49% `[1.42, 2.15]` | 0.6 / 0.5 |
| 0.98 | **not resolved** | **not resolved** `[1.5, ~8]` | — / 0.20 |

The r_e=0 row is a control: with no heralds the two decoders agree within their CIs, and they do. The
honest story is the *opposite* of what I first wrote. Where the fit resolves, with the budget held
constant the **blind** threshold barely moves (1.44% → 1.49% from r_e=0 to 0.5) — erasure conversion on
its own is close to a wash. The gain is almost entirely **herald-conditioning**: 1.38% → 2.32%. At
r_e=0.98 the collapse fit resolves *neither* threshold (herald rails ν to its bound; blind is bistable),
so I quote neither, rather than the 6.27% I once did. The robust high-r_e evidence is the ablation
(`scripts/ablation_table.py`, identical shots, no fit): the herald advantage grows to ~1000× at d=11 — a
low-statistics lower bound, but unambiguous in direction. That ablation is the thing worth building; the
threshold table is honest about where it can and cannot resolve a number.
