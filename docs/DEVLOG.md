# Development log

A narrative of the decisions and the bugs, for a reader who won't open the source. Every claim points
at something committed — a test, a CSV, a figure, or the hand-verification worksheet.

## Why this problem

I wanted to answer one question: if a fraction of the gate errors are *heralded* — the hardware flags
which qubit was disturbed and when — how much does that buy a surface-code memory, and how much of the
benefit is the physics versus the decoder? The hardware context is current: ¹⁷¹Yb neutral-atom Rydberg
gates and dual-rail superconducting cavities both natively flag leakage and photon loss. Erasure
conversion turns those flags into `HERALDED_ERASE` events. The passive benefit — a depolarized qubit is
cheaper to correct than the Pauli budget it replaced — is free. The interesting part is the
herald-conditioned decoder that reads the herald bits per shot and zeroes the corresponding
matching-graph edges. That part is not available off the shelf: PyMatching decodes a fixed graph, so the
per-shot reweighting and the DEM bookkeeping behind it were the real work.

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

## The result

From the real sweeps in `data/`, plotted in `figures/threshold_panels.png` (d≥7 collapse fits):

| r_e | herald | blind |
|---|---|---|
| 0 | 1.38% | 1.44% |
| 0.5 | 2.30% | 1.73% |
| 0.98 | 6.27% | 2.29% |

The r_e=0 row is a control: with no heralds the two decoders must agree, and they do. The threshold
rises 1.4% → 6.3%, a ≈4.5× gain. Erasure conversion alone — the blind decoder at r_e=0.98 — accounts for
1.4% → 2.3%; the herald-conditioned decoder does the rest, 2.3% → 6.3%. The decoder is more than half the
total gain, which is the whole reason it was worth building.
