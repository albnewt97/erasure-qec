"""Herald-conditioned decoder correctness (PLAN.md §8, M6 gate).

The headline test forces a specific erasure pattern with probability-1 probe
erasures and shows the herald-aware decoder correcting a pattern the blind
decoder provably fails:

  Erase (1,3) and (1,5) — the two OFF-row qubits of the x=1 vertical column —
  every shot. When both suffer an X-like error (1/4 of shots), the syndrome is
  the single detector D(2,2,1): the two X's meet at the Z-check (0,4) (fires
  twice, cancels) and terminate on the top boundary, leaving only (2,2)'s
  round-1 detector. That syndrome is EXACTLY the signature of a single X on
  (1,1) — which lies ON logical_z_support(3) and would flip Z-bar. Only the
  herald bits distinguish the two:

  - blind (static weights, erasures folded in at their small unconditional
    marginal): the one-edge L0-flipping explanation via (1,1) is cheaper than
    the two erased edges -> predicts a flip -> WRONG (truth: no flip; neither
    erased qubit is on the logical row).
  - herald-aware: both erased edges drop to weight 0 -> the two-edge no-flip
    route costs 0 -> CORRECT.

All samplers are seeded: every assertion below is deterministic and was
verified against a live run before being frozen.
"""

import numpy as np
import pytest
import stim

from erasure_qec.circuits.builder import build
from erasure_qec.config import NoiseParams
from erasure_qec.decoding.dem_partition import partition_dem
from erasure_qec.decoding.herald_matching import (
    BlindMatchingDecoder,
    HeraldMatchingDecoder,
)
from erasure_qec.noise.injector import (
    BiasedErasureInjector,
    NullInjector,
    PauliOnlyInjector,
)

# The forced pattern: probability-1 erasures on the two off-row column qubits.
PROBES = [(1 + 3j, 0), (1 + 5j, 0)]


def _sampling_circuit() -> stim.Circuit:
    """Heralds forced every shot; the probes' Paulis are the only randomness."""
    return build(3, 2, NullInjector(), probe_erasures=PROBES, probe_q=1.0)


def _decoding_circuit() -> stim.Circuit:
    """Same structure, decoder-side noise model: background Pauli noise gives
    every graph edge a positive weight; probes at LOW q (0.02) so the blind
    decoder's static weights treat the erased locations as unlikely."""
    return build(
        3, 2, PauliOnlyInjector(NoiseParams(p=1e-3)), probe_erasures=PROBES, probe_q=0.02
    )


def _det_at(circuit: stim.Circuit, x: float, y: float, t: float) -> int:
    coords = circuit.get_detector_coordinates()
    (idx,) = [i for i, c in coords.items() if c == [x, y, t]]
    return int(idx)


def test_forced_erasure_pattern_blind_fails_herald_corrects() -> None:
    """§8 correctness gate: the herald-aware decoder corrects a forced-erasure
    pattern the blind decoder fails."""
    sampling = _sampling_circuit()
    decoding = _decoding_circuit()
    # Same detector layout: only the noise/probability arguments differ.
    assert sampling.num_detectors == decoding.num_detectors == 18
    assert sampling.get_detector_coordinates() == decoding.get_detector_coordinates()

    herald = HeraldMatchingDecoder.from_circuit(decoding)
    blind = BlindMatchingDecoder(decoding)

    n = 2048
    dets, obs = sampling.compile_detector_sampler(seed=2026).sample(
        n, separate_observables=True
    )
    # Neither erased qubit is on the logical row -> the truth is NEVER a flip.
    assert not obs.any()

    herald_pred = herald.decode_batch(dets)
    blind_pred = blind.decode_batch(dets)
    herald_fails = (herald_pred != obs).any(axis=1)
    blind_fails = (blind_pred != obs).any(axis=1)

    # Herald-aware decodes every shot correctly; blind fails on ~1/4 of shots
    # (both probes X-like). Verified live: 0 vs 475 of 2048 with this seed.
    assert herald_fails.sum() == 0
    assert 0.15 * n < blind_fails.sum() < 0.35 * n

    # The specific corrected pattern, pinned on the first blind failure:
    # both heralds fired, syndrome == {D(2,2,1)} on the Z-check side, and the
    # herald-aware decoder corrects exactly where blind fails.
    i = int(np.flatnonzero(blind_fails)[0])
    d_z22_r1 = _det_at(sampling, 2, 2, 1)
    d_z04_r1 = _det_at(sampling, 0, 4, 1)
    assert dets[i, 4] and dets[i, 5]  # herald detectors (1,3,0,1) and (1,5,0,1)
    assert dets[i, d_z22_r1] and not dets[i, d_z04_r1]  # the ambiguous syndrome
    assert blind_pred[i].any() and not herald_pred[i].any()  # blind flips, herald doesn't
    assert not obs[i].any()  # truth: no flip


def test_fast_path_matches_static_matcher_when_no_heralds_exist() -> None:
    """With no heralds anywhere, every shot takes the fast path and the decoder
    must agree exactly with the plain static matcher on the same DEM."""
    circuit = build(3, 2, PauliOnlyInjector(NoiseParams(p=5e-3)))
    decoder = HeraldMatchingDecoder.from_circuit(circuit)
    static = BlindMatchingDecoder(circuit)  # full DEM == dem_pauli here
    assert partition_dem(circuit).herald_indices.tolist() == []

    dets, _ = circuit.compile_detector_sampler(seed=7).sample(512, separate_observables=True)
    assert (decoder.decode_batch(dets) == static.decode_batch(dets)).all()


def test_conditional_only_edges_decode_pure_probe_circuit() -> None:
    """Decoder compiled from the probe circuit itself: dem_pauli has ZERO error
    mechanisms, so the base graph has no edges and every heralded shot is
    decoded purely on conditionally-added weight-0 edges (§8: "adding the edge
    to the graph if it exists only conditionally"). With heralds identifying
    every fault location, decoding is perfect."""
    sampling = _sampling_circuit()
    partition = partition_dem(sampling)
    assert partition.dem_pauli.num_errors == 0

    decoder = HeraldMatchingDecoder(partition)
    dets, obs = sampling.compile_detector_sampler(seed=11).sample(
        512, separate_observables=True
    )
    preds = decoder.decode_batch(dets)
    assert (preds == obs).all()


def test_two_tier_dispatch_batch_equals_per_shot_decoding() -> None:
    """Mixed batch (some shots heralded, some not) must produce identical
    results whether decoded together or one shot at a time — i.e. the fast/slow
    dispatch is transparent."""
    circuit = build(3, 2, BiasedErasureInjector(NoiseParams(p=0.02, r_e=0.5)))
    decoder = HeraldMatchingDecoder.from_circuit(circuit)
    partition = partition_dem(circuit)

    n = 256
    dets, _ = circuit.compile_detector_sampler(seed=13).sample(n, separate_observables=True)
    heralded = dets[:, partition.herald_indices].any(axis=1)
    # The sample genuinely exercises both tiers (verified live: 73 of 256 heralded).
    assert 0 < int(heralded.sum()) < n

    batch = decoder.decode_batch(dets)
    singles = np.vstack([decoder.decode_batch(dets[i : i + 1]) for i in range(n)])
    assert (batch == singles).all()


def test_herald_bits_alone_predict_identity() -> None:
    """The I-branch: heralds fired but no syndrome -> no correction."""
    decoder = HeraldMatchingDecoder.from_circuit(_decoding_circuit())
    shot = np.zeros((1, decoder.num_detectors), dtype=bool)
    shot[0, 4] = True  # herald (1,3,0,1)
    shot[0, 5] = True  # herald (1,5,0,1)
    assert not decoder.decode_batch(shot).any()


def test_ablation_herald_beats_blind_on_biased_erasure_circuit() -> None:
    """Smoke ablation on a real biased-erasure circuit: with r_e=0.9 the herald
    information is most of the error budget, so the herald-aware decoder must
    make strictly fewer logical errors than the blind one on identical shots.
    Verified live with this seed: 371 vs 544 failures of 2048."""
    circuit = build(3, 3, BiasedErasureInjector(NoiseParams(p=0.03, r_e=0.9)))
    herald = HeraldMatchingDecoder.from_circuit(circuit)
    blind = BlindMatchingDecoder(circuit)

    dets, obs = circuit.compile_detector_sampler(seed=17).sample(
        2048, separate_observables=True
    )
    herald_fails = int((herald.decode_batch(dets) != obs).any(axis=1).sum())
    blind_fails = int((blind.decode_batch(dets) != obs).any(axis=1).sum())
    assert herald_fails < blind_fails


def test_rejects_wrong_shape() -> None:
    decoder = HeraldMatchingDecoder.from_circuit(_decoding_circuit())
    with pytest.raises(ValueError, match="shape"):
        decoder.decode_batch(np.zeros((4, decoder.num_detectors + 1), dtype=bool))
