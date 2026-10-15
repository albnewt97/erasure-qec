"""Wu et al. decoding circuit, tagged partition and decoders (P4 gate, T1-T4).

Seeded stim samplers reproduce only on the same machine and stim version (see
tests/test_herald_decoder.py), so every assertion here is structural or a band
sized from binomial standard errors; draw-dependent counts in comments carry
the machine and versions they were measured on.
"""

import itertools
import math
from collections.abc import Mapping

import numpy as np
import numpy.typing as npt
import pymatching
import pytest
import stim

from erasure_qec.circuits.builder import build
from erasure_qec.decoding.dem_partition import _split_components
from erasure_qec.decoding.herald_matching import BlindMatchingDecoder, HeraldMatchingDecoder
from erasure_qec.decoding.sinter_adapter import WU_DECODERS, wu_task
from erasure_qec.decoding.tagged_partition import partition_tagged, partition_tagged_dem
from erasure_qec.noise.wu import WuDecodingInjector, WuNoiseParams, WuSamplingInjector

SEED = 2026

_Edge = tuple[tuple[int, ...], tuple[int, ...]]


def _circuits(
    d: int, params: WuNoiseParams, overrides: Mapping[int, float] | None = None
) -> tuple[stim.Circuit, stim.Circuit]:
    o = dict(overrides or {})
    sampling = build(d, d, WuSamplingInjector(d, params, location_overrides=o))
    decoding = build(d, d, WuDecodingInjector(d, params, location_overrides=o))
    return sampling, decoding


def _decoding_dem(circuit: stim.Circuit) -> stim.DetectorErrorModel:
    return circuit.detector_error_model(decompose_errors=True, flatten_loops=True).flattened()


def _pack(dets: npt.NDArray[np.bool_]) -> npt.NDArray[np.uint8]:
    return np.packbits(dets, axis=1, bitorder="little")


def _unpack_obs(packed: npt.NDArray[np.uint8], num_obs: int) -> npt.NDArray[np.bool_]:
    return np.unpackbits(packed, axis=1, count=num_obs, bitorder="little").astype(bool)


# --- Amendments A and B: detector layout and shared location enumeration. ---


@pytest.mark.parametrize(
    ("params", "overrides"),
    [
        (WuNoiseParams(0.02, 0.5), {}),
        (WuNoiseParams(0.0, 0.0), {3: 0.5}),
        (WuNoiseParams(0.01, 0.0), {3: 0.5, 7: 0.0}),
    ],
    ids=["no_overrides", "pe0_one_override", "pe0_override_with_zero"],
)
def test_sampling_and_decoding_detector_layouts_match(
    params: WuNoiseParams, overrides: dict[int, float]
) -> None:
    sampling, decoding = _circuits(3, params, overrides)
    assert sampling.num_detectors == decoding.num_detectors
    assert sampling.get_detector_coordinates() == decoding.get_detector_coordinates()
    herald_coords = [c for c in decoding.get_detector_coordinates().values() if len(c) == 4]
    n_pairs = sum(
        len(inst.targets_copy()) // 2
        for inst in decoding.flattened()
        if isinstance(inst, stim.CircuitInstruction) and inst.name == "CX"
    )
    assert len(herald_coords) == n_pairs  # every pair flagged, active or not


@pytest.mark.parametrize("d", [3, 5])
def test_sampling_and_decoding_enumerate_locations_identically(d: int) -> None:
    params = WuNoiseParams(0.02, 0.5)
    sampling_inj = WuSamplingInjector(d, params)
    decoding_inj = WuDecodingInjector(d, params)
    build(d, d, sampling_inj)
    build(d, d, decoding_inj)
    assert sampling_inj.locations
    assert sampling_inj.locations == decoding_inj.locations


def test_partition_tagged_rejects_mismatched_layouts() -> None:
    sampling, _ = _circuits(3, WuNoiseParams(0.02, 0.5))
    _, decoding = _circuits(3, WuNoiseParams(0.02, 0.0))  # no flags -> fewer detectors
    with pytest.raises(ValueError, match="detector counts differ"):
        partition_tagged(sampling, decoding)


def test_partition_tagged_dem_rejects_malformed_tags() -> None:
    base = "detector(1, 1, 0, 1) D0\ndetector(2, 2, 0) D1\n"
    with pytest.raises(ValueError, match="unrecognised"):
        partition_tagged_dem(stim.DetectorErrorModel(base + "error[foo](0.1) D1"))
    with pytest.raises(ValueError, match="out of range"):
        partition_tagged_dem(stim.DetectorErrorModel(base + "error[wu_e1:Xa](0.1) D1"))
    with pytest.raises(ValueError, match="touches a herald"):
        partition_tagged_dem(stim.DetectorErrorModel(base + "error(0.1) D0 D1"))


# --- T1: herald_table[k] equals the edges of location k's generators alone. ---


def test_t1_herald_table_matches_isolated_location_edges() -> None:
    d = 3
    sampling, decoding = _circuits(d, WuNoiseParams(0.02, 1.0))
    partition = partition_tagged(sampling, decoding)
    n_locations = len(partition.herald_indices)
    chosen = np.random.default_rng(SEED).choice(n_locations, size=8, replace=False)

    for k in (int(k) for k in chosen):
        # Only location k's four generators, each at 0.1; all other noise zero.
        isolated = build(
            d, d, WuDecodingInjector(d, WuNoiseParams(0.0, 0.0), location_overrides={k: 0.2})
        )
        edges: set[_Edge] = set()
        for instr in _decoding_dem(isolated):
            if instr.type != "error":
                continue
            assert instr.tag.startswith(f"wu_e{k}:")
            for comp in _split_components(instr.targets_copy()):
                dets = tuple(sorted(t.val for t in comp if t.is_relative_detector_id()))
                obs = tuple(sorted({t.val for t in comp if t.is_logical_observable_id()}))
                if dets or obs:
                    edges.add((dets, obs))
        table = partition.herald_table[int(partition.herald_indices[k])]
        assert edges
        assert edges == {(e.dets, e.obs_mask) for e in table}


# --- T2: forced erasures at two locations (mirrors the M6 forced-erasure test). ---


def _first_discriminating_pair(
    d: int, n: int
) -> tuple[int, int, npt.NDArray[np.bool_], npt.NDArray[np.bool_]] | None:
    """First (k1, k2) in lexicographic order where blind fails on >= 10% of shots.

    Returns the pair with the herald and blind per-shot failure vectors, or
    ``None`` if no pair qualifies. Decoders are compiled from a decoding
    circuit at small, mis-specified ``p_e`` (0.005) as in M6.
    """
    decoding_params = WuNoiseParams(0.01, 0.5)  # p_e = 0.005, p_p = 0.005
    _, decoding = _circuits(d, decoding_params)
    dem = _decoding_dem(decoding)
    partition = partition_tagged_dem(dem)
    herald = HeraldMatchingDecoder(partition)
    blind = BlindMatchingDecoder.from_dem(dem, partition.herald_indices)

    for k1, k2 in itertools.combinations(range(len(partition.herald_indices)), 2):
        sampling, _ = _circuits(d, WuNoiseParams(0.0, 0.0), {k1: 1.0, k2: 1.0})
        assert sampling.get_detector_coordinates() == decoding.get_detector_coordinates()
        dets, obs = sampling.compile_detector_sampler(seed=SEED).sample(
            n, separate_observables=True
        )
        blind_fails = (blind.decode_batch(dets) != obs).any(axis=1)
        if blind_fails.sum() >= 0.1 * n:
            herald_fails = (herald.decode_batch(dets) != obs).any(axis=1)
            assert dets[:, partition.herald_indices].sum(axis=1).tolist() == [2] * n
            return k1, k2, herald_fails, blind_fails
    return None


def test_t2_forced_erasure_pair_blind_fails_herald_corrects() -> None:
    n = 2048
    found = _first_discriminating_pair(3, n) or _first_discriminating_pair(5, n)
    assert found is not None, "no location pair at d=3 or d=5 makes blind fail on >= 10%"
    _, _, herald_fails, blind_fails = found
    # Measured 2026-10-15 on arm64 / stim 1.16.0 / pymatching 2.4.0: first pair
    # (0, 2) at d=3, blind 514/2048, herald 0. Exhaustive d=3 scan: 845 of 2556
    # pairs give blind >= 10%; herald fails 0 shots on all 2556 pairs.
    assert blind_fails.sum() >= 0.1 * n
    assert herald_fails.sum() == 0


# --- T3: R_e = 0 -> wu_herald_mwpm equals plain PyMatching. ---


def test_t3_re_zero_herald_equals_plain_pymatching() -> None:
    d, n = 5, 5000
    task = wu_task(d, WuNoiseParams(0.01, 0.0), {"d": d})
    assert task.circuit is not None and task.detector_error_model is not None
    compiled = WU_DECODERS["wu_herald_mwpm"].compile_decoder_for_dem(dem=task.detector_error_model)
    plain = pymatching.Matching.from_detector_error_model(
        task.circuit.detector_error_model(decompose_errors=True)
    )
    dets = task.circuit.compile_detector_sampler(seed=SEED).sample(n)
    wu_pred = _unpack_obs(
        compiled.decode_shots_bit_packed(bit_packed_detection_event_data=_pack(dets)),
        task.circuit.num_observables,
    )
    plain_pred = np.asarray(plain.decode_batch(dets.astype(np.uint8)), dtype=bool)
    assert wu_pred.shape == plain_pred.shape == (n, 1)
    assert np.array_equal(wu_pred, plain_pred)


# --- T4: herald beats blind at d=5, p=0.03, R_e=0.98. ---


def test_t4_herald_beats_blind_by_five_standard_errors() -> None:
    d, n = 5, 20000
    task = wu_task(d, WuNoiseParams(0.03, 0.98), {"d": d})
    assert task.circuit is not None and task.detector_error_model is not None
    dets, obs = task.circuit.compile_detector_sampler(seed=SEED).sample(
        n, separate_observables=True
    )
    packed = _pack(dets)
    fails = {}
    for name in ("wu_herald_mwpm", "wu_blind_mwpm"):
        compiled = WU_DECODERS[name].compile_decoder_for_dem(dem=task.detector_error_model)
        pred = _unpack_obs(
            compiled.decode_shots_bit_packed(bit_packed_detection_event_data=packed),
            task.circuit.num_observables,
        )
        fails[name] = int((pred != obs).any(axis=1).sum())
    ler_h = fails["wu_herald_mwpm"] / n
    ler_b = fails["wu_blind_mwpm"] / n
    se = math.sqrt(ler_h * (1 - ler_h) / n + ler_b * (1 - ler_b) / n)
    assert ler_b - ler_h > 5 * se, (ler_h, ler_b, se)
