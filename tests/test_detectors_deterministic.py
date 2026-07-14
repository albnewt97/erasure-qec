"""Determinism + structural checks for the noiseless builder (PLAN.md §3.5)."""

from collections.abc import Callable

import pytest
import stim

from erasure_qec.circuits.builder import build, num_detectors
from erasure_qec.noise.injector import NoiseInjector, NullInjector

DISTANCES = [3, 5, 7]

NoiseFactory = Callable[[float], NoiseInjector]


def _cases() -> list[tuple[int, int]]:
    """(distance, rounds) pairs: T=1 edge, a small bulk, and T=d."""
    cases: list[tuple[int, int]] = []
    for d in DISTANCES:
        for t in sorted({1, 3, d}):
            cases.append((d, t))
    return cases


CASES = _cases()


@pytest.mark.parametrize("d,rounds", CASES)
def test_builds_and_compiles(d: int, rounds: int) -> None:
    circuit = build(d, rounds, NullInjector())
    # Must compile a detector sampler without error.
    circuit.compile_detector_sampler()


@pytest.mark.parametrize("d,rounds", CASES)
def test_all_detectors_deterministically_zero(d: int, rounds: int) -> None:
    """The gate: 1000 shots of the noiseless circuit fire no detectors."""
    circuit = build(d, rounds, NullInjector())
    sampler = circuit.compile_detector_sampler()
    dets, obs = sampler.sample(1000, separate_observables=True)
    assert not dets.any(), "noiseless circuit produced a detection event"
    assert not obs.any(), "noiseless circuit produced a logical observable flip"


@pytest.mark.parametrize("d,rounds", CASES)
def test_detector_count_closed_form(d: int, rounds: int) -> None:
    circuit = build(d, rounds, NullInjector())
    n_anc = d**2 - 1
    expected = n_anc // 2 + (rounds - 1) * n_anc + n_anc // 2
    assert circuit.num_detectors == expected
    assert num_detectors(d, rounds) == expected


@pytest.mark.parametrize("d,rounds", CASES)
def test_structural_diff_against_generated(d: int, rounds: int) -> None:
    """Detector count, observable count, and SGE length match Stim's generator."""
    mine = build(d, rounds, NullInjector())
    reference = stim.Circuit.generated(
        "surface_code:rotated_memory_z", distance=d, rounds=rounds
    )
    assert mine.num_detectors == reference.num_detectors
    assert mine.num_observables == reference.num_observables == 1


@pytest.mark.parametrize("d", DISTANCES)
def test_shortest_graphlike_error_matches_generated(
    d: int, make_noise: NoiseFactory
) -> None:
    """With matching uniform noise, both circuits have graphlike distance d."""
    p = 1e-3
    mine = build(d, d, make_noise(p))
    reference = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        distance=d,
        rounds=d,
        after_clifford_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p,
    )
    mine_len = len(mine.shortest_graphlike_error())
    ref_len = len(reference.shortest_graphlike_error())
    assert mine_len == ref_len == d


def test_num_detectors_requires_positive_rounds() -> None:
    with pytest.raises(ValueError):
        num_detectors(3, 0)
    with pytest.raises(ValueError):
        build(3, 0, NullInjector())
