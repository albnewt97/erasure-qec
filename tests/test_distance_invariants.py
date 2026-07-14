"""Distance & hook invariants for the CX schedule (PLAN.md §4).

With a small uniform depolarizing sprinkle, the hook-safe schedule must have
graphlike distance exactly ``d`` and no undetectable logical error shorter
than ``d``. The deliberately broken schedule must degrade the distance,
proving the scheduling logic is load-bearing.
"""

from collections.abc import Callable

import pytest

from erasure_qec.circuits.builder import build
from erasure_qec.circuits.scheduling import BROKEN_SCHEDULE, HOOK_SAFE_SCHEDULE
from erasure_qec.noise.injector import NoiseInjector

DISTANCES = [3, 5, 7]
NOISE_P = 1e-3

# Observed shortest graphlike error under the broken schedule: a single
# X-ancilla hook fault spans two lattice steps of the vertical logical, so the
# distance halves to ceil((d+1)/2). Recorded here as a tight regression pin.
BROKEN_SGE_LENGTH = {3: 2, 5: 3, 7: 4}

NoiseFactory = Callable[[float], NoiseInjector]


@pytest.mark.parametrize("d", DISTANCES)
def test_correct_schedule_graphlike_distance_is_d(
    d: int, make_noise: NoiseFactory
) -> None:
    circuit = build(d, d, make_noise(NOISE_P), schedule=HOOK_SAFE_SCHEDULE)
    assert len(circuit.shortest_graphlike_error()) == d


# The exhaustive undetectable-error search is fast at d=3,5 but ~30s at d=7,
# so d=7 is gated behind the `slow` marker (d<7 still runs every CI pass).
_SULE_CASES = [3, 5, pytest.param(7, marks=pytest.mark.slow)]


@pytest.mark.parametrize("d", _SULE_CASES)
def test_correct_schedule_no_undetectable_error_shorter_than_d(
    d: int, make_noise: NoiseFactory
) -> None:
    circuit = build(d, d, make_noise(NOISE_P), schedule=HOOK_SAFE_SCHEDULE)
    error = circuit.search_for_undetectable_logical_errors(
        dont_explore_detection_event_sets_with_size_above=4,
        dont_explore_edges_with_degree_above=4,
        dont_explore_edges_increasing_symptom_degree=False,
    )
    assert len(error) == d


@pytest.mark.parametrize("d", DISTANCES)
def test_broken_schedule_reduces_distance(d: int, make_noise: NoiseFactory) -> None:
    circuit = build(d, d, make_noise(NOISE_P), schedule=BROKEN_SCHEDULE)
    observed = len(circuit.shortest_graphlike_error())
    assert observed < d, f"broken schedule failed to reduce distance (got {observed})"
    assert observed == BROKEN_SGE_LENGTH[d]


@pytest.mark.parametrize("d", DISTANCES)
def test_broken_schedule_is_strictly_worse_than_correct(
    d: int, make_noise: NoiseFactory
) -> None:
    correct = build(d, d, make_noise(NOISE_P), schedule=HOOK_SAFE_SCHEDULE)
    broken = build(d, d, make_noise(NOISE_P), schedule=BROKEN_SCHEDULE)
    assert len(broken.shortest_graphlike_error()) < len(
        correct.shortest_graphlike_error()
    )
