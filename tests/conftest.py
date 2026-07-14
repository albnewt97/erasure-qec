"""Shared pytest fixtures for the erasure_qec test suite."""

from collections.abc import Callable, Sequence

import pytest
import stim


class UniformCircuitNoise:
    """Standard circuit-level depolarizing noise implementing ``NoiseInjector``.

    Sprinkles ``DEPOLARIZE2`` after every CX layer and ``X_ERROR`` around
    resets/measurements, giving the noiseless builder graphlike error
    mechanisms so distance invariants can be probed (PLAN.md §4).
    """

    def __init__(self, p: float) -> None:
        self.p = p

    def after_reset(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        if qubits:
            circuit.append("X_ERROR", list(qubits), self.p)

    def before_measure(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        if qubits:
            circuit.append("X_ERROR", list(qubits), self.p)

    def on_two_qubit_gate(
        self, circuit: stim.Circuit, pairs: Sequence[tuple[int, int]]
    ) -> list[int]:
        for a, b in pairs:
            circuit.append("DEPOLARIZE2", [a, b], self.p)
        return []

    def on_idle(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        return None


@pytest.fixture
def make_noise() -> Callable[[float], UniformCircuitNoise]:
    """Factory fixture: ``make_noise(p)`` -> a ``UniformCircuitNoise`` injector."""
    return UniformCircuitNoise
