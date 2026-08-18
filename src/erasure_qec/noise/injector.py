"""Noise injection interface (PLAN.md §3, §5).

The builder emits the bare Clifford + measurement structure and calls a
``NoiseInjector`` at each physically-motivated fault location. This keeps
``builder.py`` free of noise-model details: swapping the injector swaps the
noise model without touching circuit construction.

``on_two_qubit_gate`` returns the qubit indices that received a
``HERALDED_ERASE`` channel, **in the exact order their herald bits were
appended to the measurement record** — the builder uses this to immediately
emit the matching sentinel ``DETECTOR(x, y, t, 1)`` per herald bit (§3.4).
Injectors that never herald (``NullInjector``, ``PauliOnlyInjector``) always
return an empty list.
"""

from collections.abc import Sequence
from typing import Protocol

import stim

from erasure_qec.config import NoiseParams
from erasure_qec.noise.model import channel_rates


class NoiseInjector(Protocol):
    """Hook set the builder calls at each fault location in a round.

    Every method receives the in-progress ``stim.Circuit`` and appends its
    noise instructions in place. Qubits are given as integer indices matching
    the circuit's ``QUBIT_COORDS`` assignment.
    """

    def after_reset(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        """Append reset noise (e.g. ``X_ERROR``) right after an ``R`` on ``qubits``."""
        ...

    def before_measure(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        """Append measurement noise (e.g. ``X_ERROR``) right before an ``M``/``MR``."""
        ...

    def on_two_qubit_gate(
        self, circuit: stim.Circuit, pairs: Sequence[tuple[int, int]]
    ) -> list[int]:
        """Append gate noise for the ``(control, target)`` pairs of one CX layer.

        Returns the qubit indices that received a ``HERALDED_ERASE`` channel,
        in herald-record order (empty if none).
        """
        ...

    def on_idle(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        """Append idle noise (e.g. ``DEPOLARIZE1``) for qubits idle during a TICK."""
        ...


class NullInjector:
    """The noiseless identity injector: every hook is a no-op."""

    def after_reset(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        return None

    def before_measure(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        return None

    def on_two_qubit_gate(
        self, circuit: stim.Circuit, pairs: Sequence[tuple[int, int]]
    ) -> list[int]:
        return []

    def on_idle(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        return None


class ErasureInjector:
    """Full erasure noise model (PLAN.md §5).

    NB: ``HERALDED_ERASE`` replaces the qubit with the maximally-mixed ``I/2`` --
    an *unbiased* erasure (all four Paulis equally likely, conditioned on the
    herald). This does not model the Z-biased leakage of Sahay et al.; the
    threshold advantage here comes from herald-conditioned decoding, not bias.

    Per two-qubit gate: ``DEPOLARIZE2(p * (1 - r_e))`` residual Pauli noise,
    plus ``HERALDED_ERASE`` independently on each qubit at the rate set by
    :func:`~erasure_qec.noise.model.channel_rates` (``(2/3) * p * r_e``, chosen
    so the per-gate non-identity error budget is ``p`` for every ``r_e``).
    ``X_ERROR`` around every reset/measurement, ``DEPOLARIZE1`` on idle qubits.
    Every rate that resolves to exactly 0 is skipped entirely (no
    zero-probability instructions are ever appended), so ``p=0`` reproduces
    the noiseless circuit byte-for-byte.
    """

    def __init__(self, params: NoiseParams) -> None:
        self._rates = channel_rates(params)

    def after_reset(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        if qubits and self._rates.reset_flip > 0.0:
            circuit.append("X_ERROR", list(qubits), self._rates.reset_flip)

    def before_measure(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        if qubits and self._rates.meas_flip > 0.0:
            circuit.append("X_ERROR", list(qubits), self._rates.meas_flip)

    def on_two_qubit_gate(
        self, circuit: stim.Circuit, pairs: Sequence[tuple[int, int]]
    ) -> list[int]:
        if not pairs:
            return []
        flattened = [q for pair in pairs for q in pair]
        if self._rates.depolarize2 > 0.0:
            circuit.append("DEPOLARIZE2", flattened, self._rates.depolarize2)
        if self._rates.herald <= 0.0:
            return []
        circuit.append("HERALDED_ERASE", flattened, self._rates.herald)
        return flattened

    def on_idle(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        if qubits and self._rates.idle_depolarize > 0.0:
            circuit.append("DEPOLARIZE1", list(qubits), self._rates.idle_depolarize)


class PauliOnlyInjector:
    """Residual-Pauli-only noise: the ``r_e = 0`` special case of §5.

    Equivalent to ``ErasureInjector`` with the erasure fraction forced
    to zero: ``DEPOLARIZE2(p)`` per two-qubit gate, no ``HERALDED_ERASE`` and
    no herald detectors at all. Useful as the non-erasure baseline (e.g. the
    ``baseline_pauli`` experiment config) and for the §4 distance invariants.
    """

    def __init__(self, params: NoiseParams) -> None:
        self._delegate = ErasureInjector(
            NoiseParams(
                p=params.p,
                r_e=0.0,
                p_meas=params.p_meas,
                p_reset=params.p_reset,
                p_idle=params.p_idle,
            )
        )

    def after_reset(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        self._delegate.after_reset(circuit, qubits)

    def before_measure(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        self._delegate.before_measure(circuit, qubits)

    def on_two_qubit_gate(
        self, circuit: stim.Circuit, pairs: Sequence[tuple[int, int]]
    ) -> list[int]:
        return self._delegate.on_two_qubit_gate(circuit, pairs)

    def on_idle(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        self._delegate.on_idle(circuit, qubits)
