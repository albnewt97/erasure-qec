"""Wu et al. (2022) two-qubit-gate erasure noise: sampling circuit only.

Reference: Wu, Kolkowitz, Puri, Thompson, Nat. Commun. 13, 4657 (2022),
arXiv:2201.03540, Sec. III and Methods VII.A. Per two-qubit gate, exactly one
of the following happens:

- with probability ``p_p = p (1 - R_e)``: a Pauli uniform over the 15
  non-identity elements of ``{I,X,Y,Z}^2``;
- with probability ``p_e = p R_e``: an erasure of *both* qubits, i.e. one
  heralded event after which each qubit is in ``I/2``, sampled as a Pauli
  uniform over all 16 elements of ``{I,X,Y,Z}^2``;
- otherwise nothing.

Ancilla initialization and measurement each flip with probability ``p_m``
(``X_ERROR``). Idle errors and single-qubit gate errors are neglected.

The sampler is exclusive: Pauli and erasure branches share one
``CORRELATED_ERROR`` / ``ELSE_CORRELATED_ERROR`` chain per gate, so they never
co-occur on one gate. The decoding circuit (P4) treats Pauli and erasure as
independent channels; that is a decoder approximation at ``O(p_p * p_e)``.

The herald is a flag qubit (index ``num_code_qubits(d)`` = ``2 d^2 - 1``, no
``QUBIT_COORDS``) that is reset, flipped by the erasure branches, and measured
once per gate pair, giving one herald record bit per pair. Because the herald
and the Pauli are one correlated mechanism, stim's ``decompose_errors``
either merges the herald with a syndrome detector (``R_e = 0.98``) or fails
(``R_e = 1.0``). Never call ``detector_error_model()`` on circuits built with
:class:`WuSamplingInjector`; decode with a separate decoding circuit.

Assumptions not fixed by Wu et al.:

- Data-qubit initialization is noiseless (Wu does not specify it), and the
  final data measurement is perfect (Wu: d rounds then one round of perfect
  measurements). ``p_m`` applies to ancillas only.
- When ``p_p == 0`` the 15 Pauli branches are omitted from the chain; they
  have probability 0, so the sampled distribution is unchanged.
- At ``p_e == 0`` the gate noise is a single ``DEPOLARIZE2(p_p)`` over the
  whole CX layer, which makes the circuit identical to ``ErasureInjector`` at
  ``r_e = 0`` with no measurement, reset or idle noise.
"""

import itertools
from collections.abc import Sequence
from dataclasses import dataclass

import stim

from erasure_qec.circuits.layout import num_code_qubits
from erasure_qec.config import _check_probability

_PAULIS = "IXYZ"
_PAULI_TARGET = {"X": stim.target_x, "Y": stim.target_y, "Z": stim.target_z}


@dataclass(frozen=True)
class WuNoiseParams:
    """Wu et al. noise parameters.

    ``p`` is the two-qubit-gate error probability, ``R_e`` the erasure
    fraction, ``p_m`` the ancilla init/measurement flip probability. With
    ``p_m_follows_p`` the effective ``p_m`` is ``p``.
    """

    p: float
    R_e: float
    p_m: float = 0.0
    p_m_follows_p: bool = False

    def __post_init__(self) -> None:
        _check_probability(self.p, "p")
        _check_probability(self.R_e, "R_e")
        _check_probability(self.p_m, "p_m")
        if self.p_e + self.p_p > 1.0:
            raise ValueError(f"p_e + p_p must be <= 1, got {self.p_e + self.p_p}")

    @property
    def p_e(self) -> float:
        """Per-gate pair-erasure probability ``p * R_e``."""
        return self.p * self.R_e

    @property
    def p_p(self) -> float:
        """Per-gate non-identity Pauli probability ``p * (1 - R_e)``."""
        return self.p * (1.0 - self.R_e)

    @property
    def meas(self) -> float:
        """Effective ancilla init/measurement flip probability."""
        return self.p if self.p_m_follows_p else self.p_m


def append_wu_gate_chain(
    circuit: stim.Circuit, a: int, b: int, flag: int, p_e: float, p_p: float
) -> None:
    """Append ``R flag``, one exclusive error chain on ``(a, b)``, then ``M flag``.

    The chain has 16 erasure branches (probability ``p_e / 16`` each, targets
    ``X(flag)`` plus the non-identity factors of each ``(P_a, P_b)`` in
    ``itertools.product("IXYZ", "IXYZ")`` order), followed by 15 non-identity
    Pauli branches without the flag (``p_p / 15`` each; omitted if
    ``p_p == 0``). The first branch is ``CORRELATED_ERROR(q)``, each later one
    ``ELSE_CORRELATED_ERROR(q / (1 - cumulative))``, so branch marginals are
    exactly ``q``. Appends exactly one measurement record bit (the herald).
    Requires ``p_e > 0`` and ``p_e + p_p <= 1``.
    """
    branches: list[tuple[float, list[stim.GateTarget]]] = []
    for pa, pb in itertools.product(_PAULIS, _PAULIS):
        branches.append((p_e / 16, [stim.target_x(flag), *_pauli_targets(pa, a, pb, b)]))
    if p_p > 0.0:
        for pa, pb in itertools.product(_PAULIS, _PAULIS):
            if pa == "I" and pb == "I":
                continue
            branches.append((p_p / 15, _pauli_targets(pa, a, pb, b)))

    circuit.append("R", [flag])
    used = 0.0
    for k, (q, targets) in enumerate(branches):
        # min() only absorbs float round-off when p_e + p_p == 1.
        conditional = min(1.0, q / (1.0 - used))
        circuit.append(
            "CORRELATED_ERROR" if k == 0 else "ELSE_CORRELATED_ERROR", targets, conditional
        )
        used += q
    circuit.append("M", [flag])


def _pauli_targets(pa: str, a: int, pb: str, b: int) -> list[stim.GateTarget]:
    """Pauli targets for ``P_a(a) P_b(b)``, skipping identity factors."""
    targets: list[stim.GateTarget] = []
    if pa != "I":
        targets.append(_PAULI_TARGET[pa](a))
    if pb != "I":
        targets.append(_PAULI_TARGET[pb](b))
    return targets


class WuSamplingInjector:
    """``NoiseInjector`` emitting the Wu et al. model as a sampling circuit.

    Sampling only: do not build a detector error model from its circuits (see
    module docstring). ``locations`` records ``(a, b, k)`` for every gate pair
    that received an erasure chain, where ``k`` counts pairs in emission order
    (equal to the pair's herald-bit index among all herald bits).
    """

    def __init__(self, d: int, params: WuNoiseParams) -> None:
        self.params = params
        self.flag = num_code_qubits(d)
        self.locations: list[tuple[int, int, int]] = []

    def after_reset(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        if qubits and self.params.meas > 0.0:
            circuit.append("X_ERROR", list(qubits), self.params.meas)

    def before_measure(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        if qubits and self.params.meas > 0.0:
            circuit.append("X_ERROR", list(qubits), self.params.meas)

    def after_initial_data_reset(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        return None

    def before_final_data_measure(self, circuit: stim.Circuit, qubits: Sequence[int]) -> None:
        return None

    def on_two_qubit_gate(
        self, circuit: stim.Circuit, pairs: Sequence[tuple[int, int]]
    ) -> list[int]:
        if not pairs:
            return []
        p_e, p_p = self.params.p_e, self.params.p_p
        if p_e == 0.0:
            if p_p > 0.0:
                circuit.append("DEPOLARIZE2", [q for pair in pairs for q in pair], p_p)
            return []
        for a, b in pairs:
            append_wu_gate_chain(circuit, a, b, self.flag, p_e, p_p)
            self.locations.append((a, b, len(self.locations)))
        return [a for a, _ in pairs]

    def on_idle(self, circuit: stim.Circuit, qubits: Sequence[int]) -> list[int]:
        return []
