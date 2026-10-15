"""Wu et al. (2022) two-qubit-gate erasure noise: sampling and decoding circuits.

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

Sampling circuit (:class:`WuSamplingInjector`)
----------------------------------------------
The sampler is exclusive: Pauli and erasure branches share one
``CORRELATED_ERROR`` / ``ELSE_CORRELATED_ERROR`` chain per gate, so they never
co-occur on one gate. The herald is a flag qubit (index ``num_code_qubits(d)``
= ``2 d^2 - 1``, no ``QUBIT_COORDS``) that is reset, flipped by the erasure
branches, and measured once per gate pair, giving one herald record bit per
pair. Because the herald and the Pauli are one correlated mechanism, stim's
``decompose_errors`` either merges the herald with a syndrome detector
(``R_e = 0.98``) or fails (``R_e = 1.0``). Never call
``detector_error_model()`` on circuits built with :class:`WuSamplingInjector`.

Decoding circuit (:class:`WuDecodingInjector`)
----------------------------------------------
Same builder call and same detector layout, but at each erasure location ``k``
the erasure is emitted as four independent tagged generators
``X_ERROR[wu_e{k}:Xa]``, ``Z_ERROR[wu_e{k}:Za]``, ``X_ERROR[wu_e{k}:Xb]``,
``Z_ERROR[wu_e{k}:Zb]``, each at ``p_e / 2``, followed by a noiseless
``R flag; M flag``. A Pauli uniform over ``{I,X,Y,Z}^2`` is four independent
fair bits ``(x_a, z_a, x_b, z_b)``: conditional on the herald each generator
fires with probability 1/2 (matching weight 0); unconditionally each fires
with probability ``p_e / 2`` (the blind decoder's marginal). The tags let the
decoder attribute each DEM mechanism to its herald without the herald
detector appearing in the mechanism. The decoding circuit treats Pauli
(``DEPOLARIZE2(p_p)``) and erasure as independent channels, while the sampler
makes them exclusive; that is a decoder approximation at ``O(p_p * p_e)``.

Flag layout rule (both injectors): whenever ``params.p_e > 0`` or
``location_overrides`` is non-empty, every gate pair gets ``R flag ... M flag``
and one herald detector, even if that pair's effective ``p_e`` is 0. The
detector layout therefore does not depend on which locations are active.
Location ``k`` is the global pair index in emission order, equal to the
herald index; both injectors assign it through one shared enumerator.

Both injectors emit their per-pair instructions as stim program text (one
cached template per probability, formatted per pair and parsed once per CX
layer) instead of one ``append`` call per instruction. Probabilities are
written with ``repr()``, so parsing reproduces the exact doubles;
tests/test_wu_text_build.py checks the result against the per-instruction
``append`` construction.

Assumptions not fixed by Wu et al.:

- Data-qubit initialization is noiseless (Wu does not specify it), and the
  final data measurement is perfect (Wu: d rounds then one round of perfect
  measurements). ``p_m`` applies to ancillas only.
- Zero-probability branches are omitted: the 15 Pauli branches when
  ``p_p == 0`` and the 16 erasure branches when a location's ``p_e`` is 0.
  The sampled distribution is unchanged.
- At ``p_e == 0`` with no overrides the gate noise is a single
  ``DEPOLARIZE2(p_p)`` over the whole CX layer, which makes the circuit
  identical to ``ErasureInjector`` at ``r_e = 0`` with no measurement, reset
  or idle noise.
"""

import functools
import itertools
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType

import stim

from erasure_qec.circuits.layout import num_code_qubits
from erasure_qec.config import _check_probability

_PAULIS = "IXYZ"
_NO_OVERRIDES: Mapping[int, float] = MappingProxyType({})


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
    ``itertools.product("IXYZ", "IXYZ")`` order; omitted if ``p_e == 0``),
    followed by 15 non-identity Pauli branches without the flag (``p_p / 15``
    each; omitted if ``p_p == 0``). The first branch is
    ``CORRELATED_ERROR(q)``, each later one
    ``ELSE_CORRELATED_ERROR(q / (1 - cumulative))``, so branch marginals are
    exactly ``q``. Appends exactly one measurement record bit (the herald).
    Requires ``p_e + p_p <= 1``.
    """
    circuit.append_from_stim_program_text(_chain_template(p_e, p_p).format(a=a, b=b, f=flag))


@functools.lru_cache(maxsize=256)
def _chain_template(p_e: float, p_p: float) -> str:
    """Stim text of :func:`append_wu_gate_chain` with ``{a}``, ``{b}``, ``{f}`` placeholders."""
    branches: list[tuple[float, str]] = []
    if p_e > 0.0:
        for pa, pb in itertools.product(_PAULIS, _PAULIS):
            branches.append((p_e / 16, " ".join(["X{f}", *_pauli_text(pa, pb)])))
    if p_p > 0.0:
        for pa, pb in itertools.product(_PAULIS, _PAULIS):
            if pa == "I" and pb == "I":
                continue
            branches.append((p_p / 15, " ".join(_pauli_text(pa, pb))))

    lines = ["R {f}"]
    used = 0.0
    for k, (q, targets) in enumerate(branches):
        # min() only absorbs float round-off when p_e + p_p == 1.
        conditional = min(1.0, q / (1.0 - used))
        name = "CORRELATED_ERROR" if k == 0 else "ELSE_CORRELATED_ERROR"
        lines.append(f"{name}({conditional!r}) {targets}")
        used += q
    lines.append("M {f}")
    return "\n".join(lines) + "\n"


def _pauli_text(pa: str, pb: str) -> list[str]:
    """Stim text targets for ``P_a(a) P_b(b)``, skipping identity factors."""
    targets: list[str] = []
    if pa != "I":
        targets.append(pa + "{a}")
    if pb != "I":
        targets.append(pb + "{b}")
    return targets


@functools.lru_cache(maxsize=256)
def _decoding_template(p_e: float) -> str:
    """Stim text of one decoding location: tagged generators (if ``p_e > 0``) and the flag.

    Placeholders ``{k}`` (location), ``{a}``, ``{b}``, ``{f}`` (flag).
    """
    lines: list[str] = []
    if p_e > 0.0:
        prob = repr(p_e / 2)
        for gate, qubit, side in (
            ("X_ERROR", "{a}", "Xa"),
            ("Z_ERROR", "{a}", "Za"),
            ("X_ERROR", "{b}", "Xb"),
            ("Z_ERROR", "{b}", "Zb"),
        ):
            lines.append(f"{gate}[wu_e{{k}}:{side}]({prob}) {qubit}")
    lines += ["R {f}", "M {f}"]
    return "\n".join(lines) + "\n"


class _WuInjectorBase:
    """Shared hooks, flag layout and location enumeration of the Wu injectors.

    ``location_overrides`` maps a location index to that pair's ``p_e``. It is
    a test affordance (forcing or isolating single erasure locations); the
    default is empty. ``locations`` records ``(a, b, k)`` for every pair that
    received a flag, ``k`` being the global pair index in emission order.
    """

    def __init__(
        self,
        d: int,
        params: WuNoiseParams,
        location_overrides: Mapping[int, float] = _NO_OVERRIDES,
    ) -> None:
        for k, p_e_k in location_overrides.items():
            _check_probability(p_e_k, f"location_overrides[{k}]")
            if p_e_k + params.p_p > 1.0:
                raise ValueError(f"location {k}: p_e + p_p must be <= 1")
        self.params = params
        self.location_overrides: dict[int, float] = dict(location_overrides)
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

    def on_idle(self, circuit: stim.Circuit, qubits: Sequence[int]) -> list[int]:
        return []

    def on_two_qubit_gate(
        self, circuit: stim.Circuit, pairs: Sequence[tuple[int, int]]
    ) -> list[int]:
        if not pairs:
            return []
        if self.params.p_e == 0.0 and not self.location_overrides:
            if self.params.p_p > 0.0:
                circuit.append("DEPOLARIZE2", [q for pair in pairs for q in pair], self.params.p_p)
            return []
        self._emit_heralded_layer(circuit, pairs)
        return [a for a, _ in pairs]

    def _enumerate_locations(
        self, pairs: Sequence[tuple[int, int]]
    ) -> Iterator[tuple[int, int, int, float]]:
        """Assign location indices: yield ``(k, a, b, p_e_k)`` per pair, in order.

        The only place ``k`` is assigned, so sampling and decoding circuits
        enumerate locations identically.
        """
        for a, b in pairs:
            k = len(self.locations)
            self.locations.append((a, b, k))
            yield k, a, b, self.location_overrides.get(k, self.params.p_e)

    def _emit_heralded_layer(self, circuit: stim.Circuit, pairs: Sequence[tuple[int, int]]) -> None:
        raise NotImplementedError


class WuSamplingInjector(_WuInjectorBase):
    """``NoiseInjector`` emitting the Wu et al. model as a sampling circuit.

    Sampling only: do not build a detector error model from its circuits (see
    module docstring). Each flagged pair gets :func:`append_wu_gate_chain`
    with that location's ``p_e`` and the global ``p_p``.
    """

    def _emit_heralded_layer(self, circuit: stim.Circuit, pairs: Sequence[tuple[int, int]]) -> None:
        p_p, f = self.params.p_p, self.flag
        circuit.append_from_stim_program_text(
            "".join(
                _chain_template(p_e_k, p_p).format(a=a, b=b, f=f)
                for _, a, b, p_e_k in self._enumerate_locations(pairs)
            )
        )


class WuDecodingInjector(_WuInjectorBase):
    """``NoiseInjector`` emitting the Wu et al. decoding circuit.

    Same detector layout as :class:`WuSamplingInjector` with equal arguments.
    Per CX layer: ``DEPOLARIZE2(p_p)`` on all pairs, then per pair four tagged
    generators at ``p_e_k / 2`` and a noiseless ``R flag; M flag``.
    """

    def _emit_heralded_layer(self, circuit: stim.Circuit, pairs: Sequence[tuple[int, int]]) -> None:
        if self.params.p_p > 0.0:
            circuit.append("DEPOLARIZE2", [q for pair in pairs for q in pair], self.params.p_p)
        f = self.flag
        circuit.append_from_stim_program_text(
            "".join(
                _decoding_template(p_e_k).format(k=k, a=a, b=b, f=f)
                for k, a, b, p_e_k in self._enumerate_locations(pairs)
            )
        )
