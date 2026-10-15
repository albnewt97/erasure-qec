"""Text-template construction of the Wu circuits equals per-instruction appends.

The reference injectors below keep the original construction (one
``circuit.append`` per instruction, as committed in P3/P4) inside the test.
The production injectors emit stim program text parsed once per CX layer.
Both ``str()`` (what the golden digests hash; stim prints ~6 significant
digits) and exact ``==`` (full double precision) must agree.
"""

import itertools
from collections.abc import Sequence

import pytest
import stim

from erasure_qec.circuits.builder import build
from erasure_qec.noise.wu import WuDecodingInjector, WuNoiseParams, WuSamplingInjector

_PAULI_TARGET = {"X": stim.target_x, "Y": stim.target_y, "Z": stim.target_z}


def _reference_pauli_targets(pa: str, a: int, pb: str, b: int) -> list[stim.GateTarget]:
    targets: list[stim.GateTarget] = []
    if pa != "I":
        targets.append(_PAULI_TARGET[pa](a))
    if pb != "I":
        targets.append(_PAULI_TARGET[pb](b))
    return targets


def _reference_chain(
    circuit: stim.Circuit, a: int, b: int, flag: int, p_e: float, p_p: float
) -> None:
    branches: list[tuple[float, list[stim.GateTarget]]] = []
    if p_e > 0.0:
        for pa, pb in itertools.product("IXYZ", "IXYZ"):
            branches.append(
                (p_e / 16, [stim.target_x(flag), *_reference_pauli_targets(pa, a, pb, b)])
            )
    if p_p > 0.0:
        for pa, pb in itertools.product("IXYZ", "IXYZ"):
            if pa == "I" and pb == "I":
                continue
            branches.append((p_p / 15, _reference_pauli_targets(pa, a, pb, b)))

    circuit.append("R", [flag])
    used = 0.0
    for k, (q, targets) in enumerate(branches):
        conditional = min(1.0, q / (1.0 - used))
        circuit.append(
            "CORRELATED_ERROR" if k == 0 else "ELSE_CORRELATED_ERROR", targets, conditional
        )
        used += q
    circuit.append("M", [flag])


class _ReferenceSampling(WuSamplingInjector):
    def _emit_heralded_layer(self, circuit: stim.Circuit, pairs: Sequence[tuple[int, int]]) -> None:
        for _, a, b, p_e_k in self._enumerate_locations(pairs):
            _reference_chain(circuit, a, b, self.flag, p_e_k, self.params.p_p)


class _ReferenceDecoding(WuDecodingInjector):
    def _emit_heralded_layer(self, circuit: stim.Circuit, pairs: Sequence[tuple[int, int]]) -> None:
        if self.params.p_p > 0.0:
            circuit.append("DEPOLARIZE2", [q for pair in pairs for q in pair], self.params.p_p)
        for k, a, b, p_e_k in self._enumerate_locations(pairs):
            if p_e_k > 0.0:
                for gate, qubit, side in (
                    ("X_ERROR", a, "Xa"),
                    ("Z_ERROR", a, "Za"),
                    ("X_ERROR", b, "Xb"),
                    ("Z_ERROR", b, "Zb"),
                ):
                    circuit.append(
                        stim.CircuitInstruction(gate, [qubit], [p_e_k / 2], tag=f"wu_e{k}:{side}")
                    )
            circuit.append("R", [self.flag])
            circuit.append("M", [self.flag])


P = 0.03
OVERRIDES = {0: 0.9, 4: 0.0, 9: 0.3}


@pytest.mark.parametrize("use_overrides", [False, True], ids=["no_overrides", "overrides"])
@pytest.mark.parametrize("p_m", [0.0, 0.005])
@pytest.mark.parametrize("r_e", [0.5, 0.98, 1.0])
@pytest.mark.parametrize("d", [3, 5])
def test_text_build_equals_append_build(
    d: int, r_e: float, p_m: float, use_overrides: bool
) -> None:
    params = WuNoiseParams(P, r_e, p_m=p_m)
    overrides = OVERRIDES if use_overrides else {}
    for new_cls, ref_cls in (
        (WuSamplingInjector, _ReferenceSampling),
        (WuDecodingInjector, _ReferenceDecoding),
    ):
        new = build(d, d, new_cls(d, params, location_overrides=overrides))
        ref = build(d, d, ref_cls(d, params, location_overrides=overrides))
        assert str(new) == str(ref), new_cls.__name__
        assert new == ref, new_cls.__name__


def test_text_build_parses_scientific_notation_probabilities() -> None:
    """Small p gives repr() floats like '6.25e-06'; parsing must stay exact."""
    params = WuNoiseParams(1e-4, 0.5)
    for new_cls, ref_cls in (
        (WuSamplingInjector, _ReferenceSampling),
        (WuDecodingInjector, _ReferenceDecoding),
    ):
        new = build(3, 3, new_cls(3, params))
        ref = build(3, 3, ref_cls(3, params))
        assert "e-" in str(new)
        assert new == ref
