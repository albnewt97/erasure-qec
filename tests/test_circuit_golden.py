"""Golden SHA-256 digests of builder output for the existing noise model.

These pin the circuits the committed ``data/*.csv`` sweeps were collected
from. Each case builds ``build(d, d, injector)`` and hashes ``str(circuit)``.
A failure means the existing noise model or the builder output changed, which
invalidates the committed data: revert the change, do not update the digests.

The digests are tied to the stim version locked in ``uv.lock`` (stim 1.16.0
when recorded), because ``str(stim.Circuit)`` formatting is stim's. A stim
upgrade can change the text without changing the circuit; that case needs a
separate equivalence check, not a silent digest refresh.
"""

import hashlib
from collections.abc import Callable

import pytest
import stim

from erasure_qec.circuits.builder import build
from erasure_qec.config import NoiseParams
from erasure_qec.noise.injector import ErasureInjector, PauliOnlyInjector


def _erasure(
    d: int, *, p: float, r_e: float, convert_idle: bool = False
) -> Callable[[], stim.Circuit]:
    params = NoiseParams(p=p, r_e=r_e, convert_idle=convert_idle)
    return lambda: build(d, d, ErasureInjector(params))


GOLDEN: list[tuple[str, Callable[[], stim.Circuit], str]] = [
    (
        "erasure_d3_p0.01_re0.0",
        _erasure(3, p=0.01, r_e=0.0),
        "e5e7b23401327112a7046f1bc97afda74093dd3ea347e6c3572edf161371fd34",
    ),
    (
        "erasure_d5_p0.02_re0.5",
        _erasure(5, p=0.02, r_e=0.5),
        "88915a2067815956d27150d03711c2b7e4255e209e7e7225fe26f530bb7c47ed",
    ),
    (
        "erasure_d5_p0.03_re0.98",
        _erasure(5, p=0.03, r_e=0.98),
        "c25ce4ada5be851a17dcb47d2cecdd338fa5fd84b66f73bde3244c102ca08422",
    ),
    (
        "erasure_d3_p0.01_re0.5_convert_idle",
        _erasure(3, p=0.01, r_e=0.5, convert_idle=True),
        "14cbc40fac4927df0a2eba5f896a58447de599cf15856fe7f846787064ca5e46",
    ),
    (
        "pauli_only_d3_p0.01_re0.0",
        lambda: build(3, 3, PauliOnlyInjector(NoiseParams(p=0.01, r_e=0.0))),
        "e5e7b23401327112a7046f1bc97afda74093dd3ea347e6c3572edf161371fd34",
    ),
    (
        "erasure_d3_p0.01_re0.0_probes_1+3j_1+5j",
        lambda: build(
            3,
            3,
            ErasureInjector(NoiseParams(p=0.01, r_e=0.0)),
            probe_erasures=[(complex(1, 3), 0), (complex(1, 5), 0)],
        ),
        "22f267414f0f134d2d5cc0c691ce1a395451c55d64a219ba2928cd6598c35b1c",
    ),
]


@pytest.mark.parametrize(
    ("make_circuit", "expected"),
    [(make, digest) for _, make, digest in GOLDEN],
    ids=[name for name, _, _ in GOLDEN],
)
def test_circuit_digest_matches_golden(
    make_circuit: Callable[[], stim.Circuit], expected: str
) -> None:
    """``sha256(str(circuit))`` equals the digest recorded for the committed data."""
    digest = hashlib.sha256(str(make_circuit()).encode()).hexdigest()
    assert digest == expected
