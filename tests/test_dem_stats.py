"""Circuit-wide heralded fraction of the DEM error budget (Phase 2.1).

These pin the honest "how much of the noise is erasure-converted" number, which
is NOT r_e: measurement/reset errors are never heralded, and idle errors only
when convert_idle is set, so the circuit-wide fraction sits below r_e.
"""

import pytest

from erasure_qec.analysis.dem_stats import heralded_fraction
from erasure_qec.config import NoiseParams


@pytest.mark.slow
def test_heralded_fraction_below_r_e_gate_only() -> None:
    """convert_idle=False (the committed-sweep model): only the 2q-gate budget
    is heralded, so the circuit-wide fraction is well below r_e."""
    hf_50 = heralded_fraction(NoiseParams(p=0.02, r_e=0.5), d=5)
    hf_98 = heralded_fraction(NoiseParams(p=0.02, r_e=0.98), d=5)
    assert hf_50 == pytest.approx(0.305, abs=0.01)
    assert hf_98 == pytest.approx(0.547, abs=0.01)
    # The whole point: r_e (2q-gate fraction) is NOT the circuit-wide fraction.
    assert hf_50 < 0.5
    assert hf_98 < 0.98


@pytest.mark.slow
def test_convert_idle_raises_heralded_fraction_but_stays_below_r_e() -> None:
    for r_e in (0.5, 0.98):
        off = heralded_fraction(NoiseParams(p=0.02, r_e=r_e), d=5)
        on = heralded_fraction(NoiseParams(p=0.02, r_e=r_e, convert_idle=True), d=5)
        assert on > off  # idle conversion adds heralded mass
        assert on < r_e  # still below r_e: meas/reset are never heralded
    assert heralded_fraction(
        NoiseParams(p=0.02, r_e=0.98, convert_idle=True), d=5
    ) == pytest.approx(0.720, abs=0.01)
