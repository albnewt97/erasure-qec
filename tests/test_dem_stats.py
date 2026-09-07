"""Circuit-wide heralded fraction of the DEM error budget (Phase 2.1).

These pin the honest "how much of the noise is erasure-converted" number, which
is NOT r_e: measurement/reset errors are never heralded, and idle errors only
when convert_idle is set, so the circuit-wide fraction sits below r_e.

They also pin the no-op-herald correction: HERALDED_ERASE heralds on its
identity branch, producing DEM mechanisms whose only target is a herald
detector. Those cause no syndrome and no logical error, so they belong in
neither side of the ratio.
"""

import pytest

from erasure_qec.analysis.dem_stats import dem_error_mass, heralded_fraction
from erasure_qec.config import NoiseParams


@pytest.mark.slow
def test_no_op_herald_branches_exist_and_are_excluded() -> None:
    """The identity branch of HERALDED_ERASE puts mechanisms in the DEM that
    flip a herald detector and nothing else; both sides of the ratio drop them."""
    params = NoiseParams(p=0.02, r_e=0.98)
    mass = dem_error_mass(params, d=5)

    # The no-op bucket is real, not an empty category.
    assert mass.no_op > 0.0
    assert mass.unheralded > 0.0 and mass.heralded > 0.0

    corrected = heralded_fraction(params, d=5)
    legacy = heralded_fraction(params, d=5, count_no_op_heralds=True)

    # The correction removes the no-op mass from BOTH numerator and denominator.
    assert corrected == pytest.approx(mass.heralded / (mass.unheralded + mass.heralded))
    assert legacy == pytest.approx(
        (mass.heralded + mass.no_op) / (mass.unheralded + mass.heralded + mass.no_op)
    )
    # Counting no-ops in both inflates the fraction, so the corrected value is
    # strictly lower.
    assert corrected < legacy


@pytest.mark.slow
def test_no_op_mass_vanishes_without_erasure_conversion() -> None:
    """The no-op bucket is caused by HERALDED_ERASE: with r_e=0 there are no
    erasure channels, so it is empty and the correction changes nothing."""
    params = NoiseParams(p=0.02, r_e=0.0)
    mass = dem_error_mass(params, d=5)
    assert mass.no_op == 0.0
    assert mass.heralded == 0.0
    assert heralded_fraction(params, d=5) == heralded_fraction(
        params, d=5, count_no_op_heralds=True
    )


@pytest.mark.slow
def test_heralded_fraction_below_r_e_gate_only() -> None:
    """convert_idle=False (the committed-sweep model): only the 2q-gate budget
    is heralded, so the circuit-wide fraction is well below r_e."""
    hf_50 = heralded_fraction(NoiseParams(p=0.02, r_e=0.5), d=5)
    hf_98 = heralded_fraction(NoiseParams(p=0.02, r_e=0.98), d=5)
    assert hf_50 == pytest.approx(0.232, abs=0.01)
    assert hf_98 == pytest.approx(0.454, abs=0.01)
    # The whole point: r_e (2q-gate fraction) is NOT the circuit-wide fraction.
    assert hf_50 < 0.5
    assert hf_98 < 0.98


@pytest.mark.slow
def test_legacy_no_op_counting_reproduces_audit_numbers() -> None:
    """count_no_op_heralds=True keeps the pre-correction docs/AUDIT.md values
    reproducible."""
    assert heralded_fraction(
        NoiseParams(p=0.02, r_e=0.5), d=5, count_no_op_heralds=True
    ) == pytest.approx(0.305, abs=0.01)
    assert heralded_fraction(
        NoiseParams(p=0.02, r_e=0.98), d=5, count_no_op_heralds=True
    ) == pytest.approx(0.547, abs=0.01)


@pytest.mark.slow
def test_convert_idle_raises_heralded_fraction_but_stays_below_r_e() -> None:
    for r_e in (0.5, 0.98):
        off = heralded_fraction(NoiseParams(p=0.02, r_e=r_e), d=5)
        on = heralded_fraction(NoiseParams(p=0.02, r_e=r_e, convert_idle=True), d=5)
        assert on > off  # idle conversion adds heralded mass
        assert on < r_e  # still below r_e: meas/reset are never heralded
    assert heralded_fraction(
        NoiseParams(p=0.02, r_e=0.98, convert_idle=True), d=5
    ) == pytest.approx(0.631, abs=0.01)
