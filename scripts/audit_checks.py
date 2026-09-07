"""Re-runnable reproduction of the external audit (Phase 0).

Emits one line per audit check with the measured value next to the expected
value, so a reviewer can confirm every finding. Nothing here changes state; it
only measures. Run:

    uv run python scripts/audit_checks.py

NOTE: the "EXPECTED" column is the *pre-fix* audit baseline. The Phase-1 fixes
(commits 6b6b87d noise, 0ce6ff5 fit) have since landed, so on HEAD several rows
now differ from EXPECTED **by design** — that difference is the evidence the
fixes took effect. See docs/AUDIT.md "Resolution (Phase 1)" for the post-fix
values. To reproduce the original findings, run this against the pre-fix
commit.
"""

from __future__ import annotations

import numpy as np
import stim

from erasure_qec.analysis.dem_stats import heralded_fraction
from erasure_qec.analysis.statistics import (
    Z_95,
    SweepPoint,
    load_sweep,
    per_round_p_l,
)
from erasure_qec.analysis.threshold_fit import (
    _model,
    _select_window,
    estimate_crossing,
    fit_threshold,
    per_round_estimate,
)
from erasure_qec.circuits.builder import build
from erasure_qec.config import NoiseParams
from erasure_qec.decoding.herald_matching import (
    BlindMatchingDecoder,
    HeraldMatchingDecoder,
)
from erasure_qec.noise.injector import ErasureInjector
from erasure_qec.noise.model import nonidentity_pauli_probability

BASELINE = "tests/fixtures/real_baseline_pauli.csv"
R50 = "tests/fixtures/real_erasure_r50.csv"


def _contract_dem(circuit: stim.Circuit) -> stim.DetectorErrorModel:
    return circuit.detector_error_model(
        decompose_errors=True, flatten_loops=True, approximate_disjoint_errors=True
    ).flattened()


def _herald_detectors(dem: stim.DetectorErrorModel) -> set[int]:
    coords = dem.get_detector_coordinates()
    return {i for i, c in coords.items() if len(c) > 3 and c[3] == 1.0}


def _fit_diagnostics(csv: str, decoder: str, d_min: int | None) -> dict[str, float]:
    """Run the fit and recompute chi^2/dof, window max P_L_shot, and bound-pinning."""
    points = [p for p in load_sweep(csv) if p.decoder == decoder]
    fit = fit_threshold(points, d_min=d_min)
    filtered = [p for p in points if d_min is None or p.d >= d_min]
    center = estimate_crossing(filtered)
    used: list[SweepPoint] = _select_window(filtered, center, 1.5)
    out: dict[str, float] = {
        "p_th": fit.p_th,
        "nu": fit.nu,
        "p_th_err": fit.p_th_err,
        "converged": float(fit.converged),
        "n_points": float(len(used)),
    }
    if used and fit.converged:
        p_arr = np.array([p.p for p in used])
        d_arr = np.array([p.d for p in used], dtype=float)
        y = np.array([per_round_p_l(p.p_l_shot, p.rounds) for p in used])
        # Must match _weighted_fit_arrays: per_round_estimate is a 95% Wilson
        # interval, so divide the half-width by Z_95 to recover the 1-sigma
        # weight the fitter actually uses. Using the raw half-width here scales
        # chi^2 down by Z_95**2 ~ 3.84.
        sigma = np.array(
            [max((e.high - e.low) / (2.0 * Z_95), 1e-9) for e in map(per_round_estimate, used)]
        )
        params = (fit.p_th, fit.nu, *fit.coeffs)
        resid = (y - _model((p_arr, d_arr), *params)) / sigma
        dof = max(len(used) - 5, 1)
        out["chi2_dof"] = float(np.sum(resid**2) / dof)
        out["window_max_pl_shot"] = float(max(p.p_l_shot for p in used))
        out["p_arr_max"] = float(p_arr.max())
        out["pinned_at_bound"] = float(abs(fit.p_th - p_arr.max()) < 1e-9)
    return out


def _nonidentity_pauli_prob(r_e: float, p: float = 0.02) -> float:
    """Exact per-2q-gate probability of a non-identity Pauli, as a factor of p."""
    return nonidentity_pauli_probability(NoiseParams(p=p, r_e=r_e)) / p


def _dem_heralded_mass(d: int, p: float, r_e: float) -> float:
    """The **LEGACY** heralded-mass metric, which counts no-op herald branches.

    This row of the audit table compares against a pre-fix auditor expectation,
    so it must keep the pre-fix definition to stay comparable -- but that is the
    superseded metric. ``HERALDED_ERASE`` heralds on all four Pauli branches,
    identity included, so the decomposed DEM carries one mechanism per herald
    detector that causes no syndrome and no logical error; counting those in
    both numerator and denominator inflates the fraction. The corrected metric
    is the default of :func:`heralded_fraction`, which is what ``README.md`` and
    ``scripts/heralded_fraction.py`` report (0.454 vs 0.547 here).

    This used to be a local reimplementation, which meant running this script
    appeared to contradict the README with no indication the definitions
    differed. It now calls the library with the legacy flag set explicitly, so
    the choice is visible. Verified 2026-09-06 that the two agree to 1.5e-14
    (float summation order; not bit-identical).
    """
    return heralded_fraction(NoiseParams(p=p, r_e=r_e), d=d, count_no_op_heralds=True)


def _syndrome_density(d: int, p: float, r_e: float, shots: int = 20000) -> float:
    circ = build(d, d, ErasureInjector(NoiseParams(p=p, r_e=r_e)))
    dem = _contract_dem(circ)
    heralds = _herald_detectors(dem)
    synd = [i for i in range(circ.num_detectors) if i not in heralds]
    dets = circ.compile_detector_sampler(seed=0).sample(shots)
    return float(dets[:, synd].mean())


def _herald_free_fraction(d: int, p: float, r_e: float, shots: int = 40000) -> float:
    circ = build(d, d, ErasureInjector(NoiseParams(p=p, r_e=r_e)))
    dem = _contract_dem(circ)
    heralds = sorted(_herald_detectors(dem))
    dets = circ.compile_detector_sampler(seed=0).sample(shots)
    n_free = int((~dets[:, heralds].any(axis=1)).sum())
    return n_free / shots


def _ablation_ratio(d: int, p: float, r_e: float, shots: int = 40000) -> tuple[float, int, int]:
    circ = build(d, d, ErasureInjector(NoiseParams(p=p, r_e=r_e)))
    dets, obs = circ.compile_detector_sampler(seed=0).sample(shots, separate_observables=True)
    herald = HeraldMatchingDecoder.from_circuit(circ)
    blind = BlindMatchingDecoder(circ)
    h_pred = herald.decode_batch(dets)
    b_pred = blind.decode_batch(dets)
    h_err = int((h_pred[:, 0] != obs[:, 0]).sum())
    b_err = int((b_pred[:, 0] != obs[:, 0]).sum())
    h_pl = per_round_p_l(h_err / shots, d)
    b_pl = per_round_p_l(b_err / shots, d)
    ratio = b_pl / h_pl if h_pl > 0 else float("inf")
    return ratio, h_err, b_err


def main() -> None:
    rows: list[tuple[str, str, str]] = []

    def add(check: str, measured: str, expected: str) -> None:
        rows.append((check, measured, expected))

    b7 = _fit_diagnostics(BASELINE, "herald_mwpm", 7)
    add(
        "baseline herald d_min=7",
        f"p_th={b7['p_th'] * 100:.3f}%  nu={b7['nu']:.2f}",
        "p_th=1.376%, nu=1.45",
    )
    b_none = _fit_diagnostics(BASELINE, "herald_mwpm", None)
    add(
        "baseline herald d_min=None",
        f"p_th={b_none['p_th'] * 100:.3f}%  nu={b_none['nu']:.2f}",
        "p_th=1.708%, nu=2.22",
    )

    r7 = _fit_diagnostics(R50, "herald_mwpm", 7)
    add(
        "r50 herald d_min=7",
        f"p_th={r7['p_th'] * 100:.3f}% +/- {r7['p_th_err'] * 100:.3f}%",
        "p_th=2.185% +/- 0.278% (README says 2.30% +/- 0.15%)",
    )
    r_none = _fit_diagnostics(R50, "herald_mwpm", None)
    add(
        "r50 herald d_min=None",
        f"p_th={r_none['p_th'] * 100:.3f}%  "
        f"p_arr.max={r_none.get('p_arr_max', float('nan')) * 100:.3f}%  "
        f"pinned={bool(r_none.get('pinned_at_bound', 0))}  "
        f"chi2/dof={r_none.get('chi2_dof', float('nan')):.2f}  "
        f"converged={bool(r_none['converged'])}",
        "p_th=3.162% = p_arr.max, chi2/dof~2.89, converged=True",
    )
    r_none_blind = _fit_diagnostics(R50, "blind_mwpm", None)
    add("r50 blind d_min=None", f"nu={r_none_blind['nu']:.2f}", "nu=3.75 (unphysical; expect ~1.5)")

    chi2s = [
        _fit_diagnostics(BASELINE, "herald_mwpm", 7).get("chi2_dof", float("nan")),
        _fit_diagnostics(BASELINE, "blind_mwpm", 7).get("chi2_dof", float("nan")),
        _fit_diagnostics(R50, "herald_mwpm", 7).get("chi2_dof", float("nan")),
        _fit_diagnostics(R50, "blind_mwpm", 7).get("chi2_dof", float("nan")),
    ]
    add(
        "chi2/dof of d_min=7 fits",
        "  ".join(f"{c:.2f}" for c in chi2s),
        "0.11-0.32 (over-parameterised)",
    )

    wmax = [
        _fit_diagnostics(BASELINE, "herald_mwpm", 7).get("window_max_pl_shot", float("nan")),
        _fit_diagnostics(BASELINE, "herald_mwpm", None).get("window_max_pl_shot", float("nan")),
        _fit_diagnostics(R50, "herald_mwpm", 7).get("window_max_pl_shot", float("nan")),
        _fit_diagnostics(R50, "herald_mwpm", None).get("window_max_pl_shot", float("nan")),
    ]
    add(
        "max P_L_shot in fit window",
        f"max over 4 fits = {max(wmax):.3f}",
        "0.42-0.50 (coin-flip data fitted)",
    )

    # LEGACY metric (counts no-op herald branches) -- see _dem_heralded_mass.
    # The corrected metric, which README.md quotes, is printed beside it.
    add(
        "DEM heralded mass, LEGACY incl. no-op heralds (d=5,p=0.02,r_e=0.98)",
        f"{_dem_heralded_mass(5, 0.02, 0.98):.3f}",
        "0.475 (not 0.98)",
    )
    add(
        "DEM heralded mass, CORRECTED (= README; d=5,p=0.02,r_e=0.98)",
        f"{heralded_fraction(NoiseParams(p=0.02, r_e=0.98), d=5):.3f}",
        "n/a (post-fix metric, no pre-fix expectation)",
    )

    dens = [_syndrome_density(5, 0.02, re) for re in (0.0, 0.5, 0.98)]
    add(
        "syndrome density p=0.02 r_e=0/0.5/0.98",
        "  ".join(f"{x:.4f}" for x in dens),
        "0.2149 / 0.1993 / 0.1837 (-15%)",
    )

    npr = [_nonidentity_pauli_prob(re) for re in (0.0, 0.5, 0.98)]
    add(
        "non-identity Pauli /p  r_e=0/0.5/0.98",
        "  ".join(f"{x:.3f}" for x in npr),
        "1.0 / 0.875 / 0.755",
    )

    add(
        "herald-free shot fraction (d=5,p=0.02,r_e=0.98)",
        f"{_herald_free_fraction(5, 0.02, 0.98):.5f}",
        "~13/40000 = 0.00033",
    )

    for d in (3, 5, 7):
        ratio, h_err, b_err = _ablation_ratio(d, 0.01, 0.98)
        add(
            f"ablation ratio d={d} (p=1%,r_e=0.98)",
            f"{ratio:.1f}x  (herald_err={h_err}, blind_err={b_err})",
            {3: "1.9x", 5: "6.4x", 7: "~30x (n~22)"}[d],
        )

    width = max(len(c) for c, _, _ in rows)
    print(f"{'CHECK':<{width}}  {'MEASURED':<48}  EXPECTED")
    print("-" * (width + 48 + 20))
    for check, measured, expected in rows:
        print(f"{check:<{width}}  {measured:<48}  {expected}")


if __name__ == "__main__":
    main()
