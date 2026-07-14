"""Figure generation: byte-stability and partial-data robustness (PLAN.md §10, M8 gate)."""

import hashlib
from pathlib import Path

from erasure_qec.analysis.plotting import (
    figure_hook_regression,
    render_all,
)
from erasure_qec.analysis.synthetic import AnsatzParams, write_synthetic_csv

FIXTURES = Path(__file__).resolve().parent / "fixtures"
EXPECTED_FIGURES = {
    "threshold_panels.png",
    "lambda_vs_p.png",
    "ablation.png",
    "hook_regression.png",
}


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_render_all_produces_every_figure(tmp_path: Path) -> None:
    written = render_all(FIXTURES, tmp_path)
    assert {p.name for p in written} == EXPECTED_FIGURES
    assert all(p.exists() and p.stat().st_size > 0 for p in written)


def test_figures_regenerate_byte_stable(tmp_path: Path) -> None:
    """The M8 gate: identical CSVs -> byte-identical PNGs across runs."""
    first = {p.name: _digest(p) for p in render_all(FIXTURES, tmp_path / "a")}
    second = {p.name: _digest(p) for p in render_all(FIXTURES, tmp_path / "b")}
    assert first == second


def test_hook_figure_needs_no_csv(tmp_path: Path) -> None:
    out = figure_hook_regression(tmp_path / "hook.png", distances=(3, 5))
    assert out.exists() and out.stat().st_size > 0


def test_plotting_robust_to_empty_data(tmp_path: Path) -> None:
    """With no CSVs, only the hook figure renders (no crash on partial data)."""
    written = render_all(tmp_path / "empty_data", tmp_path / "figs")
    assert {p.name for p in written} == {"hook_regression.png"}


def test_plotting_robust_to_single_p_sweep(tmp_path: Path) -> None:
    """A one-p sweep can't fit a threshold, but curves/ablation still render."""
    ansatz = AnsatzParams(p_th=0.02, nu=1.5, a=0.09, b=22.0)
    write_synthetic_csv(
        tmp_path / "data" / "erasure_r50.csv",
        ansatz,
        decoders=("herald_mwpm", "blind_mwpm"),
        distances=(3, 5, 7),
        p_values=(0.01,),  # single point
        shots=100000,
        r_e=0.5,
    )
    written = render_all(tmp_path / "data", tmp_path / "figs")
    names = {p.name for p in written}
    # Threshold/collapse can't fit on one p, but ablation and hook always render.
    assert "hook_regression.png" in names
    assert "ablation.png" in names


def test_synthetic_fixtures_are_byte_stable(tmp_path: Path) -> None:
    """The committed fixtures regenerate identically (no sampling)."""
    ansatz = AnsatzParams(p_th=0.024, nu=1.5, a=0.09, b=22.0)
    kw = dict(
        decoders=("herald_mwpm", "blind_mwpm"),
        distances=(3, 5, 7, 9, 11),
        p_values=tuple(10.0 ** (-3.0 + i * (2.0 / 12.0)) for i in range(13)),
        shots=200_000,
        r_e=0.5,
        blind_penalty=1.35,
    )
    a = write_synthetic_csv(tmp_path / "a.csv", ansatz, **kw)  # type: ignore[arg-type]
    b = write_synthetic_csv(tmp_path / "b.csv", ansatz, **kw)  # type: ignore[arg-type]
    assert a.read_text() == b.read_text()
    # And matches the committed fixture for r50.
    committed = (FIXTURES / "erasure_r50.csv").read_text()
    assert a.read_text() == committed
