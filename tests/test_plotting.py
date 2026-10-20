"""Figure generation: byte-stability and partial-data robustness.

docs/archive/PLAN_v1.md §10, M8 gate.
"""

import hashlib
import json
from pathlib import Path

from erasure_qec.analysis.plotting import (
    figure_hook_regression,
    render_all,
)
from erasure_qec.analysis.synthetic import AnsatzParams, write_synthetic_csv

FIXTURES = Path(__file__).resolve().parent / "fixtures"
# The committed synthetic fixtures are named synthetic_* (never mistaken for
# measurements); render_all keys on explicit config names, so pass them here.
SYNTHETIC = ("synthetic_baseline_pauli", "synthetic_erasure_r50", "synthetic_erasure_r98")
EXPECTED_FIGURES = {
    "threshold_panels.png",
    "lambda_vs_p.png",
    "ablation.png",
    "hook_regression.png",
}


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _ablation_json(path: Path) -> Path:
    """A minimal results/ablation.json (schema of scripts/ablation_table.py; not a measurement)."""
    cells = [
        {"p": p, "r_e": r_e, "d": d, "gated": d == 7, "ratio": 1.0 + d * r_e, "ratio_lb": 0.9 + d}
        for p in (0.005, 0.01)
        for r_e in (0.5, 0.98)
        for d in (3, 5, 7)
    ]
    path.write_text(json.dumps({"parameters": {"shots": 1000}, "cells": cells}))
    return path


def test_render_all_produces_every_figure(tmp_path: Path) -> None:
    ablation = _ablation_json(tmp_path / "ablation.json")
    written = render_all(FIXTURES, tmp_path, config_names=SYNTHETIC, ablation_json=ablation)
    assert {p.name for p in written} == EXPECTED_FIGURES
    assert all(p.exists() and p.stat().st_size > 0 for p in written)


def test_figures_regenerate_byte_stable(tmp_path: Path) -> None:
    """The M8 gate: identical inputs -> byte-identical PNGs across runs."""
    ablation = _ablation_json(tmp_path / "ablation.json")
    first = {
        p.name: _digest(p)
        for p in render_all(
            FIXTURES, tmp_path / "a", config_names=SYNTHETIC, ablation_json=ablation
        )
    }
    second = {
        p.name: _digest(p)
        for p in render_all(
            FIXTURES, tmp_path / "b", config_names=SYNTHETIC, ablation_json=ablation
        )
    }
    assert set(first) == EXPECTED_FIGURES
    assert first == second


def test_hook_figure_needs_no_csv(tmp_path: Path) -> None:
    out = figure_hook_regression(tmp_path / "hook.png", distances=(3, 5))
    assert out.exists() and out.stat().st_size > 0


def test_plotting_robust_to_empty_data(tmp_path: Path) -> None:
    """With no CSVs, only the hook figure renders (no crash on partial data)."""
    written = render_all(tmp_path / "empty_data", tmp_path / "figs")
    assert {p.name for p in written} == {"hook_regression.png"}


def test_plotting_robust_to_single_p_sweep(tmp_path: Path) -> None:
    """A one-p sweep has no crossing, but curves and the ablation still render."""
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
    ablation = _ablation_json(tmp_path / "ablation.json")
    written = render_all(tmp_path / "data", tmp_path / "figs", ablation_json=ablation)
    names = {p.name for p in written}
    # No crossing on one p, but ablation and hook always render.
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
    committed = (FIXTURES / "synthetic_erasure_r50.csv").read_text()
    assert a.read_text() == committed


def test_ablation_figure_skipped_without_results(tmp_path: Path) -> None:
    """No ablation.json: the figure is skipped and no file is written."""
    written = render_all(
        FIXTURES, tmp_path, config_names=SYNTHETIC, ablation_json=tmp_path / "missing.json"
    )
    assert "ablation.png" not in {p.name for p in written}
    assert not (tmp_path / "ablation.png").exists()
