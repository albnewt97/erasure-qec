"""Experiment-config loading: the p grid and its additive dense bands."""

from pathlib import Path

import pytest
import yaml

from erasure_qec.config import load_experiment_config

_BASE = """\
name: test_cfg
r_e: 0.5
distances: [3, 5]
p_log10_min: -3.0
p_log10_max: -1.0
p_points: 13
max_shots: 1000
max_errors: 10
decoders: [herald_mwpm, blind_mwpm]
lambda_p: 0.01
"""


def _log_grid() -> list[float]:
    return [10.0 ** (-3.0 + i * (2.0 / 12)) for i in range(13)]


def _write(tmp_path: Path, extra: str) -> Path:
    path = tmp_path / "cfg.yaml"
    path.write_text(_BASE + extra, encoding="utf-8")
    return path


def test_without_extra_p_shot_loads_as_before(tmp_path: Path) -> None:
    """No `extra_p_shot` key: p_values is the log grid merged with `extra_p` only."""
    cfg = load_experiment_config(_write(tmp_path, "extra_p: [0.02, 0.011, 0.02]\n"))
    assert cfg.p_values == tuple(sorted(set(_log_grid()) | {0.011, 0.02}))


def test_extra_p_shot_is_merged_like_extra_p(tmp_path: Path) -> None:
    """Both bands: union with the log grid, sorted, de-duplicated (overlaps included)."""
    cfg = load_experiment_config(
        _write(
            tmp_path,
            "extra_p: [0.011, 0.02]\nextra_p_shot: [0.0079, 0.011, 0.0079, 0.1]\n",
        )
    )
    expected = tuple(sorted(set(_log_grid()) | {0.011, 0.02, 0.0079}))
    assert cfg.p_values == expected
    assert list(cfg.p_values) == sorted(set(cfg.p_values))  # sorted, no duplicates
    assert 0.1 in cfg.p_values  # coincides with the log grid's last point; kept once


def test_extra_p_shot_alone(tmp_path: Path) -> None:
    cfg = load_experiment_config(_write(tmp_path, "extra_p_shot: [0.0123]\n"))
    assert cfg.p_values == tuple(sorted(set(_log_grid()) | {0.0123}))


@pytest.mark.parametrize("name", ["baseline_pauli", "erasure_r50", "erasure_r98"])
def test_committed_configs_contain_every_extra_band(name: str) -> None:
    path = Path(__file__).resolve().parent.parent / "experiments" / "configs" / f"{name}.yaml"
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    cfg = load_experiment_config(path)
    for key in ("extra_p", "extra_p_shot"):
        for p in raw.get(key, ()):
            assert float(p) in cfg.p_values
