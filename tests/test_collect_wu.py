"""experiments/collect_wu.py: smoke run, resume safety, strong_id stability, configs.

The smoke and resume tests run ``sinter.collect`` on a d = 3 config with two
p values and 200 shots per task, writing only under ``tmp_path``.
"""

from __future__ import annotations

import importlib.util
import multiprocessing
import sys
from pathlib import Path
from typing import Any

import pytest
import sinter

from erasure_qec.decoding.sinter_adapter import wu_task

_ROOT = Path(__file__).resolve().parents[1]
_SCRIPT = _ROOT / "experiments" / "collect_wu.py"
_spec = importlib.util.spec_from_file_location("collect_wu", _SCRIPT)
assert _spec is not None and _spec.loader is not None
collect_wu = importlib.util.module_from_spec(_spec)
sys.modules["collect_wu"] = collect_wu
_spec.loader.exec_module(collect_wu)

CONFIG_DIR = _ROOT / "experiments" / "configs" / "wu"
EXPECTED_COUNTS = {
    "wu_fig3_re000": 11,
    "wu_fig3_re098": 17,
    "wu_fig3b_re100": 19,
    "wu_fig3b_re050": 11,
    "wu_fig3b_re090": 11,
    "wu_si_pm005": 17,
    "wu_si_pmp": 15,
}
META_KEYS = {"convention", "config", "d", "rounds", "p", "R_e", "p_m", "p_m_follows_p"}

TINY_YAML = """\
name: wu_tiny
convention: wu2022
R_e: 0.98
p_m: 0.005
p_m_follows_p: false
distances: [3]
p_values: [0.03, 0.04]
max_shots: 200
max_errors: 1000
decoders: [wu_herald_mwpm, wu_blind_mwpm]
"""


def _tiny_config(tmp_path: Path) -> Any:
    path = tmp_path / "wu_tiny.yaml"
    path.write_text(TINY_YAML, encoding="utf-8")
    return collect_wu.load_wu_config(path)


def _shots_by_strong_id(csv_path: Path) -> dict[str, int]:
    return {s.strong_id: s.shots for s in sinter.read_stats_from_csv_files(str(csv_path))}


def test_smoke_run_writes_every_task(tmp_path: Path) -> None:
    config = _tiny_config(tmp_path)
    collect_wu.collect(config, num_workers=2, data_dir=tmp_path / "out")
    stats = sinter.read_stats_from_csv_files(str(tmp_path / "out" / "wu_tiny.csv"))
    keys = sorted((s.decoder, s.json_metadata["d"], s.json_metadata["p"]) for s in stats)
    assert keys == sorted(
        (dec, 3, p) for dec in ("wu_herald_mwpm", "wu_blind_mwpm") for p in (0.03, 0.04)
    )
    for s in stats:
        assert s.shots == 200
        assert set(s.json_metadata) == META_KEYS
        assert s.json_metadata["convention"] == "wu2022"
        assert s.json_metadata["config"] == "wu_tiny"
        assert s.json_metadata["rounds"] == 3
        assert s.json_metadata["R_e"] == 0.98
        assert s.json_metadata["p_m"] == 0.005
        assert s.json_metadata["p_m_follows_p"] is False


def test_second_run_adds_no_tasks(tmp_path: Path) -> None:
    config = _tiny_config(tmp_path)
    out = tmp_path / "out"
    csv_path = out / "wu_tiny.csv"
    collect_wu.collect(config, num_workers=2, data_dir=out)
    first = _shots_by_strong_id(csv_path)
    collect_wu.collect(config, num_workers=2, data_dir=out)
    second = _shots_by_strong_id(csv_path)
    assert len(first) == 4
    assert second == first  # same strong_ids, and no shots added to any of them


def test_build_pool_is_joined_before_collect(tmp_path: Path) -> None:
    config = _tiny_config(tmp_path)
    tasks = collect_wu.build_tasks(config, processes=2)
    assert len(tasks) == 2
    assert multiprocessing.active_children() == []


def _strong_id_in_child(config_path: str, index: int, decoder: str) -> str:
    """Build task ``index`` of the config in this (fresh) process; its sinter strong_id."""
    config = collect_wu.load_wu_config(config_path)
    task = wu_task(*collect_wu.task_specs(config)[index])
    return str(
        sinter.Task(
            circuit=task.circuit,
            decoder=decoder,
            detector_error_model=task.detector_error_model,
            json_metadata=task.json_metadata,
        ).strong_id()
    )


def test_strong_id_identical_across_spawn_processes(tmp_path: Path) -> None:
    path = tmp_path / "wu_tiny.yaml"
    path.write_text(TINY_YAML, encoding="utf-8")
    ctx = multiprocessing.get_context("spawn")
    ids = []
    for _ in range(2):  # two separate pools: two separate spawned processes
        with ctx.Pool(1) as pool:
            ids.append(pool.apply(_strong_id_in_child, (str(path), 1, "wu_herald_mwpm")))
    assert ids[0] == ids[1]


def test_p_m_follows_p_metadata_records_effective_p_m() -> None:
    config = collect_wu.load_wu_config(CONFIG_DIR / "wu_si_pmp.yaml")
    for _, params, meta in collect_wu.task_specs(config):
        assert params.p_m_follows_p
        assert meta["p_m"] == meta["p"] == params.p
        assert meta["p_m_follows_p"] is True


@pytest.mark.parametrize("name", sorted(EXPECTED_COUNTS))
def test_committed_configs_load_with_pre_registered_settings(name: str) -> None:
    config = collect_wu.load_wu_config(CONFIG_DIR / f"{name}.yaml")
    assert config.name == name
    assert len(config.p_values) == EXPECTED_COUNTS[name]
    assert all(p == round(p, 4) for p in config.p_values)
    assert config.distances == (7, 9, 13)
    assert config.max_shots == 100_000
    assert config.max_errors == 1000
    expected = (
        ("wu_herald_mwpm",) if name == "wu_fig3_re000" else ("wu_herald_mwpm", "wu_blind_mwpm")
    )
    assert config.decoders == expected


def test_config_dir_holds_exactly_the_seven_configs() -> None:
    assert {p.stem for p in CONFIG_DIR.glob("*.yaml")} == set(EXPECTED_COUNTS)


# Edge rule for scaled grids, docs/remediation/WU_DPAIR_DECISION.md (2026-10-16).
SCALED_C = {"wu_fig3b_re050": 0.021946, "wu_fig3b_re090": 0.032103}
EDGE_UPPER_F = (1.295, 1.340, 1.385)
EDGE_LOWER_F = (0.755, 0.710, 0.665)


@pytest.mark.parametrize("name", sorted(SCALED_C))
def test_scaled_edge_extensions_are_new_points(name: str) -> None:
    grid = collect_wu.load_wu_config(CONFIG_DIR / f"{name}.yaml").p_values
    c = SCALED_C[name]
    upper = [round(c * f, 4) for f in EDGE_UPPER_F]
    lower = [round(c * f, 4) for f in EDGE_LOWER_F]
    extension = upper + lower
    assert len(set(extension)) == 6
    assert not set(extension) & set(grid)
    assert min(upper) > max(grid) and max(lower) < min(grid)
