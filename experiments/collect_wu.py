"""Wu et al. (2022) validation collection driver (P5a).

Usage:
    uv run python experiments/collect_wu.py experiments/configs/wu/<name>.yaml --workers 8

Runs one config. One ``sinter.Task`` per (d, p) grid point, rounds = d, built
with :func:`~erasure_qec.decoding.sinter_adapter.wu_tasks_parallel` (sampling
circuit plus the tagged DEM of the decoding circuit). The build pool is closed
and joined before ``sinter.collect`` starts, so the build processes and the
sinter workers never hold memory at the same time. Decoders come from
``WU_DECODERS``. Collection is resumable: re-running a config appends to
``data/wu/<name>.csv`` and skips tasks that already reached ``max_shots`` or
``max_errors``.
"""

import argparse
import multiprocessing
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import sinter
import yaml

from erasure_qec.decoding.sinter_adapter import WU_DECODERS, wu_tasks_parallel
from erasure_qec.noise.wu import WuNoiseParams

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "wu"
CONVENTION = "wu2022"


@dataclass(frozen=True)
class WuConfig:
    """One Wu collection config, loaded from experiments/configs/wu/*.yaml.

    ``p_m`` is ``None`` exactly when ``p_m_follows_p`` is true (ancilla
    init/measurement flip probability equal to ``p`` at every grid point).
    """

    name: str
    R_e: float
    p_m: float | None
    p_m_follows_p: bool
    distances: tuple[int, ...]
    p_values: tuple[float, ...]
    max_shots: int
    max_errors: int
    decoders: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.p_m_follows_p != (self.p_m is None):
            raise ValueError("p_m must be null exactly when p_m_follows_p is true")
        if not self.distances or any(d < 3 or d % 2 == 0 for d in self.distances):
            raise ValueError(f"distances must be odd and >= 3, got {self.distances}")
        if not self.p_values or len(set(self.p_values)) != len(self.p_values):
            raise ValueError("p_values must be non-empty and distinct")
        if self.max_shots < 1 or self.max_errors < 1:
            raise ValueError("max_shots and max_errors must be >= 1")
        unknown = set(self.decoders) - set(WU_DECODERS)
        if not self.decoders or unknown:
            raise ValueError(f"decoders must be a non-empty subset of WU_DECODERS, got {unknown}")

    def params(self, p: float) -> WuNoiseParams:
        """Noise parameters at grid point ``p``."""
        if self.p_m is None:
            return WuNoiseParams(p, self.R_e, p_m_follows_p=True)
        return WuNoiseParams(p, self.R_e, p_m=self.p_m)


def load_wu_config(path: str | Path) -> WuConfig:
    """Parse one experiments/configs/wu/*.yaml file; requires ``convention: wu2022``.

    The optional ``extra_p_edge`` list (points appended by the edge rule,
    docs/archive/remediation/WU_DPAIR_DECISION.md) is merged with ``p_values``, sorted
    and de-duplicated, as ``extra_p`` is in
    :func:`erasure_qec.config.load_experiment_config`.
    """
    with open(path, encoding="utf-8") as f:
        raw = yaml.safe_load(f)
    if not isinstance(raw, dict):
        raise ValueError(f"{path}: expected a YAML mapping at top level")
    if raw.get("convention") != CONVENTION:
        raise ValueError(f"{path}: convention must be {CONVENTION!r}")
    base = tuple(float(p) for p in raw["p_values"])
    if len(set(base)) != len(base):
        raise ValueError(f"{path}: p_values must be distinct")
    edge = tuple(float(p) for p in raw.get("extra_p_edge", ()))
    return WuConfig(
        name=str(raw["name"]),
        R_e=float(raw["R_e"]),
        p_m=None if raw["p_m"] is None else float(raw["p_m"]),
        p_m_follows_p=bool(raw["p_m_follows_p"]),
        distances=tuple(int(d) for d in raw["distances"]),
        p_values=tuple(sorted(set(base) | set(edge))),
        max_shots=int(raw["max_shots"]),
        max_errors=int(raw["max_errors"]),
        decoders=tuple(str(dec) for dec in raw["decoders"]),
    )


def task_specs(config: WuConfig) -> list[tuple[int, WuNoiseParams, dict[str, Any]]]:
    """``(d, params, json_metadata)`` per (d, p) grid point, d-major, rounds = d.

    ``json_metadata["p_m"]`` is the effective ancilla flip probability
    (``WuNoiseParams.meas``), which equals ``p`` when ``p_m_follows_p``.
    """
    specs: list[tuple[int, WuNoiseParams, dict[str, Any]]] = []
    for d in config.distances:
        for p in config.p_values:
            params = config.params(p)
            meta = {
                "convention": CONVENTION,
                "config": config.name,
                "d": d,
                "rounds": d,
                "p": p,
                "R_e": config.R_e,
                "p_m": params.meas,
                "p_m_follows_p": config.p_m_follows_p,
            }
            specs.append((d, params, meta))
    return specs


def build_tasks(config: WuConfig, processes: int) -> list[sinter.Task]:
    """All tasks of ``config``, built in a pool that is joined before this returns."""
    tasks = wu_tasks_parallel(task_specs(config), processes=processes)
    if multiprocessing.active_children():
        raise RuntimeError("task-build pool processes are still alive")
    return tasks


def collect(
    config: WuConfig, num_workers: int, data_dir: Path = DATA_DIR
) -> list[sinter.TaskStats]:
    """Run (or resume) ``config``; results accumulate in ``<data_dir>/<name>.csv``."""
    data_dir.mkdir(parents=True, exist_ok=True)
    save_path = data_dir / f"{config.name}.csv"
    tasks = build_tasks(config, processes=num_workers)
    stats = sinter.collect(
        num_workers=num_workers,
        tasks=tasks,
        decoders=list(config.decoders),
        custom_decoders=WU_DECODERS,
        max_shots=config.max_shots,
        max_errors=config.max_errors,
        save_resume_filepath=save_path,
        print_progress=True,
    )
    print(f"wrote {save_path} ({len(stats)} task stats)")
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="experiments/configs/wu/*.yaml")
    parser.add_argument(
        "--workers",
        type=int,
        default=8,
        help="task-build pool processes and sinter workers (default: 8)",
    )
    args = parser.parse_args()
    collect(load_wu_config(args.config), num_workers=args.workers)


if __name__ == "__main__":
    main()
