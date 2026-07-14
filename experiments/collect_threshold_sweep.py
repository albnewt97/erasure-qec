"""Threshold-sweep collection driver (PLAN.md §9).

Usage:
    uv run python experiments/collect_threshold_sweep.py experiments/configs/erasure_r50.yaml
    uv run python experiments/collect_threshold_sweep.py <config.yaml> --workers 8

Builds one ``sinter.Task`` per (d, p) point — rounds T = d, biased-erasure
noise at the config's ``r_e`` — with the §6-contract DEM attached explicitly
(``decompose_errors=True, flatten_loops=True, approximate_disjoint_errors=True``;
sinter does not add these flags itself and the Design-2 partition requires
them). Collection is resumable: re-running with the same config appends to
``data/<name>.csv`` and skips already-satisfied tasks.
"""

import argparse
import os
from collections.abc import Iterator
from pathlib import Path

import sinter

from erasure_qec.circuits.builder import build
from erasure_qec.config import ExperimentConfig, NoiseParams, load_experiment_config
from erasure_qec.decoding.sinter_adapter import CUSTOM_DECODERS, contract_dem
from erasure_qec.noise.injector import BiasedErasureInjector

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def tasks_from_config(config: ExperimentConfig) -> Iterator[sinter.Task]:
    """One task per (distance, physical error rate) grid point, T = d."""
    for d in config.distances:
        for p in config.p_values:
            circuit = build(
                d, d, BiasedErasureInjector(NoiseParams(p=p, r_e=config.r_e))
            )
            yield sinter.Task(
                circuit=circuit,
                detector_error_model=contract_dem(circuit),
                json_metadata={
                    "d": d,
                    "rounds": d,
                    "p": p,
                    "r_e": config.r_e,
                    "config": config.name,
                },
            )


def collect(config: ExperimentConfig, num_workers: int) -> list[sinter.TaskStats]:
    """Run (or resume) the sweep; results accumulate in ``data/<name>.csv``."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    save_path = DATA_DIR / f"{config.name}.csv"
    stats = sinter.collect(
        num_workers=num_workers,
        tasks=list(tasks_from_config(config)),
        decoders=list(config.decoders),
        custom_decoders=CUSTOM_DECODERS,
        max_shots=config.max_shots,
        max_errors=config.max_errors,
        save_resume_filepath=save_path,
        print_progress=True,
    )
    print(f"wrote {save_path} ({len(stats)} task stats)")
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="experiments/configs/*.yaml sweep definition")
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, (os.cpu_count() or 2) - 1),
        help="sinter worker processes (default: cpu_count - 1)",
    )
    args = parser.parse_args()
    collect(load_experiment_config(args.config), num_workers=args.workers)


if __name__ == "__main__":
    main()
