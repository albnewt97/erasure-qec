"""Lambda-factor scan collection driver (PLAN.md §1, §10).

Holds the physical error rate fixed at the config's sub-threshold ``lambda_p``
and sweeps the code distance, so ``Lambda = p_L(d) / p_L(d+2)`` can be read off
directly. One ``sinter.Task`` per (d, decoder) with the §6-contract DEM
attached; resumable exactly like the threshold sweep.

Usage:
    uv run python experiments/collect_lambda_scan.py experiments/configs/erasure_r50.yaml
    uv run python experiments/collect_lambda_scan.py <config.yaml> --workers 8
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
    """One task per distance at the fixed ``lambda_p``, T = d."""
    for d in config.distances:
        circuit = build(
            d, d, BiasedErasureInjector(NoiseParams(p=config.lambda_p, r_e=config.r_e))
        )
        yield sinter.Task(
            circuit=circuit,
            detector_error_model=contract_dem(circuit),
            json_metadata={
                "d": d,
                "rounds": d,
                "p": config.lambda_p,
                "r_e": config.r_e,
                "config": config.name,
                "scan": "lambda",
            },
        )


def collect(config: ExperimentConfig, num_workers: int) -> list[sinter.TaskStats]:
    """Run (or resume) the Lambda scan into ``data/<name>_lambda.csv``."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    save_path = DATA_DIR / f"{config.name}_lambda.csv"
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
    parser.add_argument("config", type=Path, help="experiments/configs/*.yaml definition")
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
