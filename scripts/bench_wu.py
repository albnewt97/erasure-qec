"""Single-process timing of the Wu et al. sampling circuit and decoders.

For each (d, p, R_e) point: instructions per circuit (flattened), sampler
compile time, sample time for ``--shots`` shots, decoding-DEM build time,
decoder compile time, and decode throughput (shots/s) for ``wu_herald_mwpm``
and ``wu_blind_mwpm``. Sampling and decoding are timed separately because the
31-branch chain per gate pair makes the sampling circuit much larger than the
decoding circuit. One comparison row times the repo-model ``herald_mwpm``
(``ErasureInjector``) at (15, 0.04, 0.98). Rounds are d. All timings are wall
clock (``time.perf_counter``), one repetition, one process.

Run:

    uv run python scripts/bench_wu.py              # 500 shots per point
    uv run python scripts/bench_wu.py --shots 50   # quick probe
"""

from __future__ import annotations

import argparse
import os
import platform
import subprocess
import time
from collections.abc import Callable
from dataclasses import dataclass
from typing import TypeVar

import numpy as np
import pymatching
import sinter
import stim

from erasure_qec.circuits.builder import build
from erasure_qec.config import NoiseParams
from erasure_qec.decoding.sinter_adapter import CUSTOM_DECODERS, WU_DECODERS, contract_dem
from erasure_qec.noise.injector import ErasureInjector
from erasure_qec.noise.wu import WuDecodingInjector, WuNoiseParams, WuSamplingInjector

POINTS: list[tuple[int, float, float]] = [
    (11, 0.04, 0.98),
    (15, 0.04, 0.98),
    (11, 0.05, 1.0),
    (15, 0.05, 1.0),
]
REPO_POINT = (15, 0.04, 0.98)
SEED = 2026

T = TypeVar("T")


@dataclass
class Row:
    """Timings (seconds unless noted) for one benchmark point."""

    label: str
    d: int
    p: float
    r_e: float
    instr_sampling: int
    instr_decoding: int
    t_build_sampling: float
    t_build_decoding: float
    t_sampler_compile: float
    t_sample: float
    t_dem: float
    t_compile_herald: float
    herald_shots_per_s: float
    t_compile_blind: float | None
    blind_shots_per_s: float | None


COLUMNS = (
    "model d p R_e instr_samp instr_dec build_samp_s build_dec_s sampler_compile_s "
    "sample_s dem_s compile_herald_s herald_shots/s compile_blind_s blind_shots/s"
)


def _timed(fn: Callable[[], T]) -> tuple[T, float]:
    t0 = time.perf_counter()
    out = fn()
    return out, time.perf_counter() - t0


def _decode_rate(
    decoder: sinter.Decoder, dem: stim.DetectorErrorModel, packed: np.ndarray
) -> tuple[float, float]:
    """(compile seconds, decode shots/s) for one sinter decoder on bit-packed shots."""
    compiled, t_compile = _timed(lambda: decoder.compile_decoder_for_dem(dem=dem))
    _, t_decode = _timed(
        lambda: compiled.decode_shots_bit_packed(bit_packed_detection_event_data=packed)
    )
    return t_compile, packed.shape[0] / t_decode


def _sample(circuit: stim.Circuit, shots: int) -> tuple[np.ndarray, float, float]:
    sampler, t_compile = _timed(lambda: circuit.compile_detector_sampler(seed=SEED))
    packed, t_sample = _timed(lambda: sampler.sample(shots, bit_packed=True))
    return packed, t_compile, t_sample


def bench_wu(d: int, p: float, r_e: float, shots: int) -> Row:
    params = WuNoiseParams(p, r_e)
    sampling, t_bs = _timed(lambda: build(d, d, WuSamplingInjector(d, params)))
    decoding, t_bd = _timed(lambda: build(d, d, WuDecodingInjector(d, params)))
    packed, t_sc, t_s = _sample(sampling, shots)
    dem, t_dem = _timed(
        lambda: decoding.detector_error_model(decompose_errors=True, flatten_loops=True)
    )
    t_ch, rate_h = _decode_rate(WU_DECODERS["wu_herald_mwpm"], dem, packed)
    t_cb, rate_b = _decode_rate(WU_DECODERS["wu_blind_mwpm"], dem, packed)
    return Row(
        "wu", d, p, r_e, len(sampling.flattened()), len(decoding.flattened()),
        t_bs, t_bd, t_sc, t_s, t_dem, t_ch, rate_h, t_cb, rate_b,
    )  # fmt: skip


def bench_repo(d: int, p: float, r_e: float, shots: int) -> Row:
    circuit, t_b = _timed(lambda: build(d, d, ErasureInjector(NoiseParams(p=p, r_e=r_e))))
    packed, t_sc, t_s = _sample(circuit, shots)
    dem, t_dem = _timed(lambda: contract_dem(circuit))
    t_ch, rate_h = _decode_rate(CUSTOM_DECODERS["herald_mwpm"], dem, packed)
    n = len(circuit.flattened())
    return Row("repo", d, p, r_e, n, n, t_b, t_b, t_sc, t_s, t_dem, t_ch, rate_h, None, None)


def _machine_info() -> str:
    cpu = platform.processor()
    if platform.system() == "Darwin":
        cpu = subprocess.run(
            ["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True
        ).stdout.strip()
    return (
        f"# machine: {platform.system()} {platform.release()} {platform.machine()}, "
        f"{cpu}, {os.cpu_count()} logical cores\n"
        f"# python {platform.python_version()}, stim {stim.__version__}, "
        f"pymatching {pymatching.__version__}, sinter {sinter.__version__}, "
        f"numpy {np.__version__}"
    )


def _fmt_row(r: Row) -> str:
    def f(x: float | None, spec: str) -> str:
        return "-" if x is None else format(x, spec)

    return " ".join(
        [
            r.label, str(r.d), f(r.p, ".3f"), f(r.r_e, ".2f"),
            str(r.instr_sampling), str(r.instr_decoding),
            f(r.t_build_sampling, ".2f"), f(r.t_build_decoding, ".2f"),
            f(r.t_sampler_compile, ".3f"), f(r.t_sample, ".3f"), f(r.t_dem, ".2f"),
            f(r.t_compile_herald, ".2f"), f(r.herald_shots_per_s, ".1f"),
            f(r.t_compile_blind, ".2f"), f(r.blind_shots_per_s, ".1f"),
        ]
    )  # fmt: skip


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shots", type=int, default=500)
    args = parser.parse_args()

    print(_machine_info())
    print(f"# shots per point: {args.shots}; rounds = d; seed {SEED}; single process")
    print("# repo row: build_dec_s repeats build_samp_s (one circuit); no blind decoder")
    print(COLUMNS, flush=True)
    for d, p, r_e in POINTS:
        print(_fmt_row(bench_wu(d, p, r_e, args.shots)), flush=True)
    print(_fmt_row(bench_repo(*REPO_POINT, args.shots)), flush=True)


if __name__ == "__main__":
    main()
