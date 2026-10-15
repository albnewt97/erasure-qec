"""Wall-time and memory projection for the P5 Wu collection (Task 6, pre-registered).

Three subcommands, run in this order on an otherwise idle machine:

    uv run python scripts/wu_projection.py profile --out DIR
    uv run python scripts/wu_projection.py pilot --config NAME --out DIR   # once per config
    uv run python scripts/wu_projection.py project --dir DIR

``profile`` (6a): per d, one wu_fig3_re098 task at p = 0.042, built in a fresh
process: build times (sampling, decoding, DEM), ``str()`` and ``strong_id``
time, ``len(sampling)``, ``len(str())`` in MB of the sampling circuit and the
DEM, peak RSS for building one task, and held RSS per task (slope of peak RSS
over three tasks held at once). Worker memory: peak RSS of a fresh process that
receives a pickled task (as a sinter worker does), compiles the sampler and
decoders, and decodes ``PILOT_SHOTS`` shots. Also measures the effective
parallelism of 8 concurrent processes for task builds and for decoding.

``pilot`` (6c shots model): for one config, at d in {7, 9, 11, 13, 15} and at
the lowest and highest p of its grid, times the task build, sampler compile,
sampling, and decoder compile and decode of ``PILOT_SHOTS`` shots. Herald
failures are counted only at d in {13, 15} (the candidate d_hi). Pilot data are
planning-only scratch: they are not written to data/ and not used in any result.

``project`` (6c/6d): applies the pre-registered formula and decision rule.
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing
import pickle
import platform
import resource
import subprocess
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pymatching
import sinter
import stim

from erasure_qec.circuits.builder import build
from erasure_qec.decoding.sinter_adapter import WU_DECODERS
from erasure_qec.noise.wu import WuDecodingInjector, WuNoiseParams, WuSamplingInjector

MAX_SHOTS = 100_000
MAX_ERRORS = 1000
WORKERS = 8
PILOT_SHOTS = 1000
DISTANCES = (7, 9, 11, 13, 15)
PAIRS = {"(11,15)": (11, 15), "(9,13)": (9, 13)}
CROSS_CHECK_D = 7
HOURS_LIMIT = 12.0
RAM_FRACTION_LIMIT = 0.5
DEFER_ORDER = ("wu_si_pmp", "wu_si_pm005", "wu_fig3b_re090")
SEED = 2026
BOTH = ("wu_herald_mwpm", "wu_blind_mwpm")


@dataclass(frozen=True)
class Config:
    """One P5a config: noise, p grid, decoders."""

    name: str
    r_e: float
    p_m: float | None  # None: p_m follows p
    grid: tuple[float, ...]
    decoders: tuple[str, ...]
    note: str = ""


def _steps(start: float, stop: float, step: float) -> tuple[float, ...]:
    n = round((stop - start) / step)
    return tuple(round(start + i * step, 6) for i in range(n + 1))


def _scaled(c: float) -> tuple[float, ...]:
    return tuple(round(c * (0.800 + 0.045 * i), 6) for i in range(11))


CONFIGS: tuple[Config, ...] = (
    Config("wu_fig3_re000", 0.0, 0.0, _steps(0.0075, 0.0125, 0.0005), ("wu_herald_mwpm",)),
    Config("wu_fig3_re098", 0.98, 0.0, _steps(0.034, 0.050, 0.001), BOTH),
    Config("wu_fig3b_re100", 1.0, 0.0, _steps(0.042, 0.060, 0.001), BOTH),
    Config("wu_si_pm005", 0.98, 0.005, _steps(0.030, 0.046, 0.001), BOTH),
    Config("wu_si_pmp", 0.98, None, _steps(0.022, 0.036, 0.001), BOTH),
    Config("wu_fig3b_re050", 0.5, 0.0, _scaled(0.020), BOTH, "c(0.5)=0.020 planning value"),
    Config("wu_fig3b_re090", 0.9, 0.0, _scaled(0.036), BOTH, "c(0.9)=0.036 planning value"),
)
BY_NAME = {c.name: c for c in CONFIGS}


def params_for(cfg: Config, p: float) -> WuNoiseParams:
    if cfg.p_m is None:
        return WuNoiseParams(p, cfg.r_e, p_m_follows_p=True)
    return WuNoiseParams(p, cfg.r_e, p_m=cfg.p_m)


def _peak_rss_mb() -> float:
    """Peak RSS of this process in MB (macOS reports bytes, Linux kilobytes)."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / 2**20 if platform.system() == "Darwin" else peak / 2**10


def _timed(fn: Callable[[], Any]) -> tuple[Any, float]:
    t0 = time.perf_counter()
    out = fn()
    return out, time.perf_counter() - t0


def _build_task(cfg: Config, d: int, p: float) -> tuple[sinter.Task, dict[str, float]]:
    """``wu_task`` with each stage timed separately."""
    params = params_for(cfg, p)
    sampling, t_s = _timed(lambda: build(d, d, WuSamplingInjector(d, params)))
    decoding, t_d = _timed(lambda: build(d, d, WuDecodingInjector(d, params)))
    dem, t_dem = _timed(
        lambda: decoding.detector_error_model(decompose_errors=True, flatten_loops=True)
    )
    task = sinter.Task(circuit=sampling, detector_error_model=dem, json_metadata={"d": d, "p": p})
    return task, {"build_sampling_s": t_s, "build_decoding_s": t_d, "dem_s": t_dem}


def _in_fresh_process(fn: Callable[..., dict[str, Any]], *args: Any) -> dict[str, Any]:
    with multiprocessing.get_context("spawn").Pool(1, maxtasksperchild=1) as pool:
        return pool.apply(fn, args)


# --- profile (6a) ---


def _profile_build(d: int) -> dict[str, Any]:
    cfg, p = BY_NAME["wu_fig3_re098"], 0.042
    base = _peak_rss_mb()
    task, times = _build_task(cfg, d, p)
    peak_one = _peak_rss_mb()
    assert task.circuit is not None and task.detector_error_model is not None
    (s_circ, s_dem), t_str = _timed(lambda: (str(task.circuit), str(task.detector_error_model)))
    keyed = sinter.Task(
        circuit=task.circuit,
        decoder="wu_herald_mwpm",
        detector_error_model=task.detector_error_model,
        json_metadata=task.json_metadata,
    )
    _, t_id = _timed(keyed.strong_id)
    held = [task]
    peaks = [peak_one]
    for _ in range(2):
        held.append(_build_task(cfg, d, p)[0])
        peaks.append(_peak_rss_mb())
    return {
        "d": d,
        **times,
        "str_s": t_str,
        "strong_id_s": t_id,
        "len_sampling": len(task.circuit),
        "str_sampling_mb": len(s_circ) / 2**20,
        "str_dem_mb": len(s_dem) / 2**20,
        "base_rss_mb": base,
        "peak_rss_one_task_mb": peak_one,
        "held_mb_per_task": (peaks[-1] - peaks[0]) / 2,
    }


def _worker_peak(task_bytes: bytes, decoders: tuple[str, ...]) -> dict[str, Any]:
    task = pickle.loads(task_bytes)
    after_unpickle = _peak_rss_mb()
    packed = task.circuit.compile_detector_sampler(seed=SEED).sample(PILOT_SHOTS, bit_packed=True)
    for name in decoders:
        compiled = WU_DECODERS[name].compile_decoder_for_dem(dem=task.detector_error_model)
        compiled.decode_shots_bit_packed(bit_packed_detection_event_data=packed)
    return {"after_unpickle_mb": after_unpickle, "peak_mb": _peak_rss_mb()}


def _decode_alone(task_bytes: bytes) -> float:
    task = pickle.loads(task_bytes)
    packed = task.circuit.compile_detector_sampler(seed=SEED).sample(500, bit_packed=True)
    compiled = WU_DECODERS["wu_herald_mwpm"].compile_decoder_for_dem(dem=task.detector_error_model)
    _, t = _timed(lambda: compiled.decode_shots_bit_packed(bit_packed_detection_event_data=packed))
    return float(t)


def _build_alone(d: int) -> float:
    return float(_build_task(BY_NAME["wu_fig3_re098"], d, 0.042)[1]["build_sampling_s"])


def _parallel_efficiency(fn: Callable[[Any], float], arg: Any) -> dict[str, float]:
    """Mean per-job time alone vs with WORKERS concurrent jobs: efficiency = alone / concurrent."""
    ctx = multiprocessing.get_context("spawn")
    with ctx.Pool(1) as pool:
        alone = pool.apply(fn, (arg,))
    with ctx.Pool(WORKERS) as pool:
        concurrent = pool.map(fn, [arg] * WORKERS)
    mean_c = float(np.mean(concurrent))
    return {"alone_s": alone, "concurrent_mean_s": mean_c, "efficiency": alone / mean_c}


def cmd_profile(out: Path) -> None:
    builds = [_in_fresh_process(_profile_build, d) for d in DISTANCES]
    workers = []
    for name, p in (("wu_fig3_re098", 0.050), ("wu_fig3b_re100", 0.060)):
        for d in (13, 15):
            task, _ = _build_task(BY_NAME[name], d, p)
            res = _in_fresh_process(_worker_peak, pickle.dumps(task), BOTH)
            workers.append({"config": name, "d": d, "p": p, **res})
    task13, _ = _build_task(BY_NAME["wu_fig3_re098"], 13, 0.042)
    eff = {
        "build_d11": _parallel_efficiency(_build_alone, 11),
        "decode_d13": _parallel_efficiency(_decode_alone, pickle.dumps(task13)),
    }
    (out / "profile.json").write_text(
        json.dumps({"builds": builds, "workers": workers, "parallel": eff}, indent=1)
    )


# --- pilot (6c shots model) ---


def cmd_pilot(cfg: Config, out: Path) -> None:
    rows = []
    for d in DISTANCES:
        for p in (cfg.grid[0], cfg.grid[-1]):
            task, times = _build_task(cfg, d, p)
            assert task.circuit is not None and task.detector_error_model is not None
            sampler, t_sc = _timed(lambda c=task.circuit: c.compile_detector_sampler(seed=SEED))
            (dets, obs), t_s = _timed(
                lambda s=sampler: s.sample(PILOT_SHOTS, separate_observables=True, bit_packed=True)
            )
            row: dict[str, Any] = {"d": d, "p": p, **times, "sampler_compile_s": t_sc}
            row["sample_s_per_shot"] = t_s / PILOT_SHOTS
            for name in cfg.decoders:
                compiled, t_c = _timed(
                    lambda n=name, t=task: WU_DECODERS[n].compile_decoder_for_dem(
                        dem=t.detector_error_model
                    )
                )
                pred, t_dec = _timed(
                    lambda c=compiled, x=dets: c.decode_shots_bit_packed(
                        bit_packed_detection_event_data=x
                    )
                )
                row[f"{name}_compile_s"] = t_c
                row[f"{name}_s_per_shot"] = t_dec / PILOT_SHOTS
                if name == "wu_herald_mwpm" and d in (13, 15):
                    row["herald_failures"] = int((pred != obs).any(axis=1).sum())
            rows.append(row)
            print(json.dumps(row), flush=True)
    (out / f"pilot_{cfg.name}.json").write_text(json.dumps(rows, indent=1))


# --- project (6c / 6d) ---


def _lerp(x: float, x0: float, x1: float, y0: float, y1: float) -> float:
    return y0 if x1 == x0 else y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def _p_l(row: dict[str, Any]) -> tuple[float, bool]:
    """Pilot P_L at one end; zero failures -> 0.5 / PILOT_SHOTS, flagged."""
    f = row["herald_failures"]
    return (f / PILOT_SHOTS, False) if f > 0 else (0.5 / PILOT_SHOTS, True)


def shots_for(cfg: Config, rows: list[dict[str, Any]], d_hi: int, p: float) -> int:
    """min(MAX_SHOTS, ceil(MAX_ERRORS / P_L(p))), P_L log-linear between the d_hi pilots."""
    lo = next(r for r in rows if r["d"] == d_hi and r["p"] == cfg.grid[0])
    hi = next(r for r in rows if r["d"] == d_hi and r["p"] == cfg.grid[-1])
    log_pl = _lerp(p, lo["p"], hi["p"], math.log(_p_l(lo)[0]), math.log(_p_l(hi)[0]))
    return min(MAX_SHOTS, math.ceil(MAX_ERRORS / math.exp(log_pl)))


def _interp(rows: list[dict[str, Any]], d: int, p: float, key: str) -> float:
    lo, hi = sorted((r for r in rows if r["d"] == d), key=lambda r: r["p"])
    return _lerp(p, lo["p"], hi["p"], lo[key], hi[key])


@dataclass
class ConfigProjection:
    name: str
    n_tasks: int
    total_shots: int
    build_h: float
    shots_h: float
    compile_h: float
    peak_ram_gb: float
    held_gb: float
    workers_gb: float
    total_h_measured_eff: float

    @property
    def total_h(self) -> float:
        return self.build_h + self.shots_h + self.compile_h


def project_config(
    cfg: Config, rows: list[dict[str, Any]], profile: dict[str, Any], pair: tuple[int, int]
) -> ConfigProjection:
    d_lo, d_hi = pair
    builds = {b["d"]: b for b in profile["builds"]}
    worker_peak = max(w["peak_mb"] for w in profile["workers"] if w["d"] == d_hi)
    eff = profile["parallel"]["decode_d13"]["efficiency"]
    build_s = shots_cpu_s = compile_s = 0.0
    held_mb = 0.0
    total_shots = 0
    n_tasks = 0
    for d in (CROSS_CHECK_D, d_lo, d_hi):
        b = builds[d]
        held_mb += len(cfg.grid) * max(
            b["held_mb_per_task"], b["str_sampling_mb"] + b["str_dem_mb"]
        )
        for p in cfg.grid:
            n_tasks += 1
            shots = shots_for(cfg, rows, d_hi, p)
            total_shots += shots * len(cfg.decoders)
            build_s += sum(
                _interp(rows, d, p, k) for k in ("build_sampling_s", "build_decoding_s", "dem_s")
            )
            for name in cfg.decoders:
                per_shot = _interp(rows, d, p, "sample_s_per_shot") + _interp(
                    rows, d, p, f"{name}_s_per_shot"
                )
                shots_cpu_s += shots * per_shot
                compile_s += _interp(rows, d, p, f"{name}_compile_s") + _interp(
                    rows, d, p, "sampler_compile_s"
                )
    base = builds[d_hi]["base_rss_mb"]
    pool_phase = base + held_mb + WORKERS * builds[d_hi]["peak_rss_one_task_mb"]
    collect_phase = base + held_mb + WORKERS * worker_peak
    return ConfigProjection(
        cfg.name,
        n_tasks,
        total_shots,
        build_s / WORKERS / 3600,
        shots_cpu_s / WORKERS / 3600,
        compile_s / 3600,
        max(pool_phase, collect_phase) / 1024,
        held_mb / 1024,
        WORKERS * max(worker_peak, builds[d_hi]["peak_rss_one_task_mb"]) / 1024,
        (
            build_s / (WORKERS * profile["parallel"]["build_d11"]["efficiency"])
            + shots_cpu_s / (WORKERS * eff)
            + compile_s
        )
        / 3600,
    )


def decide(
    projections: dict[str, dict[str, ConfigProjection]], ram_gb: float
) -> tuple[str, list[str]]:
    """Pre-registered rule 6d: returns (d-pair, deferred configs)."""

    def ok(pair: str, deferred: list[str]) -> bool:
        runs = [p for n, p in projections[pair].items() if n not in deferred]
        return sum(p.total_h for p in runs) <= HOURS_LIMIT and all(
            p.peak_ram_gb <= RAM_FRACTION_LIMIT * ram_gb for p in runs
        )

    if ok("(11,15)", []):
        return "(11,15)", []
    if ok("(9,13)", []):
        return "(9,13)", []
    deferred: list[str] = []
    for name in DEFER_ORDER:
        deferred.append(name)
        if ok("(9,13)", deferred):
            break
    return "(9,13)", deferred


def _machine() -> tuple[str, float]:
    if platform.system() == "Darwin":
        cpu = subprocess.run(
            ["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True
        ).stdout.strip()
        ram = int(
            subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True).stdout
        )
    else:
        cpu, ram = platform.processor(), 0
    info = (
        f"# machine: {platform.system()} {platform.machine()}, {cpu}, "
        f"{multiprocessing.cpu_count()} logical cores, RAM {ram / 2**30:.1f} GiB\n"
        f"# python {platform.python_version()}, stim {stim.__version__}, "
        f"pymatching {pymatching.__version__}, sinter {sinter.__version__}"
    )
    return info, ram / 2**30


def cmd_project(directory: Path) -> None:
    info, ram_gb = _machine()
    profile = json.loads((directory / "profile.json").read_text())
    pilots = {c.name: json.loads((directory / f"pilot_{c.name}.json").read_text()) for c in CONFIGS}
    print(info)
    print(
        f"# max_shots {MAX_SHOTS}, max_errors {MAX_ERRORS}, workers {WORKERS}, "
        f"pilot shots {PILOT_SHOTS}, d = {{{CROSS_CHECK_D}, d_lo, d_hi}}"
    )

    print("\n## 6a profile: one wu_fig3_re098 task at p = 0.042 (fresh process per d)")
    print(
        "d build_samp_s build_dec_s dem_s str_s strong_id_s len_sampling str_samp_MB "
        "str_dem_MB peak_rss_one_task_MB held_MB_per_task"
    )
    for b in profile["builds"]:
        print(
            f"{b['d']} {b['build_sampling_s']:.2f} {b['build_decoding_s']:.2f} {b['dem_s']:.2f} "
            f"{b['str_s']:.2f} {b['strong_id_s']:.2f} {b['len_sampling']} "
            f"{b['str_sampling_mb']:.1f} {b['str_dem_mb']:.1f} "
            f"{b['peak_rss_one_task_mb']:.0f} {b['held_mb_per_task']:.0f}"
        )
    print(
        "worker peak RSS (pickled task, sampler + both decoders, "
        f"{PILOT_SHOTS} shots): config d p after_unpickle_MB peak_MB"
    )
    for w in profile["workers"]:
        print(f"  {w['config']} {w['d']} {w['p']} {w['after_unpickle_mb']:.0f} {w['peak_mb']:.0f}")
    for k, v in profile["parallel"].items():
        print(
            f"parallel {k}: alone {v['alone_s']:.2f}s, {WORKERS} concurrent mean "
            f"{v['concurrent_mean_s']:.2f}s, efficiency {v['efficiency']:.2f}"
        )

    print("\n## pilot P_L (herald, planning-only scratch; zero failures flagged *)")
    print("config d_hi p failures P_L")
    for c in CONFIGS:
        for d in (13, 15):
            for r in (r for r in pilots[c.name] if r["d"] == d):
                pl, flag = _p_l(r)
                print(f"{c.name} {d} {r['p']} {r['herald_failures']} {pl:.4g}{'*' if flag else ''}")

    projections: dict[str, dict[str, ConfigProjection]] = {}
    for pair_name, pair in PAIRS.items():
        projections[pair_name] = {
            c.name: project_config(c, pilots[c.name], profile, pair) for c in CONFIGS
        }
        print(f"\n## 6c projection, d-pair {pair_name} (+ d = {CROSS_CHECK_D})")
        print(
            "config tasks shots build_h shots_h compile_h total_h peak_RAM_GiB "
            "held_tasks_GiB workers_GiB total_h_at_measured_efficiency"
        )
        for p in projections[pair_name].values():
            print(
                f"{p.name} {p.n_tasks} {p.total_shots} {p.build_h:.2f} {p.shots_h:.2f} "
                f"{p.compile_h:.2f} {p.total_h:.2f} {p.peak_ram_gb:.2f} "
                f"{p.held_gb:.2f} {p.workers_gb:.2f} {p.total_h_measured_eff:.2f}"
            )
        total = sum(p.total_h for p in projections[pair_name].values())
        total_eff = sum(p.total_h_measured_eff for p in projections[pair_name].values())
        print(
            f"TOTAL {total:.2f} h (at measured efficiency {total_eff:.2f} h); max peak RAM "
            f"{max(p.peak_ram_gb for p in projections[pair_name].values()):.2f} GiB "
            f"(limit {RAM_FRACTION_LIMIT * ram_gb:.2f} GiB)"
        )

    pair, deferred = decide(projections, ram_gb)
    print(f"\n## 6d decision: d-pair {pair}; deferred: {deferred or 'none'}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("profile").add_argument("--out", type=Path, required=True)
    pilot = sub.add_parser("pilot")
    pilot.add_argument("--config", choices=sorted(BY_NAME), required=True)
    pilot.add_argument("--out", type=Path, required=True)
    sub.add_parser("project").add_argument("--dir", type=Path, required=True)
    args = parser.parse_args()
    if args.cmd == "profile":
        args.out.mkdir(parents=True, exist_ok=True)
        cmd_profile(args.out)
    elif args.cmd == "pilot":
        args.out.mkdir(parents=True, exist_ok=True)
        cmd_pilot(BY_NAME[args.config], args.out)
    else:
        cmd_project(args.dir)


if __name__ == "__main__":
    main()
