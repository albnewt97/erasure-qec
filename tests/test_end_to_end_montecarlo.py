"""End-to-end Monte-Carlo smoke + throughput (PLAN.md §9/§11, M7 gate).

Both tests are @slow (excluded from the CI fast lane; run locally with plain
``uv run pytest``).

Throughput note (M7): the pre-optimization slow path rebuilt a matcher per
shot via Python add_edge calls and measured ~69x slower at r_e=0.98 than the
r_e=0 fast path — over the 50x trigger — so the herald-signature grouping
optimization was implemented (one vectorized-from-check-matrix reweighted
matcher per distinct fired-herald signature, batched decode per group). At
d=5, p=3e-2, r_e=0.98 nearly every signature is unique (E[fired heralds]
~ 12 of 800), so per-signature matcher construction remains the floor; the
measurement below prints the achieved rates.

Outcome, measured 2026-09-06 on an Apple M2 (8 cores, macOS 14.5): the
optimization took the gap from ~69x to ~51x. It did NOT clear the 50x
trigger. The 48x once quoted in README.md was 80,000/1,650 — a ratio of two
rounded numbers, not a measured one; the measured paired ratio is 48.6-54.6x
over 13 idle runs (median 51.3).
"""

import time
from pathlib import Path

import pytest
import sinter

from erasure_qec.circuits.builder import build
from erasure_qec.config import NoiseParams
from erasure_qec.decoding.sinter_adapter import CUSTOM_DECODERS, contract_dem
from erasure_qec.noise.injector import ErasureInjector


@pytest.mark.slow
def test_throughput_fast_vs_heralded_slow_path() -> None:
    """Time decode_shots_bit_packed on 10^4 shots at d=5, p=3e-2 for r_e=0
    (all fast path) and r_e=0.98 (essentially all slow path); print shots/sec."""
    n = 10_000
    rates: dict[float, float] = {}
    for r_e in (0.0, 0.98):
        circuit = build(5, 5, ErasureInjector(NoiseParams(p=3e-2, r_e=r_e)))
        compiled = CUSTOM_DECODERS["herald_mwpm"].compile_decoder_for_dem(dem=contract_dem(circuit))
        packed_dets, _ = circuit.compile_detector_sampler(seed=1).sample(
            n, separate_observables=True, bit_packed=True
        )
        start = time.perf_counter()
        preds = compiled.decode_shots_bit_packed(bit_packed_detection_event_data=packed_dets)
        elapsed = time.perf_counter() - start
        assert preds.shape[0] == n
        rates[r_e] = n / elapsed
        print(f"\nthroughput d=5 p=3e-2 r_e={r_e}: {rates[r_e]:,.0f} shots/sec")
    ratio = rates[0.0] / rates[0.98]
    print(f"fast/slow ratio: {ratio:.1f}x")

    # Gate: the slow path must not be CATASTROPHICALLY worse than the fast path.
    # The bound is deliberately loose, because this benchmark cannot separate the
    # regression one would most like to catch from ordinary environment noise.
    # Measured on an Apple M2 (8 cores, macOS 14.5), n=10^4, 2026-09-06:
    #
    #   current impl, idle machine .................. 48.6-54.6x (13 runs, median 51.3)
    #   current impl, 4 busy cores .................. 49.7x
    #   current impl, 8 busy cores (oversubscribed) . 55.8, 62.0, 65.3x
    #   reverted to per-shot add_edge rebuild, idle . 68.9x
    #
    # The per-shot-rebuild regression (68.9x) sits only 1.35x above the median idle
    # ratio, while CPU contention alone pushes the healthy impl to 65.3x -- the slow
    # path is Python-bound matcher construction and the fast path is a single
    # batched C++ call, so the two degrade at different rates under load. That
    # leaves [65.3, 68.9] for a bound that both catches the regression and never
    # flakes: 5% wide, against idle-run scatter that is already +/-6%. There is no
    # safe choice. Catching a revert to per-shot rebuilds needs a same-run A/B
    # against a shadow implementation, not a threshold on this ratio.
    #
    # So 100x gates catastrophe only: 1.5x above the worst ratio seen under full
    # contention, so it will not flake on slow or noisy CI, and it fires if the gap
    # roughly doubles from today's. A bare revert to per-shot add_edge rebuilds
    # would NOT trip it -- that limitation is the measurement's, not the threshold's.
    #
    # That gap is now covered STRUCTURALLY instead, by
    # tests/test_herald_decoder.py::
    #   test_slow_path_builds_one_matcher_per_signature_not_per_shot,
    # which counts Matching.from_check_matrix calls. Counting is exact and
    # machine-independent, so it catches the per-shot-rebuild regression that no
    # threshold on this ratio can. A same-run A/B against a shadow slow
    # implementation was considered and rejected: it would require keeping a
    # deliberately-slow duplicate decoder in the tree forever, and the machine
    # speed it cancels is not actually what defeats this test -- the overlap
    # between "contended healthy" and "idle regressed" is.
    #
    # Re-measured 2026-09-06, 10 consecutive runs on an idle Apple M2:
    # ratio min 48.8x, median 50.8x, max 54.4x (sd 1.6).
    assert ratio < 100.0, (
        f"slow path {ratio:.1f}x slower than the fast path; reference is 48.6-54.6x "
        "on an idle Apple M2 and <=65.3x under full CPU contention"
    )


@pytest.mark.slow
def test_end_to_end_smoke_p_l_decreases_with_distance(tmp_path: Path) -> None:
    """M7 gate: 10^4-shot end-to-end run at d in {3,5} for herald_mwpm AND
    blind_mwpm, CSV written to a tmp dir, and p_L(d=5) < p_L(d=3) below threshold
    (p=1e-2, r_e=0.5; verified live: herald p_L ~ 0.049 vs ~ 0.026)."""
    n = 10_000
    tasks = []
    for d in (3, 5):
        circuit = build(d, d, ErasureInjector(NoiseParams(p=1e-2, r_e=0.5)))
        tasks.append(
            sinter.Task(
                circuit=circuit,
                detector_error_model=contract_dem(circuit),
                json_metadata={"d": d, "rounds": d, "p": 1e-2, "r_e": 0.5},
            )
        )

    csv_path = tmp_path / "smoke_e2e.csv"  # never touch the real data/ dir
    stats = sinter.collect(
        num_workers=2,
        tasks=tasks,
        decoders=["herald_mwpm", "blind_mwpm"],
        custom_decoders=CUSTOM_DECODERS,
        max_shots=n,
        max_errors=n,  # never binds: the full 10^4 shots are collected
        save_resume_filepath=csv_path,
    )
    assert csv_path.exists()

    p_l: dict[tuple[str, int], float] = {}
    for stat in stats:
        assert stat.shots >= n, f"{stat.decoder} d={stat.json_metadata['d']} incomplete"
        assert stat.errors > 0  # p=1e-2 is comfortably above the noise floor
        p_l[(stat.decoder, stat.json_metadata["d"])] = stat.errors / stat.shots
    assert set(p_l) == {(dec, d) for dec in ("herald_mwpm", "blind_mwpm") for d in (3, 5)}

    # Below threshold, distance must help; herald must also beat blind at fixed d.
    assert p_l[("herald_mwpm", 5)] < p_l[("herald_mwpm", 3)]
    assert p_l[("herald_mwpm", 3)] < p_l[("blind_mwpm", 3)]
    assert p_l[("herald_mwpm", 5)] < p_l[("blind_mwpm", 5)]
    print(f"\nsmoke p_L: {p_l}")
