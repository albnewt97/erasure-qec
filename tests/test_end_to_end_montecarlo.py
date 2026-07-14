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
"""

import time
from pathlib import Path

import pytest
import sinter

from erasure_qec.circuits.builder import build
from erasure_qec.config import NoiseParams
from erasure_qec.decoding.sinter_adapter import CUSTOM_DECODERS, contract_dem
from erasure_qec.noise.injector import BiasedErasureInjector

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


@pytest.mark.slow
def test_throughput_fast_vs_heralded_slow_path() -> None:
    """Time decode_shots_bit_packed on 10^4 shots at d=5, p=3e-2 for r_e=0
    (all fast path) and r_e=0.98 (essentially all slow path); print shots/sec."""
    n = 10_000
    rates: dict[float, float] = {}
    for r_e in (0.0, 0.98):
        circuit = build(5, 5, BiasedErasureInjector(NoiseParams(p=3e-2, r_e=r_e)))
        compiled = CUSTOM_DECODERS["herald_mwpm"].compile_decoder_for_dem(
            dem=contract_dem(circuit)
        )
        packed_dets, _ = circuit.compile_detector_sampler(seed=1).sample(
            n, separate_observables=True, bit_packed=True
        )
        start = time.perf_counter()
        preds = compiled.decode_shots_bit_packed(
            bit_packed_detection_event_data=packed_dets
        )
        elapsed = time.perf_counter() - start
        assert preds.shape[0] == n
        rates[r_e] = n / elapsed
        print(f"\nthroughput d=5 p=3e-2 r_e={r_e}: {rates[r_e]:,.0f} shots/sec")
    print(f"fast/slow ratio: {rates[0.0] / rates[0.98]:.1f}x")


@pytest.mark.slow
def test_end_to_end_smoke_p_l_decreases_with_distance() -> None:
    """M7 gate: 10^4-shot end-to-end run at d in {3,5} for herald_mwpm AND
    blind_mwpm, CSV written to data/, and p_L(d=5) < p_L(d=3) below threshold
    (p=1e-2, r_e=0.5; verified live: herald p_L ~ 0.049 vs ~ 0.026)."""
    n = 10_000
    tasks = []
    for d in (3, 5):
        circuit = build(d, d, BiasedErasureInjector(NoiseParams(p=1e-2, r_e=0.5)))
        tasks.append(
            sinter.Task(
                circuit=circuit,
                detector_error_model=contract_dem(circuit),
                json_metadata={"d": d, "rounds": d, "p": 1e-2, "r_e": 0.5},
            )
        )

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = DATA_DIR / "smoke_e2e.csv"
    csv_path.unlink(missing_ok=True)  # stale results would confuse resume
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
