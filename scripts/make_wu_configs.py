"""Write the seven P5a Wu collection configs, experiments/configs/wu/*.yaml.

Fixed grids are ``round(k * step, 4)`` for integer ``k`` in ``[k_lo, k_hi]``.
The R_e = 0.5 and 0.9 grids are ``round(c * f, 4)`` for f in ``GRID_FACTORS``,
with ``c`` the per-shot (5,7) crossing printed by scripts/wu_coarse_scan.py
(recorded in docs/remediation/baseline_outputs/wu_coarse_scan.txt); the script
asserts that those two grids equal the ones printed there. Asserts the point
count of every grid. A config that the edge rule extended (see
docs/remediation/WU_DPAIR_DECISION.md) also gets ``extra_p_edge``, the
3 points beyond the touched end; the collection loader merges it into the grid.
Run:

    uv run python scripts/make_wu_configs.py
"""

from __future__ import annotations

import ast
import re
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "experiments" / "configs" / "wu"
SCAN_OUTPUT = ROOT / "docs" / "remediation" / "baseline_outputs" / "wu_coarse_scan.txt"
DECISION = "docs/remediation/WU_DPAIR_DECISION.md"

DISTANCES = (7, 9, 13)
MAX_SHOTS = 100_000
MAX_ERRORS = 1000
HERALD = ("wu_herald_mwpm",)
BOTH = ("wu_herald_mwpm", "wu_blind_mwpm")
GRID_FACTORS = (0.800, 0.845, 0.890, 0.935, 0.980, 1.025, 1.070, 1.115, 1.160, 1.205, 1.250)
# Edge-rule continuation of GRID_FACTORS (WU_DPAIR_DECISION.md, 2026-10-16).
EDGE_UPPER_FACTORS = (1.295, 1.340, 1.385)
EDGE_LOWER_FACTORS = (0.755, 0.710, 0.665)

# Per-shot (5,7) crossings from wu_coarse_scan.txt (repr as printed there).
C_RE050 = 0.02194620108673982
C_RE090 = 0.032103278742631816

EDGE_RULE = (
    "Edge rule (recorded in P5a, applied in P5b): if the (9,13) crossing's bracket\n"
    "touches either end of this grid, append 3 points beyond that end at the same\n"
    "step and collect this config once more."
)


@dataclass(frozen=True)
class WuConfigSpec:
    """One config file: noise, grid, decoders and the header text."""

    name: str
    r_e: float
    p_m: float | None  # None: p_m follows p
    grid: tuple[float, ...]
    decoders: tuple[str, ...]
    reference: str
    notes: tuple[str, ...] = ()
    edge: tuple[float, ...] = ()  # extra_p_edge; empty unless the edge rule applied


def stepped(k_lo: int, k_hi: int, step: float) -> tuple[float, ...]:
    """``round(k * step, 4)`` for integer ``k`` from ``k_lo`` to ``k_hi`` inclusive."""
    return tuple(round(k * step, 4) for k in range(k_lo, k_hi + 1))


def scaled(c: float) -> tuple[float, ...]:
    """``round(c * f, 4)`` for f in ``GRID_FACTORS``."""
    return tuple(round(c * f, 4) for f in GRID_FACTORS)


def scaled_edge(c: float, factors: tuple[float, ...]) -> tuple[float, ...]:
    """``round(c * f, 4)`` for f in ``factors`` (an edge-rule extension)."""
    return tuple(round(c * f, 4) for f in factors)


def _scan_grids() -> dict[float, tuple[float, ...]]:
    """The grids printed in the coarse-scan output, keyed by R_e."""
    text = SCAN_OUTPUT.read_text(encoding="utf-8")
    found = re.findall(r"^R_e = ([0-9.]+):.*?^  grid = (\[.*?\])$", text, re.M | re.S)
    return {float(r_e): tuple(ast.literal_eval(grid)) for r_e, grid in found}


SCAN_NOTE = (
    "Grid: round(c * f, 4), f = 0.800 .. 1.250 step 0.045, c = {c!r}, the per-shot\n"
    "(5,7) crossing of the coarse scan (herald decoder, 3000 shots per point;\n"
    "docs/remediation/baseline_outputs/wu_coarse_scan.txt, scripts/wu_coarse_scan.py)."
)

SPECS: tuple[WuConfigSpec, ...] = (
    WuConfigSpec(
        "wu_fig3_re000",
        0.0,
        0.0,
        stepped(15, 25, 0.0005),
        HERALD,
        "Wu reference: 0.937(4)%, Fig. 3a (R_e = 0, p_m = 0).",
        (
            "Decoders: wu_herald_mwpm only. At R_e = 0 the herald decoder's predictions\n"
            "equal plain PyMatching on the same DEM (T3,\n"
            "tests/test_wu_decoder.py::test_t3_re_zero_herald_equals_plain_pymatching),\n"
            "so wu_blind_mwpm would decode identically and is omitted.",
        ),
    ),
    WuConfigSpec(
        "wu_fig3_re098",
        0.98,
        0.0,
        stepped(34, 50, 0.001),
        BOTH,
        "Wu reference: 4.15(2)%, Fig. 3a (R_e = 0.98, p_m = 0).",
    ),
    WuConfigSpec(
        "wu_fig3b_re100",
        1.0,
        0.0,
        stepped(42, 60, 0.001),
        BOTH,
        "Wu reference: 5.13%, text, Sec. III (R_e = 1, p_m = 0).",
    ),
    WuConfigSpec(
        "wu_fig3b_re050",
        0.5,
        0.0,
        scaled(C_RE050),
        BOTH,
        'Wu reference: Fig. 3b, "not tabulated" (R_e = 0.5, p_m = 0).',
        (
            SCAN_NOTE.format(c=C_RE050),
            "Edge rule applied (P5b, lower end): the (9,13) per-shot herald curves do not\n"
            "cross on p_values and ln P_L(13) - ln P_L(9) >= 0 at the lowest p.\n"
            "extra_p_edge = round(c * f, 4), f = 0.755, 0.710, 0.665.",
        ),
        scaled_edge(C_RE050, EDGE_LOWER_FACTORS),
    ),
    WuConfigSpec(
        "wu_fig3b_re090",
        0.9,
        0.0,
        scaled(C_RE090),
        BOTH,
        'Wu reference: Fig. 3b, "not tabulated" (R_e = 0.9, p_m = 0).',
        (SCAN_NOTE.format(c=C_RE090),),
    ),
    WuConfigSpec(
        "wu_si_pm005",
        0.98,
        0.005,
        stepped(30, 46, 0.001),
        BOTH,
        "Wu reference: 3.80(2)%, SI S6 (R_e = 0.98, p_m = 0.005).",
    ),
    WuConfigSpec(
        "wu_si_pmp",
        0.98,
        None,
        stepped(22, 36, 0.001),
        BOTH,
        "Wu reference: 2.85(1)%, SI S6 (R_e = 0.98, p_m = p).",
    ),
)

EXPECTED_COUNTS = {
    "wu_fig3_re000": 11,
    "wu_fig3_re098": 17,
    "wu_fig3b_re100": 19,
    "wu_fig3b_re050": 11,
    "wu_fig3b_re090": 11,
    "wu_si_pm005": 17,
    "wu_si_pmp": 15,
}


def _comment(text: str) -> list[str]:
    return [f"# {line}" if line else "#" for line in text.split("\n")]


def render(spec: WuConfigSpec) -> str:
    """The YAML text of one config, header comment first."""
    header = [
        f"Wu et al. (2022) validation config {spec.name} (convention wu2022).",
        spec.reference,
        f"d-pair (9,13) with d = 7 cross-check: {DECISION}.",
        *spec.notes,
        EDGE_RULE,
        "Generated by scripts/make_wu_configs.py; do not edit by hand.",
    ]
    lines = [line for block in header for line in _comment(block)]
    follows = spec.p_m is None
    lines += [
        f"name: {spec.name}",
        "convention: wu2022",
        f"R_e: {spec.r_e!r}",
        f"p_m: {'null' if follows else repr(spec.p_m)}",
        f"p_m_follows_p: {'true' if follows else 'false'}",
        f"distances: [{', '.join(str(d) for d in DISTANCES)}]",
        f"p_values: [{', '.join(repr(p) for p in spec.grid)}]",
        f"max_shots: {MAX_SHOTS}",
        f"max_errors: {MAX_ERRORS}",
        f"decoders: [{', '.join(spec.decoders)}]",
    ]
    if spec.edge:
        lines.append(f"extra_p_edge: [{', '.join(repr(p) for p in spec.edge)}]")
    return "\n".join(lines) + "\n"


def main() -> None:
    scan = _scan_grids()
    assert scan[0.5] == scaled(C_RE050), (scan[0.5], scaled(C_RE050))
    assert scan[0.9] == scaled(C_RE090), (scan[0.9], scaled(C_RE090))
    assert [s.name for s in SPECS] == list(EXPECTED_COUNTS)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for spec in SPECS:
        n = len(spec.grid)
        assert n == EXPECTED_COUNTS[spec.name], (spec.name, n)
        assert len(set(spec.grid)) == n and list(spec.grid) == sorted(spec.grid), spec.name
        assert len(spec.edge) in (0, 3) and not set(spec.edge) & set(spec.grid), spec.name
        path = OUT_DIR / f"{spec.name}.yaml"
        path.write_text(render(spec), encoding="utf-8")
        print(f"{spec.name}: {n} points {list(spec.grid)}")
        if spec.edge:
            print(f"{spec.name}: extra_p_edge {list(spec.edge)}")


if __name__ == "__main__":
    main()
