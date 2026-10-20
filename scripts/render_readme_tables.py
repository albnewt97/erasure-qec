"""Render the README result tables from results/*.json.

Each table has a name; README.md carries it between ``<!-- BEGIN:name -->`` and
``<!-- END:name -->`` and tests/test_readme_tables_in_sync.py checks that the
block equals :func:`render` exactly. Tables and sources:

- ``thresholds`` (results/headline.json): per (r_e, decoder) the per-shot (7,11)
  pairwise crossing with its 95% CI and bracket, n_sign_changes_unsaturated, the
  per-shot d >= 7 collapse fit, the gap (collapse - pairwise) and the per-round
  d >= 7 collapse fit (secondary);
- ``herald_vs_blind`` (results/headline.json): herald - blind difference of the
  per-shot (7,11) crossings with its 95% CI and the replicates in which either
  crossing is missing or saturated;
- ``ablation`` (results/ablation.json): each decoder's error count, p_L(blind) /
  p_L(herald) per (p, r_e, d) and each decoder's per-round Lambda_d;
- ``wu_validation`` (results/wu_validation.json): the herald (9,13) per-shot
  crossing per config against Wu et al.'s reported value.

Rates, CI endpoints, brackets and differences print as percentages (or
percentage points) to 2 decimals; the collapse-minus-pairwise gap prints to 3
decimals in pp. Every table ends with one footnote line, chosen per table
(``TABLES``). The tables with bootstrap CIs (thresholds, herald_vs_blind,
wu_validation) end with :data:`FOOTNOTE`: CI endpoints can move in the second
decimal with the NumPy version. The ablation table has no bootstrap: its bounds
are closed-form Wilson limits and its error counts come from stim sampling, so
it ends with :data:`STIM_FOOTNOTE` instead. A table whose JSON does not exist is
not rendered.

    uv run python scripts/render_readme_tables.py            # every table
    uv run python scripts/render_readme_tables.py thresholds # one table
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FOOTNOTE = "CI endpoints can move in the second decimal with the NumPy version."
STIM_FOOTNOTE = (
    "stim's seeded sampling reproduces these counts only on the same machine and stim version."
)


def _pct(x: float | None) -> str:
    """A rate as a percentage, 2 decimals; a value that rounds to zero prints as 0.00."""
    if x is None:
        return "—"
    text = f"{100.0 * x:.2f}"
    return "0.00" if text == "-0.00" else text


def _pp(x: float | None, digits: int = 2) -> str:
    """A value already in percentage points."""
    if x is None:
        return "—"
    text = f"{x:.{digits}f}"
    return text.replace("-", "", 1) if float(text) == 0.0 else text


def _ci(ci: list[float | None] | None) -> str:
    return "—" if not ci or None in ci else f"[{_pct(ci[0])}, {_pct(ci[1])}]"


def _bracket(br: list[float] | None) -> str:
    return "—" if br is None else f"[{_pct(br[0])}, {_pct(br[1])}]"


def _crossing(c: dict[str, Any]) -> str:
    return str(c["status"]) if c["status"] != "crossing" else f"{_pct(c['value'])} {_ci(c['ci95'])}"


def _fit(f: dict[str, Any]) -> str:
    text = f"{_pct(f['p_th'])} {_ci(f['p_th_ci95'])}"
    return text if f["resolved"] else f"not resolved ({f['reason']})"


def _thresholds(doc: dict[str, Any]) -> list[str]:
    d1, d2 = doc["parameters"]["headline_pair"]
    d_min = doc["parameters"]["d_min"]
    lines = [
        f"| `r_e` | decoder | per-shot ({d1},{d2}) crossing, % [95% CI] | bracket, % "
        f"| sign changes (unsaturated) | per-shot collapse `d ≥ {d_min}`, % [95% CI] "
        f"| gap collapse − pairwise (pp) | per-round collapse `d ≥ {d_min}` (secondary), "
        "% [95% CI] |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for c in doc["conditions"]:
        h = c["headline"]
        lines.append(
            f"| {c['r_e']:g} | {c['decoder']} | {_crossing(h)} | {_bracket(h['bracket'])} "
            f"| {h['n_sign_changes_unsaturated']} | {_fit(c['collapse_shot_dmin7'])} "
            f"| {_pp(c['gap_collapse_minus_pairwise_pp'], 3)} "
            f"| {_fit(c['per_round_secondary'])} |"
        )
    return lines


def _herald_vs_blind(doc: dict[str, Any]) -> list[str]:
    d1, d2 = doc["parameters"]["headline_pair"]
    headline = {(c["r_e"], c["decoder"]): c["headline"] for c in doc["conditions"]}
    lines = [
        f"| `r_e` | herald ({d1},{d2}), % | blind ({d1},{d2}), % "
        "| herald − blind (pp) [95% CI] | replicates missing or saturated | note |",
        "|---|---|---|---|---|---|",
    ]
    for diff in doc["herald_minus_blind"]:
        r_e = diff["r_e"]
        value = None if diff["value"] is None else 100.0 * diff["value"]
        lo, hi = diff["ci95"]
        ci = "—" if lo is None or hi is None else f"[{_pp(100.0 * lo)}, {_pp(100.0 * hi)}]"
        note = "control: no heralds" if diff["control"] else ""
        lines.append(
            f"| {r_e:g} | {_crossing(headline[(r_e, 'herald_mwpm')])} "
            f"| {_crossing(headline[(r_e, 'blind_mwpm')])} | {_pp(value)} {ci} "
            f"| {diff['n_missing_or_saturated']} / {diff['n_boot']} | {note} |"
        )
    return lines


def _ratio(cell: dict[str, Any]) -> str:
    if cell["gated"]:
        return "—" if cell["ratio_lb"] is None else f"> {cell['ratio_lb']:.0f}× †"
    return "—" if cell["ratio"] is None else f"{cell['ratio']:.1f}×"


def _lam(x: float | None) -> str:
    return "—" if x is None else f"{x:.2f}"


def _ablation(doc: dict[str, Any]) -> list[str]:
    params = doc["parameters"]
    lines = [
        "| `p`, % | `r_e` | `d` | herald errors | blind errors | p_L(blind) / p_L(herald) "
        "| herald Λ_d | blind Λ_d |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for c in sorted(doc["cells"], key=lambda c: (c["p"], c["r_e"], c["d"])):
        lines.append(
            f"| {_pct(c['p'])} | {c['r_e']:g} | {c['d']} | {c['h_err']:,} | {c['b_err']:,} "
            f"| {_ratio(c)} | {_lam(c['lambda_herald'])} | {_lam(c['lambda_blind'])} |"
        )
    lines += [
        "",
        f"{params['shots']:,} shots per cell, seed {params['seed']}; error counts are out of "
        f"{params['shots']:,} shots.",
        f"† Fewer than {params['min_errors']} herald errors: a ≥ 95% one-sided lower confidence "
        "bound, not a point ratio.",
        f"Λ_d = per-round p_L(d) / p_L(d+2); — means fewer than {params['min_errors']} errors "
        "at d or d+2.",
    ]
    return lines


def _wu_validation(doc: dict[str, Any]) -> list[str]:
    lines = [
        "| config | `R_e` | `p_m` | herald (9,13), % [95% CI] | sign changes (unsaturated) "
        "| Wu et al. | ours − Wu (pp) | Wu inside CI |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for name, entry in doc["configs"].items():
        row = next(r for r in entry["rows"] if r["decoder"] == "wu_herald_mwpm")
        c = row["crossings"]["9,13"]
        inside = row["wu_inside_ci"]
        lines.append(
            f"| {name} | {entry['R_e']:g} | {entry['p_m']} | {_crossing(c)} "
            f"| {c['n_sign_changes_unsaturated']} | {row['wu']['text']} "
            f"| {_pp(row['ours_minus_wu_pp'])} | {'—' if inside is None else inside} |"
        )
    return lines


# name -> (source JSON, renderer, closing footnote line)
TABLES: dict[str, tuple[str, Callable[[dict[str, Any]], list[str]], str]] = {
    "thresholds": ("headline.json", _thresholds, FOOTNOTE),
    "herald_vs_blind": ("headline.json", _herald_vs_blind, FOOTNOTE),
    "ablation": ("ablation.json", _ablation, STIM_FOOTNOTE),
    "wu_validation": ("wu_validation.json", _wu_validation, FOOTNOTE),
}


def footnote(name: str) -> str:
    """The closing footnote line of table ``name``."""
    return TABLES[name][2]


def source(name: str, results_dir: Path = RESULTS) -> Path:
    """The JSON file table ``name`` is rendered from."""
    return results_dir / TABLES[name][0]


def render(name: str, results_dir: Path = RESULTS) -> str | None:
    """Markdown of table ``name`` (no trailing newline), or None if its JSON is missing."""
    path = source(name, results_dir)
    if not path.exists():
        return None
    _, table, footnote = TABLES[name]
    return "\n".join([*table(json.loads(path.read_text())), "", footnote])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("names", nargs="*", help=f"tables to print (default: {', '.join(TABLES)})")
    args = parser.parse_args()
    unknown = sorted(set(args.names) - set(TABLES))
    if unknown:
        raise SystemExit(f"unknown table(s) {unknown}; choose from {list(TABLES)}")
    for name in args.names or list(TABLES):
        table = render(name)
        if table is None:
            print(f"{name}: not rendered, {source(name).relative_to(ROOT)} is missing\n")
            continue
        print(f"<!-- BEGIN:{name} -->\n{table}\n<!-- END:{name} -->\n")


if __name__ == "__main__":
    main()
