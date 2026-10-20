# erasure-qec v0.2.0 — release notes (draft)

Released 2026-10-20. The full list of changes, with claim and caveat IDs, is in [CHANGELOG.md](../../CHANGELOG.md).

v0.2.0 re-derives every published number from committed scripts and `results/*.json`. It changes the headline metric, withdraws several v0.1 claims, and adds a validation against Wu et al. (Nat. Commun. 13, 4657, 2022) under their noise model.

## Headline result

At `r_e = 0.98` the per-shot `d = 7` and `d = 11` curves cross at 3.84% [3.60, 4.00] for the herald-conditioned decoder and at 1.19% [1.15, 1.23] for the blind decoder (`scripts/headline_results.py`, `results/headline.json`).
The intervals are bootstrap 95% CIs that cover binomial noise only, not interpolation error inside the grid bracket.
Their endpoints can move in the second decimal with the NumPy version.
Every threshold is a Z-basis memory threshold of the rotated surface code; see [README.md](../../README.md#results) for all conditions and caveats.

## What changed

- **Headline metric.** The headline is now the per-shot (7,11) pairwise crossing, `P_L` after `T = d` rounds. The v0.1 headline, the per-round `d ≥ 7` collapse fit, is now a secondary column.
- **Dense sweeps.** The thresholds and the herald − blind difference are computed on sweeps densified around the per-shot crossings by a rule fixed before the dense data were collected.
- **Collapse fits.** At `r_e = 0.98` the per-shot collapse fit resolves both decoders, and the per-round fit resolves neither. The per-round blind verdict depends on retaining bound-pinned bootstrap replicates: excluding them resolves that fit at 21 of 30 bootstrap seeds (`results/per_round_fit_sensitivity.json`).
- **Ablation.** The paired herald-versus-blind ablation now covers `p = 0.5%` and `p = 1.0%`. Its error counts reproduce only on the same machine and stim version.
- **Wu et al. validation.** A `wu2022` noise convention follows Wu's noise model, and the comparison uses the herald (9,13) per-shot crossing ([docs/VALIDATION.md](../VALIDATION.md)). Wu's value lies inside our 95% CI for `R_e = 1` only. The differences from Wu's setup, including the distance pair and the decoder, are listed there; the effect of none of them was measured.

## Withdrawn claims

The following v0.1 statements are withdrawn; [CHANGELOG.md](../../CHANGELOG.md#withdrawn-claims) gives the reason for each.

- "blind barely moves 1.44% -> 1.49%", and "Delta = -0.83%" with its tipping-point analysis. Both were per-round fits on the pre-dense sweeps.
- The `r_e = 0.5` dagger footnote, which rested on a probe that was never committed.
- "the r_e = 0.98 threshold is not resolved", which holds for the per-round fit only.
- The comparison of the repo model with Wu's thresholds, and "R_e defined exactly as our r_e".
- "models dual-rail superconducting cavities", "each assumption is systematically optimistic", throughput numbers and timings with no committed record.

## Added

- `scripts/reproduce.sh` with three tiers: `fast` (tests, README tables, figures), `analysis` (every analysis script in a throwaway clone, compared with the committed results) and `collect` (prints the sampling commands with recorded wall times).
- `scripts/per_round_fit_sensitivity.py`: cap sensitivity and pin retention of the per-round fit.
- `CITATION.cff` and package metadata.

## Fixed

- At `r_e = 1` with no unheralded noise the herald decoder crashed on a DEM with no declared observables; no committed sweep reached that case.
- Docstrings and a config called the erasure model "biased"; it is unbiased.

## Reproducing

Install with `uv sync --group dev`, then run `scripts/reproduce.sh fast` or `scripts/reproduce.sh analysis`.
Neither modifies the work tree.
The sweeps cannot be regenerated shot for shot, because `sinter.collect` takes no seed.
`scripts/reproduce.sh collect` prints the collection commands.

## Upgrading from v0.1

v0.1 numbers should not be quoted.
The remediation ledgers that map each v0.1 claim and caveat to its v0.2.0 status are in [docs/archive/remediation/](remediation/CAVEAT_LEDGER.md).
