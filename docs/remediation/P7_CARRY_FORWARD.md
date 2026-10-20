# P7 carry-forward ledger

Items recorded in P6 for P7 (and P8) to apply. One item per line; the bracket names the
destination document(s).

1. [README, METHODS] Per-shot failure after T = d rounds is the headline metric; per-round is secondary. Do not call per-round "biased" or "wrong". At T = d the per-round crossing never lies below the per-shot crossing (test test_round_negative_wherever_shot_nonpositive); which is closer to the asymptotic threshold is an empirical question.
2. [CHANGELOG] Per-shot crossings computed on the pre-dense grid had CIs that excluded interpolation error. Densification moved r50 blind from 1.181% to 1.104%, outside its earlier 95% CI. Every pre-dense per-shot number (including "about 3.9%") is superseded.
3. [README, CHANGELOG] The README dagger footnote on the r50 blind CI upper endpoint (2.154) is withdrawn; on the dense data it is 1.800.
4. [README] r50 blind: the point estimate sits at the lower edge of its bracket and the CI extends below it; (9,11) gives 1.028%. Never quote 1.104% without its CI.
5. [README] Report the per-shot collapse fit next to each pairwise crossing, with the gap (collapse - pairwise) as a systematic column. Collapse is lower in all six conditions (dense data: -0.004 to -0.025 pp).
6. [README, VALIDATION] Do not present the baseline per-shot threshold as agreeing with Wu's 0.937%: the repo model has measurement, reset and idle noise that Wu's does not.
7. [VALIDATION] wu_si_pm005 (9,13) has n_sign_changes_unsaturated = 3; state it next to its value. Same for the (7,9) cross-checks of re098 and pmp in the shift table.
8. [METHODS] sinter pickles circuits and DEMs via text with about 6 significant digits; sampled probabilities differ from nominal by about 1e-6 relative. No result depends on this.
9. [CHANGELOG] Fixed a latent bug: at r_e = 1 with no unheralded noise, partition_flattened_dem declared 0 observables and HeraldMatchingDecoder crashed. No committed sweep reached that case.
10. [P8] .gitignore: "figures/dev/   # comment" ignores nothing (git has no inline comments).
11. [P8] scripts/ is outside mypy's scope.
12. [VALIDATION] The saturated-crossing rule was added after data; it changed no herald value (already documented; keep it).
13. [README] The Wu validation used d-pair (9,13), set by a memory rule; Wu used (11,15).
14. [README, P8] Legacy scripts (paired_separation.py and any other script whose recorded output was produced on the pre-dense CSVs) reproduce their baseline_outputs only at commit da47e3a. Move them to scripts/legacy/ with a header naming that commit, or drop them from the README reproduce table. They must not be listed as reproducing current results.
15. [README, METHODS] At r_e = 0.98 the per-shot collapse fit (d >= 7) resolves both thresholds (herald 3.83%, blind 1.18%); the per-round fit resolves neither (herald bound_pin, blind unresolved_ci). State this as an observation with no claim about which metric is correct (see item 1).
16. [README, METHODS] Lambda is reported per-round while thresholds are per-shot; label the metric at every appearance.

## Status (P7 reconciliation, gate G6)

Every item is `done` (file and section that carry it at the end of P7) or `deferred to P8`. Filled
says whether the destination already carries the item at this commit or receives it in P7 session 2.

| item | status | destination | filled |
|---|---|---|---|
| 1 | done | README.md#metric; docs/METHODS.md#ordering-of-per-shot-and-per-round-crossings | session 2 |
| 2 | done | CHANGELOG.md#fixed | session 2 |
| 3 | done | README.md#thresholds; CHANGELOG.md#withdrawn-claims | session 2 |
| 4 | done | README.md#thresholds | session 2 |
| 5 | done | README.md#thresholds (collapse and gap columns of the rendered `thresholds` table) | yes |
| 6 | done | README.md#validation-against-wu-et-al; README.md#model | session 2 |
| 7 | done | docs/VALIDATION.md#results; docs/VALIDATION.md#differences-from-wu-et-al (shift table) | session 2 |
| 8 | done | docs/METHODS.md#sinter-text-pickling | session 2 |
| 9 | done | CHANGELOG.md#fixed | session 2 |
| 10 | deferred to P8 | — | — |
| 11 | deferred to P8 | — | — |
| 12 | done | docs/VALIDATION.md#analysis-rules-added-after-data | yes |
| 13 | done | docs/VALIDATION.md#differences-from-wu-et-al (item 1); README.md#validation-against-wu-et-al | VALIDATION yes; README session 2 |
| 14 | done | README.md#reproduce (legacy scripts listed as reproducing their baseline_outputs only at da47e3a, not as reproducing current results) | session 2 |
| 15 | done | README.md#thresholds; docs/METHODS.md#collapse-fit | session 2 |
| 16 | done | README.md#thresholds; README.md#ablation; docs/METHODS.md#lambda | session 2 |
