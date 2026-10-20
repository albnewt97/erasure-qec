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
