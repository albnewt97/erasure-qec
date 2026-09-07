# Related work and future directions

Directions I would take this next, and why each one interests me, with the paper that motivates it. I
have not built on any of these yet; each citation was verified against its arXiv/DOI page.

1. Zhang, B. et al. **Logical qubits with erasure conversion using metastable neutral atoms.** Nat.
   Phys. **22**, 910 (2026). [arXiv:2506.13724](https://arxiv.org/abs/2506.13724) — posted under the
   preprint title *Leveraging erasure errors in logical qubits with metastable ¹⁷¹Yb atoms* and
   retitled at v2. This is the experimental realization of exactly the setting I simulate, which
   makes it the natural first comparison: I would benchmark my simulated sub-threshold suppression
   against their measured logical-qubit performance.

2. Chow, M. N. H. et al. **Circuit-based leakage-to-erasure conversion in a neutral-atom quantum
   processor.** PRX Quantum **5**, 040343 (2024).
   [arXiv:2405.10434](https://arxiv.org/abs/2405.10434). My model assumes an ideal in-gate herald;
   their conversion is circuit-based and imperfect, so modeling imperfect or delayed erasure conversion
   would make the noise model more faithful to what the hardware actually does.

3. Scholl, P. et al. **Erasure conversion in a high-fidelity Rydberg quantum simulator.** Nature
   **622**, 273 (2023). [arXiv:2305.03406](https://arxiv.org/abs/2305.03406). A parallel
   experimental line that would inform the realistic herald fidelities I currently take as ideal.

4. Kang, M., Campbell, W. C. & Brown, K. R. **Quantum error correction with metastable states of
   trapped ions using erasure conversion.** PRX Quantum **4**, 020358 (2023).
   [arXiv:2210.15024](https://arxiv.org/abs/2210.15024). A third hardware profile — trapped ions —
   worth adding as an experiment config alongside the neutral-atom and dual-rail cases.

5. Gu, S., Retzker, A. & Kubica, A. **Fault-tolerant quantum architectures based on erasure qubits.**
   Phys. Rev. Res. **7**, 013249 (2025). [arXiv:2312.14060](https://arxiv.org/abs/2312.14060). — and
   Gu, S., Vaknin, Y., Retzker, A. & Kubica, A. **Optimizing quantum error correction protocols with
   erasure qubits.** PRX Quantum **6**, 040354 (2025).
   [arXiv:2408.00829](https://arxiv.org/abs/2408.00829). Between them these motivate
   extending the noise model to imperfect erasure checks and false heralds, which my ideal-herald model
   omits entirely.

6. Baranes, G. et al. **Leveraging qubit loss detection in fault-tolerant quantum algorithms.** Phys.
   Rev. X **16**, 011002 (2026). [arXiv:2502.20558](https://arxiv.org/abs/2502.20558). This points past
   the memory experiment I built toward logical algorithms, where loss/erasure information can be
   carried through a whole computation rather than a single memory round.

7. Yu, C.-C. et al. **Taming Rydberg decay with measurement-based quantum computation.** Phys. Rev.
   Lett. **136**, 160601 (2026). [arXiv:2411.04664](https://arxiv.org/abs/2411.04664). A
   measurement-based alternative to the matching-based recovery I use, and an interesting contrast in
   how the same herald information gets exploited.

8. Wu, Y. & Zhong, L. **Fusion Blossom: Fast MWPM Decoders for QEC.** In *2023 IEEE International
   Conference on Quantum Computing and Engineering (QCE)*, 928–938 (2023).
   [arXiv:2305.08307](https://arxiv.org/abs/2305.08307). An alternative MWPM backend I would use to put
   a throughput number on my slow path relative to a decoder built for speed.

9. Delfosse, N. & Zémor, G. **Linear-time maximum likelihood decoding of surface codes over the
   quantum erasure channel.** Phys. Rev. Research **2**, 033042 (2020).
   [arXiv:1703.01517](https://arxiv.org/abs/1703.01517). — and Delfosse, N. & Nickerson, N. H.
   **Almost-linear time decoding algorithm for topological codes.** Quantum **5**, 595 (2021).
   [arXiv:1709.06218](https://arxiv.org/abs/1709.06218). **The largest known decoder gap in this repo.**
   I decode with MWPM and set the heralded edges' weight to 0; peeling / weighted union-find is the
   *maximum-likelihood* decoder for pure erasure, which weight-0 MWPM is not. The gap grows with `r_e`,
   so it bites hardest at exactly the `r_e = 0.98` regime my strongest claim (the ablation) lives in —
   my numbers there should be read as a lower bound on what an ML erasure decoder would deliver. This is
   also the decoder Wu et al. use, so it is a prerequisite for any genuinely like-for-like comparison to
   their thresholds; the README's "Comparison to the literature" now flags the difference and points
   here.

## Methodology

**Paired herald-vs-blind threshold sweeps.** The threshold-difference test
(`bootstrap_threshold_difference`, docs/AUDIT.md) currently uses an *unpaired*
bootstrap because `sinter` samples each decoder independently and the CSV stores
only marginal `(shots, errors)`. Making it genuinely paired is small and
concrete — it needs neither re-thinking nor a large re-collection, just a
four-integer-per-`(p, d)` schema instead of two. The pattern already exists in
`scripts/ablation_table.py`, which samples once and decodes the same `dets`
array with both decoders; the sweep collector would do the same and record, per
`(p, d)`, the 2×2 joint outcome counts:

- `n00` — both decoders correct,
- `n11` — both wrong,
- `n10` — herald wrong, blind correct,
- `n01` — herald correct, blind wrong.

That is exactly the input a paired bootstrap (multinomial-resample the four
cells, derive each decoder's marginal errors) and a McNemar test (`n10` vs
`n01`) need. `bootstrap_threshold_difference` already accepts this via its
`joint_counts` argument (exercised by the synthetic paired tests); only the
collection path in `experiments/` needs the schema change. Since the unpaired CI
is conservative, this can only tighten the interval — it would not change the
`r_e = 0.5` verdict, only sharpen it.

**Zero-denominator replicates in the Lambda bootstrap.** `lambda_factor` builds
its CI from `Λ = p_L(d) / p_L(d+2)` replicates and must do *something* with the
draws where the `d+2` code logs zero errors and the denominator vanishes. It
currently excludes them, now counted as `n_boot_dropped` and flagged by
`LambdaEstimate.ci_conditioned` (docs/AUDIT.md, "Conditioning in the Lambda
bootstrap"). Exclusion is the conservative-looking choice but it is not the
honest one: those replicates have `Λ = +∞`, so removing them pulls the upper
endpoint down, and where they are more than 2.5% of the draws the true 97.5th
percentile is unbounded. The better estimator propagates `+∞` (giving
`high = inf` once the drop fraction passes the tail mass) or reports a one-sided
lower bound on `Λ` in that regime, the way `scripts/ablation_table.py` already
does for its sub-gate cells. This was left undone deliberately: it changes a
published quantity rather than fixing bookkeeping, and the `≥ 50`-error gate
means no currently plotted point is affected (all 425 gated-in pairs drop zero
replicates). The work is to change the estimator, re-check the gate is still the
right cut, and re-state any `Λ` quoted below the gate as a bound.

**⚠ The bootstrap pin asymmetry (read before touching either bootstrap).** `_bootstrap_pipeline`
retains bound-pinned replicates and `bootstrap_threshold_difference` discards them; that
inconsistency is deliberate and load-bearing, not a bug — making the two agree moves the r_e = 0.98
blind `rel_ci_width` from 2.43–4.26 to 0.39–1.36 (astride the 1.0 guard) and flips `resolved` at 9 of
30 seeds, but **not at seed 0**, so the seed this repo declares will show the change as inert. See
docs/AUDIT.md, "The bootstrap pin asymmetry is load-bearing (2026-09-07)".
