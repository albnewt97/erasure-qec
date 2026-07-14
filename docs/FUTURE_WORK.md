# Related work and future directions

Directions I would take this next, and why each one interests me, with the paper that motivates it. I
have not built on any of these yet; each citation was verified against its arXiv/DOI page.

1. Zhang, B. et al. **Logical qubits with erasure conversion using metastable neutral atoms.** Nat.
   Phys. **22**, 910 (2026). [arXiv:2506.13724](https://arxiv.org/abs/2506.13724). This is the
   experimental realization of exactly the setting I simulate, which makes it the natural first
   comparison: I would benchmark my simulated sub-threshold suppression against their measured
   logical-qubit performance.

2. Chow, M. N. H. et al. **Circuit-based leakage-to-erasure conversion in a neutral-atom quantum
   processor.** PRX Quantum **5**, 040343 (2024).
   [arXiv:2405.10434](https://arxiv.org/abs/2405.10434). My model assumes an ideal in-gate herald;
   their conversion is circuit-based and imperfect, so modeling imperfect or delayed erasure conversion
   would make the noise model more faithful to what the hardware actually does.

3. Scholl, P. et al. **Erasure conversion in a high-fidelity Rydberg quantum simulator.** Nature
   **622**, 273 (2023). A parallel experimental line that would inform the realistic herald fidelities
   I currently take as ideal.

4. Kang, M., Campbell, W. C. & Brown, K. R. **Quantum error correction with metastable states of
   trapped ions using erasure conversion.** PRX Quantum **4**, 020358 (2023).
   [arXiv:2210.15024](https://arxiv.org/abs/2210.15024). A third hardware profile — trapped ions —
   worth adding as an experiment config alongside the neutral-atom and dual-rail cases.

5. Gu, S., Retzker, A. & Kubica, A. **Fault-tolerant quantum architectures based on erasure qubits.**
   Phys. Rev. Res. **7**, 013249 (2025). [arXiv:2312.14060](https://arxiv.org/abs/2312.14060). — and
   Gu, S., Vaknin, Y., Retzker, A. & Kubica, A. **Optimizing quantum error correction protocols with
   erasure qubits.** [arXiv:2408.00829](https://arxiv.org/abs/2408.00829). Between them these motivate
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

8. Wu, Y. & Zhong, L. **Fusion Blossom: Fast MWPM Decoders for QEC.**
   [arXiv:2305.08307](https://arxiv.org/abs/2305.08307). An alternative MWPM backend I would use to put
   a throughput number on my slow path relative to a decoder built for speed.
