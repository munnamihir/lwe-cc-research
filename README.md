# LWE-CC: Communication Complexity of Lattice-Based Cryptographic Protocols

**Author:** Mihir Bommisetty — Independent Researcher, Mechanicsburg PA  
**Started:** April 15, 2026  
**Status:** Active research — work in progress

## What This Is

This repository contains the **first computational study of the communication complexity of the Learning With Errors (LWE) function** — the mathematical core of NIST post-quantum cryptographic standards ML-KEM (FIPS 203) and ML-DSA (FIPS 204).

## The Core Problem: LWE-CC_n

**Alice's input:** secret `s ∈ ℤ_q^n`, error `e ∈ ℤ_q^m` (small)  
**Bob's input:** matrix `A ∈ ℤ_q^{m×n}`  
**Goal:** Bob computes `b = As + e mod q` via one-way communication

**Question:** What is the minimum communication needed?

## Results (April 2026 — First in Literature)

| n | q | rank(M_f) | Fooling set | Classical LB | Quantum UB (conj.) |
|---|---|-----------|-------------|--------------|---------------------|
| 2 | 5 | 31        | ≥ 104       | ≥ 7 bits     | ~2 qubits           |
| 2 | 7 | 53        | ≥ 114       | ≥ 7 bits     | ~2 qubits           |
| 3 | 5 | 126       | ≥ 173       | ≥ 8 bits     | ~2 qubits           |

**These numbers do not appear in any prior paper.**

## Proved Theorems

**Theorem 1 (trivial, rigorous):** CC(LWE-CC_n) ≥ n · log₂(q) bits.  
*Proof:* Setting A = I_n, Bob must recover s from Alice's message. This requires log₂(|ℤ_q^n|) = n log₂ q bits.

**Theorem 2 (computational):** For (n=2,q=5), (n=2,q=7), (n=3,q=5): fooling sets of size ≥ 104, 114, 173 exist, giving CC lower bounds of 7, 7, 8 bits respectively.

## Conjectures

**Conjecture 1:** rank(M_LWE-CC_n) = Θ(q^n). *Evidence:* 31≈5², 53≈7², 126≈5³.

**Conjecture 2 (Main):** Q¹(LWE-CC_n) = O(log n) and R¹(LWE-CC_n) = Ω(n^{1/3}).  
*Connection:* Reduction from VSP (Klartag-Regev 2011) via lattice embedding.

**Conjecture 3:** Q(QV-CC_{N,k}) = O(k^{1/2} log p) for the Quorum Vault threshold reconstruction problem.

## Key Finding: Hard Distribution Structure

The discrepancy LP (coordinate descent, April 2026) finds that the hard distribution μ* for LWE-CC concentrates on:
- Small secrets s (including s=0) with non-zero errors e
- This is the "error-dominated" regime, structurally similar to the hardest VSP instances

This informs the proof strategy: a Klartag-Regev style argument using concentration of measure on the discrete torus ℤ_q^n in the error-dominated regime.

## Open Problems (in order of difficulty)

1. Prove Conjecture 1 (rank = Θ(q^n)) analytically
2. Compute optimal fooling set via ILP for n=2, q=5
3. Find the exact hard distribution μ* via full LP (not coordinate descent)
4. Prove classical lower bound Ω(n^{1/3}) using hard distribution  
5. Prove or disprove quantum upper bound O(log n)

## Connection to QuantumMail-v2

This research was motivated by [QuantumMail-v2](https://github.com/munnamihir/QuantumMail-v2) — a client-side encrypted email platform using ML-KEM-768 (FIPS 203). The key encapsulation in QuantumMail-v2 is precisely LWE-CC as a one-way communication protocol. The quorum vault recovery (2-of-N threshold, PBKDF2-derived token) defines the QV-CC multiparty problem.

## Files

```
code/
  lwe_cc_exact.py     ← first computation of rank, fooling set, discrepancy
  lwe_cc_bounds.py    ← rank proof + discrepancy LP
results/
  lwe_cc_results.json ← all computed bounds
  lp_results.json     ← discrepancy LP output
notes/
  rank-conjecture.md  ← proof sketch for Conjecture 1
  hard-distribution.md ← structure of μ*
drafts/
  lwe-cc-paper-v1.tex ← arXiv draft
tools/
  lwe-cc-exact-bounds.html ← interactive research tool
  lwe-cc-full-tool.html    ← full LWE-CC + QV-CC simulator
```

## Citation

If you build on this work, please cite:
```
@misc{bommisetty2026lwecc,
  author = {Bommisetty, Mihir},
  title  = {Communication Complexity of Lattice-Based Cryptographic Protocols: LWE-CC},
  year   = {2026},
  url    = {https://github.com/munnamihir/lwe-cc-research}
}
```
