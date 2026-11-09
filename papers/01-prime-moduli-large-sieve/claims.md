# Claims and evidence — the log log Q in Wolke's prime-moduli large sieve is necessary

| # | Claim (precise) | Rigor | Digits | Novelty | Certificate | Lean statement | Witness | Notes |
| --- | --------------- | ----- | ------ | ------- | ----------- | -------------- | ------- | ----- |
| 1 | Let 0 < η ≤ 1/2, k ≥ 1 and N = ⌊k(1−η)Q/η⌋. Then there is x₀ such that a\_n = e(−n x₀), 1 ≤ n ≤ N, satisfies Σ\_{p≤Q} Σ\_{a=1}^{p−1} \|S(a/p)\|² ≥ (4/π²)·(π(Q) − π((1−η)Q) − ω(k))/φ(k)·N·Σ\|a\_n\|². If moreover φ(k) ≤ ηQ/(2 log Q), the right side is ≥ (4(1−η)/π² + o(1))·(k/φ(k))·(Q²/log Q)·Σ\|a\_n\|² | conjectural |  | plausibly new; the method is that of Baier–Lynch–Zhao 2019 and Corrigan 2026 |  |  |  | Written proof in main.md; its finite core is proved in Lean (#3). The PNT, Mertens and Jordan steps are written, not formalised. Numerical check at Q = 4000: code/01-prime-moduli-large-sieve/claim1\_check.py. |
| 2 | Corollary. For N = N(Q) with N/Q → ∞ and N ≤ Q²/(4 log Q), the optimal prime-moduli large-sieve constant satisfies Δ\_P(Q,N) ≥ (2/π² + o(1))·e^γ log log(N/Q)·Q²/log Q. For N = Q^{2−2θ}, fixed 0 < θ < 1/2, the ratio Δ\_P(Q,N)/(Q²/log(Q/√N)) is at least (2θe^γ/π² + o(1)) log log Q → ∞. So the factor log log Q in Wolke's bound Q² log log Q/log(Q/√N) (Iwaniec 2022, Thm 1) cannot be removed for general coefficients, and Wolke's bound has the right order there | conjectural |  |  |  |  |  | This answers the question in Iwaniec, Rev. Mat. Iberoam. 38 (2022), p. 2338. Novelty is plausible but not established: Wolke's Math. Z. 1970, Acta Hungar. 1971 and Bull. LMS 1974 papers, Halberstam–Richert and Opera de Cribro are unread (main.md, Scope). |
| 3 | Finite core, for naturals k, p, b. If 0 < k, 0 < p and k ∣ 1 + b·p, then a := (1 + b·p)/k satisfies a/p = b/k + 1/(k·p) in ℚ. If moreover b < k and 2 ≤ p, then 1 ≤ a < p. For 1 < k and a finite set s of naturals coprime to k, some reduced residue class c < k contains at least #s/φ(k) elements of s (as a rational bound; ⌊#s/φ(k)⌋ in the natural-number form) | proved |  |  |  | lean/DebruijnNewmanLargeSieve/PrimeSpoke.lean: DebruijnNewmanLargeSieve.PrimeSpoke.spoke\_point, spoke\_bounds, exists\_rich\_class\_rat (and exists\_rich\_class) | spoke\_example, rich\_class\_example, spoke\_hyps\_satisfiable | All seven theorems of the file depend only on propext, Classical.choice and Quot.sound (lean/AxiomCheck.lean). |

## Reasoning

### #1

- Step 1 (pigeonhole over the φ(k) reduced classes mod k) selects a class c holding at least its share of the primes p ∈ ((1−η)Q, Q] with p ∤ k.
- Take b with bc ≡ −1 (mod k). Then k ∣ 1 + bp, and a\_p := (1 + bp)/k gives the Farey point a\_p/p = b/k + 1/(kp), with 1 ≤ a\_p < p. These are exactly the three facts proved in Lean in #3.
- As p ranges over ((1−η)Q, Q], these points fill an interval of length η/(k(1−η)Q). With N = ⌊k(1−η)Q/η⌋, every point lies within 1/(2N) of the midpoint x₀. There |S(a\_p/p)| = |sin(πNθ)/sin(πθ)| ≥ 2N/π, by Jordan's inequality.
- Dropping all other Farey points only lowers the left side, since the terms are nonnegative. So LS ≥ #class · 4N²/π², and Σ|a\_n|² = N gives the stated bound. The PNT for a fixed η supplies the asymptotic form.
- The proof was reviewed adversarially by a second, independent agent of the research system.
- A finite check at Q = 4000 confirms the non-asymptotic inequality for k = 1, 2, 6, 30, 210 (`code/01-prime-moduli-large-sieve/runs/claim1_check.txt`).
- What would change the conclusion: an error in the Farey-point bookkeeping (formalised, so ruled out) or in the Jordan step (a standard inequality).

### #2

- With η = 1/2 and k the largest primorial ≤ N/Q, Mertens gives k/φ(k) = Π\_{p≤z}(1 − 1/p)^{−1} = (e^γ + o(1)) log z, and log k = ϑ(z) ~ z.
- Consecutive primorials differ by a factor p\_{next} = exp(o(z)). So log k = (1 − o(1)) log(N/Q), and k/φ(k) = (e^γ + o(1)) log log(N/Q). This needs N/Q → ∞.
- For N = Q^{2−2θ}, log(Q/√N) = θ log Q. Dividing the bound of #1 by Q²/log(Q/√N) gives (2θe^γ/π²) log log Q → ∞.
- Iwaniec's remark after his Theorem 1 [Iw22, p. 2342] names exactly this scenario ("many rationals a/q are very close to the same point b/c") as where his method loses c/φ(c). The Corollary shows the loss is intrinsic for general coefficients.
- His Theorem 2, which removes log log Q when a\_n is supported on (n, P) = 1, is consistent with this: the extremal sequence a\_n = e(−nx₀) is supported on all n.
- Novelty is the open risk: the literature listed in the Notes column has not been read.

### #3

These are elementary facts about naturals and ℚ, so they are formalised exactly:
- spoke\_point: from k·a = 1 + b·p, dividing by k·p gives a/p = b/k + 1/(kp) (field\_simp plus linear\_combination). The hypothesis 0 < p is needed (in Lean, division by 0 is 0); 0 < k already follows from k ∣ 1 + b·p.
- spoke\_bounds: 1 ≤ a because k·a ≥ 1, and a < p because k·a = 1 + bp ≤ 1 + (k−1)p < kp.
- exists\_rich\_class\_rat: Mathlib's exists\_le\_card\_fiber\_of\_nsmul\_le\_card\_of\_maps\_to over ℚ, applied to p ↦ p % k into the φ(k) reduced residues (Nat.totient\_eq\_card\_coprime), gives the exact bound #𝒫\_c ≥ #𝒫/φ(k) of Step 1. Coprimality passes to p % k by Nat.gcd\_rec. exists\_rich\_class is the natural-number form, with ⌊#s/φ(k)⌋.

Witnesses: spoke\_example applies spoke\_point and spoke\_bounds in the paper's setting (k = 6, p = 11, b = 1: 2/11 = 1/6 + 1/66, 1 ≤ 2 < 11); rich\_class\_example applies exists\_rich\_class to k = 6 and s = {5, 7, 11, 13}, where some reduced class holds at least two elements; spoke\_hyps\_satisfiable checks the hypotheses at k = 6, p = 5, b = 1.
