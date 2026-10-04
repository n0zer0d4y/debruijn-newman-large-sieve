---
title: "The log log Q in Wolke's large sieve with prime moduli is necessary"
authors: ["Nozero Day"]
note: "Produced by an autonomous multi-agent research system; see the repository README. Unrefereed."
status: draft
date: "2026-10-04"
---

# The log log Q in Wolke's large sieve with prime moduli is necessary

## Summary

**Background.** Wolke [Wo71, Wo71b] proved the large sieve inequality with prime moduli

  Σ\_{p≤Q} Σ\_{a=1}^{p−1} |S(a/p)|² ≪ (Q² log log Q / log(Q/√N)) · Σ|a\_n|²,  where S(x) = Σ\_{n≤N} a\_n e(nx),  (W)

for N ≤ Q²(log Q)^{−8}, in the form stated by Iwaniec [Iw22, Thm 1, p. 2338]. Iwaniec removed the factor log log Q for coefficients supported on integers without small prime factors [Iw22, Thm 2], and wrote [Iw22, p. 2338]:

> "We do not know if the bound (1.3) holds true without the factor log log Q for general coefficients."

Baier's survey [Ba25, §1.2] calls the factor "slightly unsatisfactory" and records no lower bound.

**The answer.** The factor cannot be removed. For every fixed θ ∈ (0, ½) and N = Q^{2−2θ}, there are coefficients with

  Σ\_{p≤Q} Σ\_{a=1}^{p−1} |S(a/p)|² ≥ (2θ e^γ/π² + o(1)) · log log Q · (Q²/log(Q/√N)) · Σ|a\_n|².

So (W) is sharp up to a constant factor in this range, and the true order of the prime-moduli large-sieve constant there is Q² log log Q/log Q.

**Mechanism.**
- This is the prime-denominator analogue of the Baier–Lynch–Zhao lower bound for square moduli [BLZ19], which Corrigan [Co26] carried over to quadratic-polynomial moduli. The method is theirs; the new point is its application to prime moduli.
- Corrigan asks a related question: whether averaging over primes removes the corresponding factor for quadratic moduli f(p) [Co26, p. 29].
- Farey fractions a/p with prime p accumulate near b/k for a primorial k. The "spoke" m = ak − bp = 1 picks out the primes in one residue class mod k.
- There are π(Q)/φ(k) of them, rather than the Q/k that integers would give, and k/φ(k) is as large as e^γ log log k (Mertens).

## Statement

**Theorem.** Let 0 < η ≤ ½ and let k ≥ 1 be an integer. Put

N := ⌊k(1 − η)Q/η⌋.

Then there is x₀ ∈ ℝ such that a\_n := e(−nx₀), 1 ≤ n ≤ N, satisfies

  Σ\_{p≤Q} Σ\_{a=1}^{p−1} |S(a/p)|² ≥ (4/π²) · (Π\_η(Q) − ω(k))/φ(k) · N · Σ\_n |a\_n|²,

where Π\_η(Q) := π(Q) − π((1 − η)Q).

This inequality holds unconditionally. It is non-trivial once φ(k) ≤ ηQ/(2 log Q), since then #𝒫\_c ≥ 1 for large Q. By the prime number theorem, Π\_η(Q) = (η + o(1))Q/log Q. Hence, under the hypothesis φ(k) ≤ ηQ/(2 log Q), which makes the o(1) uniform in k,

  Σ\_{p≤Q} Σ\_a |S(a/p)|² ≥ (4(1 − η)/π² + o(1)) · (k/φ(k)) · (Q²/log Q) · Σ|a\_n|².  (1)

**Corollary.** Let N = N(Q) with N/Q → ∞ and N ≤ Q²/(4 log Q), and let k be the largest primorial with k ≤ N/Q. Then the optimal constant Δ\_P(Q, N) in (W)'s shape satisfies

  Δ\_P(Q, N) ≥ (2/π² + o(1)) · e^γ log log(N/Q) · Q²/log Q.

(N/Q → ∞ is needed for log log(N/Q) to be positive and asymptotic. It holds throughout Q^{1+δ} ≤ N ≤ Q^{2−δ}.)

In particular, for N = Q^{2−2θ} with fixed 0 < θ < ½,

  Δ\_P(Q, N) / (Q²/log(Q/√N)) ≥ (2θ e^γ/π² + o(1)) · log log Q → ∞.

## Proof

**Step 1: a spoke of prime-denominator Farey points.** Assume k > 1; the case k = 1 is below.
- Let 𝒫 be the set of primes in ((1 − η)Q, Q] not dividing k, so #𝒫 ≥ Π\_η(Q) − ω(k).
- These primes fall into the φ(k) reduced residue classes mod k. By pigeonhole some class c has #𝒫\_c ≥ #𝒫/φ(k), where 𝒫\_c := {p ∈ 𝒫 : p ≡ c (mod k)}.
- Let b ∈ [1, k − 1] satisfy bc ≡ −1 (mod k). For p ∈ 𝒫\_c we then have k | 1 + bp. Put a\_p := (1 + bp)/k.
- Then 1 ≤ a\_p ≤ (1 + (k − 1)p)/k < p, so a\_p/p is one of the points in the sum. Moreover

  a\_p/p = b/k + 1/(kp).

- For k = 1 take b = 0, a\_p = 1, and 𝒫\_c the primes in ((1 − η)Q, Q].
- Distinct p give distinct points.

**Step 2: the points lie in a short interval.** For p ∈ ((1 − η)Q, Q], the points b/k + 1/(kp) lie in an interval of length

1/(k(1 − η)Q) − 1/(kQ) = η/(k(1 − η)Q).

Let x₀ be its midpoint. Every such point then satisfies |a\_p/p − x₀| ≤ η/(2k(1 − η)Q) ≤ 1/(2N), by the choice of N.

**Step 3: each point contributes about N².** For θ := a\_p/p − x₀ we have |Nθ| ≤ ½ and

|S(a\_p/p)| = |Σ\_{n≤N} e(nθ)| = |sin(πNθ)|/|sin(πθ)|.

Jordan's inequality gives |sin(πNθ)| ≥ 2N|θ| for |Nθ| ≤ ½, and |sin(πθ)| ≤ π|θ|. (For θ = 0 the sum is N.) Hence |S(a\_p/p)| ≥ 2N/π, and

  Σ\_{p≤Q} Σ\_a |S(a/p)|² ≥ #𝒫\_c · 4N²/π² ≥ (4/π²) · (Π\_η(Q) − ω(k))/φ(k) · N · Σ|a\_n|²,

because Σ|a\_n|² = N.

**Step 4: the asymptotic form (1).** Insert Π\_η(Q) = (η + o(1))Q/log Q, ω(k) ≤ log₂ k = O(log Q) = o(Q/log Q), and N = k(1 − η)Q/η + O(1):

(4/π²) · (ηQ/(φ(k) log Q)) · (k(1 − η)Q/η) = (4(1 − η)/π²) · (k/φ(k)) · Q²/log Q. ∎

**Proof of the Corollary.**
- Take η = ½, so the Theorem applies with N₁ = ⌊kQ⌋ ≤ N. A sequence of length N₁ extended by zeros is a sequence of length N, so Δ\_P(Q, N) ≥ Δ\_P(Q, N₁).
- The hypothesis φ(k) ≤ Q/(4 log Q) holds because φ(k) ≤ k ≤ N/Q ≤ Q/(4 log Q). This uses the Corollary's range restriction N ≤ Q²/(4 log Q), which contains Wolke's range N ≤ Q²(log Q)^{−8}.
- For the primorial k = P(z) = Π\_{p≤z} p, Mertens gives k/φ(k) = Π\_{p≤z}(1 − 1/p)^{−1} = (e^γ + o(1)) log z. Also log k = ϑ(z) = (1 + o(1))z.
- The largest primorial ≤ N/Q has log k ≥ log(N/Q) − log z\_{next} = (1 − o(1)) log(N/Q). Hence k/φ(k) = (e^γ + o(1)) log log(N/Q).
- For N = Q^{2−2θ}, log log(N/Q) = log log Q + O(1), and log(Q/√N) = θ log Q. ∎

## Remarks

0. **Iwaniec anticipated the obstruction.** In the remarks after his proof of Theorem 1 [Iw22, p. 2342] he writes: "We lost the factor log log Q in (2.7) after applying the maximum bound for c/φ(c). … Unfortunately, if H is small, then this coprimality condition does not help to recover φ(c)/c. This case occurs when many rationals a/q are very close to the same point b/c." The Theorem realises exactly that configuration, with c = k a primorial. So the loss is genuine, not an artefact of the method.
1. **Iwaniec's Theorem 2 is consistent with this.** The extremal sequence a\_n = e(−nx₀) is supported on all n ≤ N, including those with small prime factors. Iwaniec's restriction (n, P) = 1 is exactly what destroys the resonance with the modulus k: the resonance comes from the spoke structure mod k, which a coefficient sequence sifted by the primes dividing k cannot see coherently.
2. **Analogy.** The argument parallels [BLZ19]. Their lower bound for square moduli has size Q^{(log 2)/((1+ε) log log Q)}(Q³ + N), so the factor Q^ε in Zhao's conjectured bound cannot be replaced by a power of log Q. There too, a product of the first m odd primes plays the role of k, and the fractions cluster near a single point.
3. **Numerical illustration** (numeric-experiment; `code/01-prime-moduli-large-sieve/claim1_check.py`, output in `runs/claim1_check.txt`).
   - Setting: Q = 4000, η = ¼, with the best class chosen as in Step 1. Here LS := Σ\_p Σ\_a |S(a/p)|², ‖a‖² = N, and #points = Σ\_{p≤Q}(p − 1) = 1012957.
   - **The Theorem's choice:** N = ⌊k(1 − η)Q/η⌋ = 3kQ, x₀ the midpoint. The non-asymptotic inequality holds for k = 1, 2, 6, 30, 210, and every spoke point satisfies |S| ≥ 2N/π.
   - **A smaller N:** N = kQ and x₀ = b/k + 1/(kQ(1 − η/2)). The quantity (LS/‖a‖²)/(#points + N), the operator norm of this test vector relative to the classical scale "#points + N", is 0.762, 1.450, 2.130, 2.286 and 2.704 for k = 1, 2, 6, 30, 210, against k/φ(k) = 1, 2, 3, 3.75, 4.375.
   - The ratio to k/φ(k) is 0.76, 0.73, 0.71 for k = 1, 2, 6, while the chosen class still has many primes. At k = 30 and 210 only 17 and 5 primes remain in the class at this Q, and it falls to 0.61 and 0.62.

## Scope

- This is a written proof. It is elementary, using only the pigeonhole principle, the prime number theorem, Mertens' theorem and Jordan's inequality. Its arithmetic core (the spoke identity, 1 ≤ a\_p < p, and the pigeonhole step) is formalised in Lean (claims #3); the analytic steps are not.
- **Novelty: plausibly new, not established.**
  - Iwaniec states the question as open, and for Theorem 1 cites only Wolke's papers [Wo71, Wo71b].
  - Baier [Ba25] and Baier–Lynch–Zhao [BLZ19] record no prime-moduli lower bound, and neither does Corrigan [Co26].
  - Montgomery's survey [Mo78, pp. 563–564] mentions Wolke's lower-bound work in connection with the all-moduli inequality (N + Q²)Σ|a\_n|².
  - Baier's paper on sparse sets of moduli [Ba06] gives no lower bound, and Schlage-Puchta [SP09] concerns random coefficients.
- **Not yet read, and required before novelty can be asserted:**
  - Wolke [Wo70] and [Wo71b];
  - Wolke [Wo74] (closed access);
  - Elliott's endpoint example (quoted in [Wo71]);
  - Halberstam–Richert, *Sieve Methods*;
  - Friedlander–Iwaniec, *Opera de Cribro*.

## References

- [Ba06] S. Baier, On the large sieve with sparse sets of moduli, J. Ramanujan Math. Soc. 21 (2006), 279–295, arXiv:math/0508417.
- [Ba25] S. Baier, The large sieve for square moduli, revisited, Hardy-Ramanujan J. 48 (2025), doi:10.46298/hrj.2026.17906, arXiv:2503.18009.
- [BLZ19] S. Baier, S. B. Lynch, L. Zhao, A lower bound for the large sieve with square moduli, Bull. Aust. Math. Soc. 100 (2019), 225–229, doi:10.1017/S0004972719000224, arXiv:1812.05844.
- [Co26] C. C. Corrigan, A large sieve inequality for characters to quadratic moduli, Acta Arith. 222 (2026), 295–325, doi:10.4064/aa240923-21-1, arXiv:2410.05704.
- [Iw22] H. Iwaniec, The large sieve with prime moduli, Rev. Mat. Iberoam. 38 (2022), 2337–2354, doi:10.4171/RMI/1381.
- [Mo78] H. L. Montgomery, The analytic principle of the large sieve, Bull. Amer. Math. Soc. 84 (1978), 547–567, doi:10.1090/S0002-9904-1978-14497-8.
- [SP09] J.-C. Schlage-Puchta, Lower bounds for expressions of large sieve type, Arch. Math. (Brno) 45 (2009), 79–82, arXiv:1105.1307.
- [Wo70] D. Wolke, Farey-Brüche mit primem Nenner und das große Sieb, Math. Z. 114 (1970), 145–158, doi:10.1007/BF01110322.
- [Wo71] D. Wolke, Farey fractions with prime denominator and the large sieve, Mém. Soc. Math. France 25 (1971), 183–188, doi:10.24033/msmf.52.
- [Wo71b] D. Wolke, On the large sieve with primes, Acta Math. Acad. Sci. Hungar. 22 (1971), 239–247, doi:10.1007/BF01896016.
- [Wo74] D. Wolke, A lower bound for the large sieve inequality, Bull. London Math. Soc. 6 (1974), 315–318, doi:10.1112/blms/6.3.315.
