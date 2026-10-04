---
title: "The de Bruijn–Newman constant across a family of functions with the functional equation of L(s, χ₋₇)"
authors: ["Nozero Day"]
note: "Produced by an autonomous multi-agent research system; see the repository README. Unrefereed."
status: draft
date: "2026-10-04"
---

# Λ across the Theorem A′ family for q = 7

## Abstract

The real Dirichlet series with exactly the functional equation of L(s, χ₋₇) form a circle f\_θ, with the L-function at θ = 0. We study the de Bruijn–Newman constant Λ(f\_θ) along this circle. The claims are listed, with rigor grades, in claims.md.

- **Near the L-function.**
  - (i) Arb Rouché certificates give Λ > 0.30, 0.12, 0.06 and 0.02 at θ = 1/10, 1/20, 3/100 and 1/50 (claims #1, #4). The responsible zero is a close pair of zeros of L(s, χ₋₇) at heights 50.983 and 51.297, pushed off the line.
  - (ii) A two-zero model, tan θ\_c = g²A/(8B), predicts θ\_c ≈ 0.0132 and the initial depths with no fitted parameter (numeric).
  - (iii) Λ is upper semicontinuous at θ = 0, conditional on one sketched lemma (claim #2).
- **Near θ\***, where a\_θ(1) → 0:
  - (iv) deep zeros sit near σ₀ = log₂|a(2)/a(1)|, and certificates give Λ > 24.5, 30 and 44 at a(1) ≈ 9.45·10⁻⁵, 1.85·10⁻⁶ and 1.0·10⁻¹⁰ (claims #3, #5, #6);
  - (v) a Laplace-lemma argument proves Λ(f\_θ) → ∞ as θ → θ\* (Theorem U, proof-unbounded.md; a written proof, not formalised).
- **Reading.** The functional equation of L(s, χ₋₇) places no upper bound on Λ. What controls Λ is the leading coefficients, and near θ = 0 the L-function's close zero pairs. Nothing here concerns ζ, and none of these functions has an Euler product except θ = 0.

## 1. The family

- **Which sequences.** Take the real 7-periodic odd sequences a with a(0) = 0 that are DFT-self-dual: Σ\_m a(m) e^{2πimn/7} = i√7·a(n).
- **The space is two-dimensional** (papers/supporting/riemann-siegel-identity/proof-rs-general.md, Remark 5). It is spanned by:
  - the real character χ₋₇ = (·/7);
  - b(n) = √2 Re(e^{iψ/2} χ₆(n)), where χ₆ is the odd order-6 character with χ₆(3) = e^{iπ/3} and e^{iψ} = τ(χ̄₆)/(i√7).
  - Exact self-duality of b: its DFT has χ₆-coefficient e^{−iψ/2}τ(χ̄₆) = i√7 e^{iψ/2}, and χ̄₆-coefficient e^{iψ/2}τ(χ₆) = i√7 e^{−iψ/2}, since τ(χ₆) = −conj τ(χ̄₆) for an odd character.
- **The circle.** The vectors a\_θ = cos θ·χ₋₇ + sin θ·b have the same norm for every θ. Each f\_θ(s) = Σ a\_θ(n) n^{−s} satisfies the functional equation of L(s, χ₋₇).
- **Heat flow.** So H\_t(z) = (4πt)^{−1/2} ∫ e^{−w²/4t} ξ((1+w+iz)/2) dw, with ξ(s) = (7/π)^{(s+1)/2} Γ((s+1)/2) f\_θ(s), is real entire and even for every θ.
- **Why a non-real zero bounds Λ.** H\_t is real entire of order 1. For t ≥ Λ it lies in the Laguerre–Pólya class, and the backward heat operator preserves that class. So a non-real zero at t₁ gives Λ > t₁.
- **The endpoint θ = 0** is the Dirichlet L-function. This repository certifies no upper bound for Λ(L(s, χ₋₇)).

## 2. Off-line zeros (numeric-experiment)

`code/02-functional-equation-family/x7_scan.py` scans |f\_θ| on σ ∈ [0.52, 1.1] and t ≤ 300, and refines the minima. In those scans every tested θ ≠ 0 (π/2, π/4, 0.3, 0.1) had off-line zeros; for θ = 0.1, the one closest to σ = 1 is at (0.9283, 51.14). The scan outputs are not stored in this repository. The only off-line zero used below, (0.92832, 51.1378) at θ = 0.1, is recomputed in `runs/x7_track_theta0.1.txt`.

A heuristic predicts this. Joint universality of L(s, χ₋₇), L(s, χ₆) and L(s, χ̄₆) suggests that for every θ ≠ 0 the function f\_θ has zeros with σ arbitrarily close to 1.

## 3. Certified: Λ(f\_{1/10}) > 0.30

- **Tracking.**
  - The θ = 0.1 zero at (0.92832, 51.1378), i.e. z₀ = 102.2757 − 0.8566i, was tracked under the heat flow (`x7_track.py`; `runs/x7_track_theta0.1.txt`).
  - It reaches the real axis between t = 0.30 and t = 0.32. At t = 0.30 it is at 102.08498 − 0.19249i.
- **The certificate.**
  - An Arb Rouché certificate at t = 0.30 (`x7_cert.py`; `certificates/02-functional-equation-family/rouche-theta1_10-t030.json`) shows exactly one zero of H\_{0.30} within 2·10⁻⁵ of that point: |H − p| ≤ 3.17·10⁻⁵ < 6.28·10⁻⁵ ≤ |p|.
  - The coefficient balls contain the exact values; the self-duality residual of the computed balls is below 10⁻⁴⁵.
- **Conclusion.** Hence **Λ(f\_{1/10}) > 0.30**.

## 4. Near the L-function

### 4.1 The rate as θ → 0 (claims #4)

- **The certificates.** Arb Rouché certificates for the zero near height 51.14 give Λ(f\_{1/20}) > 0.12, Λ(f\_{3/100}) > 0.06 and Λ(f\_{1/50}) > 0.02, alongside Λ(f\_{1/10}) > 0.30. Each sits within one tracking step (0.02) of the numerical collision time.
- **The tracks** (`runs/x7_track_theta*.txt`).

  | θ    | initial zero (σ, t) | initial depth y₀ = 2σ − 1 | reaches the axis for t in |
  | ---- | ------------------- | ------------------------- | ------------------------- |
  | 0.1  | (0.92832, 51.1378)  | 0.857                     | (0.30, 0.32)              |
  | 0.05 | (0.77057, 51.1391)  | 0.541                     | (0.12, 0.14)              |
  | 0.03 | (0.68096, 51.1395)  | 0.362                     | (0.06, 0.08)              |
  | 0.02 | (0.61493, 51.1396)  | 0.230                     | (0.02, 0.04)              |

- **Reading.** This fits an on-line close pair of zeros of L(s, χ₋₇) near t ≈ 51.14 that collides at some θ\_c > 0 and leaves the line. So the lower bound from this zero decays as θ → 0 and vanishes at θ\_c, consistent with the semicontinuity of §4.2.
- **The pair** (numeric; `code/02-functional-equation-family/chi7_pair.py`; `runs/chi7_pair_47_55.txt`).
  - L(s, χ₋₇) has on-line zeros at 50.98310 and 51.29663, a gap of 0.3135. The neighbouring gaps are 1.747 and 1.950.
  - The mean spacing there, 2π/log(7t/2π), is about 1.55.
  - The pair's midpoint, 51.140, is the height of the off-line zeros above.
- **The pair-collision law** (numeric; `code/02-functional-equation-family/pair_law.py`; `runs/pair_law.txt`).
  - Near a close pair with z-gap g, model H^χ(z) ≈ (A/2)[(z − z\_m)² − g²/4] and H^b(z) ≈ B. Then H^χ + tan θ·H^b has a double zero at **tan θ\_c = g²A/(8B)**.
  - Beyond it, y₀² = 2(B/A)(tan θ − tan θ\_c), and the quadratic heat flow collides at t = y₀²/2.
  - For this pair, g = 0.62706 and B/A = 3.7259, and the law predicts **θ\_c = 0.013191**.
  - It then predicts y₀ = 0.225, 0.354, 0.524, 0.806 at θ = 0.02, 0.03, 0.05, 0.1, against the tracked 0.230, 0.362, 0.541, 0.857.
  - It predicts collision times 0.025, 0.063, 0.137, 0.325, against the tracked intervals above.
  - No parameter is fitted: the L-function's pair data alone predicts where the family's Λ lower bound switches on and how fast it grows.
  - Consequence (heuristic): pairs with ever smaller g²A/B would give θ\_c → 0, hence Λ(f\_θ) > 0 for every θ ≠ 0.

### 4.2 Upper semicontinuity at the L-function (claim #2; conditional)

limsup\_{θ→0} Λ(f\_θ) ≤ Λ(L(s, χ₋₇)).

The proof (proof-semicontinuity.md) has four steps:
- (0) The heat flow is linear, so H\_t^θ = cos θ·H\_t^χ + sin θ·H\_t^b.
- (1) Uniform far-field dominance: no zeros with |x| ≥ X(t₁, y₀) and y₀ ≤ |y| ≤ Y, with X uniform for small θ. This is Lemma 1′, which is only sketched.
- (2) For t₁ > Λ₀ the zeros of H^0\_{t₁} are real and simple, and Hurwitz keeps them real for small θ.
- (3) de Bruijn's strip lemma finishes.

**Consequences.**
- Under GRH for χ₋₇, Λ(f\_θ) → 0, since Λ₀ ≥ 0 by Dobner's theorem for the extended Selberg class.
- A uniform positive lower bound for Λ(f\_θ), θ ≠ 0, would disprove GRH(χ₋₇).
- So the Euler product is not what keeps Λ small at θ = 0; continuity in the family is.

## 5. Near θ\* (claims #3, #5, #6, #7)

### 5.1 The deep zeros

- **Where.** At θ\* = arctan(−1/b(1)) = 0.704061195… (Arb; `code/02-functional-equation-family/theta_star.py`, `runs/theta_star.txt`) the leading coefficient a\_θ(1) vanishes.
- **The two-term balance.** Near θ\*, f\_θ ≈ a(1) + a(2)2^{−s}. So f\_θ has zeros near σ₀ = log₂|a(2)/a(1)|, repeating in t with period 2π/log 2 (`deep_zeros.py`; `runs/deep_zeros_*.txt`).

| θ                              | a\_θ(1)     | σ₀     | lowest deep zero (σ, t) | initial depth −Im z | certified |
| ------------------------------ | ---------- | ------ | ----------------------- | ------------------- | --------- |
| 704/1000                       | 9.454·10⁻⁵ | 12.863 | (12.8674, 4.5492)       | 24.735              | Λ > 24.5  |
| 70406/10⁵                      | 1.847·10⁻⁶ | 18.541 | (18.5413, 4.5341)       | 36.083              | Λ > 30    |
| 70406119530951995496/10²⁰      | 1.0·10⁻¹⁰  | 32.713 | (32.7135, 4.5324)       | 64.427              | Λ > 44    |

**The certificates.**
- Each is an Arb Rouché certificate (`x7_cert.py` in high-precision mode) at the stated t, for a non-real zero of H\_t. (In the heat-flow tracks these zeros come from the deep zeros. Those long tracks, to t ≈ 24.7, 31.4 and 44.5, are not stored here; only the first t ≤ 1 is, in `runs/x7_track_deep_*.txt`.)
- They use a Stirling-type tail bound and the reflection ξ(s) = ξ(1 − s); see claims.md #3.
- The Stirling tail bound exceeds the numerically integrated tail by a factor of about 1.4·10⁵ at t = 24.5 (`tail_check.py`; `runs/tail_check_t2450.txt`).

### 5.2 Theorem U and the rate

- **Theorem U** (proof-unbounded.md). Λ(f\_θ) → ∞ as θ → θ\*.
- **The initial drift.** The proof predicts the initial drift |d Im z/dt| = Re φ′(s₁) − ½ log 2 of the deep zero, where φ is the log of the Γ-factor.
  - Predicted (`drift_check.py`; `runs/drift_check.txt`): 1.0144, 1.1822 and 1.4562.
  - Tracked (`runs/x7_track_deep_*.txt`): 1.0134, 1.1817 and 1.4560.
  - These agree to 0.1%.
- **Heuristic rate.** If that drift held all the way to the axis, the collision time would be about 4σ₀/log σ₀. This suggests Λ(f\_θ) ≍ log(1/|a\_θ(1)|)/log log(1/|a\_θ(1)|) as θ → θ\*. The rate is a conjecture; only unboundedness is proved.
- **At θ\* itself** Λ is finite, since all zeros of f\_{θ\*} lie in a vertical strip (its leading coefficient is a(2) ≠ 0), and de Bruijn's strip theorem applies.

## 6. Ceiling

These are non-multiplicative functions; nothing here concerns ζ. The result measures how special the L-function is within its functional-equation family.

## References

- A. Dobner, A proof of Newman's conjecture for the extended Selberg class, Acta Arith. 201 (2021), 29–62, doi:10.4064/aa200603-23-7, arXiv:2005.05142.
- L. Ferry, D. Ghisa, F. A. Muscutar, On Davenport and Heilbronn-type of functions, arXiv:1602.06328 (not read; title and abstract only).
- H. Heilbronn, On Dirichlet series which satisfy a certain functional equation, Q. J. Math. 9 (1938), 194–195, doi:10.1093/qmath/os-9.1.194 (not read; the family of Dirichlet series with a given functional equation goes back to this paper).
- D.H.J. Polymath, Effective approximation of heat flow evolution of the Riemann ξ function, and a new upper bound for the de Bruijn–Newman constant, Res. Math. Sci. 6 (2019), Paper No. 31, doi:10.1007/s40687-019-0193-1, arXiv:1904.12438.
