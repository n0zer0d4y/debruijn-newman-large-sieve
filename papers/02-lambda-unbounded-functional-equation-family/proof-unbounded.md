# Unboundedness of Λ on the χ₋₇ functional-equation circle

A written proof, reviewed adversarially by a second, independent agent of the research system. It is not formalised, so its rigor grade is conjectural.

**Theorem U.** Let a\_θ = cos θ·(·/7) + sin θ·b be the q = 7 odd DFT-self-dual family of main.md §1. Let θ\* be the angle with a\_{θ\*}(1) = 0, which exists because b(1) ≠ 0; note a\_{θ\*}(2) ≠ 0. Then Λ(f\_θ) → +∞ as θ → θ\*. In particular sup\_θ Λ(f\_θ) = +∞, although every f\_θ has the functional equation of L(s, χ₋₇).

More generally, the argument applies to any real family of self-dual Dirichlet series with a fixed functional equation of degree 1 such that a(1) → 0, a(2) stays bounded away from 0, and |a(n)| ≤ Cn^A uniformly in the family for some constants C, A.

## Notation

- **The Γ-factor.** G(s) = κ (7/π)^{(s+1)/2} Γ((s+1)/2) with a constant κ > 0, φ(s) = log G(s), and ξ\_θ(s) = G(s) f\_θ(s), so that ξ\_θ(s) = ξ\_θ(1 − s) (Theorem A′).
- **The heat flow.** For t > 0, H\_t(z) = (4πt)^{−1/2} ∫ e^{−w²/4t} ξ\_θ((1 + w + iz)/2) dw = E[ξ\_θ(s + W/2)], where s = (1 + iz)/2 and W ~ N(0, 2t).
- **Non-real zeros.** A zero z of H\_t with Re s ≠ ½ is non-real, since Im z = −(2 Re s − 1).
- **Why a non-real zero bounds Λ.** H\_t is real entire and even of order 1 (main.md #1). Evolving to larger t (H\_t solves the backward heat equation in z, and the operator e^{−λ∂\_z²} with λ > 0 maps Laguerre–Pólya functions to Laguerre–Pólya functions) preserves real-rootedness. So {t : H\_t has only real zeros} = [Λ, ∞) (de Bruijn), and a non-real zero of H\_t gives t < Λ.
- **Parameters.** Put ε = |a\_θ(1)| → 0 and σ₀ = log(|a\_θ(2)|/ε)/log 2 → ∞, with sgn = sign(a\_θ(2)/a\_θ(1)) ∈ {±1}, which is constant for θ on one side of θ\*.
- **Constants.** All constants below depend on t only (t fixed), and o(1) means as θ → θ\*. Because a\_θ(n) is continuous in θ, |a\_θ(n)| ≤ C for all n and θ, and |a\_θ(2)| ≥ c₂ > 0 near θ\*.

## Step 1: the window

- **The point s₀.** Fix T₀ = π/log 2, so that the unperturbed two-term zeros sit at Im s = (2k+1)π/log 2 when sgn = +1, or at 2kπ/log 2 when sgn = −1; choose the k with Im s₀ in [2, 12].
- **The window.** Put L = t log σ₀ + 2 and S = {s : |Re s − σ₀| ≤ L, |Im s − Im s₀| ≤ 1}.
- **Size of φ and its derivatives on S.** By Stirling for σ = Re s → ∞ uniformly in bounded Im s:
  - φ′(s) = ½ψ((s+1)/2) + ½ log(7/π) = ½ log σ + O(1);
  - φ″(s) = O(1/σ);
  - φ‴(s) = O(1/σ²).
- **Where the Gaussian mass lies.** The Laplace analysis below concentrates the mass of W where |W| ≤ W₁ := t(log σ₀ + C) + (log σ₀)². There Re(s + W/2) ≥ σ₀ − L − W₁/2 ≥ σ₀/2 for θ near θ\*.

## Step 2: the Laplace lemma

**Lemma L.** For n ∈ {1, 2}, uniformly for s ∈ S, as σ₀ → ∞:

E[G(s + W/2) n^{−(s + W/2)}] = G(s) n^{−s} exp((t/4)(φ′(s) − log n)²) · (1 + O((log σ₀)³/σ₀)).

*Proof.*
- **The integrand.** Write g(w) = φ(s + w/2) − (w/2) log n − w²/(4t). Taylor expansion gives φ(s + w/2) = φ(s) + (w/2)φ′(s) + (w²/8)φ″(s) + R₃(w), with |R₃(w)| ≤ |w|³ sup|φ‴|/48 on the segment.
- **The saddle.** It is w\* = t(φ′ − log n)/(1 − tφ″/2), which is O(t log σ₀).
- **The quadratic part.** It gives exactly the Gaussian integral (1 − tφ″/2)^{−1/2} exp((t/4)(φ′ − log n)²/(1 − tφ″/2)). Since tφ″ = O(1/σ₀) and (φ′)² = O((log σ₀)²), this equals exp((t/4)(φ′ − log n)²)(1 + O((log σ₀)²/σ₀)).
- **The cubic remainder.** On |w − w\*| ≤ (log σ₀)², R₃ = O((log σ₀)⁶/σ₀²) = o(1).
- **Absolute values versus the complex integral.** The remainder bounds integrate |integrand| along real w, while the main term is |∫ integrand|. Their ratio is about exp(t(Im α)²), with α = (φ′ − log n)/2. Since Im φ′ = O(|Im s|/σ) = O(1/σ₀) on S, this ratio is 1 + o(1). So the relative errors below are uniform on S.
- **The rest.** The mass outside is O(e^{−(log σ₀)⁴/(8t)}) relative, because g is strictly concave there, with g″ ≤ −1/(2t) + O(1/σ₀).
- **Arguments with Re(s + w/2) < σ₀/2.** Here |w| ≥ σ₀ − 2L. The Gaussian factor is ≤ e^{−w²/(4t)}, while log|ξ\_θ(s + w/2)| ≤ C(σ₀ + |w|) log(σ₀ + |w|). This uses three facts: f\_θ is bounded for Re ≥ ½ at bounded Im; the reflection ξ(s) = ξ(1 − s) covers Re < ½; and |Γ(x + iy)| ≤ Γ(x). Since w²/(4t) dominates C(σ₀ + |w|) log(σ₀ + |w|) for |w| ≥ σ₀ − 2L, the contribution is at most e^{−σ₀²/(5t) + O(σ₀ log σ₀)}, which is negligible relative to |G(s)|e^{(t/4)(φ′)²} ≥ e^{−O(σ₀ log σ₀)} on S. □

**The other terms of f\_θ.** Only an upper bound is needed for n ≥ 3, not an asymptotic. (For log n ≍ σ₀/t the saddle is at |w\*| ≍ σ₀, where the cubic remainder |w|³φ‴ ≍ σ₀ is not small, so Lemma L itself does not extend to those n.) Split at n₁ = exp(cσ₀/t), with c small.
- *Reduction to a real-variable bound.* Since |Γ(x + iy)| ≤ Γ(x) for x > 0, we have |G(σ′ + iτ)| ≤ G\_ℝ(σ′) := κ(7/π)^{(σ′+1)/2}Γ((σ′+1)/2), and |n^{−s′}| = n^{−Re s′}. So the near-region part of the n-th term is at most E[G\_ℝ(σ + W/2) n^{−σ−W/2}; σ + W/2 ≥ σ₀/2].
- *A quadratic majorant.* On x ≥ σ₀/2, φ\_ℝ := log G\_ℝ satisfies φ\_ℝ″(x) ≤ C/σ₀. Hence φ\_ℝ(x) ≤ φ\_ℝ(σ) + φ\_ℝ′(σ)(x − σ) + (C/(2σ₀))(x − σ)². The Gaussian integral of the majorant gives, for 3 ≤ n ≤ n₁, a bound G\_ℝ(σ) n^{−σ} exp((t/4)(φ\_ℝ′(σ) − log n)²(1 + O(t/σ₀)))(1 + O(t/σ₀)).
  - With log n ≤ cσ₀/t, the relative error in the exponent is e^{O(c²σ₀)}.
  - The factor n^{−σ} e^{(t/4)log²n} ≤ e^{−(σ − (t/4)log n) log n} ≤ n^{−σ₀/2} (for c small) absorbs that error.
- *For n > n₁,* the near-region part is at most |a(n)|·n^{−σ₀/2}·E|G(s + W/2)|, because n^{−Re(s + W/2)} ≤ n^{−σ₀/2} there. Since E|G(s + W/2)| ≤ E[G\_ℝ(σ + W/2)] ≍ G\_ℝ(σ)e^{(t/4)φ\_ℝ′(σ)²}, the sum over n > n₁ is ≪ n₁^{1−σ₀/2}·(that), which is negligible.
- *Beyond the near region,* the far-region bound above applies to the whole of ξ\_θ.
- Comparing with |G(s)| (|G(s)|/G\_ℝ(σ) → 1 for bounded Im s) gives the bound on E below.

On the Gaussian mass region, Re(s + W/2) ≥ σ₀/2. Summing over n ≥ 3 with |a(n)| ≤ C, for ξ\_θ = G·Σ a(n) n^{−s} (the term n = 1 has coefficient a(1), of size ε):

H\_t(z) = G(s)e^{(t/4)φ′²} [a(1) + a(2) 2^{−s} e^{(t/4)((φ′ − log 2)² − φ′²)} + E],

with |E| ≤ C' 3^{−σ} e^{(t/4)((|φ′| + log 3)² − |φ′|²)} + (reflected and tail parts) = O(3^{−σ₀ + L} σ₀^{C t}).

## Step 3: normalise and compare with the two-term model

- **The model.** Divide by a(1) G(s) e^{(t/4)φ′²}. With s = σ₀ + u and a(2)2^{−σ₀}/a(1) = sgn,

  Φ\_θ(u) := H\_t/(a(1) G e^{(t/4)φ′²}) = 1 + sgn · 2^{−u} · exp(−(t/4)(2φ′(s) log 2 − log²2)) + E/a(1).

- **The error is small.** Since ε = |a(2)| 2^{−σ₀}, we get E/a(1) = O(2^{σ₀} 3^{−σ₀} 3^{L} σ₀^{Ct}) = O((2/3)^{σ₀} σ₀^{C't}) → 0 uniformly on S.
- **Freezing φ′.** φ′(σ₀ + u) = φ′(s₁) + O(L/σ₀) on S, where s₁ = σ₀ + i Im s₀. So the exponential factor equals K := exp(−(t/4)(2φ′(s₁) log 2 − log²2)) times (1 + O(L log σ₀/σ₀)).
- **Conclusion of the step.** Put M(u) = 1 + sgn·K·2^{−u}. Then Φ\_θ(u) = M(u) + o(1) uniformly on the disk |u − u\*| ≤ ρ of Step 4, where |K·2^{−u}| = |2^{−(u − u\*)}| ≤ 2^ρ is bounded. So Lemma L's relative error times the n = 2 term is o(1) there.
  - This is not claimed on the whole window. At its left edge |K 2^{−u}| grows like σ₀^{(t/2) log 2}, and the comparison can fail for t > 2/log 2. Rouché needs only the circle.

## Step 4: the zero

- **The model zero.** M(u) = 0 ⇔ 2^{−u} = −sgn/K ⇔ u ∈ u\* + (2πi/log 2)ℤ, where u\* = (log K − log(−sgn))/log 2 for any fixed branch of the logarithm. Then:
  - Re u\* = log|K|/log 2 = −(t/4)(2 Re φ′(s₁) log 2 − log²2)/log 2 = −(t/2)(Re φ′(s₁) − (log 2)/2);
  - the zeros are simple and spaced 2π/log 2 apart vertically.
- **It lies in the window.** Re φ′(s₁) = ½ log σ₀ + O(1), so |Re u\*| ≤ (t/4) log σ₀ + O(t) < L − 1. Choose the branch with |Im(σ₀ + u\*) − Im s₀| ≤ 1 by the choice of k, since Im K is bounded.
- **Rouché.** On the circle |u − u\*| = ρ with ρ = min(1/4, π/(2 log 2)), |M(u)| ≥ m(ρ) > 0 uniformly: |M(u)| = |1 − 2^{−(u − u\*)}| has no other zeros within ρ. Since Φ\_θ − M = o(1) there, Rouché gives Φ\_θ, hence H\_t, a zero u\_θ with |u\_θ − u\*| < ρ, for θ close to θ\*.
- **It is non-real.** Its real part in s is σ₀ + Re u\_θ ≥ σ₀ − (t/4) log σ₀ − O(t) → ∞, so it is ≠ ½. Hence z\_θ is a non-real zero of H\_t, and Λ(f\_θ) > t.

Since t was arbitrary, Λ(f\_θ) → ∞. □

## Remarks

- **The rate.** The proof is qualitative, for each fixed t: it gives Λ(f\_θ) > t once θ is close enough to θ\*. A quantitative rate such as Λ ≥ c σ₀/log σ₀ would need t ≍ σ₀/log σ₀. There the cubic Laplace remainder |w|³φ‴ ~ t³(log σ₀)³/σ₀² is no longer small, so the rate is **heuristic**.
  - The heuristic predicts a collision time of about 4σ₀/log σ₀, if the initial drift held all the way to the line.
  - The certified bounds Λ > 24.5, 30, 44 of claims #3, #5, #6, at σ₀ = 12.863, 18.541, 32.713, give Λ·log σ₀/σ₀ > 4.9, 4.7, 4.7.
- **A numerical check of Lemma L.** Step 4 predicts the initial drift |d Im z/dt| = Re φ′(s₁) − ½ log 2, since Im z = −(2 Re s − 1).
  - With s₁ = σ₀ + i Im ρ at the three deep zeros ρ of claims #3, #5, #6 (σ₀ = 12.863, 18.541, 32.713), mpmath gives Re φ′(s₁) = 1.36101, 1.52873, 1.80274, so the predicted drifts are 1.0144, 1.1822, 1.4562 (`code/02-functional-equation-family/drift_check.py`).
  - The heat-flow tracks move Im z at the average rates 1.0134 (from −24.7348 to −24.2281 over t ∈ [0, 0.5]), 1.1817 (−36.0825 → −34.9008 over [0, 1]) and 1.4560 (−64.4270 → −62.9710 over [0, 1]) (`runs/x7_track_deep_*.txt`).
  - So Lemma L predicts the initial heat-flow drift to 0.1%.
- **What the proof uses.** Only three things about the family: (i) real self-dual coefficients (for the LP/heat-flow argument); (ii) a(1) → 0 with a(2) ≠ 0; (iii) boundedness of the other coefficients. No Euler product and no information about ζ enter. The result says that the functional equation of L(s, χ₋₇) puts no upper bound on Λ.
- **Points a referee should check.**
  - (a) The uniformity of Stirling on S. Im s is bounded, so this is standard.
  - (b) Treating complex φ′ in the Gaussian integral. The Gaussian integral of exp(αw + βw²) is exact for complex α and Re β < 1/(4t), and the absolute-value remark of Lemma L handles the remainders.
  - (c) The reflected part, where Re(s + w/2) < ½.
  - (d) The choice of branch for u\*.
