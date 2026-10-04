# Upper semicontinuity of Λ at the L-function

**Status: conjectural.** Steps 0 and 2–4 are complete. Step 1 rests on Lemma 1′, which is only sketched.

Lemma 1′ uses the effective approximation of papers/supporting/riemann-siegel-identity/derivation.md, whose Prop 6.1^χ is not written out (open items O1 and O5 there). For y > 1 only a sketch is given. A self-contained qualitative large-x asymptotic would remove the dependency; this theorem needs no explicit constants. The proof was reviewed adversarially by a second, independent agent of the research system.

**Setting.**
- q = 7. The coefficient vectors are a\_θ = cos θ·χ + sin θ·b, where χ = (·/7) and b is the exact odd DFT-self-dual partner (main.md §1). Put f\_θ = Σ a\_θ(n) n^{−s}.
- ξ\_θ(s) = (7/π)^{(s+1)/2} Γ((s+1)/2) f\_θ(s). H\_t^θ(z) = ∫ Φ\_θ(u) e^{tu²} e^{izu} du, equivalently the Gaussian average of ξ\_θ((1+iz)/2).
- Λ\_θ := Λ(f\_θ) and Λ₀ := Λ(L(s, χ₋₇)).

**Theorem.** limsup\_{θ→0} Λ\_θ ≤ Λ₀.

**Proof.** Fix t₁ > Λ₀ and y₀ ∈ (0, 1).

- **Step 0 (linearity).** Φ\_θ = cos θ·Φ\_χ + sin θ·Φ\_b, so H\_t^θ = cos θ·H\_t^χ + sin θ·H\_t^b for every t.

- **Step 1 (uniform far field).**
  - **The strip.**
    - For σ > 1, |f\_θ(s)| ≥ |a\_θ(1)| − Σ\_{n≥2}|a\_θ(n)| n^{−σ}. Hence all zeros of H\_0^θ lie in |Im z| ≤ Δ(θ), where Δ(θ) = 2σ₀(θ) − 1 and σ₀(θ) solves Σ\_{n≥2}|a\_θ(n)/a\_θ(1)| n^{−σ} = 1.
    - The ratios |a\_θ(n)/a\_θ(1)| are bounded for |θ| ≤ θ₁, since |a\_θ(n)| ≤ |cos θ| + |sin θ|·max|b| and a\_θ(1) → 1. So Δ(θ) ≤ Δ\* there. (At θ = 0.1 the coefficients satisfy |a\_θ(n)| ≤ 1.122 and a\_θ(1) = 0.877, so the ratios are at most 1.28.)
    - At θ = 0, σ₀ = 1.6883170 and Δ = 2.3766340.
    - De Bruijn's strip theorem then puts the zeros of H\_{t₁}^θ in |Im z| ≤ Y := (Δ\*² − 2t₁)₊^{1/2}, so Y ≤ Δ\* ≈ 2.38 for small θ₁.
  - **The final region.** Theorem A′ applies to every f\_θ (proof-rs-general.md, Remark 5). So the final-region inequality of the effective approximation (derivation.md §2.2) holds with weights w\_n := |a\_θ(n)/a\_θ(1)|:

    1 − Σ\_{n≥2} w\_n b\_n n^{−Re s\*} − |γ| Σ\_n w\_n b\_n n^{−Re w} − e\_C − e\_{AB} > 0 on x ≥ X, y ∈ [y₀, Y].

  - **Uniformity.** The left side is decreasing in each weight w\_n ≤ W(θ₁) := sup\_{|θ|≤θ₁} max\_n |a\_θ(n)/a\_θ(1)|. As x → ∞ it tends to 1, uniformly in the weights and in y ∈ [y₀, Y], since Re s\* − ½ ≥ y₀/2 + (t₁/4)·log(7x/4π) → ∞ and the error terms → 0 (Lemma 1′). So one X(t₁, y₀, θ₁) works for all |θ| ≤ θ₁.
  - **The right exponent.**
    - b\_n = e^{(t/4)log²n} grows, so the decay comes from Re s\* − (t/4)log n ≥ (1+y)/2 + (t/8)log(7x/4π) − O(1) for n ≤ N.
    - For the B-sum it comes from Re w − (t/4)log n ≥ (1−y)/2 + (t/8)log(7x/4π) − O(1), with |γ| ≪ (7x/4π)^{−y/2}. Both tend to +∞.
    - At y > 1 the B-exponent becomes positive only once log(7x/4π) > 4(y − 1)/t₁. With y = Y ≈ 2.25 and t₁ = 0.3 that is x > 3·10⁷, so X is very large. That is harmless for this qualitative theorem.
  - **Lemma 1′.** The error terms e\_C, e\_A, e\_B of the effective approximation tend to 0 as x → ∞, uniformly for y in the compact [y₀, Y] and |θ| ≤ θ₁.
    - **Why it is needed.** Y is not small (above), while Polymath's effective approximation [P19, Thm 1.3] and its χ-version (derivation.md §2) are stated with explicit constants only for 0 ≤ y ≤ 1.
    - **Why it should hold.**
      - The restriction serves the explicit constants. The qualitative decay follows from the same saddle-point analysis with Stirling's remainder on a compact y-range, and the dependence on θ is linear, through K\_h and the coefficient norms.
      - Polymath state the extension themselves for ζ. Theorem 1.3 and Corollary 1.4 hold "in slightly larger regions than (5), for instance if y is allowed to be as large as 10" [P19, §1].
      - In §9 they write: "Strictly speaking, the estimates there required y ≤ 1 rather than y ≤ 10; however, as remarked at the beginning of Section 6, all the estimates in that section would continue to hold under this weaker hypothesis if one adjusted all the numerical constants appropriately."
    - **Sketch for the Theorem A′ family.** Fix t₁ > Λ₀, a compact range 0 ≤ y ≤ Y and |θ| ≤ θ₁. Everything below is uniform in these.
      - (i) **e\_A + e\_B ≤ ε\_max·F** (derivation.md §2.3, stated for every real σ, so with no restriction on y).
        - ε\_max ≤ exp([(t²/8)(|α| + log N)² + t/4 + 1/6]/(T − 1.25)) − 1, with T = x/2 and |α| ≤ ½ log(7x/4π) + O(1). So ε\_max ≪ log²x/x → 0.
        - F is bounded: F ≤ W Σ\_{n≤N} b\_n n^{−Re s\*} + W|γ| Σ\_{n≤N} n^y b\_n n^{−Re s\* − Re κ}. With Re s\* ≥ (1+y)/2 + (t/2)(log N − O(1)) and L = log n ≤ log N, the first sum's exponent is (t/4)L² − (t/2)L log N − ((1+y)/2)L + O(tL) ≤ −(t/4)L² − ½L + O(tL).
        - The sums run over integers, i.e. over ∫e^{g(L)+L}dL. So the A-sum is at most about (4π/t)^{1/2}e^{(1−y)²/(4t)}, and the B-sum about (4π/t)^{1/2}e^{(1+Y)²/(4t)}, both independent of x.
        - So e\_A + e\_B → 0.
      - (ii) **e\_C → 0** (derivation.md §2.5).
        - e\_C ≤ Σ\_{σ∈{(1−y)/2, (1+y)/2}} e^{π²t/64} (|M₀(iT′)|/|B\_t|) K\_h·𝔼\_v[bracket].
        - By Lemma E(a) (derivation.md §3), which is valid for every y ≥ 0, log(|M₀(iT′)|/|B\_t|) ≤ −((1+y)/2)log w\_R − (t/4)(½log(7x/4π))² + O(1) → −∞.
        - K\_h ≤ K(q)·max\_n|a\_θ(n)| + |c₀(θ)|/2 (Lemma H), which is uniform in θ.
        - The Gaussian average of Lemma I's bracket is bounded uniformly for |σ| ≤ (1+Y)/2 + 12√t₁, since c ≈ w\_R → ∞; Lemma T covers the far tail.
      - (iii) **The main term tends to 1.** The n ≥ 2 A-terms and the whole B-sum tend to 0 once log(7x/4π) > 4(Y−1)/t₁, by the exponents displayed above.
      - **Conclusion.** The final-region inequality holds on x ≥ X(t₁, y₀, Y, θ₁), y ∈ [y₀, Y]. The explicit lemmas of derivation.md that assume y ≤ 1 (Lemma M) serve the certificates and are not used here.
    - **What is missing.** A written proof of the χ-version of Polymath's §6 (Prop 6.1^χ in particular) with o(1) errors, for y up to Y.
  - **Conclusion of Step 1.** H\_{t₁}^θ ≠ 0 on x ≥ X, y₀ ≤ y ≤ Y. Negative y follows because H is real, H(z̄) = conj H(z); x ≤ −X follows by evenness.

- **Step 2 (the base point).** t₁ > Λ₀, so every zero of H\_{t₁}^0 is real. It is also simple.
  - **Citation.** This is Csordas–Smith–Varga [CSV94, Cor. 1], as attributed by Rodgers–Tao [RT20]. CSV state it for ζ's Φ.
  - **A self-contained proof for any real even Φ of this type.** Suppose H\_t is real-rooted on [t₀, t₁] with t₀ ∈ (Λ₀, t₁), and that x₀ is a real zero of H\_{t₁} of multiplicity m ≥ 2.
    - Backward heat flow gives H\_{t₁−s} = e^{sD²}H\_{t₁}.
    - With x − x₀ = √s·w, we get H\_{t₁−s}(x₀ + √s w) = s^{m/2}(c·e^{D²}w^m + O(√s)), with c ≠ 0.
    - e^{D²}w^m is a rescaled Hermite polynomial at iw. Its zeros are purely imaginary, and at least two of them are non-zero.
    - By Hurwitz, H\_{t₁−s} then has non-real zeros for small s > 0, which contradicts real-rootedness at t₁ − s.

- **Step 3 (Hurwitz on the box B = [−X−1, X+1] × [−Y−1, Y+1]).**
  - **Setup.** Choose the vertical edges of B to avoid the zeros of H\_{t₁}^0; there are none on the horizontal edges, which lie off the axis. H\_{t₁}^θ → H\_{t₁}^0 uniformly on B as θ → 0, by Step 0 with |sin θ|·sup\_B |H\_{t₁}^b| → 0.
  - **Each real zero.** Take a real zero ρ of H\_{t₁}^0 and a disk D\_ρ of radius ε, with ε so small that the disks are disjoint, avoid the other zeros and lie in |y| < y₀. By Hurwitz, for small θ each D\_ρ contains exactly one zero of H\_{t₁}^θ. It is real, since H^θ is real on ℝ and a lone zero cannot have its conjugate partner inside D\_ρ.
  - **Everywhere else.** H\_{t₁}^0 is bounded away from 0 on B minus the disks, so for small θ there are no other zeros in B.

- **Step 4 (the strip lemma).** By Steps 1 and 3, for |θ| ≤ θ₂(t₁, y₀) every zero of H\_{t₁}^θ has |Im z| < y₀. So Λ\_θ ≤ t₁ + y₀²/2 (de Bruijn), and letting t₁ ↓ Λ₀ and y₀ ↓ 0 gives the theorem. ∎

**Corollaries.**
- Under GRH for χ₋₇ we have Λ₀ = 0, since Λ₀ ≥ 0 by Dobner's theorem for the extended Selberg class [Do21, Thm 2]. So Λ\_θ → 0.
- A proof that inf\_{θ≠0} Λ\_θ > 0 would disprove GRH(χ₋₇).

**Scope.** Nothing in the proof uses q = 7 beyond the existence of the family. The same argument gives upper semicontinuity of Λ at any point of any Theorem A′ family where a(1) ≠ 0.

**References.**
- [CSV94] G. Csordas, W. Smith, R. S. Varga, Lehmer pairs of zeros, the de Bruijn–Newman constant Λ, and the Riemann hypothesis, Constr. Approx. 10 (1994), 107–129, doi:10.1007/BF01205170.
- [Do21] A. Dobner, A proof of Newman's conjecture for the extended Selberg class, Acta Arith. 201 (2021), 29–62, doi:10.4064/aa200603-23-7.
- [P19] D.H.J. Polymath, Effective approximation of heat flow evolution of the Riemann ξ function, and a new upper bound for the de Bruijn–Newman constant, Res. Math. Sci. 6 (2019), Paper No. 31, doi:10.1007/s40687-019-0193-1, arXiv:1904.12438.
- [RT20] B. Rodgers, T. Tao, The de Bruijn–Newman constant is non-negative, Forum Math. Pi 8 (2020), e6, doi:10.1017/fmp.2020.6.
