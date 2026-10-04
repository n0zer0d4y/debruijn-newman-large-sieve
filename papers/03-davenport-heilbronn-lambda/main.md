---
title: "Bounds for the de Bruijn–Newman constant of the Davenport–Heilbronn function"
authors: ["Nozero Day"]
note: "Produced by an autonomous multi-agent research system; see the repository README. Unrefereed."
status: draft
date: "2026-10-04"
---

# Bounds for the de Bruijn–Newman constant of the Davenport–Heilbronn function

**Status: draft.** The computations use Arb ball arithmetic (python-flint) and are graded as certificates in `claims.md`. The upper bound is also conditional on error-term lemmas that are only sketched (§4.3).

## Abstract

Let f(s) = Σ a\_n n^{−s} be the Davenport–Heilbronn function. Its coefficients are a\_n = (1, κ, −κ, −1, 0) for n ≡ 1, …, 5 mod 5, with κ = (√(10 − 2√5) − 2)/(√5 − 1).

Let ξ(s) = (5/π)^{(s+1)/2}Γ((s + 1)/2)f(s), so that ξ(s) = ξ(1 − s), and let H₀(z) := ξ((1 + iz)/2), which is real entire and even. Let Λ\_DH be the de Bruijn–Newman constant of the heat flow H\_t.

- **Certified:** Λ\_DH > 0.23, by a Rouché certificate for a non-real zero of H\_{0.23}.
- **Conditional:** Λ\_DH ≤ 0.58.
  - The method is a Polymath-15-type effective approximation plus a rectangle winding argument at t₀ = 0.5, y₀ = 0.4.
  - It uses the Riemann–Siegel integral identity, which we show transfers to these non-multiplicative coefficients.
  - The computation is interval arithmetic, but it rests on error-term lemmas that are not yet written out in full (§4.3).

We know of no earlier quantitative bound for the de Bruijn–Newman constant of a natural L-type function whose constant is strictly positive.

## 1. Set-up

- **Coefficients.**
  - f has real coefficients and equals ((1 − iκ)/2)L(s, χ) + ((1 + iκ)/2)L(s, χ̄), for the character χ mod 5 with χ(2) = i.
  - It is not in the Selberg class, since it has no Euler product.
  - Davenport and Heilbronn showed that it has zeros off the critical line, and Spira [Sp94] computed some of them. The lowest, at 0.808517 + 85.699348i, is recomputed in `code/03-davenport-heilbronn/dh_explore.py` (`runs/dh_explore.txt`).
  - Balanzario and Sánchez-Ortiz [BS07] computed further off-line zeros by deforming (1 + 5^{1/2−s})ζ(s) into f.
- **The heat flow.** H\_t(z) = ∫Φ(u)e^{tu²}e^{izu}du, with Φ(u) ∝ e^{3u}θ(e^{4u}) and θ(y) = Σ a\_n n e^{−πn²y/5}.
- **The evaluation formula.** Since e^{tu²} = (4πt)^{−1/2}∫e^{−w²/(4t)}e^{wu}dw,

  H\_t(z) = (4πt)^{−1/2} ∫ e^{−w²/(4t)} ξ((1 + w + iz)/2) dw.

  All certified evaluations use this Gaussian average of ξ.

## 2. Why a non-real zero bounds Λ from below

H\_t is real entire of order 1 and genus ≤ 1. For t ≥ Λ it lies in the Laguerre–Pólya class, and the backward heat operator e^{−s d²/dz²} with s ≥ 0 preserves that class (Pólya; de Bruijn). This uses no positivity of Φ, and Φ indeed changes sign here. So {t : H\_t ∈ LP} = [Λ, ∞), and a non-real zero of H\_{t₁} gives Λ > t₁.

## 3. The lower bound (Rouché)

At t = 0.23, let c be a Newton approximation of the zero descending from the off-line zero 0.86953 + 240.40467i (the track in `runs/dh_track_zero240.txt`); Im c ≈ −0.035. The quantities below are rigorous Arb enclosures.

- h₀ = H\_t(c), a point evaluation: `acb.integral` plus the analytic tail bound of `cert_dh.DH.tail_bound`.
- M ≥ sup\_{|z−c|≤R}|H\_t|, from the modulus majorant (4πt)^{−1/2}∫e^{−w²/(4t)}sup|ξ| dw. The w-axis is cut into 400 cells, each with an Arb ball bound of ξ and an erf Gaussian mass.
- d, a central difference with step δ, with |H′(c) − d| ≤ M₃δ²/6 and M₃ ≤ 6M/R³.
- M₂ ≤ 2M/(R − r)².

On |z − c| = r, with p(z) = h₀ + d(z − c), |H − p| ≤ |H′(c) − d|r + M₂r²/2, which is < |d|r − |h₀| ≤ |p|. So H\_t has exactly one zero in |z − c| < r. Since r = 10⁻⁵ < |Im c|, that zero is non-real.

| t    | Im c    | r      | \|h₀\|      | \|d\|  | \|H − p\| bound | \|p\| lower bound |
| ---- | ------- | ------ | ----------- | ------ | --------------- | ----------------- |
| 0.22 | −0.1462 | 4·10⁻⁵ | ≤ 1.2·10⁻¹³ | ≥ 7.06 | 1.37·10⁻⁴       | 2.82·10⁻⁴         |
| 0.23 | −0.0349 | 10⁻⁵   | ≤ 9.2·10⁻¹⁴ | ≥ 1.71 | 9.9·10⁻⁶        | 1.71·10⁻⁵         |

(R = 0.05; M ≤ 215 and M ≤ 248 in the two runs.)

Numerically this zero reaches the real axis for t between 0.2300 and 0.2310 (`runs/dh_track_zero240.txt`). The first off-line zero, at height 85.70, reaches it for t between 0.176 and 0.178 (`runs/dh_track_zero85_gh.txt`).

## 4. The upper bound

### 4.1 The zero-free strip: Λ\_DH ≤ 1.6026

For σ > 1, |f(s)| ≥ 1 − Σ\_{n≥2}|a\_n|n^{−σ}, and

Σ\_{n≥2}|a\_n|n^{−σ} = 5^{−σ}[ζ(σ, ⅕) + ζ(σ, ⅘) + κ(ζ(σ, ⅖) + ζ(σ, ⅗))] − 1.

This sum is decreasing in σ, and Arb bisection shows it is < 1 for σ ≥ σ₀ = 1.3951362. By ξ(s) = ξ(1 − s), every zero of H₀ lies in |Im z| ≤ Δ = 2σ₀ − 1 = 1.7902723. De Bruijn's strip theorem then gives Λ\_DH ≤ Δ²/2 = 1.60254 (`certificates/03-davenport-heilbronn/upper-bound.json`).

### 4.2 Λ\_DH ≤ 0.58 (conditional)

- **Riemann–Siegel for f.**
  - The integral identity (Theorem A′ of papers/supporting/riemann-siegel-identity/proof-rs-general.md) uses the coefficients only through periodicity, parity, a₀ = 0, Σ aₙ = 0, real values, and the Gauss-sum self-duality Σ\_m a\_m e^{2πimn/5} = i√5·aₙ.
  - The self-duality holds exactly for f. With u = √(10 − 2√5) and v = √(10 + 2√5), it reduces to 2u = (√5 − 1)v and 2v = (√5 + 1)u. Squared, both sides equal 40 ∓ 8√5, and both sides are positive.
  - A numerical check of the identity itself (`rs_identity_dh.py`; `runs/rs_identity_dh.txt`) gives relative residuals ≤ 3.7·10⁻²².
  - So H\_t/B\_t = g\_N + O(e\_C + e\_A + e\_B), with g\_N = Σ\_{n≤N} aₙ bₙ n^{−s\*} + γ Σ aₙ bₙ n^{−w} and no Euler mollifier (derivation.md §2).
- **The reduction.**
  - §4.1 and de Bruijn put every zero of H\_{t₀} in |Im z| ≤ Ytop = (Δ² − 2t₀)^{1/2} = 1.48495.
  - H is real and even, so it suffices to exclude zeros with x ≥ 0 and y₀ ≤ y ≤ Ytop.
  - Then every zero has |Im z| < y₀, and de Bruijn gives Λ ≤ t₀ + y₀²/2 = 0.58.
  - No barrier, numerical verification of GRH or positivity of Φ is needed.
- **(a) The final region x ≥ X = 5·10⁴.** Plain blocks with weights |aₙ|, up to N = 10⁶, and a tail bound valid for every y (derivation.md §3, Lemma F) give |G − 1| < 1 for y₀ ≤ y ≤ Ytop, where G = H\_{t₀}/B\_{t₀}. The minimum lower bound is 0.0604, against an error bound of 0.0272 evaluated from x = X.
- **(b) The rectangle R = [0, X] × [y₀, Ytop].**
  - **Why the boundary suffices.** B\_t is analytic and nonvanishing there, so it suffices that G avoids (−∞, 0] on ∂R. Then G has winding number 0, and H\_t has no zeros in R.
  - **How each piece is certified.** Every boundary piece carries a 1-D Chebyshev certificate: a rigorous majorant on the Bernstein ellipse, node enclosures, and a grid bound (derivation.md §3).
  - **RS pieces.** On the bottom edge for x ≥ 700 and the top edge for x ≥ 400, G is bounded by the RS approximation. Lemma I is used at fixed c per block, with K\_h = 1.179 from Lemma H. Below these points Lemma I is not usable here: its tail condition fails at c = 10.5, and near x = 640 the bottom-edge error exceeds the margin.
  - **Direct pieces.** Elsewhere, including the side x = 0, G is evaluated directly. The heat integral is computed by the trapezoid rule with the Trefethen–Weideman strip bound, from exact ξ values.
  - **Size and margins.** 50,359 pieces; worst margin 0.1025, on an RS piece of the bottom edge.
  - **Certificate:** `certificates/03-davenport-heilbronn/lambda-058.json`, with runs in `code/03-davenport-heilbronn/rect_runs_050/`.
- **Rigor.** The enclosures are Arb ball arithmetic; the exact identities go through mpmath. The error terms rest on the lemmas listed in §4.3.

### 4.3 Why the upper bound is conditional

The effective approximation of derivation.md has four open items (O1–O4 there):

- Prop 6.1^χ is asserted to follow Polymath's proof line by line, but is not written out.
- Lemma T, the Gaussian tail of e\_C, is proved only for T′ ≥ 10⁶, while the RS pieces use it from T′ ≈ 200 (x ≈ 400). The code bounds that region by 10⁻³⁰⁰. Its true size is far smaller, but this is not proved.
- Lemma M's constant is asserted, not derived, and the lemma is written for y ≤ 1, whereas the final region uses y up to Ytop = 1.485.
- Lemma I is written in full for q = 4 with even c ≥ 22. The q = 5, half-integer c ≥ 10.5 form used here is given only in compact form.

The lower bound 0.23 is unaffected: it uses the Rouché route and none of these lemmas.

## 5. Ceiling and next steps

- f is not an L-function in the Selberg class, and this result says nothing about ζ or GRH. It measures, for a natural function where RH fails, how much heat flow its zeros need.
- **Better lower bounds.** Certify off-line zeros with larger σ at greater height.
- **An unconditional upper bound.** Write out the four lemmas of §4.3.
- **A better upper bound.** A larger final-region X, together with a longer rectangle, would allow smaller (t₀, y₀). The cost is compute.

## References

- [BS07] E. P. Balanzario, J. Sánchez-Ortiz, Zeros of the Davenport–Heilbronn counterexample, Math. Comp. 76 (2007), 2045–2049, doi:10.1090/S0025-5718-07-01999-0.
- [Do21] A. Dobner, A proof of Newman's conjecture for the extended Selberg class, Acta Arith. 201 (2021), 29–62, doi:10.4064/aa200603-23-7.
- [P19] D.H.J. Polymath, Effective approximation of heat flow evolution of the Riemann ξ function, and a new upper bound for the de Bruijn–Newman constant, Res. Math. Sci. 6 (2019), Paper No. 31, doi:10.1007/s40687-019-0193-1, arXiv:1904.12438.
- [RT20] B. Rodgers, T. Tao, The de Bruijn–Newman constant is non-negative, Forum Math. Pi 8 (2020), e6, doi:10.1017/fmp.2020.6.
- [Sp94] R. Spira, Some zeros of the Titchmarsh counterexample, Math. Comp. 63 (1994), 747–748.
- [TW14] L. N. Trefethen, J. A. C. Weideman, The exponentially convergent trapezoidal rule, SIAM Rev. 56 (2014), 385–458, doi:10.1137/130932132.
