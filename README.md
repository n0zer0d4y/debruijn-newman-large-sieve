# debruijn-newman-large-sieve

Three results in analytic number theory, on the large sieve with prime moduli and on the de Bruijn–Newman constant Λ. Each comes with a written argument, interval-arithmetic certificates where computation is involved, and reproduction code. The finite core of the first result is also formalised in a small Lean 4 library.

> **Status: unrefereed.** These results were produced by an autonomous multi-agent research system (see [How these results were produced](#how-these-results-were-produced)). No human expert has reviewed them yet; such review is very welcome. Literature searches (September–October 2026) found none of these results in print. That is evidence of novelty, not proof of it; see [Prior art](#prior-art).

## The three results

### 1. The log log Q in the prime-moduli large sieve is necessary

[`papers/01-prime-moduli-large-sieve/`](papers/01-prime-moduli-large-sieve/main.md)

- **The question.** Wolke proved the large sieve inequality with prime moduli up to an extra factor log log Q. Iwaniec (Rev. Mat. Iberoam. 38 (2022), 2337–2354) writes: "We do not know if the bound (1.3) holds true without the factor log log Q for general coefficients."
- **The answer.** It cannot be removed. For N = Q^{2−2θ} with 0 < θ < 1/2, the optimal constant is at least (2θe^γ/π² + o(1)) log log Q · Q²/log(Q/√N).
- **The method.** A "primorial spoke" construction: Farey fractions with prime denominators accumulate near b/k for a primorial k. The method is that of Baier–Lynch–Zhao (square moduli) and Corrigan (quadratic moduli). The new point is its application to prime moduli.
- **Rigor.** A written, elementary proof. Its finite arithmetic core (three lemmas) is machine-checked in Lean 4 ([`lean/`](lean/)). The analytic steps (Jordan's inequality, the prime number theorem, Mertens' theorem) are not formalised.

### 2. The de Bruijn–Newman constant is unbounded on a family with one fixed functional equation

[`papers/02-lambda-unbounded-functional-equation-family/`](papers/02-lambda-unbounded-functional-equation-family/main.md)

- **The family.** The real Dirichlet series with exactly the functional equation of L(s, χ₋₇) form a circle f_θ, with the L-function itself at θ = 0.
- **Theorem U** (written proof, not formalised). Λ(f_θ) → ∞ as θ → θ\*, the angle at which the first Dirichlet coefficient vanishes. So a functional equation alone places no upper bound on Λ.
- **Certified,** by interval-arithmetic Rouché certificates:
  - Λ > 24.5, 30 and 44 as a_θ(1) → 0;
  - near the L-function, Λ > 0.30, 0.12, 0.06 and 0.02 at θ = 1/10, 1/20, 3/100 and 1/50.
- **Conditional:** Λ is upper semicontinuous at the L-function. The proof is complete except for one lemma (uniform far-field dominance), which is only sketched.
- **Numerical:** a two-zero "pair-collision law" predicts, with no fitted parameter, where the off-line zeros near the L-function appear.

### 3. Bounds for the de Bruijn–Newman constant of the Davenport–Heilbronn function

[`papers/03-davenport-heilbronn-lambda/`](papers/03-davenport-heilbronn-lambda/main.md)

- **Why this function.** The Davenport–Heilbronn function satisfies a Riemann-type functional equation but has zeros off the critical line.
- **Certified:** Λ_DH > 0.23, by an interval-arithmetic Rouché certificate for a non-real zero of the heat-flowed function.
- **Conditional:** Λ_DH ≤ 0.58.
  - The method is the Polymath effective heat-flow approximation, transferred to this function through a Riemann–Siegel-type identity ([`papers/supporting/riemann-siegel-identity/`](papers/supporting/riemann-siegel-identity/)), plus a rectangle winding certificate.
  - The computation is interval arithmetic, but several error-term lemmas are not yet written out in full (open items O1–O4 of [`derivation.md`](papers/supporting/riemann-siegel-identity/derivation.md)).
- **Context.** Dobner (Acta Arith. 201 (2021)) proved Λ ≥ 0 for the extended Selberg class, which contains these functions. He gave no upper bounds or numerical values.

## How these results were produced

The research was carried out by an autonomous multi-agent research system. The computations and proofs use standard tools (python-flint/Arb, mpmath, Lean 4 with Mathlib). What the system adds is discipline about evidence:

- **Every number is tied to a logged computation.** Before a result is reported, its numbers, confidence and citations are checked against a log of the computations and literature lookups that produced them, and anything unsupported is removed.
- **Every claim carries a rigor grade** (see [`docs/verification.md`](docs/verification.md)): proved in Lean, certified by interval arithmetic, cited, or conjectural.
- **Key results are cross-checked independently.** Two agents worked on the problems; the written proofs were reviewed adversarially by the second agent, and the certificates were re-run with corrected code before they were kept.
- **Work persists across sessions** through a research memory that links each fact to the evidence behind it.

A human set the research goals and decided what to publish; the system chose the specific problems and carried out the work.

## Repository layout

```
papers/                       one folder per result
  01-.../ 02-.../ 03-.../     main.md (the paper), claims.md (each claim with its evidence and rigor grade), proofs
  supporting/                 the Riemann–Siegel-type identity and the effective approximation used by results 2 and 3
certificates/                 interval-arithmetic certificate files (JSON)
code/                         the Python scripts that produce the certificates, with the outputs of our runs
  shared-dirichlet-pipeline/  modules shared by the result-3 code
lean/                         Lake package DebruijnNewmanLargeSieve (Lean 4 + Mathlib)
docs/                         reproducing.md, verification.md
```

## Reproducing

[`docs/reproducing.md`](docs/reproducing.md) lists the exact command for every certificate. In brief:

```
pip install -r requirements.txt
cd code/02-functional-equation-family
python x7_cert.py 1 10 0.30 102.08497690108965 -0.19249494975385603 r=2e-5 delta=1e-5
```

The last command re-derives the certificate for Λ(f\_{1/10}) > 0.30 in about 30 seconds. For the Lean library:

```
cd lean && lake exe cache get && lake build && lake env lean AxiomCheck.lean
```

## How much to trust each claim

Every claim carries a rigor grade in its paper's `claims.md`; see [`docs/verification.md`](docs/verification.md).

- **proved:** a Lean theorem, kernel-checked, depending only on `propext`, `Classical.choice` and `Quot.sound`.
- **certificate / numeric-experiment:** rigorous interval (Arb) arithmetic, reproducible from `code/`.
- **reference:** a published result.
- **conjectural:** a written proof that is not formalised, or a claim resting on lemmas that are only sketched. Theorem U, the semicontinuity theorem and the bound Λ_DH ≤ 0.58 are in this class.

## Prior art

Read:

- H. Iwaniec, Rev. Mat. Iberoam. 38 (2022), doi:10.4171/RMI/1381;
- D. Wolke, Mém. Soc. Math. France 25 (1971), doi:10.24033/msmf.52;
- S. Baier, S. B. Lynch, L. Zhao, Bull. Aust. Math. Soc. 100 (2019), doi:10.1017/S0004972719000224;
- S. Baier, arXiv:2503.18009;
- C. C. Corrigan, Acta Arith. 222 (2026), doi:10.4064/aa240923-21-1. It gives the same method for quadratic moduli and no prime-moduli lower bound;
- A. Dobner, Acta Arith. 201 (2021), doi:10.4064/aa200603-23-7;
- B. Rodgers, T. Tao, Forum Math. Pi 8 (2020), doi:10.1017/fmp.2020.6;
- D.H.J. Polymath, Res. Math. Sci. 6 (2019), doi:10.1007/s40687-019-0193-1;
- J. Andrade, A. Chang, S. J. Miller, J. Number Theory 144 (2014), doi:10.1016/j.jnt.2014.04.021;
- C. M. Newman, W. Wu, Bull. Amer. Math. Soc. 57 (2020), doi:10.1090/bull/1668.

Not yet read (closed access):

- D. Wolke, "Farey-Brüche mit primem Nenner und das große Sieb", Math. Z. 114 (1970), 145–158, doi:10.1007/BF01110322;
- D. Wolke, "On the large sieve with primes", Acta Math. Acad. Sci. Hungar. 22 (1971), 239–247, doi:10.1007/BF01896016;
- D. Wolke, "A lower bound for the large sieve inequality", Bull. London Math. Soc. 6 (1974), 315–318, doi:10.1112/blms/6.3.315.

The last could bear on result 1.

## License

Code (`code/`, `lean/`): MIT, see [`LICENSE`](LICENSE). Papers and documentation: CC BY 4.0, see [`LICENSE-docs`](LICENSE-docs).

## Citation

See [`CITATION.cff`](CITATION.cff).
