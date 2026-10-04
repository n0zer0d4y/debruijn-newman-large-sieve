# Verification

## Rigor grades

Each paper folder has a `claims.md` with one row per claim: the exact statement, its rigor grade, the certificate or Lean file it rests on, and notes. Under `## Reasoning` it gives the argument for each row.

| Grade                            | Meaning                                                                                                              | How to check it                                                                           |
| -------------------------------- | -------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| proved                           | a Lean 4 theorem, kernel-checked, depending only on the standard axioms `propext`, `Classical.choice`, `Quot.sound`  | `cd lean && lake build && lake env lean AxiomCheck.lean`                                  |
| certificate / numeric-experiment | rigorous interval (Arb ball) arithmetic via python-flint                                                             | re-run the script listed in `reproducing.md`; the output must match the certificate file  |
| reference                        | a published result                                                                                                   | the cited paper                                                                           |
| conjectural                      | a written proof that is not formalised, or a claim resting on lemmas that are only sketched                          | read the proof file; independent review is welcome                                        |
| withdrawn                        | not claimed; the row is kept only so that the numbering is stable                                                    | —                                                                                         |

**On the certificates.** They are computed with Arb ball arithmetic, so every number carries a guaranteed error radius, and the tail bounds they use are proved.

The upper bound Λ\_DH ≤ 0.58 also uses analytic error terms whose lemmas are only partly written: open items O1–O4 of `papers/supporting/riemann-siegel-identity/derivation.md`. That claim is therefore graded conjectural.

The certificates have been re-run on the machine that produced them, not on independent hardware, and are not formally verified.

## What the certificates prove

- **Rouché certificates** (results 2 and 3, lower bounds; files named `rouche-*.json`).
  - Each states that a heat-flowed function H\_t has exactly one zero in a small disc, using the test |H − p| < |p| on the circle for the linear model p.
  - The disc does not meet the real axis, so the zero is non-real.
  - H\_t is real entire and even, of order 1. For t ≥ Λ it lies in the Laguerre–Pólya class, and evolving to larger t preserves real-rootedness (de Bruijn). So a non-real zero at time t gives Λ > t.
- **The upper-bound certificate** (result 3; conditional). It is a fixed-time argument, not Polymath's barrier criterion.
  - de Bruijn's strip theorem puts every zero of H\_{t₀} in |Im z| ≤ Ytop.
  - The final region excludes zeros with x ≥ X and y₀ ≤ y ≤ Ytop.
  - A winding certificate shows that G = H\_{t₀}/B\_{t₀} avoids (−∞, 0] on the boundary of [0, X] × [y₀, Ytop], so there are no zeros inside.
  - By symmetry every zero then has |Im z| < y₀, and de Bruijn gives Λ ≤ t₀ + y₀²/2.
  - No barrier, numerical verification of GRH or positivity of Φ is used.
  - The parameters and margins are in `certificates/03-davenport-heilbronn/lambda-058.json` and `code/03-davenport-heilbronn/notes.md`.

## Review

- **The written proofs** (result 1, Theorem U, the semicontinuity theorem and the Riemann–Siegel-type identity) were each reviewed adversarially by a second, independent agent of the research system.
- **The certificates** were re-run from the code in this repository, and their outputs matched the recorded files.
- **Human review.** No human expert has refereed these results yet.
