# Λ across the q = 7 odd Theorem A′ family (functional-equation twins of L(s, χ₋₇)): code notes

All runs use Python with python-flint and mpmath, in a single process. Outputs stored in `runs/` are named in the table; "not stored" marks exploratory runs whose output is not kept here.

| file            | what                                                                                        | result                                                                       |
| --------------- | ------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------- |
| x7\_cert.py      | the exact family (Gauss-sum b, Arb balls) and the Rouché certificates                       | `certificates/02-functional-equation-family/`                                |
| x7\_scan.py      | off-line zero scan of f\_θ (grid on \|f\|, mpmath Newton)                                      | not stored                                                                   |
| x7\_track.py     | heat-flow track of one zero (Gauss–Hermite, 28 nodes)                                       | `runs/x7_track_theta*.txt` (θ = 0.1, 0.05, 0.03, 0.02); `runs/x7_track_deep_*.txt` (deep zeros, t ≤ 1) |
| deep\_zeros.py   | the deep zeros near θ\* from the two-term balance                                           | `runs/deep_zeros_*.txt`                                                      |
| theta\_star.py   | θ\* in Arb, and the rational θ with a prescribed a\_θ(1)                                     | `runs/theta_star.txt`                                                        |
| chi7\_pair.py    | zeros of L(s, χ₋₇) near height 51                                                           | `runs/chi7_pair_47_55.txt`                                                   |
| pair\_law.py     | the pair-collision law                                                                      | `runs/pair_law.txt`                                                          |
| pair\_scan.py    | pair statistics along the critical line                                                     | not stored                                                                   |
| drift\_check.py  | the initial drift predicted by Lemma L against the tracks                                   | `runs/drift_check.txt`                                                       |
| tail\_check.py   | the Stirling tail bound (TAIL = 2) against a numerical tail                                 | `runs/tail_check_t2450.txt`                                                  |
| cert\_profile.py | timing profile of the certificate                                                           | not stored                                                                   |

**Notes.**

- **The Rouché radius.** It must satisfy r < |d|R²/M. The majorant M = ∫Gaussian·sup|ξ| ignores cancellation and grows quickly with t, so at larger t the certificate keeps the Newton centre at Arb precision (high-precision mode) to make |h₀| small enough.
- **Which zeros survive.** Off-line zeros at larger height move to the axis faster, so the longest-surviving zero near the L-function is a low one.
