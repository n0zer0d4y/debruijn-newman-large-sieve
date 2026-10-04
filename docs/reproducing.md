# Reproducing the results

## Environment

- **Python 3.12** with the versions in `requirements.txt` (python-flint 0.8.0, mpmath 1.3.0, numpy 2.4.2, scipy 1.17.1). These are the exact versions used for the certificates. `code/01-prime-moduli-large-sieve/claim1_check.py` also needs sympy.
- **Lean** v4.28.0 with Mathlib v4.28.0, pinned in `lean/lean-toolchain` and `lean/lakefile.toml`.
- **Windows note.** The scripts were run on Windows with CPython. They use only relative imports (`code/shared-dirichlet-pipeline`, `code/03-davenport-heilbronn`), so the platform should not matter.

## Result 1: the prime-moduli large sieve

```
cd lean
lake exe cache get        # downloads prebuilt Mathlib
lake build
lake env lean AxiomCheck.lean
```

`AxiomCheck.lean` prints the axioms each theorem depends on. Every one should be among `propext`, `Classical.choice` and `Quot.sound`, with no `sorryAx`. `scripts/setup-shared-deps.ps1` can instead point `.lake/packages` at a shared Mathlib store on Windows.

The numerical illustration of main.md Remark 3 comes from `python code/01-prime-moduli-large-sieve/claim1_check.py` (output in `code/01-prime-moduli-large-sieve/runs/claim1_check.txt`).

## Result 2: the χ₋₇ functional-equation family

All scripts are in `code/02-functional-equation-family/`. The certificate file names encode θ and t: `rouche-theta<num>_<den>-t<100·t>.json`.

| Certificate               | Command                                                                                                                                                                 |
| ------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Λ(f\_{1/10}) > 0.30       | `python x7_cert.py 1 10 0.30 102.08497690108965 -0.19249494975385603 r=2e-5 delta=1e-5`                                                                                 |
| Λ(f\_{1/20}) > 0.12       | `python x7_cert.py 1 20 0.12 102.20308 -0.18882 r=1e-5 delta=1e-6`                                                                                                      |
| Λ(f\_{3/100}) > 0.06      | `python x7_cert.py 3 100 0.06 102.24162 -0.08561 r=1e-5 delta=1e-6`                                                                                                     |
| Λ(f\_{1/50}) > 0.02       | `python x7_cert.py 1 50 0.02 102.26686 -0.11072 r=1e-5 delta=1e-6`                                                                                                      |
| Λ(f\_{704/1000}) > 24.5   | `python x7_cert.py 704 1000 24.5 4.47195 -0.92827 hp=1 R=0.6 r=1e-6 delta=1e-12 newton=1 W=120 rt=50 tail=2 cpu=20`                                                     |
| Λ(f\_{70406/100000}) > 30 | `python x7_cert.py 70406 100000 30 4.427456398001064 -1.75066014263749 hp=1 R=1.0 r=1e-6 delta=1e-12 newton=0 W=150 rt=50 tail=2 cpu=20`                                |
| Λ(f\_θ) > 44, a(1) = 1e-10 | `python x7_cert.py 70406119530951995496 100000000000000000000 44 6.101087845997802 -1.038177241043204 hp=1 R=0.9 r=1e-6 delta=1e-12 newton=1 W=220 rt=50 tail=2 cpu=15` |

- **Run times.** The first four take about 30–50 seconds each; the last three about 3–6 minutes each.
- **Supporting numerics,** with outputs in `runs/`:
  - `deep_zeros.py` (the deep zeros near θ\*);
  - `x7_track.py` (heat-flow tracks);
  - `chi7_pair.py` and `pair_law.py` (the pair-collision law);
  - `theta_star.py`;
  - `drift_check.py`;
  - `tail_check.py` (the Stirling tail bound against a direct numerical tail).

  The exact arguments are in `code/02-functional-equation-family/notes.md` and at the top of each output file.

## Result 3: Davenport–Heilbronn

Scripts are in `code/03-davenport-heilbronn/`, with the certificates in `certificates/03-davenport-heilbronn/`.

- **Lower bound.** `python dh_rouche_run.py` runs `cert_dh2.py` at t = 0.22 and 0.23 (the Rouché certificates; about 40 and 70 seconds).
- **The zero-free strip.** `python dh_sigma0_upper.py`.
- **Upper bound (conditional).** `dh_upper.py` handles the final region and error terms; `dh_rect.py` and `dh_rect_summary.py` handle the rectangle winding certificate.
  - The rectangle runs for the 0.58 bound are in `rect_runs_050/`.
  - Run from `code/03-davenport-heilbronn/`, appending `t0=0.5 y0=0.4 X=50000` to each `dh_rect.py` command. It takes about 10 minutes on 2 processes.
  - `python dh_rect.py direct bottom 0 704`, `direct top 0 404` and `direct left 0.4 -1`.
  - `python dh_rect.py rs bottom 700 25000`, `rs bottom 25000 50000`, `rs top 400 25000` and `rs top 25000 50000`. Each prints one JSON line; save it as a file in `rect_runs_050/`.
  - Then `python dh_rect_summary.py`, which re-runs part (a) and writes `lambda-058.json`.
- **Supporting numerics,** with outputs in `runs/`:
  - `dh_explore.py`;
  - `dh_track_zero240.py`;
  - `dh_track_zero85_gh.py`;
  - `rs_identity_dh.py`;
  - `c0_census.py`;
  - `dh_family_dim.py`.
- **Details.** `notes.md` in that folder lists every script.
