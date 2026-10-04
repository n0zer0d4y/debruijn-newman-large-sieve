"""Reconnaissance (numeric-experiment): the q = 7, odd-parity Theorem A' family.
Basis: the real character chi_{-7} = Legendre (n/7) and the orthogonal real DFT-self-dual odd vector b
(eigenvalue i sqrt 7), from the SVD null space as in ../03-davenport-heilbronn/dh_family_dim.py.
For f_theta = cos(theta) chi_{-7} + sin(theta) b (as Dirichlet series sum a(n) n^{-s} via Hurwitz zeta),
scan |f| on sigma in [0.52, 1.1], t in [2, TMAX], refine local minima by Newton, and report zeros off the
critical line (|sigma - 1/2| > 1e-6).  mpmath dps 20.  Usage: python x7_scan.py theta TMAX [TMIN] [SIGMA_MAX]
"""
import sys, math
import numpy as np
import mpmath as mp
from flint import acb, arb, ctx
ctx.prec = 64

mp.mp.dps = 20
q = 7


def basis():
    n = np.arange(q)
    F = np.exp(2j * np.pi * np.outer(n, n) / q)
    M = F - 1j * math.sqrt(q) * np.eye(q)
    P = np.zeros((q, q))
    for k in range(q):
        P[k, k] += 1; P[k, (-k) % q] += 1                 # odd: a(k) + a(-k) = 0
    e0 = np.zeros((1, q)); e0[0, 0] = 1
    A = np.vstack([M.real, M.imag, P, e0, np.ones((1, q))])
    _, s, vt = np.linalg.svd(A)
    null = vt[-2:]                                        # dimension 2 (dh_family_dim.py)
    chi = np.array([0, 1, 1, -1, 1, -1, -1], dtype=float)  # (n/7)
    chi_n = chi / np.linalg.norm(chi)
    assert np.linalg.norm(null.T @ (null @ chi_n) - chi_n) < 1e-9, 'chi_{-7} not in the family'
    b = null[0] - (null[0] @ chi_n) * chi_n
    if np.linalg.norm(b) < 1e-6:
        b = null[1] - (null[1] @ chi_n) * chi_n
    b = b / np.linalg.norm(b) * np.linalg.norm(chi)
    return chi, b


def make_f(theta):
    chi, b = basis()
    a = math.cos(theta) * chi + math.sin(theta) * b
    A = [mp.mpf(float(x)) for x in a]
    def f(s):
        return sum(A[r] * mp.zeta(s, mp.mpf(r) / q) for r in range(1, q)) * mp.mpf(q) ** (-s)
    return f, a


if __name__ == '__main__':
    theta = float(sys.argv[1]); TMAX = float(sys.argv[2]) if len(sys.argv) > 2 else 150
    TMIN = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
    f, a = make_f(theta)
    print('theta', theta, 'coefficients', np.round(a, 6).tolist(), flush=True)
    SMAX = float(sys.argv[4]) if len(sys.argv) > 4 else 1.1
    sig = np.linspace(0.52, SMAX, max(15, int(25 * (SMAX - 0.5)))); ts = np.arange(TMIN, TMAX, 0.1)
    AF = [arb(float(x)) for x in a]
    def fa(s):
        v = sum((AF[r] * s.zeta(acb(arb(r) / q)) for r in range(1, q)), acb(0)) * acb(q) ** (-s)
        return abs(complex(v.real.mid(), v.imag.mid()))
    G = np.array([[fa(acb(s_, t_)) for s_ in sig] for t_ in ts])
    found = []
    for i in range(1, len(ts) - 1):
        for j in range(1, len(sig) - 1):
            v = G[i, j]
            if v < G[i - 1:i + 2, j - 1:j + 2].min() + 1e-15 and v < 0.3:
                try:
                    z = mp.findroot(f, mp.mpc(sig[j], ts[i]))
                except Exception:
                    continue
                if abs(f(z)) < 1e-10 and abs(z.real - 0.5) > 1e-6 and 0.5 < z.real < SMAX + 0.1:
                    key = (round(float(z.real), 6), round(float(z.imag), 5))
                    if key not in found:
                        found.append(key)
    print('off-line zeros (sigma, t):', sorted(found, key=lambda x: x[1]), flush=True)
