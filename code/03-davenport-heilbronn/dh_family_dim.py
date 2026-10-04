"""Theorem A' family (proof-rs-general.md Remark 5), numeric-experiment: for each modulus q and parity a,
the dimension of the real space of q-periodic sequences with a(-n) = (-1)^a a(n), a(0) = 0, sum a = 0 and
DFT self-duality  sum_m a(m) e^{2 pi i m n/q} = i^a sqrt(q) a(n).  (Real Dirichlet characters and the
Davenport-Heilbronn vector are elements; a dimension >= 2 means genuinely non-character members exist.)
Computed as the null space of the real linear system, with an SVD rank tolerance of 1e-9.
"""
import numpy as np

for q in range(3, 20):
    n = np.arange(q)
    F = np.exp(2j * np.pi * np.outer(n, n) / q)
    for a in (0, 1):
        rows = []
        # DFT self-duality: F x - i^a sqrt(q) x = 0 (real and imaginary parts)
        M = F - (1j ** a) * np.sqrt(q) * np.eye(q)
        rows += [M.real, M.imag]
        P = np.zeros((q, q))
        for k in range(q):
            P[k, k] += 1; P[k, (-k) % q] -= (-1) ** a       # parity
        rows.append(P)
        e0 = np.zeros((1, q)); e0[0, 0] = 1; rows.append(e0)  # a(0) = 0
        rows.append(np.ones((1, q)))                          # sum = 0
        A = np.vstack(rows)
        s = np.linalg.svd(A, compute_uv=False)
        dim = q - int(np.sum(s > 1e-9 * s[0]))
        print('q=%2d parity=%d  dim=%d' % (q, a, dim))
