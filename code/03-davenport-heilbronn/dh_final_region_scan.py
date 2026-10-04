"""Upper bound, step 2 (float feasibility, numeric-experiment): the plain final-region lower bound
  LB = 1 - sum_{2<=n<=N} |a_n| b_n n^{-sigma} - sum_{n<=N} |a_n| |gamma| n^y b_n n^{-sigma}
for the Davenport-Heilbronn coefficients (q = 5, parity 1), using the conventions of
code/shared-dirichlet-pipeline/newman_chi.py (Theorem A holds for DH: rs_identity_dh.py).
No Euler mollifier exists (f has no Euler product).  For each (t, y), report LB at x = 2T.
"""
import sys, math
import os as _os
sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'shared-dirichlet-pipeline'))
import numpy as np
import mpmath as mp
from newman_chi import Chi

k = float((mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1))
dh = Chi(5, 1, [0.0, 1.0, k, -k, -1.0])


def plain_lb(x, t, y):
    sig, gam = dh.sig_gam(x, y, t)
    N = int(math.floor(math.sqrt(5 * x / (4 * math.pi) + 5 * t / 16)))
    n = np.arange(1, N + 1, dtype=float)
    a = np.abs(np.array([dh(int(i)) for i in n]))
    L = np.log(n)
    w = a * np.exp(t / 4 * L * L - sig * L)
    return 1 - (w[1:].sum()) - gam * (w * np.exp(y * L)).sum(), N, sig, gam


if __name__ == '__main__':
    for T in (1e3, 1e4, 1e5, 1e6):
        x = 2 * T
        row = []
        for t in (0.3, 0.4, 0.5, 0.6, 0.8):
            for y in (0.2, 0.4):
                lb, N, sig, gam = plain_lb(x, t, y)
                row.append('t=%.1f y=%.1f LB=%+.3f' % (t, y, lb))
        print('T=%.0e N=%d |' % (T, N), ' ; '.join(row), flush=True)
