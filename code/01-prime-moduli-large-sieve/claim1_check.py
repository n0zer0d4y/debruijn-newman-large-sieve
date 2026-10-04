"""(numeric) Check of the Theorem of papers/01-prime-moduli-large-sieve at Q = 4000, eta = 1/4 (main.md Remark 3).
LS = sum_{p<=Q} sum_{a=1}^{p-1} |S(a/p)|^2 with a_n = e(-n x0), |S(a/p)| = |sin(pi N th)/sin(pi th)|, th = a/p - x0.
The chosen class is the one holding the most primes of ((1-eta)Q, Q] (Step 1). Float64; a sanity check, not a proof."""
import numpy as np, math
from sympy import primerange, totient, primefactors
Q = 4000; eta = 0.25
P = np.array(list(primerange(2, Q + 1)))
A = np.concatenate([np.arange(1, p) for p in P]); Pp = np.repeat(P, P - 1)
X = A / Pp
npts = len(X)
def LS(N, x0):
    th = X - x0; th = th - np.round(th)
    s = np.sin(np.pi * th); num = np.sin(np.pi * N * th)
    v = np.where(np.abs(s) < 1e-15, float(N) ** 2, (num / np.where(s == 0, 1, s)) ** 2)
    return float(v.sum())
win = [int(p) for p in P if (1 - eta) * Q < p <= Q]
print("#points =", npts, " Pi_eta =", len(win))
for k in (1, 2, 6, 30, 210):
    phi = int(totient(k)); om = len(primefactors(k))
    if k == 1: b, cls = 0, win
    else:
        cnt = {}
        for p in win:
            if k % p: cnt.setdefault(int(p % k), []).append(int(p))
        c, cls = max(cnt.items(), key=lambda kv: len(kv[1]))
        b = (-pow(c, -1, k)) % k
    # Theorem: N = floor(k(1-eta)Q/eta), x0 = midpoint of the spoke interval
    N = int(k * (1 - eta) * Q / eta)
    x0 = b / k + 0.5 * (1 / (k * (1 - eta) * Q) + 1 / (k * Q))
    ls = LS(N, x0); rhs = (4 / math.pi ** 2) * (len(win) - om) / phi * N * N
    # each spoke point individually >= 2N/pi ?
    pts = np.array([((1 + b * p) // k) / p for p in cls]); th = pts - x0
    mn = np.min(np.abs(np.sin(np.pi * N * th) / np.sin(np.pi * th)))
    # Remark 3: N = kQ, x0 = b/k + 1/(kQ(1 - eta/2))
    N3 = k * Q; x3 = b / k + 1 / (k * Q * (1 - eta / 2))
    r3 = LS(N3, x3) / N3 / (npts + N3)
    print(f"k={k:4d} class size={len(cls):4d} N={N:8d}  LS/|a|^2={ls/N:.4e}  Thm RHS/|a|^2={rhs/N:.4e}  holds={ls >= rhs}  "
          f"min|S|/(2N/pi)={mn/(2*N/math.pi):.4f}  Remark3 ratio={r3:.3f} (k/phi={k/phi:.3f})")
