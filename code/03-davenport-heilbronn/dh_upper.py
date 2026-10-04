"""Upper bound, part (a): the final region for the Davenport-Heilbronn heat flow at t0, x >= X.

Goal: H_t0(x + iy) != 0 for x >= X and y0 <= y <= Ytop, where Ytop^2 = Delta^2 - 2 t0
(de Bruijn: every zero of H_t lies in Im^2 <= Delta^2 - 2t; Delta = 2 sigma0 - 1, claim #2).
Method: the final region (papers/supporting/riemann-siegel-identity/derivation.md §3, Lemmas E and F), with the DH coefficients.
  * Theorem A holds for DH (rs_identity_dh.py), so the RS approximation f_t = A + B, with the
    error terms e_A, e_B, e_C (Lemma I; K_h from Lemma H', c0 = 0 for q = 5), carries over.
  * There is no Euler product, hence no mollifier: |E| = 1 and plain triangle blocks are used
    from N(X) to 1e6 (cert_final.Final.pln takes weight 1 >= |a_n|, which is conservative), then
    Lemma F's uniform tail.
  * Error terms bounded uniformly for x >= X by Lemma E (cert_extra), with the y-range [y0, Ytop].
Verified arithmetic as in cert_final (Arb scalars, rig.vexp sums).  Single process.
"""
import sys, math, json, time
import os as _os
sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'shared-dirichlet-pipeline'))
import numpy as np
from flint import arb, acb, ctx
import rig, cert_final as cf, cert_barrier as cb, cert_extra as ce

ctx.prec = 96
PI = arb.pi()
K = ((10 - 2 * arb(5).sqrt()).sqrt() - 2) / (arb(5).sqrt() - 1)
KF = float(K.mid())


def tail_any_y(ch, N0, y1, y2, t):
    """For every N > N0 and y in [y1, y2] (any y2, no monotonicity needed): the n in (N0, N] part of
    the plain bound is <= e^{-t L0^2/4 + cA L0} + G e^{-t L0^2/4 + (cA + y2) L0},
    cA = (1 - y1)/2 + t eps/2, eps = q t/(16 N0^2) + 1/x(N0), G = |gamma(x(N0), y1)|.
    Proof: each part is <= N max(f(N0), f(N)) (the exponent is convex in log n); with
    sigma_N >= (1+y)/2 + (t/2)(log N - eps), both endpoint forms decrease in L_N = log N for
    L_N >= L0 once (t/2) L0 > max(1, cA + y2) (asserted), so the maximum is at L_N = L0."""
    t = arb(t); q = ch.q
    xN0 = 4 * PI * (arb(N0) ** 2 - arb(q) * t / 16) / q
    eps = arb(q) * t / (16 * arb(N0) ** 2) + 1 / xN0
    L0 = arb(N0).log()
    cA = (1 - arb(y1)) / 2 + t * eps / 2
    assert t * L0 / 2 > 1 and t * L0 / 2 > cA + arb(y2)
    _, _, G = ch.sig_gam(xN0, arb(y1), t)
    TA = (-t * L0 ** 2 / 4 + cA * L0).exp()
    TB = G * (-t * L0 ** 2 / 4 + (cA + arb(y2)) * L0).exp()
    Zt = t / 2 * arb(y2) * (1 / xN0 + 2 / xN0 ** 2) * (L0 + 1) * 3
    return cf.up(TA), cf.up(TB), cf.up(Zt)


def pln_w(F, a, b, y1, y2, t):
    """cert_final.Final.pln with the true weights |a_n| (1 or k) instead of 1:
    1 - (sum_{n<=b} |a_n| E_n - 1) - sum |a_n| g_max,n E_n - Z, all rigorous as in Final.pln
    (|a_n| rounded up by 2u)."""
    U = rig.U
    xlo, s0, cg = F.corners(a, b, y1, y2, t)
    xl = cf.down(xlo); F.check_range(xl, b)
    idx = F.MS[:F.cnt(b)]
    E, relE = cf.weighted_E(F.L, t, s0, idx)
    _, gmax = F.gm_range(cg, idx, y2 - y1, xl)
    wa = np.abs(F.CHS[:len(idx)]) * (1 + 2 * U)
    S = float(np.sum(wa * E)); SG = float(np.sum(wa * gmax * E))
    eps = (relE + rig.gamma_n(len(idx) + 6)) * 1.001
    kap = F.kappa(xlo, y2, t); Lb = F.L[idx[-1]]
    Z = kap * Lb * math.exp(kap * Lb) * SG * (1 + eps) * 1.001
    return 1 - (S * (1 + eps) - 1) - SG * (1 + eps) - Z - 4 * U


class DHChar(cf.Char):
    def __init__(self, t0):
        self.D, self.t = 'DH', t0
        self.q = 5; self.A = 1
        self.tab = [0.0, 1.0, KF, -KF, -1.0]
        self.XPL = None
        self.lq = (arb(5) / PI).log()

    def chi(self, n):
        return self.tab[n % 5]


def final_region(t0=0.5, y0=0.4, X=5.0e4, Delta=1.7902723164702195):
    T1 = time.time()
    rig._init_constants()
    ch = DHChar(t0)
    Ytop = math.sqrt(Delta ** 2 - 2 * t0) + 1e-9
    L = rig.arb_logs(1000002)
    MS = np.array([k_ for k_ in range(1, 1000002) if k_ % 5 != 0], dtype=np.int64)
    CHS = np.array([ch.chi(int(k_)) for k_ in MS], dtype=float)
    F = cf.Final(ch, L, MS, CHS)
    q, t = 5, t0
    Ns = int(math.floor(math.sqrt(q * X / (4 * math.pi) + q * t / 16)))
    ycells = []; y = y0
    while y < Ytop - 1e-12:
        h = 0.05 if y < 0.8 else 0.2
        ycells.append((y, min(y + h, Ytop))); y = min(y + h, Ytop)
    # plain blocks N in [Ns, 1e6]
    mn = 9.0; at = None; blocks = 0; a = Ns; fails = []
    stack = []
    while a < 1000001:
        b = min(max(int(a * 1.25), a + 1), 1000001)
        for (y1, y2) in ycells: stack.append((a, b, y1, y2))
        a = b + 1
    while stack:
        a, b, y1, y2 = stack.pop()
        v = pln_w(F, a, b, y1, y2, t); blocks += 1
        if v > 0.06:
            if v < mn: mn, at = v, (a, b, y1, y2)
            continue
        if b - a > 1:
            m = (a + b) // 2; stack += [(a, m, y1, y2), (m + 1, b, y1, y2)]
        elif y2 - y1 > 0.002:
            ym = (y1 + y2) / 2; stack += [(a, b, y1, ym), (a, b, ym, y2)]
        else:
            fails.append((a, b, y1, y2, v))
    # N >= 1e6 tail (Lemma F), for the last block values
    tmin = 9.0
    for (y1, y2) in ycells:
        v = pln_w(F, 1000001 - int(1000001 / 1.25), 1000001, y1, y2, t)
        TA, TB, Zt = tail_any_y(ch, 1000001, y1, y2, t)
        tmin = min(tmin, v - TA - TB - Zt)
    # error terms (Lemma E) with the DH y-range; |E| = 1
    xs = X; Tps = xs / 2                    # part (b)'s rectangle ends at x = X, so the bounds must hold from X
    cs = math.floor(math.sqrt(q * Tps / (2 * math.pi))) - 1
    Kh = cb.kh_analytic(ch)[0]
    lr = cf.up(ce.gamma_ratio_log_bound(ch, xs, y0, t).exp())
    eC = 0.0
    for sig in (arb((1 - Ytop) / 2).union(arb((1 - y0) / 2)), arb((1 + y0) / 2).union(arb((1 + Ytop) / 2))):
        eC += cf.up((PI ** 2 * arb(t) / 64).exp()) * lr * Kh * ce.gauss_uniform(ch, sig, t, cs, Tps)
    G = cf.up(ch.sig_gam(arb(xs), arb(y0), arb(t))[2])
    Fb = math.sqrt(4 * math.pi / t) * (math.exp(((1 - y0) / 2) ** 2 / t) + G * math.exp(((1 + Ytop) / 2) ** 2 / t)) + 2
    N = math.sqrt(q * xs / (4 * math.pi)); lN = math.log(N)
    eps = math.expm1(((t * t / 8) * (2 * lN + 3) ** 2 + t / 4 + 1 / 6) / (xs / 2 - 1.25)) * 1.01
    err = (eC + eps * Fb) * (1 + 1e-9)
    lbmin = min(mn, tmin)
    return dict(t0=t0, y0=y0, X=X, Ytop=Ytop, Ns=Ns, blocks=blocks, FR_min=mn, FR_min_at=at, fails=fails,
                tail_min=tmin, Kh=Kh, eC=eC, eAB=eps * Fb, err_total=err,
                ok=(not fails) and lbmin > err, Lambda_bound=t0 + y0 * y0 / 2, secs=round(time.time() - T1, 1))


if __name__ == '__main__':
    kw = {}
    for a_ in sys.argv[1:]:
        k_, v_ = a_.split('='); kw[k_] = float(v_)
    print(json.dumps(final_region(**kw), default=str))
