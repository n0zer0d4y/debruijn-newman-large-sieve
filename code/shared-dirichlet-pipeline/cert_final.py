"""Certificate, part (ii): the final region in verified arithmetic (papers/supporting/riemann-siegel-identity/derivation.md §3).

Same inequality and block structure as family_driver2.run (Polymath Lemma 8.5 with the
chi-mollifier, N in [N_s, 1e5]; plain triangle blocks on [1e5, 1e6]; uniform tail for
N >= 1e6), with:
  * sigma and |gamma| at cell corners in Arb (python-flint acb), rounded outward;
  * log m from Arb; exp by rig.vexp (no libm); every float sum carries a gamma_n bound;
  * the remainder Z of Polymath (97) bounded in every block:
      Z <= (sum_d |lam_d| d^-s0) * |kappa|max * L_b * e^{|kappa| L_b} * sum_m g_max,m E_m,
    |kappa| <= (t/2) y (1/x + 2/x^2)   (|alpha'(s)| <= 1/x + 2/x^2 for |s + a| >= x/2);
  * the N >= N0 tail proved for all N at once (see tail_bound).
Run: <Python with python-flint> cert_final.py D t0   -> JSON on stdout.
"""
import sys, json, math, time, itertools
import numpy as np
from flint import arb, acb, ctx
import rig
from family_driver2 import kron

ctx.prec = 96
U = rig.U
PI = arb.pi()


def down(a):
    return float(np.nextafter(float(a.lower()), -np.inf))


def up(a):
    return float(np.nextafter(float(a.upper()), np.inf))


class Char:
    def __init__(self, D, t0):
        self.D, self.t = D, t0
        q = self.q = abs(D); self.A = 0 if D > 0 else 1
        self.tab = [kron(D, n) for n in range(q)]
        self.XPL = 2 * max(1e8 / q, (7.5 if q % 2 == 0 else 3.75) * 1e7 / q + 200)
        self.lq = (arb(q) / PI).log()

    def chi(self, n):
        return self.tab[n % self.q]

    # --- Arb scalars -------------------------------------------------------
    def am(self, s):
        sa = s + self.A
        return self.lq / 2 + (sa / 2).log() / 2 - 1 / (2 * sa)

    def lM0(self, s):
        w = (s + self.A) / 2
        return w * self.lq + (2 * PI).sqrt().log() + (w - arb(1) / 2) * w.log() - w

    def lMt(self, s, t):
        return (t / 4) * self.am(s) ** 2 + self.lM0(s)

    def sig_gam(self, x, y, t):
        """(sigma, |gamma|) as arb balls at the exact point (x, y, t) (x, y, t arb)."""
        sp = (1 + y - acb(0, 1) * x) / 2
        sm = (1 - y + acb(0, 1) * x) / 2
        sig = (sp + (t / 2) * self.am(sp)).real
        g = (self.lMt(sm, t) - self.lMt(sp, t)).exp().abs_upper()
        gl = (self.lMt(sm, t) - self.lMt(sp, t)).exp().abs_lower()
        return sig, gl, g

    def xr(self, a, b, t):
        q = self.q
        return 4 * PI * (arb(a) ** 2 - arb(q) * t / 16) / q, 4 * PI * (arb(b + 1) ** 2 - arb(q) * t / 16) / q


def mollifier(ch, P, t):
    """lam_d = prod_{p|d} (-chi(p) b_p) in Arb; returns list of (d, lam_d as arb)."""
    out = []
    for r in range(len(P) + 1):
        for S in itertools.combinations(P, r):
            d, l = 1, arb(1)
            for p in S:
                d *= p
                l *= -ch.chi(p) * ((t / 4) * arb(p).log() ** 2).exp()
            out.append((d, l))
    return out


def weighted_E(L, t, s0, idx):
    """E_m = b_m m^{-s0} = exp(t/4 L^2 - s0 L) for the m in idx; returns (E, rel)."""
    Lm = L[idx]
    z = (t / 4) * Lm * Lm - s0 * Lm
    ez = 6 * U * ((t / 4) * Lm * Lm + s0 * Lm)
    v, rel = rig.vexp(z, ez)
    return v, float(np.max(rel))


def mpow(L, y, idx):
    """m^y = exp(y L) with its relative error bound."""
    Lm = L[idx]
    return rig.vexp(y * Lm, 3 * U * y * Lm)


class Final:
    def __init__(self, ch, L, MS, CHS):
        self.ch, self.L, self.MS, self.CHS = ch, L, MS, CHS

    def cnt(self, N):
        return int(np.searchsorted(self.MS, N, side='right'))

    def corners(self, a, b, y1, y2, t):
        ch = self.ch
        xlo, xhi = ch.xr(a, b, t)
        s0 = down(ch.sig_gam(xlo, arb(y1), arb(t))[0])
        cg = []
        for xx in (xlo, xhi):
            for yy in (y1, y2):
                _, gl, gu = ch.sig_gam(xx, arb(yy), arb(t))
                cg.append((down(gl), up(gu), yy))
        return xlo, s0, cg

    def gm_range(self, cg, idx, dy=0.0, xlo=None):
        """Per-m [g_min, g_max] over the four corners, rounded outward.  log g_m is monotone in x
        (Lemma M(b)) and within (dy^2/8) sup|d^2_y log g_m| <= dy^2/(8x) of affine in y
        (|d_y Re alpha(s+-)| <= 1/(2x) each), so the corner
        extremes are inflated by e^{dy^2/(4 x_lo)}; this covers m = N, where the sign of d_y may flip."""
        lo = None; hi = None
        for gl, gu, yy in cg:
            v, rel = mpow(self.L, yy, idx)
            vl = gl * v * (1 - rel - 4 * U); vh = gu * v * (1 + rel + 4 * U)
            lo = vl if lo is None else np.minimum(lo, vl)
            hi = vh if hi is None else np.maximum(hi, vh)
        if xlo is not None:
            f = math.exp(dy * dy / (4 * xlo)) * (1 + 4 * U)
            lo = lo / f; hi = hi * f
        return lo, hi

    def check_range(self, xl, b):
        """Hypotheses of Lemma M: x >= 1e3 and w_R > q/(2 pi) (w_R >= N >= b's lower end)."""
        assert xl >= 1e3, xl
        assert math.sqrt(self.ch.q * xl / (4 * math.pi)) > self.ch.q / (2 * math.pi)

    def kappa(self, xlo, y2, t):
        x = down(xlo)
        return t / 2 * y2 * (1 / x + 2 / x ** 2) * (1 + 8 * U)

    def blk(self, a, b, y1, y2, t, P):
        ch = self.ch
        xlo, s0, cg = self.corners(a, b, y1, y2, t)
        xl = down(xlo); self.check_range(xl, b)
        idx = self.MS[:self.cnt(b)]
        E, relE = weighted_E(self.L, t, s0, idx)
        gmin, gmax = self.gm_range(cg, idx, y2 - y1, xl)
        g0 = np.sqrt(gmin * gmax)
        w = np.maximum(gmax - g0, g0 - gmin) * (1 + 4 * U)
        chm = self.CHS[:len(idx)]
        # lam_d d^{-s0} with sign, as a double with relative error <= 2u (Arb radius checked)
        lamw = []
        for d, l in mollifier(ch, P, arb(t)):
            v = l * arb(d) ** (-arb(s0))
            assert v.rad() <= abs(v.mid()) * 2.0**-60
            lamw.append((d, float(v.mid())))
        inb = idx <= a
        m_in = idx[inb]; E_in = E[inb]; ch_in = chm[inb]; g0_in = g0[inb]; w_in = w[inb]
        ids, be, al, dv = [], [], [], []
        for d, lw in lamw:
            ids.append(d * m_in); t_ = lw * ch_in * E_in
            be.append(t_); al.append(t_ * g0_in); dv.append(abs(lw) * E_in * w_in)
        ids = np.concatenate(ids); uu, inv = np.unique(ids, return_inverse=True)
        B = np.bincount(inv, np.concatenate(be)); A = np.bincount(inv, np.concatenate(al))
        Dv = np.bincount(inv, np.concatenate(dv))
        absT = np.bincount(inv, np.abs(np.concatenate(be)) * (1 + np.concatenate([g0_in] * len(lamw))))
        k = (1 - gmin[0]) / (1 + gmin[0]) * (1 + 4 * U)
        mk = uu > 1
        terms = np.maximum(np.abs(B - A), k * np.abs(B + A)) + Dv
        Y = float(np.sum(terms[mk]))
        nterm = len(ids)
        rho = relE + 4 * U                      # per product term (lam_w: 2u, chi exact, E: relE, g0: 1u)
        errY = ((rho + rig.gamma_n(len(lamw) + 4)) * float(np.sum(absT[mk]) + np.sum(Dv[mk])) * 1.001
                + rig.gamma_n(nterm + 8) * (Y + float(np.sum(absT))) * 1.001)
        la = sum(abs(lw) for _, lw in lamw) * (1 + rig.gamma_n(len(lamw) + 2))
        tl = ~inb
        tr = la * float(np.sum(E[tl] * (1 + gmax[tl]))) * (1 + relE + rig.gamma_n(len(idx) + 4))
        kap = self.kappa(xlo, y2, t)
        Lb = self.L[idx[-1]]
        Z = la * kap * Lb * math.exp(kap * Lb) * float(np.sum(E * gmax)) * (1 + relE + rig.gamma_n(len(idx) + 8)) * 1.001
        lb = 1 - gmax[0] * (1 + 2 * U) - Y - errY - tr - Z
        return lb - 4 * U, tr, dict(errY=errY, Z=Z)

    def pln(self, a, b, y1, y2, t):
        xlo, s0, cg = self.corners(a, b, y1, y2, t)
        xl = down(xlo); self.check_range(xl, b)
        idx = self.MS[:self.cnt(b)]
        E, relE = weighted_E(self.L, t, s0, idx)
        _, gmax = self.gm_range(cg, idx, y2 - y1, xl)
        S = float(np.sum(E)); SG = float(np.sum(gmax * E))
        eps = (relE + rig.gamma_n(len(idx) + 4)) * 1.001
        kap = self.kappa(xlo, y2, t); Lb = self.L[idx[-1]]
        Z = kap * Lb * math.exp(kap * Lb) * SG * (1 + eps) * 1.001
        return 1 - (S * (1 + eps) - 1) - SG * (1 + eps) - Z - 4 * U


def tail_bound(ch, N0, y1, y2, t):
    """For every N > N0 and y in [y1, y2]: the contribution of n in (N0, N] to the plain
    bound is <= TA + TB (returned, rigorous), where, with L0 = log N0,
    c_A = (1-y1)/2 + t eps/2, c_B = (1+y2)/2 + t eps/2, eps = q t/(16 N0^2) + 1/x(N0):
      sum_{N0<n<=N} b_n n^{-sigma_N}          <= e^{-t L0^2/4 + c_A L0} / (e beta_A),
      sum_{N0<n<=N} |gamma| n^y b_n n^{-sigma_N} <= G e^{-t L0^2/4 + c_B L0} / (e beta_B),
    beta = t L0/2 - c, G = |gamma(x(N0), y1)| (|gamma| decreases in x and in y).
    Proof: sigma_N >= (1+y)/2 + (t/2)(log N - eps) (Re alpha(s+) >= (1/2)log(qx/4pi) - 1/x
    and qx/4pi >= N^2 - qt/16); the summand is decreasing on n <= N, so the sum is
    <= int_{L0}^{L_N} exp(phi(L)) dL with phi convex, <= (L_N - L0) max(endpoints)
    <= C0 v e^{-beta v} <= C0/(e beta), v = L_N - L0.
    """
    t = arb(t); q = ch.q
    xN0 = 4 * PI * (arb(N0) ** 2 - arb(q) * t / 16) / q
    eps = arb(q) * t / (16 * arb(N0) ** 2) + 1 / xN0
    L0 = arb(N0).log()
    cA = (1 - arb(y1)) / 2 + t * eps / 2
    cB = (1 + arb(y2)) / 2 + t * eps / 2
    bA = t * L0 / 2 - cA; bB = t * L0 / 2 - cB
    assert bA > 0 and bB > 0
    # the B-summand is decreasing on n <= N only if y2 < 1 - t eps
    assert arb(y2) < 1 - t * eps
    _, _, G = ch.sig_gam(xN0, arb(y1), t)
    e = arb(1).exp()
    TA = (-t * L0 ** 2 / 4 + cA * L0).exp() / (e * bA)
    TB = G * (-t * L0 ** 2 / 4 + cB * L0).exp() / (e * bB)
    # Z for N >= N0: |kappa| L_N e^{|kappa|L_N} * (sum g E) with sum g E <= 2 (head <= 1 + tails)
    Zt = t / 2 * arb(y2) * (1 / xN0 + 2 / xN0 ** 2) * (L0 + 1) * 3
    return up(TA), up(TB), up(Zt)


_CACHE = {}


def run(D, T0, log=None, ratio_m=1.02, yh=(0.01, 0.05), ratio_p=1.05, yp=(0.02, 0.1), thr=0.01, XC=None, max_fails=None):
    """ratio_m / yh: initial mollified N-block ratio and y-cell widths (below / above y = 0.3);
    ratio_p / yp: the same for the plain region; thr: acceptance threshold for a block
    (must exceed sup|E|(e_A+e_B+e_C), which cert_extra.final_errors bounds separately).
    Blocks below thr are subdivided adaptively, so coarser initial blocks only change speed."""
    T1 = time.time()
    rig._init_constants()
    ch = Char(D, T0)
    if 'L' not in _CACHE:
        _CACHE['L'] = rig.arb_logs(1000002)
    L = _CACHE['L']
    MS = np.array([k for k in range(1, 1000002) if ch.chi(k) != 0], dtype=np.int64)
    CHS = np.array([ch.chi(int(k)) for k in MS], dtype=float)
    F = Final(ch, L, MS, CHS)
    # barrier location (same deterministic rule as family_driver2; only a choice of X)
    q = ch.q
    Pp = np.array([p for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] if ch.chi(p) != 0])
    cP = np.array([ch.chi(int(p)) for p in Pp])
    if XC is None:
        xs = np.arange(ch.XPL - min(2e5, ch.XPL / 4), ch.XPL - 1.0, 0.05); val = np.zeros_like(xs)
        for p, c in zip(Pp, cP):
            val += -np.log(np.abs(1 - c * p ** (-0.6) * np.exp(-1j * xs / 2 * np.log(p))))
        XC = float(xs[int(np.argmax(val))])
    pr = [p for p in [2, 3, 5, 7, 11, 13, 17] if ch.chi(p) != 0]; ps, pb = pr[:5], pr[:3]
    t = T0; ytop = math.sqrt(1 - 2 * t)
    Ns = int(math.floor(math.sqrt(q * (XC + 0.5) / (4 * math.pi) + q * t / 16)))
    Nb = []; a = Ns
    while a < 100001:
        b = min(max(int(a * ratio_m), a + 1), 100001); Nb.append((a, b)); a = b + 1
    yc = []; y = 0.2
    while y < ytop - 1e-12:
        h = yh[0] if y < 0.3 else yh[1]; yc.append((y, min(y + h, ytop))); y = min(y + h, ytop)
    Q = [(a, b, y1, y2) for (a, b) in Nb for (y1, y2) in yc]
    mn = 9.0; mn_at = None; fails = []; nev = 0; maxZ = 0; maxerr = 0
    while Q:
        a, b, y1, y2 = Q.pop()
        v, tr, info = F.blk(a, b, y1, y2, t, ps if a < 15000 else pb); nev += 1
        maxZ = max(maxZ, info['Z']); maxerr = max(maxerr, info['errY'])
        if v > thr:
            if v < mn: mn, mn_at = v, (a, b, y1, y2)
            continue
        if nev > 20000 or (b - a <= 1 and y2 - y1 < 0.0005):
            fails.append((a, b, y1, y2, v))
            if max_fails is not None and len(fails) > max_fails:   # early exit: this t0 cannot certify
                return dict(D=D, t0=T0, XC=XC, Ns=Ns, blocks=nev, FR_min=mn, FR_min_at=mn_at, FR_fails=fails,
                            plain_min=None, tail_min=None, maxZ=maxZ, max_errY=maxerr, early_exit=True,
                            secs=round(time.time() - T1, 1))
            continue
        if b - a > 1 and tr > 0.3 * max(thr - v, 1e-3):
            mid = (a + b) // 2; Q += [(a, mid, y1, y2), (mid + 1, b, y1, y2)]
        else:
            ym = (y1 + y2) / 2; Q += [(a, b, y1, ym), (a, b, ym, y2)]
    pm = 9.0; a = 100001; last = {}
    while a < 1000001:
        b = min(int(a * ratio_p), 1000001); y = 0.2
        while y < ytop - 1e-12:
            y2 = min(y + (yp[0] if y < 0.3 else yp[1]), ytop)
            v = F.pln(a, b, y, y2, t); pm = min(pm, v)
            if b == 1000001: last[(y, y2)] = v
            y = y2
        a = b + 1
    tmin = 9.0
    for (y1, y2), v in last.items():
        TA, TB, Zt = tail_bound(ch, 1000001, y1, y2, t)
        tmin = min(tmin, v - TA - TB - Zt)
    out = dict(D=D, t0=T0, XC=XC, Ns=Ns, blocks=nev, FR_min=mn, FR_min_at=mn_at, FR_fails=fails,
               plain_min=pm, tail_min=tmin, maxZ=maxZ, max_errY=maxerr, secs=round(time.time() - T1, 1))
    return out


if __name__ == '__main__':
    D = int(sys.argv[1]); T0 = float(sys.argv[2])
    print(json.dumps(run(D, T0)))
