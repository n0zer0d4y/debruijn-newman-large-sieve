"""Certificate, part (iii): the barrier, with a rigorous K_h (Lemma H) and a rigorous
e_C (§2.7-2.8), in verified arithmetic.  Same covering as family_driver2.run, except:
  * each edge is cut into n = ceil(len/h) equal segments (the float driver's np.arange
    centres could leave up to one step uncovered at the far end of an edge);
  * N is asserted constant on the whole rectangle for all t in [0, t0];
  * sigma, kappa, gamma and the Euler factor E at every mesh point are Arb balls;
    the Dirichlet sums use rig.vexp / rig.vcis with explicit error bounds;
  * the Lipschitz constants are the float driver's termwise bounds (derivation §3.4),
    evaluated with outward rounding (the analytic slack terms are the same);
  * e_A + e_B <= eps_max * (F_A + F_B)  (Prop 6.1^chi), e_C by eC_box (Arb, per cell).
Run: <Python with python-flint> cert_barrier.py D t0 [XC]
"""
import sys, json, math, time
import numpy as np
from flint import arb, acb, ctx
import rig
from cert_final import Char, down, up, PI

ctx.prec = 96
U = rig.U
I_ = acb(0, 1)


# ---------------------------------------------------------------- K_h, analytic (Lemma H')
def kh_analytic(ch, du=0.01):
    """Lemma H' (any q).  On the line c + u e^{i pi/4} (c in Z + 1/2), with b = pi u/(sqrt2 q),
    beta = |1 - e^{-b}| and theta* = dist(1/2 + u/sqrt2, Z):
      (1/2q) sum_{n mod 2q} F_n^{-1/2} <= min( e^{b/2}(H_{q-1} + 1)/2 + min(1/(q beta), e^{b/2}/(4 theta*)),
                                             1/beta )
    Proof: drop chi(n) != 0 (all 2q residues); a_n/2 mod pi runs over pi(m + theta)/(2q), m < 2q, with
    theta = frac(c + u/sqrt2); sin x >= (2/pi) min(x, pi - x) gives sin(a_n/2) >= d_m/q with
    {d_m} = {j + theta, j + 1 - theta : 0 <= j < q}; G(d) = (beta^2 + 4 e^{-b} d^2/q^2)^{-1/2} is
    decreasing, so sum <= G(theta) + G(1-theta) + q e^{b/2} H_{q-1}; of theta, 1-theta one is >= 1/2,
    giving the e^{b/2}/2 (folded into '+1') and the min(...) term; and F >= beta^2 gives 1/beta.
    The sup over u is taken over cells of width du with outward interval bounds (theta* >= 0 when a
    half-integer crossing lies in the cell; the bound is decreasing in beta and increasing in b)."""
    q = ch.q
    H = sum(1.0 / j for j in range(1, q)) * (1 + 1e-12)
    s2 = math.sqrt(2)
    best = 0.0
    # u beyond which 1/beta <= (H+1)/2 for u > 0 (the first bound is >= (H+1)/2 e^{b/2} there)
    umax = 1.0
    while 1 / (1 - math.exp(-math.pi * umax / (s2 * q))) > (H + 1) / 2:
        umax *= 1.5
    # shared endpoints with the last one exactly umax (np.arange could stop an ulp short)
    edges = np.linspace(-umax, umax, int(math.ceil(2 * umax / du)) + 1)
    assert edges[0] == -umax and edges[-1] == umax
    lo, hi = edges[:-1], edges[1:]
    b_lo, b_hi = math.pi * lo / (s2 * q), math.pi * hi / (s2 * q)
    # beta minimum over the cell (0 if the cell contains u = 0); -expm1 avoids cancellation at small b
    beta = np.minimum(np.abs(np.expm1(-b_lo)), np.abs(np.expm1(-b_hi))) * (1 - 1e-9)
    beta = np.where((lo <= 0) & (hi >= 0), 0.0, beta)
    eb2 = np.exp(np.maximum(b_lo, b_hi) / 2) * (1 + 1e-9)
    # theta* minimum over the cell: 0 if 1/2 + u/sqrt2 crosses an integer in the cell; the rounding
    # of w = 1/2 + u/sqrt2 is <= 4u|w| <= 1e-12 for |u| <= 1e4, deflated below
    w_lo = 0.5 + lo / s2; w_hi = 0.5 + hi / s2
    cross = np.floor(w_hi) != np.floor(w_lo)
    fl = np.minimum(np.abs(w_lo - np.round(w_lo)), np.abs(w_hi - np.round(w_hi)))
    th = np.where(cross, 0.0, fl) - 1e-11
    with np.errstate(divide='ignore'):
        t1 = np.minimum(1 / (q * beta), eb2 / (4 * np.maximum(th, 0)))
        b1 = eb2 * (H + 1) / 2 + t1
        b2 = 1 / beta
    tot = np.minimum(b1, b2)
    best = float(np.max(tot))
    # |u| > umax: for u > 0 use 1/beta (decreasing); for u < 0, F >= (e^{|b|} - 1)^2
    tail = max(1 / (1 - math.exp(-math.pi * umax / (s2 * q))), 1 / (math.exp(math.pi * umax / (s2 * q)) - 1))
    c0 = sum(ch.chi(n) * (acb(0, -1) * (PI * n * n / q)).exp() for n in range(2 * q)) / (2 * q)
    return max(best, tail) * (1 + 1e-9) + up(c0.abs_upper()) / 2, up(c0.abs_upper())


# ---------------------------------------------------------------- K_h (Lemma H)
def kh_rig(ch, U_MAX=6.0, du=1e-3):
    """sup over the 2q half-integer line classes c and |u| <= U_MAX of
    (1/2q) sum_{chi(n) != 0} F_n(u)^{-1/2} + |c0|/2, with
    F = (1 - e^{-b})^2 + 4 e^{-b} sin^2(a/2), a = pi(c - n + u/sqrt2)/q, b = pi u/(sqrt2 q);
    each grid cell [u_i, u_i + du] is bounded below by interval reasoning:
    (1-e^{-b})^2 >= its minimum over the cell (0 if b changes sign), e^{-b} >= e^{-b_hi},
    sin^2(a/2) >= 0 if a zero of sin(a/2) lies in the cell, else the smaller endpoint value.
    For |u| > U_MAX the tail bound #{chi != 0}/(2q)/(1 - e^{-pi U_MAX/(sqrt2 q)}) applies."""
    q = ch.q
    ns = [n for n in range(2 * q) if ch.chi(n) != 0]
    s2 = math.sqrt(2)
    # shared endpoints, first = -U_MAX and last = U_MAX exactly (np.arange with
    # separately rounded ends could leave ulp-wide gaps and stop short of U_MAX)
    edges = np.linspace(-U_MAX, U_MAX, int(round(2 * U_MAX / du)) + 1)
    assert edges[0] == -U_MAX and edges[-1] == U_MAX
    ulo, uhi = edges[:-1], edges[1:]
    best = 0.0
    for kk in range(2 * q):
        cc = kk + 0.5
        tot = np.zeros_like(ulo)
        for n in ns:
            blo = math.pi * ulo / (s2 * q); bhi = math.pi * uhi / (s2 * q)
            # rounding of b = pi u/(sqrt2 q) propagated into vexp (<= 8u|b|)
            e_lo, re_lo = rig.vexp(-blo, 8 * U * np.abs(blo)); e_hi, re_hi = rig.vexp(-bhi, 8 * U * np.abs(bhi))
            ebmin = e_hi * (1 - re_hi - 4 * U)                            # e^{-b} >= e^{-b_hi}
            # |1 - e^{-b}| increases in |b|, so its minimum over a cell not containing b = 0
            # is at an endpoint; the vexp error is relative to e^{-b}.
            one_lo = np.abs(1 - e_lo); one_hi = np.abs(1 - e_hi)
            p1 = (np.minimum(one_lo, one_hi) - 1.01 * (re_lo * e_lo + re_hi * e_hi) - 4 * U)
            p1 = np.where((blo < 0) & (bhi > 0), 0.0, np.maximum(p1, 0.0))
            alo = math.pi * (cc - n + ulo / s2) / q; ahi = math.pi * (cc - n + uhi / s2) / q
            # zero of sin(a/2) in the cell iff a multiple of 2pi lies in [alo, ahi]
            kz = np.floor(ahi / (2 * math.pi))
            has0 = (kz * 2 * math.pi >= alo - 1e-12)
            h_lo = np.abs(alo / 2); h_hi = np.abs(ahi / 2)
            # rounding of a = pi(c - n + u/sqrt2)/q propagated into vcis (<= 8u|a|)
            c1, s1, e1 = rig.vcis(np.abs(alo / 2), 8 * U * np.abs(alo / 2) + 8 * U)
            c2, s2_, e2 = rig.vcis(np.abs(ahi / 2), 8 * U * np.abs(ahi / 2) + 8 * U)
            smin = np.minimum(np.abs(s1) - e1, np.abs(s2_) - e2)
            smin = np.where(has0, 0.0, np.maximum(smin, 0.0))
            F = p1 * p1 + 4 * ebmin * smin * smin
            F = F * (1 - 1e-10)
            tot += 1 / np.sqrt(F) * (1 + 4 * U)
        best = max(best, float(tot.max()) * (1 + rig.gamma_n(len(ns) + 2)) / (2 * q))
    tail = len(ns) / (2 * q) / (1 - math.exp(-math.pi * U_MAX / (s2 * q))) * (1 + 1e-12)
    c0 = sum(ch.chi(n) * acb(0, -1).__mul__(PI * n * n / q).exp() for n in range(2 * q)) / (2 * q)
    return max(best, tail) + up(c0.abs_upper()) / 2, up(c0.abs_upper())


def rsqrt(t):
    """sqrt of a nonnegative quantity enclosed by the ball t (which may touch 0)."""
    return arb(0).union(t.upper().sqrt()) if t.lower() <= 0 else t.sqrt()


# ---------------------------------------------------------------- e_C (§2.7 Step 5, §2.8)
def lemI(sig, Tp, c, q, rho=arb(1) / 4):
    Qs = arb(4) / q; s2 = arb(2).sqrt()
    lam = q * Tp / (2 * PI * c * c)
    l0 = Qs * PI * c * (lam - 1) / (2 * s2); l = l0 - sig / (s2 * c)
    k = Qs * PI * (1 + lam) / 4
    G0 = (PI / k).sqrt() * (l * l / (4 * k)).exp()
    A = Qs * PI * lam / (6 * c * (1 - rho)) + abs(sig) / (3 * c ** 3 * (1 - rho))
    kp = k - A * rho * c
    assert kp > k / 2
    # + Lemma I's Tail term (<= 1e-10000 for c >= 1000; 1e-300 is a generous stand-in)
    assert c > 1000
    return (G0 + A * (l * l / (4 * kp)).exp() * 4 * (1 / kp ** 2 + abs(l / (2 * kp)) ** 3 * (PI / kp).sqrt())
            + arb(10) ** -300)


def f_int(sv, Tp, c, de, q):
    return ((1 + de / c) ** sv * lemI(sv, Tp, c, q)
            * (1 / (6 * (Tp - arb('0.66'))) + abs(sv) / Tp ** 2 + sv * sv / (4 * (Tp - 1))).exp())


def tail_v(sig, t, Tp, c, de, q, V=12, lam=None, l0=None):
    """2/sqrt(pi) int_V^inf of the majorant exp(K0 + K1|s_v| + K2 s_v^2) e^{-v^2};
    derivation: log f <= K0 + K1|s| + K2 s^2 with
      K0 = 1 + l0^2/k + log(sqrt(pi/k) + 4 max(A0, A1) (4/k^2 + beta')),  K1 = 5 + 2/c,
      K2 = 1/Tp + 1/(2 c^2 k), valid while kp >= k/2 (checked at |s| = Tp/4);
    |s_v| <= |sig| + sqrt(t) v and s_v^2 <= 2 sig^2 + 2 t v^2; |s_v| > Tp/4 is Lemma T."""
    Qs = arb(4) / q; s2 = arb(2).sqrt(); rho = arb(1) / 4
    # lam, l0 may be passed as balls (uniform over delta, eps); c and Tp are then the
    # smallest values, where A0, A1, K1, K2 are largest
    if lam is None:
        lam = q * Tp / (2 * PI * c * c)
    if l0 is None:
        l0 = Qs * PI * c * (lam - 1) / (2 * s2)
    k = Qs * PI * (1 + lam) / 4
    A0 = Qs * PI * lam / (6 * c * (1 - rho)); A1 = 1 / (3 * c ** 3 * (1 - rho))
    Amax = A0 + A1 * Tp / 4
    assert k - Amax * rho * c > k / 2
    # sums instead of max(): Python's max on overlapping balls may return the smaller one
    bet = ((abs(l0) + 1 / (s2 * c)) / k) ** 3 * (2 * PI / k).sqrt()
    K0 = 1 + l0 * l0 / k + ((PI / k).sqrt() + 4 * (A0 + A1) * (4 / k ** 2 + bet)).log()
    K1 = 5 + 2 / c; K2 = 1 / Tp + 1 / (2 * c * c * k)
    a = 1 - 2 * K2 * t; b = K1 * rsqrt(t)
    assert 2 * a * V > b
    pre = (K0 + K1 * abs(sig) + 2 * K2 * sig * sig).exp()
    return 2 / PI.sqrt() * pre * (b * V - a * V * V).exp() / (2 * a * V - b)


def gauss_avg(sig, t, Tp, c, de, q, V=12, nseg=48):
    edges = [arb(-V) + arb(2 * V) * j / nseg for j in range(nseg + 1)]
    tot = arb(0)
    for j in range(nseg):
        v = arb(edges[j]).union(edges[j + 1])
        w = (edges[j + 1].erf() - edges[j].erf()) / 2
        tot += f_int(sig + rsqrt(t) * v, Tp, c, de, q) * w
    return tot + tail_v(sig, t, Tp, c, de, q, V) + arb(10) ** -300   # + Lemma T region


def eC_box(ch, x1, x2, y1, y2, t1, t2, Kh, nx=4, ny=4):
    """Upper bound for e_C over [x1,x2] x [y1,y2] x [t1,t2] (sub-boxes as Arb balls)."""
    q = ch.q; best = 0.0
    for i in range(nx):
        xb = arb(x1 + (x2 - x1) * i / nx).union(arb(x1 + (x2 - x1) * (i + 1) / nx))
        for j in range(ny):
            yb = arb(y1 + (y2 - y1) * j / ny).union(arb(y1 + (y2 - y1) * (j + 1) / ny))
            tb = arb(t1).union(arb(t2))
            Tp = xb / 2 + PI * tb / 8
            wR = (q * Tp / (2 * PI)).sqrt()
            fl, fh = math.floor(float(wR.lower())), math.floor(float(wR.upper()))
            assert fl == fh, 'c changes inside a sub-box'
            c = arb(fl) + arb(1) / 2; de = wR - c
            sp = (1 + yb - I_ * xb) / 2
            ratio = ((ch.lM0(I_ * Tp)).real - (ch.lMt(sp, tb)).real).exp()
            tot = arb(0)
            for sig in ((1 - yb) / 2, (1 + yb) / 2):
                tot += (PI ** 2 * tb / 64).exp() * ratio * gauss_avg(sig, tb, Tp, c, de, q) * Kh
            best = max(best, up(tot))
    return best


# ---------------------------------------------------------------- the barrier sums
class Barrier:
    def __init__(self, ch, XC, T0):
        self.ch, self.XC, self.T0 = ch, XC, T0
        q = ch.q
        self.N = int(math.floor(math.sqrt(q * XC / (4 * math.pi))))
        n_lo = math.floor(math.sqrt(q * (XC - 0.5) / (4 * math.pi)))
        n_hi = math.floor(math.sqrt(q * (XC + 0.5) / (4 * math.pi) + q * T0 / 16))
        assert n_lo == n_hi == self.N, (n_lo, n_hi, self.N)
        # Lemma M hypotheses on the barrier rectangle (asserted)
        assert XC - 0.5 >= 1e3 and math.sqrt(q * (XC - 0.5) / (4 * math.pi)) > q / (2 * math.pi)
        nn = [k for k in range(1, self.N + 1) if ch.chi(k) != 0]
        self.nn = np.array(nn, dtype=float)
        L = rig.arb_logs(self.N)
        self.LN = L[np.array(nn)]
        self.CN = np.array([ch.chi(k) for k in nn], dtype=float)
        P = [p for p in [2, 3, 5, 7, 11, 13, 17, 19, 23] if ch.chi(p) != 0][:7]
        self.P = P

    def point(self, x, y, t):
        """Arb data at one mesh point: (sr, si, kr, ki, gam, E) as balls."""
        ch = self.ch
        xb, yb, tb = arb(x), arb(y), arb(t)
        sp = (1 + yb - I_ * xb) / 2; sm = (1 - yb + I_ * xb) / 2
        ss = sp + (tb / 2) * ch.am(sp)
        kap = (tb / 2) * (ch.am(sm) - ch.am((1 + yb + I_ * xb) / 2))
        gam = (ch.lMt(sm, tb) - ch.lMt(sp, tb)).exp()
        E = acb(1)
        for p in self.P:
            E *= 1 - ch.chi(p) * ((tb / 4) * arb(p).log() ** 2).exp() * (-ss * arb(p).log()).exp()
        return ss, kap, gam, E

    def Ef(self, X, Y, t):
        """Enclosures of E f_t at the points (X_i, Y_i): returns (v, err) arrays."""
        LN, CN = self.LN, self.CN
        bl = (t / 4) * LN * LN
        out = np.empty(len(X), dtype=complex); err = np.empty(len(X))
        for i, (x, y) in enumerate(zip(X, Y)):
            ss, kap, gam, E = self.point(x, y, t)
            sr, si = float(ss.real.mid()), float(ss.imag.mid())
            esr = up(ss.real.rad()) + U * abs(sr); esi = up(ss.imag.rad()) + U * abs(si)
            kr, ki = float(kap.real.mid()), float(kap.imag.mid())
            ekr = up(kap.real.rad()) + U * abs(kr); eki = up(kap.imag.rad()) + U * abs(ki)
            # A: sum chi b_n n^{-ss} = sum chi e^{bl - sr L} e^{-i si L}
            zA = bl - sr * LN; ezA = 6 * U * (np.abs(bl) + np.abs(sr * LN)) + esr * LN
            mA, rA = rig.vexp(zA, ezA)
            thA = -si * LN; assert np.all(thA >= 0)          # si < 0
            cA, sA, eA_ = rig.vcis(thA, 3 * U * thA + esi * LN)
            # e^{-i si L} = e^{+i thA} = cos + i sin
            A = np.sum(CN * mA * (cA + 1j * sA))
            magA = float(np.sum(mA))
            errA = float(np.sum(mA * (rA + 1.5 * eA_ + 4 * U))) + rig.gamma_n(len(LN) + 4) * magA * 2
            # B: sum chi b_n n^{-(conj ss - y + kap)} = sum chi e^{bl - (sr - y + kr)L} e^{i(si - ki) L}
            zB = bl - (sr - y + kr) * LN
            ezB = 6 * U * (np.abs(bl) + np.abs((sr - y + kr) * LN)) + (esr + ekr + U * y) * LN
            mB, rB = rig.vexp(zB, ezB)
            thB = -(si - ki) * LN; assert np.all(thB >= 0)
            cB, sB, eB_ = rig.vcis(thB, 3 * U * thB + (esi + eki) * LN)
            B = np.sum(CN * mB * (cB - 1j * sB))                  # e^{i(si-ki)L} = e^{-i thB}
            magB = float(np.sum(mB))
            errB = float(np.sum(mB * (rB + 1.5 * eB_ + 4 * U))) + rig.gamma_n(len(LN) + 4) * magB * 2
            g = complex(float(gam.real.mid()), float(gam.imag.mid()))
            eg = up(gam.real.rad()) + up(gam.imag.rad()) + 2 * U * abs(g)
            e = complex(float(E.real.mid()), float(E.imag.mid()))
            ee = up(E.real.rad()) + up(E.imag.rad()) + 2 * U * abs(e)
            f = A + g * B
            ef = errA + (abs(g) + eg) * errB + eg * (magB + errB) + 8 * U * (magA + abs(g) * magB)
            v = e * f
            ev = (abs(e) + ee) * ef + ee * abs(f) + 8 * U * abs(e) * (magA + abs(g) * magB)
            out[i] = v; err[i] = ev * (1 + 1e-9)
        return out, err

    def gam_abs(self, x, y, t):
        ch = self.ch
        sp = (1 + arb(y) - I_ * arb(x)) / 2; sm = (1 - arb(y) + I_ * arb(x)) / 2
        return up((ch.lMt(sm, arb(t)) - ch.lMt(sp, arb(t))).exp().abs_upper())

    def lip(self, yl, tl, th):
        """Upper bounds (Ls, Lt, Em, F) for the cell: the termwise bounds of derivation §3.4,
        with sigma taken at (XC - 1/2, yl, tl) (sigma increases in x, y, t), b at th, and
        |gamma| n^y <= 1.01 max(|gamma(yl)| n^yl, |gamma(1)| n) (log|gamma| is linear in t,
        and log-affine in y and constant in x up to O(1/x) = 4e-8 relative)."""
        ch, XC = self.ch, self.XC
        LN, nn = self.LN, self.nn
        sp = (1 + arb(yl) - I_ * arb(XC - 0.5)) / 2
        sr = down((sp + (arb(tl) / 2) * ch.am(sp)).real)
        spc = (1 + arb(yl) - I_ * arb(XC)) / 2; smc = (1 - arb(yl) + I_ * arb(XC)) / 2
        ap = up(ch.am(spc).abs_upper()); amm = up(ch.am(smc).abs_upper())
        dif = up((ch.am(smc) ** 2 - ch.am(spc) ** 2).abs_upper())
        gl = max(self.gam_abs(XC, yl, tl), self.gam_abs(XC, yl, th)) * 1.01
        g1 = max(self.gam_abs(XC, 1.0, tl), self.gam_abs(XC, 1.0, th)) * 1.01
        bw, rb = rig.vexp((th / 4) * LN * LN - sr * LN, 6 * U * ((th / 4) * LN * LN + sr * LN))
        up_ = 1 + float(np.max(rb)) + 8 * U
        wA = bw * up_
        pl, rp = rig.vexp(yl * LN, 3 * U * yl * LN)
        wB = wA * np.maximum(gl * pl * (1 + float(np.max(rp)) + 4 * U), g1 * nn)
        Gxy = 0.5 * (ap + amm) + 1; Gt = dif / 4 + 1; aa = max(ap, amm)
        s = 1 + rig.gamma_n(len(LN) + 8)
        F = (wA.sum() + wB.sum()) * s
        Fs = ((wA * 0.51 * LN).sum() + (wB * (1.51 * LN + Gxy)).sum()) * s
        Ft = ((wA * (LN ** 2 / 4 + aa * LN / 2 + 0.1)).sum() + (wB * (LN ** 2 / 4 + aa * LN / 2 + 0.1 + Gt)).sum()) * s
        LP = np.array([float(arb(p).log().mid()) for p in self.P])
        bp, rbp = rig.vexp((th / 4) * LP ** 2 - sr * LP, 6 * U * ((th / 4) * LP ** 2 + sr * LP))
        bp = bp * (1 + float(np.max(rbp)) + 4 * U)
        assert np.all(bp < 1)
        Em = float(np.prod(1 + bp)) * (1 + 16 * U)
        Es = Em * float(np.sum(bp * 0.51 * LP / (1 - bp))) * (1 + 32 * U)
        Et = Em * float(np.sum(bp * (LP ** 2 / 4 + aa * LP / 2) / (1 - bp))) * (1 + 32 * U)
        return (Es * F + Em * Fs) * (1 + 8 * U), (Et * F + Em * Ft) * (1 + 8 * U), Em, F, aa

    def eps_AB(self, yl, tl, th, aa, F):
        """e_A + e_B <= eps_max * F (Prop 6.1^chi), eps = exp([(t^2/8)(|alpha| + log N)^2 + t/4 + 1/6]/(T - 1.25)) - 1."""
        T = (self.XC - 0.5) / 2
        lN = math.log(self.N) * (1 + 4 * U)
        z = ((th ** 2 / 8) * (aa + 1 + lN) ** 2 + th / 4 + 1 / 6) / (T - 1.25)
        return math.expm1(z) * 1.01 * F


def edge_points(a, b, h):
    n = int(math.ceil((b - a) / h)); hh = (b - a) / n
    return a + hh * (np.arange(n) + 0.5), hh


def run(D, T0, XC=None, tmax_secs=1e9, log=None):
    T1 = time.time()
    rig._init_constants()
    ch = Char(D, T0)
    if XC is None:
        from cert_final import run as _  # noqa: F401  (XC rule lives in cert_final.run)
        q = ch.q
        Pp = np.array([p for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] if ch.chi(p) != 0])
        cP = np.array([ch.chi(int(p)) for p in Pp])
        xs = np.arange(ch.XPL - 2e5, ch.XPL - 1.0, 0.05); val = np.zeros_like(xs)
        for p, c in zip(Pp, cP):
            val += -np.log(np.abs(1 - c * p ** (-0.6) * np.exp(-1j * xs / 2 * np.log(p))))
        XC = float(xs[int(np.argmax(val))])
    Kh, c0 = kh_rig(ch)
    B = Barrier(ch, XC, T0)
    yl_ = lambda t: math.sqrt(0.04 + 2 * (T0 - t))
    t = 0.0; cells = 0; worst = 9.0; mind = 9.0; Emax = 0.0; ok = True; npts = 0; eCmax = 0.0
    while t < T0 - 1e-12:
        tau = 0.02
        while True:
            th = min(t + tau, T0); yl = yl_(th); Ls, Lt, Em, F, aa = B.lip(yl, t, th)
            if Lt * (th - t) / 2 < 0.3: break
            tau /= 1.5
        eC = eC_box(ch, XC - 0.5, XC + 0.5, yl, 1.0, t, th, arb(Kh))
        eAB = B.eps_AB(yl, t, th, aa, F)
        eCmax = max(eCmax, eC)
        Emax = max(Emax, Em); err = Em * (eC + eAB)
        h = min(0.05, 0.6 / Ls)
        xs, hx = edge_points(XC - 0.5, XC + 0.5, h); ys, hy = edge_points(yl, 1.0, h)
        X = np.concatenate([xs, xs, np.full(len(ys), XC - 0.5), np.full(len(ys), XC + 0.5)])
        Y = np.concatenate([np.full(len(xs), yl), np.full(len(xs), 1.0), ys, ys])
        v, ev = B.Ef(X, Y, (t + th) / 2); npts += len(X)
        d = np.where(v.real >= 0, np.abs(v), np.abs(v.imag)) - ev - 4 * U * np.abs(v)
        mind = min(mind, float(d.min()))
        mg = float(d.min()) - err - Ls * max(hx, hy) / 2 - Lt * (th - t) / 2
        worst = min(worst, mg); cells += 1
        if mg <= 0: ok = False; break
        t = th
        if time.time() - T1 > tmax_secs: ok = False; break
    return dict(D=D, t0=T0, XC=XC, N=B.N, Kh_rig=Kh, c0=c0, eC_max=eCmax, barrier_ok=ok, barrier_to=t,
                cells=cells, points=npts, barrier_worst=worst, mind=mind, supE=Emax, secs=round(time.time() - T1, 1))


if __name__ == '__main__':
    D = int(sys.argv[1]); T0 = float(sys.argv[2]); XC = float(sys.argv[3]) if len(sys.argv) > 3 else None
    print(json.dumps(run(D, T0, XC)))
