"""Upper bound, part (b): H_t0 has no zeros in R = [0, X] x [y0, Ytop] (t0 = 0.5, y0 = 0.4, X = 5e4,
Ytop^2 = Delta^2 - 2 t0), for the Davenport-Heilbronn heat flow.  (Conditional on the open items O1-O4 of papers/supporting/riemann-siegel-identity/derivation.md.)

Argument.  B_t(z) = exp(lMt((1 + y - ix)/2, t)) (cert_final.Char conventions, q = 5, A = 1) is analytic and
nonvanishing on a neighbourhood of R (Re((s + 1)/2) > 0 there).  G := H_t/B_t is analytic on R.  If
G avoids (-inf, 0] on the boundary of R, its winding number around R is 0, so G, hence H_t, has no zeros in R.
The boundary is covered by:
  (RS)  bottom and top edges, x >= 700 (bottom) and x >= 400 (top) for the 0.58 run:  G = g_N + O(e_C + e_AB),
        g_N = A(s*) + gamma B(w) with
        N = floor(sqrt(5x/4pi + 5t/16)) constant on the block [x_N, x_{N+1}] (then c = N + 1/2 in Lemma I);
        e_C from Lemma I via gauss_fixed, with the tail q-scaled (Qs = 4/q; Lemma I's tail needs Qs pi U/4 > a,
        which fails at c = 10.5, and near x = 640 the RS error exceeds the bottom-edge margin, hence the RS starts
        at 700 and 400), the ratio |M0(iT')|/|B_t| in Arb
        over the block, e_AB from Prop 6.1^chi (the formula of dh_upper.final_region).  Endpoint x_{N+1}
        of block N is covered by block N + 1 (where N is correct).
  (D)   bottom edge x in [0, 704], top edge x in [0, 404] (0.58 run), and the left side x = 0: G by the direct Arb heat integral
        (cert_dh.DH.xi, acb.integral, rigorous tail).
  (a)   the right side x = X: |G - 1| < 1 by the final region of part (a), so Re G > 0.
Each edge piece: Chebyshev interpolant p of the analytic function (g_N or G) at J+1 2nd-kind points,
|f - p| <= 4 M rho^-J/(rho - 1) + Lam_J max(node error) + float error (Trefethen ATAP Thm 8.2), M = sup|f| on
the Bernstein ellipse (Arb), and dist(p, (-inf, 0]) on a fine theta-grid minus a Markov-Lipschitz step.
The certificate for the piece is  dist(p) - e_grid - |f - p| - err > 0  (err = e_C + e_AB on RS pieces).
Usage: python dh_rect.py rs <edge> <xa> <xb>  |  python dh_rect.py direct <edge> <a> <b>
(edge in bottom, top, left; for 'left' a, b are y-limits).  Output: one JSON line.
"""
import sys, json, math, time
import os as _os
sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'shared-dirichlet-pipeline'))
import numpy as np
from flint import arb, acb, ctx
import rig, cert_final as cf, cert_barrier as cb, cert_barrier2 as cb2, cert_extra as ce
import cert_dh
from dh_upper import DHChar, K as KAPPA, KF

U = rig.U
I_ = acb(0, 1)
T0, Y0, XR = 0.5, 0.4, 5.0e4
DELTA = 1.7902723164702195                      # claim #2 (certificate): zeros of H_0 in |Im z| <= DELTA
YTOP = math.sqrt(DELTA ** 2 - 2 * T0) + 1e-9
PI = arb.pi()
KERR = float((KAPPA - KF).abs_upper()) + 1e-300   # |kappa - float(kappa)|


def cheb1(F):
    """values at 2nd-kind points u_j = cos(pi j/J) -> Chebyshev coefficients."""
    J = len(F) - 1
    j = np.arange(J + 1)
    M = np.cos(np.pi * np.outer(j, j) / J) * (2.0 / J)
    M[:, 0] *= 0.5; M[:, -1] *= 0.5; M[0, :] *= 0.5; M[-1, :] *= 0.5
    return M @ F


def cheb_lb(F, evmax, M, rho, nu=2000):
    """lower bound for dist(f, (-inf, 0]) on [-1, 1], f analytic in E_rho with |f| <= M there,
    F = float node values with |F_j - f(u_j)| <= evmax."""
    J = len(F) - 1
    C = cheb1(F)
    e_int = 4 * M * rho ** -J / (rho - 1)
    e_node = cb2.lebesgue(J) * evmax
    g = rig.gamma_n(J + 6)
    e_float = (g * (2.0 / J) * float(np.abs(F).sum()) * (J + 1) + g * float(np.abs(C).sum())) * 2
    th = np.linspace(0, math.pi, nu)
    P = np.cos(np.outer(th, np.arange(J + 1))) @ C
    d = np.where(P.real >= 0, np.abs(P), np.abs(P.imag))
    Lu = float(np.sum(np.abs(C) * np.arange(J + 1)))
    e_grid = Lu * (math.pi / (nu - 1)) / 2
    dmin = float(d.min())
    return dmin - e_grid - e_int - e_node - e_float, dict(dmin=dmin, e_grid=e_grid, e_int=e_int, e_node=e_node,
                                                            e_float=e_float, M=M)


def ellipse(c0, h, rho):
    a = h * (rho + 1 / rho) / 2; b = h * (rho - 1 / rho) / 2
    return arb(c0 - a).union(arb(c0 + a)), arb(-b).union(arb(b))


# ------------------------------------------------------------------ (RS) pieces
class RSN(cb.Barrier):
    """The Barrier sums with a prescribed N and no Euler factor (E = 1: f has no Euler product)."""
    def __init__(self, ch, N):
        self.ch, self.N = ch, N
        nn = [k for k in range(1, N + 1) if ch.chi(k) != 0]
        self.nn = np.array(nn, dtype=float)
        L = rig.arb_logs(N)
        self.LN = L[np.array(nn)]
        self.CN = np.array([ch.chi(k) for k in nn], dtype=float)
        self.WN = np.abs(self.CN) * (1 + 2 * U) + KERR
        self.P = []

    def majorant(self, xr, xi, y, t):
        """sup |g_N| over x in xr + i xi (z = x + iy), t fixed: sum |a_n| b_n n^{-Re s*} + |gamma| sum |a_n| b_n n^{-Re w}."""
        ch = self.ch
        x = acb(xr, xi); yb = arb(y); tb = arb(t)
        sp = (1 + yb - I_ * x) / 2; sm = (1 - yb + I_ * x) / 2
        ss = sp + (tb / 2) * ch.am(sp); w = sm + (tb / 2) * ch.am(sm)
        slo = cf.down(ss.real); wlo = cf.down(w.real)
        ghi = cf.up((ch.lMt(sm, tb) - ch.lMt(sp, tb)).real.exp())
        L = self.LN; th = t * (1 + 2 * U)
        a, ra = rig.vexp(th / 4 * L * L - slo * L, 6 * U * (th / 4 * L * L + abs(slo) * L))
        b, rb = rig.vexp(th / 4 * L * L - wlo * L, 6 * U * (th / 4 * L * L + abs(wlo) * L))
        s = 1 + rig.gamma_n(len(L) + 4)
        MA = float(np.sum(self.WN * a)) * (1 + float(np.max(ra))) * s
        MB = float(np.sum(self.WN * b)) * (1 + float(np.max(rb))) * s
        return (MA + ghi * MB) * (1 + 1e-12)


def x_of_N(N, t=T0, q=5):
    return 4 * math.pi * (N * N - q * t / 16) / q


def lemI_tail_fixed(sv, l0, c, q):
    """Lemma I's Tail for one fixed c (cert_extra.lemI_tail without the uniformity in c).
    On w = c + u e^{i pi/4}, |e^{i pi w^2/q}| = e^{-(pi/q)(sqrt2 c u + u^2)}, so Step 3 gives the Gaussian
    e^{-(pi/(2q)) u^2} = e^{-Qs pi u^2/8}, Qs = 4/q (the q = 4 derivation is the case Qs = 1):
    Tail <= 2 int_U^inf e^{0.35|sigma| + a u - Qs pi u^2/8} du <= 2 e^{0.35|sigma|} e^{aU - Qs pi U^2/8}/(Qs pi U/4 - a),
    U = c/4, a = |sigma|/c + |l0|, since the exponent is concave with slope a - Qs pi u/4 < 0 on u >= U."""
    Qs = arb(4) / q
    U_ = arb(c) / 4
    s = abs(sv).upper(); la = abs(l0).upper()
    a = s / arb(c) + la
    assert Qs * PI * U_ / 4 > a
    return 2 * (arb('0.35') * s).exp() * (a * U_ - Qs * PI * U_ * U_ / 8).exp() / (Qs * PI * U_ / 4 - a)


def gauss_fixed(ch, sig, t, c, Tps, dlo=-0.5, dhi=0.5):
    """cert_extra.gauss_uniform with c fixed (eps = 1/c a point), for every delta in [dlo, dhi] and T' >= Tps."""
    q = ch.q; Qs = arb(4) / q; s2 = arb(2).sqrt(); rho = arb(1) / 4
    eps = 1 / arb(c); de = arb(dlo).union(arb(dhi)); iT = arb(0).union(1 / arb(Tps))
    t = arb(t)

    def f(sv):
        lam = (1 + de * eps) ** 2
        l0 = Qs * PI * de * (2 + de * eps) / (2 * s2)
        l = l0 - sv * eps / s2
        k = Qs * PI * (1 + lam) / 4
        G0 = (PI / k).sqrt() * (l * l / (4 * k)).exp()
        A = Qs * PI * lam * eps / (6 * (1 - rho)) + abs(sv) * eps ** 3 / (3 * (1 - rho))
        kp = k - Qs * PI * lam * rho / (6 * (1 - rho)) - abs(sv) * rho * eps ** 2 / (3 * (1 - rho))
        assert kp > k / 2
        lem = G0 + A * (l * l / (4 * kp)).exp() * 4 * (1 / kp ** 2 + abs(l / (2 * kp)) ** 3 * (PI / kp).sqrt()) + lemI_tail_fixed(sv, l0, c, q)
        st = (iT / (6 * (1 - arb('0.66') * iT)) + abs(sv) * iT ** 2 + sv * sv * iT / (4 * (1 - iT))).exp()
        return (1 + de * eps) ** sv * lem * st

    V = 12; nseg = 48; tot = arb(0)
    edges = [arb(-V) + arb(2 * V) * j / nseg for j in range(nseg + 1)]
    for j in range(nseg):
        v = edges[j].union(edges[j + 1]); w = (edges[j + 1].erf() - edges[j].erf()) / 2
        tot += f(sig + cb.rsqrt(t) * v) * w
    lam_b = (1 + de * eps) ** 2
    l0_b = Qs * PI * de * (2 + de * eps) / (2 * s2)
    tv = cf.up(cb.tail_v(sig, t, arb(Tps), arb(c), de, q, lam=lam_b, l0=l0_b))
    return cf.up(tot) + tv + 1e-300


def block_err(ch, B, x1, x2, y, t, Kh):
    """e_C + e_AB, uniformly on the sub-block x in [x1, x2] (N constant, c = N + 1/2), at height y.
    delta = w_R - c is increasing in x, so delta lies in [w_R(x1) - c, w_R(x2) - c] (outward, within [-1/2, 1/2])."""
    q = 5
    xb = arb(x1).union(arb(x2)); tb = arb(t)
    Tp = xb / 2 + PI * tb / 8
    Tps = cf.down(Tp)
    N = B.N
    cs = N + 0.5
    w1 = (q * (arb(x1) / 2 + PI * tb / 8) / (2 * PI)).sqrt(); w2 = (q * (arb(x2) / 2 + PI * tb / 8) / (2 * PI)).sqrt()
    assert N >= 10 and w1.lower() >= N - 1e-12 and w2.upper() <= N + 1 + 1e-12, 'Lemma I needs c >= 10.5; N constant'
    dlo = max(-0.5, cf.down(w1 - cs)); dhi = min(0.5, cf.up(w2 - cs))
    # |M0(iT')|/|B_t| <= exp(Lemma E bound at x1): the bound decreases in x (cert_extra.gamma_ratio_log_bound,
    # as in part (a)).  (Evaluating the ratio on the x-ball directly loses a factor ~e^{pi dx/4} to the
    # dependency between the two log-Gamma terms.)
    lr = cf.up(ce.gamma_ratio_log_bound(ch, x1, y, t).exp())
    eC = 0.0
    for sig in (arb((1 - y) / 2), arb((1 + y) / 2)):
        eC += cf.up((PI ** 2 * tb / 64).exp()) * lr * Kh * gauss_fixed(ch, sig, t, cs, Tps, dlo, dhi)
    # e_A + e_B <= eps_max * F (Prop 6.1^chi), F = sup over the sub-block of
    # sum |a_n| b_n n^{-Re s*} + |gamma| sum |a_n| b_n n^{-Re w} (RSN.majorant on real x-sub-balls of width <= 1/4);
    # eps_max as in dh_upper.final_region (decreasing in x on x >= x1)
    m = max(1, int(math.ceil((x2 - x1) * 4)))
    F = max(B.majorant(arb(x1 + (x2 - x1) * i / m).union(arb(x1 + (x2 - x1) * (i + 1) / m)), arb(0), y, t)
            for i in range(m))
    lN = math.log(math.sqrt(q * x1 / (4 * math.pi)))
    eps = math.expm1(((t * t / 8) * (2 * lN + 3) ** 2 + t / 4 + 1 / 6) / (x1 / 2 - 1.25)) * 1.01
    return (eC + eps * F) * (1 + 1e-9), eC, eps * F


def x_of_w(w, t=T0, q=5):
    return 4 * math.pi * (w * w - q * t / 16) / q


def rs_edge(edge, xa, xb, t=T0, ell=2.0, J=20, rho=4.0, nb=5, tmax=1e9):
    """Certify dist(G, (-inf, 0]) > 0 on the edge for x in [xa, xb]; each N-block is cut into nb sub-blocks
    of equal delta-width (the Lemma I bound is evaluated per sub-block).  Reports failing pieces and 'cut',
    the smallest x beyond which every piece passed."""
    T1 = time.time()
    rig._init_constants()
    ch = DHChar(t)
    Kh = ce.kh_best(ch)                          # Lemma H grid bound (q = 5 <= 29)
    y = Y0 if edge == 'bottom' else YTOP
    N = int(math.floor(math.sqrt(5 * xa / (4 * math.pi) + 5 * t / 16)))
    worst = math.inf; worst_at = None; pieces = 0; fails = []; errmax = 0.0; x = xa; last_fail = None
    while x < xb:
        x2N = min(x_of_N(N + 1, t), xb)
        if x2N <= x: N += 1; continue
        B = RSN(ch, N)
        cuts = sorted(set([x, x2N] + [x_of_w(N + j / nb, t) for j in range(1, nb) if x < x_of_w(N + j / nb, t) < x2N]))
        for s1, s2 in zip(cuts[:-1], cuts[1:]):
            err, eC, eAB = block_err(ch, B, s1, s2, y, t, Kh)
            errmax = max(errmax, err)
            n = max(1, int(math.ceil((s2 - s1) / ell)))
            stack = [(s1 + (s2 - s1) * i / n, s1 + (s2 - s1) * (i + 1) / n, 0) for i in range(n)]
            while stack:
                a, b, dep = stack.pop()
                c0, h = (a + b) / 2, (b - a) / 2
                uu = np.cos(np.pi * np.arange(J + 1) / J)
                X = c0 + h * uu; X[0] = b; X[-1] = a
                v, ev = B.Ef(X, np.full(J + 1, y), t)
                xr, xi = ellipse(c0, h, rho)
                M = B.majorant(xr, xi, y, t)
                lb, info = cheb_lb(v, float(ev.max()), M, rho)
                m = lb - err; pieces += 1
                if not (m > 0):                   # also catches nan
                    if dep < 3 and lb - info['dmin'] < -0.1 * err: stack += [(a, c0, dep + 1), (c0, b, dep + 1)]
                    else:
                        fails.append((a, b, round(m, 4))); last_fail = b if last_fail is None else max(last_fail, b)
                    continue
                if m < worst: worst, worst_at = m, dict(a=a, b=b, N=N, err=err, eC=eC, eAB=eAB, **info)
        x = x2N; N += 1
        if time.time() - T1 > tmax: break
    return dict(part='rs', edge=edge, y=y, xa=xa, xb=xb, reached=x, t0=t, Kh=Kh, pieces=pieces, nfails=len(fails),
                fails=fails[:40], cut=last_fail, worst_margin=worst, worst_at=worst_at, err_max=errmax,
                ok=bool((not fails) and x >= xb and worst > 0), secs=round(time.time() - T1, 1))


# ------------------------------------------------------------------ (D) pieces
PREC_D = 80; RT_D = 24
AS, HS = 3.0, 0.6                                 # trapezoid strip half-width a and step h (H_trap)


def lMt_arb(s, t):
    lq = (arb(5) / PI).log(); sa = s + 1
    am = lq / 2 + (sa / 2).log() / 2 - 1 / (2 * sa)
    w = sa / 2
    return (t / 4) * am ** 2 + w * lq + (2 * PI).sqrt().log() + (w - arb(1) / 2) * w.log() - w


def Wof(X0, t):
    return max(20.0, math.ceil(math.sqrt(8 * t * (math.pi * X0 / 8 + 45))))


def H_direct(dh, t, z, W):
    tt = arb(t)
    g = lambda w, analytic: (-(w * w) / (4 * tt)).exp() * dh.xi((1 + w + I_ * z) / 2)
    val = acb.integral(g, -W, W, rel_tol=arb(2) ** -RT_D)
    return val / (4 * PI * tt).sqrt() + dh.tail_bound(tt, z, W)


def H_trap(dh, t, z, W, h, a, Mstrip):
    """H_t(z) by the trapezoid rule on w in [-W, W] with step h, with a rigorous error ball:
      F(w) = (4 pi t)^{-1/2} e^{-w^2/4t} xi((1 + w + iz)/2) is entire and -> 0 in every horizontal strip, so
      |int F - h sum_{k in Z} F(kh)| <= 2 M/(e^{2 pi a/h} - 1)  (Trefethen-Weideman, SIAM Rev. 56 (2014), Thm 5.1),
      M = sup_{|b| < a} int |F(w + ib)| dw <= e^{a^2/4t} sup_{|x' - Re z| <= a} (4 pi t)^{-1/2} int e^{-w^2/4t}
      |xi((1 + w - Im z)/2 + i x'/2)| dw =: e^{a^2/4t} Mstrip  (Mstrip from H_sup on that segment);
      h sum_{|kh| > W} |F(kh)| <= int_{|w| >= W} m(w) dw for the decreasing majorant m of cert_dh.DH.tail_bound,
      which is the tail_bound itself (normalised by (4 pi t)^{-1/2}, like F)."""
    tt = arb(t)
    s0 = (1 + I_ * z) / 2
    # grid offset: the trapezoid bound holds for the shifted grid {kh + off}; near x = 0 choose off so that no node
    # s_k = s0 + w_k/2 comes within 0.05 of s = 1 or s = 0 (the Hurwitz poles, also after the reflection in xi_eval)
    off = 0.0
    if abs(float(z.real.mid())) < 0.2:
        y_ = float(z.imag.mid())
        for j in range(64):
            off = h * j / 64
            if all(abs(((c - off) / h) - round((c - off) / h)) * h > 0.1 for c in (1 + y_, y_ - 1)): break
        else:
            raise AssertionError('no admissible grid offset')
    K_ = int(math.floor(W / h)) + 1
    tot = acb(0)
    for k in range(-K_, K_ + 1):
        wf = k * h + off
        if abs(wf) > W: continue
        w = arb(k) * arb(h) + arb(off)
        tot += (-(w * w) / (4 * tt)).exp() * xi_eval(dh, s0 + w / 2)
    val = tot * arb(h) / (4 * PI * tt).sqrt()
    disc = 2 * (arb(a) ** 2 / (4 * tt)).exp() * arb(Mstrip) / ((2 * PI * arb(a) / arb(h)).exp() - 1)
    # excluded nodes have |w| > W, so each h|F(w_k)| <= int_{w_k - h}^{w_k} m with m decreasing from W - h
    tail = dh.tail_bound(tt, z, W - h)
    er = disc.upper() + tail.real.upper()
    return val + acb(arb(0, er), arb(0, er))


def xi_disk_bound(dh, r=0.45, n=96):
    """sup_{|s - 1| <= r} |xi(s)|, by the maximum modulus principle: the circle |s - 1| = r covered by n balls
    (each at distance >= r - 2 pi r/n > 0.4 from the Hurwitz poles' common point s = 1; xi itself is entire,
    since sum a_r = 0)."""
    key = ('disk', r, n)
    if key in dh.__dict__: return dh.__dict__[key]
    best = arb(0); rr = arb(r) * (arb.pi() / n) * arb('1.01')
    for k in range(n):
        th = 2 * arb.pi() * k / n
        c = acb(1 + arb(r) * th.cos(), arb(r) * th.sin())
        b = acb(c.real + arb(0, rr.upper()), c.imag + arb(0, rr.upper()))
        v = dh.xi(b).abs_upper()
        if v.upper() > best.upper(): best = arb(v.upper())
    dh.__dict__[key] = arb(best.upper())
    return dh.__dict__[key]


def xi_eval(dh, s):
    """xi(s), evaluated at 1 - s when Re s < 1/2 (xi(s) = xi(1 - s)), which keeps Gamma((s + 1)/2) away from its
    poles s = -1, -3, ... (xi is entire; only the factorised evaluation is singular there)."""
    return dh.xi(s) if float(s.real.mid()) >= 0.5 else dh.xi(1 - s)


def xi_sup(dh, s, depth=0):
    """upper bound for |xi| on the complex box s.  A box straddling Re s = 1/2 is split there; a box with
    Re s < 1/2 is reflected (xi(s) = xi(1 - s)); a box meeting |s - 1| < 0.3 is split down to half-widths <= 0.05
    (then contained in |s - 1| <= 0.45) and bounded by xi_disk_bound (maximum modulus)."""
    re, im = s.real, s.imag
    if re.lower() < 0.5 - 1e-6 and re.upper() > 0.5 + 1e-6:
        lo = arb(re.lower()).union(arb(0.5)); hi = arb(0.5).union(arb(re.upper()))
        a1 = xi_sup(dh, acb(lo, im), depth + 1); a2 = xi_sup(dh, acb(hi, im), depth + 1)
        return arb(a1.upper()) if a1.upper() > a2.upper() else arb(a2.upper())
    if re.upper() <= 0.5 + 1e-6:
        s = 1 - s; re, im = s.real, s.imag
    rx, ry = float(re.rad()), float(im.rad()); mx, my = float(re.mid()), float(im.mid())
    dx = max(0.0, abs(mx - 1) - rx); dy = max(0.0, abs(my) - ry)
    if math.hypot(dx, dy) > 0.3:
        return dh.xi(s).abs_upper()
    if max(rx, ry) <= 0.05:
        return xi_disk_bound(dh)
    best = arb(0)
    for sx in ((-1, 1) if rx > 0.05 else (0,)):
        for sy in ((-1, 1) if ry > 0.05 else (0,)):
            nr = arb(mx + sx * rx / 2, rx / 2 * 1.0000001) if sx else re
            ni = arb(my + sy * ry / 2, ry / 2 * 1.0000001) if sy else im
            v = xi_sup(dh, acb(nr, ni), depth + 1)
            if v.upper() > best.upper(): best = arb(v.upper())
    return best


def H_sup(dh, t, zr, zi, W, cells=160):
    """sup_{z in zr + i zi} (4 pi t)^{-1/2} int e^{-w^2/4t} |xi((1 + w + iz)/2)| dw (+ the tail), in Arb."""
    tt = arb(t)
    edges = [-W + 2 * W * k / cells for k in range(cells + 1)]
    tot = arb(0)
    for k in range(cells):
        wb = arb(edges[k]).union(arb(edges[k + 1]))
        s = acb((1 + wb - zi) / 2, zr / 2)
        mass = ((arb(edges[k + 1]) / (2 * tt.sqrt())).erf() - (arb(edges[k]) / (2 * tt.sqrt())).erf()) / 2
        tot += xi_sup(dh, s) * mass
    return tot + dh.tail_bound(tt, acb(zr, zi), W).real.upper() * 2


def G_direct(dh, t, z, W):
    sp = (1 + z.imag - I_ * z.real) / 2
    return H_direct(dh, t, z, W) / (dh.norm * lMt_arb(sp, arb(t)).exp())


def G_majorant(dh, t, zr, zi, W, sub=(4, 2)):
    """sup |G| over the z-box zr + i zi: sup|H| / inf|B_t|, cut into sub[0] x sub[1] boxes."""
    tt = arb(t); best = 0.0
    xs = [zr.lower() + (zr.upper() - zr.lower()) * k / sub[0] for k in range(sub[0] + 1)]
    ys = [zi.lower() + (zi.upper() - zi.lower()) * k / sub[1] for k in range(sub[1] + 1)]
    for i in range(sub[0]):
        for j in range(sub[1]):
            zre = arb(xs[i]).union(arb(xs[i + 1])); zim = arb(ys[j]).union(arb(ys[j + 1]))
            Hs = H_sup(dh, t, zre, zim, W)
            sp = (1 + zim - I_ * zre) / 2
            Bl = (lMt_arb(sp, tt).real).exp().lower() * dh.norm
            assert Bl > 0
            best = max(best, cf.up(Hs / Bl))
    return best


def direct_edge(edge, a0, b0, t=T0, ell=4.0, J=24, rho=2.0, tmax=1e9):
    """edge 'bottom'/'top': x in [a0, b0] at y = Y0/YTOP;  edge 'left': y in [a0, b0] at x = 0."""
    T1 = time.time()
    ctx.prec = PREC_D
    worst = math.inf; worst_at = None; pieces = 0; fails = []; evals = 0
    n = max(1, int(math.ceil((b0 - a0) / ell)))
    stack = [(a0 + (b0 - a0) * i / n, a0 + (b0 - a0) * (i + 1) / n, 0) for i in range(n)][::-1]
    while stack:
        a, b, dep = stack.pop()
        c0, h = (a + b) / 2, (b - a) / 2
        uu = np.cos(np.pi * np.arange(J + 1) / J)
        S = c0 + h * uu; S[0] = b; S[-1] = a
        X0 = c0 if edge != 'left' else 0.0
        dh = cert_dh.DH(round(X0))
        # the ellipse box (in z) and the tail width
        er, ei = ellipse(c0, h, rho)
        if edge == 'left':
            zr, zi = ei, er                       # z = i*y(u): Re z = -h Im u, Im z = c0 + h Re u
        else:
            yv = Y0 if edge == 'bottom' else YTOP
            zr, zi = er, ei + yv
        # the trapezoid strip (half-width AS in Re z around every node of the piece)
        if edge == 'left':
            sr, si = arb(-AS).union(arb(AS)), arb(a).union(arb(b))
        else:
            sr, si = arb(a - AS).union(arb(b + AS)), arb(Y0 if edge == 'bottom' else YTOP)
        W = Wof(max(abs(float(zr.lower())), abs(float(zr.upper())), abs(float(sr.lower())), abs(float(sr.upper()))), t)
        Ms = H_sup(dh, t, sr, si, W)
        F = np.empty(J + 1, dtype=complex); evmax = 0.0
        for j, sv in enumerate(S):
            z = acb(0, sv) if edge == 'left' else acb(sv, Y0 if edge == 'bottom' else YTOP)
            Hz = H_trap(dh, t, z, W, HS, AS, Ms)
            spz = (1 + z.imag - I_ * z.real) / 2
            G = Hz / (dh.norm * lMt_arb(spz, arb(t)).exp())
            F[j] = complex(float(G.real.mid()), float(G.imag.mid()))
            evmax = max(evmax, cf.up(G.real.rad()) + cf.up(G.imag.rad()) + 2 * U * abs(F[j]))
        evals += J + 1
        M = G_majorant(dh, t, zr, zi, W, sub=(4, 2))
        lb, info = cheb_lb(F, evmax, M, rho)
        pieces += 1
        if not (lb > 0):                      # also catches nan
            if dep < 4: stack += [(c0, b, dep + 1), (a, c0, dep + 1)]
            else: fails.append((a, b, lb))
            continue
        if lb < worst: worst, worst_at = lb, dict(a=a, b=b, **info)
        if time.time() - T1 > tmax: break
    return dict(part='direct', edge=edge, a=a0, b=b0, t0=t, pieces=pieces, evals=evals, fails=fails,
                pending=len(stack), worst_margin=worst, worst_at=worst_at,
                ok=bool((not fails) and not stack and worst > 0), secs=round(time.time() - T1, 1))


def set_params(t0, y0, X):
    """Reset the rectangle parameters (the module constants read by rs_edge/direct_edge)."""
    global T0, Y0, XR, YTOP
    T0, Y0, XR = t0, y0, X
    YTOP = math.sqrt(DELTA ** 2 - 2 * T0) + 1e-9


if __name__ == '__main__':
    part, edge, a, b = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
    kw = {}; prm = {}
    for s_ in sys.argv[5:]:
        k_, v_ = s_.split('=')
        if k_ in ('t0', 'y0', 'X'): prm[k_] = float(v_)
        else: kw[k_] = float(v_) if k_ in ('ell', 'rho', 'tmax') else int(v_)
    if prm:
        set_params(prm.get('t0', T0), prm.get('y0', Y0), prm.get('X', XR))
    if edge == 'left' and b < 0: b = YTOP              # b = -1: up to Ytop
    kw['t'] = T0
    out = rs_edge(edge, a, b, **kw) if part == 'rs' else direct_edge(edge, a, b, **kw)
    out.update(y0=Y0, Ytop=YTOP, X=XR)
    print(json.dumps(out, default=str))
