"""Certificate, remaining pieces (papers/supporting/riemann-siegel-identity/derivation.md §3):

(c) Positivity used in Polymath §3.
    For every D except -19: theta_chi(y) > 0 for all y > 0, where theta_chi(y) =
    sum_n chi(n) n^a e^{-pi n^2 y/q} (a = parity); then Phi(u) ∝ e^{(1+2a)u} theta(e^{4u}) > 0.
      y in [1, 20]: Arb interval bisection (terms n <= M plus a rigorous tail);
      y >= 20:      theta >= e^{-pi y/q}(1 - sum_{n>=2} n e^{-pi(n^2-1)y/q}) > 0, the sum
                    being decreasing in y and checked at y = 20;
      y < 1:        theta(1/y) = y^{a+1/2} theta(y) (real primitive chi has root number 1).
    For D = -19 (theta_1 changes sign): H_t(iy) = int_0^inf Phi(u) e^{t u^2} cosh(yu) du > 0
    on [0, t0] x [0, 1], by acb.integral on [0, 1.6] with (t, y) as balls, plus the u > 1.6 tail.

(e) Final-region error terms.  For all x >= x_s = XC + 1/2, y in [0.2, sqrt(1-2t)], t = t0:
      |E| (e_A + e_B + e_C) <= la_max * (eAB_sup + eC_sup),
    with la_max = prod_{p in P5} (1 + b_p p^{-s_min}) (s_min = sigma at the corner),
    eC_sup from the monotone majorant of the gamma-ratio (derivation §6, Lemma E) and a
    uniform Arb bound of the Gaussian average over eps = 1/c in [0, 1/c_s], delta in [-1/2, 1/2],
    1/T' in [0, 1/T'_s]; eAB_sup = eps(x_s) * F with F the crude bound of Lemma E.
"""
import math, sys, json
import numpy as np
from flint import arb, acb, ctx
from cert_final import Char, down, up, PI
import cert_barrier as cb

ctx.prec = 96
I_ = acb(0, 1)


# ------------------------------------------------------------------ (c) theta > 0
def theta_ball(ch, y, M=None):
    if M is None:
        # terms with pi n^2 y/q > 750 are below e^-750: M = sqrt(750 q/(pi y_lo)) + 2
        M = int(math.sqrt(750 * ch.q / (math.pi * float(y.lower())))) + 2
    q, a = ch.q, ch.A
    s = arb(0)
    for n in range(1, M + 1):
        c = ch.chi(n)
        if c:
            s += c * arb(n) ** a * (-PI * n * n * y / q).exp()
    # tail n > M: sum n e^{-pi n^2 y/q} <= e^{-pi M^2 y/q} * sum_{k>=1} (M+k) e^{-pi(2Mk)y/q}
    r = (-2 * PI * M * y.lower() / q).exp()
    tail = (-PI * M * M * y.lower() / q).exp() * (M + 1) * r / (1 - r) ** 2
    return s + arb(0, 1) * tail if False else s.union(s + tail).union(s - tail)


def dominance_Y(q):
    """Smallest Y (doubling from 20) with sum_{n>=2} n e^{-pi (n^2 - 1) Y/q} < 1/2 (checked in Arb);
    the sum decreases in Y, so theta(y) > 0 for all y >= Y."""
    Y = 20.0
    M = 400
    while True:
        s = arb(0); y = arb(Y)
        for n in range(2, M):
            s += arb(n) * (-PI * (n * n - 1) * y / q).exp()
        # n >= M: sum_k (M+k) e^{-pi((M+k)^2-1)Y/q} <= e^{-pi(M^2-1)Y/q} (M/(1-r) + r/(1-r)^2), r = e^{-2 pi M Y/q}
        r = (-2 * PI * M * y / q).exp()
        s += (-PI * (M * M - 1) * y / q).exp() * (M / (1 - r) + r / (1 - r) ** 2)
        if s < arb(1) / 2:
            return Y, up(s)
        Y *= 2


def theta_positive(ch, max_cells=400):
    q = ch.q
    Y, sY = dominance_Y(q)
    stack = [(1.0, Y)]; cells = 0; mn = 9.0
    while stack:
        a, b = stack.pop()
        v = theta_ball(ch, arb(a).union(arb(b)))
        lo = float(v.lower()); cells += 1
        if lo > 0:
            # normalise by the n = 1 term to report a scale-free margin
            mn = min(mn, lo / float((-PI * a / q).exp().upper()))
            continue
        if b - a < 1e-6 or cells > max_cells:
            return dict(ok=False, at=(a, b), cells=cells)
        m = (a + b) / 2; stack += [(a, m), (m, b)]
    return dict(ok=True, cells=cells, min_rel=mn, Y=Y, tailY=sY)


def H_iy_positive(ch, T0, max_boxes=400):
    """H_t(iy) = int_0^inf Phi(u) e^{t u^2} cosh(yu) du > 0 on [0, T0] x [0, 1], where
    Phi(u) = e^{(1+2a)u} theta_chi(e^{4u}) (up to a positive constant).  Adaptive bisection over
    (t, y) boxes given as Arb balls; each box: acb.integral on [0, U] plus a rigorous tail.
    U = log(300 q)/4, so that pi e^{4U}/q = 300 pi.  Tail (v = u - U >= 0): with
    e^{4u} >= e^{4U}(1 + 4v + 8v^2), |Phi| e^{tu^2} cosh(yu) <= sum_n n^a C_n exp(A_n v - B_n v^2),
    C_n = e^{(1+2a+y)U + tU^2} e^{-pi n^2 e^{4U}/q}, A_n = 1+2a+y+2tU - 4 pi n^2 e^{4U}/q,
    B_n = 8 pi n^2 e^{4U}/q - t > 0, and int_0^inf e^{Av - Bv^2} dv <= sqrt(pi/B) e^{max(A,0)^2/(4B)}.
    Terms with pi n^2 e^{4u}/q > 800 on [0, U] are dropped only in the tail sum (bounded the same way)."""
    q, a = ch.q, ch.A
    U = math.log(300 * q) / 4
    Nmax = int(math.sqrt(800 * q / math.pi)) + 2          # at u = 0 the n-th term is <= e^{-800} beyond
    ns = [(n, ch.chi(n)) for n in range(1, Nmax + 1) if ch.chi(n) != 0]

    def phi(u):
        e4 = (4 * u).exp()
        s = acb(0)
        for n, c in ns:
            s += c * arb(n) ** a * (-PI * n * n * e4 / q).exp()
        return ((1 + 2 * a) * u).exp() * s

    def tail(tb, yb):
        t = arb(tb.upper()); y = arb(yb.upper()); Ub = arb(U); e4U = (4 * Ub).exp()
        tot = arb(0)
        for n in range(1, 4 * Nmax):
            C = ((1 + 2 * a + y) * Ub + t * Ub * Ub).exp() * (-PI * n * n * e4U / q).exp()
            A = 1 + 2 * a + y + 2 * t * Ub - 4 * PI * n * n * e4U / q
            B = 8 * PI * n * n * e4U / q - t
            Ap = A if A > 0 else arb(0)
            tot += arb(n) ** a * C * (PI / B).sqrt() * (Ap * Ap / (4 * B)).exp()
        # n >= 4 Nmax: the terms are < e^{-800 * 16}, far below the arb(10)^-300 slack
        # the terms n > Nmax dropped from phi on [0, U]: |sum| <= sum_{n>Nmax} n e^{-pi n^2/q}
        # (e^{4u} >= 1), times int_0^U e^{(1+2a)u} e^{t u^2} cosh(yu) du <= U e^{(3 + 1 + t U) U}
        r = (-2 * PI * Nmax / q).exp()
        drop = (-PI * Nmax * Nmax / q).exp() * (Nmax / (1 - r) + r / (1 - r) ** 2)
        drop *= Ub * ((4 + t * Ub) * Ub).exp()
        return tot + drop + arb(10) ** -300

    mn = 9.0; boxes = 0
    stack = [(0.0, T0, 0.0, 1.0)]
    while stack:
        tlo, thi, ylo, yhi = stack.pop()
        tb = arb(tlo).union(arb(thi)); yb = arb(ylo).union(arb(yhi))
        f = lambda u, analytic: phi(u) * (tb * u * u).exp() * (yb * u).cosh()
        val = acb.integral(f, 0, U).real
        lo = float((val - tail(tb, yb)).lower()); boxes += 1
        if lo > 0:
            mn = min(mn, lo); continue
        if boxes > max_boxes:
            return dict(ok=False, at=(tlo, thi, ylo, yhi), lo=lo, boxes=boxes)
        if (thi - tlo) / T0 >= (yhi - ylo):
            m = (tlo + thi) / 2; stack += [(tlo, m, ylo, yhi), (m, thi, ylo, yhi)]
        else:
            m = (ylo + yhi) / 2; stack += [(tlo, thi, ylo, m), (tlo, thi, m, yhi)]
    return dict(ok=True, boxes=boxes, min_lower=mn, U=U, Nmax=Nmax)


# ------------------------------------------------------------------ (e) final-region errors
def gamma_ratio_log_bound(ch, x, y, t):
    """Lemma E: log(|M0(iT')|/|M_t(s+)|) <= -(1+y)/2 log w_R(T') + (1+y)/(4T') + (pi/4 + 1/T') pi t/8
       - (t/4)[(1/2 log(qx/4pi) - 1/x)^2 - (pi/4 + 1/x)^2]; decreasing in x and in y."""
    q = ch.q; x = arb(x); y = arb(y); t = arb(t)
    Tp = x / 2 + PI * t / 8
    wR = (q * Tp / (2 * PI)).sqrt()
    ra = (q * x / (4 * PI)).log() / 2 - 1 / x
    assert ra > 0
    return (-(1 + y) / 2 * wR.log() + (1 + y) / (4 * Tp) + (PI / 4 + 1 / Tp) * PI * t / 8
            - (t / 4) * (ra ** 2 - (PI / 4 + 1 / x) ** 2))


def lemI_tail(sv, l0, eps, q):
    """Lemma I's Tail := int_{|u| > rho c} e^{|sigma| m(u) + l0 u - Qs pi u^2/8} du (rho = 1/4, Qs = 4/q: on the line
    c + u e^{i pi/4}, |e^{i pi w^2/q}| = e^{-(pi/q)(sqrt2 c u + u^2)}), for every c >= 1/eps_upper (eps a ball in
    [0, 1/c_s]).  m(u) <= |u|/c + 0.35, so with U = c/4, a = |sigma|/c + |l0|:
    Tail <= 2 e^{0.35|sigma|} e^{aU - Qs pi U^2/8}/(Qs pi U/4 - a)  (needs Qs pi U/4 > a).
    The bound decreases in U once Qs pi U/4 > 2a, so U = c_s/4 is the worst case (checked below)."""
    Qs = arb(4) / q
    cmin = 1 / eps.upper()
    U = cmin / 4
    s = abs(sv).upper(); la = abs(l0).upper()
    a = s / cmin + la
    assert Qs * PI * U / 4 > 2 * a
    return 2 * (arb('0.35') * s).exp() * (a * U - Qs * PI * U * U / 8).exp() / (Qs * PI * U / 4 - a)


def gauss_uniform(ch, sig, t, cs, Tps):
    """Uniform bound for the Gaussian average of f over c >= cs, |delta| <= 1/2, T' >= Tps,
    written in eps = 1/c, delta, iT = 1/T' (derivation §6, Lemma E)."""
    q = ch.q; Qs = arb(4) / q; s2 = arb(2).sqrt(); rho = arb(1) / 4
    eps = arb(0).union(1 / arb(cs)); de = arb(-0.5).union(arb(0.5)); iT = arb(0).union(1 / arb(Tps))
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
        lem = G0 + A * (l * l / (4 * kp)).exp() * 4 * (1 / kp ** 2 + abs(l / (2 * kp)) ** 3 * (PI / kp).sqrt()) + lemI_tail(sv, l0, eps, q)
        st = (iT / (6 * (1 - arb('0.66') * iT)) + abs(sv) * iT ** 2 + sv * sv * iT / (4 * (1 - iT))).exp()
        return (1 + de * eps) ** sv * lem * st

    V = 12; nseg = 48; tot = arb(0)
    edges = [arb(-V) + arb(2 * V) * j / nseg for j in range(nseg + 1)]
    for j in range(nseg):
        v = edges[j].union(edges[j + 1]); w = (edges[j + 1].erf() - edges[j].erf()) / 2
        tot += f(sig + cb.rsqrt(t) * v) * w
    # |v| > 12: the tail majorant with lam and l0 as balls over (delta, eps), and c = cs,
    # Tp = Tps (A0, A1, K1, K2 are largest there; Tp/c^2 is ~ 2 pi lam/q, uniformly bounded)
    lam_b = (1 + de * eps) ** 2
    l0_b = Qs * PI * de * (2 + de * eps) / (2 * s2)
    tv = up(cb.tail_v(sig, t, arb(Tps), arb(cs), de, q, lam=lam_b, l0=l0_b))
    return up(tot) + tv + 1e-300


def kh_best(ch):
    """Grid K_h (Lemma H) for q <= 29, analytic (Lemma H') above; both rigorous upper bounds."""
    return cb.kh_rig(ch)[0] if ch.q <= 29 else cb.kh_analytic(ch)[0]


def final_errors(ch, XC, T0, P5, Kh=None):
    t = T0; ytop = math.sqrt(1 - 2 * t)
    xs = XC + 0.5
    Tps = xs / 2; cs = math.floor(math.sqrt(ch.q * Tps / (2 * math.pi))) - 1
    if Kh is None:
        Kh = kh_best(ch)
    lr = up(gamma_ratio_log_bound(ch, xs, 0.2, t).exp())
    eC = 0.0
    for sig in (arb((1 - ytop) / 2).union(arb(0.4)), arb(0.6).union(arb((1 + ytop) / 2))):
        eC += up((PI ** 2 * arb(t) / 64).exp()) * lr * Kh * gauss_uniform(ch, sig, t, cs, Tps)
    # la_max: sigma at the corner (x_s, 0.2, t) is the minimum over the final region
    s_min = down(ch.sig_gam(arb(xs), arb(0.2), arb(t))[0])
    la = 1.0
    for p in P5:
        la *= 1 + up(((arb(t) / 4) * arb(p).log() ** 2 - s_min * arb(p).log()).exp())
    # e_A + e_B: eps(x_s) * F; F <= sqrt(4 pi/t)(e^{0.16/t} + G e^{((1+ytop)/2)^2/t}) + 2 (Lemma E)
    G = up(ch.sig_gam(arb(xs), arb(0.2), arb(t))[2])
    F = math.sqrt(4 * math.pi / t) * (math.exp(0.16 / t) + G * math.exp(((1 + ytop) / 2) ** 2 / t)) + 2
    N = math.sqrt(ch.q * xs / (4 * math.pi)); lN = math.log(N)
    eps = math.expm1(((t * t / 8) * (2 * lN + 3) ** 2 + t / 4 + 1 / 6) / (xs / 2 - 1.25)) * 1.01
    return dict(Kh=Kh, gamma_ratio=lr, eC_sup=eC, la_max=la * (1 + 1e-12), eAB_sup=eps * F,
                total=la * (eC + eps * F) * (1 + 1e-9))


if __name__ == '__main__':
    import rig
    rig._init_constants()
    D = int(sys.argv[1]); T0 = float(sys.argv[2]); XC = float(sys.argv[3])
    ch = Char(D, T0)
    pr = [p for p in [2, 3, 5, 7, 11, 13, 17] if ch.chi(p) != 0]
    out = dict(D=D)
    out['positivity'] = H_iy_positive(ch, T0) if D == -19 else theta_positive(ch)
    out['final_errors'] = final_errors(ch, XC, T0, pr[:5])
    print(json.dumps(out))
