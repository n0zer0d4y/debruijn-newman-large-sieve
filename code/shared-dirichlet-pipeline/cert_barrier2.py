"""Certificate, part (iii), fast version: the barrier by tensor Chebyshev interpolation.

For a t-box [ta, tb] the rectangle is R' = [XC - 1/2, XC + 1/2] x [yl(tb), 1] (it contains R_t for
every t in the box).  On each edge, with u in [-1, 1] the edge parameter and v in [-1, 1] the rescaled
time, g(u, v) = E_t f_t(z(u), t(v)) = E(s*) (A(s*) + gamma B(w)), where s* = s+ + (t/2) alpha(s+) and
w = s- + (t/2) alpha(s-) (= conj(s*) - y + kappa for real data), is analytic in (u, v).

  p = tensor Chebyshev interpolant at (J+1) x (K+1) Chebyshev points of the 2nd kind, node values
      from cert_barrier.Barrier.Ef (verified, with error ev);
  |g - p| <= 4 M1 rho1^-J/(rho1 - 1) + Lam_J 4 M2 rho2^-K/(rho2 - 1) + Lam_J Lam_K max(ev) + e_float,
      (Trefethen, ATAP Thm 8.2 in each variable; Lam_n <= (2/pi) log(n+1) + 1 the Lebesgue constant),
      with M1 = sup|g| on E_rho1 x [-1,1], M2 = sup|g| on [-1,1] x E_rho2, bounded termwise in Arb:
      |A| <= sum b_n(Re t_hi) n^{-Re s*_lo}, |B| <= sum b_n n^{-Re w_lo}, |gamma| <= exp(Re(lMt(s-) - lMt(s+)))_hi,
      |E| <= prod (1 + b_p p^{-Re s*_lo});
  dist(p, (-inf, 0]) is bounded below on a grid in the angles (u = cos theta), with |d/dtheta cos(j theta)| <= j
      for the Lipschitz step.
Then dist(E H_t/B_t, (-inf,0]) >= dist(p) - |g - p| - sup|E| (e_C + e_A + e_B) > 0 on the boundary of R' for
every t in the box, so the winding number is 0 and H_t has no zeros in R' (derivation §3.4 / §6.1).
"""
import sys, json, math, time
import numpy as np
from flint import arb, acb, ctx
import rig
from cert_final import Char, down, up, PI
import cert_barrier as cb

ctx.prec = 96
U = rig.U
I_ = acb(0, 1)


def cheb_pts(n):
    return np.cos(np.pi * np.arange(n + 1) / n)          # 2nd kind, from +1 down to -1


def cheb_coeffs(F):
    """F[j, k] = g(u_j, v_k) at 2nd-kind points -> coefficients c[j, k] of sum c T_j(u) T_k(v)."""
    def dct1(n):
        j = np.arange(n + 1)
        M = np.cos(np.pi * np.outer(j, j) / n) * (2.0 / n)
        M[:, 0] *= 0.5; M[:, -1] *= 0.5
        M[0, :] *= 0.5; M[-1, :] *= 0.5
        return M
    J, K = F.shape[0] - 1, F.shape[1] - 1
    return dct1(J) @ F @ dct1(K).T


def lebesgue(n):
    return 2 / math.pi * math.log(n + 1) + 1


def ellipse_box(center, half, rho):
    """Real and imaginary ranges of center + half * E_rho (E_rho the Bernstein ellipse)."""
    a = half * (rho + 1 / rho) / 2; b = half * (rho - 1 / rho) / 2
    return arb(center - a).union(arb(center + a)), arb(-b).union(arb(b))


class Box:
    def __init__(self, B, ch, XC, T0):
        self.B, self.ch, self.XC, self.T0 = B, ch, XC, T0
        self.LN = B.LN
        self.LP = np.array([float(arb(p).log().mid()) for p in B.P])

    def majorant(self, xr, xi, yr, yi, tr, ti):
        """sup |g| over x in xr + i xi, y in yr + i yi, t in tr + i ti (all arb intervals)."""
        ch = self.ch
        x = acb(xr, xi); y = acb(yr, yi); t = acb(tr, ti)
        sp = (1 + y - I_ * x) / 2; sm = (1 - y + I_ * x) / 2
        ss = sp + (t / 2) * ch.am(sp)
        w = sm + (t / 2) * ch.am(sm)
        slo = down(ss.real); wlo = down(w.real)
        glog = (ch.lMt(sm, t) - ch.lMt(sp, t)).real
        ghi = up(glog.exp())
        thi = up(tr)
        L = self.LN
        a, ra = rig.vexp(thi / 4 * L * L - slo * L, 6 * U * (abs(thi) / 4 * L * L + abs(slo) * L))
        b, rb = rig.vexp(thi / 4 * L * L - wlo * L, 6 * U * (abs(thi) / 4 * L * L + abs(wlo) * L))
        s = 1 + rig.gamma_n(len(L) + 4)
        MA = float(np.sum(a)) * (1 + float(np.max(ra))) * s
        MB = float(np.sum(b)) * (1 + float(np.max(rb))) * s
        LP = self.LP
        e, re = rig.vexp(thi / 4 * LP * LP - slo * LP, 6 * U * (abs(thi) / 4 * LP * LP + abs(slo) * LP))
        ME = float(np.prod(1 + e * (1 + float(np.max(re))))) * (1 + 16 * U)
        return ME * (MA + ghi * MB) * (1 + 1e-12)

    def edge_geom(self, edge, yl, lo=-1.0, hi=1.0):
        """Edge as z(u) = z0 + half*u on the sub-interval [lo, hi] of the full edge parameter."""
        XC = self.XC
        if edge in ('bottom', 'top'):
            y = yl if edge == 'bottom' else 1.0
            a, b = XC - 0.5, XC + 0.5
            c0 = (a + b) / 2 + (lo + hi) / 2 * (b - a) / 2; h = (hi - lo) / 2 * (b - a) / 2
            return 'x', c0, h, y
        x = XC - 0.5 if edge == 'left' else XC + 0.5
        a, b = yl, 1.0
        c0 = (a + b) / 2 + (lo + hi) / 2 * (b - a) / 2; h = (hi - lo) / 2 * (b - a) / 2
        return 'y', c0, h, x

    def check(self, edge, yl, ta, tb, J=20, K=10, rho1=4.0, rho2=4.0, lo=-1.0, hi=1.0, grid=(600, 120)):
        kind, c0, h, other = self.edge_geom(edge, yl, lo, hi)
        uu = cheb_pts(J); vv = cheb_pts(K)
        tc, th_ = (ta + tb) / 2, (tb - ta) / 2
        F = np.empty((J + 1, K + 1), dtype=complex); evmax = 0.0
        for k, v in enumerate(vv):
            t = tc + th_ * v
            zs = c0 + h * uu
            if kind == 'x':
                val, ev = self.B.Ef(zs, np.full(J + 1, other), t)
            else:
                val, ev = self.B.Ef(np.full(J + 1, other), zs, t)
            F[:, k] = val; evmax = max(evmax, float(ev.max()))
        C = cheb_coeffs(F)
        # majorants on the two complexified domains
        zr1, zi1 = ellipse_box(c0, h, rho1)
        tr1 = arb(ta).union(arb(tb)); ti1 = arb(0)
        zr2 = arb(c0 - h).union(arb(c0 + h)); zi2 = arb(0)
        tr2, ti2 = ellipse_box(tc, th_, rho2)
        if kind == 'x':
            M1 = self.majorant(zr1, zi1, arb(other), arb(0), tr1, ti1)
            M2 = self.majorant(zr2, zi2, arb(other), arb(0), tr2, ti2)
        else:
            M1 = self.majorant(arb(other), arb(0), zr1, zi1, tr1, ti1)
            M2 = self.majorant(arb(other), arb(0), zr2, zi2, tr2, ti2)
        LJ, LK = lebesgue(J), lebesgue(K)
        e_int = 4 * M1 * rho1 ** -J / (rho1 - 1) + LJ * 4 * M2 * rho2 ** -K / (rho2 - 1)
        e_node = LJ * LK * evmax
        # float error of the coefficients and of evaluating p (derived, not heuristic):
        # each c_jk is a (J+1)(K+1)-term sum with factors |cos| <= 1 and weights <= 2/J, 2/K, so
        # |dc_jk| <= gamma_{J+K+6} * (4/(J K)) * sum|F|; evaluating sum c T T adds gamma_{J+K+6} sum|c|
        gJK = rig.gamma_n(J + K + 6)
        e_float = (gJK * 4.0 / (J * K) * float(np.abs(F).sum()) * (J + 1) * (K + 1)
                   + gJK * float(np.abs(C).sum())) * 2
        # dist(p) on a grid + Markov-Lipschitz step
        nu, nv = grid
        tu = np.linspace(0, math.pi, nu); tv = np.linspace(0, math.pi, nv)   # u = cos(theta)
        Tu = np.cos(np.outer(tu, np.arange(J + 1))); Tv = np.cos(np.outer(tv, np.arange(K + 1)))
        P = Tu @ C @ Tv.T
        d = np.where(P.real >= 0, np.abs(P), np.abs(P.imag))
        absC = np.abs(C)
        jj = np.arange(J + 1)[:, None]; kk = np.arange(K + 1)[None, :]
        # |dp/dtheta_u| <= sum j |c|, |dp/dtheta_v| <= sum k |c|  (d/dtheta cos(j theta) = -j sin)
        Lu = float(np.sum(absC * jj)); Lv = float(np.sum(absC * kk))
        e_grid = Lu * (math.pi / (nu - 1)) / 2 + Lv * (math.pi / (nv - 1)) / 2
        dmin = float(d.min())
        lb = dmin - e_grid - e_int - e_node - e_float
        return lb, dict(dmin=dmin, e_grid=e_grid, e_int=e_int, e_node=e_node, M1=M1, M2=M2)


def run(D, T0, XC, tau0=0.04, J=20, K=10):
    T1 = time.time()
    rig._init_constants()
    ch = Char(D, T0)
    import cert_extra
    Kh = cert_extra.kh_best(ch); c0 = None
    B = cb.Barrier(ch, XC, T0)
    bx = Box(B, ch, XC, T0)
    yl_ = lambda t: math.sqrt(0.04 + 2 * (T0 - t))
    t = 0.0; worst = 9.0; boxes = 0; evals = 0; eCmax = 0.0; Emax = 0.0; detail = None
    while t < T0 - 1e-12:
        tau = min(tau0, T0 - t)
        while True:
            tb = min(t + tau, T0); yl = yl_(tb)
            Ls, Lt, Em, Fm, aa = B.lip(yl, t, tb)
            eC = cb.eC_box(ch, XC - 0.5, XC + 0.5, yl, 1.0, t, tb, arb(Kh))
            err = Em * (eC + B.eps_AB(yl, t, tb, aa, Fm))
            ok = True; mg = 9.0
            for edge in ('bottom', 'top', 'left', 'right'):
                lb, info = bx.check(edge, yl, t, tb, J, K)
                boxes += 1; evals += (J + 1) * (K + 1)
                m = lb - err
                if m < mg: mg, det = m, (edge, info)
                if m <= 0: ok = False; break
            if ok or tau < 1e-4: break
            tau /= 2
        if not ok:
            return dict(D=D, t0=T0, barrier_ok=False, barrier_to=t, fail=det, secs=round(time.time() - T1, 1))
        worst = min(worst, mg); eCmax = max(eCmax, eC); Emax = max(Emax, Em)
        if detail is None or mg <= worst: detail = det
        t = tb
    return dict(D=D, t0=T0, XC=XC, N=B.N, Kh_rig=Kh, c0=c0, method='chebyshev', J=J, K=K, eC_max=eCmax,
                barrier_ok=True, barrier_to=t, boxes=boxes, evals=evals, barrier_worst=worst,
                worst_detail=detail, supE=Emax, secs=round(time.time() - T1, 1))


if __name__ == '__main__':
    D = int(sys.argv[1]); T0 = float(sys.argv[2]); XC = float(sys.argv[3])
    print(json.dumps(run(D, T0, XC), default=str))
