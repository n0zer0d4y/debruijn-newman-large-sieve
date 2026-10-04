"""Rouche certificate for a non-real zero of H_t for f_theta in the q = 7 odd Theorem A' family.

Exact family (Arb balls, so the certificate covers the exact coefficients):
  chi(n) = (n/7) (odd real character, tau = i sqrt 7);  chi6(3^k) = e^{i pi k/3} (odd, order 6);
  e^{i psi} = tau(conj chi6)/(i sqrt 7),  b(n) = sqrt2 Re(e^{i psi/2} chi6(n)) * sgn,
  so that sum_m b(m) e^{2 pi i m n/7} = i sqrt 7 b(n) (checked below), b is orthogonal to chi and |b| = |chi|;
  sgn is chosen to match x7_scan.basis() (the float SVD direction used in the scans).
  a_theta = cos(theta) chi + sin(theta) b  (theta given as an exact rational).
xi(s) = (7/pi)^{(s+1)/2} Gamma((s+1)/2) f(s) * norm satisfies xi(s) = xi(1-s) (Theorem A' / DFT self-duality),
so H_t(z) = (4 pi t)^{-1/2} int e^{-w^2/4t} xi((1 + w + iz)/2) dw is real entire and even, and a non-real zero
of H_{t1} gives Lambda(f_theta) > t1 (the same argument as papers/03-davenport-heilbronn-lambda #1).
The Rouche step is cert_dh2.certify's, run with this xi and a tail bound adapted to q = 7:
  |f(s)| <= zeta(sigma') C_a for sigma' >= 2, with C_a = sum_r |a_r|, and (7/pi) in place of (5/pi).
Usage: python x7_cert.py theta_num theta_den t cx cy [R= r= delta=]
"""
import sys, json, math, time
import os as _os
sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '03-davenport-heilbronn'))
from flint import arb, acb, ctx
import cert_dh, cert_dh2

ctx.prec = 160
PI = arb.pi()
I_ = acb(0, 1)
Q = 7


def family(theta):
    chi = [0, 1, 1, -1, 1, -1, -1]
    w = (I_ * PI / 3).exp()
    chi6 = [acb(0)] * Q
    g = 1
    for k in range(6):
        chi6[g] = w ** k
        g = (g * 3) % Q
    tau_bar = sum(chi6[n].conjugate() * (2 * PI * I_ * n / Q).exp() for n in range(1, Q))
    e_psi = tau_bar / (I_ * arb(Q).sqrt())
    half = (e_psi.log() / 2).exp()                         # e^{i psi/2}, principal branch
    b = [arb(0)] + [arb(2).sqrt() * (half * chi6[n]).real for n in range(1, Q)]
    # sign convention of x7_scan.basis(): b(1) < 0
    if b[1] > 0:
        b = [-x for x in b]
    a = [arb(theta).cos() * chi[n] + arb(theta).sin() * b[n] for n in range(Q)]
    # checks: DFT self-duality of b, orthogonality, norm
    res = max(abs(sum(b[m] * (2 * PI * I_ * m * n / Q).exp() for m in range(Q)) - I_ * arb(Q).sqrt() * b[n]).upper()
              for n in range(Q))
    orth = abs(sum(b[n] * chi[n] for n in range(Q))).upper()
    nb = sum(b[n] ** 2 for n in range(Q)) - 6
    return a, b, dict(selfdual_residual=float(res), orth=float(orth), norm_minus_6=float(abs(nb).upper()))


class Fam(cert_dh.DH):
    def __init__(self, X0, a, W=20.0):
        super().__init__(X0)
        self.a = a
        self.C = sum(abs(x) for x in a[1:])
        self.W = W

    def H(self, t, z, W=None):
        return super().H(t, z, W=self.W if W is None else W)

    def f(self, s):
        tot = acb(0)
        for r in range(1, Q):
            tot += self.a[r] * s.zeta(acb(arb(r) / Q))
        return tot * acb(Q) ** (-s)

    REFLECT = False                                         # True: evaluate xi(1 - s) when Re s < 1/2 (exact identity)

    def xi(self, s):
        if self.REFLECT and s.real.mid() < 0.5:
            s = 1 - s                                       # xi(s) = xi(1 - s): real DFT-self-dual a, root number 1
        return (acb(Q) / PI) ** ((s + 1) / 2) * ((s + 1) / 2).gamma() * self.f(s) * self.norm

    TAIL = 1                                                # 1: the cert_dh-style bound; 2: the Stirling bound below

    def tail_bound(self, t, z, W):
        """As cert_dh.DH.tail_bound, with (7/pi) and C_a = sum |a_r| >= max |a_n| in place of (5/pi) and 1 + kappa.
        TAIL = 2: for |w| >= W put x = (|w| + B)/4 >= (sigma' + 1)/2, B = |Im z| + 3, sigma' = max(Re s, 1 - Re s) >= 2.
        By xi(s) = xi(1-s), |Gamma(u + iv)| <= Gamma(u), |f| <= C_a zeta(2) and Stirling
        log Gamma(x) <= (x - 1/2) log x - x + log(2 pi)/2 + 1/(12x) (x >= 1), the integrand is <= e^{g(w)} with
        g(w) = -w^2/(4t) + log(norm C_a zeta(2)) + x log(7/pi) + (x - 1/2) log x - x + log(2 pi)/2 + 1/(12x).
        g' <= -w/(2t) + log(7x/pi)/4 and g'' <= -1/(2t) + 5/(48x) < 0 once x > 5t/24, so for g'(W) < 0,
        int_{|w|>W} e^{g}/sqrt(4 pi t) <= 2 e^{g(W)}/(|g'(W)| sqrt(4 pi t))."""
        if self.TAIL == 2:
            B = abs(z.imag).upper() + 3
            t = arb(t)
            x = (arb(W) + B) / 4
            assert x >= 1 and arb(W) >= abs(z.imag).upper() + 5    # sigma' >= (W - |Im z| - 1)/2 >= 2, as zeta(2) needs
            assert x > 5 * t / 24
            gp = -arb(W) / (2 * t) + (7 * x / PI).log() / 4
            assert gp < 0
            G = (self.norm * self.C * PI ** 2 / 6).log() + x * (arb(7) / PI).log() + (x - 0.5) * x.log() - x \
                + (2 * PI).log() / 2 + 1 / (12 * x)
            gW = -(arb(W) ** 2) / (4 * t) + G
            tb = 2 * gW.exp() / (-gp) / (4 * PI * t).sqrt()
            return acb(arb(0, tb.upper()), arb(0, tb.upper()))
        t = arb(t)
        B = abs(z.imag).upper() + 3
        x = (arb(W) + B) / 4
        cW = x * (Q * x / PI).log()
        assert cW < arb(W) ** 2 / (8 * t)
        assert arb(W) / (4 * t) > ((Q * x / PI).log() + 1) / 4
        # divided by (4 pi t)^{1/2}: H normalises the integral
        tb = 2 * self.norm * self.C * (PI ** 2 / 6) * (4 * t / W) * (-(arb(W) ** 2) / (8 * t)).exp() / (4 * PI * t).sqrt()
        return acb(arb(0, tb.upper()), arb(0, tb.upper()))


def certify(theta, t, cx, cy, R=0.05, r=2e-4, delta=1e-4, newton=6, W=20.0):
    a, b, chk = family(theta)
    dh = Fam(round(cx), a, W)
    orig, orig_def = cert_dh2.DH, cert_dh2.majorant.__defaults__
    cert_dh2.DH = lambda X0: dh                             # reuse cert_dh2.certify with this xi
    cert_dh2.majorant.__defaults__ = (float(W), int(20 * W))  # (W, cells): same cell width as W = 20, 400 cells
    try:
        out = cert_dh2.certify(t, cx, cy, R=R, r=r, delta=delta, newton=int(newton))
    finally:
        cert_dh2.DH = orig; cert_dh2.majorant.__defaults__ = orig_def
    out['W'] = W
    out['family_checks'] = chk
    out['coefficients_mid'] = [float(x.mid()) for x in a]
    return out


def certify_hp(theta, t, cx, cy, R=0.2, r=4e-15, delta=1e-20, newton=4, W=66.0, rt=110, tail=1, cpu=20):
    """cert_dh2.certify with (i) the Newton centre kept as an Arb midpoint at the working precision (no float
    rounding), (ii) H by acb.integral with rel_tol 2^-rt, and (iii) the majorant radius padded by 1e-9 relative, so
    that its box around the float image of c contains the disk |z - c| <= R about the exact c.  The Rouche test is
    unchanged: |H'(c) - d| r + M2 r^2/2 < |d| r - |h0|, with M3 <= 6M/R^3 and M2 <= 2M/(R - r)^2."""
    T0 = time.time()
    a, b, chk = family(theta)
    dh = Fam(round(cx), a, W)
    dh.TAIL = int(tail)
    dh.REFLECT = int(tail) == 2                                # the new mode also evaluates xi by reflection
    tt = arb(t)

    def Hq(z):
        g = lambda w, analytic: (-(w * w) / (4 * tt)).exp() * dh.xi((1 + w + I_ * z) / 2)
        val = acb.integral(g, -dh.W, dh.W, rel_tol=arb(2) ** -rt)
        return val / (4 * PI * tt).sqrt() + dh.tail_bound(tt, z, dh.W)

    c = acb(cx, cy); dl = arb(delta)
    for _ in range(int(newton)):
        h = Hq(c); d = (Hq(c + dl) - Hq(c - dl)) / (2 * dl)
        c = c - h / d
        c = acb(c.real.mid(), c.imag.mid())
    h0 = Hq(c); d = (Hq(c + dl) - Hq(c - dl)) / (2 * dl)
    M = cert_dh2.majorant(dh, t, c, R * (1 + 1e-9), W=float(W), cells=int(cpu * W))
    M3 = 6 * M / arb(R) ** 3; M2 = 2 * M / (arb(R) - r) ** 2
    lhs = M3 * dl ** 2 / 6 * r + M2 * arb(r) ** 2 / 2
    rhs = d.abs_lower() * r - h0.abs_upper()
    ok = bool(lhs < rhs) and (r < -float(c.imag.mid()))
    return dict(t=t, center=[float(c.real.mid()), float(c.imag.mid())], r=r, R=R, delta=delta, W=W, rel_tol_bits=rt, tail=int(tail), cells=int(cpu * W),
                H_center_abs_upper=float(h0.abs_upper()), d_abs_lower=float(d.abs_lower()), M=float(M.upper()),
                lhs=float(lhs.upper()), rhs=float(rhs.lower()), ok=ok, family_checks=chk,
                coefficients_mid=[float(x.mid()) for x in a], secs=round(time.time() - T0, 1))


if __name__ == '__main__':
    th = arb(int(sys.argv[1])) / int(sys.argv[2])
    t = float(sys.argv[3]); cx = float(sys.argv[4]); cy = float(sys.argv[5])
    kw = {}
    for s_ in sys.argv[6:]:
        k_, v_ = s_.split('='); kw[k_] = float(v_)
    T = time.time()
    hp = kw.pop('hp', 0)
    out = certify_hp(th, t, cx, cy, **kw) if hp else certify(th, t, cx, cy, **kw)
    out['theta'] = '%s/%s' % (sys.argv[1], sys.argv[2]); out['secs_total'] = round(time.time() - T, 1)
    if out.get('ok'):
        out['conclusion'] = 'H_t has a non-real zero for f_theta; hence Lambda(f_theta) > %g' % t
    print(json.dumps(out, default=str))
