"""a certified strictly positive lower bound for the de Bruijn-Newman constant of the
Davenport-Heilbronn function.

f(s) = sum a_n n^{-s}, a_n = (1, k, -k, -1, 0) for n = 1..5 mod 5, k = (sqrt(10 - 2 sqrt5) - 2)/(sqrt5 - 1)
(Davenport-Heilbronn; Balanzario-Sanchez-Ortiz, Math. Comp. 76 (2007) 2045, eq. (1)).
xi(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) f(s) satisfies xi(s) = xi(1-s) and is real on Re s = 1/2.
H_0(z) := xi((1 + iz)/2) is real entire, even, of order 1; Phi(u) = e^{3u} theta(e^{4u}) up to a constant.
H_t(z) = int Phi(u) e^{t u^2} e^{izu} du = (4 pi t)^{-1/2} int e^{-w^2/(4t)} H_0(z - iw) dw.

Claim: if H_{t1} has a non-real zero, then Lambda_DH > t1.  (For t >= Lambda, H_t is in the
Laguerre-Polya class; the backward heat operator exp(-s d^2/dz^2) preserves LP for s >= 0, which needs
no positivity of Phi.)  The certificate: a rigorous winding number 1 of H_{t1} around a box in Im z < 0.

Rigor: Arb (python-flint).  H_t(z) = I(t, z) + tail, I = acb.integral over w in [-W, W] (z may be a ball:
the enclosure then covers every z in it), tail bounded analytically (tail_bound).  The winding number:
the box boundary is cut into segments S_k; for each, the enclosure B_k of H(S_k) must exclude 0 with
angular half-width < pi/2, so the continuous argument along S_k stays in one sector and
Delta_k = arg g(end) - arg g(start) is the principal difference.  W = sum Delta_k / 2 pi.
"""
import sys, json, math, time
from flint import arb, acb, ctx

ctx.prec = 160
PI = arb.pi()
K = ((10 - 2 * arb(5).sqrt()).sqrt() - 2) / (arb(5).sqrt() - 1)
A = [None, arb(1), K, -K, arb(-1)]           # a_r for r = 1..4 (a_5 = 0)
I_ = acb(0, 1)


def f(s):
    tot = acb(0)
    for r in range(1, 5):
        tot += A[r] * s.zeta(acb(arb(r) / 5))
    return tot * acb(5) ** (-s)


class DH:
    def __init__(self, X0):
        self.norm = (PI * arb(X0) / 8).exp()        # a positive constant: keeps |xi| = O(1) near x = X0

    def xi(self, s):
        return (acb(5) / PI) ** ((s + 1) / 2) * ((s + 1) / 2).gamma() * f(s) * self.norm

    def H(self, t, z, W=20.0):
        t = arb(t)
        g = lambda w, analytic: (-(w * w) / (4 * t)).exp() * self.xi((1 + w + I_ * z) / 2)
        val = acb.integral(g, -W, W, rel_tol=arb(2) ** -60)
        return val / (4 * PI * t).sqrt() + self.tail_bound(t, z, W)

    def tail_bound(self, t, z, W):
        """|w| > W: sigma' = max(Re s, 1 - Re s) >= (W - |Im z| - 1)/2 >= 2, and by xi(s) = xi(1 - s),
        |xi(s)| <= norm (5/pi)^{(s'+1)/2} Gamma((s'+1)/2) zeta(s') (1 + k) (|Gamma(x+iy)| <= Gamma(x),
        |f| <= zeta(sigma)(1 + k)).  With x = (sigma'+1)/2 <= (|w| + B)/4, B = |Im z| + 3, and
        log Gamma(x) <= x log x, the integrand is <= norm (1+k) zeta(2) exp(-w^2/(4t) + c(w)), where
        c(w) = x log(5 x/pi) is increasing and convex; we bound c(w) <= c(W) + c'(W)... simply by
        c(w) <= w^2/(8t) for w >= W (checked at W; the gap w^2/(8t) - c(w) is increasing there), so the
        tail is <= 2 norm (1+k) zeta(2) int_W^inf e^{-w^2/(8t)} dw <= 2 norm (1+k) zeta(2) (4t/W) e^{-W^2/(8t)},
        divided by (4 pi t)^{1/2} because H normalises the integral."""
        t = arb(t)
        B = abs(z.imag).upper() + 3
        x = (arb(W) + B) / 4
        cW = x * (5 * x / PI).log()
        assert cW < arb(W) ** 2 / (8 * t)
        # gap increasing: d/dw [w^2/(8t)] = w/(4t) >= c'(w) = (log(5x/pi) + 1)/4 for w >= W
        assert arb(W) / (4 * t) > ((5 * x / PI).log() + 1) / 4
        tb = 2 * self.norm * (1 + K) * (PI ** 2 / 6) * (4 * t / W) * (-(arb(W) ** 2) / (8 * t)).exp() / (4 * PI * t).sqrt()
        return acb(arb(0, tb.upper()), arb(0, tb.upper()))     # the ball |Re|, |Im| <= tb


def winding(dh, t, c, h, nseg=8):
    """Winding number of H_t around the square |Re(z-c)|, |Im(z-c)| <= h (counter-clockwise)."""
    corners = [c + acb(h, -h), c + acb(h, h), c + acb(-h, h), c + acb(-h, -h)]
    pts = []
    for i in range(4):
        a, b = corners[i], corners[(i + 1) % 4]
        for j in range(nseg):
            pts.append((a + (b - a) * j / nseg, a + (b - a) * (j + 1) / nseg))
    total = 0.0; worst = 0.0; minabs = 9e9
    vals = {}
    def pv(z):
        key = (float(z.real.mid()), float(z.imag.mid()))
        if key not in vals: vals[key] = dh.H(t, z)
        return vals[key]
    for (za, zb) in pts:
        seg = acb(za.real.union(zb.real), za.imag.union(zb.imag))
        Bk = dh.H(t, seg)
        r = float(Bk.rad()); m = abs(complex(Bk.mid()))
        if not (m > 1.5 * r):
            return dict(ok=False, reason='segment enclosure too wide or contains 0', seg=[str(za), str(zb)], mid=m, rad=r)
        half = math.asin(min(1.0, r / m)) if m > 0 else math.pi
        worst = max(worst, half); minabs = min(minabs, m - r)
        ga, gb = complex(pv(za).mid()), complex(pv(zb).mid())
        d = math.atan2((gb / ga).imag, (gb / ga).real)
        total += d
    w = total / (2 * math.pi)
    return dict(ok=abs(w - round(w)) < 0.25, winding=round(w), winding_raw=w, worst_half_angle=worst,
                min_abs_lower=minabs, segments=len(pts))


if __name__ == '__main__':
    t1 = float(sys.argv[1]); cx = float(sys.argv[2]); cy = float(sys.argv[3]); h = float(sys.argv[4])
    nseg = int(sys.argv[5]) if len(sys.argv) > 5 else 8
    T = time.time()
    dh = DH(round(cx))
    c = acb(cx, cy)
    out = dict(t=t1, center=[cx, cy], half_width=h, box_Im_range=[cy - h, cy + h])
    assert cy + h < 0, 'box must lie strictly below the real axis'
    out.update(winding(dh, t1, c, h, nseg))
    out['secs'] = round(time.time() - T, 1)
    print(json.dumps(out))
