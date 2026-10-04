"""Certificate by Rouche (replaces the segment-ball winding of cert_dh.py, whose z-ball enclosures of
the oscillatory w-integral were too wide: radius 3.9 against a midpoint of 0.29).

Given t and a centre c (Im c < 0), with all quantities rigorous (Arb):
  h0 = H_t(c)                          point evaluation (acb.integral at a point z, plus tail);
  M  >= sup_{|z-c| <= R} |H_t(z)|      modulus majorant: (4 pi t)^{-1/2} int e^{-w^2/(4t)} sup|xi| dw, with the
                                       w-axis cut into cells, sup|xi| over each (s-ball) by Arb, the Gaussian mass
                                       of each cell by erf, plus the analytic tail of cert_dh.DH.tail_bound;
  d  = (H(c + delta) - H(c - delta))/(2 delta), with |H'(c) - d| <= M3 delta^2/6, M3 <= 6M/R^3 (Cauchy);
  M2 <= 2M/(R - r)^2 on |z - c| <= r.
On |z - c| = r, with p(z) = h0 + d (z - c):
  |H - p| <= |H'(c) - d| r + M2 r^2/2  <  |d| r - |h0| <= |p|
gives (Rouche) exactly one zero of H_t in |z - c| < r.  If r < |Im c|, that zero is non-real, so
H_t is not in the Laguerre-Polya class and Lambda_DH > t.
"""
import sys, json, math, time
from flint import arb, acb, ctx
from cert_dh import DH, K, PI, I_

ctx.prec = 160


def point(dh, t, z):
    return dh.H(t, z)


def majorant(dh, t, c, R, W=20.0, cells=400):
    """sup_{|z - c| <= R} |H_t(z)|."""
    t = arb(t)
    # box [c - R, c + R]^2 built in Arb from the exact double centre, with R padded by 1e-9 R so that the
    # float centre and float R cannot leave a sliver of the disc |z - c| <= R uncovered
    Rp = arb(R) * (1 + arb(10) ** -9)
    cr, ci = arb(float(c.real.mid())), arb(float(c.imag.mid()))
    zre = (cr - Rp).union(cr + Rp)
    zim = (ci - Rp).union(ci + Rp)
    tot = arb(0)
    edges = [-W + 2 * W * j / cells for j in range(cells + 1)]
    for j in range(cells):
        a, b = edges[j], edges[j + 1]
        wb = arb(a).union(arb(b))
        # s = (1 + w + i z)/2 = (1 + w - Im z)/2 + i Re z / 2
        s = acb((1 + wb - zim) / 2, zre / 2)
        sup = dh.xi(s).abs_upper()
        mass = ((arb(b) / (2 * t.sqrt())).erf() - (arb(a) / (2 * t.sqrt())).erf()) / 2   # int e^{-w^2/4t}/sqrt(4 pi t)
        tot += sup * mass
    tail = dh.tail_bound(t, acb(zre, zim), W)
    return tot + tail.real.upper() * 2


def certify(t, cx, cy, R=0.05, r=2e-4, delta=1e-4, newton=6):
    T = time.time()
    dh = DH(round(cx))
    c = acb(cx, cy)
    # Newton on rigorous point values (the iteration itself need not be rigorous)
    for _ in range(newton):
        h = point(dh, t, c)
        d = (point(dh, t, c + delta) - point(dh, t, c - delta)) / (2 * delta)
        step = h / d
        c = acb(float((c - step).real.mid()), float((c - step).imag.mid()))
        if abs(complex(step.mid())) < 1e-13: break
    h0 = point(dh, t, c)
    d = (point(dh, t, c + delta) - point(dh, t, c - delta)) / (2 * delta)
    M = majorant(dh, t, c, R)
    M3 = 6 * M / arb(R) ** 3
    M2 = 2 * M / (arb(R) - r) ** 2
    err_d = M3 * arb(delta) ** 2 / 6
    lhs = err_d * r + M2 * arb(r) ** 2 / 2          # bound for |H - p| on the circle
    rhs = d.abs_lower() * r - h0.abs_upper()        # lower bound for |p| on the circle
    ok = bool(lhs < rhs) and (r < -float(c.imag.mid()))
    return dict(t=t, center=[float(c.real.mid()), float(c.imag.mid())], r=r, R=R,
                H_center_abs_upper=float(h0.abs_upper()), d_abs_lower=float(d.abs_lower()),
                M=float(M.upper()), lhs=float(lhs.upper()), rhs=float(rhs.lower()), ok=ok,
                conclusion=('H_t has exactly one zero in |z - c| < r, non-real; hence Lambda_DH > %g' % t) if ok else None,
                secs=round(time.time() - T, 1))


if __name__ == '__main__':
    t = float(sys.argv[1]); cx = float(sys.argv[2]); cy = float(sys.argv[3])
    kw = {}
    for a in sys.argv[4:]:
        k_, v = a.split('='); kw[k_] = float(v)
    print(json.dumps(certify(t, cx, cy, **kw)))
