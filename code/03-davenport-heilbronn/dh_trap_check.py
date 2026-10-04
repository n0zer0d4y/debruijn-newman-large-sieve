"""Upper bound, part (b): cross-check of the trapezoid H_t (dh_rect.H_trap, rigorous error ball) against the adaptive
acb.integral H_t (dh_rect.H_direct) at a few points, with timings, at two precisions.  The two enclosures must
overlap; the trapezoid is the fast route for the direct Chebyshev nodes.
"""
import time, math
import dh_rect as d
from flint import arb, acb, ctx

t = 0.6
for prec in (80, 112):
    ctx.prec = prec
    for (x, y) in ((279.5, 0.4), (100.3, 0.4), (485.0, 0.4), (20.0, 1.0), (0.0, 0.9)):
        dh = d.cert_dh.DH(round(x)); W = d.Wof(x, t)
        a, h = 3.0, 0.6
        T = time.time()
        Ms = d.H_sup(dh, t, arb(x - a).union(arb(x + a)), arb(y), W)
        t_sup = time.time() - T
        T = time.time(); v1 = d.H_trap(dh, t, acb(x, y), W, h, a, Ms); t1 = time.time() - T
        T = time.time(); v2 = d.H_direct(dh, t, acb(x, y), W); t2 = time.time() - T
        diff = v1 - v2
        print('prec %d x=%g y=%g  trap %s (%.2fs, sup %.2fs, M %.3g)  integral %s (%.2fs)  overlap %s  |diff| %s' % (
            prec, x, y, v1.str(6), t1, t_sup, float(Ms.upper()), v2.str(6), t2, v1.overlaps(v2), diff.abs_upper().str(3)), flush=True)
