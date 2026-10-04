"""Upper bound, part (b): timing and sanity of the direct Arb evaluation of H_t (cert_dh.DH.H) for x < 280,
and the ratio G = H_t/B_t (B_t = exp(lMt(s+, t)) of cert_final.Char conventions, A = 1, q = 5),
at a few (prec, W, rel_tol) settings.  The tail bound of cert_dh needs W^2/(8t) > pi X0/8 + margin.
"""
import sys, time, math
import os as _os
sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'shared-dirichlet-pipeline'))
from flint import arb, acb, ctx
import cert_dh


def lMt(s, t):
    PI = arb.pi(); lq = (arb(5) / PI).log(); sa = s + 1
    am = lq / 2 + (sa / 2).log() / 2 - 1 / (2 * sa)
    w = sa / 2
    lM0 = w * lq + (2 * PI).sqrt().log() + (w - arb(1) / 2) * w.log() - w
    return (t / 4) * am ** 2 + lM0


def H(dh, t, z, W, rt):
    tt = arb(t)
    g = lambda w, analytic: (-(w * w) / (4 * tt)).exp() * dh.xi((1 + w + acb(0, 1) * z) / 2)
    val = acb.integral(g, -W, W, rel_tol=arb(2) ** -rt)
    return val / (4 * arb.pi() * tt).sqrt() + dh.tail_bound(tt, z, W)


t = 0.6
for (x, y) in ((279.5, 0.4), (100.3, 0.4), (20.0, 1.0)):
    X0 = round(x)
    W = max(20.0, math.ceil(math.sqrt(8 * t * (math.pi * X0 / 8 + 45))))
    for prec, rt in ((160, 60), (96, 30), (80, 24)):
        ctx.prec = prec
        dh = cert_dh.DH(X0)
        z = acb(x, y)
        T = time.time(); v = H(dh, t, z, W, rt); dt = time.time() - T
        sp = (1 + arb(y) - acb(0, 1) * arb(x)) / 2
        G = v / (dh.norm * lMt(sp, arb(t)).exp())
        print('x=%g y=%g W=%g prec=%d rt=%d  %.2fs  H=%s  G=%s' % (x, y, W, prec, rt, dt, v.str(5), G.str(5)), flush=True)
