"""theta* (a_theta(1) = 0) in exact Arb, and rational theta = p/10^k with a prescribed small a_theta(1) > 0.
a_theta(1) = cos theta + sin theta b(1), so theta* = atan(-1/b(1)) (b(1) < 0) and a'(theta*) = -1/cos(theta*)... computed.
Usage: python theta_star.py eps1 [eps2 ...]   (targets for a_theta(1))"""
import sys
sys.path.insert(0, __file__.rsplit('\\', 1)[0])
from flint import arb
from x7_cert import family
a0, b, chk = family(arb(0))
b1 = b[1]
ts = (-1 / b1).atan()
print('b(1) =', b1.str(25), ' theta* =', ts.str(30))
for e in sys.argv[1:]:
    e = arb(e)
    # a(theta) = cos + sin b1; derivative -sin + cos b1 < 0 near theta*
    d = -ts.sin() + ts.cos() * b1
    th = ts - e / d * -1 if False else ts + e / d            # first-order step toward a = e (d < 0, so theta < theta*)
    k = 20
    p = int((th * 10 ** k).mid().floor().unique_fmpz())
    thr = arb(p) / 10 ** k
    a1 = thr.cos() + thr.sin() * b1
    print('target %s: theta = %d/10^%d, a(1) = %s' % (e.str(3), p, k, a1.str(6)))
