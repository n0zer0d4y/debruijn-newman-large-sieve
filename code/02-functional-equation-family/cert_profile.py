"""Profile the pieces of x7_cert.certify_hp at large t: one rigorous point value of H_t (acb.integral) and the majorant.
Usage: python cert_profile.py t cx cy W rt cells R"""
import sys, time
sys.path.insert(0, __file__.rsplit('\\', 1)[0])
from flint import arb, acb
from x7_cert import family, Fam, PI, I_
import cert_dh2

t, cx, cy, W, rt, cells, R = sys.argv[1:8]
t = float(t); cx = float(cx); cy = float(cy); W = float(W); rt = int(rt); cells = int(cells); R = float(R)
a, b, chk = family(arb(704) / 1000)
dh = Fam(round(cx), a, W); dh.TAIL = 2; dh.REFLECT = True
tt = arb(t); z = acb(cx, cy)
T = time.time()
g = lambda w, analytic: (-(w * w) / (4 * tt)).exp() * dh.xi((1 + w + I_ * z) / 2)
val = acb.integral(g, -dh.W, dh.W, rel_tol=arb(2) ** -rt)
h = val / (4 * PI * tt).sqrt() + dh.tail_bound(tt, z, dh.W)
print('H point value: %s  (%.1fs)' % (h.str(radius=True)[:80], time.time() - T), flush=True)
T = time.time()
M = cert_dh2.majorant(dh, t, z, R, W=W, cells=cells)
print('majorant M <= %s  (%.1fs, %d cells)' % (M.str()[:40], time.time() - T, cells), flush=True)
