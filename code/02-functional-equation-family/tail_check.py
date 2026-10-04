"""Sanity check of Fam.tail_bound TAIL = 2 (x7_cert.py) against TAIL = 1 and a direct numerical tail integral.
Usage: python tail_check.py t W cx cy"""
import sys
sys.path.insert(0, __file__.rsplit('\\', 1)[0])
from flint import arb, acb
import mpmath as mp
import x7_cert
from x7_cert import family, Fam, PI, I_

t, W, cx, cy = map(float, sys.argv[1:5])
a, b, chk = family(arb(704) / 1000)
dh = Fam(round(cx), a, W)
dh.REFLECT = True
z = acb(cx, cy)
out = {}
for mode in (1, 2):
    dh.TAIL = mode
    try:
        out[mode] = float(dh.tail_bound(arb(t), z, W).real.rad())
    except AssertionError as e:
        out[mode] = 'assert failed'
print('tail bound TAIL=1:', out[1], '  TAIL=2:', out[2])
# direct numerical tail: int_{|w| > W} e^{-w^2/4t} |xi((1 + w + i z)/2)| dw / sqrt(4 pi t), mpmath quad on [W, W + 60]
mp.mp.dps = 20
def integrand(w):
    s = acb((1 + w - cy) / 2, cx / 2)
    v = dh.xi(s)
    return mp.e ** (-(w * w) / (4 * t)) * float(v.abs_upper())
tot = 0
for w0 in range(int(W), int(W) + 60, 4):
    tot += mp.quad(lambda w: integrand(float(w)), [w0, w0 + 4])
for w0 in range(-int(W) - 60, -int(W), 4):
    tot += mp.quad(lambda w: integrand(float(w)), [w0, w0 + 4])
print('numerical tail (both sides, W..W+60):', mp.nstr(tot / mp.sqrt(4 * mp.pi * t), 6))
