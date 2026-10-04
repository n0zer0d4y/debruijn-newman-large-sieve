"""(numeric) deep zeros of f_theta near theta* (a_theta(1) -> 0) from the two-term balance
a(1) + a(2) 2^{-s} = 0, i.e. sigma0 = log|a(2)/a(1)|/log 2, t_k = (arg(-a(1)/a(2)) + 2 pi k)/log 2 ... refined by
findroot on the full f_theta (x7_scan.make_f).  Prints a(1), sigma0 and the refined zeros with t < TMAX.
Usage: python deep_zeros.py theta [TMAX]
"""
import sys, math
import mpmath as mp
sys.path.insert(0, __file__.rsplit('\\', 1)[0])
from x7_scan import make_f

mp.mp.dps = 30
theta = float(sys.argv[1]); TMAX = float(sys.argv[2]) if len(sys.argv) > 2 else 40
f, a = make_f(theta)
a1, a2 = mp.mpf(float(a[1])), mp.mpf(float(a[2]))
print('a(1) = %s, a(2) = %s, a(3) = %s' % (mp.nstr(a1, 8), mp.nstr(a2, 8), mp.nstr(mp.mpf(float(a[3])), 8)))
r = -a1 / a2                                   # 2^{-s} = r
sig0 = -mp.log(abs(r)) / mp.log(2)
ph = 0 if r > 0 else mp.pi                     # arg(2^{-s}) = -t log 2 = arg r (mod 2 pi)
print('two-term balance sigma0 = %s' % mp.nstr(sig0, 8))
k = 0; found = []
while True:
    t0 = (2 * mp.pi * k - ph) / mp.log(2)
    if t0 > TMAX:
        break
    if t0 > 0.5:
        try:
            z = mp.findroot(f, mp.mpc(sig0, t0))
            found.append(z)
            print('k=%d  start (%s, %s)  zero (%s, %s)  |f| = %s' % (k, mp.nstr(sig0, 6), mp.nstr(t0, 6),
                  mp.nstr(z.real, 10), mp.nstr(z.imag, 10), mp.nstr(abs(f(z)), 3)))
        except Exception as e:
            print('k=%d no convergence: %s' % (k, str(e)[:80]))
    k += 1
