"""(numeric-experiment) heat-flow track of an off-line zero of f_theta (q = 7, odd Theorem A' family,
x7_scan.make_f), in the variable z of H_t(z) = E_v[xi((1 + 2 sqrt(t) v + iz)/2)] (Gauss-Hermite, 28 nodes),
xi(s) = (7/pi)^{(s+1)/2} Gamma((s+1)/2) f_theta(s) (times a constant). The same scheme as
../03-davenport-heilbronn/dh_track_zero240.py.  It prints (t, x, y) until the zero reaches the real axis.
Usage: python x7_track.py theta sigma0 t0 [tmax] [tstart zx zy [dt]]
"""
import sys, time
import numpy as np
import mpmath as mp
sys.path.insert(0, __file__.rsplit('\\', 1)[0])
from x7_scan import make_f

mp.mp.dps = 18
theta, s0, g0 = float(sys.argv[1]), sys.argv[2], sys.argv[3]
TMAX = float(sys.argv[4]) if len(sys.argv) > 4 else 0.8
f, a = make_f(theta)
rho = mp.findroot(f, mp.mpc(s0, g0))
print('rho =', rho, '|f(rho)| =', abs(f(rho)), flush=True)
X0 = 2 * rho.imag


def xi_n(s):
    return (mp.mpf(7) / mp.pi) ** ((s + 1) / 2) * mp.gamma((s + 1) / 2) * f(s) * mp.e ** (mp.pi * X0 / 8)


import os
v, wt = np.polynomial.hermite.hermgauss(int(os.environ.get('GH_NODES', '28')))   # GH_NODES overrides 28


def H(t, z):
    st = mp.sqrt(t)
    return sum(mp.mpf(float(wj)) * xi_n((1 + 2 * st * mp.mpf(float(vj)) + 1j * z) / 2) for vj, wj in zip(v, wt)) / mp.sqrt(mp.pi)


z = 2 * (rho - mp.mpf(1) / 2) / 1j
t = 0.0; dt = 0.02
if len(sys.argv) > 7:                                   # restart: tstart zx zy [dt]
    t = float(sys.argv[5]); z = mp.mpc(sys.argv[6], sys.argv[7])
    dt = float(sys.argv[8]) if len(sys.argv) > 8 else 0.1
    z = mp.findroot(lambda zz: H(t, zz), z, tol=1e-9)
print('z0 =', z, 't0 =', t, flush=True)
T1 = time.time()
while t < TMAX and time.time() - T1 < 250:
    t2 = t + dt
    try:
        z2 = mp.findroot(lambda zz: H(t2, zz), z, tol=1e-9)
    except Exception as e:
        if dt > 0.0011:
            dt /= 4; continue
        print('fail at', t2, e); break
    if abs(z2.imag) < 0.01:
        print('reached axis between t=%.4f and %.4f: z=%s' % (t, t2, complex(z2))); break
    z, t = z2, t2
    print('t=%.4f  x=%.5f  y=%.5f  (%.0fs)' % (t, float(z.real), float(z.imag), time.time() - T1), flush=True)
