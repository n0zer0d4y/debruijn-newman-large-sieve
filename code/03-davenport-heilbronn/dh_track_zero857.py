# (numeric-experiment) Davenport-Heilbronn: see code/03-davenport-heilbronn/notes.md.
import mpmath as mp, time, numpy as np
mp.mp.dps = 18
k = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
A = [0, 1, k, -k, -1]
def f(s):
    return sum(A[r] * mp.zeta(s, mp.mpf(r)/5) for r in range(1, 5)) * mp.mpf(5)**(-s)
v, wt = np.polynomial.hermite.hermgauss(28)
T1 = time.time()
for guess in [('0.85610', '857.2958')]:
    rho = mp.findroot(f, mp.mpc(*guess))
    X0 = 2 * rho.imag
    xi_n = lambda s: (mp.mpf(5)/mp.pi)**((s+1)/2) * mp.gamma((s+1)/2) * f(s) * mp.e**(mp.pi * X0 / 8)
    H = lambda t, z: sum(mp.mpf(float(wj)) * xi_n((1 + 2*mp.sqrt(t)*mp.mpf(float(vj)) + 1j*z)/2) for vj, wj in zip(v, wt)) / mp.sqrt(mp.pi)
    z = 2*(rho - mp.mpf(1)/2)/1j
    print('rho', complex(rho), 'z0', complex(z), flush=True)
    t = 0.0; dt = 0.04
    while t < 0.5 and time.time() - T1 < 260:
        t2 = t + dt
        try:
            z2 = mp.findroot(lambda zz: H(t2, zz), z, tol=1e-9)
        except Exception as e:
            if dt > 0.0011: dt /= 4; continue
            print('fail', t2, e); break
        if abs(z2.imag) < 0.01:
            print('axis between t=%.4f and %.4f' % (t, t2)); break
        z, t = z2, t2
        print('t=%.4f y=%.5f (%.0fs)' % (t, float(z.imag), time.time()-T1), flush=True)
        if abs(z.imag) < 0.25: dt = min(dt, 0.005)
        if abs(z.imag) < 0.08: dt = min(dt, 0.001)
