# (numeric-experiment) Davenport-Heilbronn: see code/03-davenport-heilbronn/notes.md.
import mpmath as mp, time
mp.mp.dps = 18
k = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
A = [0, 1, k, -k, -1]
def f(s):
    return sum(A[r] * mp.zeta(s, mp.mpf(r)/5) for r in range(1, 5)) * mp.mpf(5)**(-s)
X0 = mp.mpf(171)
def xi_n(s):
    return (mp.mpf(5)/mp.pi)**((s+1)/2) * mp.gamma((s+1)/2) * f(s) * mp.e**(mp.pi * X0 / 8)
# Gauss-Hermite: int e^{-w^2/(4t)} F(w) dw / sqrt(4 pi t) = (1/sqrt(pi)) sum_j w_j F(2 sqrt(t) v_j)
import numpy as np
v, wt = np.polynomial.hermite.hermgauss(24)
def H(t, z):
    st = mp.sqrt(t)
    return sum(mp.mpf(float(wj)) * xi_n((1 + 2*st*mp.mpf(float(vj)) + 1j*z)/2) for vj, wj in zip(v, wt)) / mp.sqrt(mp.pi)
z = mp.mpc('171.3986969707551843', '-0.6170343649132747')
T1 = time.time()
# accuracy check vs previous adaptive value at t = 0.01: zero x=171.389367 y=-0.598013
z = mp.findroot(lambda zz: H(0.01, zz), z, tol=1e-10)
print('t=0.01 GH zero', complex(z), round(time.time()-T1,1), 's', flush=True)
t = 0.01; dt = 0.01
while t < 0.3 and time.time() - T1 < 250:
    t2 = t + dt
    try:
        z2 = mp.findroot(lambda zz: H(t2, zz), z, tol=1e-10)
    except Exception as e:
        print('fail at', t2, e); break
    if abs(z2.imag) < 0.02:
        print('near axis at t=%.3f: z=%s' % (t2, complex(z2))); z = z2; t = t2; break
    z, t = z2, t2
    print('t=%.3f  x=%.6f  y=%.6f  (%.0fs)' % (t, float(z.real), float(z.imag), time.time()-T1), flush=True)
    if abs(z.imag) < 0.15: dt = 0.002
