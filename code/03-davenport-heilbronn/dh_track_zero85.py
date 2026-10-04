# (numeric-experiment) Davenport-Heilbronn: see code/03-davenport-heilbronn/notes.md.
import mpmath as mp, time
mp.mp.dps = 20
k = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
A = [0, 1, k, -k, -1]
def f(s):
    return sum(A[r] * mp.zeta(s, mp.mpf(r)/5) for r in range(1, 5)) * mp.mpf(5)**(-s)
X0 = mp.mpf(171)
def xi_n(s):
    # xi normalised by a fixed factor so values are O(1) near x = 171: divide by |Gamma((1 + i X0/2)/2 ...)| scale
    return (mp.mpf(5)/mp.pi)**((s+1)/2) * mp.gamma((s+1)/2) * f(s) * mp.e**(mp.pi * X0 / 8)
def H(t, z):
    if t == 0: return xi_n((1 + 1j*z)/2)
    W = 9 * mp.sqrt(t)
    g = lambda w: mp.e**(-w*w/(4*t)) * xi_n((1 + w + 1j*z)/2)
    return mp.quad(g, mp.linspace(-W, W, 13)) / mp.sqrt(4*mp.pi*t)
z = mp.mpc('171.3986969707551843', '-0.6170343649132747')
print('H_0(z0) ~', abs(H(0, z)), ' scale |H_0(z0+1)| ~', abs(H(0, z + 1)))
T1 = time.time(); out = []
for t in [0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18]:
    try:
        z = mp.findroot(lambda zz: H(t, zz), z, tol=1e-12)
        out.append((t, float(z.real), float(z.imag)))
        print('t=%.2f  zero x=%.6f  y=%.6f  (%.0fs)' % (t, float(z.real), float(z.imag), time.time()-T1), flush=True)
    except Exception as e:
        print('t=%.2f failed: %s' % (t, e)); break
    if time.time() - T1 > 240: break
