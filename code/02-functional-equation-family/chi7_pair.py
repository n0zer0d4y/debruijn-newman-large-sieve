"""Check for claim 02-#4 (numeric): zeros of L(s, chi_{-7}) on the critical line near t = 51.14, via sign changes of the
Hardy-type function Z(t) = e^{i theta(t)} L(1/2+it) (real for the odd real character), mpmath dps 25.
Usage: python chi7_pair.py T1 T2 step"""
import sys
import mpmath as mp
mp.mp.dps = 25
chi = [0, 1, 1, -1, 1, -1, -1]
def L(s):
    return sum(chi[r] * mp.zeta(s, mp.mpf(r) / 7) for r in range(1, 7)) * mp.mpf(7) ** (-s)
def xi(t):                       # completed, real on the line for this odd character (root number 1)
    s = mp.mpf(1) / 2 + 1j * t
    v = (mp.mpf(7) / mp.pi) ** ((s + 1) / 2) * mp.gamma((s + 1) / 2) * L(s)
    return v
T1, T2, h = map(float, sys.argv[1:4])
prev = None; t = T1; zs = []
while t <= T2:
    v = xi(t)
    re = mp.re(v * mp.exp(mp.pi * t / 4))
    if prev is not None and mp.sign(re) != mp.sign(prev[1]):
        z = mp.findroot(lambda u: mp.re(xi(u) * mp.exp(mp.pi * u / 4)), (prev[0] + t) / 2)
        zs.append(z)
    prev = (t, re); t += h
print('imag-part check at T1:', mp.nstr(abs(mp.im(xi(T1))) / abs(xi(T1)), 3))
print('zeros:', [mp.nstr(z, 12) for z in zs])
print('gaps:', [mp.nstr(zs[i + 1] - zs[i], 6) for i in range(len(zs) - 1)])
