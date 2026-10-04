"""first-order pair-collision law (numeric check).  Near a close pair z1 < z2 of real zeros of H^chi
(H(z) = xi((1+iz)/2), chi = chi_{-7}), model H^chi(z) ~ (A/2)[(z-zm)^2 - g^2/4], H^b(z) ~ B (the DFT-self-dual
partner b of x7_cert.family).  Then H^chi + tan(theta) H^b has a double real zero at
    tan(theta_c) = g^2 A / (8 B)            (needs A B > 0 for theta > 0),
and for theta > theta_c an off-line pair with Im z = y0, y0^2 = 2 (B/A)(tan theta - tan theta_c); under the heat flow
dH/dt = -H'' the quadratic model collides at t = y0^2/2.
Usage: python pair_law.py gamma1 gamma2   (ordinates of the pair of L(s, chi_{-7}))
"""
import sys
import mpmath as mp
sys.path.insert(0, __file__.rsplit('\\', 1)[0])
from x7_scan import basis

mp.mp.dps = 30
g1, g2 = mp.mpf(sys.argv[1]), mp.mpf(sys.argv[2])
chi = [0, 1, 1, -1, 1, -1, -1]
_, bvec = basis()                                   # float SVD partner (sign b(1) < 0), unit norm as in scans
b = [mp.mpf(float(x)) for x in bvec]
nc = mp.sqrt(sum(mp.mpf(c) ** 2 for c in chi)); nb = mp.sqrt(sum(x ** 2 for x in b))
chi_u = [mp.mpf(c) / nc for c in chi]; b_u = [x / nb for x in b]


def xi(a, s):
    L = sum(a[r] * mp.zeta(s, mp.mpf(r) / 7) for r in range(1, 7)) * mp.mpf(7) ** (-s)
    return (mp.mpf(7) / mp.pi) ** ((s + 1) / 2) * mp.gamma((s + 1) / 2) * L


def H(a, z):
    return xi(a, (1 + 1j * z) / 2)


z1, z2 = 2 * g1, 2 * g2
zm = (z1 + z2) / 2; gz = z2 - z1
A = mp.re(mp.diff(lambda z: H(chi_u, z), zm, 2))
B = mp.re(H(b_u, zm))
print('check H^chi(z1), H^chi(z2) relative:', mp.nstr(abs(H(chi_u, z1)) / abs(A), 3), mp.nstr(abs(H(chi_u, z2)) / abs(A), 3))
print('imag parts (should be ~0):', mp.nstr(mp.im(H(chi_u, zm)) / abs(A), 3), mp.nstr(mp.im(H(b_u, zm)) / abs(A), 3))
print('A = H^chi\'\'(zm) =', mp.nstr(A, 8), '  B = H^b(zm) =', mp.nstr(B, 8), '  B/A =', mp.nstr(B / A, 8))
tc = mp.atan(gz ** 2 * A / (8 * B))
print('z-gap g =', mp.nstr(gz, 8), '  predicted theta_c =', mp.nstr(tc, 8))
for th in ['0.02', '0.03', '0.05', '0.1']:
    th = mp.mpf(th)
    y02 = 2 * (B / A) * (mp.tan(th) - mp.tan(tc))
    print('theta %s: predicted y0 = %s, t_c(quadratic) = %s' % (th, mp.nstr(mp.sqrt(y02), 5) if y02 > 0 else '-', mp.nstr(y02 / 2, 5)))
