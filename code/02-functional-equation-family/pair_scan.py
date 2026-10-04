"""pair-collision angles for every consecutive zero pair of L(s, chi_{-7}) up to height T (numeric).
For each pair (g1, g2): z-gap g = 2(g2 - g1), A = H^chi''(zm), B = H^b(zm) (unit-normalised chi_{-7} and its
DFT-self-dual partner b, as pair_law.py), and the local-quadratic collision angle tan(theta_c) = g^2 A/(8B).
A positive theta_c means the pair leaves the line for theta slightly above theta_c; a negative one, below.
Prints the pairs with the smallest |theta_c| and summary statistics.
Usage: python pair_scan.py T [step]
"""
import sys
import mpmath as mp
sys.path.insert(0, __file__.rsplit('\\', 1)[0])
from x7_scan import basis

mp.mp.dps = 20
T = float(sys.argv[1]); h = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
chi = [0, 1, 1, -1, 1, -1, -1]
_, bvec = basis()
b = [mp.mpf(float(x)) for x in bvec]
nc = mp.sqrt(sum(mp.mpf(c) ** 2 for c in chi)); nb = mp.sqrt(sum(x ** 2 for x in b))
chi_u = [mp.mpf(c) / nc for c in chi]; b_u = [x / nb for x in b]


def xi(a, s):
    L = sum(a[r] * mp.zeta(s, mp.mpf(r) / 7) for r in range(1, 7)) * mp.mpf(7) ** (-s)
    return (mp.mpf(7) / mp.pi) ** ((s + 1) / 2) * mp.gamma((s + 1) / 2) * L


def Zn(a, t):                       # xi on the line, scaled by e^{pi t/4} to remove the Gamma decay (real)
    return mp.re(xi(a, mp.mpf(1) / 2 + 1j * t) * mp.exp(mp.pi * t / 4))


zs = []; t = 1.0; prev = Zn(chi_u, t)
while t < T:
    t2 = t + h; v = Zn(chi_u, t2)
    if mp.sign(v) != mp.sign(prev):
        zs.append(mp.findroot(lambda u: Zn(chi_u, u), (t + t2) / 2))
    t, prev = t2, v
rows = []
for i in range(len(zs) - 1):
    g1, g2 = zs[i], zs[i + 1]
    tm = (g1 + g2) / 2
    # in the z variable (z = 2t), with the same e^{pi t/4} scaling for both functions (ratio-invariant)
    A = mp.diff(lambda u: Zn(chi_u, u), tm, 2) / 4          # d^2/dz^2 = (1/4) d^2/dt^2
    B = Zn(b_u, tm)
    gz = 2 * (g2 - g1)
    tc = mp.atan(gz ** 2 * A / (8 * B))
    rows.append((abs(tc), tc, g1, g2, g2 - g1, B / A))
print('zeros found:', len(zs), ' pairs:', len(rows))
rows.sort(key=lambda r: r[0])
print('smallest |theta_c|:')
for r in rows[:8]:
    print('  theta_c = %s  pair %s / %s  gap %s  B/A %s' % tuple(mp.nstr(x, 6) for x in r[1:]))
check = [r for r in rows if abs(r[2] - 50.983) < 0.01]
print('check pair at 50.983:', [mp.nstr(r[1], 6) for r in check])
neg = [r for r in rows if r[1] < 0]
print('pairs with theta_c < 0 (collide for theta < 0):', len(neg), 'of', len(rows))
if neg:
    m = min(neg, key=lambda r: r[0])
    print('  smallest negative: theta_c = %s at pair %s / %s' % (mp.nstr(m[1], 6), mp.nstr(m[2], 8), mp.nstr(m[3], 8)))
print('B/A quantiles:', [mp.nstr(sorted(r[5] for r in rows)[int(q * (len(rows) - 1))], 4) for q in (0, 0.1, 0.5, 0.9, 1)])
