"""Upper bound, step 1: does Theorem A (papers/supporting/riemann-siegel-identity/proof-rs-general.md) hold with the
Davenport-Heilbronn coefficients in place of a real character?

Theorem A (for real primitive chi mod q, parity a):
  L(s) = int_{c SW} w^{-s} Kc(w) dw + (-1)^a X(s) int_{c SE} w^{s-1} conj(Kc(-conj w)) dw,
  Kc = K + (c0/2) e^{i pi w^2/q},  K(w) = e^{i pi w^2/q} sum_{j<2q} c_j e^{pi i j w/q}/(e^{2 pi i w} - 1),
  e^{i pi n^2/q} sum_j c_j e^{pi i j n/q} = a_n (n mod 2q),  c0 = (1/2q) sum_{n mod 2q} a_n e^{-i pi n^2/q},
  X(s) = (pi/q)^{s-1/2} Gamma((1 - s + a)/2)/Gamma((s + a)/2);
  l_c = {c + u w}, w = e^{i pi/4}, traversed u: +inf -> -inf ('c SW'); {c + u conj(w)}, u: -inf -> +inf ('c SE').
Here a_n = (1, k, -k, -1, 0) for n = 1..5 mod 5 (real, period 5), parity a = 1 (Gamma((s+1)/2)).
Prints the relative residual at several s (on and off the line, small and moderate |Im s|),
for the symmetrised kernel and for K alone.  mpmath, a single process, about a minute.
"""
import mpmath as mp

mp.mp.dps = 30
q, a = 5, 1
k = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
tab = [mp.mpf(0), mp.mpf(1), k, -k, mp.mpf(-1)]
an = lambda n: tab[n % q]


def f(s):
    return sum(tab[r] * mp.zeta(s, mp.mpf(r) / q) for r in range(1, q)) * mp.mpf(q) ** (-s)


M = 2 * q
Amat = mp.matrix(M, M); b = mp.matrix(M, 1)
for n in range(M):
    for j in range(M):
        Amat[n, j] = mp.exp(1j * mp.pi * n * n / q) * mp.exp(mp.pi * 1j * j * n / q)
    b[n] = an(n)
cj = mp.lu_solve(Amat, b)
c0 = sum(an(n) * mp.exp(-1j * mp.pi * n * n / q) for n in range(M)) / M
print('c0 =', mp.nstr(c0, 12))


def K(w):
    return mp.exp(1j * mp.pi * w * w / q) * sum(cj[j] * mp.exp(mp.pi * 1j * j * w / q) for j in range(M)) / (mp.exp(2j * mp.pi * w) - 1)


def Kc(w, sym=True):
    return K(w) + ((c0 / 2) * mp.exp(1j * mp.pi * w * w / q) if sym else 0)


om = mp.exp(1j * mp.pi / 4)


def rhs(s, c=mp.mpf('0.5'), sym=True, U=9):
    # c SW: w = c + u om, u from +U to -U  ->  int = -int_{-U}^{U} ... om du
    I1 = -mp.quad(lambda u: (c + u * om) ** (-s) * Kc(c + u * om, sym) * om, [-U, -2, 0, 2, U])
    omb = mp.conj(om)
    I2 = mp.quad(lambda u: (c + u * omb) ** (s - 1) * mp.conj(Kc(-mp.conj(c + u * omb), sym)) * omb, [-U, -2, 0, 2, U])
    X = (mp.pi / q) ** (s - mp.mpf(1) / 2) * mp.gamma((1 - s + a) / 2) / mp.gamma((s + a) / 2)
    return I1 + (-1) ** a * X * I2


for s in (mp.mpc(0.5, 3), mp.mpc(0.3, 7), mp.mpc(1.7, -2), mp.mpc(0.8085, 12)):
    lhs = f(s)
    for sym in (True, False):
        r = rhs(s, sym=sym)
        print('s =', mp.nstr(s, 6), 'sym' if sym else 'K  ', 'rel residual', mp.nstr(abs(r - lhs) / abs(lhs), 4))
