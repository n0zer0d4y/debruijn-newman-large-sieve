# (numeric-experiment) Davenport-Heilbronn: see code/03-davenport-heilbronn/notes.md.
import mpmath as mp
mp.mp.dps = 25
k = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
A = [0, 1, k, -k, -1]           # a_n for n mod 5 (real Davenport-Heilbronn coefficients)
chi = [0, 1, 1j, -1j, -1]
def L(s, c):
    return sum(c[r] * mp.zeta(s, mp.mpf(r)/5) for r in range(1, 5)) * mp.mpf(5)**(-s)
def f(s):
    return L(s, A)
# check a_n = Re((1 - i k) chi(n)) and f = ((1-ik)/2) L(chi) + ((1+ik)/2) L(conj chi)
chib = [0, 1, -1j, 1j, -1]
s0 = mp.mpc(0.3, 7)
print('combination check', abs(f(s0) - ((1-1j*k)/2*L(s0, chi) + (1+1j*k)/2*L(s0, chib))))
# functional equation: xi(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) f(s), xi(s) = xi(1-s)
def xi(s):
    return (mp.mpf(5)/mp.pi)**((s+1)/2) * mp.gamma((s+1)/2) * f(s)
for s in (mp.mpc(0.3, 7), mp.mpc(0.7, 20)):
    print('FE residual', abs(xi(s) - xi(1-s)) / abs(xi(s)))
# the known off-line zero
rho = mp.findroot(f, mp.mpc(0.8085, 85.6993))
print('zero rho =', rho, '|f| =', abs(f(rho)))
# heat flow: Phi(u) = e^{3u} theta(e^{4u}), theta(y) = sum a_n n e^{-pi n^2 y/5}; H_t(z) = 2 int_0^inf Phi e^{t u^2} cos(z u) du
def theta(y):
    s = mp.mpf(0); n = 1
    while True:
        term = A[n % 5] * n * mp.e**(-mp.pi * n * n * y / 5)
        s += term
        if n > 5 and abs(mp.e**(-mp.pi * n * n * y / 5)) * n < mp.mpf(10)**-30: break
        n += 1
    return s
def H(t, z):
    g = lambda u: mp.e**(3*u) * theta(mp.e**(4*u)) * mp.e**(t*u*u) * mp.cos(z*u)
    return 2 * mp.quad(g, mp.linspace(0, 2.2, 60))
z0 = 2*(rho - mp.mpf(1)/2) / 1j
print('z of rho:', z0)
# H_0 should vanish there, relative to its size
print('|H_0(z0)| =', abs(H(0, z0)), ' |H_0(z0 + 0.5)| =', abs(H(0, z0 + 0.5)))
