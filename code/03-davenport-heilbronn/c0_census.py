"""(numeric) Census of c0 = (1/2q) sum_{n mod 2q} chi_D(n) e^{-i pi n^2/q} over fundamental discriminants 0 < |D| <= 400
(papers/supporting/riemann-siegel-identity/proof-rs-general.md, Lemma G3 and section 9).
chi_D is the Kronecker symbol (D/.), q = |D|. Lemma G3 predicts c0 = 0 for odd chi (D < 0) and for odd q.
Usage: python c0_census.py
"""
import cmath, math


def kron(D, n):
    """Kronecker symbol (D/n) for n >= 0."""
    if n == 0:
        return 1 if abs(D) == 1 else 0
    res = 1
    while n % 2 == 0:
        n //= 2
        if D % 2 == 0:
            return 0
        res *= 1 if D % 8 in (1, 7) else -1
    a, m = D % n, n
    while a:                                    # Jacobi symbol (a/m), m odd
        while a % 2 == 0:
            a //= 2
            if m % 8 in (3, 5):
                res = -res
        a, m = m, a
        if a % 4 == 3 and m % 4 == 3:
            res = -res
        a %= m
    return res if m == 1 else 0


def fundamental(D):
    if D in (0, 1):
        return False
    if D % 4 == 1:
        return all(D % (p * p) for p in range(2, int(abs(D) ** 0.5) + 1))
    if D % 4 == 0:
        m = D // 4
        return m % 4 in (2, 3) and all(m % (p * p) for p in range(2, int(abs(m) ** 0.5) + 2))
    return False


Ds = [D for D in range(-400, 401) if fundamental(D)]
nz = []
for D in Ds:
    q = abs(D)
    c0 = sum(kron(D, n) * cmath.exp(-1j * math.pi * n * n / q) for n in range(2 * q)) / (2 * q)
    if abs(c0) > 1e-9:
        nz.append((D, round(abs(c0), 6), round(cmath.phase(c0) / math.pi, 6)))
print('fundamental discriminants with 0 < |D| <= 400:', len(Ds))
print('c0 != 0 for (D, |c0|, arg c0 / pi):', nz)
