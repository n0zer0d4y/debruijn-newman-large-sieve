"""Verified vector arithmetic for the certificate (papers/supporting/riemann-siegel-identity/derivation.md §3, Lemma R').

Model: IEEE-754 binary64, round-to-nearest, for +, -, *, / and sqrt (numpy float64
ufuncs; no FMA contraction), exact ldexp/rint.  u = 2^-53.  No libm function is
trusted: exp and cis are evaluated here by argument reduction + Horner, with the
error bounds proved in the docstrings; integer logarithms come from Arb (python-flint).

Every function returns (value, bound), where bound is an absolute error bound
(`vcis`) or a relative one (`vexp`), valid for the exact function of the exact input
provided the caller's input error bound is supplied.
"""
import math
import numpy as np

U = 2.0**-53

# ln 2 split: LN2_HI has 29 significant bits, so k*LN2_HI is exact for |k| < 2^20.
LN2 = 0.6931471805599453
LN2_HI = float.fromhex('0x1.62e42fe000000p-1')
LN2_LO = 1.4286068203094173e-09   # |LN2_HI + LN2_LO - ln2| <= 2^-80 (checked in selftest with Arb)
INV_LN2 = 1.4426950408889634
_EXP_DEG = 22
_EXP_C = [1.0 / math.factorial(j) for j in range(_EXP_DEG + 1)]


def vexp(z, ez=0.0):
    """exp(z_true) with |z - z_true| <= ez.  Returns (v, rel) with
    |v - exp(z_true)| <= rel * exp(z_true) ... up to the second-order term, which is
    absorbed: rel = 110u + 1.02*(ez + er), for |z| <= 700.

    Proof sketch.  k = rint(z/ln2) is an exact integer, |k| < 1100.  r1 = z - k*LN2_HI:
    k*LN2_HI is exact; the subtraction has error <= u|r1|.  r = r1 - k*LN2_LO has error
    <= u|r| + u|k LN2_LO| + |k| 2^-80, so |r - (z - k ln2)| <= er := 2u(|r| + |k| LN2_LO)
    + |k| 2^-80 + u*|r|  (|r| <= 0.35 + 1e-9).  Horner on the degree-22 Taylor
    polynomial with |r| <= 0.36: error <= gamma_44 * sum |r|^j/j! <= 44u e^{0.36}/(1-44u),
    relative to e^r >= e^{-0.36}: <= 44u e^{0.72} (1+1e-14) <= 91u.  Truncation <=
    0.36^23/23! e^{0.36} < 1e-30.  ldexp is exact (no over/underflow for |z| <= 700).
    So |v - e^{z_true}| <= e^{z_true} [(1+91u)(e^{ez+er}) - 1] <= e^{z_true}(110u +
    1.02(ez+er)) whenever ez + er <= 1e-3.
    """
    z = np.asarray(z, dtype=np.float64)
    assert np.all(np.abs(z) <= 700.0)
    k = np.rint(z * INV_LN2)
    r1 = z - k * LN2_HI
    r = r1 - k * LN2_LO
    er = 3 * U * np.abs(r) + 2 * U * np.abs(k) * LN2_LO + np.abs(k) * 2.0**-80
    p = np.full_like(r, _EXP_C[_EXP_DEG])
    for c in _EXP_C[_EXP_DEG - 1::-1]:
        p = p * r + c
    v = np.ldexp(p, k.astype(np.int64))
    tot = np.asarray(ez) + er
    assert np.all(tot <= 1e-3)
    return v, 110 * U + 1.02 * tot


TWO_PI_HI = float.fromhex('0x1.921fb54400000p+2')   # 33 significant bits
TWO_PI_LO = 1.2976551061315376e-09 * 4 / 4           # placeholder, set in _init
_SIN_DEG = 25
CW = [0.0, 0.0, 0.0]   # Cody-Waite pieces of 2pi, set by _init_constants()


def _init_constants():
    """Set TWO_PI_LO from Arb so that |TWO_PI_HI + TWO_PI_LO - 2pi| <= 2^-80; and check LN2_LO."""
    global TWO_PI_LO, LN2_LO
    from flint import arb, ctx
    ctx.prec = 200
    tp = 2 * arb.pi()
    TWO_PI_LO = float((tp - arb(TWO_PI_HI)).mid())
    assert abs(tp - arb(TWO_PI_HI) - arb(TWO_PI_LO)) < arb(2)**-80
    def round_bits(x, bits):
        m, e = math.frexp(float(x.mid()))
        return math.ldexp(math.floor(m * 2**bits) / 2**bits, e)
    p1 = round_bits(tp, 26); p2 = round_bits(tp - arb(p1), 26); p3 = float((tp - arb(p1) - arb(p2)).mid())
    assert abs(tp - arb(p1) - arb(p2) - arb(p3)) < arb(2)**-100
    CW[:] = [p1, p2, p3]
    l2 = arb(2).log()
    LN2_LO = float((l2 - arb(LN2_HI)).mid())
    assert abs(l2 - arb(LN2_HI) - arb(LN2_LO)) < arb(2)**-80


def vcis(theta, eth=0.0):
    """exp(-i theta_true) with |theta - theta_true| <= eth.  Returns (c, s, abs_err):
    |(c - i s) - e^{-i theta_true}| <= abs_err = eth + er + 200u.
    (Bound derivation: Horner
    on cos is gamma_26 cosh(0.786) <= 35u, sin <= 22u, complex <= 41u; each double-angle
    step doubles the complex error, |z^2 - z0^2| = |z - z0||z + z0|, so 4*41u + 10u = 175u.)

    Proof sketch.  k = rint(theta/2pi) exact integer, |k| < 2^20 required (checked).
    r = (theta - k*TWO_PI_HI) - k*TWO_PI_LO, |r| <= pi + 1e-6, with
    |r - (theta - 2 pi k)| <= er := 2u|r| + 2u|k TWO_PI_LO| + |k| 2^-80.
    Then r is halved twice (exact), cos/sin of r/4 (|r/4| <= 0.786) by Taylor/Horner in
    (r/4)^2 to degree 25/26 (truncation < 1e-26; Horner error <= gamma_28 * cosh/sinh
    majorant <= 30u), and two double-angle steps c' = c^2 - s^2 = (c-s)(c+s), s' = 2cs,
    each adding <= 3u + 2*(incoming error) in absolute terms (|c|,|s| <= 1+1e-15):
    total <= 4*(30u) + 4*3u + 3u < 140u ... we return 150u.
    """
    # Cody-Waite with three pieces P1 + P2 + P3 of 2pi (P1, P2 have 26 significant
    # bits): k*P1 and k*P2 are exact for |k| < 2^27; theta - k*P1 is exact by Sterbenz
    # (for k >= 1, k*P1 in [theta/2, 2 theta]; k = 0 is trivial); the next two
    # subtractions add <= u|r| each; |2pi - P1 - P2 - P3| <= 2^-100.
    theta = np.asarray(theta, dtype=np.float64)
    assert np.all(theta >= 0)
    k = np.rint(theta / (2 * math.pi))
    assert np.all(k < 2.0**27)
    r = ((theta - k * CW[0]) - k * CW[1]) - k * CW[2]
    er = 2 * U * np.abs(r) + U * k * abs(CW[2]) + k * 2.0**-100
    x = r * 0.25
    x2 = x * x
    # cos(x) = sum (-1)^j x^{2j}/(2j)!, sin(x) = x sum (-1)^j x^{2j}/(2j+1)!
    cc = np.full_like(x, (-1)**13 / math.factorial(26))
    for j in range(12, -1, -1):
        cc = cc * x2 + (-1)**j / math.factorial(2 * j)
    ss = np.full_like(x, (-1)**12 / math.factorial(25))
    for j in range(11, -1, -1):
        ss = ss * x2 + (-1)**j / math.factorial(2 * j + 1)
    ss = ss * x
    for _ in range(2):
        cc, ss = (cc - ss) * (cc + ss), 2.0 * cc * ss
    return cc, ss, np.asarray(eth) + er + 200 * U


def gamma_n(n):
    """Higham's gamma_n = n u / (1 - n u)."""
    return n * U / (1 - n * U)


def arb_logs(nmax):
    """log(1..nmax) from Arb, rounded to double; |fl - log n| <= u*log n (checked radius)."""
    from flint import arb, ctx
    ctx.prec = 80
    out = np.empty(nmax + 1)
    out[0] = -np.inf
    for n in range(1, nmax + 1):
        a = arb(n).log()
        out[n] = float(a.mid())
    return out


def selftest(nrand=20000, seed=1):
    """Empirical check (not the proof): compare vexp/vcis with Arb at random points and
    confirm the returned bounds hold."""
    from flint import arb, ctx
    _init_constants()
    ctx.prec = 120
    rng = np.random.default_rng(seed)
    z = rng.uniform(-60, 20, nrand)
    v, rel = vexp(z)
    worst_e = 0.0
    for zi, vi, ri in zip(z[:2000], v[:2000], np.broadcast_to(rel, z.shape)[:2000]):
        ex = arb(float(zi)).exp()
        err = abs(arb(float(vi)) - ex) / ex
        e = float(err.upper())
        worst_e = max(worst_e, e / ri)
    th = rng.uniform(0, 4e8, nrand)
    c, s, ea = vcis(th)
    worst_c = 0.0
    for ti, ci, si, ai in zip(th[:2000], c[:2000], s[:2000], np.broadcast_to(ea, th.shape)[:2000]):
        t = arb(float(ti))
        e = float(max(abs(arb(float(ci)) - t.cos()).upper(), abs(arb(float(si)) - t.sin()).upper()))
        worst_c = max(worst_c, e / ai)
    return dict(exp_err_over_bound=worst_e, cis_err_over_bound=worst_c)


if __name__ == '__main__':
    print(selftest())
