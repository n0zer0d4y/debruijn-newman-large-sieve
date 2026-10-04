"""Upper bound, part (b), float reconnaissance (numeric-experiment, not a certificate): the RS main term
g = A(s*) + gamma B(w) of the Davenport-Heilbronn heat flow at t0 along the rectangle edges
y = y0 and y = Ytop, x in [X1, X2], and along the vertical sides x = X1, X2 (y in [y0, Ytop]).
Reports min |g| per edge and the float winding number of B_t g around [X1, X2] x [y0, Ytop]
(arg B_t from its analytic Stirling logarithm lMt; arg g by accumulating principal steps).
Conventions as in cert_barrier / cert_barrier2 (E = 1: no Euler product for DH).
"""
import sys, math, json
import numpy as np

q, A_, t0 = 5, 1, 0.6
k = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
TAB = [0.0, 1.0, k, -k, -1.0]
LQ = math.log(q / math.pi)


def am(s):
    sa = s + A_
    return LQ / 2 + np.log(sa / 2) / 2 - 1 / (2 * sa)


def lM0(s):
    w = (s + A_) / 2
    return w * LQ + math.log(math.sqrt(2 * math.pi)) + (w - 0.5) * np.log(w) - w


def lMt(s, t):
    return (t / 4) * am(s) ** 2 + lM0(s)


def g_and_logB(x, y, t=t0):
    x = np.asarray(x, float); y = np.broadcast_to(np.asarray(y, float), x.shape)
    sp = (1 + y - 1j * x) / 2; sm = (1 - y + 1j * x) / 2
    ss = sp + (t / 2) * am(sp); w = sm + (t / 2) * am(sm)
    gam = np.exp(lMt(sm, t) - lMt(sp, t))
    N = np.floor(np.sqrt(q * x / (4 * math.pi) + q * t / 16)).astype(int)
    g = np.zeros(x.shape, complex)
    for n in range(1, int(N.max()) + 1):
        a = TAB[n % 5]
        if a == 0: continue
        L = math.log(n); b = math.exp(t / 4 * L * L)
        m = N >= n
        g += np.where(m, a * b * (np.exp(-ss * L) + gam * np.exp(-w * L)), 0)
    return g, lMt(sp, t)


def arg_change(v):
    d = np.angle(v[1:] / v[:-1])
    return float(np.sum(d)), float(np.max(np.abs(d)))


def edge(xs, ys):
    g, lB = g_and_logB(xs, ys)
    dg, mx = arg_change(g)
    dB = float(lB[-1].imag - lB[0].imag)
    return dict(min_abs_g=float(np.abs(g).min()), at=float(xs[np.argmin(np.abs(g))]) if np.ndim(xs) else None,
                darg_g=dg, darg_B=dB, max_step=mx)


if __name__ == '__main__':
    X1 = float(sys.argv[1]) if len(sys.argv) > 1 else 280.0
    X2 = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0e4
    if len(sys.argv) > 3: t0 = float(sys.argv[3])
    y0 = float(sys.argv[4]) if len(sys.argv) > 4 else 0.4
    Ytop = math.sqrt(1.7902723164702195 ** 2 - 2 * t0)
    G_ = lambda X, Y: g_and_logB(X, Y, t0)
    h = float(sys.argv[5]) if len(sys.argv) > 5 else 0.01
    xs = np.linspace(X1, X2, int((X2 - X1) / h) + 1)
    ys = np.linspace(y0, Ytop, 400)
    out = {}
    tot = 0.0
    # counter-clockwise: bottom (x up), right (y up), top (x down), left (y down)
    parts = [('bottom', xs, np.full_like(xs, y0)), ('right', np.full_like(ys, X2), ys),
             ('top', xs[::-1], np.full_like(xs, Ytop)), ('left', np.full_like(ys, X1), ys[::-1])]
    for name, X, Y in parts:
        g, lB = G_(X, Y)
        dg, mx = arg_change(g)
        dB = float(lB[-1].imag - lB[0].imag)
        out[name] = dict(min_abs_g=float(np.abs(g).min()), at=[float(X[np.argmin(np.abs(g))]), float(Y[np.argmin(np.abs(g))])],
                         darg_g=dg, darg_B=dB, max_step=mx)
        tot += dg + dB
    out['winding'] = tot / (2 * math.pi)
    # bottom-edge profile: the 10 smallest local minima of |g|
    gb, _ = G_(xs, np.full_like(xs, y0))
    a = np.abs(gb); loc = np.where((a[1:-1] < a[:-2]) & (a[1:-1] < a[2:]))[0] + 1
    sel = loc[np.argsort(a[loc])[:10]]
    out['bottom_smallest_minima'] = [[float(xs[i]), float(a[i])] for i in sel]
    out['t0'], out['y0'] = t0, y0
    print(json.dumps(out))
