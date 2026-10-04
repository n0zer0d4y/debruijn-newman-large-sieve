# Driver for real primitive Dirichlet characters (not used by the certificates of this repository; cert_final imports kron() from here): run(D, T0) = barrier location + final region + rigorous K_h (Lemma H) + e_C + barrier.
# Paste into a local Python run followed by e.g. `for D,T0 in [(24,0.23)]: print(run(D,T0))`. Float64/mpmath (numeric-experiment).
import os, itertools, math, time
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import mpmath as mp, numpy as np
from math import gcd
mp.mp.dps = 20
def kron(D, n):
    if n == 0 or gcd(abs(D), n) != 1: return 0
    res = 1; m = n
    while m % 2 == 0:
        m //= 2; res *= 1 if D % 8 in (1, 7) else -1
    a = D % m; r = 1; mm = m
    while a != 0:
        while a % 2 == 0:
            a //= 2
            if mm % 8 in (3, 5): r = -r
        a, mm = mm, a
        if a % 4 == 3 and mm % 4 == 3: r = -r
        a %= mm
    return res * (r if mm == 1 else 0)
def run(D, T0, tb=260, t_start=None):
    T_all = time.time(); q = abs(D); A_ = 0 if D > 0 else 1; TAB = [kron(D, n) for n in range(q)]; chi = lambda n: TAB[n % q]
    XPL = 2 * max(1e8 / q, (7.5 if q % 2 == 0 else 3.75) * 1e7 / q + 200)
    am = lambda s: mp.log(q / mp.pi) / 2 + mp.log((s + A_) / 2) / 2 - 1 / (2 * (s + A_))
    def lM0(s):
        w = (s + A_) / 2; return (s + A_) / 2 * mp.log(q / mp.pi) + mp.log(mp.sqrt(2 * mp.pi)) + (w - mp.mpf(1) / 2) * mp.log(w) - w
    lMt = lambda s, t: (t / 4) * am(s) ** 2 + lM0(s)
    def sg(x, y, t):
        x, y, t = mp.mpf(x), mp.mpf(y), mp.mpf(t); sp = (1 + y - 1j * x) / 2; sm = (1 - y + 1j * x) / 2
        return float(mp.re(sp + (t / 2) * am(sp))), float(abs(mp.exp(lMt(sm, t) - lMt(sp, t))))
    xr = lambda a, b, t: (4 * math.pi * (a ** 2 - q * t / 16) / q, 4 * math.pi * ((b + 1) ** 2 - q * t / 16) / q)
    Nof = lambda x, t: int(math.floor(math.sqrt(q * x / (4 * math.pi) + q * t / 16)))
    MS = np.array([k for k in range(1, 1000002) if chi(k) != 0], dtype=np.int64); CHS = np.array([chi(int(k)) for k in MS], dtype=float)
    cnt = lambda N: int(np.searchsorted(MS, N, side='right'))
    def moll(P, t):
        b = lambda n: math.exp(t / 4 * math.log(n) ** 2); lam = {}
        for r in range(len(P) + 1):
            for S in itertools.combinations(P, r):
                d, l = 1, 1.0
                for p in S: d *= p; l *= -chi(p) * b(p)
                lam[d] = l
        return lam
    def blk(a, b, y1, y2, t, P):
        xlo, xhi = xr(a, b, t); s0, _ = sg(xlo, y1, t); cor = [sg(xx, yy, t)[1] for xx in (xlo, xhi) for yy in (y1, y2)]
        lam = moll(P, t); m = MS[:cnt(b)]; lm = np.log(m.astype(float)); bm = np.exp(t / 4 * lm ** 2); chm = CHS[:len(m)]
        gs = np.array([g * np.exp(yy * lm) for g, yy in zip(cor, (y1, y2, y1, y2))]); gmin, gmax = gs.min(0), gs.max(0); g0 = np.sqrt(gmin * gmax); w = np.maximum(gmax - g0, g0 - gmin)
        inb = m <= a; idx, be, al, dv = [], [], [], []
        for d, ld in lam.items():
            idx.append(d * m[inb]); be.append(ld * chm[inb] * bm[inb]); al.append(ld * chm[inb] * bm[inb] * g0[inb]); dv.append(abs(ld) * bm[inb] * w[inb])
        idx = np.concatenate(idx); u, inv = np.unique(idx, return_inverse=True)
        B = np.bincount(inv, np.concatenate(be)); A = np.bincount(inv, np.concatenate(al)); Dv = np.bincount(inv, np.concatenate(dv))
        k = (1 - gmin[0]) / (1 + gmin[0]); wt = np.exp(-s0 * np.log(u.astype(float))); mk = u > 1
        Y = float(np.sum((np.maximum(np.abs(B - A), k * np.abs(B + A)) + Dv)[mk] * wt[mk]))
        la = sum(abs(ld) * d ** (-s0) for d, ld in lam.items()); tl = ~inb
        tr = la * float(np.sum(bm[tl] * (1 + gmax[tl]) * np.exp(-s0 * lm[tl])))
        return 1 - gmax[0] - Y - tr, tr
    def pln(a, b, y1, y2, t):
        xlo, xhi = xr(a, b, t); s0, _ = sg(xlo, y1, t); cor = [(sg(xx, yy, t)[1], yy) for xx in (xlo, xhi) for yy in (y1, y2)]
        n = MS[:cnt(b)].astype(float); ln = np.log(n); wt = np.exp(t / 4 * ln ** 2 - s0 * ln); gm = np.max([g * np.exp(yy * ln) for g, yy in cor], axis=0)
        return 1 - (np.sum(wt) - 1) - np.sum(gm * wt)
    # barrier location
    Pp = np.array([p for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] if chi(p) != 0]); cP = np.array([chi(int(p)) for p in Pp])
    xs = np.arange(XPL - 2e5, XPL - 1.0, 0.05); val = np.zeros_like(xs)
    for p, c in zip(Pp, cP): val += -np.log(np.abs(1 - c * p ** (-0.6) * np.exp(-1j * xs / 2 * np.log(p))))
    XC = float(xs[int(np.argmax(val))])
    pr = [p for p in [2, 3, 5, 7, 11, 13, 17] if chi(p) != 0]; ps, pb = pr[:5], pr[:3]
    # final region
    t = T0; ytop = math.sqrt(1 - 2 * t); Ns = Nof(XC + 0.5, t); Nb = []; a = Ns
    while a < 100001:
        b = min(max(int(a * 1.02), a + 1), 100001); Nb.append((a, b)); a = b + 1
    yc = []; y = 0.2
    while y < ytop - 1e-12:
        h = 0.01 if y < 0.3 else 0.05; yc.append((y, min(y + h, ytop))); y = min(y + h, ytop)
    Q = [(a, b, y1, y2) for (a, b) in Nb for (y1, y2) in yc]; mn = 9; fails = 0; nev = 0
    while Q:
        a, b, y1, y2 = Q.pop(); v, tr = blk(a, b, y1, y2, t, ps if a < 15000 else pb); nev += 1
        if v > 0.01: mn = min(mn, v); continue
        if nev > 9000 or (b - a <= 1 and y2 - y1 < 0.0005): fails += 1; continue
        if b - a > 1 and tr > 0.3 * max(0.01 - v, 1e-3):
            mid = (a + b) // 2; Q += [(a, mid, y1, y2), (mid + 1, b, y1, y2)]
        else:
            ym = (y1 + y2) / 2; Q += [(a, b, y1, ym), (a, b, ym, y2)]
    pm = 9; a = 100001
    while a < 1000001:
        b = min(int(a * 1.05), 1000001); y = 0.2
        while y < ytop - 1e-12:
            y2 = min(y + (0.02 if y < 0.3 else 0.1), ytop); pm = min(pm, pln(a, b, y, y2, t)); y = y2
        a = b + 1
    N0 = 1000001; xlo, _ = xr(N0, N0, t); s0, g0_ = sg(xlo, 0.2, t); L0 = math.log(N0)
    tails = max(math.exp((1 - s0) * L0 + t / 4 * L0 ** 2), g0_ * math.exp((1 - (s0 - 0.2)) * L0 + t / 4 * L0 ** 2))
    # rigorous K_h (Lemma H) with symmetrisation
    chi2 = [chi(n) for n in range(2 * q)]; c0 = abs(sum(chi2[n] * complex(math.cos(math.pi * n * n / q), -math.sin(math.pi * n * n / q)) for n in range(2 * q)) / (2 * q))
    ns = [n for n in range(2 * q) if chi2[n] != 0]; s2 = math.sqrt(2); us = np.arange(-6, 6.00002, 2e-5); best = 0; dbest = 0
    for kk in range(2 * q):
        cc = kk + 0.5; tot = np.zeros_like(us); dtot = np.zeros_like(us)
        for n in ns:
            aa = math.pi * (cc - n + us / s2) / q; bb = math.pi * us / (s2 * q); eb = np.exp(-bb); F = (1 - eb) ** 2 + 4 * eb * np.sin(aa / 2) ** 2
            dF = (math.pi / (s2 * q)) * (2 * np.abs(1 - eb) * eb + 4 * eb * np.sin(aa / 2) ** 2 + 2 * eb * np.abs(np.sin(aa)))
            tot += F ** -0.5; dtot += 0.5 * F ** -1.5 * dF
        best = max(best, float(tot.max())); dbest = max(dbest, float(dtot.max()))
    Kh = max(best / (2 * q) + 2 * dbest / (2 * q) * 1e-5, len(ns) / (2 * q) / (1 - math.exp(-6 * math.pi / (s2 * q)))) + c0 / 2
    def lemI(sig, Tp, c, rho=0.25):
        Qs = 4.0 / q; lam = q * Tp / (2 * math.pi * c * c); l0 = Qs * math.pi * c * (lam - 1) / (2 * s2); l = l0 - sig / (s2 * c)
        k = Qs * math.pi * (1 + lam) / 4; G0 = math.sqrt(math.pi / k) * math.exp(l * l / (4 * k)); A = Qs * math.pi * lam / (6 * c * (1 - rho)) + abs(sig) / (3 * c ** 3 * (1 - rho)); kp = k - A * rho * c
        return G0 + A * math.exp(l * l / (4 * kp)) * 4 * (1 / kp ** 2 + abs(l / (2 * kp)) ** 3 * math.sqrt(math.pi / kp))
    gx, gw = np.polynomial.hermite.hermgauss(60)
    def eC(x, y, t):
        Tp = x / 2 + math.pi * t / 8; wR = math.sqrt(q * Tp / (2 * math.pi)); c = math.floor(wR) + 0.5; de = wR - c
        ratio = float(mp.exp(mp.re(lM0(1j * mp.mpf(Tp))) - mp.re(lMt((1 + mp.mpf(y) - 1j * mp.mpf(x)) / 2, t)))); tot = 0
        for sig in [(1 - y) / 2, (1 + y) / 2]:
            f = lambda sv: (1 + de / c) ** sv * lemI(sv, Tp, c) * math.exp(1 / (6 * (Tp - 0.66)) + abs(sv) / Tp ** 2 + sv * sv / (4 * (Tp - 1)))
            I = float(np.dot(gw, [f(sig + math.sqrt(t) * v) for v in gx]) / math.sqrt(math.pi)) if t > 0 else f(sig)
            tot += math.exp(math.pi ** 2 * t / 64) * ratio * I * Kh
        return tot
    eCm = max(eC(XC, math.sqrt(0.04 + 2 * T0), 0.0), eC(XC, 1.0, 0.0), eC(XC, 0.2, T0))
    # barrier
    N = int(math.floor(math.sqrt(q * XC / (4 * math.pi)))); Nck = (Nof(XC + 0.5, T0), int(math.floor(math.sqrt(q * (XC - 0.5) / (4 * math.pi)))))
    P = np.array([p for p in [2, 3, 5, 7, 11, 13, 17, 19, 23] if chi(p) != 0][:7], dtype=float); CP = np.array([chi(int(p)) for p in P], dtype=float); LP = np.log(P)
    nn = np.array([k for k in range(1, N + 1) if chi(k) != 0], dtype=float); LN = np.log(nn); CN = np.array([chi(int(k)) for k in nn], dtype=float)
    al_ = lambda s: np.log(q / np.pi) / 2 + np.log((s + A_) / 2) / 2 - 1 / (2 * (s + A_))
    def lm0(s):
        w = (s + A_) / 2; return (s + A_) / 2 * np.log(q / np.pi) + 0.5 * np.log(2 * np.pi) + (w - 0.5) * np.log(w) - w
    lmt = lambda s, t: (t / 4) * al_(s) ** 2 + lm0(s)
    def Ef(x, y, t):
        sp = (1 + y - 1j * x) / 2; sm = (1 - y + 1j * x) / 2; gam = np.exp(lmt(sm, t) - lmt(sp, t)); ss = sp + (t / 2) * al_(sp); kap = (t / 2) * (al_(sm) - al_((1 + y + 1j * x) / 2))
        sr, si = ss.real, ss.imag; b = CN * np.exp(t / 4 * LN ** 2); out = np.empty(len(x), dtype=complex)
        for i0 in range(0, len(x), 1024):
            sl = slice(i0, i0 + 1024); ph = np.exp(-1j * np.mod(np.outer(si[sl], LN), 2 * np.pi))
            A = (b * np.exp(-np.outer(sr[sl], LN)) * ph).sum(1); Bs = (b * np.exp(np.outer(y[sl] - sr[sl], LN) - np.outer(kap[sl], LN)) * np.conj(ph)).sum(1)
            E = np.prod(1 - CP * np.exp(t / 4 * LP ** 2 - np.outer(sr[sl], LP) - 1j * np.mod(np.outer(si[sl], LP), 2 * np.pi)), axis=1); out[sl] = E * (A + gam[sl] * Bs)
        return out
    def ga(y, t):
        sp = (1 + y - 1j * XC) / 2; sm = (1 - y + 1j * XC) / 2; return abs(np.exp(lmt(sm, t) - lmt(sp, t)))
    def lip(yl, tl, th):
        sp = (1 + yl - 1j * XC) / 2; sm = (1 - yl + 1j * XC) / 2; ap, amm = al_(sp), al_(sm); sr = (sp + (tl / 2) * ap).real
        gl = max(ga(yl, tl), ga(yl, th)) * 1.01; g1 = max(ga(1.0, tl), ga(1.0, th)) * 1.01
        b = np.exp(th / 4 * LN ** 2); wA = b * np.exp(-sr * LN); wB = b * np.exp(-sr * LN) * np.maximum(gl * np.exp(yl * LN), g1 * nn)
        Gxy = 0.5 * (abs(ap) + abs(amm)) + 1; Gt = abs(amm ** 2 - ap ** 2) / 4 + 1; aa = max(abs(ap), abs(amm))
        F = wA.sum() + wB.sum(); Fs = (wA * 0.51 * LN).sum() + (wB * (1.51 * LN + Gxy)).sum(); Ft = (wA * (LN ** 2 / 4 + aa * LN / 2 + 0.1)).sum() + (wB * (LN ** 2 / 4 + aa * LN / 2 + 0.1 + Gt)).sum()
        bp = np.exp(th / 4 * LP ** 2 - sr * LP); Em = np.prod(1 + bp); Es = Em * np.sum(bp * 0.51 * LP / (1 - bp)); Et = Em * np.sum(bp * (LP ** 2 / 4 + aa * LP / 2) / (1 - bp))
        return Es * F + Em * Fs, Et * F + Em * Ft, Em
    yl_ = lambda t: math.sqrt(0.04 + 2 * (T0 - t))
    t = 0.0; cells = 0; worst = 9; mind = 9; Emax = 0; ok = True
    while t < T0 - 1e-12:
        tau = 0.02
        while True:
            th = min(t + tau, T0); yl = yl_(th); Ls, Lt, Em = lip(yl, t, th)
            if Lt * (th - t) / 2 < 0.3: break
            tau /= 1.5
        Emax = max(Emax, Em); err = Em * eCm * 1.2 + 1e-6
        h = min(0.05, 0.6 / Ls); xs = np.arange(XC - 0.5 + h / 2, XC + 0.5, h); ys = np.arange(yl + h / 2, 1.0, h)
        X = np.concatenate([xs, xs, np.full(len(ys), XC - 0.5), np.full(len(ys), XC + 0.5)]); Y = np.concatenate([np.full(len(xs), yl), np.full(len(xs), 1.0), ys, ys])
        v = Ef(X, Y, (t + th) / 2); d = np.where(v.real >= 0, np.abs(v), np.abs(v.imag)); mind = min(mind, float(d.min()))
        mg = d.min() - err - Ls * h / 2 - Lt * (th - t) / 2; worst = min(worst, float(mg)); cells += 1
        if mg <= 0: ok = False; break
        t = th
        if time.time() - T_all > tb: ok = False; break
    return dict(D=D, t0=T0, Lambda=round(T0 + 0.02, 3), XC=XC, N=N, Nck=Nck, FR_min=round(float(mn), 4), FR_fails=fails, plain_min=round(float(pm), 4), tail=round(tails, 5), Kh_rig=round(Kh, 4), c0=round(c0, 4), eCmax=f'{eCm:.2e}', barrier_ok=ok, barrier_to=round(t, 4), cells=cells, barrier_worst=round(worst, 4), mind=round(mind, 4), supE=round(float(Emax), 3), secs=round(time.time() - T_all, 1))
