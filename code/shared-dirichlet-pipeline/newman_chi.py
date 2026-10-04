"""General pipeline for a real primitive character chi mod q, parity a (used here only by code/03-davenport-heilbronn/dh_final_region_scan.py).

Conventions: s+ = (1+y-ix)/2, s- = (1-y+ix)/2, s* = s+ + (t/2) alpha(s+), sigma = Re s*,
gamma = M_t(s-)/M_t(s+), w_R^2 = q T'/(2 pi) = q x/(4 pi) + q t/16, N = floor(w_R),
RS contour through c = N + 1/2 (half-integer), delta = w_R - c in [-1/2, 1/2).
Kernel K = e^{i pi w^2/q} h(w), h(w) = sum_j c_j e^{pi i j w/q}/(e^{2 pi i w}-1) (period 2q in w).
Float64 / mpmath exploration grade (numeric-experiment).
"""
import itertools
import math

import mpmath as mp
import numpy as np

mp.mp.dps = 20


class Chi:
    def __init__(self, q, a, table):
        self.q, self.a, self.table = q, a, table

    def __call__(self, n):
        return self.table[n % self.q]

    # ---- gamma-factor model (derivation §2.1) ----
    def alpha(self, s):
        q, a = self.q, self.a
        return mp.log(q / mp.pi) / 2 + mp.log((s + a) / 2) / 2 - 1 / (2 * (s + a))

    def logM0(self, s):
        q, a = self.q, self.a
        w = (s + a) / 2
        return (s + a) / 2 * mp.log(q / mp.pi) + mp.log(mp.sqrt(2 * mp.pi)) + (w - mp.mpf(1) / 2) * mp.log(w) - w

    def logMt(self, s, t):
        return (t / 4) * self.alpha(s) ** 2 + self.logM0(s)

    def sig_gam(self, x, y, t):
        x, y, t = mp.mpf(x), mp.mpf(y), mp.mpf(t)
        sp = (1 + y - 1j * x) / 2
        sm = (1 - y + 1j * x) / 2
        return float(mp.re(sp + (t / 2) * self.alpha(sp))), float(abs(mp.exp(self.logMt(sm, t) - self.logMt(sp, t))))

    def x_range(self, Nlo, Nhi, t):
        q = self.q
        return 4 * math.pi * (Nlo ** 2 - q * t / 16) / q, 4 * math.pi * ((Nhi + 1) ** 2 - q * t / 16) / q

    # ---- kernel (derivation §2.6a) ----
    def kernel_coeffs(self):
        q = self.q
        M = 2 * q
        A = mp.matrix(M, M)
        b = mp.matrix(M, 1)
        for n in range(M):
            for j in range(M):
                A[n, j] = mp.exp(1j * mp.pi * n * n / q) * mp.exp(mp.pi * 1j * j * n / q)
            b[n] = self(n)
        return mp.lu_solve(A, b)

    def h_sup(self, U0=3.0, du=0.002):
        """sup over c in (Z + 1/2) and u of |h(c + u e^{i pi/4})|: grid on |u| <= U0 (+ Lipschitz slack),
        and the analytic tail bound sum|c_j| / (1 - e^{-sqrt2 pi U0}) for |u| > U0 (derivation §4.1)."""
        q = self.q
        c = [complex(v) for v in self.kernel_coeffs()]
        d = complex(mp.exp(1j * mp.pi / 4))
        us = np.arange(-U0, U0 + du, du)
        best = 0.0
        for cc in [k + 0.5 for k in range(2 * q)]:
            w = cc + us * d
            num = sum(cj * np.exp(1j * np.pi * j * w / q) for j, cj in enumerate(c))
            hv = np.abs(num / (np.exp(2j * np.pi * w) - 1))
            best = max(best, float(hv.max()))
        # |u| > U0: Im w = u/sqrt2. u > U0: |e^{pi i j w/q}| = e^{-pi j u/(sqrt2 q)}, |e^{2 pi i w}-1| >= 1-e^{-sqrt2 pi u};
        # u < -U0: divide by e^{sqrt2 pi |u|}: |c_j| e^{-pi (2q-j)|u|/(sqrt2 q)} / (1 - e^{-sqrt2 pi |u|}); both decrease in |u|.
        den = 1 - math.exp(-math.sqrt(2) * math.pi * U0)
        tail_pos = sum(abs(cj) * math.exp(-math.pi * j * U0 / (math.sqrt(2) * q)) for j, cj in enumerate(c)) / den
        tail_neg = sum(abs(cj) * math.exp(-math.pi * (2 * q - j) * U0 / (math.sqrt(2) * q)) for j, cj in enumerate(c)) / den
        return best, max(tail_pos, tail_neg)

    def N_of(self, x, t):
        return int(math.floor(math.sqrt(self.q * x / (4 * math.pi) + self.q * t / 16)))

    # ---- Lemma I_q (derivation §2.7 scaled by Qs = 4/q) ----
    def lemmaI(self, sig, Tp, c, rho=0.25):
        Qs = 4.0 / self.q
        lam = self.q * Tp / (2 * math.pi * c * c)
        l0 = Qs * math.pi * c * (lam - 1) / (2 * math.sqrt(2))
        l = l0 - sig / (math.sqrt(2) * c)
        k = Qs * math.pi * (1 + lam) / 4
        G0 = math.sqrt(math.pi / k) * math.exp(l * l / (4 * k))
        A = Qs * math.pi * lam / (6 * c * (1 - rho)) + abs(sig) / (3 * c ** 3 * (1 - rho))
        kp = k - A * rho * c
        E = A * math.exp(l * l / (4 * kp)) * 4 * (1 / kp ** 2 + abs(l / (2 * kp)) ** 3 * math.sqrt(math.pi / kp))
        return G0 + E

    def eC(self, x, y, t, Kh):
        gx, gw = np.polynomial.hermite.hermgauss(60)
        Tp = x / 2 + math.pi * t / 8
        wR = math.sqrt(self.q * Tp / (2 * math.pi))
        c = math.floor(wR) + 0.5
        delta = wR - c
        ratio = float(mp.exp(mp.re(self.logM0(1j * mp.mpf(Tp))) - mp.re(self.logMt((1 + mp.mpf(y) - 1j * mp.mpf(x)) / 2, t))))
        tot = 0.0
        for sig in [(1 - y) / 2, (1 + y) / 2]:
            f = lambda sv: (1 + delta / c) ** sv * self.lemmaI(sv, Tp, c) * math.exp(1 / (6 * (Tp - 0.66)) + abs(sv) / Tp ** 2 + sv * sv / (4 * (Tp - 1)))
            I = float(np.dot(gw, [f(sig + math.sqrt(t) * v) for v in gx]) / math.sqrt(math.pi)) if t > 0 else f(sig)
            tot += math.exp(math.pi ** 2 * t / 64) * ratio * I * Kh
        return tot

    # ---- final region (derivation §3.3 logic) ----
    def mollifier(self, primes, t):
        b = lambda n: math.exp(t / 4 * math.log(n) ** 2)
        lam = {}
        for r in range(len(primes) + 1):
            for S in itertools.combinations(primes, r):
                d, l = 1, 1.0
                for p in S:
                    d *= p
                    l *= -self(p) * b(p)
                lam[d] = l
        return lam

    def block_lb(self, Nlo, Nhi, ylo, yhi, t, primes):
        xlo, xhi = self.x_range(Nlo, Nhi, t)
        sig_min, _ = self.sig_gam(xlo, ylo, t)
        corners = [self.sig_gam(xx, yy, t)[1] for xx in (xlo, xhi) for yy in (ylo, yhi)]
        lam = self.mollifier(primes, t)
        m = np.array([k for k in range(1, Nhi + 1) if self(k) != 0], dtype=np.int64)
        lm = np.log(m.astype(float))
        bm = np.exp(t / 4 * lm ** 2)
        chm = np.array([self(int(k)) for k in m], dtype=float)
        gs = np.array([g * np.exp(yy * lm) for g, yy in zip(corners, (ylo, yhi, ylo, yhi))])
        gmin, gmax = gs.min(0), gs.max(0)
        g0 = np.sqrt(gmin * gmax)
        w = np.maximum(gmax - g0, g0 - gmin)
        inb = m <= Nlo
        idx, bet, alp, dev = [], [], [], []
        for d, ld in lam.items():
            idx.append(d * m[inb]); bet.append(ld * chm[inb] * bm[inb])
            alp.append(ld * chm[inb] * bm[inb] * g0[inb]); dev.append(abs(ld) * bm[inb] * w[inb])
        idx = np.concatenate(idx)
        uniq, inv = np.unique(idx, return_inverse=True)
        B = np.bincount(inv, np.concatenate(bet)); A = np.bincount(inv, np.concatenate(alp)); Dv = np.bincount(inv, np.concatenate(dev))
        a1min, a1max = gmin[0], gmax[0]
        k = (1 - a1min) / (1 + a1min)
        wt = np.exp(-sig_min * np.log(uniq.astype(float)))
        mask = uniq > 1
        Y = float(np.sum((np.maximum(np.abs(B - A), k * np.abs(B + A)) + Dv)[mask] * wt[mask]))
        lam_abs = sum(abs(ld) * d ** (-sig_min) for d, ld in lam.items())
        tail = ~inb
        trunc = lam_abs * float(np.sum(bm[tail] * (1 + gmax[tail]) * np.exp(-sig_min * lm[tail])))
        return 1 - a1max - Y - trunc, trunc

    def plain_block(self, Nlo, Nhi, ylo, yhi, t):
        xlo, xhi = self.x_range(Nlo, Nhi, t)
        sig, _ = self.sig_gam(xlo, ylo, t)
        cor = [(self.sig_gam(xx, yy, t)[1], yy) for xx in (xlo, xhi) for yy in (ylo, yhi)]
        n = np.array([k for k in range(1, Nhi + 1) if self(k) != 0], dtype=float)
        ln = np.log(n)
        wt = np.exp(t / 4 * ln ** 2 - sig * ln)
        gmax = np.max([g * np.exp(yy * ln) for g, yy in cor], axis=0)
        return 1 - (np.sum(wt) - 1) - np.sum(gmax * wt)
