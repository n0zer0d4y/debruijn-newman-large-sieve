# (numeric-experiment) Davenport-Heilbronn: see code/03-davenport-heilbronn/notes.md.
import mpmath as mp, time
mp.mp.dps = 18
k = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
A = [0, 1, k, -k, -1]
def f(s):
    return sum(A[r] * mp.zeta(s, mp.mpf(r)/5) for r in range(1, 5)) * mp.mpf(5)**(-s)
T = time.time()
for g in [(0.8625, 1354.8), (0.8875, 1464.5), (0.9125, 2442.5)]:
    best = None
    for ds in (-0.15, 0.0, 0.15):
        for dt in (-0.4, 0.0, 0.4):
            try:
                r = mp.findroot(f, mp.mpc(g[0]+ds, g[1]+dt), tol=1e-12)
                if abs(r.imag - g[1]) < 3 and abs(f(r)) < 1e-10:
                    if best is None or r.real > best.real: best = r
            except Exception: pass
    print(g, '->', None if best is None else (float(best.real), float(best.imag)), round(time.time()-T), 's', flush=True)
