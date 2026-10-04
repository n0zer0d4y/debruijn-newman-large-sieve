# (numeric-experiment) Davenport-Heilbronn: see code/03-davenport-heilbronn/notes.md.
import subprocess
import sys as _sys
PY = _sys.executable   # any Python with python-flint, mpmath, numpy
code = r'''
from flint import arb, ctx
ctx.prec = 128
K = ((10 - 2*arb(5).sqrt()).sqrt() - 2) / (arb(5).sqrt() - 1)
print('kappa =', K)
def S(s):
    s = arb(s)
    return arb(5)**(-s) * (s.zeta(arb(1)/5) + s.zeta(arb(4)/5) + K*(s.zeta(arb(2)/5) + s.zeta(arb(3)/5))) - 1
lo, hi = arb(1.05), arb(3)
for _ in range(60):
    m = (lo + hi) / 2
    if S(m) < 1: hi = m
    else: lo = m
print('sigma0 upper =', hi.upper(), ' S(sigma0) =', S(hi.upper()))
D = 2*hi.upper() - 1
print('Delta =', D, ' Lambda_DH <= Delta^2/2 =', (D*D/2).upper())
'''
r = subprocess.run([PY, '-c', code], capture_output=True, text=True, timeout=120)
print(r.stdout, r.stderr[-500:])
