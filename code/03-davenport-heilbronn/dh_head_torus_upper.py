# (numeric-experiment) Davenport-Heilbronn: see code/03-davenport-heilbronn/notes.md.
import numpy as np, math, os
os.environ['OPENBLAS_NUM_THREADS']='1'
k = (math.sqrt(10-2*math.sqrt(5))-2)/(math.sqrt(5)-1)
A = {1:1.0, 2:k, 3:-k, 4:-1.0, 0:0.0}
NMAX = 200000
n = np.arange(1, NMAX+1)
a = np.array([A[i % 5] for i in n])
def smooth23(m):
    while m % 2 == 0: m //= 2
    while m % 3 == 0: m //= 3
    return m == 1
sm = np.array([smooth23(int(i)) for i in range(1, 5001)] + [False]*(NMAX-5000))
g = np.linspace(0, 2*np.pi, 721)
U, V = np.meshgrid(np.exp(1j*g), np.exp(1j*g), indexing='ij')
idx = [i for i in range(1, 5001) if sm[i-1] and a[i-1] != 0]
exps = []
for i in idx:
    e2 = e3 = 0; m = i
    while m % 2 == 0: m//=2; e2+=1
    while m % 3 == 0: m//=3; e3+=1
    exps.append((i, e2, e3))
def bound(sig):
    P = np.zeros_like(U)
    for i, e2, e3 in exps:
        P += a[i-1] * i**(-sig) * U**e2 * V**e3
    tail = float(np.sum(np.abs(a[~sm]) * n[~sm]**(-sig))) + (NMAX**(1-sig))/(sig-1)
    return float(np.abs(P).min()), tail
for sig in [1.40, 1.30, 1.25, 1.20, 1.15, 1.10]:
    m, tl = bound(sig)
    print('sigma', sig, 'min|P|', round(m,4), 'tail', round(tl,4), 'zero-free' if m > tl else 'not shown')
