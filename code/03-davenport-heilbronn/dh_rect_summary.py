"""assemble the Lambda_DH <= t0 + y0^2/2 certificate (default: the 0.58 run, t0 = 0.5, y0 = 0.4, X = 5e4) from
  part (a): dh_upper.final_region()  (re-run here; H_t0 != 0 for x >= X, y0 <= y <= Ytop, with |G - 1| < 1);
  part (b): the edge runs of dh_rect.py in rect_runs_050/*.json (G = H/B_t avoids (-inf, 0] on each boundary piece).
Checks every run's ok flag, the coverage of the four edges of R = [0, X] x [y0, Ytop], and the common
parameters, then writes certificates/03-davenport-heilbronn/lambda-058.json.  Single process.
"""
import json, glob, os, sys, math, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dh_upper, dh_rect

prm = dict(a.split('=') for a in sys.argv[1:])
T0_, Y0_, X_ = float(prm.get('t0', 0.5)), float(prm.get('y0', 0.4)), float(prm.get('X', 5e4))
RUNDIR = prm.get('runs', 'rect_runs_050'); OUT = prm.get('out', 'lambda-058.json')
dh_rect.set_params(T0_, Y0_, X_)
runs = {}
for f in sorted(glob.glob(os.path.join(HERE, RUNDIR, '*.json'))):
    d = json.load(open(f)); runs[os.path.basename(f)] = d


def cover(intervals, a, b):
    """True iff the union of closed intervals contains [a, b]."""
    x = a
    for lo, hi in sorted(intervals):
        if lo > x + 1e-12: return False
        x = max(x, hi)
        if x >= b: return True
    return x >= b


edges = {'bottom': [], 'top': [], 'left': []}
bad = []
for name, d in runs.items():
    if not d.get('ok'): bad.append(name); continue
    if d['part'] == 'rs':
        edges[d['edge']].append((d['xa'], d['xb']))
    else:
        edges[d['edge']].append((d['a'], d['b']))
T = time.time()
fa = dh_upper.final_region(t0=T0_, y0=Y0_, X=X_)
X, Y0, YT = dh_rect.XR, dh_rect.Y0, dh_rect.YTOP
checks = dict(
    runs_ok=not bad, bad_runs=bad,
    bottom_covered=cover(edges['bottom'], 0.0, X), top_covered=cover(edges['top'], 0.0, X),
    left_covered=cover(edges['left'], Y0, YT),
    right_side_part_a=bool(fa['ok']) and fa['X'] <= X and fa['y0'] <= Y0 and fa['Ytop'] >= YT - 1e-12,
    same_t0=all(abs(d['t0'] - dh_rect.T0) < 1e-15 for d in runs.values()) and abs(fa['t0'] - dh_rect.T0) < 1e-15,
    same_rect=all(abs(d.get('y0', Y0) - Y0) < 1e-15 and abs(d.get('X', X) - X) < 1e-9 for d in runs.values()),
)
ok = all(v for k, v in checks.items() if k != 'bad_runs')
out = dict(
    claim='H_%g of the Davenport-Heilbronn function has no zeros with Im z >= %g; hence Lambda_DH <= %g' % (T0_, Y0_, T0_ + Y0_ ** 2 / 2),
    t0=dh_rect.T0, y0=Y0, Ytop=YT, X=X, Delta=dh_rect.DELTA,
    argument=['zeros of H_t0 lie in Im^2 <= Delta^2 - 2 t0 (de Bruijn; Delta from claim #2)',
              'x >= X: final region (part a): |G - 1| <= 1 - FR_min + err < 1',
              'R = [0, X] x [y0, Ytop]: G = H/B_t analytic, B_t nonvanishing; G avoids (-inf, 0] on dR (part b), '
              'so winding 0 and no zeros in R',
              'H even and real: zeros symmetric under z -> -conj z and z -> conj z',
              'so every zero of H_t0 has |Im z| < y0; de Bruijn: Lambda <= t0 + y0^2/2'],
    Lambda_bound=dh_rect.T0 + Y0 ** 2 / 2, checks=checks, ok=ok,
    status=('conditional: the computation is Arb-rigorous, but the effective-approximation error terms rest on '
            'lemmas not yet written out in full (papers/supporting/riemann-siegel-identity/derivation.md, open items O1-O4)'),
    part_a={k: fa[k] for k in ('FR_min', 'tail_min', 'eC', 'eAB', 'err_total', 'Kh', 'blocks', 'ok', 'secs')},
    part_b={n: {k: d.get(k) for k in ('part', 'edge', 'xa', 'xb', 'a', 'b', 'y', 'pieces', 'evals', 'worst_margin',
                                         'err_max', 'Kh', 'ok', 'secs') if k in d} for n, d in runs.items()},
    secs_part_a_rerun=round(time.time() - T, 1),
)
dst = os.path.join(HERE, '..', '..', 'certificates', '03-davenport-heilbronn', OUT)
json.dump(out, open(dst, 'w'), indent=1, default=str)
print(json.dumps(dict(ok=ok, checks=checks, Lambda_bound=out['Lambda_bound'], part_a=out['part_a']), default=str))
