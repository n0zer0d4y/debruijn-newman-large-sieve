"""float feasibility (numeric-experiment): which (t0, y0, X) could push Lambda_DH below 0.68?
For each (t0, y0) and X, the plain final-region lower bound LB (dh_final_region_scan.plain_lb) at the worst of
a few x in [X, 3X] and y in {y0, Ytop}, plus the float RS-error size at X (newman_chi.Chi.eC with the grid K_h).
Prints Lambda = t0 + y0^2/2 for the feasible cells (LB - eC > 0.05).
"""
import sys, math
import os as _os
sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'shared-dirichlet-pipeline'))
import numpy as np
from dh_final_region_scan import plain_lb, dh

DELTA = 1.7902723164702195
KH = 1.1793613456444694          # kh_rig for q = 5 (dh_rect run)
rows = []
for t0 in (0.35, 0.4, 0.45, 0.5, 0.55):
    Yt = math.sqrt(DELTA ** 2 - 2 * t0)
    for y0 in (0.2, 0.25, 0.3, 0.35, 0.4):
        for X in (2e4, 5e4, 1e5, 2e5, 5e5, 1e6):
            lbs = []
            for x in (X, 1.3 * X, 1.7 * X, 2.5 * X):
                for y in (y0, Yt):
                    lbs.append(plain_lb(x, t0, y)[0])
            lb = min(lbs)
            ec = dh.eC(X, y0, t0, KH)
            if lb - ec > 0.05:
                rows.append((t0 + y0 * y0 / 2, t0, y0, X, round(lb, 3), round(ec, 4)))
                break
rows.sort()
for r in rows[:15]:
    print('Lambda %.4f  t0 %.2f y0 %.2f  X %.0e  LB %.3f  eC %.4f' % r, flush=True)
