"""(numeric) The initial heat-flow drift predicted by Lemma L of proof-unbounded.md, against the tracked drift.

For the deep zero near the two-term balance sigma0 = log2|a(2)/a(1)|, Step 4 of the proof predicts
|d Im z/dt| = Re phi'(s1) - (log 2)/2 at t = 0, with s1 = sigma0 + i Im(rho) and
phi'(s) = (1/2) psi((s+1)/2) + (1/2) log(7/pi).
The tracked values come from runs/x7_track_deep_*.txt (x7_track.py): the change of Im z over the first time step range.
Usage: python drift_check.py
"""
import mpmath as mp

mp.mp.dps = 20
# (theta label, sigma0 from runs/deep_zeros_*.txt, Im rho, Im z at t = 0, t1, Im z at t1) from runs/x7_track_deep_*.txt
DATA = [('704/1000', '12.862998', '4.54922536305453443', '-24.7348271551807644', 0.5, '-24.22814'),
        ('70406/10^5', '18.540798', '4.53405811825604376', '-36.0825132705847808', 1.0, '-34.90077'),
        ('a(1) = 1e-10', '32.713497', '4.53236556875250106', '-64.4269975424084191', 1.0, '-62.97095')]
for lab, s0, im, y0, t1, y1 in DATA:
    s1 = mp.mpc(s0, im)
    dphi = mp.digamma((s1 + 1) / 2) / 2 + mp.log(mp.mpf(7) / mp.pi) / 2
    pred = mp.re(dphi) - mp.log(2) / 2
    tracked = (mp.mpf(y1) - mp.mpf(y0)) / t1
    print('%-13s sigma0 = %s  Re phi\'(s1) = %s  predicted drift = %s  tracked = %s  ratio = %s'
          % (lab, s0, mp.nstr(mp.re(dphi), 6), mp.nstr(pred, 5), mp.nstr(tracked, 5), mp.nstr(tracked / pred, 5)))
