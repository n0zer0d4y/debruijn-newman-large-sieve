# (numeric-experiment) Davenport-Heilbronn: see code/03-davenport-heilbronn/notes.md.
import subprocess, json
import sys as _sys
PY = _sys.executable   # any Python with python-flint, mpmath, numpy
import os as _os2
HERE = _os2.path.dirname(_os2.path.abspath(__file__))
for args in (['0.22', '480.6959411390449', '-0.14620587779148794', 'r=4e-5', 'delta=2e-5'],
             ['0.23', '480.69056', '-0.03486', 'r=1e-5', 'delta=5e-6']):
    r = subprocess.run([PY, 'cert_dh2.py'] + args, cwd=HERE, capture_output=True, text=True, timeout=130)
    print('ARGS', args); print(r.stdout.strip()); print(r.stderr.strip()[-400:])
