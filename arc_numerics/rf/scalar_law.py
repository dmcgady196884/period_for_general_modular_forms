"""Item D, sharpened: is hat r_{f_z} = C * lambda(z) * V for a FIXED V in W?

cm_vs_generic_k18.py put all eight hat r_{f_z} on one complex line in W to 7.7e-7, and the
row norms obey |hat r_{f_z}| ~ e^{-4 pi Im z} to 0.05% (slope -12.562 against 4 pi = 12.5664,
on a baseline Im z = 1.3 -> 2.179, with no power of y riding along: a y^18 would have moved
the slope to -1.97).  A holomorphic scalar with |.| ~ e^{-4 pi Im z} has q_z^2 leading, so the
candidates are products of two weight-(k-2)-ish cusp forms.

The test needs no PSLQ and no basis interpretation: if lambda(z) is one of the candidates,
then c_0 / lambda(z) is constant across all eight z in MODULUS AND PHASE.  c_0 is read from
the same fixed SVD basis for every row, so its phase is comparable between rows.

  g h    g = Delta E_4 in M_16 (the numerator of f_z), h = Delta E_6 spanning S_18
  g^2    weight 32
  D^2    Delta^2, weight 24
  g D    weight 28

If g h wins, the answer to "how do the W^pm coefficients of g j'/(j - j(z)) depend on z" is:
they do not, except through the single scalar g(z) h(z) -- with h the cusp form whose periods
span the cuspidal part of W.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, E4, E6, Delta                                    # noqa: E402

mp.mp.dps = 25

# (label, z, c_0) with c_0 transcribed from logs/cm18.txt -- 10 significant digits
ROWS = [
    ("CM  d=7",    (1 + I * mp.sqrt(7)) / 2,      -4.111937873e+13 - 2.924308739e+13j),
    ("CM  d=8",    I * mp.sqrt(2),                -1.297501082e+13 - 9.227538141e+12j),
    ("CM  d=11",   (1 + I * mp.sqrt(11)) / 2,     -6.053634848e+11 - 4.305203744e+11j),
    ("CM  d=19",   (1 + I * mp.sqrt(19)) / 2,     -866410352.9 - 616170964.3j),
    ("ctrl 1.3i",  mp.mpf('1.3') * I,             -5.438713663e+13 - 3.867896642e+13j),
    ("ctrl .5+1.4i", mp.mpf('0.5') + mp.mpf('1.4') * I,
     -1.55788669e+13 - 1.107931471e+13j),
    ("gen  .2+1.3i", mp.mpf('0.2') + mp.mpf('1.3') * I,
     6.692557152e+13 - 9.446932464e+11j),
    ("gen  .35+1.5i", mp.mpf('0.35') + mp.mpf('1.5') * I,
     -1.621309781e+12 + 5.184916551e+12j),
]

g = lambda t: Delta(t) * E4(t)
h = lambda t: Delta(t) * E6(t)

CAND = [("g h", lambda t: g(t) * h(t)),
        ("g^2", lambda t: g(t)**2),
        ("D^2", lambda t: Delta(t)**2),
        ("g D", lambda t: g(t) * Delta(t))]

for nm, lam in CAND:
    print("=" * 78, flush=True)
    print("lambda(z) = %s" % nm, flush=True)
    vals = []
    for lbl, z, c0 in ROWS:
        r = mp.mpc(c0) / lam(z)
        vals.append(r)
        print("   %-15s c_0/lambda = %-34s |.| = %-14s arg = %s"
              % (lbl, mp.nstr(r, 10), mp.nstr(abs(r), 9),
                 mp.nstr(mp.arg(r), 9)), flush=True)
    m = sum(vals) / len(vals)
    sp = max(abs(v - m) for v in vals) / abs(m)
    print("   ---> relative spread about the mean: %s" % mp.nstr(sp, 5), flush=True)
