"""Is the short/long connector ambiguity at rho exactly one wall-crossing?

basepoint_shift_rho.py finds L* finite and flat in eps at a pole on SL2(Z).rho, but the two
routings of the T-segment connector around rho give different (individually flat) values.  The
two paths differ by one loop about rho, so lem:wall_arc predicts

    short - long = 2 pi i e^{-i pi s/2} Res_{tau=rho}( f ktil )

with NO tau^{s-1} term: the connector lives on the T-segment only, Gamma_S being shared.  If it
matches, the ambiguity is an ordinary wall-crossing already described by eq:wall_arc, not a new
phenomenon.  Measured at s=7, k=4, f = E_4/j (double pole at rho):
    short - long = 0.00313285273687852494 + 0.00325281242127539154 j

Res is taken as a circle about rho, radius swept to show it is a genuine residue (r-independent)
and not a quadrature accident.  Note j has a TRIPLE zero at rho, so f = E_4/j has a pole of order
3m - a = 3 - 1 = 2 there; the circle must stay inside the disc where rho is the only pole, and the
nearest other point of SL2(Z).rho is rho+1 at distance 1.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, jay, ktil                             # noqa: E402

mp.mp.dps = 40
RHO = mp.e**(2 * I * pi / 3)
MEASURED = {7: mp.mpf('0.00313285273687852494') + I * mp.mpf('0.00325281242127539154')}


def residue(g, c, r, n=48):
    ps = [2 * pi * mp.mpf(j) / n for j in range(n + 1)]
    tot = sum(mp.quad(lambda p: g(c + r * mp.e**(I * p)) * I * r * mp.e**(I * p), [u, v])
              for u, v in zip(ps[:-1], ps[1:]))
    return tot / (2 * I * pi)


f = lambda t: E4(t) / jay(t)
for s in (mp.mpf(7), mp.mpf(11)):
    print("s = %s" % mp.nstr(s, 3), flush=True)
    for r in (mp.mpf('0.08'), mp.mpf('0.04'), mp.mpf('0.02')):
        R = residue(lambda t: f(t) * ktil(t, s, 4), RHO, r)
        pred = 2 * I * pi * mp.e**(-I * pi * s / 2) * R
        line = "   r=%-6s  2 pi i e^{-i pi s/2} Res = %s" % (mp.nstr(r, 3), mp.nstr(pred, 18))
        if int(s) in MEASURED:
            line += "   |diff| = %s" % mp.nstr(abs(pred - MEASURED[int(s)]), 6)
        print(line, flush=True)
