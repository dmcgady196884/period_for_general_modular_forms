"""eq:geoIN alone at several precisions, no quadrature.  Seconds, not half an hour.

repro_arcside.py (dps 40, NMAX 120) gets F(-pi) = the contour value
    0.1916488053709276 + 0.2724964207185487 i
arc_pole_side_test.py (dps 60, NMAX 200) reports F(-pi) wrong by 0.5 at the same z, s, N
with what looks like identical code.  Since the formula needs no contour integral, just
evaluate it at both settings and compare against that known target.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi                                            # noqa: E402

KK = 12
TARGET = "0.1916488053709276 + 0.2724964207185487"


def F(dps, nmax, argm1):
    mp.mp.dps = dps
    S = mp.mpc('2.4', '0.6')
    RHO = mp.e**(2 * I * pi / 3)
    TAU0 = mp.e**(I * pi / 3)
    Z = mp.e**(I * mp.mpf('1.4'))
    e2, Gs = mp.e**(2 * I * pi * S), mp.gamma(S)
    acc = mp.mpc(0)
    for n in range(1, nmax + 1):
        xr, xt = 2 * I * pi * n * RHO, 2 * I * pi * n * TAU0
        pre = mp.e**(-S * (mp.log(2 * pi * n) + I * pi / 2))
        Sn = pre * (e2 * mp.gammainc(S, xr, mp.inf) + (1 - e2) * Gs
                    - mp.gammainc(S, xt, mp.inf))
        # (-2 pi n)^{-w} written as ONE explicit exponential.  Never build the negative
        # real first and then raise it to a complex power: mp.e**(L +- i pi) carries a
        # roundoff-sized imaginary part, and **w then takes mpmath's PRINCIPAL branch of
        # that, so the intended branch is decided by the sign of the roundoff and flips
        # with dps and n.  That is what made every earlier branch verdict noise.
        lg = mp.log(2 * pi * n) + I * argm1
        Tn = (mp.gammainc(S, xt, mp.inf) * mp.e**(-S * lg)
              + I**KK * mp.gammainc(KK - S, xt, mp.inf) * mp.e**(-(KK - S) * lg))
        acc += mp.e**(2 * I * pi * n * Z) * (mp.e**(-I * pi * S / 2) * Sn + Tn)
    G0 = (mp.e**(-I * pi * S / 2) * (TAU0**S - RHO**S) / S
          - (TAU0 / I)**S / S - I**KK * (TAU0 / I)**(KK - S) / (KK - S))
    return -G0 - acc


print("target (the contour value) = %s i" % TARGET, flush=True)
for dps in (30, 40, 50, 60, 80):
    for nmax in (120, 200):
        fm = F(dps, nmax, -pi)
        fp = F(dps, nmax, pi)
        print("  dps %-4d NMAX %-4d  F(-pi) = %-42s  F(+pi) = %s"
              % (dps, nmax, mp.nstr(fm, 16), mp.nstr(fp, 16)), flush=True)
