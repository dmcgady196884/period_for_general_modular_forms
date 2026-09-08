"""Minimal reproduction: arc_pole_side_test.py's inside-dip contour vs both branches.

branch_showdown.py says eq:geoG needs arg(-2 pi n) = -pi, and that at z = e^{1.4i},
s = 2.4+0.6i the two branches differ by 2.339 relative.  arc_pole_side_test.py used +pi
at that same z and s and reported |inside - F| / |out| = 2.19e-54.  Both cannot hold.
Raw numbers here, no ratios, so whichever is wrong is visible.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, ktil                                      # noqa: E402

mp.mp.dps = 40
NMAX = 120
KK = 12
S = mp.mpc('2.4', '0.6')
AMP = mp.mpf('0.15')
WID = mp.mpf('0.15')
RHO = mp.e**(2 * I * pi / 3)
TAU0 = mp.e**(I * pi / 3)
TH = mp.mpf('1.4')
Z = mp.e**(I * TH)


def block0(tau):
    w = mp.e**(2 * I * pi * (tau - Z))
    return w / (1 - w)


def bump(th):
    u = (th - TH) / WID
    return mp.mpf(0) if abs(u) >= 1 else mp.e * mp.e**(-1 / (1 - u**2))


def dbump(th):
    u = (th - TH) / WID
    if abs(u) >= 1:
        return mp.mpf(0)
    return bump(th) * (-2 * u / (1 - u**2)**2) / WID


def polar_int(g, sign, panels=48):
    def integrand(th):
        r = mp.e**(sign * AMP * bump(th))
        return g(r * mp.e**(I * th)) * r * (sign * AMP * dbump(th) + I) * mp.e**(I * th)
    a, b = pi / 3, 2 * pi / 3
    ns = sorted(set([a, b, TH, TH - WID, TH + WID]
                    + [a + (b - a) * m / panels for m in range(panels + 1)]))
    return -sum(mp.quad(integrand, [u, v]) for u, v in zip(ns[:-1], ns[1:]))


def chord_int(g, panels=60):
    A, B = RHO, TAU0
    return sum(mp.quad(lambda x: g(A + (B - A) * x) * (B - A),
                       [mp.mpf(m) / panels, mp.mpf(m + 1) / panels])
               for m in range(panels))


def Lstar_polar(sign):
    return mp.e**(-I * pi * S / 2) * (
        polar_int(lambda t: block0(t) * t**(S - 1), sign)
        + polar_int(lambda t: block0(t) * ktil(t, S, KK), sign))


def Lstar_chord():
    return mp.e**(-I * pi * S / 2) * (
        chord_int(lambda t: block0(t) * t**(S - 1))
        + chord_int(lambda t: block0(t) * ktil(t, S, KK)))


def Gtab(argm1):
    e2, Gs = mp.e**(2 * I * pi * S), mp.gamma(S)
    out = []
    for n in range(1, NMAX + 1):
        xr, xt = 2 * I * pi * n * RHO, 2 * I * pi * n * TAU0
        pre = mp.e**(-S * (mp.log(2 * pi * n) + I * pi / 2))
        Sn = pre * (e2 * mp.gammainc(S, xr, mp.inf) + (1 - e2) * Gs
                    - mp.gammainc(S, xt, mp.inf))
        base = mp.e**(mp.log(2 * pi * n) + I * argm1)
        Tn = (mp.gammainc(S, xt, mp.inf) / base**S
              + I**KK * mp.gammainc(KK - S, xt, mp.inf) / base**(KK - S))
        out.append(mp.e**(-I * pi * S / 2) * Sn + Tn)
    return out


def G0():
    return (mp.e**(-I * pi * S / 2) * (TAU0**S - RHO**S) / S
            - (TAU0 / I)**S / S - I**KK * (TAU0 / I)**(KK - S) / (KK - S))


def geoIN(G):
    return -G0() - sum(mp.e**(2 * I * pi * n * Z) * G[n - 1]
                       for n in range(1, NMAX + 1))


print("s = %s   z = e^{1.4 i}   N = 0" % mp.nstr(S, 8), flush=True)
print("%-6s %-6s %-16s %-16s %-16s %s"
      % ("dps", "NMAX", "|ins-chord|", "|ins-F(+pi)|", "|ins-F(-pi)|", "|F(+pi)-F(-pi)|"),
      flush=True)

for dps in (40, 60):
    for nmax in (120, 200):
        mp.mp.dps = dps
        globals()['NMAX'] = nmax
        ins = Lstar_polar(-1)
        cho = Lstar_chord()
        fp, fm = geoIN(Gtab(pi)), geoIN(Gtab(-pi))
        sc = abs(ins)
        print("%-6d %-6d %-16s %-16s %-16s %s"
              % (dps, nmax,
                 mp.nstr(abs(ins - cho) / sc, 6),
                 mp.nstr(abs(ins - fp) / sc, 6),
                 mp.nstr(abs(ins - fm) / sc, 6),
                 mp.nstr(abs(fp - fm) / sc, 6)), flush=True)
