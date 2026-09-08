"""Which branch of (-2 pi n)^s does eq:geoG take?  Direct, one script, no inference.

geoG_branch_check.py (mode-by-mode, T-part alone) says arg(-2 pi n) = -pi.
arc_pole_side_test.py (assembled, arg = +pi) reproduces its contour to 1e-54.
Both cannot be right.  This compares the ASSEMBLED eq:geoIN under both branches against
quadrature, at one complex s, for a pole above the arc and a pole on it.

Note the two candidates differ by a constant factor e^{-2 pi i s} on the T-part alone,
of modulus e^{2 pi Im s} = 43.4 at s = 2.4 + 0.6i, so the wrong one cannot hide.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, ktil                                      # noqa: E402

mp.mp.dps = 40
NMAX = 120
KK = 12
S = mp.mpc('2.4', '0.6')
RHO = mp.e**(2 * I * pi / 3)
TAU0 = mp.e**(I * pi / 3)


def block0(tau, z):
    w = mp.e**(2 * I * pi * (tau - z))
    return w / (1 - w)


def chord_int(g, panels=60):
    A, B = RHO, TAU0
    return sum(mp.quad(lambda x: g(A + (B - A) * x) * (B - A),
                       [mp.mpf(m) / panels, mp.mpf(m + 1) / panels])
               for m in range(panels))


def arc_int(g, panels=60):
    a, b = pi / 3, 2 * pi / 3
    th = [a + (b - a) * m / panels for m in range(panels + 1)]
    return -sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
                for u, v in zip(th[:-1], th[1:]))


def Lstar(z, integrator):
    g = lambda t: block0(t, z)
    return mp.e**(-I * pi * S / 2) * (integrator(lambda t: g(t) * t**(S - 1))
                                      + integrator(lambda t: g(t) * ktil(t, S, KK)))


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


def G0(argm1):
    base = mp.e**(I * argm1 - I * pi / 2)       # (tau0/i) carries no (-1); kept explicit
    del base
    return (mp.e**(-I * pi * S / 2) * (TAU0**S - RHO**S) / S
            - (TAU0 / I)**S / S - I**KK * (TAU0 / I)**(KK - S) / (KK - S))


def geoIN(z, G, argm1):
    return -G0(argm1) - sum(mp.e**(2 * I * pi * n * z) * G[n - 1]
                            for n in range(1, NMAX + 1))


Gp, Gm = Gtab(pi), Gtab(-pi)
print("s = %s,  |e^{-2 pi i s}| = %s,  NMAX = %d, dps = %d"
      % (mp.nstr(S, 8), mp.nstr(abs(mp.e**(-2 * I * pi * S)), 8), NMAX, mp.mp.dps),
      flush=True)
print("%-22s %-14s %-14s %-14s %s"
      % ("z", "|F(+pi)-F(-pi)|", "arg=+pi", "arg=-pi", "verdict"), flush=True)

# NOTE: the chord is the right target for BOTH poles.  For z above the arc nothing is
# crossed, and for z ON the arc the inside-dipping contour deforms to the chord without
# crossing, so chord quadrature is the bare eq:geoIN value in either case.  Never call
# arc_int for a pole on |tau| = 1: that runs the contour through the pole.
for lbl, z in (("0.3+1.4i (above arc)", mp.mpc('0.3', '1.4')),
               ("e^{1.4i} (ON the arc)", mp.e**(I * mp.mpf('1.4')))):
    ch = Lstar(z, chord_int)
    fp = geoIN(z, Gp, pi)
    fm = geoIN(z, Gm, -pi)
    ep = abs(fp - ch) / abs(ch)
    em = abs(fm - ch) / abs(ch)
    print("%-22s %-14s %-14s %-14s %s"
          % (lbl, mp.nstr(abs(fp - fm) / abs(ch), 4), mp.nstr(ep, 4), mp.nstr(em, 4),
             "+pi" if ep < em else "-pi"), flush=True)
