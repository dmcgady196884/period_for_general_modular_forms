"""Which branch of (-2 pi n)^s does the T-part of eq:geoG take?

lens_residue_check.py reproduces the chord to 1e-100 at s = 3 but misses by 0.74 at
s = 2.4 + 0.6i, with the crossing-term columns unaffected.  That points at a branch
choice invisible at integer s.  The S-part of eq:geoG carries
    e^{-i pi s/2} (2 pi i n)^{-s} = (2 pi n)^{-s} e^{-i pi s},
i.e. arg(-1) = -pi, while the T-part is written with (-2 pi n)^{-s} whose branch is set
by lem:Tclosed-Mk.  The two candidates differ by e^{2 pi i s}.

GROUND TRUTH.  L^*(q^{-n}, s) in the geodesic class, by quadrature on the chord:
    e^{-i pi s/2} int_rho^{tau0} e^{-2 pi i n tau} tau^{s-1} dtau        (S-segment)
  + e^{-i pi s/2} int_rho^{tau0} e^{-2 pi i n tau} ktil(tau,s) dtau      (T-segment)
Both segments are the chord; q^{-n} is entire, so no residue enters and the chord IS the
arc value here.  This is exactly G^geo_{s,k}(-n) of eq:geoG.

Reported per (s, n): the S-part alone against its quadrature, then the T-part with
arg(-1) = +pi and with arg(-1) = -pi against its quadrature.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, ktil                                      # noqa: E402

mp.mp.dps = 40
KK = 12
RHO = mp.e**(2 * I * pi / 3)
TAU0 = mp.e**(I * pi / 3)
SVALS = [mp.mpf(3), mp.mpf(7), mp.mpc('2.4', '0.6'), mp.mpc('-0.3', '2.2')]


def chord(g, panels=40):
    A, B = RHO, TAU0
    return sum(mp.quad(lambda x: g(A + (B - A) * x) * (B - A),
                       [mp.mpf(m) / panels, mp.mpf(m + 1) / panels])
               for m in range(panels))


def S_quad(s, n):
    return mp.e**(-I * pi * s / 2) * chord(
        lambda t: mp.e**(-2 * I * pi * n * t) * t**(s - 1))


def T_quad(s, n):
    return mp.e**(-I * pi * s / 2) * chord(
        lambda t: mp.e**(-2 * I * pi * n * t) * ktil(t, s, KK))


def S_form(s, n):
    """eq:geoSupper, with the eq:sheet factor on the rho endpoint"""
    e2 = mp.e**(2 * I * pi * s)
    pre = mp.e**(-s * (mp.log(2 * pi * n) + I * pi / 2))      # (2 pi i n)^{-s}
    return mp.e**(-I * pi * s / 2) * pre * (
        e2 * mp.gammainc(s, 2 * I * pi * n * RHO, mp.inf)
        + (1 - e2) * mp.gamma(s)
        - mp.gammainc(s, 2 * I * pi * n * TAU0, mp.inf))


def T_form(s, n, argm1):
    """last two terms of eq:geoG, with arg(-1) = argm1 in (-2 pi n)^s"""
    xt = 2 * I * pi * n * TAU0
    base = mp.e**(mp.log(2 * pi * n) + I * argm1)
    return (mp.gammainc(s, xt, mp.inf) / base**s
            + I**KK * mp.gammainc(KK - s, xt, mp.inf) / base**(KK - s))


print("dps = %d, k = %d" % (mp.mp.dps, KK), flush=True)
for s in SVALS:
    print("\ns = %s" % mp.nstr(s, 8), flush=True)
    for n in (1, 2, 3):
        sq, tq = S_quad(s, n), T_quad(s, n)
        sf = S_form(s, n)
        tp = T_form(s, n, pi)
        tm = T_form(s, n, -pi)
        print("  n=%d  S: %-11s   T(arg=+pi): %-11s  T(arg=-pi): %s"
              % (n,
                 mp.nstr(abs(sf - sq) / abs(sq), 4),
                 mp.nstr(abs(tp - tq) / abs(tq), 4),
                 mp.nstr(abs(tm - tq) / abs(tq), 4)), flush=True)
