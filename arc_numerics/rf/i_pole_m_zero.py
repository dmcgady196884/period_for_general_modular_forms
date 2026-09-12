"""At i, is the W-selecting configuration m = 0, i.e. the MEAN of the two detours?

At rho the (1+U+U^2) defect was measured to be exactly D = -A m, where m is the connector's
winding about rho's image DOWNSTAIRS in SL2(Z)\\H (one upstairs turn = n = 3 downstairs loops,
the valence formula's 1/3).  Four windings gave |A|, 2|A|, 5|A|, 4|A| at m = -1, 2, 5, -4:
exact.  The endpoints S tau_0 and T^{-1} tau_0 are one stabiliser step apart, so paths between
them realise only m = 2 mod 3 and m = 0 is unreachable -- hence no contour is in W there.

At i the same count predicts m in 1 + 2 Z (odd), so again never 0, and m = 0 is the MIDPOINT of
m = +-1 -- the mean of the inside and outside detours.  detour_vs_retract.py already found that
mean to be finite, exactly radius-independent, and exactly REAL.  This tests whether it is also
the W-selecting configuration, which is the claim that unifies i (1/2) and rho (1/3) under
"evaluate the affine-in-m family at m = 0", i.e. add 1/n of the elliptic residue.

Setup is deliberately the UNSHIFTED base point tau_0 = rho+1, both segments = gamma^arc: the pole
at i is INTERIOR to the arc, so no base-point shift and no connector are needed.  The two
detours of radius r about i are m = -1 and m = +1; their mean is m = 0.

PREDICTIONS.
 (1) Each single side FAILS, and fails on (1+S): S acts by (r,theta) -> (1/r, pi-theta) and so
     carries an inside detour at i to an outside one, meaning neither side is S-symmetric.  This
     is the one place (1+S) discriminates -- at rho it never did, since the spiral was
     S-symmetric by construction.
 (2) The MEAN passes BOTH relations, ~1e-20 or better.  S swaps the two detours, so the mean is
     the unique S-invariant combination.
 (3) Each side's defect is r-independent, and the two are equal and opposite about the mean
     (D = -A m at m = -+1), so |D_in| = |D_out|.
 (4) Phi differs between sides by 2 pi i Res_i f, and the mean's Phi is their average.

f = Delta/(j-1728) = Delta^2/E_6^2, k = 12: E_6 has a simple zero at i, so P = 2.  Radii 0.1 and
0.05 check r-independence; the nearest other point of SL2(Z).i is i +- 1, at distance 1.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, Delta, jay, ktil, E12,                   # noqa: E402
                    rvec, inW, NU, check_orientation)

mp.mp.dps = 30
K = 12
N = K - 2


def arcseg(g, th0, th1, n=20):
    xs = [th0 + (th1 - th0) * mp.mpf(j) / n for j in range(n + 1)]
    return sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
               for u, v in zip(xs[:-1], xs[1:]))


def detour_int(g, r, side, n=24):
    """gamma^arc from rho to rho+1 (theta DECREASING), detoured about i at radius r.

    |e^{i th} - i| = r at th = pi/2 +- al, al = 2 arcsin(r/2), and
        e^{i(pi/2+al)} - i = -2 sin(al/2) e^{i al/2}   (arg pi + al/2)
        e^{i(pi/2-al)} - i =  2 sin(al/2) e^{-i al/2}  (arg  - al/2)
    so psi runs pi+al/2 -> -al/2 through pi/2 (above i, |tau|>1, OUTSIDE) or
    pi+al/2 -> 2pi-al/2 through 3pi/2 (below i, |tau|<1, INSIDE).
    """
    al = 2 * mp.asin(r / 2)
    psi0 = pi + al / 2
    psi1 = -al / 2 if side == 'out' else 2 * pi - al / 2
    ps = [psi0 + (psi1 - psi0) * mp.mpf(j) / n for j in range(n + 1)]
    circ = sum(mp.quad(lambda p: g(I + r * mp.e**(I * p)) * I * r * mp.e**(I * p), [u, v])
               for u, v in zip(ps[:-1], ps[1:]))
    return arcseg(g, 2 * pi / 3, pi / 2 + al) + circ + arcseg(g, pi / 2 - al, pi / 3)


def rhat(g, r, side, rE):
    rf = []
    for l in range(N + 1):
        s = mp.mpf(l + 1)
        L = detour_int(lambda t: g(t) * t**(s - 1), r, side) + \
            detour_int(lambda t: g(t) * ktil(t, s, K), r, side)
        rf.append((2 * pi * I)**(N + 1) * (-1)**l * mp.binomial(N, l) * L)
    Ph = detour_int(g, r, side)
    return [a - Ph * b for a, b in zip(rf, rE)], Ph


print("guard: Phi(E_12) plain arc = %s" % mp.nstr(check_orientation(), 12), flush=True)
f = lambda t: Delta(t) / (jay(t) - 1728)
rE = rvec(E12, K)

for r in (mp.mpf('0.1'), mp.mpf('0.05')):
    print("=" * 76, flush=True)
    print("radius r = %s" % mp.nstr(r, 4), flush=True)
    got = {}
    for side in ('in', 'out'):
        v, Ph = rhat(f, r, side, rE)
        got[side] = (v, Ph)
        s_, u_ = inW(v, K)
        print("   %-4s (m=%+d)  Phi = %-24s |(1+S)| = %-12s |(1+U+U^2)| = %-12s |defect| = %s"
              % (side, -1 if side == 'in' else 1, mp.nstr(Ph, 10), mp.nstr(s_, 6),
                 mp.nstr(u_, 6), mp.nstr(max(abs(x) for x in NU(v, K)), 10)), flush=True)
    mean = [(a + b) / 2 for a, b in zip(got['in'][0], got['out'][0])]
    Phm = (got['in'][1] + got['out'][1]) / 2
    s_, u_ = inW(mean, K)
    print("   MEAN (m= 0)  Phi = %-24s |(1+S)| = %-12s |(1+U+U^2)| = %-12s |defect| = %s"
          % (mp.nstr(Phm, 10), mp.nstr(s_, 6), mp.nstr(u_, 6),
             mp.nstr(max(abs(x) for x in NU(mean, K)), 10)), flush=True)
    print("   Phi(in) - Phi(out) = %s" % mp.nstr(got['in'][1] - got['out'][1], 12), flush=True)
