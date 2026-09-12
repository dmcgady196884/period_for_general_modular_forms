"""Does the contour integral through an elliptic pole actually blow up?  No -- move the CONTOUR.

The delta^{1-P} divergence measured by finite_part_scan.py is a property of retracting the
POLE, not of L*.  f is modular, so a pole cannot be moved alone: the stabiliser of an elliptic
point of order n acts as rotation by 2 pi / n about it, and the contour runs through the fixed
point, so every displacement sends one image to each side.  The contour is pinched and the
integral scales as epsilon^{a + 1 - 2m} = delta^{1-P}.

Section 4.1 never had this problem because it moved the CONTOUR: the reference path sits at
i(1+delta) - 1 -> i(1+delta), lifted off the elliptic point, with f untouched.  On the arc the
same move is a detour of radius r around i, inside (|tau| < 1) or outside (|tau| > 1).

Claims tested, on f = E_4^3/(j - 1728) = E_4^3 Delta / E_6^2, k = 12, double pole at i:
  (1) each detour is FINITE and independent of r  (homotopy invariance: no other pole of f is
      enclosed -- the nearest points of SL2(Z).i are i +- 1, at distance 1);
  (2) inside - outside = 2 pi i Res_{tau=i} f K, evaluated independently as a small circle.

Geometry.  |e^{i theta} - i| = r at theta = pi/2 +- alpha, alpha = 2 arcsin(r/2), and
    e^{i(pi/2 + alpha)} - i = -2 sin(alpha/2) e^{i alpha/2}      (argument pi + alpha/2)
    e^{i(pi/2 - alpha)} - i =  2 sin(alpha/2) e^{-i alpha/2}     (argument   - alpha/2)
so the detour runs psi: pi + alpha/2 -> -alpha/2 through psi = pi/2 (above i, |tau| > 1, OUTSIDE)
or psi: pi + alpha/2 -> 2 pi - alpha/2 through psi = 3 pi/2 (below i, |tau| < 1, INSIDE).
Traversal is rho -> rho+1, i.e. theta decreasing, so the detour is entered at pi/2 + alpha.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, E6, Delta, ktil                       # noqa: E402

mp.mp.dps = 30
KK = 12


def f(t):
    return E4(t)**3 * Delta(t) / E6(t)**2


def integrand(t, s):
    return f(t) * (t**(s - 1) + ktil(t, s, KK))


def arcpiece(th0, th1, s, n=24):
    """Integral along |tau| = 1 from angle th0 to th1, subdivided into n panels."""
    xs = [th0 + (th1 - th0) * mp.mpf(j) / n for j in range(n + 1)]
    return sum(mp.quad(lambda x: integrand(mp.e**(I * x), s) * I * mp.e**(I * x), [u, v])
               for u, v in zip(xs[:-1], xs[1:]))


def detour(r, s, side, n=24):
    """Full arc rho -> rho+1 with a radius-r circular detour around i on the given side."""
    al = 2 * mp.asin(r / 2)
    psi0 = pi + al / 2
    psi1 = -al / 2 if side == 'out' else 2 * pi - al / 2
    ps = [psi0 + (psi1 - psi0) * mp.mpf(j) / n for j in range(n + 1)]
    arc_in = arcpiece(2 * pi / 3, pi / 2 + al, s)
    arc_out = arcpiece(pi / 2 - al, pi / 3, s)
    circ = sum(mp.quad(lambda p: integrand(I + r * mp.e**(I * p), s) * I * r * mp.e**(I * p),
                       [u, v]) for u, v in zip(ps[:-1], ps[1:]))
    return mp.e**(-I * pi * s / 2) * (arc_in + circ + arc_out)


def residue(r, s, n=32):
    """Res_{tau=i} f K, as a CCW circle of radius r about i."""
    ps = [2 * pi * mp.mpf(j) / n for j in range(n + 1)]
    c = sum(mp.quad(lambda p: integrand(I + r * mp.e**(I * p), s) * I * r * mp.e**(I * p),
                    [u, v]) for u, v in zip(ps[:-1], ps[1:]))
    return c / (2 * I * pi)


for s in (mp.mpf(7), mp.mpf(11)):
    print("=" * 78, flush=True)
    print("s = %s   f = E_4^3/(j-1728), k = 12, double pole at i" % mp.nstr(s, 3), flush=True)
    vals = {}
    for side in ('in', 'out'):
        for r in (mp.mpf('0.12'), mp.mpf('0.06'), mp.mpf('0.02')):
            vals[(side, r)] = detour(r, s, side)
            print("  %-3s r=%-6s L* = %s" % (side, mp.nstr(r, 3),
                                             mp.nstr(vals[(side, r)], 20)), flush=True)
    rs = [mp.mpf('0.12'), mp.mpf('0.06'), mp.mpf('0.02')]
    for side in ('in', 'out'):
        sp = max(abs(vals[(side, a)] - vals[(side, b)]) for a in rs for b in rs)
        print("  spread over r, %-3s : %s" % (side, mp.nstr(sp, 6)), flush=True)
    d = vals[('in', rs[-1])] - vals[('out', rs[-1])]
    R = residue(mp.mpf('0.07'), s)
    print("  in - out            = %s" % mp.nstr(d, 20), flush=True)
    print("  2 pi i Res_i(f K)   = %s" % mp.nstr(2 * I * pi * R * mp.e**(-I * pi * s / 2), 20),
          flush=True)
    print("  difference          = %s" % mp.nstr(abs(d - 2 * I * pi * R * mp.e**(-I * pi * s / 2)), 6),
          flush=True)
