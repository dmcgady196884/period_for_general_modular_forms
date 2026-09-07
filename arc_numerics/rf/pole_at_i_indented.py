"""SUPERSEDED by i_limit.py.  Kept because its FAILURE is the instructive part.

The defect this script measures is NOT a property of f.  Indenting past i on one side by
a fixed h passes around BOTH members of the pinching stabiliser pair {z, Sz} as soon as
h > delta -- a different homotopy class, and not the continuation of the arc.  That is
the whole source of the factor-3 asymmetry and |(1+U+U^2)| = 2.67 reported below.
Deform and rescale instead (i_limit.py): no obstruction at any pole order.

Its one sound result is the k=18 control, where the average lands in W to 2e-30 / 2e-26.
That works because a SIMPLE pole at i is not a merger -- there is nothing to pinch.

Original docstring follows.
---------------------------------------------------------------------------
Poles at i, done with INDENTED contours instead of symmetric excision.

Supersedes the excision analysis in pole_at_i.py, which measured an artifact.

The point.  L^* is homotopy-invariant (prop:indepHomotopy), and a pole sitting ON the
contour is resolved by choosing a side.  BOTH one-sided contours are finite at every
pole order -- the contour simply does not touch the pole -- and by lem:wall they differ
by a full residue.  So no divergence here can be intrinsic.

What symmetric excision does wrong.  Deleting |theta - pi/2| < delta and letting
delta -> 0 is NOT the average of the two indented contours once the pole is worse than
simple.  Indenting adds a small semicircle whose own contribution diverges like
1/delta and exactly cancels the divergence of the excised piece; excision drops that
counterterm.  For a SIMPLE pole the semicircle contributes only the (finite) half
residue and the two prescriptions agree -- which is why the k=18 control passed and
concealed the error.

So the object of interest is
    hat r_f := (r_f^{above} + r_f^{below}) / 2 ,
finite at every order, and S-fixed because S exchanges the two sides.

Contours.  Excise |theta - pi/2| < delta from gamma^arc and rejoin the endpoints
e^{i(pi/2 + delta)} -> e^{i(pi/2 - delta)} by a polyline through i(1 +- h).  "Above"
(+h) is the |tau| > 1 side.  Homotopy invariance says each side's answer must not
depend on delta or h, which is the built-in self-check below.

Cases (both at dim S_k = 1, dim W = 3):
    k = 18   E4^6/E6      simple pole   -- must reproduce the PV answer
    k = 20   E4^11/E6^2   double pole   -- the 4 | k case that excision mangled
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common import mp, I, pi, E4, E6, ktil, check_orientation, inW

from validate import path_int

mp.mp.dps = 30

print("Phi(E_12) =", mp.nstr(check_orientation(), 12), " (guard passed)\n", flush=True)


def _nodes_toward(a, b, depth):
    out, d = [a], b - a
    for _ in range(depth):
        d /= 2
        out.append(b - d)
    return sorted(set(out + [b]))


def arc_indented(g, delta, h, side, depth=8):
    """gamma^arc with |theta - pi/2| < delta replaced by a detour through i(1 + side*h).

    Orientation rho -> rho+1, i.e. theta decreasing, so the detour runs from
    e^{i(pi/2 + delta)} to e^{i(pi/2 - delta)}.
    """
    tot = mp.mpc(0)
    for a, b, toward in ((pi / 3, pi / 2 - delta, True),
                         (pi / 2 + delta, 2 * pi / 3, False)):
        th = (_nodes_toward(a, b, depth) if toward
              else [a + b - x for x in reversed(_nodes_toward(a, b, depth))])
        th = sorted(set(th))
        tot += sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
                   for u, v in zip(th[:-1], th[1:]))
    arcpart = -tot
    p1 = mp.e**(I * (pi / 2 + delta))
    p2 = mp.e**(I * (pi / 2 - delta))
    apex = I * (1 + side * h)
    detour = path_int(g, [p1, apex, p2], nsub=24)
    return arcpart + detour


def rvec_indented(g, k, delta, h, side, depth=8):
    n = k - 2
    out = []
    for l in range(n + 1):
        s = mp.mpf(l + 1)
        L = (arc_indented(lambda t: g(t) * t**(s - 1), delta, h, side, depth)
             + arc_indented(lambda t: g(t) * ktil(t, s, k), delta, h, side, depth))
        out.append((2 * pi * I)**(n + 1) * (-1)**l * mp.binomial(n, l) * L)
    return out


CASES = [
    ("k=18  E4^6/E6     SIMPLE", 18, lambda t: E4(t)**6 / E6(t),
     lambda t: E4(t)**3 * E6(t)),
    ("k=20  E4^11/E6^2  DOUBLE", 20, lambda t: E4(t)**11 / E6(t)**2,
     lambda t: E4(t)**5),
]
if len(sys.argv) > 1:
    want = {int(a) for a in sys.argv[1:]}
    CASES = [c for c in CASES if c[1] in want]

GEOM = [(mp.mpf('0.10'), mp.mpf('0.10')),
        (mp.mpf('0.05'), mp.mpf('0.05')),
        (mp.mpf('0.05'), mp.mpf('0.12'))]     # same class, different shape

for name, k, f, ref in CASES:
    print("=" * 74)
    print(name)
    print("=" * 74)
    avgs = []
    for delta, h in GEOM:
        row = {}
        for side, tag in ((+1, "above"), (-1, "below")):
            rE = rvec_indented(ref, k, delta, h, side)
            rf = rvec_indented(f, k, delta, h, side)
            Phi = arc_indented(f, delta, h, side)
            row[tag] = [a - Phi * b for a, b in zip(rf, rE)]
        avg = [(a + b) / 2 for a, b in zip(row["above"], row["below"])]
        avgs.append(avg)
        sa, ua = inW(row["above"], k)
        sb, ub = inW(row["below"], k)
        sm, um = inW(avg, k)
        print("  delta=%s h=%s" % (mp.nstr(delta, 3), mp.nstr(h, 3)))
        print("     above   |tilde r| = %-16s in W: %s / %s"
              % (mp.nstr(max(abs(x) for x in row["above"]), 9),
                 mp.nstr(sa, 4), mp.nstr(ua, 4)))
        print("     below   |tilde r| = %-16s in W: %s / %s"
              % (mp.nstr(max(abs(x) for x in row["below"]), 9),
                 mp.nstr(sb, 4), mp.nstr(ub, 4)))
        print("     AVERAGE |tilde r| = %-16s in W: %s / %s"
              % (mp.nstr(max(abs(x) for x in avg), 9),
                 mp.nstr(sm, 4), mp.nstr(um, 4)), flush=True)
    print("  homotopy self-check -- the average must not depend on (delta, h):")
    base = avgs[0]
    nb = max(abs(x) for x in base)
    for m in range(1, len(avgs)):
        print("     geom %d vs geom 0:  max|d| / |avg| = %s"
              % (m, mp.nstr(max(abs(a - b) for a, b in zip(avgs[m], base)) / nb, 6)))
    print()
