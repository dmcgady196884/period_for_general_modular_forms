"""The pedestrian on-arc case: does the shifted harness handle a pole ON gamma^arc?

shifted_harness_validation.py checks poles OFF the contour.  This is the other pedestrian case,
and the one lem:arcon actually governs: x = j(z) in (0,1728) has TWO preimages on the arc, z and
Sz with arg z + arg Sz = pi, so f = Delta/(j-x) has an S-pair of simple poles sitting on the
reference contour.

The shifted spiral resolves them for free.  Its radius is (1+eps)^{(pi/2 - theta)/(pi/6)}, so it
runs OUTSIDE |tau| = 1 for theta < pi/2 and INSIDE for theta > pi/2 -- and S acts by
(r,theta) -> (1/r, pi-theta), so the two members of the S-pair are passed on opposite sides in
exactly the S-symmetric way lem:arcon's (1+S) argument needs.  No indentation is chosen by hand;
it falls out of displacing the base point.  (The PLAIN harness cannot run these cases at all: the
arc goes straight through both poles.)

Cases, k=12:
  75 deg   x = j(e^{5 i pi/12}), DAM's angle, on the arc (the arc is 60..120 deg); S-partner at
           105 deg.
  x = 500  the configuration already measured by arc_eight_classes.py part C, whose two diagonal
           pairs gave |(1+U+U^2)| = 3.25e-19 and |hat r| = 312512.435.  Recorded here as a VALUE
           cross-check, with the caveat that arc_eight_classes.py may not share common.rvec's
           normalisation -- if the defect matches and the magnitude does not, suspect the
           convention, not the contour.

PREDICTIONS.
 (1) |(1+S)| ~ 0 and |(1+U+U^2)| ~ 0 for BOTH routings and both cases.  The connector sits at
     radius eps about rho, where f has NO pole, so it encloses nothing.
 (2) short - long identically 0, for the same reason (cf. the pole-free control in
     basepoint_shift_rho.py, which gave 0.0 + 0.0j exactly).
 (3) Phi = 0: f = Delta/(j-x) starts at q^1, and Phi(f) = c_f(0).
If (1) fails here, the harness is broken and the rho verdict is void.  If it holds here AND in
shifted_harness_validation.py, the harness is sound and the rho defects are real.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, Delta, jay, ktil, E12, RHO,              # noqa: E402
                    rvec, inW, check_orientation)

mp.mp.dps = 30
TAU0 = mp.e**(I * pi / 3)
EPS = mp.mpf('1e-3')
K = 12


def clustered(n, depth):
    xs = [mp.mpf(j) / n for j in range(n + 1)]
    d = mp.mpf(1) / n
    for _ in range(depth):
        d /= 2
        xs += [d, 1 - d]
    return sorted(set(xs))


def integrate(path, g, n=20, depth=13):
    ts = clustered(n, depth)
    return sum(mp.quad(lambda t: g(path(t)[0]) * path(t)[1], [u, v])
               for u, v in zip(ts[:-1], ts[1:]))


def spiral(eps):
    L = mp.log(1 + eps)

    def p(t):
        th = 2 * pi / 3 + t * (pi / 3 - 2 * pi / 3)
        tau = mp.e**(L * (pi / 2 - th) / (pi / 6) + I * th)
        return tau, tau * (-L / (pi / 6) + I) * (pi / 3 - 2 * pi / 3)
    return p


def connector(eps, longway):
    a = (TAU0 * (1 + eps) - 1) - RHO
    b = RHO / (1 + eps) - RHO
    la, lb = mp.log(a), mp.log(b)
    d = mp.im(lb - la)
    while d > pi:
        d -= 2 * pi
    while d <= -pi:
        d += 2 * pi
    if longway:
        d += 2 * pi if d < 0 else -2 * pi
    lb = mp.mpc(mp.re(lb), mp.im(la) + d)

    def p(t):
        w = la + t * (lb - la)
        return RHO + mp.e**w, mp.e**w * (lb - la)
    return p


def shifted_tilde_r(g, k, eps, longway, rE):
    S, C = spiral(eps), connector(eps, longway)
    n = k - 2
    rf = []
    for l in range(n + 1):
        s = mp.mpf(l + 1)
        LS = integrate(S, lambda t: g(t) * t**(s - 1))
        LT = integrate(C, lambda t: g(t) * ktil(t, s, k)) + \
            integrate(S, lambda t: g(t) * ktil(t, s, k))
        rf.append((2 * pi * I)**(n + 1) * (-1)**l * mp.binomial(n, l) * (LS + LT))
    Ph = integrate(C, g) + integrate(S, g)
    return [a - Ph * b for a, b in zip(rf, rE)], Ph


print("guard: Phi(E_12) plain arc = %s" % mp.nstr(check_orientation(), 12), flush=True)
rE = rvec(E12, K)

Z75 = mp.e**(5 * I * pi / 12)
print("75 deg: z = %s, |z| = %s, j(z) = %s"
      % (mp.nstr(Z75, 10), mp.nstr(abs(Z75), 6), mp.nstr(jay(Z75), 12)), flush=True)

CASES = [
    ("75 deg   x = j(e^{5 i pi/12}), S-partner at 105 deg", jay(Z75)),
    ("x = 500  arc_eight_classes part C: defect 3.25e-19, |hat r| 312512.435", mp.mpf(500)),
]

for lbl, x in CASES:
    print("=" * 76, flush=True)
    print("%s\n   x = %s" % (lbl, mp.nstr(x, 12)), flush=True)
    f = lambda t, x=x: Delta(t) / (jay(t) - x)
    store = {}
    for lw in (False, True):
        tr, Ph = shifted_tilde_r(f, K, EPS, lw, rE)
        s_, u_ = inW(tr, K)
        store["long" if lw else "short"] = tr
        print("   %-5s Phi = %-20s |(1+S)| = %-12s |(1+U+U^2)| = %-12s |hat r| = %s"
              % ("long" if lw else "short", mp.nstr(Ph, 8), mp.nstr(s_, 6), mp.nstr(u_, 6),
                 mp.nstr(max(abs(v) for v in tr), 12)), flush=True)
    d = max(abs(a - b) for a, b in zip(store['short'], store['long']))
    print("   short - long (absolute) = %s" % mp.nstr(d, 6), flush=True)
