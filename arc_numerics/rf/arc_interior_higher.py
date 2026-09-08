"""Arc-interior poles, continued: higher order, higher weight, and PV vs deformation.

Companion to arc_interior_poles.py, which covers the simple pole at k=12.  Run parts
selectively:  python arc_interior_higher.py 2        (double pole, k=12)
              python arc_interior_higher.py 24       (simple pole, k=24, dim W = 5)
              python arc_interior_higher.py pv       (PV vs deformation average, k=12)
with no argument running all three.  Part 24 is the slow one: n = 22 means 23
coefficients, so ~184 deformed-arc quadratures.

What each part is for.

 PART 2 (double pole, k=12, f = Delta/(j - j0)^2).  Off the elliptic points a deformed
   pole leaves the arc on ONE side and nothing is pinched, so def:arcsplit's kappa should
   be 0 at every pole order.  If |tilde r| stays bounded here, that claim survives at
   P = 2; if it grows, the footnote on def:arcsplit is wrong and kappa > 0 generically.

 PART 24 (simple pole, k=24, f = Delta^2/(j - j0)).  dim S_24 = 2 so dim W = 5.  Both
   deformations should still land in W, and the point is that they do so in a 5-dimensional
   space rather than a 3-dimensional one, where "lands in W" is a much weaker accident.

 PART PV (k=12).  Is the average of the two S-symmetric deformations equal to the average
   of the two deformation limits (j0 -> j0 +- i eps)?  If yes, the two prescriptions agree
   away from the elliptic points and def:arcsplit really is one object with the elliptic
   case as its degeneration.  If no, the generic case is not understood either.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common import (mp, I, pi, E4, E6, Delta, jay, ktil, check_orientation, inW, E12,
                    arcint, rvec)

mp.mp.dps = 25
J0 = mp.mpf(500)
ETA = mp.mpf('0.15')

WHICH = set(sys.argv[1:]) or {"2", "24", "pv"}


def deformed_arc(g, eta, panels=40):
    def integrand(th):
        u = th - pi / 2
        r = mp.e**(eta * mp.sin(6 * u))
        drdth = r * eta * 6 * mp.cos(6 * u)
        tau = r * mp.e**(I * th)
        return g(tau) * (drdth + I * r) * mp.e**(I * th)
    a, b = pi / 3, 2 * pi / 3
    th = [a + (b - a) * m / panels for m in range(panels + 1)]
    return -sum(mp.quad(integrand, [u, v]) for u, v in zip(th[:-1], th[1:]))


def rvec_def(g, eta, k, panels=40):
    n = k - 2
    out = []
    for l in range(n + 1):
        s = mp.mpf(l + 1)
        out.append((2 * pi * I)**(n + 1) * (-1)**l * mp.binomial(n, l)
                   * (deformed_arc(lambda t: g(t) * t**(s - 1), eta, panels)
                      + deformed_arc(lambda t: g(t) * ktil(t, s, k), eta, panels)))
    return out


def tilde_def(g, ref, eta, k, panels=40):
    Phi = deformed_arc(g, eta, panels)
    rR = rvec_def(ref, eta, k, panels)
    return [a - Phi * b for a, b in zip(rvec_def(g, eta, k, panels), rR)], Phi


def report(nm, t, ph, k):
    s_, u_ = inW(t, k)
    print("  %-8s Phi = %-26s |tilde r| = %-16s in W: %s / %s"
          % (nm, mp.nstr(ph, 10), mp.nstr(max(abs(x) for x in t), 9),
             mp.nstr(s_, 4), mp.nstr(u_, 4)), flush=True)


print("Phi(E_12) =", mp.nstr(check_orientation(), 10), flush=True)

if "2" in WHICH:
    print()
    print("=" * 74)
    print("PART 2: DOUBLE pole on the arc, k = 12, f = Delta/(j - j0)^2")
    print("=" * 74)
    f2 = lambda t: Delta(t) / (jay(t) - J0)**2
    tA, phA = tilde_def(f2, E12, ETA, 12)
    tB, phB = tilde_def(f2, E12, -ETA, 12)
    report("eta>0", tA, phA, 12)
    report("eta<0", tB, phB, 12)
    D = [a - b for a, b in zip(tA, tB)]
    print("  |difference| / |tilde r| = %s   difference in W: %s / %s"
          % (mp.nstr(max(abs(x) for x in D) / max(abs(x) for x in tA), 6),
             *[mp.nstr(x, 5) for x in inW(D, 12)]))
    print("  NB kappa = 0 survives iff |tilde r| is bounded, i.e. comparable to the")
    print("     simple-pole case rather than blowing up.")

if "24" in WHICH:
    print()
    print("=" * 74)
    print("PART 24: simple pole on the arc, k = 24, dim S_24 = 2 so dim W = 5")
    print("=" * 74)
    f24 = lambda t: Delta(t)**2 / (jay(t) - J0)
    ref24 = lambda t: E4(t)**6                    # constant term 1, so Phi = 1
    print("  Phi(E4^6) =", mp.nstr(deformed_arc(ref24, ETA), 10), " (should be +1)",
          flush=True)
    tA, phA = tilde_def(f24, ref24, ETA, 24, panels=28)
    tB, phB = tilde_def(f24, ref24, -ETA, 24, panels=28)
    report("eta>0", tA, phA, 24)
    report("eta<0", tB, phB, 24)
    D = [a - b for a, b in zip(tA, tB)]
    print("  |difference| / |tilde r| = %s   difference in W: %s / %s"
          % (mp.nstr(max(abs(x) for x in D) / max(abs(x) for x in tA), 6),
             *[mp.nstr(x, 5) for x in inW(D, 24)]))

if "pv" in WHICH:
    print()
    print("=" * 74)
    print("PART PV: average of the two deformations vs average of the two limits, k = 12")
    print("=" * 74)
    f = lambda t: Delta(t) / (jay(t) - J0)
    tA, _ = tilde_def(f, E12, ETA, 12)
    tB, _ = tilde_def(f, E12, -ETA, 12)
    PV = [(a + b) / 2 for a, b in zip(tA, tB)]
    print("  PV (average of the two S-symmetric deformations): in W: %s / %s"
          % tuple(mp.nstr(x, 5) for x in inW(PV, 12)), flush=True)
    rE = rvec(E12, 12)
    for e in ('1e-2', '1e-3', '1e-4'):
        lims = []
        for sgn in (1, -1):
            j0 = J0 + sgn * I * mp.mpf(e)
            g = lambda t, j0=j0: Delta(t) / (jay(t) - j0)
            Phi = arcint(g)
            lims.append([a - Phi * b for a, b in zip(rvec(g, 12), rE)])
        avg = [(a + b) / 2 for a, b in zip(*lims)]
        d = max(abs(a - b) for a, b in zip(avg, PV)) / max(abs(x) for x in PV)
        print("  |eps|=%-5s  max|avg(limits) - PV| / |PV| = %s   avg in W: %s / %s"
              % (e, mp.nstr(d, 6), *[mp.nstr(x, 4) for x in inW(avg, 12)]), flush=True)
