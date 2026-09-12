"""Retraction prescription, with the quadrature actually resolving the pole.

retraction_check.py used common.arcint, whose mesh clusters only at the arc's ENDPOINTS.
The displaced poles sit in the MIDDLE of the arc at distance ~ |Im x| / |j'(z)| ~ |Im x|/5000
from it, so that run never resolved them and both signs returned the same smoothed value.
arc_eight_classes.py records the same trap: the mesh must be clustered at z and Sz below
eps/5000, and two clustering depths must be run so that unresolved quadrature shows up as a
discrepancy rather than as a plausible wrong number.

Here the theta-mesh is clustered geometrically at BOTH pole angles as well as the endpoints,
to depth DEPTH (spacing ~ 0.26/2^DEPTH), and every quantity is computed at two depths.  If the
two depths disagree, the row is not converged and must not be read.

Checked, for f = Delta/(j - x), k = 12:
  Q2  does lim_{delta->0} L*(f_delta, s) exist with NO rescaling, along the retraction?
  Q3  do the two one-sided j-shifts x0 +- i eps give the SAME L* and Phi, or different?
      arc_eight_classes.py's Proposition B says DIFFERENT (they are the two diagonal
      classes A and B); the unresolved run said the same.  This settles it.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, Delta, jay, ktil                          # noqa: E402

mp.mp.dps = 30
KK = 12
X0 = mp.mpf(500)
YW = mp.mpf('0.95')
SVALS = [mp.mpf(3), mp.mpf(7)]

zp = mp.findroot(lambda t: jay(t) - X0, mp.e**(I * mp.mpf('1.29')))
zp = zp / abs(zp)
zm = -1 / zp
PHIP, PHIM = mp.arg(zp), mp.arg(zm)
W = I * YW


def arc_cl(g, depth):
    """arc integral, rho -> rho+1, mesh clustered at the endpoints AND both pole angles"""
    a, b = pi / 3, 2 * pi / 3
    ns = [a, b, (a + b) / 2]
    for base in (a, b, PHIP, PHIM):
        d = (b - a) / 4
        for _ in range(depth):
            d /= 2
            ns += [base - d, base + d]
        ns.append(base)
    ns = sorted(set(x for x in ns if a <= x <= b))
    return -sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
                for u, v in zip(ns[:-1], ns[1:]))


def Lc(x, s, depth):
    f = lambda t: Delta(t) / (jay(t) - x)
    return mp.e**(-I * pi * s / 2) * (arc_cl(lambda t: f(t) * t**(s - 1), depth)
                                      + arc_cl(lambda t: f(t) * ktil(t, s, KK), depth))


def Phic(x, depth):
    return arc_cl(lambda t: Delta(t) / (jay(t) - x), depth)


print("dps=%d  x0=%s  poles at arg %s and %s"
      % (mp.mp.dps, mp.nstr(X0, 6), mp.nstr(PHIP, 8), mp.nstr(PHIM, 8)), flush=True)
print("|j'(z)| ~ %s, so a pole with Im x = e sits ~ e/|j'| off the arc\n"
      % mp.nstr(abs((jay(zp * mp.e**(I * mp.mpf('1e-8'))) - X0) / (zp * (mp.e**(I * mp.mpf('1e-8')) - 1))), 6),
      flush=True)

print("Q3: the two one-sided j-shifts, resolved.  Two depths per entry.", flush=True)
for es, dep in (('1e-4', 30), ('1e-5', 34)):
    eps = mp.mpf(es)
    for s in SVALS:
        a1 = Lc(X0 + I * eps, s, dep)
        a2 = Lc(X0 + I * eps, s, dep + 4)
        b1 = Lc(X0 - I * eps, s, dep)
        b2 = Lc(X0 - I * eps, s, dep + 4)
        conv = max(abs(a1 - a2), abs(b1 - b2)) / abs(a2)
        print("   eps=%-6s s=%-4s depth%d-vs-%d: %-11s   |L*_+ - L*_-|/|L*| = %s"
              % (es, mp.nstr(s, 3), dep, dep + 4, mp.nstr(conv, 4),
                 mp.nstr(abs(a2 - b2) / abs(a2), 6)), flush=True)
    p1, p2 = Phic(X0 + I * eps, dep + 4), Phic(X0 - I * eps, dep + 4)
    print("   eps=%-6s Phi_+ = %-30s Phi_- = %-30s rel diff = %s"
          % (es, mp.nstr(p1, 12), mp.nstr(p2, 12),
             mp.nstr(abs(p1 - p2) / abs(p1), 6)), flush=True)

print("\nQ2: retraction limit, resolved.", flush=True)
for s in SVALS:
    print("   s = %s" % mp.nstr(s, 4), flush=True)
    prev = None
    for ds, dep in (('1e-2', 22), ('1e-3', 26), ('1e-4', 30), ('1e-5', 34)):
        d = mp.mpf(ds)
        xd = jay(zp - d * (zp - W))
        v1 = Lc(xd, s, dep)
        v2 = Lc(xd, s, dep + 4)
        ch = "" if prev is None else "  change %s" % mp.nstr(abs(v2 - prev), 6)
        print("     delta=%-6s L* = %-34s depth-conv %-11s%s"
              % (ds, mp.nstr(v2, 14), mp.nstr(abs(v1 - v2) / abs(v2), 4), ch), flush=True)
        prev = v2
