"""The delta-retraction prescription: does the limit exist, and which diagonal is it?

DAM's proposal replaces contour indentation by a canonical deformation of f.  For a pole at
z = e^{i phi} on the arc with pi/3 <= phi <= pi/2, and w = i y with sqrt3/2 < y < 1, set

    z_w(delta) := z - delta (z - w),

a straight ray from z (delta = 0) into the interior of the unit disk.  Modularity carries the
partner to S(z_w(delta)), which goes OUTSIDE the disk; the retracted pole itself lands in the
lens.  Restricting phi to [pi/3, pi/2] is not a convention -- the other half of the arc is
forced by S.

The family is realised simply, since moving the pole is moving its j-value:

    f_delta := g / ( j - x(delta) )^m ,     x(delta) := j( z_w(delta) ),

so the whole prescription is a j-shift with a SPECIFIED direction.  Two questions:

  Q1  Which diagonal?  Equivalent to the sign of Im x(delta).  Note the retraction moves the
      POLE, while def:georef's "+i0" rule moved the CONTOUR, and a pole approaching from one
      side leaves the contour indented AWAY from it -- so the two rules may well name
      OPPOSITE classes.  Measured here rather than argued.
  Q2  Does lim_{delta->0} L*(f_delta,s) exist with no rescaling at a non-elliptic on-arc pole?
      prop:arcside says yes for the j-shift family; this checks it for the retraction.

Reported: sign of Im x(delta), convergence of L* along the retraction, and the two one-sided
j-shift limits x0 +- i eps for comparison.  Poles are off the arc for delta > 0, so plain
arcint applies throughout and no indented contour is ever needed.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, Delta, jay, ktil, arcint                   # noqa: E402

mp.mp.dps = 25
KK = 12
X0 = mp.mpf(500)
YW = mp.mpf('0.95')
SVALS = [mp.mpf(3), mp.mpf(7)]
DELTAS = ['1e-2', '1e-3', '1e-4', '1e-5', '1e-6']


def L(x, s):
    f = lambda t: Delta(t) / (jay(t) - x)
    return mp.e**(-I * pi * s / 2) * (arcint(lambda t: f(t) * t**(s - 1))
                                      + arcint(lambda t: f(t) * ktil(t, s, KK)))


# the pole of Delta/(j - X0) on the arc with pi/3 <= phi <= pi/2
zp = mp.findroot(lambda t: jay(t) - X0, mp.e**(I * mp.mpf('1.29')))
zp = zp / abs(zp)
w = I * YW
print("dps=%d  x0=%s  w=%s" % (mp.mp.dps, mp.nstr(X0, 6), mp.nstr(w, 6)), flush=True)
print("pole z+ = %s,  arg = %s  (pi/3 = %s, pi/2 = %s)"
      % (mp.nstr(zp, 10), mp.nstr(mp.arg(zp), 8),
         mp.nstr(pi / 3, 6), mp.nstr(pi / 2, 6)), flush=True)
print("partner S z+ = %s,  arg = %s\n"
      % (mp.nstr(-1 / zp, 10), mp.nstr(mp.arg(-1 / zp), 8)), flush=True)

print("Q1: where does the retraction send j?", flush=True)
for ds in DELTAS:
    d = mp.mpf(ds)
    zd = zp - d * (zp - w)
    xd = jay(zd)
    print("   delta=%-7s |z_w| = %-14s  x(delta) = %-30s  Im x = %s"
          % (ds, mp.nstr(abs(zd), 10), mp.nstr(xd, 12),
             "POSITIVE" if xd.imag > 0 else "NEGATIVE"), flush=True)

print("\nQ2: convergence of L* along the retraction, vs the two one-sided j-shifts",
      flush=True)
for s in SVALS:
    print("  s = %s" % mp.nstr(s, 4), flush=True)
    prev = None
    for ds in DELTAS:
        d = mp.mpf(ds)
        xd = jay(zp - d * (zp - w))
        v = L(xd, s)
        ch = "" if prev is None else "  change %s" % mp.nstr(abs(v - prev), 6)
        print("    retract delta=%-7s L* = %-34s%s" % (ds, mp.nstr(v, 14), ch),
              flush=True)
        prev = v
    for sgn, lbl in ((1, "x0 + i eps"), (-1, "x0 - i eps")):
        vv = L(X0 + sgn * I * mp.mpf('1e-6'), s)
        print("    %-12s eps=1e-6   L* = %-34s  |diff from retraction| = %s"
              % (lbl, mp.nstr(vv, 14), mp.nstr(abs(vv - prev), 6)), flush=True)
