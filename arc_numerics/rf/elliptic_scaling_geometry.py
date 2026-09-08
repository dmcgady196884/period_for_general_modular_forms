"""Geometry of the rescaled contour C in the proof of lem:arcell.

The proof substitutes tau = tau_e + delta u with c delta^n = eta, where
    j(tau) - j(tau_e) = c (tau - tau_e)^n (1 + O(tau - tau_e)),  c != 0.
After rescaling the n nearby poles sit at u^n = 1, and the limit integral is over C, the
image of the contour.  The proof is unconditional ONLY IF C misses those poles.  Since
u-space is tau-space rotated by -arg(delta) = -(arg eta - arg c)/n, this is a statement
about the arc's tangent direction at tau_e measured against arg c.

What is computed here:
  c at each elliptic point, by Richardson-extrapolating (j - j(tau_e))/(tau - tau_e)^n
    along several directions of approach (the limit must be direction-independent, which
    is itself a check that the zero really has order exactly n);
  the arc's tangent direction at tau_e;
  the resulting direction of C in the u-plane, for the admissible eta of def:arcsplit;
  the angular separation between C and the nearest pole ray arg u = 2 pi l / n.

At i (n=2) the arc passes THROUGH, so C is a full line and the poles are at u = +-1.
At rho (n=3) the arc ENDS, so C is a half-line, and rho+1 supplies a second one; the two
are S-related.  Both are reported.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, jay

mp.mp.dps = 40

def leading_coeff(tau_e, jval, n, dirs, hs=('1e-6','1e-7','1e-8','1e-9')):
    """(j - jval)/(tau - tau_e)^n as tau -> tau_e along several directions"""
    out = {}
    for d in dirs:
        vals = []
        for hs_ in hs:
            h = mp.mpf(hs_)
            t = tau_e + h*d
            vals.append((jay(t) - jval)/(h*d)**n)
        out[d] = vals
    return out

print("="*76)
print("tau_e = i,  n = 2,  j(i) = 1728")
print("="*76)
dirs2 = [mp.mpf(1), I, mp.e**(I*pi/4)]
res = leading_coeff(I, mp.mpf(1728), 2, dirs2)
for d, vals in res.items():
    print("  approach along %-22s : %s" % (mp.nstr(d,6), ", ".join(mp.nstr(v,10) for v in vals)))
c_i = res[dirs2[0]][-1]
print("  c_i  =", mp.nstr(c_i, 14), "   arg c_i =", mp.nstr(mp.arg(c_i), 10))
print("  |c_i| =", mp.nstr(abs(c_i), 12))

print()
print("="*76)
print("tau_e = rho = e^{2 pi i/3},  n = 3,  j(rho) = 0")
print("="*76)
rho = mp.e**(2*I*pi/3)
dirs3 = [mp.mpf(1), I, mp.e**(I*pi/4)]
res3 = leading_coeff(rho, mp.mpf(0), 3, dirs3)
for d, vals in res3.items():
    print("  approach along %-22s : %s" % (mp.nstr(d,6), ", ".join(mp.nstr(v,10) for v in vals)))
c_r = res3[dirs3[0]][-1]
print("  c_rho  =", mp.nstr(c_r, 14), "   arg c_rho =", mp.nstr(mp.arg(c_r), 10))

print()
print("="*76)
print("tangent of the arc, and the direction of C")
print("="*76)
# arc: tau(theta) = e^{i theta}, dtau/dtheta = i e^{i theta}
for name, tau_e, n, cval, theta_e, adm in (
        ("i",     I,     2, c_i, pi/2,     "eta > 0  (arg eta = 0)"),
        ("rho",   rho,   3, c_r, 2*pi/3,   "arg eta = pi  (eta < 0; eta > 0 puts poles ON the arc)"),
        ("rho+1", mp.e**(I*pi/3), 3, None, pi/3, "same eta as at rho")):
    tang = I*mp.e**(I*theta_e)
    print("  %-6s tangent dtau/dtheta = %s   arg = %s"
          % (name, mp.nstr(tang, 10), mp.nstr(mp.arg(tang), 8)))

print()
for name, n, cval, theta_e, argeta in (("i", 2, c_i, pi/2, mp.mpf(0)),
                                       ("rho", 3, c_r, 2*pi/3, pi)):
    argdelta = (argeta - mp.arg(cval))/n
    tang = I*mp.e**(I*theta_e)
    argC = mp.arg(tang) - argdelta
    # reduce mod 2pi/n and compare against pole rays arg u = 2 pi l / n
    sep = min(abs(mp.fmod(argC - 2*pi*l/n + pi, 2*pi) - pi) for l in range(n))
    print("  %-5s n=%d  arg c = %-12s arg delta = %-12s arg C = %-12s"
          % (name, n, mp.nstr(mp.arg(cval),8), mp.nstr(argdelta,8), mp.nstr(argC,8)))
    print("        angular separation of C from the nearest pole ray arg u = 2 pi l/n : %s rad"
          % mp.nstr(sep, 8))
    print("        (0 would mean C runs straight into a pole; pi/n = %s is maximal)"
          % mp.nstr(pi/n, 8))

# ---------------------------------------------------------------------------
# CORRECTED: the arc is traversed with theta DECREASING (2pi/3 -> pi/3).  At i the
# contour passes THROUGH, so C is a full line and the sign is irrelevant.  At rho and
# rho+1 the contour ENDS, so C is a RAY and the sign matters: the ray leaves rho in the
# direction -dtau/dtheta, and at rho+1 it points back along +dtau/dtheta.
print()
print("="*76)
print("CORRECTED: C as a ray, with the arc traversed theta-decreasing")
print("="*76)
rho = mp.e**(2*I*pi/3)
def sep_from_poles(argC, n):
    return min(abs(mp.fmod(argC - 2*pi*l/n + pi, 2*pi) - pi) for l in range(n))

for name, theta_e, n, cval, argeta, sgn, kind in (
        ("i",     pi/2,   2, c_i, mp.mpf(0), +1, "full line"),
        ("rho",   2*pi/3, 3, c_r, pi,        -1, "ray leaving rho"),
        ("rho+1", pi/3,   3, c_r, pi,        +1, "ray leaving rho+1 back along the arc")):
    argdelta = (argeta - mp.arg(cval))/n
    ray = sgn * I*mp.e**(I*theta_e)
    argC = mp.arg(ray) - argdelta
    s = sep_from_poles(argC, n)
    print("  %-6s (%s)" % (name, kind))
    print("     arg(ray in tau) = %-12s arg delta = %-12s arg C = %-12s"
          % (mp.nstr(mp.arg(ray),8), mp.nstr(argdelta,8), mp.nstr(argC,8)))
    print("     separation from nearest pole ray = %-14s   maximum pi/n = %s"
          % (mp.nstr(s,10), mp.nstr(pi/n,10)))
    print("     ratio to maximum = %s" % mp.nstr(s/(pi/n), 10))
