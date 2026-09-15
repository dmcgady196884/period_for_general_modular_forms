"""The W^pm localisation of Xi for f_z = Delta E_4 J'/(J-J(z)) in F_18: the point [mu_-:mu_+].

From Wpm18_projection.py, Xi_{f;z} = g(z)[ mu_-(z) W_- + mu_+(z) W_+ + mu_0(z) p_0 ] with
g = Delta E_4 and, after dividing by (2 pi i)^{n+2},

    mu_-  supported on EVEN m, anti-palindromic (coeff at 16-m = -coeff at m),
    mu_+  supported on ODD  m, palindromic      (coeff at 16-m = +coeff at m),

both with rational coefficients in the E_18 reference.  The v_0 = p_0 direction is set aside here.

WHAT THE LOCALISATION IS.  The scalar g(z) and the overall normalisation of each basis vector both
cancel in the RATIO, so the position of Xi inside the W^pm plane is the single rational function

    R(z) := mu_+(z) / mu_-(z) ,

a point of P^1 depending on z alone -- no periods, no transcendence, no dependence on f beyond the
fact that its pole is simple.

TWO INVARIANCES, both predicted before computing:
  * eq:thetafe gives Xi^o(-1/z) = -z^{-n} Xi^o(z) coordinate-wise, so BOTH mu_pm scale by the same
    factor and R(-1/z) = R(z).  The W^pm-localisation depends only on the S-ORBIT of the pole.
  * mu_- is even and mu_+ odd, so R(-z) = -R(z).
Anti-palindromic of even degree forces mu_-(+-1) = 0, so R has poles at z = +-1; palindromic mu_+
of odd degree forces mu_+(-1) = 0 as well.  Roots and the reduced form are reported.

Pure arithmetic on exact rationals -- no quadrature, nothing reused but the coefficient table.
"""
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I                                                   # noqa: E402

mp.mp.dps = 30
N = 16

MU_M = {0: F(-12000, 43867), 2: F(16, 3), 4: F(-50, 3), 6: F(52, 3),
        10: F(-52, 3), 12: F(50, 3), 14: F(-16, 3), 16: F(12000, 43867)}
MU_P = {1: F(-16), 3: F(308, 3), 5: F(-182), 7: F(286, 3),
        9: F(286, 3), 11: F(-182), 13: F(308, 3), 15: F(-16)}


def ev(c, z):
    return sum(mp.mpf(v.numerator) / v.denominator * z**m for m, v in c.items())


def clear(c):
    """scale to coprime integers; returns (coeffs, the scale removed)."""
    from math import gcd, lcm
    L = 1
    for v in c.values():
        L = lcm(L, v.denominator)
    ints = {m: int(v * L) for m, v in c.items()}
    g = 0
    for v in ints.values():
        g = gcd(g, abs(v))
    return {m: v // g for m, v in ints.items()}, F(g, L)


print("PART 1  the two polynomials, cleared to coprime integers", flush=True)
for nm, c in (("mu_-", MU_M), ("mu_+", MU_P)):
    ints, sc = clear(c)
    terms = "  ".join("%+d z^%d" % (v, m) for m, v in sorted(ints.items()))
    print("   %s = (%s) * [ %s ]" % (nm, sc, terms), flush=True)

print("\nPART 2  symmetry checks", flush=True)
for nm, c, sgn in (("mu_-", MU_M, -1), ("mu_+", MU_P, +1)):
    ok = all(c.get(N - m, F(0)) == sgn * v for m, v in c.items())
    par = "even" if all(m % 2 == 0 for m in c) else "odd"
    print("   %s: %s in z, coeff(16-m) = %s coeff(m) -> %s"
          % (nm, par, "+" if sgn > 0 else "-", "yes" if ok else "NO"), flush=True)
for z in (mp.mpf(1), mp.mpf(-1)):
    print("   mu_-(%s) = %-12s   mu_+(%s) = %s"
          % (mp.nstr(z, 2), mp.nstr(ev(MU_M, z), 6), mp.nstr(z, 2), mp.nstr(ev(MU_P, z), 6)),
          flush=True)

print("\nPART 3  R(z) = mu_+/mu_-, and its invariances", flush=True)
ZS = [("0.31+1.27i", mp.mpf('0.31') + I * mp.mpf('1.27')),
      ("i sqrt2        (d=8)", I * mp.sqrt(2)),
      ("(1+i sqrt7)/2  (d=7)", (1 + I * mp.sqrt(7)) / 2),
      ("(1+i sqrt11)/2 (d=11)", (1 + I * mp.sqrt(11)) / 2),
      ("(1+i sqrt19)/2 (d=19)", (1 + I * mp.sqrt(19)) / 2),
      ("(1+i sqrt43)/2 (d=43)", (1 + I * mp.sqrt(43)) / 2),
      ("(1+i sqrt163)/2(d=163)", (1 + I * mp.sqrt(163)) / 2),
      ("i              (d=4)", I),
      ("e^{2 pi i/3}+1 (d=3)", mp.mpf('0.5') + I * mp.sqrt(3) / 2)]
print("   %-24s %-34s %-12s %s" % ("z", "R(z)", "R(-1/z)/R(z)", "R(-z)/R(z)"), flush=True)
for lbl, z in ZS:
    d = ev(MU_M, z)
    if abs(d) < mp.mpf('1e-25'):
        print("   %-24s mu_- vanishes" % lbl, flush=True)
        continue
    R = ev(MU_P, z) / d
    Rs = ev(MU_P, -1 / z) / ev(MU_M, -1 / z)
    Rn = ev(MU_P, -z) / ev(MU_M, -z)
    print("   %-24s %-34s %-12s %s"
          % (lbl, mp.nstr(R, 16), mp.nstr(Rs / R, 6), mp.nstr(Rn / R, 6)), flush=True)

print("\nPART 4  is R(z) algebraic at the CM points?  (R is rational in z, so this is automatic;",
      flush=True)
print("        printed to show the actual values rather than to test anything)", flush=True)
for lbl, z in ZS[1:7]:
    d = ev(MU_M, z)
    if abs(d) < mp.mpf('1e-25'):
        continue
    R = ev(MU_P, z) / d
    print("   %-24s R = %s" % (lbl, mp.nstr(R, 20)), flush=True)
