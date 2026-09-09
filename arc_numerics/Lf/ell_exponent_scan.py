"""What is the true exponent kappa of def:arcsplit at an elliptic point on the arc?

lem:ellLstar claims kappa = (P-1)/n.  model_integral.py showed that at tau_e = i (n = 2,
interior of the arc, so C is a full LINE) the model integral carries a factor
(1 + (-1)^a) and therefore VANISHES for odd a -- equivalently odd P, equivalently
k = 2 mod 4.  That does not say the contribution vanishes; it says the rescaling
eta^{(P-1)/2} is too aggressive, so the true kappa is SMALLER and def:arcsplit's "least
exponent for which the limit exists" is not (P-1)/n there.

Measured here: |eta^kappa L*(f_eta,s)| down a sequence of eta, for several kappa.  The
correct kappa is the column that settles to a finite nonzero constant; too-large kappa
drives it to 0, too-small kappa blows up.

  TEST   g = E_6, m = 2.  ord_i E_6 = 1, so a = 1 (ODD), P = 2m - a = 3, weight k = 6,
         and 6 = 2 mod 4, matching lem:ellLaurent's parity for an odd-order pole at i.
         lem:ellLstar predicts kappa = 1; the parity result says that is too big.

  CONTROL  g = E_4^3, m = 1.  E_4(i) != 0 so a = 0 (EVEN), P = 2, weight k = 12, 4 | 12.
         Here the model integral does not vanish and kappa = (P-1)/n = 1/2 should hold.
         If the control does not flatten at 1/2 the method is wrong, not the lemma.

eta is taken purely imaginary so that 1728 + eta avoids j(gamma^arc) = [0,1728].
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, E6, Delta, jay, ktil                  # noqa: E402

mp.mp.dps = 30
SVALS = [mp.mpf(3), mp.mpf(7)]
ETAS = ['1e-2', '1e-3', '1e-4', '1e-5']


def kappa_grid(n):
    """multiples of 1/n bracketing the candidates, plus 0"""
    return [mp.mpf(r) / n for r in range(0, 3 * n + 1)]

# The CONTROL has P = 2, so its predicted (P-1)/n is ALSO 1/2 -- it validates the method
# but cannot tell (P-1)/n apart from a constant 1/2.  DISCRIM has even a at m = 2, where
# (P-1)/n = 3/2 while a constant-1/2 alternative would still say 1/2.
# (label, g, m, k, j(tau_e), n).  RHO sits at the arc's two ENDPOINTS, where j = 0 and
# the elliptic order is 3; the others sit at i, where j = 1728 and the order is 2.
CASES = [("TEST    g=E_6,   m=2, a=1, P=3, k=6 ", lambda t: E6(t), 2, 6, 1728, 2),
         ("CONTROL g=E_4^3, m=1, a=0, P=2, k=12", lambda t: E4(t)**3, 1, 12, 1728, 2),
         ("DISCRIM g=E_4^3, m=2, a=0, P=4, k=12", lambda t: E4(t)**3, 2, 12, 1728, 2),
         # RHOA has a = 0, where (P-1)/n and (P-1-a)/n coincide -- it checks the method at
         # n = 3 but does not discriminate.  RHOB has a = 1: (P-1)/n = 1/3 if the parity
         # shift is special to i, versus 0 if it also operates at rho.
         ("RHOA    g=Delta, m=1, a=0, P=3, k=12", lambda t: Delta(t), 1, 12, 0, 3),
         ("RHOB    g=E_4,   m=1, a=1, P=2, k=4 ", lambda t: E4(t), 1, 4, 0, 3),
         # RHOC has a = n-1 = 2 at rho, the rho analogue of odd a at i: there the relative
         # cube root is w^{a+1} = w^0 = 1 AND modularity forces k = 2 mod 6, which makes
         # K(rho) = K(tau_0), so the two rays cancel and the leading term dies.  Unified
         # rule: the leading model term vanishes iff a = n-1.  Predicts (P-2)/n = 2/3, not
         # the (P-1)/n = 1 that eq:ellkappa currently claims at rho.
         ("RHOC    g=E_4^2, m=2, a=2, P=4, k=8 ", lambda t: E4(t)**2, 2, 8, 0, 3)]
if len(sys.argv) > 1:
    CASES = [c for c in CASES if c[0].split()[0] == sys.argv[1]]


def arc_int(g, depth=22):
    """int over gamma^arc, rho -> rho+1 (theta DECREASING), clustered hard at theta=pi/2
    where the displaced poles crowd the contour, and at both endpoints."""
    a, b = pi / 3, 2 * pi / 3
    mid = (a + b) / 2
    ns = [a, b, mid]
    for base in (a, b, mid):
        d = (b - a) / 4
        for _ in range(depth):
            d /= 2
            ns += [base - d, base + d]
    ns = sorted(set(x for x in ns if a <= x <= b))
    return -sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
                for u, v in zip(ns[:-1], ns[1:]))


def Lstar(f, s, k):
    return mp.e**(-I * pi * s / 2) * (arc_int(lambda t: f(t) * t**(s - 1))
                                      + arc_int(lambda t: f(t) * ktil(t, s, k)))


print("dps = %d   eta purely imaginary" % mp.mp.dps, flush=True)
for lbl, g, m, k, je, n in CASES:
    KAPPAS = kappa_grid(n)
    for s in SVALS:
        print("\n%s\n%s   n = %d   s = %s\n%s"
              % ("=" * 78, lbl, n, mp.nstr(s, 6), "=" * 78), flush=True)
        print("  eta      " + "".join("k=%-12s" % mp.nstr(kp, 4) for kp in KAPPAS),
              flush=True)
        for es in ETAS:
            eta = I * mp.mpf(es)
            f = lambda t, eta=eta, g=g, m=m, je=je: g(t) / (jay(t) - je - eta)**m
            L = Lstar(f, s, k)
            row = "".join("%-14s" % mp.nstr(abs(eta**kp * L), 7) for kp in KAPPAS)
            print("  %-9s%s" % (es, row), flush=True)
