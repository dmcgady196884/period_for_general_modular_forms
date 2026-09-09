"""Hole 3: does def:arcsplit's limit depend on the deforming family?

def:arcsplit sets L*(f,s) := lim eta^kappa L*(f_eta,s) over SOME admissible family, kappa
the least exponent for which the limit exists.  def:ellfam picks one family; the question
is whether any other gives the same answer.

Given f_0 with pole order P at tau_e, the representation g/(j - j_e)^m with 0 <= a < n is
UNIQUE (m = ceil(P/n), a = mn - P, g = f_0 (j-j_e)^m), so def:ellfam is canonical -- but it
is still a choice.  Two admissible families with the SAME f_0:

    A:  f_eta = g / (j - j_e - eta)^m           (def:ellfam; one displaced orbit)
    B:  f_eta = g / ((j - j_e)^m - eta)         (m displaced orbits, at j - j_e = eta^{1/m})

Both are holomorphic in eta, both have f_0 = g/(j-j_e)^m, and for eta off the real axis
neither has a pole on SL2(Z).gamma^arc.

PREDICTION.  Near tau_e, A scales with c delta^n = eta, so delta ~ eta^{1/n}; B has
c^m delta^{mn} = eta, so delta ~ eta^{1/(mn)}.  The integrand carries delta^{-(P-1)} in both
cases, so

    kappa_A = (P-1)/n        while      kappa_B = (P-1)/(mn),

different whenever m > 1.  If that is what the scan shows, family-independence is FALSE at
the elliptic points and the honest fix is to name def:ellfam as the definition rather than
claim independence.

Run at the DISCRIM configuration, where A is already verified: tau_e = i, n = 2, g = E_4^3,
a = 0, m = 2, P = 4, k = 12.  kappa_A = 3/2 (measured, flat to 8 digits); kappa_B = 3/4.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, jay, ktil                             # noqa: E402

mp.mp.dps = 30
KK, NLOC, MM, JE = 12, 2, 2, mp.mpf(1728)
SVALS = [mp.mpf(3), mp.mpf(7)]
KAPPAS = [mp.mpf(r) / 4 for r in range(0, 9)]
ETAS = ['1e-2', '1e-3', '1e-4', '1e-5']


def arc_int(g, depth=22):
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


def fam_A(t, eta):
    return E4(t)**3 / (jay(t) - JE - eta)**MM


def fam_B(t, eta):
    return E4(t)**3 / ((jay(t) - JE)**MM - eta)


print("dps = %d   tau_e = i, n = %d, m = %d, P = %d, k = %d"
      % (mp.mp.dps, NLOC, MM, MM * NLOC, KK), flush=True)
print("predicted: kappa_A = %s, kappa_B = %s\n"
      % (mp.nstr(mp.mpf(MM * NLOC - 1) / NLOC, 4),
         mp.nstr(mp.mpf(MM * NLOC - 1) / (MM * NLOC), 4)), flush=True)

for name, fam in (("A  g/(j-j_e-eta)^m ", fam_A), ("B  g/((j-j_e)^m-eta)", fam_B)):
    for s in SVALS:
        print("%s\nfamily %s   s = %s\n%s"
              % ("=" * 78, name, mp.nstr(s, 4), "=" * 78), flush=True)
        print("  eta      " + "".join("k=%-11s" % mp.nstr(kp, 4) for kp in KAPPAS),
              flush=True)
        for es in ETAS:
            eta = I * mp.mpf(es)
            L = Lstar(lambda t, eta=eta: fam(t, eta), s, KK)
            row = "".join("%-13s" % mp.nstr(abs(eta**kp * L), 6) for kp in KAPPAS)
            print("  %-9s%s" % (es, row), flush=True)
