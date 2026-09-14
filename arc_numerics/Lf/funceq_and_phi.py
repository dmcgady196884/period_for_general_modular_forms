"""Items B and C of the open-questions list: the functional equation, and Phi on M_k.

B.  L*(f,s) = i^k L*(f,k-s) for ALL f in F_k in the geodesic class.  This should be a lemma,
not an experiment.  Substituting tau -> -1/tau with f|_k S = f,
    int_{gamma^S} f tau^{s-1} dtau = (-1)^{s-1} int_{S gamma^S} f tau^{k-s-1} dtau,
and gamma^arc is S-stable with reversed orientation, so L*_S(f,s) = i^k L*_S(f,k-s).  The
T-segment needs NO path property: for even k the kernel obeys
    ktil(tau, k-s) = -e^{i pi (k-s-1)} ktil(tau, s),
giving L*_T(f,s) = i^{-k} L*_T(f,k-s), and i^k = i^{-k} = (-1)^{k/2} for even k.  So both
halves carry the same factor and the full L* satisfies the equation.  Checked here on
MEROMORPHIC f, where it has never been tested; the suite only covers f in M_k^!.

C.  Phi(g) = c_g(0) for HOLOMORPHIC g, since "no pole above the contour in the strip" is
vacuous.  DAM's observation: Phi(E_16) = Phi(E_4^4) = Phi(E_6^2 E_4) = 1, all three having
constant term 1 -- so dim S_k > 0 does not make the NORMALISATION ambiguous, only the
reference form.  The admissible reference forms are exactly {g in M_k : c_g(0) = 1}, and two
choices shift hat r_f by Phi(f) r_{g-g'} with g-g' in S_k.  Tested at k = 16 and k = 28.

PREDICTIONS.
 (1) functional equation at ~quadrature precision for every meromorphic f tried, at generic
     complex s as well as real s (a real-s-only check could hide an i^k vs i^{-k} error,
     since both are +-1 there).
 (2) Phi = 1 for every holomorphic g with c_g(0) = 1; Phi = 0 for every cusp form.
 (3) the SAME at k = 28, where dim S_28 = 2, so two independent cusp directions vanish.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, E6, Delta, jay, ktil, arcint, check_orientation  # noqa

mp.mp.dps = 30


def Lstar(f, k, s, depth=12):
    r"""def:Lint INCLUDING the e^{-i pi s/2} prefactor.  Omitting it compares the raw
    Lambda = I + J, which obeys Lambda(s) = -e^{i pi (s-1)} Lambda(k-s) instead -- at s=3,
    k=12 that factor is -1 against the expected +1, giving exactly 2.0 for every f."""
    g = lambda t: f(t) * (t**(s - 1) + ktil(t, s, k))
    return mp.e**(-I * pi * s / 2) * arcint(g, depth)


print("guard: Phi(E_12) = %s" % mp.nstr(check_orientation(), 12), flush=True)

print("=" * 78, flush=True)
print("B.  functional equation  L*(f,s) = i^k L*(f,k-s)", flush=True)
CASES = [
    ("k=12  Delta/(j-500)      poles ON the arc", lambda t: Delta(t) / (jay(t) - 500), 12),
    ("k=12  Delta/(j-j(2i))    pole above the arc",
     lambda t: Delta(t) / (jay(t) - jay(2 * I)), 12),
    ("k=12  E_4^3              holomorphic control", lambda t: E4(t)**3, 12),
    ("k=16  E_4^4/(j-500)      other weight", lambda t: E4(t)**4 / (jay(t) - 500), 16),
]
for lbl, f, k in CASES:
    for s in (mp.mpf(3), mp.mpf(7), mp.mpf('2.4') + mp.mpf('0.6') * I):
        a, b = Lstar(f, k, s), Lstar(f, k, k - s)
        ik = I**k
        num = abs(a - ik * b)
        den = max(abs(a), mp.mpf('1e-300'))
        print("   %-38s s=%-12s |L*(s) - i^k L*(k-s)|/|L*(s)| = %s"
              % (lbl, mp.nstr(s, 4), mp.nstr(num / den, 6)), flush=True)

print("=" * 78, flush=True)
print("C.  Phi(g) = c_g(0) on M_k;  admissible reference forms are c_g(0) = 1", flush=True)
HOL = [
    ("k=16  E_4^4        c(0)=1", lambda t: E4(t)**4, 1),
    ("k=16  E_6^2 E_4    c(0)=1", lambda t: E6(t)**2 * E4(t), 1),
    ("k=16  Delta E_4    c(0)=0  (cusp)", lambda t: Delta(t) * E4(t), 0),
    ("k=28  E_4^7        c(0)=1", lambda t: E4(t)**7, 1),
    ("k=28  E_4^4 E_6^2  c(0)=1", lambda t: E4(t)**4 * E6(t)**2, 1),
    ("k=28  E_4 E_6^4    c(0)=1", lambda t: E4(t) * E6(t)**4, 1),
    ("k=28  Delta E_4^4  c(0)=0  (cusp)", lambda t: Delta(t) * E4(t)**4, 0),
    ("k=28  Delta^2 E_4  c(0)=0  (cusp)", lambda t: Delta(t)**2 * E4(t), 0),
]
for lbl, g, want in HOL:
    v = arcint(g, 12)
    print("   %-34s Phi = %-26s  target %d   err %s"
          % (lbl, mp.nstr(v, 14), want, mp.nstr(abs(v - want), 6)), flush=True)

print("=" * 78, flush=True)
print("   the reference-form ambiguity is Phi(f) r(S_k): it VANISHES when Phi(f) = 0.",
      flush=True)
for lbl, f in [("k=12  Delta/(j-j(2i))   pole ABOVE the arc, c_f(0)=0",
                lambda t: Delta(t) / (jay(t) - jay(2 * I))),
               ("k=12  Delta/(j-500)     poles ON the arc, c_f(0)=0",
                lambda t: Delta(t) / (jay(t) - 500)),
               ("k=12  E_4^3             holomorphic, c(0)=1", lambda t: E4(t)**3)]:
    print("   %-46s Phi(f) = %s" % (lbl, mp.nstr(arcint(f, 12), 12)), flush=True)
print("   (Delta/j is omitted: its poles ARE rho and rho+1, the arc's own endpoints, so the", flush=True)
print("    plain arcint integrates through them and returns garbage -- 1.1e62.  That case", flush=True)
print("    needs the shifted contour of def:arcshift, read at W = 0.)", flush=True)
