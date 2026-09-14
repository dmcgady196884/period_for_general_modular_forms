"""Item D, answered: hat r_{f_z} as an explicit power series in 1/j(z).

f_z = g j'/(j - j(z)) depends on z ONLY through x = j(z).  On the arc |j(tau)| <= j(i) = 1728,
so for every z with |j(z)| > 1728 the geometric expansion converges uniformly on the contour:

    j'/(j - x) = -(j'/x) sum_{m>=0} (j/x)^m   ==>   hat r_{f_z} = sum_{m>=1} x^{-m} v_m ,
    v_m = -hat r_{g j^{m-1} j'} .

With g = Delta E_4 and j' = -2 pi i E_4^2 E_6/Delta,

    g j' = -2 pi i E_4^3 E_6 = -2 pi i * (the reference form),

so v_1 = 0 IDENTICALLY: hat r annihilates the reference form because Phi(E_4^3 E_6) = 1.  That,
and not any CM arithmetic, is why cm_vs_generic_k18.py measured |hat r_{f_z}| ~ e^{-4 pi Im z}
(i.e. |q|^2) instead of |q|.  Note this is REFERENCE-DEPENDENT: replacing E_4^3 E_6 by
E_4^3 E_6 + c Delta E_6 (equally admissible, Phi = 1) restores a nonzero v_1 proportional to
r(Delta E_6) and the law reverts to |q|.  So this family singles out ref = g j'/c_0(g j').

PREDICTIONS, stated before the run:
  (1) v_m = 2 pi i hat r_{h_m}, h_m = E_4^3 E_6 j^{m-1} in M_18^! -- each lies in W.
  (2) the truncation residual at order M falls by a factor 1728/|x| per term, since 1728 is
      where the geometric series is slowest on the contour.  Concretely, for the M = 2
      truncation: d=19 (|x| = 884736) ~ 2e-3, d=11 (32768) ~ 5e-2, d=7 (3375) ~ 5e-1.
      These match the crude scale test already done by hand on |c_0 x^2|: 0.2%, 4%, 31%.
  (3) d=7 sits at |x|/1728 = 1.95, barely inside the disc of convergence, so its series decays
      by only a factor two per term and five terms will NOT suffice there.  That is a feature:
      it shows the expansion is the real mechanism rather than a fit.
  (4) the v_m are nearly parallel -- this is what makes the eight hat r_{f_z} sit on one line
      in P(W) to 7.7e-7.  Same hierarchy rank_of_rSkbang.py found at k=12 (rank 3 at 1e-18,
      2 at 1e-12, 1 at 1e-8).  Pairwise angles are printed.

Everything here is holomorphic on H (poles only at the cusp), so plain arcint applies and
TRAP 4 does not arise.  depth 10 / dps 30 is justified: direction_floor.py showed a 4x mesh
refinement moves the answer by 1.1e-28.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, E4, E6, Delta, jay, ktil, arcint, inW)      # noqa: E402

mp.mp.dps = 30
K, N, DEPTH, MMAX = 18, 16, 10, 5
g = lambda t: Delta(t) * E4(t)
ref = lambda t: E4(t)**3 * E6(t)


def rv(h):
    out = []
    for l in range(N + 1):
        s = mp.mpf(l + 1)
        L = mp.e**(-I * pi * s / 2) * arcint(
            lambda t: h(t) * (t**(s - 1) + ktil(t, s, K)), DEPTH)
        out.append((2 * pi * I)**(N + 1) * (-1)**l * mp.binomial(N, l) * I**(l + 1) * L)
    return out


def nrm(v):
    return max(abs(x) for x in v)


rE = rv(ref)
PhE = arcint(ref, DEPTH)
print("Phi(E_4^3 E_6) = %s   (target 1)" % mp.nstr(PhE, 14), flush=True)


def hatr(h):
    r = rv(h)
    P = arcint(h, DEPTH)
    return [a - P * b for a, b in zip(r, rE)], P


print("\nthe fixed vectors v_m = 2 pi i hat r_{E_4^3 E_6 j^{m-1}}", flush=True)
V = {}
for m in range(1, MMAX + 1):
    h = (lambda t, m=m: ref(t) * jay(t)**(m - 1))
    hr, P = hatr(h)
    V[m] = [2 * pi * I * x for x in hr]
    s_, u_ = inW(V[m], K) if nrm(V[m]) > 0 else (mp.mpf(0), mp.mpf(0))
    print("   m=%d  Phi = %-24s |v_m| = %-14s  inW %s / %s"
          % (m, mp.nstr(P, 10), mp.nstr(nrm(V[m]), 8), mp.nstr(s_, 4), mp.nstr(u_, 4)),
          flush=True)
print("   |v_1|/|v_2| = %s   (prediction: zero)"
      % mp.nstr(nrm(V[1]) / nrm(V[2]), 6), flush=True)

print("\npairwise directions of v_2..v_5 (1 - |<a,b>|/(|a||b|), Hermitian)", flush=True)
for a in range(2, MMAX + 1):
    for b in range(a + 1, MMAX + 1):
        ip = abs(sum(mp.conj(V[a][l]) * V[b][l] for l in range(N + 1)))
        na = mp.sqrt(sum(abs(V[a][l])**2 for l in range(N + 1)))
        nb = mp.sqrt(sum(abs(V[b][l])**2 for l in range(N + 1)))
        print("   v_%d . v_%d :  1 - cos = %s" % (a, b, mp.nstr(1 - ip / (na * nb), 6)),
              flush=True)

ZS = [("CM d=7 ", (1 + I * mp.sqrt(7)) / 2), ("CM d=11", (1 + I * mp.sqrt(11)) / 2),
      ("CM d=19", (1 + I * mp.sqrt(19)) / 2)]
print("\ntruncation test:  hat r_{f_z} =? sum_{m=2}^{M} v_m x^{-m}", flush=True)
for lbl, z in ZS:
    x = jay(z)
    f = lambda t: g(t) * (-2 * I * pi * jay(t) * E6(t) / E4(t)) / (jay(t) - x)
    v, _ = hatr(f)
    print("   %s  x = j(z) = %-22s |x|/1728 = %s"
          % (lbl, mp.nstr(x, 10), mp.nstr(abs(x) / 1728, 7)), flush=True)
    for M in range(2, MMAX + 1):
        approx = [sum(V[m][l] / x**m for m in range(2, M + 1)) for l in range(N + 1)]
        res = nrm([a - b for a, b in zip(v, approx)]) / nrm(v)
        print("        M=%d  relative residual = %s" % (M, mp.nstr(res, 6)), flush=True)
    print("        predicted per-term decay 1728/|x| = %s"
          % mp.nstr(1728 / abs(x), 6), flush=True)
