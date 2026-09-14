"""How big is r(S_k^!) inside W?  The reference-form ambiguity, measured rather than inferred.

The admissible reference forms in def:rfhat are {g : Phi(g) = 1}, and Phi(g) = c_g(0) for any g
with no poles in H -- weakly holomorphic included, since int_{tau_0-1}^{tau_0} q^n dtau = 0 for
negative n as well.  So the admissible set is {g in M_k^! : c_g(0) = 1} and the ambiguity in
hat r_f is Phi(f) * r(S_k^!), NOT Phi(f) * r(S_k) as a restriction to holomorphic g would give.

DAM's argument for why that is still not all of W: modulo the Bol image D^{k-1}(M^!_{2-k}), which
r annihilates, S_k^! is spanned by Delta_k and its Guerzhoy companion hat Delta_k, so
    dim r(S_k^!) = 2 dim S_k,
against dim W = 2 dim S_k + 1.  Codimension one, the missing direction being the Eichler-Shimura
C summand, i.e. p_0.  Consequence: hat r_f is canonical modulo r(S_k^!), with its p_0 component
the one surviving invariant -- and fully canonical when Phi(f) = 0, which is exactly the class
F_k^circ of lem:rfWmero.

THIS TEST converts that from a literature inference into a measurement.  At k = 12, dim W = 3, so
the prediction is rank 2, with every p_0 coordinate vanishing.

CONSTRUCTION of elements of S_12^!.  Any weight-12 weakly holomorphic form minus c_0 times E_4^3
lies in S_12^! (E_4^3 = Delta j has c_0 = 1).  And c_0 is READ OFF by arcint, since Phi(g) =
c_g(0) for g holomorphic on H -- no q-expansion arithmetic needed.  Take Delta itself, then
Delta j^m - c_0(Delta j^m) E_4^3 for m = 2,3,4.  (m = 1 is degenerate: Delta j = E_4^3 exactly,
so the combination vanishes identically.)

These f are holomorphic on H -- their only pole is at the cusp -- so plain arc quadrature is
correct for them and TRAP 4 does not arise.  And Phi(f) = c_f(0) = 0, so hat r_f = r_f with no
subtraction, which the run verifies rather than assumes.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, E4, Delta, jay, rvec, inW, decompose12,          # noqa: E402
                    arcint, rank, check_orientation)

mp.mp.dps = 30
K = 12

print("guard: Phi(E_12) = %s" % mp.nstr(check_orientation(), 12), flush=True)

E43 = lambda t: E4(t)**3
c0_E43 = arcint(E43, 11)
print("c_0(E_4^3) = Phi(E_4^3) = %s   (target 1)" % mp.nstr(c0_E43, 12), flush=True)

FORMS = [("Delta", Delta)]
for m in (2, 3, 4):
    raw = (lambda t, m=m: Delta(t) * jay(t)**m)
    c = arcint(raw, 11)
    FORMS.append(("Delta j^%d - %s E_4^3" % (m, mp.nstr(c, 8)),
                  (lambda t, raw=raw, c=c: raw(t) - c * E43(t))))

rows = []
print("=" * 82, flush=True)
for lbl, f in FORMS:
    Ph = arcint(f, 11)
    v = rvec(f, K, 11)
    a, al, be = decompose12(v)
    s_, u_ = inW(v, K)
    rows.append([a, al, be])
    print("%-34s Phi = %-14s |(1+S)| = %-10s |(1+U+U^2)| = %s"
          % (lbl, mp.nstr(Ph, 6), mp.nstr(s_, 4), mp.nstr(u_, 4)), flush=True)
    print("     W_- = %-24s W_+ = %-24s p_0 = %s"
          % (mp.nstr(a, 12), mp.nstr(al, 12), mp.nstr(be, 12)), flush=True)

M = mp.matrix(len(rows), 3)
for i, r in enumerate(rows):
    for j in range(3):
        M[i, j] = r[j]
print("=" * 82, flush=True)
sc = max(abs(M[i, j]) for i in range(M.rows) for j in range(3))
for tol in ('1e-18', '1e-12', '1e-8'):
    print("rank of the %d x 3 coordinate matrix (tol %s): %d"
          % (M.rows, tol, rank(M, mp.mpf(tol))), flush=True)
print("largest |p_0| coordinate among these f: %s   (scale %s)"
      % (mp.nstr(max(abs(r[2]) for r in rows), 8), mp.nstr(sc, 8)), flush=True)
print("PREDICTION: rank 2, and every p_0 coordinate zero.", flush=True)
