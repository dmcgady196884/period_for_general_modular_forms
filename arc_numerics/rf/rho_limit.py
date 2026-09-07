"""Poles at rho: fixing the correction by CONTINUITY rather than by choosing w.

lem:georho leaves tilde r_f outside W with defect -Q_f, and lem:kersum only says a
correction w EXISTS (unique mod W, hence not pinned).  Continuity pins it.

Setup.  f_{j0} = Delta/(j - j0), k = 12.  j has a TRIPLE zero at rho, so for j0 != 0
the form has three SIMPLE poles z, Uz, U^2z at distance |j0|^{1/3} from rho, merging
into the triple pole of Delta/j = Delta^2/E4^3 at j0 = 0.  Modularity forces those
three to be a full U-orbit -- an asymmetric split is impossible -- so the only freedom
is the path in j0.  For j0 != 0 the poles are off the arc, so thm:geoperiod applies
with NO correction, and

    tilde r_{f_{j0}} = j0^{-2/3} ( c v + O(j0^{1/3}) ),   v in W.

The cube root is an overall scalar, so the LINE C.v is free of ray/branch ambiguity.

What this script establishes, in order:
 (A) The quadrature is not the limitation: depth 10 vs 13 agree at ~1e-41, so all
     residuals below are genuine j0-truncation.
 (B) Neville extrapolation in u = j0^{1/3} of the direction normalised by its l=5
     entry (which quotients out the divergent modulus AND the cube-root branch).
     NB the odd entries land on 2/21, -25/42, 1 -- but that is FORCED, not a discovery:
     at k=12 the odd part of W is 1-dimensional, so every element of W has odd part
     proportional to W_-.  The real content is the two ratios alpha/a, beta/a.
 (C) Canonicality under change of family.  B and C differ from A by HOLOMORPHIC forms;
     since tilde r is linear, that adds a bounded vector against a leading term
     diverging like j0^{-2/3}, so the direction must converge to the same limit.
     D is the control: E4^3 has a triple zero at rho that cancels the pole, so
     E4^3/(j - j0) must stay BOUNDED.  If D diverges the whole picture is wrong.

Caveat this script does NOT settle: the limit is not intrinsic.  It moves with the
reference Eisenstein series -- see rf/eisenstein_dependence.py.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common import (mp, I, pi, E4, Delta, jay, check_orientation, arcint,
                    rvec, E12, inW, neville, decompose12)

mp.mp.dps = 40
KK, NN = 12, 10
PSI = pi                      # ray, kept off the cut j(arc) = [0, 1728]

print("Phi(E_12) =", mp.nstr(check_orientation(), 12), " (guard passed)", flush=True)

_cache = {}


def tr_of(f, depth):
    if ('E12', depth) not in _cache:
        _cache[('E12', depth)] = rvec(E12, KK, depth)
    rE = _cache[('E12', depth)]
    Phi = arcint(f, depth)
    return [a - Phi * b for a, b in zip(rvec(f, KK, depth), rE)], Phi


famA = lambda j0: (lambda t: Delta(t) / (jay(t) - j0))
famB = lambda j0: (lambda t: (Delta(t) + 3 * E4(t)**3) / (jay(t) - j0))
famC = lambda j0: (lambda t: Delta(t) / (jay(t) - j0) + 7 * Delta(t))
famD = lambda j0: (lambda t: E4(t)**3 / (jay(t) - j0))

print()
print("=" * 76)
print("(A) quadrature floor: same j0, depth 10 vs depth 13")
print("=" * 76)
for e in ('1e-3', '1e-5'):
    j0 = mp.mpf(e) * mp.e**(I * PSI)
    t10, _ = tr_of(famA(j0), 10)
    t13, _ = tr_of(famA(j0), 13)
    n10 = [x / t10[5] for x in t10]
    n13 = [x / t13[5] for x in t13]
    print("  |j0|=%-5s max|normalised(d10) - normalised(d13)| = %s"
          % (e, mp.nstr(max(abs(a - b) for a, b in zip(n10, n13)), 6)), flush=True)

print()
print("=" * 76)
print("(B) extrapolation in u = j0^{1/3}")
print("=" * 76)
DEPTH = 13
us, vs = [], []
for e in ('1e-2', '1e-3', '1e-4', '1e-5', '1e-6'):
    j0 = mp.mpf(e) * mp.e**(I * PSI)
    tr, _ = tr_of(famA(j0), DEPTH)
    s_, u_ = inW(tr, KK)
    us.append(j0**(mp.mpf(1) / 3))
    vs.append([x / tr[5] for x in tr])
    print("  |j0|=%-5s |tilde r| = %-16s in W: %s / %s"
          % (e, mp.nstr(max(abs(x) for x in tr), 8), mp.nstr(s_, 4), mp.nstr(u_, 4)),
          flush=True)

vlim, shifts = [], []
for l in range(NN + 1):
    val, sh = neville(us, [vs[m][l] for m in range(len(us))])
    vlim.append(val)
    shifts.append(sh)

print()
print("  extrapolated v (v_5 = 1), with per-coefficient error bar:")
for l in range(NN + 1):
    print("    l=%2d  %-46s [shift %s]"
          % (l, mp.nstr(vlim[l], 20), mp.nstr(shifts[l], 3)))

print()
print("  odd entries vs (2/21, -25/42, 1) -- forced by dim W^odd = 1, so this is a")
print("  check on the extrapolation landing in W, not a result:")
targ = [mp.mpf(2) / 21, mp.mpf(-25) / 42, mp.mpf(1), mp.mpf(-25) / 42, mp.mpf(2) / 21]
for idx, l in enumerate((1, 3, 5, 7, 9)):
    print("    l=%2d  diff = %s" % (l, mp.nstr(abs(vlim[l] - targ[idx]), 4)))

a_, al_, be_ = decompose12(vlim)
print()
print("  THE CONTENT -- v = a W_- + alpha W_+ + beta p_0 (read-off, disjoint supports):")
print("    a       =", mp.nstr(a_, 20), "  (= 1/42 exactly)")
print("    alpha/a =", mp.nstr(al_ / a_, 20))
print("    beta/a  =", mp.nstr(be_ / a_, 20))
w0, w2, w4 = vlim[0] / I, vlim[2] / I, vlim[4] / I
print("    w4/w2   =", mp.nstr(w4 / w2, 18), "  (W_+ has (1,-3,3,-1), so -3 expected)")
print("    w2/w0   =", mp.nstr(w2 / w0, 18))
for nm, x in (("w2/w0", w2 / w0), ("w4/w2", w4 / w2), ("w0", w0)):
    try:
        rel = mp.pslq([mp.re(x), mp.mpf(1)], maxcoeff=10**6, maxsteps=10**5)
    except Exception:
        rel = None
    print("    PSLQ %-6s rational? %s" % (nm, rel))

print()
print("=" * 76)
print("(C) canonicality under change of family, at |j0| = 1e-4")
print("=" * 76)
j0 = mp.mpf('1e-4') * mp.e**(I * PSI)
base = None
for nm, fam in (("A  Delta/(j-j0)           ", famA),
                ("B  (Delta+3E4^3)/(j-j0)   ", famB),
                ("C  Delta/(j-j0) + 7Delta  ", famC),
                ("D  E4^3/(j-j0)  [CONTROL] ", famD)):
    tr, _ = tr_of(fam(j0), DEPTH)
    nv = [x / tr[5] for x in tr]
    if base is None:
        base, extra = nv, "<- reference"
    else:
        extra = "max|dir - dirA| = %s" % mp.nstr(
            max(abs(a - b) for a, b in zip(nv, base)), 6)
    print("  %s |tilde r| = %-18s %s"
          % (nm, mp.nstr(max(abs(x) for x in tr), 9), extra), flush=True)
print()
print("  D must be BOUNDED (its numerator's triple zero at rho kills the pole);")
print("  A, B, C all diverge like j0^{-2/3}.  D's direction comparison is meaningless.")
