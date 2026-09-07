"""SUPERSEDED by i_limit.py.  Kept as a record and as corroboration.

There is NO obstruction at i at any pole order.  This script evaluates ON the contour,
which is the wrong thing to do at an elliptic point: the pole there is the merger of a
pinching stabiliser orbit, and any on-contour prescription is a degenerate limit.
Deform and rescale instead (i_limit.py): tilde r ~ eta^{-1/2} and the rescaled limit is
in W to 3e-25, ray-independent to 8e-25.

What survives from this script: the divergence it measures is real, and its extracted
A = lim eps * tilde r IS that same limit seen through a cruder regulator (A in W to
5.8e-5, below its own 1.7e-4 extraction error), while the Hadamard finite part misses W
by O(1).  So symmetric excision points the right way; it simply drops the indenting
semicircle's counterterm.

Original docstring follows.
---------------------------------------------------------------------------
Poles at i: the simple case converges, the 4 | k double case does NOT.

i is the fixed point of the symmetry the contour relies on ($S$ acts on the arc by
theta -> pi - theta), so oddness of log r forces every S-symmetric path through i and
no deformation avoids it.  The prescription is a principal value: excise
|theta - pi/2| < eps symmetrically and let eps -> 0.

Why that works for a SIMPLE pole and cannot for a double one.  Both segments are the
same arc here, so the effective kernel is the sum

    K(tau, s) = tau^{s-1} + ktil(tau, s, k).

Near i write u = theta - pi/2, so tau - i ~ -u.  A simple pole contributes ~ 1/u,
which is ODD, and symmetric excision cancels it.  A double pole contributes ~ 1/u^2,
which is EVEN, and symmetric excision does NOT: the integral over |u| > eps behaves
like 2 a_{-2} K(i,s) / eps.  That vanishes only if K(i,s) = 0, and it does not --
|K(i,s)| runs from 4.0 to 2344 at k = 20.

lem:ellLaurent permits the even order exactly when 4 | k: at mu = P = 2 the relation
reads r*(2) = i^k r*(2), so even-order poles are forbidden for k = 2 mod 4 and allowed
for 4 | k.  So the divergence is not an artifact of a bad test form; it is the generic
situation at 4 | k.

Test forms (both at dim S_k = 1, so dim W = 3 and the check is not degenerate):
    k = 18   E4^6 / E6       simple pole at i   (control: must converge)
    k = 20   E4^11 / E6^2    double pole at i   (must diverge like 1/eps)

Reference form: any holomorphic weight-k form with constant term 1 has Phi = 1, so
E4^3 E6 at k = 18 and E4^5 at k = 20.  These differ from E_k by a cusp form, which
shifts tilde r within W and so cannot affect a W-membership verdict.

Last block extracts the Hadamard finite part by fitting A/eps + B + C eps + ... and
asks whether B lies in W.  If it does, def:FP is the right prescription at 4 | k and
the PV simply is not.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common import mp, I, pi, E4, E6, ktil, check_orientation, inW, arcint, E12

mp.mp.dps = 30

print("Phi(E_12) =", mp.nstr(check_orientation(), 12), " (guard passed)\n", flush=True)


def _nodes_toward(a, b, depth):
    """nodes on [a,b] clustered geometrically toward b (the excision edge)"""
    out, d = [a], b - a
    for _ in range(depth):
        d /= 2
        out.append(b - d)
    return sorted(set(out + [b]))


def arc_excised(g, eps, depth=9):
    """int over gamma^arc with |theta - pi/2| < eps removed, oriented rho -> rho+1.

    Panels cluster into the excision edge: the pole sits eps away from it, so uniform
    panels of width ~(pi/6 - eps)/N are wider than the distance to the singularity for
    any practical N.  Getting this wrong makes the k=20 numbers pure noise.
    """
    tot = mp.mpc(0)
    for a, b, toward_b in ((pi / 3, pi / 2 - eps, True),
                           (pi / 2 + eps, 2 * pi / 3, False)):
        th = (_nodes_toward(a, b, depth) if toward_b
              else [a + b - x for x in reversed(_nodes_toward(a, b, depth))])
        th = sorted(set(th))
        tot += sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
                   for u, v in zip(th[:-1], th[1:]))
    return -tot


def rvec_excised(g, k, eps, depth=9):
    n = k - 2
    out = []
    for l in range(n + 1):
        s = mp.mpf(l + 1)
        L = (arc_excised(lambda t: g(t) * t**(s - 1), eps, depth)
             + arc_excised(lambda t: g(t) * ktil(t, s, k), eps, depth))
        out.append((2 * pi * I)**(n + 1) * (-1)**l * mp.binomial(n, l) * L)
    return out


CASES = [
    ("k=18  E4^6/E6     SIMPLE", 18, lambda t: E4(t)**6 / E6(t),
     lambda t: E4(t)**3 * E6(t)),
    ("k=20  E4^11/E6^2  DOUBLE", 20, lambda t: E4(t)**11 / E6(t)**2,
     lambda t: E4(t)**5),
]

EPS = [mp.mpf(x) for x in ('0.08', '0.04', '0.02', '0.01', '0.005')]

if len(sys.argv) > 1:                      # e.g. `python pole_at_i.py 20`
    want = {int(a) for a in sys.argv[1:]}
    CASES = [c for c in CASES if c[1] in want]

print("=" * 74)
print("RESOLUTION CHECK -- depth 8 vs depth 11 at the tightest eps, double pole.")
print("If these disagree, the sweep below is noise rather than physics.")
print("=" * 74)
for name, k, f, ref in CASES:
    if k != 20:
        continue
    eps = EPS[-1]
    a8, a11 = arc_excised(f, eps, 8), arc_excised(f, eps, 11)
    print("  eps=%s  Phi(d8)=%s  Phi(d11)=%s  rel diff=%s"
          % (mp.nstr(eps, 3), mp.nstr(a8, 12), mp.nstr(a11, 12),
             mp.nstr(abs(a8 - a11) / abs(a11), 4)), flush=True)
print()

store = {}
for name, k, f, ref in CASES:
    print("=" * 74)
    print(name, "   (dim W = 3)")
    print("=" * 74)
    rows = []
    for eps in EPS:
        rE = rvec_excised(ref, k, eps)
        rf = rvec_excised(f, k, eps)
        Phi = arc_excised(f, eps)
        tr = [a - Phi * b for a, b in zip(rf, rE)]
        sc = max(abs(x) for x in tr)
        s_, u_ = inW(tr, k)
        rows.append((eps, tr, sc))
        print("  eps=%-7s |Phi|=%-14s |tilde r|=%-16s in W: %s / %s"
              % (mp.nstr(eps, 3), mp.nstr(abs(Phi), 8), mp.nstr(sc, 9),
                 mp.nstr(s_, 4), mp.nstr(u_, 4)), flush=True)
    store[k] = rows
    print("  |tilde r| ratios between successive eps (2.0 each = pure 1/eps blowup):")
    print("   ", ", ".join(mp.nstr(rows[m + 1][2] / rows[m][2], 6)
                           for m in range(len(rows) - 1)))
    print()

print("=" * 74)
print("Hadamard finite part:  fit tilde r = A/eps + B + C eps + D eps^2 and test B")
print("=" * 74)
for name, k, f, ref in CASES:
    rows = store[k]
    n = k - 2
    M = mp.matrix(len(rows), 4)
    for m, (eps, tr, sc) in enumerate(rows):
        M[m, 0] = 1 / eps
        M[m, 1] = mp.mpf(1)
        M[m, 2] = eps
        M[m, 3] = eps**2
    A, B = [], []
    for l in range(n + 1):
        y = mp.matrix([rows[m][1][l] for m in range(len(rows))])
        cof = mp.lu_solve(M.T * M, M.T * y)
        A.append(cof[0])
        B.append(cof[1])
    print("  %s" % name)
    print("     |A| (1/eps coefficient) = %-16s  A in W:  %s / %s"
          % (mp.nstr(max(abs(x) for x in A), 8),
             *[mp.nstr(x, 5) for x in inW(A, k)]))
    print("     |B| (finite part)       = %-16s  B in W:  %s / %s"
          % (mp.nstr(max(abs(x) for x in B), 8),
             *[mp.nstr(x, 5) for x in inW(B, k)]), flush=True)

    # independent route to A: Neville-extrapolate eps * tilde r to eps = 0
    from common import neville
    xs = [rows[m][0] for m in range(len(rows))]
    Aex, sh = [], []
    for l in range(n + 1):
        val, s2 = neville(xs, [rows[m][0] * rows[m][1][l] for m in range(len(rows))])
        Aex.append(val)
        sh.append(s2)
    print("     eps*tilde r extrapolated to eps=0:  |A'| = %s  [shift %s]"
          % (mp.nstr(max(abs(x) for x in Aex), 8), mp.nstr(max(sh), 3)))
    print("     A' in W:  %s / %s" % tuple(mp.nstr(x, 5) for x in inW(Aex, k)))
    print("     max|A - A'| / |A| = %s"
          % mp.nstr(max(abs(a - b) for a, b in zip(A, Aex))
                    / max(abs(x) for x in A), 5), flush=True)
