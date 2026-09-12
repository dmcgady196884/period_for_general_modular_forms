"""Is the delta^0 term of L*(f_{delta,w},s) contaminated by log delta, at i and at rho?

fit_sanity.py shows the check used earlier -- "c_0 drifts by -c_log log 2 when the grid is
halved" -- is an exact identity of the interpolation scheme, true whether or not the function
has a log.  It is not evidence.  The discriminant is how the FITTED c_log responds to halving
the whole grid:

    c_log(grid/2) / c_log(grid)  ->  1      genuine log delta (coefficient is grid-independent)
                                 ->  2^-3   no log; c_log is only absorbing the first omitted
                                            power delta^3, which shrinks with the grid

(fit_sanity.py measures 0.1249 and 0.99998 on synthetic functions of each kind.)  The ratio is
only sharp once the grid is small enough that delta^3 has died: there the synthetic genuine-log
case reads 0.271 at 2^-5 but 0.992 at 2^-8, so the grid start is swept.

PREDICTION.  A log sits at order delta^0 exactly where the local model integral
int u^b du/(u^n - 1)^mu turns logarithmic, i.e. at the Taylor coefficient g_{P-1} of the regular
part, and it survives only if that coefficient survives symmetrisation over the order-n
stabiliser.  At i (n=2) it should not.  At rho (n=3) the relative root is omega^{b+1} = omega^{3m}
= 1, so the coefficient is proportional to K(rho) - K(tau_0), which vanishes only for
k = 2 mod 6; here k = 4.  So: ratio -> 2^-3 at i, ratio -> 1 at rho.

Both cases are m=1, P=2, the cheapest pinch and the one where the delta^-1 cancellation is mildest
(at delta = 2^-15 the leading term is only ~4 orders above c_0, so 30 dps is not strained).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, jay, ktil                             # noqa: E402

mp.mp.dps = 30
YW = mp.mpf('0.95')
RHO1 = mp.e**(I * pi / 3)
EXPS = [-1, 0, 1, 2]

CASES = [
    ("i    P=2 k=12 (E_4^3, m=1)  [control: expect 2^-3]", I, lambda t: E4(t)**3, 1, 12),
    ("rho  P=2 k=4  (E_4,   m=1)  [expect 1]", RHO1, lambda t: E4(t), 1, 4),
]


def arc_cl(g, depth):
    a, b = pi / 3, 2 * pi / 3
    ns = [a, b]
    for base in (a, b, pi / 2):
        d = (b - a) / 4
        for _ in range(depth):
            d /= 2
            ns += [base - d, base + d]
        ns.append(base)
    ns = sorted(set(x for x in ns if a <= x <= b))
    return -sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
                for u, v in zip(ns[:-1], ns[1:]))


def Lstar(gf, m, x, s, k, depth):
    f = lambda t: gf(t) / (jay(t) - x)**m
    return mp.e**(-I * pi * s / 2) * (arc_cl(lambda t: f(t) * t**(s - 1), depth)
                                      + arc_cl(lambda t: f(t) * ktil(t, s, k), depth))


def solve(ds, vals, withlog):
    cols = list(EXPS) + (['log'] if withlog else [])
    A = mp.matrix(len(cols), len(cols))
    b = mp.matrix(len(cols), 1)
    for i in range(len(cols)):
        for jx, c in enumerate(cols):
            A[i, jx] = mp.log(ds[i]) if c == 'log' else ds[i]**c
        b[i] = vals[i]
    x = mp.lu_solve(A, b)
    return {c: x[jx] for jx, c in enumerate(cols)}


for lbl, z, gf, m, k in CASES:
    print("=" * 78, flush=True)
    print(lbl, flush=True)
    for s in (mp.mpf(7), mp.mpf(11)):
        print("  s = %s" % mp.nstr(s, 3), flush=True)
        for e0 in (5, 8, 11):
            need = len(EXPS) + 1
            g1 = [mp.mpf(2)**(-e) for e in range(e0, e0 + need)]
            g2 = [mp.mpf(2)**(-e) for e in range(e0 + 1, e0 + 1 + need)]
            res = []
            for ds in (g1, g2):
                v1 = [Lstar(gf, m, jay(z - d * (z - I * YW)), s, k, 26) for d in ds]
                v2 = [Lstar(gf, m, jay(z - d * (z - I * YW)), s, k, 30) for d in ds]
                conv = max(abs(a - b) / abs(b) for a, b in zip(v1, v2))
                res.append((conv, solve(ds, v2, False), solve(ds, v2, True)))
            (c1, n1, w1), (c2, n2, w2) = res
            r = w2['log'] / w1['log'] if w1['log'] != 0 else mp.mpf(0)
            print("    2^-%-3d c_log %-22s -> %-22s  RATIO %s"
                  % (e0, mp.nstr(w1['log'], 8), mp.nstr(w2['log'], 8), mp.nstr(r, 8)),
                  flush=True)
            print("           c_0 no-log %-24s c_0 w/log %-24s conv %s/%s"
                  % (mp.nstr(n1[0], 14), mp.nstr(w1[0], 14),
                     mp.nstr(c1, 3), mp.nstr(c2, 3)), flush=True)
