"""Finite part at BOTH elliptic points and several pole orders: are there log delta terms?

finite_part.py settled one case (P = 2 at i): no log, c_0 stable to 13 digits under both a
grid shift and the addition of a log term.  That is one point of one order at one of the two
elliptic points, and rho is the case where a log HAS already appeared once (the c_0 log J
divergence of the block mode sum), so it must be checked before "no logs" goes into a
definition.

Cases.  At i, j - 1728 has a double zero, so P = 2m - a with a = ord_i g:
    P=2  g=E_4^3, m=1, k=12      P=3  g=E_6,  m=2, k=6      P=4  g=E_4^3, m=2, k=12
At rho, j has a triple zero, so P = 3m - a with a = ord_rho g:
    P=2  g=E_4,   m=1, k=4       P=3  g=Delta, m=1, k=12    P=4  g=E_4^2, m=2, k=8

Retraction: z_w(delta) = z - delta(z - w) with w = i y.  At i that runs straight down; at
rho we retract the endpoint rho+1 = e^{i pi/3}, which lies in the range [pi/3, pi/2] of
def:retract, and modularity carries rho along with it.  In both cases the deformed form is
g/(j - x(delta))^m with x(delta) = j(z_w(delta)), so nothing but w is chosen.

Ansatz: L*(delta) = sum_{r=-(P-1)}^{2} c_r delta^r  (+ c_L log delta), solved on a square
grid, then repeated on a grid shifted by one octave.  Every value is computed at two mesh
depths.  s = 7 and 11 throughout: s = 3 is unusable at k = 8 (structural zero) and at k = 6
(s = k/2 with i^k = -1).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, E6, Delta, jay, ktil                  # noqa: E402

mp.mp.dps = 30
YW = mp.mpf('0.95')
RHO1 = mp.e**(I * pi / 3)
SVALS = [mp.mpf(7), mp.mpf(11)]

CASES = [
    ("i    P=2 k=12 (E_4^3, m=1)", I, lambda t: E4(t)**3, 1, 2, 12),
    ("i    P=3 k=6  (E_6,   m=2)", I, lambda t: E6(t), 2, 3, 6),
    ("i    P=4 k=12 (E_4^3, m=2)", I, lambda t: E4(t)**3, 2, 4, 12),
    ("rho  P=2 k=4  (E_4,   m=1)", RHO1, lambda t: E4(t), 1, 2, 4),
    ("rho  P=3 k=12 (Delta, m=1)", RHO1, lambda t: Delta(t), 1, 3, 12),
    ("rho  P=4 k=8  (E_4^2, m=2)", RHO1, lambda t: E4(t)**2, 2, 4, 8),
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


def solve(ds, vals, exps, withlog):
    cols = list(exps) + (['log'] if withlog else [])
    A = mp.matrix(len(cols), len(cols))
    b = mp.matrix(len(cols), 1)
    for i in range(len(cols)):
        for jx, c in enumerate(cols):
            A[i, jx] = mp.log(ds[i]) if c == 'log' else ds[i]**c
        b[i] = vals[i]
    x = mp.lu_solve(A, b)
    return {c: x[jx] for jx, c in enumerate(cols)}


for lbl, z, gf, m, P, k in CASES:
    exps = list(range(-(P - 1), 3))
    need = len(exps) + 1
    g1 = [mp.mpf(2)**(-e) for e in range(5, 5 + need)]
    g2 = [mp.mpf(2)**(-e) for e in range(6, 6 + need)]
    print("=" * 78, flush=True)
    print("%s   exponents %d..2" % (lbl, -(P - 1)), flush=True)
    for s in SVALS:
        out = []
        for ds in (g1, g2):
            v1 = [Lstar(gf, m, jay(z - d * (z - I * YW)), s, k, 26) for d in ds]
            v2 = [Lstar(gf, m, jay(z - d * (z - I * YW)), s, k, 30) for d in ds]
            conv = max(abs(a - b) / abs(b) for a, b in zip(v1, v2))
            nl = solve(ds, v2, exps, False)
            wl = solve(ds, v2, exps, True)
            out.append((conv, nl, wl))
        (c1, n1, w1), (c2, n2, w2) = out
        print("  s=%-4s depth-conv %s / %s" % (mp.nstr(s, 3), mp.nstr(c1, 3),
                                               mp.nstr(c2, 3)), flush=True)
        print("     c_lead(%d) = %s" % (-(P - 1), mp.nstr(n1[-(P - 1)], 12)), flush=True)
        print("     c_0 grid1 = %-26s  grid2 = %s"
              % (mp.nstr(n1[0], 16), mp.nstr(n2[0], 16)), flush=True)
        print("     c_0 w/log = %-26s  c_log = %s"
              % (mp.nstr(w1[0], 16), mp.nstr(w1['log'], 8)), flush=True)
