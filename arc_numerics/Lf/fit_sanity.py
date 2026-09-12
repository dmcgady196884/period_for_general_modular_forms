"""Is "drift == -c_log log 2" evidence of a log, or an identity of the fitting scheme?

finite_part_scan.py extracts c_0 three ways on a geometric grid d_e = 2^-e:
    no-log  fit on d_1..d_N        (N = number of power columns)
    no-log  fit on d_2..d_{N+1}    ("grid shifted one octave")
    with-log fit on d_1..d_{N+1}
all by exact interpolation (square systems).  The first two differ by one Neville step
from the third, so a relation between their disagreement and c_log may hold whatever the
function is.  Feed the scheme a POLYNOMIAL-PLUS-POLE with no log at all and see.

Second question, the one that actually discriminates: how does the fitted c_log respond to
halving the whole grid?  A genuine log has a grid-independent coefficient.  A c_log that is
only absorbing the first omitted power delta^3 must shrink by 2^-3.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp                                                   # noqa: E402

mp.mp.dps = 40
EXPS = [-1, 0, 1, 2]


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


def report(tag, F, e0):
    need = len(EXPS) + 1
    g1 = [mp.mpf(2)**(-e) for e in range(e0, e0 + need)]
    g2 = [mp.mpf(2)**(-e) for e in range(e0 + 1, e0 + 1 + need)]
    n1, n2 = solve(g1, [F(d) for d in g1], EXPS, False), solve(g2, [F(d) for d in g2], EXPS, False)
    w1, w2 = solve(g1, [F(d) for d in g1], EXPS, True), solve(g2, [F(d) for d in g2], EXPS, True)
    drift = n2[0] - n1[0]
    pred = -w1['log'] * mp.log(2)
    print("%s   (grid starts at 2^-%d)" % (tag, e0), flush=True)
    print("   c_0 no-log  grid1 = %s   grid2 = %s"
          % (mp.nstr(n1[0], 16), mp.nstr(n2[0], 16)), flush=True)
    print("   drift = %-24s  -c_log log2 = %-24s  ratio = %s"
          % (mp.nstr(drift, 10), mp.nstr(pred, 10),
             mp.nstr(drift / pred, 10) if pred != 0 else 'n/a'), flush=True)
    print("   c_log grid1 = %-22s grid2 = %-22s ratio = %s"
          % (mp.nstr(w1['log'], 10), mp.nstr(w2['log'], 10),
             mp.nstr(w2['log'] / w1['log'], 8) if w1['log'] != 0 else 'n/a'), flush=True)


# (a) no log anywhere; the only thing outside the ansatz is delta^3 (and delta^4).
print("=" * 78, flush=True)
print("NO LOG.  F = 1/d + 2 + 3d + 4d^2 + 5d^3 + 6d^4", flush=True)
F = lambda d: 1 / d + 2 + 3 * d + 4 * d**2 + 5 * d**3 + 6 * d**4
for e0 in (5, 8, 11):
    report("   ", F, e0)

# (b) same, plus a log of known size.
print("=" * 78, flush=True)
print("WITH LOG.  F = 1/d + 2 + 3d + 4d^2 + 5d^3 + 6d^4 + (1e-5) log d", flush=True)
G = lambda d: F(d) + mp.mpf('1e-5') * mp.log(d)
for e0 in (5, 8, 11):
    report("   ", G, e0)
