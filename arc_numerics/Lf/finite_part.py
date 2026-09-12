"""The delta^0 term of L*(f_{delta,w},s) at a pinched (elliptic) pole, and whether logs occur.

WHY.  A pole approaching the contour from ONE side gives a finite limit: the antiderivative
(u - i delta)^{1-m}/(1-m) is single-valued along the path, so the integral is a difference of
endpoint values, bounded as delta -> 0 for every m >= 2 (m = 1 gives the log, hence the
Sokhotski-Plemelj pair PV +- i pi Res).  The divergence is caused by PINCHING, which happens
only at an elliptic point: retracting i downward puts a pole at i(1 - e) and, by modularity,
its partner at i/(1 - e), one on each side of the contour.

Scaling L* by delta^kappa is therefore not a definition -- it is not linear in f and would
annihilate the contribution of every pole away from the pinch.  The finite part does neither.
Wanted: the delta^0 coefficient of

    L*(f_{delta,w}, s) = c_{-1}/delta + c_0 + c_1 delta + ...     (+ c_L log delta ?)

SIMPLEST PINCH.  f = E_4^3/(j - 1728), k = 12: a pole of order P = m n - a = 2 at i.  The
divergence is then a simple pole in delta, so c_0 is a clean extraction and the log question
is sharp.  Retraction z_w(delta) = i(1 - delta(1-y)) runs straight down the imaginary axis;
x(delta) = j(z_w(delta)) is real and > 1728, hence off SL2(Z).gamma^arc as required.

The displaced pole sits (1-y) delta below i and its partner the same distance above, both at
the arc's midpoint, so the theta-mesh is clustered at pi/2; every value is computed at two
depths so an unresolved row announces itself.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, jay, ktil                             # noqa: E402

mp.mp.dps = 30
KK = 12
YW = mp.mpf('0.95')
EXPS = [-1, 0, 1, 2, 3]


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


def Lstar(x, s, depth):
    f = lambda t: E4(t)**3 / (jay(t) - x)
    return mp.e**(-I * pi * s / 2) * (arc_cl(lambda t: f(t) * t**(s - 1), depth)
                                      + arc_cl(lambda t: f(t) * ktil(t, s, KK), depth))


def solve(ds, vals, withlog):
    cols = list(EXPS) + (['log'] if withlog else [])
    A = mp.matrix(len(ds), len(cols))
    b = mp.matrix(len(ds), 1)
    for i, d in enumerate(ds):
        for jx, c in enumerate(cols):
            A[i, jx] = mp.log(d) if c == 'log' else d**c
        b[i] = vals[i]
    x = mp.lu_solve(A, b)
    return {c: x[jx] for jx, c in enumerate(cols)}


DS = [mp.mpf(2)**(-e) for e in range(6, 12)]
print("dps=%d  f = E_4^3/(j-1728), k=12, pole order 2 at i" % mp.mp.dps, flush=True)
print("retraction i -> i(1 - delta(1-y)), y = %s\n" % mp.nstr(YW, 4), flush=True)

for s in (mp.mpf(3), mp.mpf(7)):
    v1 = [Lstar(jay(I * (1 - d * (1 - YW))), s, 24) for d in DS]
    v2 = [Lstar(jay(I * (1 - d * (1 - YW))), s, 28) for d in DS]
    conv = max(abs(a - b) / abs(b) for a, b in zip(v1, v2))
    a = solve(DS, v2, False)
    bfit = solve(DS[:-1] + [DS[-1]], v2, True) if False else solve(DS, v2, False)
    aw = solve(DS[:len(EXPS) + 1], v2[:len(EXPS) + 1], True)
    print("s = %s   depth24-vs-28: %s" % (mp.nstr(s, 4), mp.nstr(conv, 4)), flush=True)
    print("   c_-1 = %s" % mp.nstr(a[-1], 16), flush=True)
    print("   c_0  = %-30s  (no log in the ansatz)" % mp.nstr(a[0], 16), flush=True)
    print("   c_0  = %-30s  c_log = %s   (log allowed)"
          % (mp.nstr(aw[0], 16), mp.nstr(aw['log'], 8)), flush=True)
    sub = solve(DS[1:], v2[1:], False)
    print("   c_0  = %-30s  (no log, grid shifted one octave)"
          % mp.nstr(sub[0], 16), flush=True)
