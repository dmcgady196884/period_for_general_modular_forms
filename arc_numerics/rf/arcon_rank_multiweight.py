"""Hole 4: does lem:arcon's "unless" ever bite, at weights other than 12?

lem:arcon ends: the defect eq:arcondefect is non-zero UNLESS D_E := r_{E_k}|(1+U+U^2) lies
in the span of {(X - p_j Y)^n}.  The Vandermonde argument in the proof is general and shows
the residue terms are independent of each other; what has only ever been checked at k = 12
is whether D_E can actually fall into that span.  If it can, the lemma is still true as
stated but its conclusion is empty for those configurations.

Setup: f with poles on gamma^arc has exactly TWO of them, at z and Sz with
arg z + arg Sz = pi, since j is 2-to-1 from the arc onto [0,1728].  So m = 2 and the span is
{(X-zY)^n, (X-SzY)^n}.  Measured here, over several weights and the whole cut:

    dist = || D_E - proj_span D_E || / || D_E ||,

so dist = 0 would mean the "unless" bites.  Basis convention matches common.rvec:
component l is the coefficient of X^{n-l} Y^l, hence (X-pY)^n has v[l] = C(n,l) (-p)^l.

E_k at general weight from its q-series, E_k = 1 - (2k/B_k) sum sigma_{k-1}(m) q^m; on the
arc |q| <= e^{-pi sqrt3} = 0.0043, so a couple of dozen terms is ample.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, jay, rvec, NU, dim_Sk                     # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "..", "end_to_end_numerical_testing"))
from validate import _sigma                                            # noqa: E402

mp.mp.dps = 25
NQ = 30
WEIGHTS = [12, 14, 16, 18, 20, 22, 26]
XS = ['1', '50', '300', '864', '1200', '1700', '1727']


def Ek(k):
    coef = -2 * k / mp.bernoulli(k)
    def f(t):
        q = mp.e**(2 * I * pi * t)
        return 1 + coef * sum(_sigma(k - 1, m) * q**m for m in range(1, NQ + 1))
    return f


def polyvec(p, n):
    return [mp.binomial(n, l) * (-p)**l for l in range(n + 1)]


def dist_to_span(d, vs):
    """|| d - proj_span d || / || d || , Hermitian inner product"""
    G = [[sum(mp.conj(a) * b for a, b in zip(u, v)) for v in vs] for u in vs]
    rhs = [sum(mp.conj(u) * a for u, a in zip(v, d)) for v in vs]
    M = mp.matrix(len(vs), len(vs))
    for r in range(len(vs)):
        for c in range(len(vs)):
            M[r, c] = G[r][c]
    b = mp.matrix(rhs)
    co = mp.lu_solve(M, b)
    res = [a - sum(co[j] * vs[j][idx] for j in range(len(vs)))
           for idx, a in enumerate(d)]
    nrm = lambda v: mp.sqrt(sum(abs(x)**2 for x in v))
    return nrm(res) / nrm(d)


print("dps = %d, NQ = %d\n" % (mp.mp.dps, NQ), flush=True)
print("  k   dimS_k   " + "".join("x=%-11s" % x for x in XS), flush=True)
for k in WEIGHTS:
    n = k - 2
    rE = rvec(Ek(k), k)
    DE = NU(rE, k)
    if max(abs(x) for x in DE) < mp.mpf('1e-18') * max(abs(x) for x in rE):
        print("  %-3d %-8d  D_E is ZERO -- span condition vacuous" % (k, dim_Sk(k)),
              flush=True)
        continue
    row = ""
    for xs in XS:
        x = mp.mpf(xs)
        z = mp.findroot(lambda t: jay(t) - x, mp.e**(I * mp.mpf('1.85')))
        z = z / abs(z)                       # pin it to the unit arc
        sz = -1 / z
        d = dist_to_span(DE, [polyvec(z, n), polyvec(sz, n)])
        row += "%-13s" % mp.nstr(d, 6)
    print("  %-3d %-8d  %s" % (k, dim_Sk(k), row), flush=True)
