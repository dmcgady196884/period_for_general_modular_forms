"""Parity split of W, computed without the per-row normalisation that corrupted theta_k18's PART 6.

theta_k18_DeltaE4.py reported dim W^even = 2 AND dim W^odd = 2 with dim W = 3, which is impossible,
and in the same run reported rank 4 for the 19 coefficients c_m.  Both are the same artifact:
each row was normalised by its OWN max, so a row that is essentially zero (c_8, and the odd part of
a basis vector that is nearly even) gets its roundoff amplified to O(1) and contributes a spurious
direction.  Normalise by a GLOBAL scale and drop rows below it instead.

W IS eps-stable, so the split is genuine.  With eps: P(X,Y) -> P(-X,Y),

    eps S eps^-1 = -S = S in PSL_2(Z)   (n even, so -I acts trivially on V_n),
    eps U eps^-1 = S U^-1 S^-1 ,

and on the subspace where P|(1+S) = 0 one has P|_S = -P, so P|(1+ eps U eps^-1 + (...)^2) = 0 is
equivalent to P|(1+U+U^2) = 0.  Hence eps maps W to W, and W = W^even + W^odd with eps acting as
(-1)^l on the X^{n-l}Y^l coefficient.

PREDICTION.  p_0 = X^n - Y^n has l = 0 and l = n, both even, so p_0 lies in W^even.  The cuspidal
part contributes r^+ (even) and r^- (odd), dim S_k of each.  So

    dim W^even = dim S_k + 1 ,   dim W^odd = dim S_k ,

i.e. (2,1) at k = 18 and at k = 12, and (3,2) at k = 24 where dim S_24 = 2.  No reference form and
no quadrature enter -- W is cut out of V_n by the S- and U-relations alone.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, NS, NU, opmat, nullspace, rank, dim_Sk, slash       # noqa: E402

mp.mp.dps = 40


def Wbasis(k):
    n = k - 2
    A, B = opmat(lambda v: NS(v, k), k), opmat(lambda v: NU(v, k), k)
    M = mp.matrix(2 * (n + 1), n + 1)
    for r in range(n + 1):
        for c in range(n + 1):
            M[r, c], M[r + n + 1, c] = A[r, c], B[r, c]
    return nullspace(M, mp.mpf('1e-25'))


def split_rank(vs, n, odd):
    """rank of the parity projections, normalised by a GLOBAL scale (never per-row)."""
    proj = [[x if (l % 2 == (1 if odd else 0)) else mp.mpc(0) for l, x in enumerate(v)]
            for v in vs]
    sc = max((max(abs(x) for x in p) for p in proj), default=mp.mpf(0))
    if sc == 0:
        return 0
    keep = [p for p in proj if max(abs(x) for x in p) > sc * mp.mpf('1e-20')]
    if not keep:
        return 0
    M = mp.matrix(len(keep), n + 1)
    for a, p in enumerate(keep):
        for l in range(n + 1):
            M[a, l] = p[l] / sc
    return rank(M, mp.mpf('1e-18'))


print("%-5s %-8s %-8s %-12s %-12s %s" % ("k", "dim S_k", "dim W", "W^even", "W^odd", "eps-stable?"),
      flush=True)
for k in (12, 16, 18, 20, 24, 26):
    n = k - 2
    WB = Wbasis(k)
    de = split_rank(WB, n, False)
    do = split_rank(WB, n, True)
    # eps-stability: project each basis vector and check the projection is still in W
    worst = mp.mpf(0)
    for v in WB:
        for od in (False, True):
            p = [x if (l % 2 == (1 if od else 0)) else mp.mpc(0) for l, x in enumerate(v)]
            s = max(abs(x) for x in p)
            if s < mp.mpf('1e-20'):
                continue
            d = max(max(abs(x) for x in NS(p, k)), max(abs(x) for x in NU(p, k))) / s
            worst = max(worst, d)
    print("%-5d %-8d %-8d %-12d %-12d %s"
          % (k, dim_Sk(k), len(WB), de, do, mp.nstr(worst, 4)), flush=True)
print("\npredicted: dim W^even = dim S_k + 1,  dim W^odd = dim S_k", flush=True)

# p_0 = X^n - Y^n really is in W and even
print("\np_0 = X^n - Y^n:", flush=True)
for k in (12, 18):
    n = k - 2
    p0 = [mp.mpc(0)] * (n + 1)
    p0[0], p0[n] = mp.mpc(1), mp.mpc(-1)
    s = max(max(abs(x) for x in NS(p0, k)), max(abs(x) for x in NU(p0, k)))
    print("   k=%d  in W? %s   supported at l = 0, %d (both even)"
          % (k, mp.nstr(s, 4), n), flush=True)
