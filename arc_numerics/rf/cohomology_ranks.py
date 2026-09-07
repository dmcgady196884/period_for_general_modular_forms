"""lem:kersum and lem:Qinv -- why a pole at rho is never an obstruction.

lem:georho leaves tilde r_f outside W with defect -Q_f.  A correction w exists iff
Q_f lies in im[(1+U+U^2)|ker(1+S)].  Two legs:

 LEG 1 (lem:kersum).  im[(1+U+U^2)|ker(1+S)] = V_n^U, the FULL U-fixed subspace.
     Since 1+S and 1+U+U^2 are self-adjoint for eq:inner and im(1+g) = V_n^g for
     finite-order g, this reduces to ker(1+S) + ker(1+U+U^2) = V_n, dually to
     V_n^S cap V_n^U = 0.  That intersection is V_n^{SL2(Z)} because <S,U> = SL2(Z):
     T-invariance forces c Y^n, which S sends to c X^n, so c = 0.
     Checked here: rank == dim V_n^U at every weight, and the fixed spaces meet in 0.
     NB dim ker(1+S) is n/2 for n = 0 mod 4 and n/2 + 1 for n = 2 mod 4, so the rank
     identity is NOT a restatement of one dimension formula.

 LEG 2 (lem:Qinv).  Q_f = (2 pi i)^{n+2} Res_{tau=p}[f (X - tau Y)^n], and the form
     omega_f = f (X - tau Y)^n dtau satisfies gamma^* omega_f = omega_f | gamma^{-1}.
     With U p = p for p = rho+1 this gives Q_f | U = Q_f, i.e. Q_f in V_n^U.  That is
     lem:ellLaurent transcribed to V_n.
     Checked here by computing the residue as a circle integral and slashing by U,
     against a control of generic (non-modular) Laurent weights through the same
     B_alpha assembly.

Together: Q_f in V_n^U = im, so w exists.  Expected output: LEG 1 matches at every
weight; LEG 2 gives ~1e-28 for genuine modular forms and O(1) for the control.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common import (mp, I, pi, E4, E6, Delta, slash, Sm, Um, TAU0,
                    opmat, nullspace, rank, dim_Sk, NS, NU)

mp.mp.dps = 30


def stack(A, B):
    M = mp.matrix(A.rows + B.rows, A.cols)
    for i in range(A.rows):
        for j in range(A.cols):
            M[i, j] = A[i, j]
    for i in range(B.rows):
        for j in range(B.cols):
            M[A.rows + i, j] = B[i, j]
    return M


def ident(n):
    M = mp.matrix(n + 1, n + 1)
    for i in range(n + 1):
        M[i, i] = mp.mpc(1)
    return M


print("=" * 78)
print("LEG 1 (lem:kersum):  rank == dim V_n^U ?   and   V_n^S cap V_n^U = 0 ?")
print("=" * 78)
print("   k    n | dimker(1+S)  dim W  2dimS_k+1 | rank  dim V_n^U | V^S cap V^U")
print("-" * 78)
ok = True
for k in range(4, 29, 2):
    n = k - 2
    Smat = opmat(lambda v, k=k: slash(v, k - 2, Sm), k)
    Umat = opmat(lambda v, k=k: slash(v, k - 2, Um), k)
    A = opmat(lambda v, k=k: NS(v, k), k)          # 1 + S
    NUm = opmat(lambda v, k=k: NU(v, k), k)        # 1 + U + U^2
    Id = ident(n)
    dim_kerS = (n + 1) - rank(A)
    dim_W = (n + 1) - rank(stack(A, NUm))
    dim_VU = rank(NUm)
    basis = nullspace(A)
    P = mp.matrix(n + 1, max(len(basis), 1))
    for c, v in enumerate(basis):
        col = NU(v, k)
        for r in range(n + 1):
            P[r, c] = col[r]
    rk = rank(P) if basis else 0
    dim_fix = (n + 1) - rank(stack(Id - Smat, Id - Umat))
    good = (rk == dim_VU) and (dim_fix == 0) and (dim_W == 2 * dim_Sk(k) + 1)
    ok &= good
    print("  %3d  %3d |     %3d       %3d      %3d     |  %3d      %3d      |     %d %s"
          % (k, n, dim_kerS, dim_W, 2 * dim_Sk(k) + 1, rk, dim_VU, dim_fix,
             "" if good else "  <-- MISMATCH"))
print("-" * 78)
print("LEG 1:", "PASS" if ok else "FAIL")

print()
print("=" * 78)
print("LEG 2 (lem:Qinv):  is Res_p[f (X - tau Y)^n] fixed by U ?")
print("=" * 78)
p = TAU0
print("p = %s   |U p - p| = %s" % (mp.nstr(p, 8), mp.nstr(abs((p - 1) / p - p), 5)))

CASES = [("k=12  Delta^2/E4^3  P=3", 12, lambda t: Delta(t)**2 / E4(t)**3),
         ("k=16  Delta^2/E4^2  P=2", 16, lambda t: Delta(t)**2 / E4(t)**2),
         ("k=20  E6^4/E4       P=1", 20, lambda t: E6(t)**4 / E4(t)),
         ("k=24  Delta^3/E4^3  P=3", 24, lambda t: Delta(t)**3 / E4(t)**3)]

for name, k, f in CASES:
    n, rad = k - 2, mp.mpf('0.12')

    def coef(l, f=f, rad=rad):
        g = lambda th: (f(p + rad * mp.e**(I * th)) * (p + rad * mp.e**(I * th))**l
                        * rad * I * mp.e**(I * th))
        return mp.quad(g, [0, pi / 2, pi, 3 * pi / 2, 2 * pi]) / (2 * pi * I)

    R = [(-1)**l * mp.binomial(n, l) * coef(l) for l in range(n + 1)]
    sc = max(abs(x) for x in R)
    dev = max(abs(a - b) for a, b in zip(R, slash(R, n, Um))) / sc
    print("  %-24s |R|U - R| / |R| = %s" % (name, mp.nstr(dev, 6)), flush=True)

print()
print("  control: same B_alpha assembly, GENERIC non-modular weights (k=12, P=3)")
n = 10
for trial in range(3):
    ws = [mp.mpc(mp.cos(3 * trial + a), mp.sin(2 * trial - a)) for a in range(1, 4)]
    Rg = [mp.mpc(0)] * (n + 1)
    for al in range(1, 4):
        m = n - al + 1
        pref = mp.binomial(n, al - 1) * (-1)**(al - 1)
        for t in range(m + 1):
            Rg[t + al - 1] += ws[al - 1] * pref * mp.binomial(m, t) * (-p)**t
    scg = max(abs(x) for x in Rg)
    print("     trial %d:  |R|U - R| / |R| = %s"
          % (trial + 1, mp.nstr(max(abs(a - b) for a, b in
                                    zip(Rg, slash(Rg, n, Um))) / scg, 6)))
