"""Is Theta(z) -- the residue polynomial per unit residue -- a Laurent polynomial in z?

Xi_{f;z} = a_z Theta(z) at a simple pole, with (eq:respolyexp)

  Theta(z)_l = (2 pi i)^{n+2}(-1)^l C(n,l)[ z^l - z^{n-l} + ktil(z,l+1) - z^n ktil(-1/z,l+1) ]
               - 2 pi i (1 - z^n) (r_Ek)_l .

ktil(w,l+1) = -B_{l+1}(w+1)/(l+1) + (-1)^l B_{n-l+1}(w+1)/(n-l+1) is a POLYNOMIAL of degree
max(l+1, n-l+1) <= n+1.  So ktil(z,.) contributes up to z^{n+1}, and z^n ktil(-1/z,.) contributes
down to z^{n - max(l+1,n-l+1)}, whose minimum over l is z^{-1}, attained only at l = 0 and l = n.
Hence the PREDICTION

    Theta(z) = sum_{m=-1}^{n+1} c_m z^m ,      c_m in V_n ,

i.e. z Theta(z) is a polynomial of degree <= n+2, with a simple pole at z=0 and nothing else.

THE STRUCTURAL POINT.  Theta(z) lies in W for every z (prop:respolyW), the monomials z^m are
linearly independent as functions, and W is a linear subspace.  Therefore EVERY c_m lies in W.  At
k=12 that is 13 vectors inside a 3-dimensional space: ten linear relations, and the W-coordinates
of Theta(z) are three explicit Laurent polynomials in z.  That is the answer to "does the
localisation within W have any pretty dependence on z" -- if it survives.

TRANSCENDENCE SPLIT.  r_Ek enters only through -2 pi i (1-z^n) r_Ek, i.e. only at m = 0 and m = n.
So L(z) := [Theta(z) + 2 pi i (1-z^n) r_Ek] / (2 pi i)^{n+2} should have RATIONAL coefficients.
Tested by PSLQ on each coefficient of each monomial.

METHOD.  Interpolate z Theta(z) from n+3 sample points, then VERIFY at fresh points not used in
the fit -- a fit that is not tested out of sample proves nothing.  Everything is the closed form,
so this is exact algebra at dps 50; no quadrature anywhere.

ALSO ANSWERED HERE (question 1).  Is Theta(z) ever zero, and how much does its direction in P(W)
move with z?  For hat r_{f_z} the direction was nearly constant (3e-6 spread) because the v_m were
almost parallel -- that killed the CM test.  Theta is a different object and need not behave that
way; the coefficient rank below decides it.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, E4, E6, ktil, slash, rvec, NS, NU, inW,     # noqa: E402
                    opmat, nullspace, arcint, rank)

mp.mp.dps = 50
K, N = 12, 10
Sm = (0, -1, 1, 0)
E43 = lambda t: E4(t)**3
rE = rvec(E43, K, 11)
print("reference E_4^3  Phi = %s" % mp.nstr(arcint(E43, 11), 12), flush=True)


def Theta(z):
    P = [mp.binomial(N, l) * (-z)**l for l in range(N + 1)]
    PS = slash(P, N, Sm)
    out = []
    for l in range(N + 1):
        kz = ktil(z, mp.mpf(l + 1), K)
        kw = ktil(-1 / z, mp.mpf(l + 1), K)
        br = (P[l] - PS[l]) + mp.binomial(N, l) * (-1)**l * (kz - z**N * kw)
        out.append((2 * pi * I)**(N + 1) * br)
    return [2 * pi * I * (a - (1 - z**N) * b) for a, b in zip(out, rE)]


# ---------------------------------------------------------------- PART 1: Laurent fit
MS = list(range(-1, N + 2))                       # m = -1 .. n+1, so n+3 = 13 unknowns
FIT = [mp.mpf('0.31') + mp.mpf('1.27') * I, mp.mpf('-0.4') + mp.mpf('1.1') * I,
       mp.mpf('0.7') + mp.mpf('1.9') * I, I * mp.mpf('1.3'), mp.mpf('0.5') + I * mp.mpf('1.4'),
       mp.mpf('-0.25') + I * mp.mpf('1.6'), mp.mpf('0.15') + I * mp.mpf('2.1'),
       mp.mpf('0.9') + I * mp.mpf('1.05'), mp.mpf('-0.6') + I * mp.mpf('1.35'),
       mp.mpf('0.05') + I * mp.mpf('1.15'), mp.mpf('0.44') + I * mp.mpf('1.72'),
       mp.mpf('-0.8') + I * mp.mpf('2.4'), mp.mpf('0.33') + I * mp.mpf('1.05')]
assert len(FIT) == len(MS)

V = mp.matrix(len(MS), len(MS))
for a, z in enumerate(FIT):
    for b, m in enumerate(MS):
        V[a, b] = z**m
TH = [Theta(z) for z in FIT]
C = []                                            # C[b] = c_{MS[b]} as a V_n vector
for l in range(N + 1):
    rhs = mp.matrix([TH[a][l] for a in range(len(MS))])
    C.append(mp.lu_solve(V, rhs))
c = [[C[l][b] for l in range(N + 1)] for b in range(len(MS))]

print("\nPART 1  out-of-sample verification of the Laurent form", flush=True)
for z in (mp.mpf('0.62') + I * mp.mpf('1.11'), mp.mpf('-0.37') + I * mp.mpf('1.83'),
          I * mp.mpf('2.6'), mp.mpf('0.5') + I * mp.sqrt(3) / 2):
    got = Theta(z)
    pred = [sum(c[b][l] * z**MS[b] for b in range(len(MS))) for l in range(N + 1)]
    e = max(abs(x - y) for x, y in zip(got, pred)) / max(abs(x) for x in got)
    print("   z = %-28s relative error %s" % (mp.nstr(z, 8), mp.nstr(e, 6)), flush=True)

# ---------------------------------------------------------------- PART 2: c_m in W?
print("\nPART 2  is every coefficient c_m itself in W?", flush=True)
for b, m in enumerate(MS):
    nb = max(abs(x) for x in c[b])
    if nb < mp.mpf('1e-30'):
        print("   m = %-3d  c_m = 0" % m, flush=True)
        continue
    s_, u_ = inW(c[b], K)
    print("   m = %-3d  |c_m| = %-14s |(1+S)| = %-11s |(1+U+U^2)| = %s"
          % (m, mp.nstr(nb, 8), mp.nstr(s_, 4), mp.nstr(u_, 4)), flush=True)

# ---------------------------------------------------------------- PART 3: rank
A, B_ = opmat(lambda v: NS(v, K), K), opmat(lambda v: NU(v, K), K)
MM = mp.matrix(2 * (N + 1), N + 1)
for r_ in range(N + 1):
    for cc in range(N + 1):
        MM[r_, cc], MM[r_ + N + 1, cc] = A[r_, cc], B_[r_, cc]
WB = nullspace(MM, mp.mpf('1e-30'))
print("\nPART 3  dim W = %d;  rank of the %d coefficient vectors"
      % (len(WB), len(MS)), flush=True)
G = mp.matrix(len(MS), N + 1)
for b in range(len(MS)):
    nb = max(abs(x) for x in c[b]) or mp.mpf(1)
    for l in range(N + 1):
        G[b, l] = c[b][l] / nb
for tol in ('1e-30', '1e-20', '1e-10'):
    print("   rank (tol %s) = %d" % (tol, rank(G, mp.mpf(tol))), flush=True)

# ---------------------------------------------------------------- PART 4: rationality
print("\nPART 4  L(z) := [Theta + 2 pi i (1-z^n) r_Ek]/(2 pi i)^{n+2}: rational coefficients?",
      flush=True)
for b, m in enumerate(MS):
    lv = []
    for l in range(N + 1):
        corr = 2 * pi * I * ((1 if m == 0 else 0) - (1 if m == N else 0)) * rE[l]
        lv.append((c[b][l] + corr) / (2 * pi * I)**(N + 2))
    nb = max(abs(x) for x in lv)
    if nb < mp.mpf('1e-30'):
        print("   m = %-3d  L_m = 0" % m, flush=True)
        continue
    hits = []
    for l in (0, 1, N // 2, N):
        x = lv[l]
        if abs(x) < mp.mpf('1e-30'):
            hits.append("%d:0" % l)
            continue
        # L_m is expected REAL rational; report the imaginary part as the guard
        pr, pim = mp.re(x), mp.im(x)
        q = mp.pslq([pr, mp.mpf(1)], tol=mp.mpf('1e-35'), maxcoeff=10**12, maxsteps=10**5)
        tag = ("%s/%s" % (-q[1], q[0])) if q else "irrational?"
        if abs(pim) > mp.mpf('1e-30') * max(abs(pr), mp.mpf(1)):
            tag += "+%si" % mp.nstr(pim, 4)
        hits.append("%d:%s" % (l, tag))
    print("   m = %-3d  |L_m| = %-14s  %s" % (m, mp.nstr(nb, 8), "  ".join(hits)), flush=True)

# ---------------------------------------------------------------- PART 5: direction
print("\nPART 4b  functional equation  Theta(-1/z) = -z^{-n} Theta(z)", flush=True)
# gamma_o(Sz) = c_{Sz} - c_z = -gamma_o(z), so Xi_{f;Sz} = -Xi_{f;z}; with Xi = a_z Theta(z)
# and a_{Sz} = z^n a_z this forces the displayed identity, hence c_{n-m} = (-1)^{m+1} c_m.
for z in (mp.mpf('0.31') + I * mp.mpf('1.27'), I * mp.sqrt(2),
          mp.mpf('-0.37') + I * mp.mpf('1.83')):
    lhs, rhs = Theta(-1 / z), [-x / z**N for x in Theta(z)]
    print("   z = %-24s relative error %s"
          % (mp.nstr(z, 8), mp.nstr(max(abs(a - b) for a, b in zip(lhs, rhs))
                                    / max(abs(x) for x in rhs), 6)), flush=True)
for b, m in enumerate(MS):
    mm = N - m
    if mm not in MS or mm < m:
        continue
    bb = MS.index(mm)
    sgn = (-1)**(m + 1)
    d = max(abs(c[bb][l] - sgn * c[b][l]) for l in range(N + 1))
    sc_ = max(abs(x) for x in c[bb]) or mp.mpf(1)
    print("   c_%d = %sc_%d :  %s" % (mm, "+" if sgn > 0 else "-", m, mp.nstr(d / sc_, 6)),
          flush=True)

print("\nPART 5  does the direction of Theta(z) in P(W) move with z?", flush=True)
ZS = [("0.31+1.27i", FIT[0]), ("i sqrt2", I * mp.sqrt(2)),
      ("(1+i sqrt7)/2", (1 + I * mp.sqrt(7)) / 2), ("(1+i sqrt11)/2", (1 + I * mp.sqrt(11)) / 2),
      ("2.6i", I * mp.mpf('2.6'))]
base = None
for lbl, z in ZS:
    v = Theta(z)
    nv = mp.sqrt(sum(abs(x)**2 for x in v))
    u = [x / nv for x in v]
    if base is None:
        base = u
        print("   %-16s |Theta| = %-14s (reference direction)"
              % (lbl, mp.nstr(nv, 8)), flush=True)
        continue
    ip = abs(sum(mp.conj(base[l]) * u[l] for l in range(N + 1)))
    print("   %-16s |Theta| = %-14s 1 - |<.,.>| = %s"
          % (lbl, mp.nstr(nv, 8), mp.nstr(1 - ip, 6)), flush=True)
