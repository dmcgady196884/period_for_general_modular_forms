"""DAM's example, run properly: f_z = Delta E_4 J'/(J - J(z)) in F_18, and the image of Xi.

g = Delta E_4 in M_16 and J'/(J - J(z)) has weight 2, so f_z in F_18.  J = j - 744, so J' = j' =
-2 pi i (J + 744) E_6/E_4, and Res_z f_z = g(z) (Lemma lem:resequiv; already confirmed to 1e-8 via
Phi(f_z) = 2 pi i g(z) in logs/cm18.txt).  Hence the whole object FACTORS,

    Xi_{f_z; z} = g(z) * Theta_18(z) ,

with g(z) = Delta(z)E_4(z) carrying every transcendental ingredient and Theta_18 -- universal,
depending on k alone -- carrying the W-geometry.  prop:thetalaurent says

    Theta(z) = sum_{m=-1}^{n+1} c_m z^m ,  c_m in W ,  c_{n-m} = (-1)^{m+1} c_m ,

which at k=18 (n=16) is 19 coefficient vectors inside a 3-dimensional W.  theta_laurent.py checked
this at k=12 only; this run checks the weight DAM actually asked about.

THE IMAGE OF Xi.  Xi is linear in f, so {Xi_f : f in F_k} is a linear SUBSPACE of W.  If the
Theta(z) span W as z varies then that subspace is all of W, i.e. Xi is surjective: every element of
W is the residue polynomial of some meromorphic form.  Equivalently the curve z -> Theta(z) lies in
no proper subspace.  PART 4 measures the rank directly, over CM points and generic points both,
rather than inferring it from the k=12 coefficient rank.

PART 5 gives the explicit W-coordinates of Xi at the class-number-one CM points, which is the
"pretty analytic form of the result of the residue at z" in the only sense that is basis-independent:
magnitude |g(z)| against e^{-4 pi Im z}, and direction in P(W).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, E4, E6, Delta, ktil, slash, rvec, inW,      # noqa: E402
                    opmat, NS, NU, nullspace, arcint, rank)

mp.mp.dps = 50
K, N = 18, 16
Sm = (0, -1, 1, 0)
g = lambda t: Delta(t) * E4(t)                     # M_16
ref = lambda t: E4(t)**3 * E6(t)                   # M_18, c_0 = 1
print("reference E_4^3E_6  Phi = %s" % mp.nstr(arcint(ref, 11), 12), flush=True)
rE = rvec(ref, K, 11)


def Theta(z):
    P = [mp.binomial(N, l) * (-z)**l for l in range(N + 1)]
    PS = slash(P, N, Sm)
    out = []
    for l in range(N + 1):
        br = (P[l] - PS[l]) + mp.binomial(N, l) * (-1)**l * (
            ktil(z, mp.mpf(l + 1), K) - z**N * ktil(-1 / z, mp.mpf(l + 1), K))
        out.append((2 * pi * I)**(N + 1) * br)
    return [2 * pi * I * (a - (1 - z**N) * b) for a, b in zip(out, rE)]


MS = list(range(-1, N + 2))                        # 19 unknowns
FIT = [mp.mpf(str(re_)) + I * mp.mpf(str(im_)) for re_, im_ in
       [('0.31', '1.27'), ('-0.4', '1.1'), ('0.7', '1.9'), ('0.0', '1.3'), ('0.5', '1.4'),
        ('-0.25', '1.6'), ('0.15', '2.1'), ('0.9', '1.05'), ('-0.6', '1.35'), ('0.05', '1.15'),
        ('0.44', '1.72'), ('-0.8', '2.4'), ('0.33', '1.05'), ('0.22', '2.8'), ('-0.15', '1.45'),
        ('0.66', '1.25'), ('-0.55', '1.95'), ('0.12', '1.62'), ('0.85', '2.2')]]
assert len(FIT) == len(MS)

V = mp.matrix(len(MS), len(MS))
for a, z in enumerate(FIT):
    for b, m in enumerate(MS):
        V[a, b] = z**m
TH = [Theta(z) for z in FIT]
c = []
cols = []
for l in range(N + 1):
    cols.append(mp.lu_solve(V, mp.matrix([TH[a][l] for a in range(len(MS))])))
c = [[cols[l][b] for l in range(N + 1)] for b in range(len(MS))]

print("\nPART 1  out-of-sample check of the Laurent form at k=18", flush=True)
for z in (mp.mpf('0.62') + I * mp.mpf('1.11'), I * mp.sqrt(2), (1 + I * mp.sqrt(7)) / 2):
    got = Theta(z)
    pred = [sum(c[b][l] * z**MS[b] for b in range(len(MS))) for l in range(N + 1)]
    print("   z = %-26s relative error %s"
          % (mp.nstr(z, 8), mp.nstr(max(abs(x - y) for x, y in zip(got, pred))
                                    / max(abs(x) for x in got), 6)), flush=True)

print("\nPART 2  every c_m in W?  and the functional equation c_{n-m} = (-1)^{m+1} c_m", flush=True)
bad = 0
for b, m in enumerate(MS):
    nb = max(abs(x) for x in c[b])
    s_, u_ = inW(c[b], K) if nb > mp.mpf('1e-30') else (mp.mpf(0), mp.mpf(0))
    if max(s_, u_) > mp.mpf('1e-25'):
        bad += 1
    mm = N - m
    fe = ""
    if mm in MS and mm >= m:
        bb = MS.index(mm)
        sgn = (-1)**(m + 1)
        fe = "   c_%d=%sc_%d: %s" % (mm, "+" if sgn > 0 else "-", m,
                                     mp.nstr(max(abs(c[bb][l] - sgn * c[b][l])
                                                 for l in range(N + 1))
                                             / (max(abs(x) for x in c[bb]) or mp.mpf(1)), 4))
    print("   m=%-3d |c_m|=%-14s (1+S)=%-10s (1+U+U^2)=%-10s%s"
          % (m, mp.nstr(nb, 8), mp.nstr(s_, 4), mp.nstr(u_, 4), fe), flush=True)
print("   coefficients failing W-membership at 1e-25: %d" % bad, flush=True)

print("\nPART 3  functional equation directly:  Theta(-1/z) = -z^{-n} Theta(z)", flush=True)
for z in (FIT[0], I * mp.sqrt(2), (1 + I * mp.sqrt(11)) / 2):
    lhs, rhs = Theta(-1 / z), [-x / z**N for x in Theta(z)]
    print("   z = %-26s relative error %s"
          % (mp.nstr(z, 8), mp.nstr(max(abs(a - b) for a, b in zip(lhs, rhs))
                                    / max(abs(x) for x in rhs), 6)), flush=True)

A_, B_ = opmat(lambda v: NS(v, K), K), opmat(lambda v: NU(v, K), K)
MM = mp.matrix(2 * (N + 1), N + 1)
for r_ in range(N + 1):
    for cc in range(N + 1):
        MM[r_, cc], MM[r_ + N + 1, cc] = A_[r_, cc], B_[r_, cc]
WB = nullspace(MM, mp.mpf('1e-30'))

ZS = [("d=7   (1+isq7)/2", (1 + I * mp.sqrt(7)) / 2), ("d=8   i sqrt2", I * mp.sqrt(2)),
      ("d=11  (1+isq11)/2", (1 + I * mp.sqrt(11)) / 2),
      ("d=19  (1+isq19)/2", (1 + I * mp.sqrt(19)) / 2),
      ("d=43  (1+isq43)/2", (1 + I * mp.sqrt(43)) / 2),
      ("generic 0.2+1.3i", mp.mpf('0.2') + I * mp.mpf('1.3')),
      ("generic 0.35+1.5i", mp.mpf('0.35') + I * mp.mpf('1.5'))]

print("\nPART 4  IMAGE of Xi: do the Theta(z) span W?  dim W = %d" % len(WB), flush=True)
G = mp.matrix(len(ZS), N + 1)
for a, (_, z) in enumerate(ZS):
    v = Theta(z)
    nv = max(abs(x) for x in v)
    for l in range(N + 1):
        G[a, l] = v[l] / nv
for tol in ('1e-30', '1e-20', '1e-10', '1e-6'):
    print("   rank of the %d Theta(z) (tol %s) = %d" % (len(ZS), tol, rank(G, mp.mpf(tol))),
          flush=True)
Gc = mp.matrix(len(MS), N + 1)
for b in range(len(MS)):
    nb = max(abs(x) for x in c[b]) or mp.mpf(1)
    for l in range(N + 1):
        Gc[b, l] = c[b][l] / nb
print("   rank of the %d coefficients c_m (tol 1e-20) = %d" % (len(MS), rank(Gc, mp.mpf('1e-20'))),
      flush=True)

print("\nPART 6  localisation in W by parity.  (P|eps)(X,Y) = P(-X,Y) acts as (-1)^l on the", flush=True)
print("        X^{n-l}Y^l coefficient, so W splits into even-l and odd-l parts.  dim S_18 = 1,", flush=True)
print("        so dim W^odd should be 1 and the odd part of Theta is ONE scalar Laurent poly.", flush=True)


def par(v, odd):
    return [x if (l % 2 == (1 if odd else 0)) else mp.mpc(0) for l, x in enumerate(v)]


for tag, od in (("even", False), ("odd", True)):
    Gp = mp.matrix(len(WB), N + 1)
    for a in range(len(WB)):
        pv = par(WB[a], od)
        nb = max(abs(x) for x in pv) or mp.mpf(1)
        for l in range(N + 1):
            Gp[a, l] = pv[l] / nb
    print("   dim W^%-5s = %d" % (tag, rank(Gp, mp.mpf('1e-20'))), flush=True)

wodd = None
for a in range(len(WB)):
    pv = par(WB[a], True)
    if max(abs(x) for x in pv) > mp.mpf('1e-20'):
        wodd = pv
        break
l0 = max(range(N + 1), key=lambda l: abs(wodd[l]))
print("   odd basis vector normalised at l = %d" % l0, flush=True)
lam = [par(c[b], True)[l0] / wodd[l0] for b in range(len(MS))]
print("   lambda_-(z) = sum_m lambda_m z^m,  lambda_m / (2 pi i)^{n+2}, PSLQ'd:", flush=True)
for b, m in enumerate(MS):
    x = lam[b] / (2 * pi * I)**(N + 2)
    if abs(x) < mp.mpf('1e-30'):
        print("      m=%-3d  0" % m, flush=True)
        continue
    pr = mp.re(x)
    q = mp.pslq([pr, mp.mpf(1)], tol=mp.mpf('1e-30'), maxcoeff=10**14, maxsteps=10**5)
    tag = ("%s/%s" % (-q[1], q[0])) if q else mp.nstr(pr, 14)
    im = mp.im(x)
    ex = ("  +%si" % mp.nstr(im, 4)) if abs(im) > mp.mpf('1e-25') * max(abs(pr), mp.mpf(1)) else ""
    print("      m=%-3d  %s%s" % (m, tag, ex), flush=True)
# residual: is Theta^odd really rank one, i.e. exactly lambda_-(z) * wodd?
z = FIT[3]
tv = par(Theta(z), True)
pred = [sum(lam[b] * z**MS[b] for b in range(len(MS))) * x for x in wodd]
print("   Theta^odd(z) = lambda_-(z) w_odd at z=%s : rel %s"
      % (mp.nstr(z, 6), mp.nstr(max(abs(a - b) for a, b in zip(tv, pred))
                                / max(abs(x) for x in tv), 6)), flush=True)

print("\nPART 5  Xi_{f_z;z} = g(z) Theta(z) for f_z = Delta E_4 J'/(J-J(z))", flush=True)
print("   %-20s %-22s %-14s %s" % ("z", "g(z) = Delta E_4 (z)", "|Xi|", "|g| vs e^{-4 pi Im z}"),
      flush=True)
for lbl, z in ZS:
    gz = g(z)
    v = Theta(z)
    xi = [gz * x for x in v]
    s_, u_ = inW(xi, K)
    pred = mp.e**(-4 * pi * mp.im(z))
    print("   %-20s %-22s %-14s ratio %s   inW %s/%s"
          % (lbl, mp.nstr(gz, 8), mp.nstr(max(abs(x) for x in xi), 8),
             mp.nstr(abs(gz) / pred, 8), mp.nstr(s_, 3), mp.nstr(u_, 3)), flush=True)
