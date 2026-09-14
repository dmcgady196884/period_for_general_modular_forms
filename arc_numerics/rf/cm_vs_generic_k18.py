"""Item D at k=18: are the W-coordinates of hat r_{f_z} arithmetic at CM points?

FORM.  f_z = g j'/(j - j(z)) with g = Delta E_4 in M_16, so f_z has weight 16+2 = 18
(j'/(j-j(z)) has weight 2).  dim S_18 = 1, so dim W = 2 dim S_k + 1 = 3 and there is a genuine
W^+ / W^- split -- unlike Delta j'/(j-j(z)), which sits at k=14 where dim W = 1.  Taking g
CUSPIDAL also gives c_{f_z}(0) = 0, since j'/(j-j(z)) = -2 pi i sum_{n>=0} j_n(z) q^n starts at
j_0 = 1 and Delta E_4 starts at q.

ANALYTIC INPUT (derived, and the run must reproduce it):
  * Res_p f_z = g(p) at EVERY pole p in SL2(Z).z.  For weight k and a simple pole of residue a
    at z, modularity gives Res_{gamma z} f = a (cz+d)^{k-2}, and g in M_{k-2} transforms with
    exactly that factor.  This is why g must have weight k-2, i.e. why this f and not
    Delta/(j-x).
  * At an ELLIPTIC point the residue is n_p g(p) with n_p the stabiliser order: j - 1728 has a
    double zero at i against a simple zero of j', giving 2/(tau - i); at rho, triple against
    double, giving 3/(tau - rho).  Same 1/n as def:arczero, here as a multiplicity.  Hence
    d = 3, 4 are the hard cases (poles ON the arc, needing the W = 0 prescription) and are left
    out of this run.
  * |hat r_{f_z}| ~ |g(z)| ~ e^{-2 pi Im z}.  The k=12 run already shows this: 2.18e11 at
    z=1.2i against 1.40e9 at z=2i is a factor 156, versus e^{2 pi (0.8)} = 153.
  * z in F lies ABOVE the arc, not in the lens, so eps_z = 0 and no crossing term enters.

THE TEST.  By Chowla-Selberg g(alpha_d)/Omega_d^{16} is algebraic, so hat r is
(algebraic) x Omega_d^{16} x (transcendental content of the L* machinery).  The scale-free
coordinate RATIOS kill the common factor, so they are algebraic iff that remaining content
factors out uniformly across the three coordinates -- not guaranteed.

CONTROLS ARE THE POINT.  Every h(d)=1 CM point has Re z in {0, 1/2}, hence lies on a
reflection-fixed locus, hence has real coordinates automatically.  Reality is therefore a
REFLECTION signal carrying no CM information, and without non-CM points on the same loci
"CM points give real coordinates" would look like a finding and be vacuous.  So: CM points,
non-CM points on the same loci, and genuinely generic points, in one run.

W-basis is built from nullspace(stacked (1+S), (1+U+U^2)) at n = 16 and its dimension asserted;
common.decompose12 is k=12-specific and its read-off convention is unverified, so it is not used.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, E4, E6, Delta, jay, ktil, arcint, rvec,   # noqa: E402
                    NS, NU, opmat, nullspace, inW, check_orientation)

mp.mp.dps = 30
K, N = 18, 16
DEPTH = 10
g = lambda t: Delta(t) * E4(t)                     # M_16, cuspidal
ref = lambda t: E4(t)**3 * E6(t)                   # M_18, c_0 = 1


def Wbasis():
    A, B = opmat(lambda v: NS(v, K), K), opmat(lambda v: NU(v, K), K)
    M = mp.matrix(2 * (N + 1), N + 1)
    for r in range(N + 1):
        for c in range(N + 1):
            M[r, c], M[r + N + 1, c] = A[r, c], B[r, c]
    return nullspace(M, mp.mpf('1e-18'))


def coords(v, B):
    G = mp.matrix(len(B), len(B))
    rhs = mp.matrix(len(B), 1)
    for a in range(len(B)):
        for b in range(len(B)):
            G[a, b] = sum(mp.conj(B[a][l]) * B[b][l] for l in range(N + 1))
        rhs[a] = sum(mp.conj(B[a][l]) * v[l] for l in range(N + 1))
    c = mp.lu_solve(G, rhs)
    rec = [sum(c[b] * B[b][l] for b in range(len(B))) for l in range(N + 1)]
    err = max(abs(a - b) for a, b in zip(v, rec)) / max(abs(a) for a in v)
    return [c[b] for b in range(len(B))], err


rE = rvec(ref, K, DEPTH)


def rhat(z):
    x = jay(z)
    f = lambda t: g(t) * (-2 * I * pi * jay(t) * E6(t) / E4(t)) / (jay(t) - x)
    L = lambda s: mp.e**(-I * pi * s / 2) * arcint(
        lambda t: f(t) * (t**(s - 1) + ktil(t, s, K)), DEPTH)
    rf = [(2 * pi * I)**(N + 1) * (-1)**l * mp.binomial(N, l) * I**(l + 1)
          * L(mp.mpf(l + 1)) for l in range(N + 1)]
    Ph = arcint(f, DEPTH)
    return [a - Ph * b for a, b in zip(rf, rE)], Ph, f


def try_pslq(x, tag):
    for part, nm in ((mp.re(x), "Re"), (mp.im(x), "Im")):
        if abs(part) < mp.mpf('1e-20'):
            continue
        for D in (1, 2, 3):
            v = [part**e for e in range(D + 1)]
            r = mp.pslq(v, tol=mp.mpf('1e-22'), maxcoeff=10**6, maxsteps=20000)
            if r:
                print("        PSLQ %s %s deg%d: %s" % (tag, nm, D, r), flush=True)
                break


print("guard: Phi(E_12) = %s" % mp.nstr(check_orientation(), 12), flush=True)
B = Wbasis()
print("dim W at k=18: %d  (expected 3 = 2 dim S_18 + 1)" % len(B), flush=True)

s7 = (1 + I * mp.sqrt(7)) / 2
s11 = (1 + I * mp.sqrt(11)) / 2
s19 = (1 + I * mp.sqrt(19)) / 2
ZS = [("CM  d=7   (1+i sqrt7)/2", s7), ("CM  d=8   i sqrt2", I * mp.sqrt(2)),
      ("CM  d=11  (1+i sqrt11)/2", s11), ("CM  d=19  (1+i sqrt19)/2", s19),
      ("ctrl  imag axis  1.3i", mp.mpf('1.3') * I),
      ("ctrl  vert edge  0.5+1.4i", mp.mpf('0.5') + mp.mpf('1.4') * I),
      ("generic  0.2+1.3i", mp.mpf('0.2') + mp.mpf('1.3') * I),
      ("generic  0.35+1.5i", mp.mpf('0.35') + mp.mpf('1.5') * I)]

for lbl, z in ZS:
    v, Ph, f = rhat(z)
    c, err = coords(v, B)
    s_, u_ = inW(v, K)
    gz = g(z)
    print("=" * 80, flush=True)
    print("%-30s Im z = %-10s  g(z) = %s" % (lbl, mp.nstr(mp.im(z), 6), mp.nstr(gz, 8)),
          flush=True)
    print("   in W: %-10s %-10s   basis reconstruction err = %-10s   Phi = %s"
          % (mp.nstr(s_, 4), mp.nstr(u_, 4), mp.nstr(err, 4), mp.nstr(Ph, 8)), flush=True)
    print("   coords          %s" % "  ".join(mp.nstr(x, 10) for x in c), flush=True)
    if abs(c[2]) > mp.mpf('1e-60'):
        r1, r2 = c[0] / c[2], c[1] / c[2]
        print("   ratios c0/c2 = %-30s c1/c2 = %s" % (mp.nstr(r1, 14), mp.nstr(r2, 14)),
              flush=True)
        print("   reality: |Im|/|.| = %-14s %-14s"
              % (mp.nstr(abs(mp.im(r1)) / abs(r1), 5), mp.nstr(abs(mp.im(r2)) / abs(r2), 5)),
              flush=True)
        try_pslq(r1, "c0/c2")
        try_pslq(r2, "c1/c2")
