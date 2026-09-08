"""Conjecture C / C' and Proposition B: the eight admissible pairs at a pole on the arc.

SETUP.  k = 12, n = 10, f = Delta/(j - x) with x = 500 in (0,1728).  On gamma^arc the
Hauptmodul j is real and runs 0 -> 1728 -> 0, so it is 2-to-1 onto [0,1728] except at i.
Hence f has exactly TWO poles ON gamma^arc (not on its Gamma-orbit, which carries all of
them), both simple, at z and Sz with arg z + arg Sz = pi.

CLASSES.  Among pairs admissible per def:segments we restrict to those whose S-segment is
S-symmetric, S gamma^S = gamma^S as sets with reversed orientation; this is a restriction
we impose, not part of def:segments, and it is what makes the S^2-loop degenerate.  In
polar form S acts by (r,theta) -> (1/r, pi-theta), so S-symmetry means log r is ODD about
pi/2, and a side at z FORCES the opposite side at Sz.  gamma^T carries no such condition.
Classifying by sides at z and Sz only, with no winding about other poles:

    sigma in {+-1}   side of gamma^S at z (+1 = |tau| > 1); side at Sz is -sigma
    (a,b) in {+-1}^2 sides of gamma^T at z and at Sz, independent
                                                              --> eight pairs

CONJECTURE C.   hat r_f|(1+S) = 0 AND hat r_f|(1+U+U^2) = 0  <=>  (a,b) = (sigma,-sigma).
    NB the (1+S) half is NOT assumed: r_f is assembled from L^*, which contains both
    segment integrals, so r_f|(1+S) need not be blind to gamma^T.  Only the two diagonal
    pairs have been checked before.

CONJECTURE C'.  For (a,b) != (sigma,-sigma) the defect is accounted for by the crossings:
    r_T-residues at the poles where gamma^T and gamma^S differ, PLUS the induced shift
    -(Delta Phi) r_{E_k}.  Only gamma^T moves between these classes, so only the ktil
    kernel contributes residues -- not tau^{s-1}.  The Phi term matters: Phi is the
    T-segment integral, it changes with (a,b), and r_{E_k}|(1+U+U^2) is NOT zero.  Omitting
    it is what left a 3.6e-4 residual in the earlier wall-crossing fit.

PROPOSITION B.  With both segments on gamma^arc and j0 = x +- i eps, the limits as
    eps -> 0 are the two DIAGONAL pairs, one for each sign.
    Resolution: |z_eps - z| ~ eps/|j'(z)| with |j'(z)| ~ 5e3, so the mesh must be clustered
    at z and Sz below eps/5000.  Each eps is run at two clustering depths so that
    unresolved quadrature shows up as a discrepancy instead of a plausible wrong number.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common import mp, I, pi, Delta, jay, ktil, check_orientation, inW, E12, rvec

mp.mp.dps = 25
KK, NN = 12, 10
X0 = mp.mpf(500)
AMP = mp.mpf('0.15')
WID = mp.mpf('0.15')
PART = (sys.argv[1] if len(sys.argv) > 1 else "all")

print("Phi(E_12) on the plain arc =", mp.nstr(check_orientation(), 10), flush=True)
zp = mp.findroot(lambda z: jay(z) - X0, mp.e**(I * mp.mpf('1.85')))
zm = -1 / zp
PHIP, PHIM = mp.arg(zp), mp.arg(zm)
print("poles: arg z = %s, arg Sz = %s, sum = %s"
      % (mp.nstr(PHIP, 8), mp.nstr(PHIM, 8), mp.nstr(PHIP + PHIM, 8)), flush=True)


def bump(th, c):
    u = (th - c) / WID
    return mp.mpf(0) if abs(u) >= 1 else mp.e * mp.e**(-1 / (1 - u**2))


def dbump(th, c):
    u = (th - c) / WID
    if abs(u) >= 1:
        return mp.mpf(0)
    return bump(th, c) * (-2 * u / (1 - u**2)**2) / WID


def polar_int(g, logr, dlogr, panels=64):
    def integrand(th):
        r = mp.e**logr(th)
        return g(r * mp.e**(I * th)) * r * (dlogr(th) + I) * mp.e**(I * th)
    a, b = pi / 3, 2 * pi / 3
    nodes = sorted(set([a, b, PHIP, PHIM]
                       + [a + (b - a) * m / panels for m in range(panels + 1)]))
    return -sum(mp.quad(integrand, [u, v]) for u, v in zip(nodes[:-1], nodes[1:]))


f = lambda t: Delta(t) / (jay(t) - X0)
rE = rvec(E12, KK)                       # E_12 holomorphic: contour-independent
dE = None


def hat_r(sig, a, b):
    lS = lambda th: sig * AMP * mp.sin(6 * (th - pi / 2))
    dS = lambda th: sig * AMP * 6 * mp.cos(6 * (th - pi / 2))
    lT = lambda th: a * AMP * bump(th, PHIP) + b * AMP * bump(th, PHIM)
    dT = lambda th: a * AMP * dbump(th, PHIP) + b * AMP * dbump(th, PHIM)
    Phi = polar_int(f, lT, dT)
    out = []
    for l in range(NN + 1):
        s = mp.mpf(l + 1)
        L = (polar_int(lambda t: f(t) * t**(s - 1), lS, dS)
             + polar_int(lambda t: f(t) * ktil(t, s, KK), lT, dT))
        out.append((2 * pi * I)**(NN + 1) * (-1)**l * mp.binomial(NN, l) * L)
    return [x - Phi * y for x, y in zip(out, rE)], Phi


store, phis = {}, {}
if PART in ("all", "C"):
    print()
    print("=" * 78)
    print("CONJECTURE C: eight pairs.  diagonal means (a,b) = (sigma,-sigma).")
    print("=" * 78)
    print(" sig   a     b    diag   |hat r|          (1+S)         (1+U+U^2)")
    for sig in (1, -1):
        for a in (1, -1):
            for b in (1, -1):
                t, Ph = hat_r(sig, a, b)
                s_, u_ = inW(t, KK)
                store[(sig, a, b)] = t
                phis[(sig, a, b)] = Ph
                print("  %+d   %+d    %+d    %-5s  %-16s %-13s %s"
                      % (sig, a, b, "YES" if (a == sig and b == -sig) else "no",
                         mp.nstr(max(abs(x) for x in t), 9),
                         mp.nstr(s_, 4), mp.nstr(u_, 4)), flush=True)

if PART in ("all", "C") and store:
    print()
    print("=" * 78)
    print("CONJECTURE C': are the six defects r_T-residues plus the Phi shift?")
    print("=" * 78)
    dEU = None

    def NU(v):
        from common import slash, Um
        p = slash(v, NN, Um)
        q = slash(p, NN, Um)
        return [x + y + z for x, y, z in zip(v, p, q)]

    def rT_poly(p, mult):
        """(2 pi i)^{n+2} * mult * Res_p[f ktil] assembled into V_n"""
        h = mp.mpf('1e-12')
        Res = Delta(p) / ((jay(p + h) - jay(p - h)) / (2 * h))
        out = []
        for l in range(NN + 1):
            s = mp.mpf(l + 1)
            out.append((2 * pi * I)**(NN + 1) * (-1)**l * mp.binomial(NN, l)
                       * 2 * pi * I * mult * Res * ktil(p, s, KK))
        return out

    ref = store[(1, 1, -1)]
    refPhi = phis[(1, 1, -1)]
    for key, t in store.items():
        sig, a, b = key
        if key == (1, 1, -1):
            continue
        G = [x - y for x, y in zip(t, ref)]
        nG = max(abs(x) for x in G)
        best = None
        for mz in (-2, -1, 0, 1, 2):
            for mm in (-2, -1, 0, 1, 2):
                P = [p + q for p, q in zip(rT_poly(zp, mz), rT_poly(zm, mm))]
                P = [p - (phis[key] - refPhi) * e for p, e in zip(P, rE)]
                rel = max(abs(x - y) for x, y in zip(G, P)) / nG
                if best is None or rel < best[0]:
                    best = (rel, mz, mm)
        print("  (%+d,%+d,%+d) - (+1,+1,-1):  |G| = %-16s best m_z=%+d m_Sz=%+d at %s"
              % (sig, a, b, mp.nstr(nG, 9), best[1], best[2], mp.nstr(best[0], 6)),
              flush=True)

if PART in ("all", "B"):
    print()
    print("=" * 78)
    print("PROPOSITION B: pole motion on the plain arc, two clustering depths per eps")
    print("=" * 78)

    def nodes_at(depth):
        a, b = pi / 3, 2 * pi / 3
        ns = [a, b, (a + b) / 2]
        for c in (PHIP, PHIM):
            d = mp.mpf('0.2')
            for _ in range(depth):
                d /= 2
                ns += [c - d, c + d]
            ns.append(c)
        return sorted(set([x for x in ns if a <= x <= b]))

    def plain(g, depth):
        ns = nodes_at(depth)
        return -sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
                    for u, v in zip(ns[:-1], ns[1:]))

    if not store:
        for sig, a, b in ((1, 1, -1), (-1, -1, 1)):
            store[(sig, a, b)], phis[(sig, a, b)] = hat_r(sig, a, b)
    dA, dB = store[(1, 1, -1)], store[(-1, -1, 1)]
    print("  diagonals: |A| = %s   |B| = %s"
          % (mp.nstr(max(abs(x) for x in dA), 9), mp.nstr(max(abs(x) for x in dB), 9)))
    for sgn, lbl in ((1, "+"), (-1, "-")):
        for e in ('1e-2', '1e-3', '1e-4'):
            row = []
            for depth in (18, 26):
                j0 = X0 + sgn * I * mp.mpf(e)
                g = lambda t, j0=j0: Delta(t) / (jay(t) - j0)
                Ph = plain(g, depth)
                rf = []
                for l in range(NN + 1):
                    s = mp.mpf(l + 1)
                    rf.append((2 * pi * I)**(NN + 1) * (-1)**l * mp.binomial(NN, l)
                              * (plain(lambda t: g(t) * t**(s - 1), depth)
                                 + plain(lambda t: g(t) * ktil(t, s, KK), depth)))
                row.append([x - Ph * y for x, y in zip(rf, rE)])
            conv = max(abs(x - y) for x, y in zip(*row)) / max(abs(x) for x in row[1])
            tr = row[1]
            qA = max(abs(x - y) for x, y in zip(tr, dA)) / max(abs(x) for x in dA)
            qB = max(abs(x - y) for x, y in zip(tr, dB)) / max(abs(x) for x in dB)
            print("  %s eps=%-5s depth18-vs-26: %-11s inW %s/%s  to A: %-11s to B: %s"
                  % (lbl, e, mp.nstr(conv, 4),
                     *[mp.nstr(x, 3) for x in inW(tr, KK)],
                     mp.nstr(qA, 5), mp.nstr(qB, 5)), flush=True)
