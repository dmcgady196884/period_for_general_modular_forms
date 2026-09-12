"""On-arc poles, with the mesh clustered WHERE THE POLES ARE. Two jobs.

TRAP THIS FIXES.  shifted_onarc_75.py clustered its quadrature mesh at the PARAMETER endpoints
t = 0, 1 -- the arc's ends -- while an on-arc pole at angle phi sits at
t = (2pi/3 - phi)/(pi/3), i.e. mid-interval (t ~ 0.235 and 0.765 for x = 500).  At eps = 1e-3 the
spiral passes 5e-4 from the pole against a panel width of 0.05, so the quadrature never saw it and
returned a smooth wrong number that still looked like it was in W (2e-20).  This is verbatim the
trap recorded for lem:arcdeform ("the pole sat 2e-8 from a mesh of width 0.26 ... returned a
smooth converged number that matched neither diagonal").  x500_discrepancy.py, clustering at the
pole angles, gets 212068.947897 for the spiral where shifted_onarc_75.py reported 312510.344.
NB the ELLIPTIC work is unaffected: there the poles ARE at t = 0, 1, so endpoint clustering was
exactly right.

JOB 1 -- decompose the wiggle/spiral difference at x = 500.  The two contours flip sides at BOTH
poles but in OPPOSITE senses (wiggle is inside at 1.2932 and outside at 1.8484; the spiral is the
reverse), so a common winding sign cannot be right -- and indeed X = +1 at both missed by relative
1.55, with Delta Phi = 5.15e-7 - 2.62e-6 i complex while 2 pi i sum_p Res_p f = 5.836e-6 is real.
Here the residues are taken at z and Sz SEPARATELY and all four (X_z, X_Sz) in {+-1}^2 are tested:
    Delta Phi   = 2 pi i [ X_z a_z + X_Sz a_Sz ],
    Delta r[l]  = (2 pi i)^{n+1} (-1)^l C(n,l) * 2 pi i *
                  sum_p X_p [ Res_p(f tau^{s-1}) + Res_p(f ktil(.,l+1)) ],
    Delta hat r = Delta r - Delta Phi * r_{E_k}.

JOB 2 -- redo 75 deg (DAM's angle) resolved, since that claim rested on the broken mesh.

PREDICTIONS.
 (1) (X_z, X_Sz) = (+1,-1) or (-1,+1) reproduces the measured difference to harness precision;
     the two common-sign combinations do not.
 (2) At 75 deg both contours are in W at ~1e-28 and DIFFER -- lem:arcon's free side choice at a
     non-elliptic on-arc pole, where the trivial stabiliser gives no 2pi/n to weight and no
     endpoint pinning to define m against.  The earlier |hat r| = 286185.562742 is suspect.
 (3) The wiggle reproduces 312512.434941 again (it did with pole-clustering, and its poles sit
     0.16 off the contour so it was never mesh-sensitive).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, Delta, jay, ktil, E12,                   # noqa: E402
                    rvec, inW, check_orientation)

mp.mp.dps = 30
K, N = 12, 10
EPS = mp.mpf('1e-1')


def nodes(a, b, n, depth, centres):
    xs = [a + (b - a) * mp.mpf(j) / n for j in range(n + 1)]
    for c in centres:
        d = abs(b - a) / n
        for _ in range(depth):
            d /= 2
            xs += [c - d, c + d]
        xs.append(c)
    return sorted(set(x for x in xs if min(a, b) <= x <= max(a, b)), reverse=True)


def wig(th):
    return mp.e**(mp.mpf('0.15') * mp.sin(6 * (th - pi / 2)) + I * th)


def dwig(th):
    return wig(th) * (mp.mpf('0.15') * 6 * mp.cos(6 * (th - pi / 2)) + I)


def spi(th):
    L = mp.log(1 + EPS)
    return mp.e**(L * (pi / 2 - th) / (pi / 6) + I * th)


def dspi(th):
    L = mp.log(1 + EPS)
    return spi(th) * (-L / (pi / 6) + I)


def pint(g, which, th):
    P, dP = (wig, dwig) if which == 'wiggle' else (spi, dspi)
    ns = nodes(2 * pi / 3, pi / 3, 24, 12, th)
    return sum(mp.quad(lambda x: g(P(x)) * dP(x), [u, v]) for u, v in zip(ns[:-1], ns[1:]))


def hatr(f, which, th, rE):
    out = []
    for l in range(N + 1):
        s = mp.mpf(l + 1)
        L = pint(lambda t: f(t) * t**(s - 1), which, th) + \
            pint(lambda t: f(t) * ktil(t, s, K), which, th)
        out.append((2 * pi * I)**(N + 1) * (-1)**l * mp.binomial(N, l) * L)
    Ph = pint(f, which, th)
    return [a - Ph * b for a, b in zip(out, rE)], Ph


def residue(g, c, r=mp.mpf('0.05'), n=48):
    ps = [2 * pi * mp.mpf(j) / n for j in range(n + 1)]
    return sum(mp.quad(lambda p: g(c + r * mp.e**(I * p)) * I * r * mp.e**(I * p), [u, v])
               for u, v in zip(ps[:-1], ps[1:])) / (2 * I * pi)


def theta_of(x):
    lo, hi = pi / 3, pi / 2
    for _ in range(200):
        mid = (lo + hi) / 2
        if mp.re(jay(mp.e**(I * mid))) < x:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


print("guard: Phi(E_12) = %s" % mp.nstr(check_orientation(), 12), flush=True)
rE = rvec(E12, K)

CASES = []
t500 = theta_of(mp.mpf(500))
CASES.append(("x = 500", mp.mpf(500), t500))
t75 = 5 * pi / 12
CASES.append(("75 deg", mp.re(jay(mp.e**(I * t75))), t75))

for lbl, x, t1 in CASES:
    print("=" * 78, flush=True)
    z, Sz = mp.e**(I * t1), mp.e**(I * (pi - t1))
    th = sorted([t1, pi - t1])
    print("%s : x = %s, arg z = %s, arg Sz = %s"
          % (lbl, mp.nstr(x, 12), mp.nstr(t1, 10), mp.nstr(pi - t1, 10)), flush=True)
    f = lambda t, x=x: Delta(t) / (jay(t) - x)
    got = {}
    for which in ('wiggle', 'spiral'):
        v, Ph = hatr(f, which, th, rE)
        got[which] = (v, Ph)
        s_, u_ = inW(v, K)
        print("   %-7s |hat r| = %-18s Phi = %-24s |(1+S)| = %-11s |(1+U+U^2)| = %s"
              % (which, mp.nstr(max(abs(y) for y in v), 12), mp.nstr(Ph, 10),
                 mp.nstr(s_, 5), mp.nstr(u_, 5)), flush=True)
    meas = [a - b for a, b in zip(got['wiggle'][0], got['spiral'][0])]
    az, aSz = residue(f, z), residue(f, Sz)
    print("   a_z = %-26s a_Sz = %s" % (mp.nstr(az, 10), mp.nstr(aSz, 10)), flush=True)
    print("   Phi(wiggle)-Phi(spiral) = %s"
          % mp.nstr(got['wiggle'][1] - got['spiral'][1], 12), flush=True)
    RS = []
    for l in range(N + 1):
        s = mp.mpf(l + 1)
        RS.append((residue(lambda t, s=s: f(t) * (t**(s - 1) + ktil(t, s, K)), z),
                   residue(lambda t, s=s: f(t) * (t**(s - 1) + ktil(t, s, K)), Sz)))
    for Xz in (1, -1):
        for Xs in (1, -1):
            dPh = 2 * I * pi * (Xz * az + Xs * aSz)
            pred = [(2 * pi * I)**(N + 1) * (-1)**l * mp.binomial(N, l) * 2 * I * pi
                    * (Xz * RS[l][0] + Xs * RS[l][1]) for l in range(N + 1)]
            dh = [a - dPh * b for a, b in zip(pred, rE)]
            rel = max(abs(a - b) for a, b in zip(meas, dh)) / max(abs(a) for a in meas)
            print("   (X_z,X_Sz) = (%+d,%+d) : relative = %s" % (Xz, Xs, mp.nstr(rel, 8)),
                  flush=True)
