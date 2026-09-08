"""Poles ON the arc, away from the elliptic points: which element of W, and why.

DAM's point: the choice of deformation is NOT cosmetic.  Which residues are enclosed is
exactly what fixing a reference class is for.  Both S-symmetric deformations put
tilde r_f in W, but at DIFFERENT points of W, and the difference is the content.

Setup.  j maps gamma^arc onto [0, 1728], so a real j0 in (0,1728) puts poles ON the arc,
at e^{i phi} and its S-image e^{i(pi - phi)} (lem:arcpairs / the S-pair prose).  Take
f = Delta/(j - j0) at k = 12, where dim S_12 = 1 and dim W = 3, so W is not degenerate.

Contours.  tau(theta) = r(theta) e^{i theta} with log r = eta sin(6(theta - pi/2)), which
is odd about pi/2 and hence S-equivariant.  eta > 0 and eta < 0 are the two admissible
choices: each passes outside one member of the S-pair and inside the other.

Questions.
 (1) Do both land in W, and how far apart are they?
 (2) Does their difference match lem:wall_arc -- i.e. 2 pi i times the residues of the two
     crossed poles, with winding +-1?  This is the real check: both being in W is
     automatic once each is, but matching the residue prediction is not.
 (3) Does the DEFORMATION prescription (def:arcsplit, j0 -> j0 +- i eps) reproduce the
     corresponding one-sided contour as eps -> 0?  Off the elliptic points nothing is
     pinched, so kappa should be 0 and the limit should be finite and equal to a
     one-sided value.
 (4) Same at higher pole order, f = Delta/(j - j0)^2.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common import (mp, I, pi, E4, Delta, jay, ktil, check_orientation, inW, E12,
                    decompose12, W12_minus, W12_plus, P0_12)

mp.mp.dps = 25
KK, NN = 12, 10
J0 = mp.mpf(500)                      # real, in (0,1728): poles land on the arc
ETA = mp.mpf('0.15')


def deformed_arc(g, eta, panels=40):
    """int over tau(theta) = r e^{i theta}, log r = eta sin(6(theta-pi/2)), rho -> rho+1."""
    def integrand(th):
        u = th - pi / 2
        r = mp.e**(eta * mp.sin(6 * u))
        drdth = r * eta * 6 * mp.cos(6 * u)
        tau = r * mp.e**(I * th)
        dtau = (drdth + I * r) * mp.e**(I * th)
        return g(tau) * dtau
    a, b = pi / 3, 2 * pi / 3
    th = [a + (b - a) * m / panels for m in range(panels + 1)]
    return -sum(mp.quad(integrand, [u, v]) for u, v in zip(th[:-1], th[1:]))


def rvec_def(g, eta):
    out = []
    for l in range(NN + 1):
        s = mp.mpf(l + 1)
        out.append((2 * pi * I)**(NN + 1) * (-1)**l * mp.binomial(NN, l)
                   * (deformed_arc(lambda t: g(t) * t**(s - 1), eta)
                      + deformed_arc(lambda t: g(t) * ktil(t, s, KK), eta)))
    return out


def tilde(g, eta):
    Phi = deformed_arc(g, eta)
    rE = rvec_def(E12, eta)
    return [a - Phi * b for a, b in zip(rvec_def(g, eta), rE)], Phi


print("Phi(E_12) on the undeformed arc =", mp.nstr(check_orientation(), 10), flush=True)

# locate the S-pair
zp = mp.findroot(lambda z: jay(z) - J0, mp.e**(I * mp.mpf('1.85')))
zm = -1 / zp
print("S-pair on the arc:  z = %s  (|z| = %s),   Sz = %s  (|Sz| = %s)"
      % (mp.nstr(zp, 10), mp.nstr(abs(zp), 6), mp.nstr(zm, 10), mp.nstr(abs(zm), 6)))
print("   arg z = %s,  arg Sz = %s,  sum = pi ? %s"
      % (mp.nstr(mp.arg(zp), 8), mp.nstr(mp.arg(zm), 8),
         mp.nstr(mp.arg(zp) + mp.arg(zm), 8)), flush=True)

f = lambda t: Delta(t) / (jay(t) - J0)

print()
print("=" * 74)
print("(1) the two admissible S-symmetric deformations")
print("=" * 74)
tA, phiA = tilde(f, ETA)
tB, phiB = tilde(f, -ETA)
for nm, t, ph in (("eta>0", tA, phiA), ("eta<0", tB, phiB)):
    s_, u_ = inW(t, KK)
    print("  %s  Phi = %-24s |tilde r| = %-16s in W: %s / %s"
          % (nm, mp.nstr(ph, 10), mp.nstr(max(abs(x) for x in t), 9),
             mp.nstr(s_, 4), mp.nstr(u_, 4)), flush=True)
D = [a - b for a, b in zip(tA, tB)]
print("  |difference| / |tilde r| = %s"
      % mp.nstr(max(abs(x) for x in D) / max(abs(x) for x in tA), 6))
print("  difference in W: %s / %s" % tuple(mp.nstr(x, 5) for x in inW(D, KK)))
a_, al_, be_ = decompose12(D)
print("  difference in the W-basis:  a = %s,  alpha = %s,  beta = %s"
      % (mp.nstr(a_, 12), mp.nstr(al_, 12), mp.nstr(be_, 12)), flush=True)

print()
print("=" * 74)
print("(2) does the difference match lem:wall_arc?")
print("=" * 74)
# residue of f at a simple pole p:  Delta(p) / j'(p)
def jprime(z, h=mp.mpf('1e-10')):
    return (jay(z + h) - jay(z - h)) / (2 * h)

for p, nm in ((zp, "z"), (zm, "Sz")):
    print("  Res_%s f = %s" % (nm, mp.nstr(Delta(p) / jprime(p), 12)))

def wall_pred(Xp, Xm):
    """2 pi i [ X_p (r_S + r_T)(z) + X_m (r_S + r_T)(Sz) ], assembled into V_n"""
    out = []
    for l in range(NN + 1):
        s = mp.mpf(l + 1)
        tot = mp.mpc(0)
        for p, X in ((zp, Xp), (zm, Xm)):
            R = Delta(p) / jprime(p)
            tot += X * R * (p**(s - 1) + ktil(p, s, KK))
        out.append((2 * pi * I)**(NN + 1) * (-1)**l * mp.binomial(NN, l) * 2 * pi * I * tot)
    return out

nD = max(abs(x) for x in D)
best = None
for Xp in (-1, 0, 1):
    for Xm in (-1, 0, 1):
        if Xp == 0 and Xm == 0:
            continue
        P = wall_pred(Xp, Xm)
        rel = max(abs(a - b) for a, b in zip(D, P)) / nD
        if best is None or rel < best[0]:
            best = (rel, Xp, Xm)
        print("  X_z=%+d X_Sz=%+d :  max|D - pred| / |D| = %s" % (Xp, Xm, mp.nstr(rel, 6)),
              flush=True)
print("  best: X_z=%+d X_Sz=%+d at %s" % (best[1], best[2], mp.nstr(best[0], 6)))

print()
print("=" * 74)
print("(3) deformation prescription: j0 -> j0 + i eps, does it reproduce a one-sided answer?")
print("=" * 74)
from common import arcint, rvec
for sgn, lbl in ((1, "+i eps"), (-1, "-i eps")):
    prev = None
    for e in ('1e-2', '1e-3', '1e-4'):
        j0 = J0 + sgn * I * mp.mpf(e)
        g = lambda t, j0=j0: Delta(t) / (jay(t) - j0)
        Phi = arcint(g)
        tr = [a - Phi * b for a, b in zip(rvec(g, KK), rvec(E12, KK))]
        s_, u_ = inW(tr, KK)
        dA = max(abs(a - b) for a, b in zip(tr, tA)) / max(abs(x) for x in tA)
        dB = max(abs(a - b) for a, b in zip(tr, tB)) / max(abs(x) for x in tB)
        print("  %s |eps|=%-5s |tilde r|=%-14s inW %s/%s   dist to eta>0: %-11s to eta<0: %s"
              % (lbl, e, mp.nstr(max(abs(x) for x in tr), 8), mp.nstr(s_, 3),
                 mp.nstr(u_, 3), mp.nstr(dA, 5), mp.nstr(dB, 5)), flush=True)
