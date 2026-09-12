"""The one untested regime: the connector ENCLOSING a pole.

All five validation rows (7i, danger zone, below arc, 75 deg, x=500) have short == long, because
none of those f has a pole at SL2(Z).rho and so the connector encloses nothing.  Only a pole AT
rho exercises that regime, and there the only evidence so far is Phi^long - Phi^short = 2 pi i
a_rho -- one check, and one that tests winding bookkeeping rather than the r_f assembly.

TEST.  Predict the whole difference VECTOR from independently-computed residues.  lem:wall_arc
with X_S = 0 (Gamma_S is shared) and X_T = -1 gives, in common.rvec's RAW convention (no
e^{-i pi s/2}, since def:rf's i^{l+1} cancels it),

    Delta L*_raw(l+1) = -2 pi i Res_rho( f ktil(.,l+1) ),
    Delta Phi         = -2 pi i a_rho,          a_rho := Res_rho f,

so with hat r = r_f - Phi r_{E_k} and r_f[l] = (2 pi i)^{n+1} (-1)^l C(n,l) L*_raw(l+1),

    hat r^short[l] - hat r^long[l]
        = -(2 pi i)^{n+2} (-1)^l C(n,l) Res_rho( f ktil(.,l+1) )  +  2 pi i a_rho r_{E_k}[l].

Every residue on the right is computed by a small circle about rho -- nothing from the contour
machinery -- so this is an INDEPENDENT prediction, not a self-consistency check.  Radii are swept
so a genuine residue is radius-independent.

SECONDARY.  Report the defect vectors themselves and test whether they lie along
(X - rho Y)^n (coefficients C(n,l)(-rho)^l) and D_E = r_{E_k}|(1+U+U^2), which is the closed form
eq:arcondefect asserts.  If the defect is that shape it is a corrector built from the single
number a_rho, in the manner of thm:rtildeW.

f = Delta/j, k=12, pole of order P=3 at rho (a_rho = 1.57562462756e-6 i, measured, purely
imaginary).  NOTE the pole is NOT simple, so Res_rho(f ktil) involves ktil and its first two
derivatives at rho; the prediction above does not care, but the (X - rho Y)^n shape might.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, Delta, jay, ktil, E12, RHO,              # noqa: E402
                    rvec, NS, NU, check_orientation)

mp.mp.dps = 30
TAU0 = mp.e**(I * pi / 3)
EPS = mp.mpf('1e-3')
K = 12
N = K - 2


def clustered(n, depth):
    xs = [mp.mpf(j) / n for j in range(n + 1)]
    d = mp.mpf(1) / n
    for _ in range(depth):
        d /= 2
        xs += [d, 1 - d]
    return sorted(set(xs))


def integrate(path, g, n=20, depth=13):
    ts = clustered(n, depth)
    return sum(mp.quad(lambda t: g(path(t)[0]) * path(t)[1], [u, v])
               for u, v in zip(ts[:-1], ts[1:]))


def spiral(eps):
    L = mp.log(1 + eps)

    def p(t):
        th = 2 * pi / 3 + t * (pi / 3 - 2 * pi / 3)
        tau = mp.e**(L * (pi / 2 - th) / (pi / 6) + I * th)
        return tau, tau * (-L / (pi / 6) + I) * (pi / 3 - 2 * pi / 3)
    return p


def connector(eps, longway):
    a = (TAU0 * (1 + eps) - 1) - RHO
    b = RHO / (1 + eps) - RHO
    la, lb = mp.log(a), mp.log(b)
    d = mp.im(lb - la)
    while d > pi:
        d -= 2 * pi
    while d <= -pi:
        d += 2 * pi
    if longway:
        d += 2 * pi if d < 0 else -2 * pi
    lb = mp.mpc(mp.re(lb), mp.im(la) + d)

    def p(t):
        w = la + t * (lb - la)
        return RHO + mp.e**w, mp.e**w * (lb - la)
    return p


def shifted(g, eps, longway, rE):
    S, C = spiral(eps), connector(eps, longway)
    rf = []
    for l in range(N + 1):
        s = mp.mpf(l + 1)
        LS = integrate(S, lambda t: g(t) * t**(s - 1))
        LT = integrate(C, lambda t: g(t) * ktil(t, s, K)) + \
            integrate(S, lambda t: g(t) * ktil(t, s, K))
        rf.append((2 * pi * I)**(N + 1) * (-1)**l * mp.binomial(N, l) * (LS + LT))
        print("      l=%-3d" % l, flush=True)
    Ph = integrate(C, g) + integrate(S, g)
    return [a - Ph * b for a, b in zip(rf, rE)], Ph


def residue(g, r, n=60):
    ps = [2 * pi * mp.mpf(j) / n for j in range(n + 1)]
    tot = sum(mp.quad(lambda p: g(RHO + r * mp.e**(I * p)) * I * r * mp.e**(I * p), [u, v])
              for u, v in zip(ps[:-1], ps[1:]))
    return tot / (2 * I * pi)


def align(u, v):
    """relative distance from u to the line spanned by v"""
    nv = sum(abs(x)**2 for x in v)
    c = sum(x * mp.conj(y) for x, y in zip(u, v)) / nv
    num = max(abs(a - c * b) for a, b in zip(u, v))
    return num / max(abs(a) for a in u)


print("guard: Phi(E_12) = %s" % mp.nstr(check_orientation(), 12), flush=True)
f = lambda t: Delta(t) / jay(t)
rE = rvec(E12, K)

print("residues about rho (two radii, must agree):", flush=True)
a_rho = [residue(f, r) for r in (mp.mpf('0.06'), mp.mpf('0.03'))]
print("   a_rho = %s   %s" % (mp.nstr(a_rho[0], 12), mp.nstr(a_rho[1], 12)), flush=True)
RK = []
for l in range(N + 1):
    s = mp.mpf(l + 1)
    v = [residue(lambda t, s=s: f(t) * ktil(t, s, K), r) for r in (mp.mpf('0.06'), mp.mpf('0.03'))]
    RK.append(v[0])
    print("   l=%-3d Res(f ktil) = %-28s  radius drift %s"
          % (l, mp.nstr(v[0], 12), mp.nstr(abs(v[0] - v[1]), 4)), flush=True)

print("short routing:", flush=True)
rs, Phs = shifted(f, EPS, False, rE)
print("long routing:", flush=True)
rl, Phl = shifted(f, EPS, True, rE)

print("=" * 76, flush=True)
print("Phi^short - Phi^long = %-26s   -2 pi i a_rho = %s"
      % (mp.nstr(Phs - Phl, 12), mp.nstr(-2 * pi * I * a_rho[0], 12)), flush=True)

pred = [-(2 * pi * I)**(N + 2) * (-1)**l * mp.binomial(N, l) * RK[l]
        + 2 * pi * I * a_rho[0] * rE[l] for l in range(N + 1)]
meas = [a - b for a, b in zip(rs, rl)]
rel = max(abs(a - b) for a, b in zip(meas, pred)) / max(abs(a) for a in meas)
print("difference vector: predicted vs measured, relative = %s" % mp.nstr(rel, 8), flush=True)
for l in range(N + 1):
    print("   l=%-3d meas %-30s pred %s"
          % (l, mp.nstr(meas[l], 10), mp.nstr(pred[l], 10)), flush=True)

print("=" * 76, flush=True)
XY = [mp.binomial(N, l) * (-RHO)**l for l in range(N + 1)]
DE = NU(rE, K)
for tag, v in (("short", rs), ("long", rl)):
    d = NU(v, K)
    print("%-5s defect: |.| = %-16s  dist to (X-rho Y)^n = %-14s  dist to D_E = %s"
          % (tag, mp.nstr(max(abs(x) for x in d), 10), mp.nstr(align(d, XY), 6),
             mp.nstr(align(d, DE), 6)), flush=True)
print("diff vector: dist to (X-rho Y)^n = %-14s  dist to D_E = %s"
      % (mp.nstr(align(meas, XY), 6), mp.nstr(align(meas, DE), 6)), flush=True)
