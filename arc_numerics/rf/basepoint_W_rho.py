"""Does W-membership select the connector routing at rho?

basepoint_shift_rho.py: displacing tau_0 off rho+1 makes L* finite and exactly flat in eps even
when f has a pole at SL2(Z).rho, but hat gamma^T then needs a connector around rho and the two
routings (short 2pi/3, long 4pi/3) give different L*, differing by one lem:wall_arc crossing
(X_S = 0, X_T = -1).  Reality does NOT select between them (short is real in all six L* rows, long
in three).  res_at_rho.py: Res_rho f != 0 at every P tested -- 1.52e-4 i (P=2), 1.58e-6 i (P=3),
-3.26e-7 i (P=4), radius-independent -- so eq:arcondefect's defect, which is built from
a_p = Res_{tau=p} f and NOT from Res(f ktil), is nonzero for whichever routing winds.

(An earlier argument that the stabiliser forces Res_rho f = 0 unless P = 1 mod 3 was WRONG: V acts
linearly on w = (tau-rho)/(tau-rhobar), not on u = tau-rho, so the Laurent-support constraint
j = 0 mod n holds in w and mixes orders in u.  Nothing is forced to vanish.)

PREDICTIONS.
 (1) |(1+S)| ~ 0 for BOTH routings.  Gamma_S is S-symmetric by construction (log r odd about
     theta = pi/2, and S acts by (r,theta) -> (1/r, pi-theta)), and lem:arcon's (1+S) argument uses
     only that -- no property of hat gamma^T enters.  So (1+S) cannot discriminate.
 (2) |(1+U+U^2)| ~ 0 for EXACTLY ONE routing, O(1) for the other.  U = TS mixes the segments, and
     by eq:arcondefect the defect is -(2 pi i)^{n+2} nu a_rho (X - rho Y)^n - 2 pi i X_T a_rho D_E,
     non-zero unless D_E lies in span{(X - rho Y)^n} -- which the rank scans exclude at k=12.
 (3) The winner is SHORT.  Weakly held: short is the minimal detour and the one that is real in all
     six L* rows, but I cannot derive it, and "neither" is the outcome that would sink the
     base-point route.
 (4) tilde_r^short != tilde_r^long, since a_rho != 0.

Conventions.  rvec() drops e^{-i pi s/2} because def:rf's i^{l+1} cancels it at s = l+1, so the
segment integrals here are taken RAW.  r_{E_k} is the plain-arc value from common.rvec: E_12 is
holomorphic, and basepoint_shift_rho.py's control showed the shifted class reproduces plain_arc
to 20 digits with short - long identically zero.  eps = 1e-3: flat there, and better conditioned
than 1e-5.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, E4, Delta, jay, ktil, E12, RHO,          # noqa: E402
                    rvec, inW, check_orientation)

mp.mp.dps = 30
TAU0 = mp.e**(I * pi / 3)
EPS = mp.mpf('1e-3')


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


def rvec_shifted(g, k, eps, longway):
    """r_f in the shifted class.  RAW segment integrals, matching common.rvec."""
    S, C = spiral(eps), connector(eps, longway)
    n = k - 2
    out = []
    for l in range(n + 1):
        s = mp.mpf(l + 1)
        LS = integrate(S, lambda t: g(t) * t**(s - 1))
        LT = integrate(C, lambda t: g(t) * ktil(t, s, k)) + \
            integrate(S, lambda t: g(t) * ktil(t, s, k))
        out.append((2 * pi * I)**(n + 1) * (-1)**l * mp.binomial(n, l) * (LS + LT))
        print("      l=%-3d done" % l, flush=True)
    return out


def Phi_shifted(g, eps, longway):
    S, C = spiral(eps), connector(eps, longway)
    return integrate(C, g) + integrate(S, g)


print("orientation guard, Phi(E_12) on the plain arc = %s"
      % mp.nstr(check_orientation(), 12), flush=True)

CASES = [
    ("k=12  Delta/j    P=3   (dim S_12 = 1, dim W = 3)", lambda t: Delta(t) / jay(t), 12),
    ("k=8   E_4^2/j^2  P=4   (dim S_8 = 0)", lambda t: E4(t)**2 / jay(t)**2, 8),
]

for lbl, f, k in CASES:
    print("=" * 76, flush=True)
    print(lbl, flush=True)
    rE = rvec(E12, k) if k == 12 else rvec(lambda t: E4(t)**2, k)
    res = {}
    for lw in (False, True):
        tag = "long" if lw else "short"
        print("   %s routing:" % tag, flush=True)
        rf = rvec_shifted(f, k, EPS, lw)
        Ph = Phi_shifted(f, EPS, lw)
        tr = [a - Ph * b for a, b in zip(rf, rE)]
        s_def, u_def = inW(tr, k)
        res[tag] = tr
        print("      Phi = %s" % mp.nstr(Ph, 12), flush=True)
        print("      |(1+S)|       = %s" % mp.nstr(s_def, 8), flush=True)
        print("      |(1+U+U^2)|   = %s" % mp.nstr(u_def, 8), flush=True)
    d = max(abs(a - b) for a, b in zip(res['short'], res['long']))
    sc = max(abs(x) for x in res['short'])
    print("   |tilde_r^short - tilde_r^long| / scale = %s" % mp.nstr(d / sc, 8), flush=True)
