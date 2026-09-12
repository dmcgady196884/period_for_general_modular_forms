"""Does the connector winding X_T = +1 put hat r_f back in W at rho?

connector_enclosing_check.py measured the two (1+U+U^2) defect VECTORS for f = Delta/j, k=12,
pole of order 3 at rho:
    short (X_T = 0):  |.| = 501060.8453    long (X_T = -1):  |.| = 1002121.691
with IDENTICAL distances to (X - rho Y)^n (0.875822) and to D_E (0.516774), i.e. the two vectors
are parallel and the long is exactly twice the short.  Writing D(X_T) = A + X_T B, short gives A
and long gives A - B = 2A, so B = -A and

    D(X_T) = A (1 - X_T),     which VANISHES at X_T = +1.

That is an integer, so a routing in W should exist -- one unit of winding the OTHER way from
short.  My earlier claim that no winding could work (from eq:arcondefect needing D_E in
span{(X - rho Y)^n}) does not apply: this pole is of order 3, not simple, and the defect lies
along neither line.

The connector is a log-spiral about rho from T^{-1}tau_0 to S tau_0 whose angular travel d is
taken modulo 2 pi; `turns` shifts it by whole turns.  short is d0 in (-pi, pi], long is the
routing connector() produced by pushing d0 across zero, measured to be X_T = -1 (via
Phi^short - Phi^long = -2 pi i a_rho = +9.89994150954e-6).  So X_T = +1 is d0 + 2 pi.

PREDICTIONS.
 (1) |(1+U+U^2)| ~ 1e-20, the harness floor established by the five validation cases -- NOT zero,
     but 25 orders below the 501060 it replaces.
 (2) |(1+S)| ~ 1e-30 as always (Gamma_S is untouched and S-symmetric).
 (3) Phi = +9.89994150954e-6: Phi falls by 2 pi |a_rho| per unit DECREASE in X_T (Phi(0) ~ 0,
     Phi(-1) = -9.89994e-6), so Phi(+1) = +9.89994e-6.
 (4) X_T = +2 should give |D| = 501060.8453 again, by D(X_T) = A(1 - X_T) at X_T = 2 giving -A.
     Included as the falsifier: if the affine law is real, this row is forced.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, Delta, jay, ktil, E12, RHO,              # noqa: E402
                    rvec, inW, NU, check_orientation)

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


def connector(eps, turns):
    """turns = 0 is the short routing; turns shifts the angular travel by whole turns."""
    a = (TAU0 * (1 + eps) - 1) - RHO
    b = RHO / (1 + eps) - RHO
    la, lb = mp.log(a), mp.log(b)
    d = mp.im(lb - la)
    while d > pi:
        d -= 2 * pi
    while d <= -pi:
        d += 2 * pi
    d += 2 * pi * turns
    lb = mp.mpc(mp.re(lb), mp.im(la) + d)

    def p(t):
        w = la + t * (lb - la)
        return RHO + mp.e**w, mp.e**w * (lb - la)
    return p


def run(g, eps, turns, rE):
    S, C = spiral(eps), connector(eps, turns)
    rf = []
    for l in range(N + 1):
        s = mp.mpf(l + 1)
        LS = integrate(S, lambda t: g(t) * t**(s - 1))
        LT = integrate(C, lambda t: g(t) * ktil(t, s, K)) + \
            integrate(S, lambda t: g(t) * ktil(t, s, K))
        rf.append((2 * pi * I)**(N + 1) * (-1)**l * mp.binomial(N, l) * (LS + LT))
    Ph = integrate(C, g) + integrate(S, g)
    return [a - Ph * b for a, b in zip(rf, rE)], Ph


print("guard: Phi(E_12) = %s" % mp.nstr(check_orientation(), 12), flush=True)
f = lambda t: Delta(t) / jay(t)
rE = rvec(E12, K)

# connector() with turns=0 reproduces the earlier "short"; the earlier "long" pushed d across
# zero, which for this geometry is turns = -1.  Sweep both signs and one extra turn.
for turns in (0, -1, 1, 2):
    tr, Ph = run(f, EPS, turns, rE)
    s_, u_ = inW(tr, K)
    d = NU(tr, K)
    print("turns=%-3d  Phi = %-26s |(1+S)| = %-12s |(1+U+U^2)| = %-12s |defect| = %s"
          % (turns, mp.nstr(Ph, 12), mp.nstr(s_, 6), mp.nstr(u_, 6),
             mp.nstr(max(abs(x) for x in d), 12)), flush=True)
