"""Is the rho m=0 residual numerical conditioning, or a real nonzero defect?

m_zero_Pk_sweep.py gave these |(1+U+U^2)| at m = 0:
    (i,   P=2, k=12) unshifted detour   2.4e-28
    (i,   P=4, k=12) unshifted detour   9.3e-26
    (rho, P=2, k=4 ) shifted spiral     8.7e-12
    (rho, P=4, k=8 ) shifted spiral     1.9e-2      <-- the outlier
The i control at P=4 is clean, so P mod n is NOT the cause: the prescription does not care about
the pole order.  What the residual tracks is the GEOMETRY.  In the shifted construction the
connector runs at distance eps from a pole of order P, so the near-pole integrand reaches
eps^{-P} = 10^{3P} at eps = 1e-3: six more orders of cancellation at P=4 than at P=2, against ten
observed.  Right direction, right rough size.

THE TEST.  L* is EXACTLY eps-independent (proved -- the tau_0-independence identity -- and measured
flat to 20 digits over four decades in basepoint_shift_rho.py), so eps = 0.1 is as legitimate as
eps = 1e-3 and far better conditioned: eps^{-4} falls from 1e12 to 1e4.

PREDICTIONS.  If the residual is conditioning, it falls steeply with eps:
 (1) (rho, P=4, k=8) at eps = 1e-1 should reach ~1e-11 or better, a gain of ~9 orders on the
     1.9e-2 at eps = 1e-3; at eps = 1e-2 something intermediate, ~1e-8.
 (2) (rho, P=2, k=4) at eps = 1e-1 should improve on 8.7e-12 by ~6 orders, to ~1e-18, approaching
     the harness floor the off-contour validations sat at.
 (3) The realisable windings' defects |D(-1)| and |D(2)| must be eps-INDEPENDENT (they are
     genuine contour integrals of an eps-flat quantity), and their ratio must stay exactly 2.
     That is the internal control: if THOSE drift with eps, the eps-flatness claim is wrong and
     something worse is going on than conditioning.
If instead the m=0 residual sits at 1.9e-2 for every eps, it is a real defect and the prescription
fails in the shifted geometry -- which, given the i control, would mean the two elliptic points
genuinely differ.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, E4, jay, ktil, RHO,                      # noqa: E402
                    rvec, inW, NU, check_orientation)

mp.mp.dps = 30
TAU0 = mp.e**(I * pi / 3)


def clustered(n, depth):
    xs = [mp.mpf(j) / n for j in range(n + 1)]
    d = mp.mpf(1) / n
    for _ in range(depth):
        d /= 2
        xs += [d, 1 - d]
    return sorted(set(xs))


def integ(path, g, n=20, depth=13):
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


def rhat(g, k, eps, turns, rE):
    S, C = spiral(eps), connector(eps, turns)
    n = k - 2
    rf = []
    for l in range(n + 1):
        s = mp.mpf(l + 1)
        LS = integ(S, lambda t: g(t) * t**(s - 1))
        LT = integ(C, lambda t: g(t) * ktil(t, s, k)) + integ(S, lambda t: g(t) * ktil(t, s, k))
        rf.append((2 * pi * I)**(n + 1) * (-1)**l * mp.binomial(n, l) * (LS + LT))
    Ph = integ(C, g) + integ(S, g)
    return [a - Ph * b for a, b in zip(rf, rE)], Ph


CASES = [
    ("rho P=4  E_4^2/j^2  k=8 ", lambda t: E4(t)**2 / jay(t)**2, 8, lambda t: E4(t)**2),
    ("rho P=2  E_4/j      k=4 ", lambda t: E4(t) / jay(t), 4, E4),
]

print("guard: Phi(E_12) plain arc = %s" % mp.nstr(check_orientation(), 12), flush=True)

for lbl, f, k, ref in CASES:
    print("=" * 78, flush=True)
    print(lbl, flush=True)
    rE = rvec(ref, k)
    for eps in (mp.mpf('1e-1'), mp.mpf('1e-2'), mp.mpf('1e-3')):
        d = {}
        for turns in (0, 1):
            v, Ph = rhat(f, k, eps, turns, rE)
            d[turns] = (v, Ph)
        v0 = [(2 * a + b) / 3 for a, b in zip(d[0][0], d[1][0])]
        n1 = max(abs(x) for x in NU(d[0][0], k))
        n2 = max(abs(x) for x in NU(d[1][0], k))
        print("   eps=%-6s |D(-1)| = %-16s |D(2)| = %-16s ratio = %-12s "
              "m=0: |(1+U+U^2)| = %-12s |defect| = %s"
              % (mp.nstr(eps, 3), mp.nstr(n1, 12), mp.nstr(n2, 12), mp.nstr(n2 / n1, 10),
                 mp.nstr(inW(v0, k)[1], 6), mp.nstr(max(abs(x) for x in NU(v0, k)), 6)),
              flush=True)
