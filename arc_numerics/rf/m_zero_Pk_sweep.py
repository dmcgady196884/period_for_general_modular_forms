"""Is D(m=0) = 0 general, or an artefact of the two (P,k) already checked?

That the (1+U+U^2) defect is affine in the quotient winding m is a THEOREM -- lem:wall_arc is
linear in the residues and hat r is assembled linearly from L* and Phi.  The content of the m=0
prescription is that the CONSTANT TERM vanishes, i.e. the winding-free routing is off by exactly
one unit.  So far verified at (i, P=2) directly (rf/i_pole_m_zero.py: sides fail with identical
raw defects 2744182.69, mean in W at 8.4e-31 / 2.4e-28) and at (rho, P=3) via the four-winding
law D = -A m, A = 501060.84529 (rf/winding_plus_one.py).  P mod n has mattered at every other
stage of this problem, so two points is not a rule.

Cases.  At rho, P = 3m - a; at i, P = 2m - a.
  rho P=2  E_4/j          k=4   (a = ord_rho E_4 = 1, m = 1)
  rho P=4  E_4^2/j^2      k=8   (a = 2, m = 2)          [P = 1 mod 3, the degenerate class]
  i   P=4  E_4^3/(j-1728)^2  k=12 (a = 0, m = 2)

Two geometries, because the two elliptic points sit differently on gamma^arc:
  rho -- pole at the PINNED ENDPOINTS, so the base point must be displaced and hat gamma^T needs
         a connector.  m = -1 + 3*turns, so turns 0 and 1 give m = -1 and 2, and
             hat r(m=0) = (2/3) hat r(-1) + (1/3) hat r(2).
  i   -- pole INTERIOR, so the UNSHIFTED base point works and the two detours are m = -+1, with
             hat r(m=0) = (1/2) hat r(-1) + (1/2) hat r(+1).

PREDICTIONS.
 (1) |(1+U+U^2)| at m=0 down at the harness floor (1e-17..1e-20 for the shifted geometry at rho,
     1e-28 for the unshifted detour geometry at i), against O(1) for every realisable winding.
 (2) At i, BOTH single sides fail and the failure shows in (1+S) too, since S carries an inside
     detour to an outside one.  At rho, (1+S) ~ 1e-26 for every routing (the spiral is
     S-symmetric by construction) and only U discriminates.
 (3) At i the two single sides have IDENTICAL raw defect magnitudes, since D = -A m at m = -+1.
 (4) rho P=4 is the P = 1 mod 3 class, where the OLD delta^kappa analysis degenerated
     (c_lead vanished to 1e-19/1e-20).  Nothing in the m=0 argument cares about P mod n, so a
     failure HERE and nowhere else would say the constant term's vanishing is P-dependent -- the
     single most informative way this could break.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, E4, Delta, jay, ktil, E12, RHO,          # noqa: E402
                    rvec, inW, NU, check_orientation)

mp.mp.dps = 30
TAU0 = mp.e**(I * pi / 3)
EPS = mp.mpf('1e-3')


# ---------------------------------------------------------------- shifted geometry (rho)
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


def rhat_rho(g, k, turns, rE):
    S, C = spiral(EPS), connector(EPS, turns)
    n = k - 2
    rf = []
    for l in range(n + 1):
        s = mp.mpf(l + 1)
        LS = integ(S, lambda t: g(t) * t**(s - 1))
        LT = integ(C, lambda t: g(t) * ktil(t, s, k)) + integ(S, lambda t: g(t) * ktil(t, s, k))
        rf.append((2 * pi * I)**(n + 1) * (-1)**l * mp.binomial(n, l) * (LS + LT))
    Ph = integ(C, g) + integ(S, g)
    return [a - Ph * b for a, b in zip(rf, rE)], Ph


# ---------------------------------------------------------------- unshifted geometry (i)
def arcseg(g, th0, th1, n=20):
    xs = [th0 + (th1 - th0) * mp.mpf(j) / n for j in range(n + 1)]
    return sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
               for u, v in zip(xs[:-1], xs[1:]))


def detour_int(g, r, side, n=24):
    al = 2 * mp.asin(r / 2)
    psi0 = pi + al / 2
    psi1 = -al / 2 if side == 'out' else 2 * pi - al / 2
    ps = [psi0 + (psi1 - psi0) * mp.mpf(j) / n for j in range(n + 1)]
    circ = sum(mp.quad(lambda p: g(I + r * mp.e**(I * p)) * I * r * mp.e**(I * p), [u, v])
               for u, v in zip(ps[:-1], ps[1:]))
    return arcseg(g, 2 * pi / 3, pi / 2 + al) + circ + arcseg(g, pi / 2 - al, pi / 3)


def rhat_i(g, k, r, side, rE):
    n = k - 2
    rf = []
    for l in range(n + 1):
        s = mp.mpf(l + 1)
        L = detour_int(lambda t: g(t) * t**(s - 1), r, side) + \
            detour_int(lambda t: g(t) * ktil(t, s, k), r, side)
        rf.append((2 * pi * I)**(n + 1) * (-1)**l * mp.binomial(n, l) * L)
    Ph = detour_int(g, r, side)
    return [a - Ph * b for a, b in zip(rf, rE)], Ph


def show(tag, v, k, Ph):
    s_, u_ = inW(v, k)
    print("   %-14s Phi = %-24s |(1+S)| = %-12s |(1+U+U^2)| = %-12s |defect| = %s"
          % (tag, mp.nstr(Ph, 10), mp.nstr(s_, 6), mp.nstr(u_, 6),
             mp.nstr(max(abs(x) for x in NU(v, k)), 10)), flush=True)


print("guard: Phi(E_12) plain arc = %s" % mp.nstr(check_orientation(), 12), flush=True)

E8 = lambda t: E4(t)**2

RHO_CASES = [
    ("rho P=2  E_4/j      k=4 ", lambda t: E4(t) / jay(t), 4, E4),
    ("rho P=4  E_4^2/j^2  k=8 ", lambda t: E4(t)**2 / jay(t)**2, 8, E8),
]

for lbl, f, k, ref in RHO_CASES:
    print("=" * 78, flush=True)
    print(lbl, flush=True)
    rE = rvec(ref, k)
    got = {}
    for turns in (0, 1):
        v, Ph = rhat_rho(f, k, turns, rE)
        got[turns] = (v, Ph)
        show("turns=%d (m=%+d)" % (turns, -1 + 3 * turns), v, k, Ph)
    v0 = [(2 * a + b) / 3 for a, b in zip(got[0][0], got[1][0])]
    Ph0 = (2 * got[0][1] + got[1][1]) / 3
    show("m=0 (2/3,1/3)", v0, k, Ph0)

print("=" * 78, flush=True)
print("i   P=4  E_4^3/(j-1728)^2  k=12", flush=True)
fi = lambda t: E4(t)**3 / (jay(t) - 1728)**2
rE = rvec(E12, 12)
got = {}
for side in ('in', 'out'):
    v, Ph = rhat_i(fi, 12, mp.mpf('0.1'), side, rE)
    got[side] = (v, Ph)
    show("%s (m=%+d)" % (side, -1 if side == 'in' else 1), v, 12, Ph)
v0 = [(a + b) / 2 for a, b in zip(got['in'][0], got['out'][0])]
show("m=0 (mean)", v0, 12, (got['in'][1] + got['out'][1]) / 2)
