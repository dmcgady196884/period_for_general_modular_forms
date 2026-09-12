"""Shift tau_0 off rho+1 and the divergence at a pole on SL2(Z).rho goes away.

endpoint_vs_interior.py showed the plain arc integral DIVERGES when f has a pole at rho: the
contour terminates at the pole and no detour helps (c_-1 = 1.0e-4 + 1.2e-4 i at k=4, P=2).
DAM's move: do at rho what Section 4.1 did at i -- displace the BASE POINT and take a limit.

WHY IT SHOULD WORK, and why the prediction is flatness rather than mere convergence.  Moving
tau_0 -> tau_0' changes the two segments by

    Delta L*_S = -int_{tau_0}^{tau_0'} f tau^{s-1} dtau + int_{S tau_0}^{S tau_0'} f tau^{s-1} dtau

and tau -> S tau with f|_k S = f sends the second integrand to (-1)^{s-1} f tau^{k-s-1} dtau.
On the T-side, T-periodicity of f and the shift identity
    ktil(tau) - ktil(tau-1) = -tau^{s-1} + e^{i pi (s-1)} tau^{k-s-1}
give Delta L*_T = int_{tau_0}^{tau_0'} f [tau^{s-1} - e^{i pi(s-1)} tau^{k-s-1}] dtau.  Adding,
the bracket is -tau^{s-1} + (-1)^{s-1} tau^{k-s-1} + tau^{s-1} - e^{i pi(s-1)} tau^{k-s-1} = 0,
since (-1)^{s-1} = e^{i pi (s-1)}.  So L* is EXACTLY tau_0-independent whenever the connecting
path misses the poles.  eps = 0 is the one degenerate base point, where both segments terminate
ON the pole.

GEOMETRY.  tau_0 = (1+eps) e^{i pi/3}, so S tau_0 = (1+eps)^{-1} e^{2 i pi/3} and
T^{-1} tau_0 = rho + eps e^{i pi/3}: both at distance ~eps from rho, in directions 2 pi/3 apart.
That is forced -- V := T^{-1} S fixes rho with V'(rho) = rho, and V S tau_0 = T^{-1} tau_0
exactly in PSL2 -- and it is the whole point of displacing the base point: at eps = 0 the two
endpoints collide and carry K_S + K_T together, while at eps != 0 they sit on different points
of the stabiliser orbit and pick up different powers of omega = e^{2 pi i/3}.

Paths, both reducing to gamma^arc (traversed rho -> rho+1, i.e. theta DECREASING) at eps = 0:
  Gamma_S  log-spiral theta: 2pi/3 -> pi/3, log r linear and ODD about pi/2, hence S-symmetric
           (S acts by (r,theta) -> (1/r, pi-theta)), which is what lem:arcon's (1+S) needs.
  Gamma_T  log-spiral connector about rho from T^{-1}tau_0 to S tau_0, then Gamma_S.  The
           connector's two routings (short 2pi/3, long 4pi/3) differ by one loop about rho and
           so by 2 pi i e^{-i pi s/2} Res_rho(f ktil) -- the reopened Figure-1(b) ambiguity.

GUARDS.  Phi(E_12) = +1 (orientation trap: getting theta's direction backwards flips r_f, Phi
and r_{E_k} together and reads as the theorem failing), and a pole-free control checked against
a direct plain-arc integral.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, E6, Delta, jay, ktil                  # noqa: E402

mp.mp.dps = 40
RHO = mp.e**(2 * I * pi / 3)
TAU0 = mp.e**(I * pi / 3)
EPS = [mp.mpf(10)**(-e) for e in (2, 3, 4, 5)]
E12 = lambda t: E4(t)**3 - mp.mpf(432000) / 691 * Delta(t)


def clustered(n, depth):
    """n+1 nodes on [0,1], bisecting towards BOTH ends depth times (poles sit there)."""
    xs = [mp.mpf(j) / n for j in range(n + 1)]
    d = mp.mpf(1) / n
    for _ in range(depth):
        d /= 2
        xs += [d, 1 - d]
    return sorted(set(xs))


def integrate(path, g, n=24, depth=14):
    ts = clustered(n, depth)
    tot = 0
    for u, v in zip(ts[:-1], ts[1:]):
        tot += mp.quad(lambda t: g(path(t)[0]) * path(t)[1], [u, v])
    return tot


def spiral(eps):
    """theta 2pi/3 -> pi/3, log r odd about pi/2.  Returns (tau, dtau/dt)."""
    L = mp.log(1 + eps)

    def p(t):
        th = 2 * pi / 3 + t * (pi / 3 - 2 * pi / 3)
        lr = L * (pi / 2 - th) / (pi / 6)
        tau = mp.e**(lr + I * th)
        dth = pi / 3 - 2 * pi / 3
        return tau, tau * (-L / (pi / 6) + I) * dth
    return p


def connector(eps, longway):
    """log-spiral about rho from T^{-1}tau_0 to S tau_0."""
    a = (TAU0 * (1 + eps) - 1) - RHO            # T^{-1} tau_0 - rho
    b = RHO / (1 + eps) - RHO                   # S tau_0 - rho
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


def Lstar(f, s, k, eps, longway):
    S = spiral(eps)
    C = connector(eps, longway)
    ls = integrate(S, lambda t: f(t) * t**(s - 1))
    lt = integrate(C, lambda t: f(t) * ktil(t, s, k)) + \
        integrate(S, lambda t: f(t) * ktil(t, s, k))
    return mp.e**(-I * pi * s / 2) * (ls + lt)


def Phi(f, eps, longway):
    S, C = spiral(eps), connector(eps, longway)
    return integrate(C, f) + integrate(S, f)


def plain_arc(f, s, k):
    xs = [pi / 3 + (pi / 3) * mp.mpf(j) / 32 for j in range(33)]
    g = lambda t: f(t) * (t**(s - 1) + ktil(t, s, k))
    v = sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, w])
            for u, w in zip(xs[:-1], xs[1:]))
    return mp.e**(-I * pi * s / 2) * (-v)


print("GUARD  Phi(E_12) should be +1", flush=True)
for eps in EPS[:3]:
    print("   eps=%-8s Phi = %s" % (mp.nstr(eps, 3), mp.nstr(Phi(E12, eps, False), 18)),
          flush=True)

CASES = [
    ("control  E_4^3 (no poles),  k=12", lambda t: E4(t)**3, 12, 0),
    ("rho P=2  E_4/j,             k=4 ", lambda t: E4(t) / jay(t), 4, 2),
    ("rho P=3  Delta/j,           k=12", lambda t: Delta(t) / jay(t), 12, 3),
    ("rho P=4  E_4^2/j^2,         k=8 ", lambda t: E4(t)**2 / jay(t)**2, 8, 4),
]

for lbl, f, k, P in CASES:
    print("=" * 78, flush=True)
    print(lbl, flush=True)
    for s in (mp.mpf(7), mp.mpf(11)):
        if P == 0:
            print("   s=%-4s plain arc  %s" % (mp.nstr(s, 3), mp.nstr(plain_arc(f, s, k), 20)),
                  flush=True)
        vals = {}
        for lw in (False, True):
            for eps in EPS:
                vals[(lw, eps)] = Lstar(f, s, k, eps, lw)
                print("   s=%-4s %-5s eps=%-8s L* = %s"
                      % (mp.nstr(s, 3), "long" if lw else "short", mp.nstr(eps, 3),
                         mp.nstr(vals[(lw, eps)], 20)), flush=True)
        for lw in (False, True):
            sp = max(abs(vals[(lw, a)] - vals[(lw, b)]) for a in EPS for b in EPS)
            print("      spread over eps (%s): %s"
                  % ("long" if lw else "short", mp.nstr(sp, 6)), flush=True)
        d = vals[(False, EPS[-1])] - vals[(True, EPS[-1])]
        print("      short - long = %s" % mp.nstr(d, 18), flush=True)
        print("      Phi = %s" % mp.nstr(Phi(f, EPS[-1], False), 18), flush=True)
