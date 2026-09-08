"""A pole at exactly Im z = sqrt3/2: does the DEFINING series reproduce the contour?

For z above the arc one inverts Li_{-N}.  At Im z = sqrt3/2 the pole lies BELOW the arc
(at Re z the arc sits at sqrt(1 - Re z^2) > sqrt3/2), so the defining series is the right
one:
    I^geo_N(s,z) = sum_{j>=1} j^N e^{-2 pi i j z} G^geo_{s,k}(j),
absolutely convergent for Im z < sqrt3/2 and exactly MARGINAL at sqrt3/2, because the arc
attains that height -- at its two endpoints only.

PREDICTION.  Endpoint integration by parts gives G^geo(j) ~ e^{2 pi i j rho} C / j, and
e^{2 pi i j tau_0} = e^{2 pi i j rho} by T-periodicity, so the j-th term behaves like
C j^{N-1} e^{2 pi i j (rho - z)} with rho - z REAL (both at height sqrt3/2).  Hence:
  N = 0, z != rho : partial sums converge like a Li_1, error ~ 1/J.
  N = 0, z  = rho : argument is 1, Li_1(1) diverges logarithmically -- partial sums grow
                    like log J and never settle.
Only N = 0 is tested by naive summation; for N >= 1 the terms are O(j^{N-1}) and do not
even tend to zero, which is precisely why the gamma unfolding (peeling ell = 0..N) is
needed rather than direct summation.

G^geo(m) at m > 0 needs no sheet factor: arg(-2 pi i m) = -pi/2 and arg(tau) in
[pi/3, 2pi/3], so arg(-2 pi i m tau) lies in [-pi/6, pi/6], principal throughout.
Powers are written as single explicit exponentials (see TRAP 3 in common.py).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, ktil                                      # noqa: E402

mp.mp.dps = 30
KK = 12
SVALS = [mp.mpf(3), mp.mpc('2.4', '0.6')]
RHO = mp.e**(2 * I * pi / 3)
TAU0 = mp.e**(I * pi / 3)
SQ32 = mp.sqrt(3) / 2

# CONTROL FIRST.  At Im z = 0.7 < sqrt3/2 the pole is strictly below the arc everywhere,
# the defining series converges ABSOLUTELY at rate e^{-2 pi (sqrt3/2 - 0.7)} = 0.35, and
# term-by-term integration is unimpeachable.  If the sum fails to reproduce the quadrature
# THERE, the bug is in Ggeo_pos and not in the threshold analysis.
POLES = [("z = 0.2 + 0.70i      (CONTROL)", mp.mpf('0.2') + I * mp.mpf('0.7')),
         ("z = 0.2 + i sqrt3/2  (off rho)", mp.mpf('0.2') + I * SQ32),
         ("z = rho              (at rho) ", RHO)]


def block0(tau, z):
    w = mp.e**(2 * I * pi * (tau - z))
    return w / (1 - w)


def arc_int(g, z, depth=12):
    """clustered at both endpoints, where the arc comes down to height sqrt3/2"""
    a, b = pi / 3, 2 * pi / 3
    ns = [a, b, (a + b) / 2]
    for base in (a, b):
        d = (b - a) / 4
        for _ in range(depth):
            d /= 2
            ns += [base - d, base + d]
    ns = sorted(set(x for x in ns if a <= x <= b))
    return -sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
                for u, v in zip(ns[:-1], ns[1:]))


def Igeo_quad(s, z):
    g = lambda t: block0(t, z)
    return mp.e**(-I * pi * s / 2) * (arc_int(lambda t: g(t) * t**(s - 1), z)
                                      + arc_int(lambda t: g(t) * ktil(t, s, KK), z))


def Ggeo_pos(s, m):
    """eq:geoG at m > 0.  All arguments principal; no eq:sheet factor needed."""
    xr, xt = -2 * I * pi * m * RHO, -2 * I * pi * m * TAU0
    # (-2 pi i m)^{-s}, arg(-2 pi i m) = -pi/2, as one exponential
    pre = mp.e**(-s * (mp.log(2 * pi * m) - I * pi / 2))
    Sn = pre * (mp.gammainc(s, xr, mp.inf) - mp.gammainc(s, xt, mp.inf))
    Tn = (mp.gammainc(s, xt, mp.inf) * mp.e**(-s * mp.log(2 * pi * m))
          + I**KK * mp.gammainc(KK - s, xt, mp.inf)
          * mp.e**(-(KK - s) * mp.log(2 * pi * m)))
    return mp.e**(-I * pi * s / 2) * Sn + Tn


for s in SVALS:
    print("\n%s\ns = %s,  k = %d,  dps = %d\n%s"
          % ("=" * 74, mp.nstr(s, 8), KK, mp.mp.dps, "=" * 74), flush=True)
    for lbl, z in POLES:
        x = TAU0 - z
        print("  %s   tau_0 - z = %s" % (lbl, mp.nstr(x, 10)), flush=True)
        target = Igeo_quad(s, z)
        print("    arc quadrature = %s" % mp.nstr(target, 14), flush=True)
        acc = mp.mpc(0)
        nxt = 250
        for j in range(1, 4001):
            acc += mp.e**(-2 * I * pi * j * z) * Ggeo_pos(s, j)
            if j == nxt:
                err = abs(acc - target) / abs(target)
                print("      J = %-6d partial = %-34s rel.err = %-12s J*err = %s"
                      % (j, mp.nstr(acc, 12), mp.nstr(err, 6),
                         mp.nstr(j * err, 6)), flush=True)
                nxt *= 2
