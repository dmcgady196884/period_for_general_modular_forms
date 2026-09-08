"""The Pi + R decomposition at the threshold line Im z = sqrt3/2.

With E_j := e^{2 pi i j rho} = e^{2 pi i j tau_0} (equal, by T-periodicity),

    G^geo_{s,k}(j) ~ E_j sum_{l>=0} c_l(s,k) / (2 pi j)^{l+1},

    c_l = (s-1)_(l) [ -i i^l e^{-i pi s/2} (tau_0^{s-1-l} - rho^{s-1-l})
                      + (-i tau_0)^{s-1-l} ]
          + i^k (k-s-1)_(l) (-i tau_0)^{k-s-1-l},

(s-1)_(l) the FALLING factorial.  The first bracket is the gamma^* / S-segment endpoint
pair, obtained by repeated integration by parts of int_rho^tau_0 e^{2 pi i j tau}
tau^{s-1} dtau; the last term is the pair of incomplete gammas.  At tau_0 = i the S
segment degenerates and c_l collapses to (s-1)_(l) + i^k (k-s-1)_(l).

Then, for Im z = sqrt3/2 and any L >= N,

    I^geo_N(s,z) = sum_{l=0}^{L} c_l/(2 pi)^{l+1} Li_{l+1-N}(e^{2 pi i (tau_0 - z)})
                   + sum_{j>=1} j^N e^{-2 pi i j z} [ G^geo(j) - E_j sum_{l<=L} c_l/(2 pi j)^{l+1} ],

the second sum absolutely convergent with terms O(j^{N-L-2}).  Taking L > N is only a
numerical convenience: the extra Li's converge absolutely anyway.

TESTS
  A  order-by-order: | G/E_j - sum_{l<=L} c_l/(2 pi j)^{l+1} | * (2 pi j)^{L+2} tends to a
     constant, for L = 0,1,2 -- i.e. the expansion is right term by term.
  B  Pi + R against arc quadrature, N = 0,1,2, at a threshold pole off rho.
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
Z = mp.mpf('0.2') + I * SQ32            # on the threshold, away from rho
JMAX = 3000


def falling(a, l):
    p = mp.mpf(1)
    for j in range(l):
        p *= (a - j)
    return p


def c_l(s, l):
    """coefficient of (2 pi j)^{-l-1} in G^geo(j)/E_j"""
    mit = -I * TAU0
    Spart = (-I * I**l * mp.e**(-I * pi * s / 2)
             * (TAU0**(s - 1 - l) - RHO**(s - 1 - l)))
    return (falling(s - 1, l) * (Spart + mit**(s - 1 - l))
            + I**KK * falling(KK - s - 1, l) * mit**(KK - s - 1 - l))


def Ggeo_pos(s, m):
    xr, xt = -2 * I * pi * m * RHO, -2 * I * pi * m * TAU0
    pre = mp.e**(-s * (mp.log(2 * pi * m) - I * pi / 2))
    Sn = pre * (mp.gammainc(s, xr, mp.inf) - mp.gammainc(s, xt, mp.inf))
    Tn = (mp.gammainc(s, xt, mp.inf) * mp.e**(-s * mp.log(2 * pi * m))
          + I**KK * mp.gammainc(KK - s, xt, mp.inf)
          * mp.e**(-(KK - s) * mp.log(2 * pi * m)))
    return mp.e**(-I * pi * s / 2) * Sn + Tn


def arc_int(g, depth=12):
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


def block(tau, N):
    w = mp.e**(2 * I * pi * (tau - Z))
    return {0: w / (1 - w), 1: w / (1 - w)**2, 2: w * (1 + w) / (1 - w)**3}[N]


def Igeo_quad(s, N):
    g = lambda t: block(t, N)
    return mp.e**(-I * pi * s / 2) * (arc_int(lambda t: g(t) * t**(s - 1))
                                      + arc_int(lambda t: g(t) * ktil(t, s, KK)))


for s in SVALS:
    print("\n%s\ns = %s,  k = %d,  z = 0.2 + i sqrt3/2\n%s"
          % ("=" * 76, mp.nstr(s, 8), KK, "=" * 76), flush=True)

    print("  TEST A: residual after subtracting l <= L, scaled by (2 pi j)^{L+2}",
          flush=True)
    for L in (0, 1, 2):
        cs = [c_l(s, l) for l in range(L + 1)]
        row = []
        for j in (200, 400, 800, 1600):
            Ej = mp.e**(2 * I * pi * j * RHO)
            resid = Ggeo_pos(s, j) / Ej - sum(
                cs[l] / (2 * pi * j)**(l + 1) for l in range(L + 1))
            row.append(abs(resid) * (2 * pi * j)**(L + 2))
        print("    L=%d   " % L + "  ".join(mp.nstr(v, 8) for v in row), flush=True)

    print("  TEST B: Pi + R against arc quadrature", flush=True)
    for N in (0, 1, 2):
        L = N + 3
        cs = [c_l(s, l) for l in range(L + 1)]
        w = mp.e**(2 * I * pi * (TAU0 - Z))
        Pi = sum(cs[l] / (2 * pi)**(l + 1) * mp.polylog(l + 1 - N, w)
                 for l in range(L + 1))
        R = mp.mpc(0)
        for j in range(1, JMAX + 1):
            Ej = mp.e**(2 * I * pi * j * RHO)
            tail = Ggeo_pos(s, j) - Ej * sum(
                cs[l] / (2 * pi * j)**(l + 1) for l in range(L + 1))
            R += mp.mpf(j)**N * mp.e**(-2 * I * pi * j * Z) * tail
        tgt = Igeo_quad(s, N)
        print("    N=%d  quad = %-34s  Pi+R = %-34s  rel = %s"
              % (N, mp.nstr(tgt, 12), mp.nstr(Pi + R, 12),
                 mp.nstr(abs(Pi + R - tgt) / abs(tgt), 6)), flush=True)
