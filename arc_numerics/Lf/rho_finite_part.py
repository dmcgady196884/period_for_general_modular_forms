"""Hole 2: the finite part of I^geo_N at z = rho, where the mode sum diverges.

At z = rho the block lambda_{N,rho} has poles at BOTH arc endpoints (rho and tau_0 are the
same point of the quotient).  Its contour integral converges when K(rho) = K(tau_0), i.e.
when k = 2 mod 6 (eq:P1cancel) -- which P = 1 forces -- yet its mode sum diverges like
c_0 log J.  There is no contradiction: at z = rho the defining series does not converge
uniformly at the endpoints, so term-by-term integration is invalid.  The regularisation
therefore cannot be read off the series; it must be matched to something convergent.

lem:geoedge is that something.  It is verified on the whole threshold line Im z = sqrt3/2
EXCEPT z = rho, and reads

    I^geo_N(s,z) = sum_{l<=L} c_l/(2 pi)^{l+1} Li_{l+1-N}(e^{2 pi i (tau_0 - z)})
                   + R^geo_{N,L}(s,z).

Let z -> rho ALONG that line.  Then e^{2 pi i (tau_0 - z)} -> 1 and the l = N term carries
Li_1(1), which diverges.  If I^geo_N(s,z) itself stays finite, R^geo must diverge oppositely
and the finite part is the limit of the sum.  Printed below, separately:

    A = the Li-sum, B = R^geo, A+B = I^geo, and (A+B) extrapolated as z -> rho.

k = 8 (= 2 mod 6, as P = 1 requires) and s = 7; s = 3 is unusable at weight 8, where
L*(f,3) vanishes identically.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, ktil                                      # noqa: E402

mp.mp.dps = 30
KK = 8
NN = 0
LL = NN + 3
NMAX = 4000
RHO = mp.e**(2 * I * pi / 3)
TAU0 = mp.e**(I * pi / 3)
SQ32 = mp.sqrt(3) / 2


def falling(a, l):
    p = mp.mpf(1)
    for j in range(l):
        p *= (a - j)
    return p


def c_l(s, l):
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


def pieces(s, z, cs, G):
    w = mp.e**(2 * I * pi * (TAU0 - z))
    A = sum(cs[l] / (2 * pi)**(l + 1) * mp.polylog(l + 1 - NN, w)
            for l in range(LL + 1))
    B = mp.mpc(0)
    for n in range(1, NMAX + 1):
        tail = G[n - 1] - mp.e**(2 * I * pi * n * RHO) * sum(
            cs[l] / (2 * pi * n)**(l + 1) for l in range(LL + 1))
        B += mp.mpf(n)**NN * mp.e**(-2 * I * pi * n * z) * tail
    return A, B


for sv in (7, 11):
    s = mp.mpf(sv)
    cs = [c_l(s, l) for l in range(LL + 1)]
    G = [Ggeo_pos(s, n) for n in range(1, NMAX + 1)]
    print("\n%s\nk = %d, N = %d, s = %s   (c_0 = %s)\n%s"
          % ("=" * 76, KK, NN, mp.nstr(s, 6), mp.nstr(cs[0], 10), "=" * 76), flush=True)
    print("  eps        |Li-sum A|       |R^geo B|        A+B", flush=True)
    prev = None
    for es in ('1e-1', '1e-2', '1e-3', '1e-4', '1e-5'):
        z = RHO + mp.mpf(es)          # along the threshold line Im z = sqrt3/2
        A, B = pieces(s, z, cs, G)
        d = "" if prev is None else "   change %s" % mp.nstr(abs(A + B - prev), 6)
        print("  %-10s %-16s %-16s %s%s"
              % (es, mp.nstr(abs(A), 10), mp.nstr(abs(B), 10),
                 mp.nstr(A + B, 14), d), flush=True)
        prev = A + B
