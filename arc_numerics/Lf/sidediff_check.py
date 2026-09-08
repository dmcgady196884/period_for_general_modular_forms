"""eq:sidediff of prop:arcside, at f-level, with the sign.

f = Delta/(j - x), k = 12, x in (0,1728), so j is 2-to-1 onto [0,1728] along the arc and f
has exactly two simple poles ON gamma^arc, at z and Sz with arg z + arg Sz = pi.

CLAIM.  Indenting outward at z rather than inward changes L^*(f,s) by
    -2 pi i e^{-i pi s/2} [ r_S(f,z) + r_T(f,z) ],
    r_S = Res_{tau=z} f(tau) tau^{s-1},   r_T = Res_{tau=z} f(tau) ktil(tau,s).
The sign comes from the arc running rho -> rho+1, i.e. theta DECREASING, so the outward
radial direction e^{i theta} = i (-i e^{i theta}) lies to the LEFT of travel; the sliver
between the two indentations is then to the RIGHT of the outward one and the loop winds
NEGATIVELY about z.

The second pole Sz is indented on a FIXED side in both runs, so it cancels in the
difference and only z is probed.

Residues exactly, not by finite difference: dj/dtau = 2 pi i (-E_4^2 E_6 / Delta), so
    Res_{tau=p} f g = Delta(p) g(p) / j'(p) = -Delta(p)^2 g(p) / (2 pi i E_4(p)^2 E_6(p)).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, E6, Delta, jay, ktil                  # noqa: E402

mp.mp.dps = 30
KK = 12
X0 = mp.mpf(500)
AMP = mp.mpf('0.15')
WID = mp.mpf('0.15')
SVALS = [mp.mpf(3), mp.mpf(7), mp.mpc('2.4', '0.6')]

ZP = mp.findroot(lambda t: jay(t) - X0, mp.e**(I * mp.mpf('1.85')))
ZM = -1 / ZP
THP, THM = mp.arg(ZP), mp.arg(ZM)


def f(t):
    return Delta(t) / (jay(t) - X0)


def res(p, g):
    """Res_{tau=p} f(tau) g(tau), exact via dj/dtau = 2 pi i (-E_4^2 E_6/Delta)"""
    return -Delta(p)**2 * g(p) / (2 * I * pi * E4(p)**2 * E6(p))


def bump(th, c):
    u = (th - c) / WID
    return mp.mpf(0) if abs(u) >= 1 else mp.e * mp.e**(-1 / (1 - u**2))


def dbump(th, c):
    u = (th - c) / WID
    if abs(u) >= 1:
        return mp.mpf(0)
    return bump(th, c) * (-2 * u / (1 - u**2)**2) / WID


def polar_int(g, sz, panels=64):
    """sz = sign of the indentation at z; the pole at Sz is always indented outward."""
    def lr(th):
        return AMP * (sz * bump(th, THP) + bump(th, THM))

    def dlr(th):
        return AMP * (sz * dbump(th, THP) + dbump(th, THM))

    def integrand(th):
        r = mp.e**(lr(th))
        return g(r * mp.e**(I * th)) * r * (dlr(th) + I) * mp.e**(I * th)
    a, b = pi / 3, 2 * pi / 3
    ns = sorted(set([a, b, THP, THM, THP - WID, THP + WID, THM - WID, THM + WID]
                    + [a + (b - a) * m / panels for m in range(panels + 1)]))
    ns = [x for x in ns if a <= x <= b]
    return -sum(mp.quad(integrand, [u, v]) for u, v in zip(ns[:-1], ns[1:]))


def Lstar(s, sz):
    return mp.e**(-I * pi * s / 2) * (
        polar_int(lambda t: f(t) * t**(s - 1), sz)
        + polar_int(lambda t: f(t) * ktil(t, s, KK), sz))


print("k = %d, x = %s, dps = %d" % (KK, mp.nstr(X0, 6), mp.mp.dps), flush=True)
print("poles on the arc: arg z = %s, arg Sz = %s, sum = %s"
      % (mp.nstr(THP, 8), mp.nstr(THM, 8), mp.nstr(THP + THM, 8)), flush=True)
print("%-14s %-26s %-26s %s" % ("s", "|L*_+ - L*_-|", "|prediction|", "rel. error"),
      flush=True)

for s in SVALS:
    Lp = Lstar(s, +1)
    Lm = Lstar(s, -1)
    rS = res(ZP, lambda t: t**(s - 1))
    rT = res(ZP, lambda t: ktil(t, s, KK))
    pred = -2 * I * pi * mp.e**(-I * pi * s / 2) * (rS + rT)
    d = Lp - Lm
    print("%-14s %-26s %-26s %s"
          % (mp.nstr(s, 8), mp.nstr(abs(d), 14), mp.nstr(abs(pred), 14),
             mp.nstr(abs(d - pred) / abs(pred), 6)), flush=True)
