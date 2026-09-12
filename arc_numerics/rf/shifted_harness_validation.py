"""Is the shifted-base-point harness sound?  Validate it where the answer is already known.

basepoint_W_rho.py reported |(1+U+U^2)| = O(1) for BOTH connector routings at rho and I read that
as mathematics.  It may be the construction: that run was the FIRST use of the spiral+connector
paths to build r_f, and nothing had checked that the harness reproduces a known answer.  The
control in basepoint_shift_rho.py validated L* for a pole-free form; it never validated hat r_f.

DAM's test.  Take f in F_k whose poles are nowhere near rho and nowhere on gamma^arc, where
lem:arcoff already guarantees hat r_f in W, and push them through BOTH harnesses:
  PLAIN    tau_0 = rho+1, both segments = gamma^arc          (common.tilde_r, already trusted)
  SHIFTED  tau_0 = (1+eps) e^{i pi/3}, Gamma_S = spiral,
           Gamma_T = connector + spiral                       (the machinery under suspicion)
For these f the two classes are homotopic -- no pole lies between them -- so they must agree.

Cases (all k=12, f = Delta/(j-x), so a SOLITARY simple pole orbit, and c_f(0)=0):
  7i           x = j(7i), huge real.  Far above the contour; the cleanest possible case.
  danger zone  tau = 0.2+0.99i: |tau| > 1 and Im tau < 1, i.e. just ABOVE the arc and below the
               height of i -- the horocyclic region that is dangerous for tau_0 = i, harmless here.
  below arc    tau = e^{i pi/4}: theta = 45 deg, BELOW the arc (Im = 0.707 < sqrt3/2), so x =
               j(e^{i pi/4}) is real and NEGATIVE and the pole orbit sits on the vertical edges of
               F.  (DAM proposed this as an on-arc point; it is not one -- gamma^arc is
               theta in [60,120] deg.  Kept because it is still off the contour and so still a
               lem:arcoff case.)

PREDICTIONS.
 (1) PLAIN: |(1+S)| and |(1+U+U^2)| both ~1e-20 or smaller for all three.  This is lem:arcoff and
     is really a check on common.py, which has been exercised before.
 (2) SHIFTED, short AND long: the SAME, ~1e-20.  With no pole between the two classes the
     connector routing cannot matter, and short - long must vanish identically (as it did for the
     pole-free control in basepoint_shift_rho.py).
 (3) tilde_r agrees between harnesses to quadrature precision.
If (2) fails, the spiral+connector harness is broken and the rho verdict is VOID.
If (2) holds, the rho defects are real and the base-point shift genuinely breaks W-membership.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, Delta, jay, ktil, E12, RHO,              # noqa: E402
                    rvec, tilde_r, inW, check_orientation)

mp.mp.dps = 30
TAU0 = mp.e**(I * pi / 3)
EPS = mp.mpf('1e-3')
K = 12


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


def shifted_tilde_r(g, k, eps, longway, rE):
    S, C = spiral(eps), connector(eps, longway)
    n = k - 2
    rf = []
    for l in range(n + 1):
        s = mp.mpf(l + 1)
        LS = integrate(S, lambda t: g(t) * t**(s - 1))
        LT = integrate(C, lambda t: g(t) * ktil(t, s, k)) + \
            integrate(S, lambda t: g(t) * ktil(t, s, k))
        rf.append((2 * pi * I)**(n + 1) * (-1)**l * mp.binomial(n, l) * (LS + LT))
    Ph = integrate(C, g) + integrate(S, g)
    return [a - Ph * b for a, b in zip(rf, rE)], Ph


print("guard: Phi(E_12) plain arc = %s" % mp.nstr(check_orientation(), 12), flush=True)
rE = rvec(E12, K)

CASES = [
    ("7i            pole far above the contour", jay(7 * I)),
    ("danger zone   tau = 0.2+0.99i (|tau|>1, Im<1)", jay(mp.mpf('0.2') + mp.mpf('0.99') * I)),
    ("below arc     tau = e^{i pi/4} (Im=0.707)", jay(mp.e**(I * pi / 4))),
]

for lbl, x in CASES:
    print("=" * 76, flush=True)
    print("%s      x = j = %s" % (lbl, mp.nstr(x, 10)), flush=True)
    f = lambda t, x=x: Delta(t) / (jay(t) - x)
    tp, Php = tilde_r(f, K)
    sp, up = inW(tp, K)
    print("   PLAIN    Phi = %-22s |(1+S)| = %-12s |(1+U+U^2)| = %s"
          % (mp.nstr(Php, 8), mp.nstr(sp, 6), mp.nstr(up, 6)), flush=True)
    store = {}
    for lw in (False, True):
        tr, Ph = shifted_tilde_r(f, K, EPS, lw, rE)
        s_, u_ = inW(tr, K)
        store["long" if lw else "short"] = tr
        print("   SHIFTED %-5s Phi = %-22s |(1+S)| = %-12s |(1+U+U^2)| = %s"
              % ("long" if lw else "short", mp.nstr(Ph, 8), mp.nstr(s_, 6),
                 mp.nstr(u_, 6)), flush=True)
    sc = max(abs(v) for v in tp)
    d_sl = max(abs(a - b) for a, b in zip(store['short'], store['long'])) / sc
    d_hp = max(abs(a - b) for a, b in zip(store['short'], tp)) / sc
    print("   short-long / scale = %-14s   shifted-plain / scale = %s"
          % (mp.nstr(d_sl, 6), mp.nstr(d_hp, 6)), flush=True)
