"""A block pole ON gamma^arc: is it a PV, a residue, or nothing?

DAM's question.  The claim tested here is that it is a SIDE CHOICE, and that the two
sides differ by exactly the crossing term X_N of eq:geocross:

    contour dips INSIDE the unit disk near z  ->  z sits above the contour, is not swept
                                                  when we push down to the chord
                                              ->  eq:geoIN bare, NO residue
    contour bulges OUTSIDE the unit disk      ->  z lies between contour and chord
                                              ->  eq:geoIN + X_N, FULL residue
    mean of the two                           ->  eq:geoIN + X_N/2

Note the mode sum is not in trouble on the arc: a pole at |z| = 1 with z != rho, rho+1
has sqrt3/2 < Im z <= 1, so the inversion still converges on the chord.  The only datum
lem:geopolylog is missing there is the side.

CAVEAT on the word "PV".  Each indented path is a finite number fixed by its homotopy
class, and the two classes differ by 2 pi i Res whatever the pole order.  That is NOT
symmetric excision: excising a disc of radius r about a pole of order N+1 and letting
r -> 0 diverges like r^{-N} for N >= 1.  The mean of the two side-choices is finite for
every N; it agrees with the Cauchy principal value only at N = 0.

CONTOURS.  Radial perturbation of the arc, as in rf/arc_eight_classes.py: in polar form
tau = e^{log r(theta) + i theta} with log r a bump of amplitude +-AMP centred at arg z.
Only r moves, so arg tau stays in [pi/3, 2pi/3] and tau^{s-1} keeps its branch.

CHECKS, per (s, N):
    P1  inside            vs  eq:geoIN
    P2  outside           vs  eq:geoIN + X_N
    P3  mean              vs  eq:geoIN + X_N/2
    P4  outside - inside  vs  X_N            <- formula-free, no truncation floor
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, ktil                                      # noqa: E402

mp.mp.dps = 60
NMAX = 200
KK = 12
CS = len(sys.argv) > 1 and sys.argv[1] == "cs"          # complex s only, for a fast re-run
SVALS = [mp.mpc('2.4', '0.6')] if CS else [mp.mpf(3), mp.mpc('2.4', '0.6')]
NLIST = [0, 1, 2]
AMP = mp.mpf('0.15')
WID = mp.mpf('0.15')

RHO = mp.e**(2 * I * pi / 3)
TAU0 = mp.e**(I * pi / 3)

# theta = pi/2 is the elliptic point i, fixed by S.  Included because an S-SYMMETRIC
# gamma^S cannot pass on a definite side of it: the side at Sz is forced opposite to the
# side at z, and Si = i.  At block level no such constraint is imposed, so both sides are
# legitimate here and the test below should behave exactly as at a generic arc pole; the
# obstruction is a statement about assembling f, not about lambda_{N,z}.
POLES = [("generic arc pole, theta = 1.4", mp.mpf('1.4')),
         ("elliptic point i, theta = pi/2", pi / 2)]


def block(tau, z, N):
    w = mp.e**(2 * I * pi * (tau - z))
    return {0: w / (1 - w), 1: w / (1 - w)**2, 2: w * (1 + w) / (1 - w)**3}[N]


# ------------------------------------------------------- bump-perturbed arc contour
def bump(th, c):
    u = (th - c) / WID
    return mp.mpf(0) if abs(u) >= 1 else mp.e * mp.e**(-1 / (1 - u**2))


def dbump(th, c):
    u = (th - c) / WID
    if abs(u) >= 1:
        return mp.mpf(0)
    return bump(th, c) * (-2 * u / (1 - u**2)**2) / WID


def polar_int(g, sign, c, panels=48):
    """int over tau = exp(log r + i theta), log r = sign * AMP * bump, rho -> rho+1.
    sign = +1 bulges OUTSIDE the unit circle, sign = -1 dips INSIDE."""
    def integrand(th):
        r = mp.e**(sign * AMP * bump(th, c))
        return g(r * mp.e**(I * th)) * r * (sign * AMP * dbump(th, c) + I) * mp.e**(I * th)
    a, b = pi / 3, 2 * pi / 3
    ns = sorted(set([a, b, c, c - WID, c + WID]
                    + [a + (b - a) * m / panels for m in range(panels + 1)]))
    return -sum(mp.quad(integrand, [u, v]) for u, v in zip(ns[:-1], ns[1:]))


def Lstar(z, N, s, sign, c):
    g = lambda t: block(t, z, N)
    return mp.e**(-I * pi * s / 2) * (polar_int(lambda t: g(t) * t**(s - 1), sign, c)
                                      + polar_int(lambda t: g(t) * ktil(t, s, KK), sign, c))


# ------------------------------------------------------------------ eq:geoIN, eq:geocross
def Gtab(s, nmax):
    e2, Gs = mp.e**(2 * I * pi * s), mp.gamma(s)
    out = []
    for n in range(1, nmax + 1):
        xr, xt = 2 * I * pi * n * RHO, 2 * I * pi * n * TAU0
        pre = mp.e**(-s * (mp.log(2 * pi * n) + I * pi / 2))
        Sn = pre * (e2 * mp.gammainc(s, xr, mp.inf) + (1 - e2) * Gs
                    - mp.gammainc(s, xt, mp.inf))
        # (-2 pi n)^{-w} as ONE explicit exponential, arg(-2 pi n) = +pi.  Never build the
        # negative real and then raise it to a complex power -- see Lf/F_only_check.py.
        lg = mp.log(2 * pi * n) + I * pi
        Tn = (mp.gammainc(s, xt, mp.inf) * mp.e**(-s * lg)
              + I**KK * mp.gammainc(KK - s, xt, mp.inf) * mp.e**(-(KK - s) * lg))
        out.append(mp.e**(-I * pi * s / 2) * Sn + Tn)
    return out


def G0(s):
    return (mp.e**(-I * pi * s / 2) * (TAU0**s - RHO**s) / s
            - (TAU0 / I)**s / s - I**KK * (TAU0 / I)**(KK - s) / (KK - s))


def geoIN(s, z, N, G):
    acc = sum(mp.mpf(n)**N * mp.e**(2 * I * pi * n * z) * G[n - 1]
              for n in range(1, len(G) + 1))
    return (-G0(s) if N == 0 else mp.mpc(0)) + (-1)**(N + 1) * acc


def prod(terms):
    p = mp.mpf(1)
    for t in terms:
        p *= t
    return p


def X(s, z, N):
    """eq:geocross"""
    pw = (-1)**N * z**(s - 1 - N) * prod(s - 1 - j for j in range(N))
    z1 = mp.zeta(1 - s + N, z + 1) * prod(1 - s + j for j in range(N))
    z2 = (mp.e**(I * pi * (s - 1)) * mp.zeta(1 - KK + s + N, z + 1)
          * prod(1 - KK + s + j for j in range(N)))
    return mp.e**(-I * pi * s / 2) / (2 * I * pi)**N * (pw + z1 - z2)


print("dps = %d  NMAX = %d  k = %d" % (mp.mp.dps, NMAX, KK), flush=True)
for lbl, th in POLES:
    z = mp.e**(I * th)
    print("  %s:  |z| = %s,  Im z = %s  (sqrt3/2 = %s)"
          % (lbl, mp.nstr(abs(z), 12), mp.nstr(z.imag, 10),
             mp.nstr(mp.sqrt(3) / 2, 10)), flush=True)

for s in SVALS:
    G = Gtab(s, NMAX)
    print("\n%s\ns = %s\n%s" % ("=" * 74, mp.nstr(s, 8), "=" * 74), flush=True)
    for lbl, th in POLES:
        z = mp.e**(I * th)
        for N in NLIST:
            ins = Lstar(z, N, s, -1, th)
            out = Lstar(z, N, s, +1, th)
            mean = (ins + out) / 2
            F = geoIN(s, z, N, G)
            Xn = X(s, z, N)
            sc = abs(out)
            print("  %-32s N=%d  P1 %-11s P2 %-11s P3 %-11s P4 %s"
                  % (lbl, N,
                     mp.nstr(abs(ins - F) / sc, 4),
                     mp.nstr(abs(out - (F + Xn)) / sc, 4),
                     mp.nstr(abs(mean - (F + Xn / 2)) / sc, 4),
                     mp.nstr(abs((out - ins) - Xn) / sc, 4)), flush=True)
