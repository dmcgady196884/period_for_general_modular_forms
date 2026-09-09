"""END-TO-END: the assembled L*(f,s) on the arc against direct quadrature.

Every ingredient of sec:arcL has been checked in isolation.  This checks the ASSEMBLY --
lem:geotail + lem:geowhole + lem:geopolylog together, with the eps_z lens bookkeeping --
which is where a sign or a missed crossing would hide.

    L*(f,s) = c_Rf(0) G^geo(0) + sum_{n != 0} c_Rf(n) G^geo(n)                 [lem:geotail]
              + sum_p sum_{m=1}^{P_p} r*_{f,tau_p}(m) (-2 pi i)^m/(m-1)!
                                       I^geo_{m-1}(s,tau_p)                    [lem:geowhole]

    I^geo_N(s,z) = -delta_{N,0} G^geo(0) + (-1)^{N+1} sum_{n>=1} n^N e^{2 pi i n z} G^geo(-n)
                   + eps_z X_N(s,z)                                            [lem:geopolylog]

with eps_z = 1 exactly when z lies in the lens between the arc and the chord.

GROUND TRUTH.  Both segments are the same arc here, so
    L*(f,s) = e^{-i pi s/2} int_{gamma^arc} f(tau) [ tau^{s-1} + ktil(tau,s) ] dtau,
by quadrature, which never sees any of the above.

TEST FORMS.  f = Delta / prod_i (j - j(z_i))^{m_i} at k = 12, with the z_i placed by hand:
  A  one pole ABOVE the arc, outside the lens      -> mode sum only, eps_z = 0
  B  one pole IN the lens                          -> mode sum + X_N, eps_z = 1
  C  both at once                                  -> two projected poles, both branches
  D  a DOUBLE pole above the arc                   -> exercises the m-sum of lem:geowhole
Each z_i is chosen in the standard fundamental domain, so the rest of its SL2(Z)-orbit lies
below sqrt3/2 and survives into c_Rf; that is what makes the tail non-trivial.

c_Rf(n) is read off Rf = f - sum_p Phat_{tau_p} f by the trapezoid rule along the horizontal
line Im tau = sqrt3/2, which is spectrally accurate since Rf is periodic and holomorphic
there (lem:projspace_arc puts its poles strictly below).  Delta/(j-x)^m vanishes to order m
at the cusp and the blocks carry only positive modes, so c_Rf(n) = 0 for n <= 0.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, Delta, jay, ktil, arcint                   # noqa: E402

mp.mp.dps = 30
KK = 12
SVALS = [mp.mpf(3), mp.mpf(7), mp.mpc('2.4', '0.6')]
NMODE = 70           # modes of c_Rf kept
NSAMP = 256          # trapezoid samples (>> 2*NMODE, so no aliasing)
NPOLY = 220          # terms of the block mode sums
SQ32 = mp.sqrt(3) / 2
RHO = mp.e**(2 * I * pi / 3)
TAU0 = mp.e**(I * pi / 3)

CASES = [
    ("A  one pole above the arc      ", [(mp.mpc('0.30', '1.40'), 1)]),
    ("B  one pole in the lens        ", [(mp.mpc('0.20', '0.93'), 1)]),
    ("C  one above + one in the lens ", [(mp.mpc('0.30', '1.40'), 1),
                                         (mp.mpc('0.20', '0.93'), 1)]),
    ("D  a double pole above the arc ", [(mp.mpc('0.30', '1.40'), 2)]),
]


def in_lens(z):
    return abs(z) < 1 and z.imag > SQ32 and abs(z.real) <= mp.mpf(1) / 2


def Tred(z):
    """translate into |Re| <= 1/2"""
    return z - mp.floor(z.real + mp.mpf(1) / 2)


def to_project(z):
    """Every T-inequivalent pole of f in the SL2(Z)-orbit of z with Im >= sqrt3/2.

    def:geoproj projects all of them, not just the one placed.  A pole inside the unit
    disk is never alone: |z| < 1 sends Im(Sz) = Im z/|z|^2 ABOVE Im z, so a lens pole
    always drags its S-partner above the threshold too.  Deeper orbit words drop well
    below sqrt3/2 (checked), so {z, Sz} is the whole list.
    """
    out = []
    for w in (Tred(z), Tred(-1 / z)):
        if w.imag >= SQ32 and not any(abs(w - u) < mp.mpf('1e-18') for u in out):
            out.append(w)
    for w in (Tred(-1 / (z + 1)), Tred(-1 / (z - 1))):
        if w.imag >= SQ32:
            raise RuntimeError("deeper orbit member above threshold: %s" % mp.nstr(w, 8))
    return out


def make_f(poles):
    xs = [(jay(z), m) for z, m in poles]
    def f(t):
        d = mp.mpf(1)
        for x, m in xs:
            d *= (jay(t) - x)**m
        return Delta(t) / d
    return f


def block(tau, z, N):
    w = mp.e**(2 * I * pi * (tau - z))
    return {0: w / (1 - w), 1: w / (1 - w)**2, 2: w * (1 + w) / (1 - w)**3}[N]


def proj(tau, z, weights):
    """Phat_z f = sum_m r*(m) (-2 pi i)^m/(m-1)! Li_{1-m}(e^{2 pi i (tau - z)})"""
    return sum(r * (-2 * I * pi)**m / mp.factorial(m - 1) * block(tau, z, m - 1)
               for m, r in weights)


def pole_weights(f, z, P, rad=mp.mpf('0.02')):
    """r*(m) = (1/2 pi i) oint f(tau) (tau - z)^{m-1} dtau,  eq:laurent_arc"""
    out = []
    for m in range(1, P + 1):
        g = lambda th: (f(z + rad * mp.e**(I * th)) * (rad * mp.e**(I * th))**(m - 1)
                        * I * rad * mp.e**(I * th))
        v = sum(mp.quad(g, [2 * pi * a / 24, 2 * pi * (a + 1) / 24]) for a in range(24))
        out.append((m, v / (2 * I * pi)))
    return out


def Ggeo(s, n):
    """eq:geoG at any n != 0; (2 pi n)^w for n<0 uses arg = +pi, as one exponential."""
    xr, xt = -2 * I * pi * n * RHO, -2 * I * pi * n * TAU0
    lg = mp.log(2 * pi * abs(n)) + (I * pi if n < 0 else 0)
    pre = mp.e**(-s * (mp.log(2 * pi * abs(n)) - I * pi / 2 * mp.sign(n)))
    if n > 0:
        Sn = pre * (mp.gammainc(s, xr, mp.inf) - mp.gammainc(s, xt, mp.inf))
    else:
        e2 = mp.e**(2 * I * pi * s)
        pre = mp.e**(-s * (mp.log(2 * pi * abs(n)) + I * pi / 2))
        Sn = pre * (e2 * mp.gammainc(s, xr, mp.inf) + (1 - e2) * mp.gamma(s)
                    - mp.gammainc(s, xt, mp.inf))
    Tn = (mp.gammainc(s, xt, mp.inf) * mp.e**(-s * lg)
          + I**KK * mp.gammainc(KK - s, xt, mp.inf) * mp.e**(-(KK - s) * lg))
    return mp.e**(-I * pi * s / 2) * Sn + Tn


def G0(s):
    return (mp.e**(-I * pi * s / 2) * (TAU0**s - RHO**s) / s
            - (TAU0 / I)**s / s - I**KK * (TAU0 / I)**(KK - s) / (KK - s))


def prod(terms):
    p = mp.mpf(1)
    for t in terms:
        p *= t
    return p


def Xcross(s, z, N):
    """eq:geocross"""
    pw = (-1)**N * z**(s - 1 - N) * prod(s - 1 - j for j in range(N))
    z1 = mp.zeta(1 - s + N, z + 1) * prod(1 - s + j for j in range(N))
    z2 = (mp.e**(I * pi * (s - 1)) * mp.zeta(1 - KK + s + N, z + 1)
          * prod(1 - KK + s + j for j in range(N)))
    return mp.e**(-I * pi * s / 2) / (2 * I * pi)**N * (pw + z1 - z2)


def Igeo(s, z, N, Gneg):
    acc = sum(mp.mpf(n)**N * mp.e**(2 * I * pi * n * z) * Gneg[n - 1]
              for n in range(1, NPOLY + 1))
    out = (-G0(s) if N == 0 else mp.mpc(0)) + (-1)**(N + 1) * acc
    if in_lens(z):
        out += Xcross(s, z, N)
    return out


def cRf_modes(f, plist):
    """trapezoid along Im tau = sqrt3/2; returns c_Rf(n) for n = 1..NMODE"""
    y = SQ32
    vals = []
    for r in range(NSAMP):
        t = mp.mpf(r) / NSAMP + I * y
        vals.append(f(t) - sum(proj(t, z, w) for z, w in plist))
    out = []
    for n in range(1, NMODE + 1):
        acc = sum(v * mp.e**(-2 * I * pi * n * mp.mpf(r) / NSAMP)
                  for r, v in enumerate(vals)) / NSAMP
        out.append(acc * mp.e**(2 * pi * n * y))
    return out


print("dps=%d  k=%d  NMODE=%d  NSAMP=%d  NPOLY=%d\n"
      % (mp.mp.dps, KK, NMODE, NSAMP, NPOLY), flush=True)

for lbl, poles in CASES:
    f = make_f(poles)
    plist = []
    for z0, P in poles:
        for z in to_project(z0):
            plist.append((z, pole_weights(f, z, P)))
    print("%s\n%s" % ("=" * 78, lbl), flush=True)
    for z, _ in plist:
        print("    projected z = %-22s |z| = %-14s lens: %s"
              % (mp.nstr(z, 8), mp.nstr(abs(z), 8), in_lens(z)), flush=True)
    cR = cRf_modes(f, plist)
    print("    |c_Rf(1)| = %-14s |c_Rf(%d)| = %s"
          % (mp.nstr(abs(cR[0]), 8), NMODE, mp.nstr(abs(cR[-1]), 8)), flush=True)
    for s in SVALS:
        Gpos = [Ggeo(s, n) for n in range(1, NMODE + 1)]
        Gneg = [Ggeo(s, -n) for n in range(1, NPOLY + 1)]
        tail = sum(cR[n - 1] * Gpos[n - 1] for n in range(1, NMODE + 1))
        polar = mp.mpc(0)
        for (z, weights) in plist:
            polar += sum(r * (-2 * I * pi)**m / mp.factorial(m - 1)
                         * Igeo(s, z, m - 1, Gneg) for m, r in weights)
        got = tail + polar
        want = mp.e**(-I * pi * s / 2) * (arcint(lambda t: f(t) * t**(s - 1))
                                          + arcint(lambda t: f(t) * ktil(t, s, KK)))
        print("      s=%-12s assembled %-26s quad %-26s rel %s"
              % (mp.nstr(s, 6), mp.nstr(got, 14), mp.nstr(want, 14),
                 mp.nstr(abs(got - want) / abs(want), 6)), flush=True)
