"""Is a pole in the lens an obstruction to lem:geopolylog, or just a residue term?

DAM's reading, which this script tests: there is no obstruction.  Deform the arc down
to the chord, pick up the residue you cross, integrate.  eq:geoIN then returns the
CHORD value and the arc value differs from it by one explicit term.

THE LENS.  Lambda = { |tau| < 1, Im tau > sqrt3/2, |Re tau| <= 1/2 }, bounded above by
gamma^arc and below by the horizontal chord rho -> rho+1.  Lambda lies strictly inside
the strip |Re tau| <= 1/2, so of the block's poles at z + Z at most ONE, namely z
itself, ever lies in Lambda.  (This is also why the T-translation z -> z+1 "rescue" of
the previous session was a no-op: it moves the named pole out of the strip but leaves a
translate behind in Lambda.  Same function, same poles.)

THE PRINCIPAL PART.  With u = tau - z and v = 2 pi i u + t,
    sum_N Li_{-N}(e^{2 pi i u}) t^N / N!  =  e^v/(1 - e^v)  =  -sum_{m>=0} B_m (-v)^m / (m! v),
and only the m = 0 term -1/v carries negative powers of u:
    -1/(2 pi i u + t) = -sum_{r>=0} (-1)^r t^r / (2 pi i u)^{r+1},
whose t^N/N! coefficient is N! (-1)^{N+1} / (2 pi i u)^{N+1}.  So Li_{-N}(e^{2 pi i u})
has principal part EXACTLY
    (-1)^{N+1} N! / (2 pi i u)^{N+1},
a single term with no lower-order singular part.  Hence for a kernel kappa holomorphic
at z,
    Res_{tau = z} [ lambda_{N,z} kappa ] = (-1)^{N+1} kappa^{(N)}(z) / (2 pi i)^{N+1}.
Checked below against a small-circle contour integral.

THE SIGN.  gamma^arc and the chord both run rho -> rho+1, the arc above.  The closed
loop (arc forward, chord backward) bounds Lambda CLOCKWISE, so arc - chord = -2 pi i Res
and therefore

    I^geo_N(s,z)  =  [ RHS of eq:geoIN ]
                     + (-1)^N (2 pi i)^{-N} e^{-i pi s/2}
                       [ d^N/dz^N (z^{s-1})  +  d^N/dz^N ktil(z,s) ]        (z in Lambda)

which is the crossing cost already stated in the prose above lem:geopolylog.

WHAT IS CHECKED, per (pole, s, N):
    A   chord quadrature   vs  eq:geoIN                 -- formula computes the CHORD
    B   arc - chord        vs  the crossing term        -- the residue accounts for it
    C   arc quadrature     vs  eq:geoIN + crossing      -- the extended lemma
    D   crossing term      vs  small-circle residue     -- the principal part above
CONTROL: z = 0.3 + 1.4i is outside Lambda, where the crossing term must be absent and
A, B, C must all hold with it set to zero.

Usage:  lens_residue_check.py [quick]      quick = dps 30, NMAX 60, s = 3 only
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, ktil                                      # noqa: E402

QUICK = len(sys.argv) > 1 and sys.argv[1] == "quick"
_CS = len(sys.argv) > 1 and sys.argv[1] == "cs"
# dps 60 for the complex-s re-run: A and C are limited by the NMAX truncation (1e-54 down
# to 1e-37), not by precision, and complex-s Hurwitz zeta in the quadrature makes dps 100
# cost ~20 min per row for no gain.  B and D still clear the dps-60 floor by 40 orders.
mp.mp.dps = 30 if QUICK else 60 if _CS else 100
NMAX = 60 if QUICK else 220
KK = 12
CS = len(sys.argv) > 1 and sys.argv[1] == "cs"          # complex s only, for a fast re-run
SVALS = ([mp.mpf(3)] if QUICK else
         [mp.mpc('2.4', '0.6')] if CS else [mp.mpf(3), mp.mpc('2.4', '0.6')])
NLIST = [0] if QUICK else [0, 1, 2]

RHO = mp.e**(2 * I * pi / 3)
TAU0 = mp.e**(I * pi / 3)
SQ32 = mp.sqrt(3) / 2

POLES = [("LENS    z = 0.95i     ", mp.mpc('0', '0.95'), True),
         ("LENS    z = 0.2+0.93i ", mp.mpc('0.2', '0.93'), True),
         ("CONTROL z = 0.3+1.4i  ", mp.mpc('0.3', '1.4'), False)]


# ------------------------------------------------------------------ the blocks
def block(tau, z, N):
    """lambda_{N,z} = Li_{-N}(e^{2 pi i (tau - z)}), exact rational form in w."""
    w = mp.e**(2 * I * pi * (tau - z))
    if N == 0:
        return w / (1 - w)
    if N == 1:
        return w / (1 - w)**2
    if N == 2:
        return w * (1 + w) / (1 - w)**3
    raise ValueError(N)


# ------------------------------------------------------------------ quadrature
def _cluster(a, b, c, depth):
    """nodes on [a,b] clustered geometrically at the endpoints and at c"""
    ns = [a, b, (a + b) / 2]
    for base in (a, b, c):
        d = (b - a) / 4
        for _ in range(depth):
            d /= 2
            ns += [base - d, base + d]
        ns.append(base)
    return sorted(set(x for x in ns if a <= x <= b))


def arc_int(g, z, depth=11):
    """int over gamma^arc, rho -> rho+1 (theta DECREASING); clustered at arg z"""
    a, b = pi / 3, 2 * pi / 3
    ns = _cluster(a, b, mp.arg(z), depth)
    return -sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
                for u, v in zip(ns[:-1], ns[1:]))


def chord_int(g, z, depth=11):
    """int over the horizontal chord rho -> rho+1 at height sqrt3/2, clustered at Re z"""
    A, B = RHO, TAU0
    ns = _cluster(mp.mpf(0), mp.mpf(1), (z.real - A.real) / (B - A).real, depth)
    return sum(mp.quad(lambda x: g(A + (B - A) * x) * (B - A), [u, v])
               for u, v in zip(ns[:-1], ns[1:]))


def Lstar(z, N, s, integrator):
    g = lambda t: block(t, z, N)
    return mp.e**(-I * pi * s / 2) * (integrator(lambda t: g(t) * t**(s - 1), z)
                                      + integrator(lambda t: g(t) * ktil(t, s, KK), z))


# ------------------------------------------------------------------ eq:geoIN
def Gtab(s, nmax):
    """G^geo_{s,k}(-n) for n = 1..nmax, via eq:geoSupper for the S-part (the gamma^*
    series of eq:geoG loses ~2.7n digits and is useless at these n)."""
    e2 = mp.e**(2 * I * pi * s)
    Gs = mp.gamma(s)
    out = []
    for n in range(1, nmax + 1):
        xr, xt = 2 * I * pi * n * RHO, 2 * I * pi * n * TAU0
        # (2 pi i n)^{-s} with arg(2 pi i n) = pi/2, matching arg(x_tau0) = 5 pi/6
        pre = mp.e**(-s * (mp.log(2 * pi * n) + I * pi / 2))
        Sn = pre * (e2 * mp.gammainc(s, xr, mp.inf) + (1 - e2) * Gs
                    - mp.gammainc(s, xt, mp.inf))
        # (-2 pi n)^{-w} as ONE explicit exponential with arg(-2 pi n) = +pi, the branch
        # for which arg(-2 pi n) + arg(tau_0/i) = 5 pi/6.  Do NOT write
        # base = mp.e**(log(2 pi n) + i pi) and then base**w: that negative real carries a
        # roundoff-sized imaginary part and **w takes mpmath's principal branch of it, so
        # the branch ends up decided by the sign of the roundoff and flips with dps and n.
        # See Lf/F_only_check.py.
        lg = mp.log(2 * pi * n) + I * pi
        Tn = (mp.gammainc(s, xt, mp.inf) * mp.e**(-s * lg)
              + I**KK * mp.gammainc(KK - s, xt, mp.inf) * mp.e**(-(KK - s) * lg))
        out.append(mp.e**(-I * pi * s / 2) * Sn + Tn)
    return out


def G0(s):
    return (mp.e**(-I * pi * s / 2) * (TAU0**s - RHO**s) / s
            - (TAU0 / I)**s / s - I**KK * (TAU0 / I)**(KK - s) / (KK - s))


def geoIN(s, z, N, G, nmax):
    tot = -G0(s) if N == 0 else mp.mpc(0)
    acc = mp.mpc(0)
    for n in range(1, nmax + 1):
        acc += mp.mpf(n)**N * mp.e**(2 * I * pi * n * z) * G[n - 1]
    return tot + (-1)**(N + 1) * acc


# ------------------------------------------------------------------ the crossing
def rising(a, m):
    p = mp.mpf(1)
    for j in range(m):
        p *= (a + j)
    return p


def dN_pow(z, s, N):
    """d^N/dz^N z^{s-1}"""
    c = mp.mpf(1)
    for j in range(N):
        c *= (s - 1 - j)
    return c * z**(s - 1 - N)


def dN_ktil(z, s, N):
    """d^N/dz^N ktil(z,s):  d_a^N zeta(w,a) = (-1)^N (w)_N zeta(w+N,a),  a = z+1"""
    return (-1)**N * (rising(1 - s, N) * mp.zeta(1 - s + N, z + 1)
                      - mp.e**(I * pi * (s - 1))
                      * rising(1 - KK + s, N) * mp.zeta(1 - KK + s + N, z + 1))


def crossing(z, s, N):
    """(-1)^N (2 pi i)^{-N} e^{-i pi s/2} [ d^N z^{s-1} + d^N ktil ]"""
    return ((-1)**N / (2 * I * pi)**N * mp.e**(-I * pi * s / 2)
            * (dN_pow(z, s, N) + dN_ktil(z, s, N)))


def circle_residue(z, s, N, r=mp.mpf('0.01'), m=48):
    """-2 pi i * (sum of residues at z), by direct small-circle quadrature.
    Sign: the arc-minus-chord loop bounds the lens CLOCKWISE."""
    def g(th):
        t = z + r * mp.e**(I * th)
        return (block(t, z, N) * (t**(s - 1) + ktil(t, s, KK))
                * I * r * mp.e**(I * th))
    val = sum(mp.quad(g, [2 * pi * a / m, 2 * pi * (a + 1) / m]) for a in range(m))
    return -mp.e**(-I * pi * s / 2) * val


# ------------------------------------------------------------------ run
print("dps = %d   NMAX = %d   k = %d" % (mp.mp.dps, NMAX, KK), flush=True)
print("sqrt3/2 = %s" % mp.nstr(SQ32, 12), flush=True)
for lbl, z, inlens in POLES:
    print("   %s |z| = %-14s Im z = %-10s in lens: %s"
          % (lbl, mp.nstr(abs(z), 10), mp.nstr(z.imag, 8), inlens), flush=True)

for s in SVALS:
    t0 = time.time()
    G = Gtab(s, NMAX)
    print("\n%s\ns = %s   (G table: %.1f s)\n%s"
          % ("=" * 78, mp.nstr(s, 8), time.time() - t0, "=" * 78), flush=True)
    for lbl, z, inlens in POLES:
        for N in NLIST:
            arc = Lstar(z, N, s, arc_int)
            cho = Lstar(z, N, s, chord_int)
            fml = geoIN(s, z, N, G, NMAX)
            X = crossing(z, s, N) if inlens else mp.mpc(0)
            circ = circle_residue(z, s, N) if inlens else mp.mpc(0)
            sc = abs(arc)
            A = abs(cho - fml) / sc
            B = abs((arc - cho) - X) / sc
            C = abs(arc - (fml + X)) / sc
            D = abs(X - circ) / (abs(X) if inlens else mp.mpf(1))
            print("  %s N=%d  |arc|=%-12s A %-10s B %-10s C %-10s D %s"
                  % (lbl, N, mp.nstr(sc, 8), mp.nstr(A, 4), mp.nstr(B, 4),
                     mp.nstr(C, 4), mp.nstr(D, 4)), flush=True)
