"""Shared machinery for the geodesic reference class, tau_0 = rho + 1.

This is the ONLY place the arc contour is defined.  Every script in arc_numerics/
imports from here, so the two traps below cannot be re-introduced independently in
eight different files (which is how the second one survived for a whole session).

TRAP 1 -- ORIENTATION.  gamma^arc runs rho -> rho+1, i.e. theta DECREASING from
2pi/3 to pi/3.  Integrating the other way flips r_f, Phi and r_{E_k} all at once,
which turns  r_f - Phi r_{E_k}  into  -(r_f + Phi r_{E_k})  -- a wrong combination
that still looks plausible and shows up only as a spurious (1+U+U^2) defect of
order 1.  `check_orientation()` asserts Phi(E_12) = +1 and every script calls it.

TRAP 2 -- WHICH Phi.  def:Phi is the integral over the FULL arc.  The delta-class
chord integral (tau_0(delta) - 1 -> tau_0(delta)) is a DIFFERENT number.  Using the
chord makes Phi(f) = 0 for a form with a pole at rho and nothing lands in W.  Use
`arcint` for def:Phi; `chordint` exists only for the eq:geodelta computations that
genuinely need it.

TRAP 3 -- BRANCH OF A NEGATIVE REAL RAISED TO A COMPLEX POWER.  eq:geoG contains
(-2 pi n)^w.  Writing  base = mp.e**(mp.log(2*pi*n) + I*pi)  and then  base**w  DESTROYS
the branch: that base is a negative real only up to a roundoff-sized imaginary part, and
**w takes mpmath's principal branch of whatever it actually is, so the result is decided
by the SIGN OF THE ROUNDOFF and flips with dps and with n.  It is stable-looking at
integer w (both branches agree there) and silently erratic at complex w.  Always write
one explicit exponential instead:

    lg = mp.log(2*pi*n) + I*pi          # arg(-2 pi n) = +pi
    term = mp.gammainc(w, xt, mp.inf) * mp.e**(-w * lg)

The correct branch is arg(-2 pi n) = +pi, the one making
arg(-2 pi n) + arg(tau_0/i) = 5 pi/6 = arg(2 pi i n tau_0).  Diagnosed in
Lf/F_only_check.py; the S-part was never affected because it was already written as a
single exponential.

TRAP 4 -- CLUSTERING AT THE ENDPOINTS WHEN THE POLE IS NOT AT AN ENDPOINT.  arc_nodes()
clusters at theta = pi/3, 2pi/3 because that is where the ELLIPTIC poles (rho, rho+1) sit.
An on-arc pole at angle phi is at t = (2pi/3 - phi)/(pi/3), i.e. MID-interval -- t ~ 0.235
and 0.765 for x = 500 -- and an endpoint-clustered mesh steps straight over it.  At
eps = 1e-3 the shifted spiral passes 5e-4 from such a pole against a panel width 0.05: the
quadrature never sees it and returns a SMOOTH WRONG NUMBER that can still look like it is
in W (2e-20 was reported).  Resolved value for x = 500 with pole-centred clustering is
|hat r| = 212068.947897 against 312510.344 from the broken mesh -- 32% out, with no
convergence warning.  This has now bitten twice: once in the lem:arcdeform test (pole
2e-8 from a mesh of width 0.26, read as two prescriptions disagreeing) and once in
rf/shifted_onarc_75.py.  ALWAYS pass the pole angles as extra clustering centres; see
rf/onarc_resolved.py for the pattern.  Note the hand-indented wiggle contour
(log r = 0.15 sin(6(theta - pi/2))) is NOT mesh-sensitive here, since it stands 0.16 off
the poles, and it reproduces arc_eight_classes.py's |hat r| = 312512.435 to all digits.

Quadrature: nodes cluster geometrically at both endpoints, since that is where the
ELLIPTIC poles approach the contour.  depth=10 already agrees with depth=13 to 5e-41, so
depth 10 is ample; the cost is linear in depth.  For any pole elsewhere on the arc this
default is WRONG -- see TRAP 4.

L^*(f,s) here is RAW QUADRATURE of the two contour integrals of def:Lint.  There is
no closed form for L^* in this class yet (the thm:mero analogue is unwritten), so the
period-polynomial results in rf/ are logically independent of that machinery and
cannot inherit an error from it.  The only special function used is the kernel
ktil = zeta(1-s, tau+1) - e^{i pi (s-1)} zeta(1-(k-s), tau+1), inside the integrand.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "end_to_end_numerical_testing"))

from validate import mp, I, pi, E4, E6, Delta, jay, ktil, slash   # noqa: E402

# ---------------------------------------------------------------- group elements
Sm = (0, -1, 1, 0)
Um = (1, -1, 1, 0)          # U = TS,  U^3 = -I,  U fixes rho+1

RHO = mp.e**(2 * I * pi / 3)        # e^{2 i pi/3}, left endpoint of the arc
TAU0 = mp.e**(I * pi / 3)           # rho + 1, the base point, fixed by U

C691 = mp.mpf(432000) / 691         # E_4^3 = E_12 + C691 * Delta


def E12(t):
    return E4(t)**3 - C691 * Delta(t)


def E43(t):
    return E4(t)**3


# ---------------------------------------------------------------- the arc contour
def arc_nodes(depth):
    """theta nodes on [pi/3, 2pi/3], clustered geometrically toward both ends"""
    a, b = pi / 3, 2 * pi / 3
    mid = (a + b) / 2
    left, right, d = [], [], mid - a
    for _ in range(depth):
        d /= 2
        left.append(a + d)
        right.append(b - d)
    return sorted(set([a] + left + [mid] + right + [b]))


def arcint(g, depth=10):
    """int_{gamma^arc} g(tau) dtau, oriented rho -> rho+1 (theta decreasing).

    This is def:Phi when g = f.  See TRAP 1.
    """
    th = arc_nodes(depth)
    return -sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
                for u, v in zip(th[:-1], th[1:]))


def subarcint(g, delta, depth=10):
    """int over the eq:geodelta sub-arc, from e^{i(2pi/3 - delta)} to e^{i(pi/3 + delta)}"""
    a, b = 2 * pi / 3 - delta, pi / 3 + delta
    th = [a + (b - a) * m / (2 ** depth) for m in range(2 ** depth + 1)]
    return sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
               for u, v in zip(th[:-1], th[1:]))


def chordint(g, delta, nsub=16):
    """int over the straight chord tau_0(delta) - 1 -> tau_0(delta).

    Horizontal at height sin(pi/3 + delta), hence ABOVE the poles at rho, rho+1.
    NOT def:Phi.  See TRAP 2.
    """
    from validate import path_int
    t0 = mp.e**(I * (pi / 3 + delta))
    return path_int(g, [t0 - 1, t0], nsub=nsub)


def check_orientation(depth=10, tol=mp.mpf('1e-20')):
    """Phi(E_12) must be +1.  Raises if the arc is being traversed backwards."""
    v = arcint(E12, depth)
    if abs(v - 1) > tol:
        raise AssertionError(
            "arc orientation/normalisation wrong: Phi(E_12) = %s, expected +1. "
            "See TRAP 1 in arc_numerics/common.py." % mp.nstr(v, 12))
    return v


# ---------------------------------------------------------------- period polynomial
def rvec(g, k=12, depth=10):
    """r_f of def:rf, in the geodesic class, by quadrature.

    r_f = (2 pi i)^{n+1} sum_l (-1)^l C(n,l) L^*(f, l+1) X^{n-l} Y^l, and def:rf's
    i^{l+1} cancels def:Lint's e^{-i pi s/2} at s = l+1, so neither appears.
    Both segments are the full arc (they coincide in this class).
    """
    n = k - 2
    out = []
    for l in range(n + 1):
        s = mp.mpf(l + 1)
        L = (arcint(lambda t: g(t) * t**(s - 1), depth)
             + arcint(lambda t: g(t) * ktil(t, s, k), depth))
        out.append((2 * pi * I)**(n + 1) * (-1)**l * mp.binomial(n, l) * L)
    return out


def tilde_r(g, k=12, depth=10, ref=None):
    """thm:geoperiod's tilde r_f = r_f - Phi(f) r_{E_k}.  Returns (tilde_r, Phi).

    `ref` defaults to E_12 at k = 12.  Passing E43 instead is what
    rf/eisenstein_dependence.py varies -- the choice is NOT inert.
    """
    if ref is None:
        ref = E12
    Phi = arcint(g, depth)
    rg, rr = rvec(g, k, depth), rvec(ref, k, depth)
    return [a - Phi * b for a, b in zip(rg, rr)], Phi


# ---------------------------------------------------------------- V_n operators
def NS(v, k=12):
    return [a + b for a, b in zip(v, slash(v, k - 2, Sm))]


def NU(v, k=12):
    a = slash(v, k - 2, Um)
    b = slash(a, k - 2, Um)
    return [x + y + z for x, y, z in zip(v, a, b)]


def inW(v, k=12):
    """(relative |(1+S)| , relative |(1+U+U^2)|)"""
    sc = max(abs(x) for x in v)
    return (max(abs(x) for x in NS(v, k)) / sc,
            max(abs(x) for x in NU(v, k)) / sc)


def inner(P, Q, k=12):
    """eq:inner"""
    n = k - 2
    return sum((-1)**a * P[a] * Q[n - a] / mp.binomial(n, a) for a in range(n + 1))


def opmat(op, k=12):
    """matrix of a linear map on V_n, in the monomial basis"""
    n = k - 2
    M = mp.matrix(n + 1, n + 1)
    for l in range(n + 1):
        e = [mp.mpc(0)] * (n + 1)
        e[l] = mp.mpc(1)
        col = op(e)
        for r in range(n + 1):
            M[r, l] = col[r]
    return M


def nullspace(M, tol=mp.mpf('1e-18')):
    n = M.cols
    _, sv, Vh = mp.svd_c(M)
    top = max(abs(sv[i]) for i in range(len(sv)))
    out = []
    for i in range(n):
        s = abs(sv[i]) / top if i < len(sv) else mp.mpf(0)
        if i >= len(sv) or s <= tol:
            out.append([mp.conj(Vh[i, j]) for j in range(n)])
    return out


def rank(M, tol=mp.mpf('1e-18')):
    sv = mp.svd_c(M, compute_uv=False)
    top = max(abs(sv[i]) for i in range(len(sv)))
    return 0 if top == 0 else sum(1 for i in range(len(sv))
                                  if abs(sv[i]) / top > tol)


def dim_Sk(k):
    if k < 12 or k % 2:
        return 0
    return k // 12 - (1 if k % 12 == 2 else 0)


# ---------------------------------------------------------------- def:Wpm at k = 12
# Derived from def:Wpm (Bernoulli formula) and cross-checked against the classical
# period polynomials of Delta.  Supports are DISJOINT -- W_- on odd l, W_+ on
# l = 2,4,6,8 (def:Wpm's c p_0 subtraction clears l = 0,10), p_0 on l = 0,10 -- so
# a decomposition in this basis is a read-off, not a solve.
W12_minus = [0, 4, 0, -25, 0, 42, 0, -25, 0, 4, 0]
W12_plus = [0, 0, 1, 0, -3, 0, 3, 0, -1, 0, 0]
P0_12 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]


def decompose12(v):
    """v = a W_- + alpha W_+ + beta p_0 at k = 12, by read-off. Returns (a, alpha, beta)."""
    return v[5] / 42, v[2], v[0]


# ---------------------------------------------------------------- extrapolation
def neville(xs, ys):
    """Extrapolate y(x) to x = 0. Returns (value, |last-order shift|) as an error bar."""
    n = len(xs)
    T = [list(ys)]
    for kk in range(1, n):
        row = []
        for i in range(n - kk):
            row.append((T[kk - 1][i] * (0 - xs[i + kk]) - T[kk - 1][i + 1] * (0 - xs[i]))
                       / (xs[i] - xs[i + kk]))
        T.append(row)
    return T[n - 1][0], abs(T[n - 1][0] - T[n - 2][0])
