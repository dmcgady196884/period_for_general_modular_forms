#!/usr/bin/env python3
r"""validate.py -- end-to-end numerical validation of finite_contour_cocycles_short.tex

GROUND TRUTH IS ALWAYS RAW QUADRATURE of a defining contour integral.  Every closed form in
the note is checked against that, never against another closed form.  Each check carries the
\label it validates, so a failure names a line in the .tex rather than a formula in Python.

Companion prose: validation_notes.tex (same directory), section by section.

    python3 validate.py            # fast tier: coarse dps, one case per regime
    python3 validate.py --slow     # slow tier: high precision, full cross-product
    python3 validate.py --layer A  # restrict to one layer

Layers (see validation_notes.tex):
    A  kernel and unfolding          <- upstream of everything; done
    B  tau_0-independence            <- done
    C  thm:weakL                     <- done
    D  periods and quasi-periods     <- pending
    E  meromorphic / section 4       <- pending
    F  external anchors (1806, Hurwitz class numbers)   <- pending
    G  suite hygiene                 <- pending
"""
import sys
import mpmath as mp

# ----------------------------------------------------------------------------- config
SLOW = "--slow" in sys.argv
ONLY_LAYER = None
if "--layer" in sys.argv:
    ONLY_LAYER = sys.argv[sys.argv.index("--layer") + 1].upper()

mp.mp.dps = 30 if SLOW else 15
TOL = mp.mpf(10)**(-(mp.mp.dps - 6))
NT = 80 if SLOW else 40          # q-series truncation for E4/E6
NMAX = 40 if SLOW else 18        # mode-sum truncation

I, pi = mp.j, mp.pi
FL = dict(flush=True)


# ------------------------------------------------------------------- modular forms
def _sigma(a, m):
    return sum(d**a for d in range(1, m + 1) if m % d == 0)


_E4C = [240 * _sigma(3, m) for m in range(NT, 0, -1)] + [mp.mpf(1)]
_E6C = [-504 * _sigma(5, m) for m in range(NT, 0, -1)] + [mp.mpf(1)]


def q(t):
    return mp.e**(2 * pi * I * t)


def E4(t):
    return mp.polyval(_E4C, q(t))


def E6(t):
    return mp.polyval(_E6C, q(t))


def Delta(t):
    return (E4(t)**3 - E6(t)**2) / 1728


def invDelta(t):
    return 1 / Delta(t)


# --- exact integer q-expansions -------------------------------------------------
# A validation suite's ground truth must be EXACT.  Extracting c_f(n) by a Fourier
# integral is catastrophically ill-conditioned for forms with a pole at the cusp:
# for 1/Delta the integrand runs as e^{2 pi n y} while the answer is only e^{4 pi sqrt n},
# so c(18) at y=1.3 asks for a 1e23 number out of 1e67-sized values.  Use integer
# series arithmetic instead.
def _smul(a, b, N):
    out = [0] * (N + 1)
    for i, ai in enumerate(a):
        if ai == 0 or i > N:
            continue
        for j, bj in enumerate(b):
            if i + j > N:
                break
            out[i + j] += ai * bj
    return out


def _sinv(a, N):
    """1/a for a[0] == 1."""
    out = [0] * (N + 1)
    out[0] = 1
    for n in range(1, N + 1):
        out[n] = -sum(a[k] * out[n - k] for k in range(1, min(n, len(a) - 1) + 1))
    return out


def _eta_prod(N):
    r"""\prod_{n\ge1}(1-q^n) truncated at q^N."""
    out = [0] * (N + 1)
    out[0] = 1
    for n in range(1, N + 1):
        fac = [0] * (N + 1)
        fac[0] = 1
        if n <= N:
            fac[n] = -1
        out = _smul(out, fac, N)
    return out


_NSER = NT + 4
_ETA24 = _eta_prod(_NSER)
for _ in range(3):                      # ()^24 = (((x^2)^2)^2)^3
    _ETA24 = _smul(_ETA24, _ETA24, _NSER)
_ETA24 = _smul(_smul(_ETA24, _ETA24, _NSER), _ETA24, _NSER)
_INV_ETA24 = _sinv(_ETA24, _NSER)


def cf_Delta(n):
    r"""Ramanujan tau: \Delta = q\prod(1-q^n)^{24}, so c(n) = [q^{n-1}]\prod(1-q^n)^{24}."""
    return mp.mpf(_ETA24[n - 1]) if 1 <= n <= _NSER else mp.mpf(0)


def cf_invDelta(n):
    r"""1/\Delta = q^{-1}\prod(1-q^n)^{-24}, so c(n) = [q^{n+1}]\prod(1-q^n)^{-24}."""
    return mp.mpf(_INV_ETA24[n + 1]) if -1 <= n <= _NSER - 1 else mp.mpf(0)


def cf_E4(n):
    return mp.mpf(1) if n == 0 else (mp.mpf(240 * _sigma(3, n)) if n > 0 else mp.mpf(0))


def cf_E6(n):
    return mp.mpf(1) if n == 0 else (mp.mpf(-504 * _sigma(5, n)) if n > 0 else mp.mpf(0))


_E4S = [1] + [240 * _sigma(3, m) for m in range(1, _NSER + 1)]
_E4CUBE = _smul(_smul(_E4S, _E4S, _NSER), _E4S, _NSER)
_JS = _smul(_E4CUBE, _INV_ETA24, _NSER)          # j = E4^3/Delta = q^{-1} E4^3 prod(1-q^n)^{-24}


def cf_j(n):
    r"""j = E_4^3/\Delta = q^{-1} + 744 + 196884q + ... ; weight 0."""
    return mp.mpf(_JS[n + 1]) if -1 <= n <= _NSER - 1 else mp.mpf(0)


def jay(t):
    return E4(t)**3 / Delta(t)


def residue_at(g, s0, r=mp.mpf('0.3'), M=None):
    r"""Res_{s=s0} g(s) = (1/2 pi i)\oint g ds, by the trapezoid rule on |s-s0|=r.

    Spectrally accurate for g meromorphic with a single pole inside, and far more stable than
    extrapolating s*g(s) -> Res, which only converges linearly in s."""
    M = M or (24 if SLOW else 12)
    tot = mp.mpc(0)
    for jj in range(M):
        th = 2 * pi * jj / M
        z = mp.e**(I * th)
        tot += g(s0 + r * z) * z
    return tot * r / M


# ------------------------------------------------------------------------- kernels
def ktil(tau, s, k):
    r"""Hurwitz T-kernel, Definition~\ref{def:kernel}."""
    return mp.zeta(1 - s, tau + 1) - mp.e**(I * pi * (s - 1)) * mp.zeta(1 - (k - s), tau + 1)


def LT_raw(f, k, s, tau0):
    r"""L^*_T by quadrature: e^{-i pi s/2} \int_{tau_0-1}^{tau_0} f ktil dtau.  GROUND TRUTH."""
    g = lambda x: f(tau0 - 1 + x) * ktil(tau0 - 1 + x, s, k)
    return mp.e**(-I * pi * s / 2) * mp.quad(g, [0, mp.mpf('0.25'), mp.mpf('0.5'),
                                                 mp.mpf('0.75'), 1])


def _gamma_cont(a, g_principal):
    r"""Gamma(a,z) continued once anticlockwise across the negative-real cut.

    Gamma(a,z) = Gamma(a) - gamma(a,z) and gamma(a,z) = z^a (entire), so the only
    multivaluedness is z^a -> e^{2 pi i a} z^a, giving
        Gamma_cont = Gamma(a)(1 - e^{2 pi i a}) + e^{2 pi i a} Gamma_principal.
    Needed for n<0 once Re(tau_0)<0: there z = -2 pi i n tau_0 has Im(z) carrying the sign of
    Re(tau_0), so the cut is crossed exactly as Re(tau_0) changes sign.  The integral varies
    continuously across that; the principal branch jumps.  This is what prop:HGG's phrase
    'the identity continues to n<0' requires -- NOT the literal principal value."""
    e = mp.e**(2 * pi * I * a)
    return mp.gamma(a) * (1 - e) + e * g_principal


def LT_closed(cf, k, s, tau0, nmax=None, arg_sign=+1):
    r"""eq:T-segment-closed-Mk.  cf(n) supplies the q-coefficients.

    For n<0 the symbol (2 pi n)^s needs a branch.  The phase is applied EXPLICITLY, never
    folded into the base as (2 pi |n| e^{i pi})^s: numerically e^{i pi} = -1 + eps*i with the
    sign of eps a rounding artifact, so the base-folded form silently flips branch with dps.
    arg_sign=+1 is the convention pinned by quadrature (see validation_notes.tex, A3).
    """
    nmax = nmax or NMAX
    tot = -cf(0) * ((tau0 / I)**s / s + I**k * (tau0 / I)**(k - s) / (k - s))
    for n in range(-nmax, nmax + 1):
        if n == 0:
            continue
        c = cf(n)
        if c == 0:
            continue
        z = -2 * pi * I * n * tau0
        den_s, den_ks = (2 * pi * abs(n))**s, (2 * pi * abs(n))**(k - s)
        if n < 0:
            den_s *= mp.e**(I * pi * s * arg_sign)
            den_ks *= mp.e**(I * pi * (k - s) * arg_sign)
        gs, gks = mp.gammainc(s, z, mp.inf), mp.gammainc(k - s, z, mp.inf)
        if n < 0 and mp.re(tau0) < 0:
            gs, gks = _gamma_cont(s, gs), _gamma_cont(k - s, gks)
        tot += c * (gs / den_s + I**k * gks / den_ks)
    return tot


# ------------------------------------------------------------------------ harness
CHECKS = []


def check(label, layer, tier="fast", desc="", tol=None):
    def deco(fn):
        CHECKS.append(dict(label=label, layer=layer, tier=tier, desc=desc, fn=fn, tol=tol))
        return fn
    return deco


def rel(got, want):
    r"""True relative error against a nonzero target; absolute error against an exact zero.

    An earlier version floored the denominator at 1 unconditionally.  That silently turns
    every small-magnitude comparison into an absolute one: the lem:wall checks below have
    jumps of order 1e-9 before normalisation, and a sign error there produced |got-want| ~
    1e-9, which 'passed'.  Only want == 0 exactly (a target that IS zero, e.g. the k=0
    residue) may be compared absolutely."""
    d = abs(got - want)
    return d if want == 0 else d / abs(want)


# ============================================================ LAYER A: kernel/unfolding
@check("def:kernel", "A", desc=r"$\delta_T\tk$ shift identity")
def a1_kernel_shift():
    r"""\delta_T ktil(tau,s) = -(tau+1)^{s-1} + e^{i pi (s-1)} (tau+1)^{k-1-s}."""
    for k in (4, 12, 18):
        for tau in (I, mp.mpf('0.3') + mp.mpf('1.4') * I, mp.mpf('-0.4') + 2 * I):
            for s in (mp.mpc('2.3', '0.4'), mp.mpc('0.5', '7.1'), mp.mpf(3)):
                got = ktil(tau + 1, s, k) - ktil(tau, s, k)
                want = (-(tau + 1)**(s - 1)
                        + mp.e**(I * pi * (s - 1)) * (tau + 1)**(k - 1 - s))
                yield ("k=%d tau=%s s=%s" % (k, mp.nstr(tau, 3), mp.nstr(s, 3)), got, want)


@check("eq:hurwitz-antideriv", "A", desc=r"Hurwitz antiderivative and shift")
def a2_hurwitz():
    r"""d_a[zeta(sigma-1,a)/(1-sigma)] = zeta(sigma,a);  zeta(sigma-1,a+1)-zeta(sigma-1,a) = -a^{1-sigma};
    and the definite integral \int_{i-1}^{i} zeta(sigma,tau+1) dtau = i^{1-sigma}/(sigma-1)."""
    for sigma in (mp.mpc('1.7', '0.3'), mp.mpc('-2.2', '1.1'), mp.mpf('3.5')):
        for a in (I, mp.mpf('0.4') + mp.mpf('1.3') * I):
            d = mp.diff(lambda x: mp.zeta(sigma - 1, x) / (1 - sigma), a)
            yield ("deriv sigma=%s a=%s" % (mp.nstr(sigma, 3), mp.nstr(a, 3)),
                   d, mp.zeta(sigma, a))
            got = mp.zeta(sigma - 1, a + 1) - mp.zeta(sigma - 1, a)
            yield ("shift sigma=%s a=%s" % (mp.nstr(sigma, 3), mp.nstr(a, 3)),
                   got, -a**(1 - sigma))
        got = mp.quad(lambda t: mp.zeta(sigma, I - 1 + t + 1), [0, 1])
        yield ("integral sigma=%s" % mp.nstr(sigma, 3), got, I**(1 - sigma) / (sigma - 1))


@check("eq:T-segment-closed-Mk", "A", desc=r"general-$\tau_0$ unfolding, positive modes")
def a3a_unfold_pos():
    r"""Delta (k=12): only n>=1, so no branch ambiguity.  Closed form vs raw quadrature."""
    k = 12
    cf = cf_Delta
    taus = [I, mp.mpf('0.3') + mp.mpf('1.3') * I]
    if SLOW:
        taus += [mp.mpf('-0.45') + mp.mpf('1.05') * I, 2 * I]
    for tau0 in taus:
        for s in (mp.mpc('2.3', '0.4'), mp.mpc('6', '3')):
            yield ("tau0=%s s=%s" % (mp.nstr(tau0, 3), mp.nstr(s, 3)),
                   LT_closed(cf, k, s, tau0), LT_raw(Delta, k, s, tau0))


@check("eq:T-segment-closed-Mk", "A", desc=r"general-$\tau_0$ unfolding, $c_f(0)\neq0$")
def a3b_unfold_const():
    r"""E4 (k=4): c_f(0)=1, so this exercises the -c_f(0)((tau_0/i)^s/s + ...) term."""
    k = 4
    cf = cf_E4
    for tau0 in (I, mp.mpf('0.3') + mp.mpf('1.3') * I):
        for s in (mp.mpc('2.3', '0.4'),):
            yield ("tau0=%s s=%s" % (mp.nstr(tau0, 3), mp.nstr(s, 3)),
                   LT_closed(cf, k, s, tau0), LT_raw(E4, k, s, tau0))


@check("eq:T-segment-closed-Mk", "A", desc=r"negative modes, ${\rm Re}\,\tau_0\ge0$")
def a3c_unfold_neg():
    r"""1/Delta (k=-12) has n=-1, so this exercises the branch of (2 pi n)^s and of
    Gamma(s,-2 pi i n tau_0) -- the 'incomplete-gamma unfolding ambiguity'.  On-axis tau_0
    (the note's actual contour, Definition~\ref{def:null_homotopy}) puts -2 pi i n tau_0
    exactly on the negative real axis, i.e. squarely on the cut."""
    k = -12
    taus = [I, mp.mpf('1.35') * I, mp.mpf('0.3') + mp.mpf('1.3') * I]
    if SLOW:
        taus += [mp.mpf('0.45') + mp.mpf('1.05') * I, 2 * I]
    for tau0 in taus:
        for s in (mp.mpc('2.3', '0.4'),):
            yield ("tau0=%s s=%s" % (mp.nstr(tau0, 3), mp.nstr(s, 3)),
                   LT_closed(cf_invDelta, k, s, tau0), LT_raw(invDelta, k, s, tau0))


@check("eq:T-segment-closed-Mk", "A", desc=r"negative modes, ${\rm Re}\,\tau_0<0$ (continued $\Gamma$)")
def a3d_unfold_neg_left():
    r"""For n<0 and Re(tau_0)<0 the point z = -2 pi i n tau_0 sits below the negative-real cut
    of Gamma, having crossed it as Re(tau_0) changed sign.  The closed form then requires the
    CONTINUED Gamma of _gamma_cont, not the principal value; with that, eq:T-segment-closed-Mk
    holds here exactly, as prop:HGG's continuation argument says it should."""
    k = -12
    taus = [mp.mpf('-0.35') + mp.mpf('1.45') * I, mp.mpf('-0.2') + mp.mpf('1.6') * I]
    if SLOW:
        taus += [mp.mpf('-0.48') + mp.mpf('1.1') * I]
    for tau0 in taus:
        yield ("tau0=%s" % mp.nstr(tau0, 3),
               LT_closed(cf_invDelta, k, mp.mpc('2.3', '0.4'), tau0),
               LT_raw(invDelta, k, mp.mpc('2.3', '0.4'), tau0))


# -------------------------------------------------------- paths, residues, full L^*
def path_int(g, pts, nsub=None):
    r"""\int g dtau along the polyline through pts (a list of points of \HH).

    Each leg is split into nsub panels before quadrature.  The wall-crossing paths of
    lem:wall run within ~0.2 of a simple pole, where a single panel per leg reaches only
    ~1e-8 -- enough to confirm a residue's sign, not enough to test an identity."""
    nsub = nsub or (12 if SLOW else 8)
    tot = mp.mpc(0)
    for a, b in zip(pts[:-1], pts[1:]):
        for m in range(nsub):
            u0, u1 = mp.mpf(m) / nsub, mp.mpf(m + 1) / nsub
            p, r = a + u0 * (b - a), a + u1 * (b - a)
            tot += (r - p) * mp.quad(lambda u: g(p + u * (r - p)),
                                     [0, mp.mpf('0.5'), 1])
    return tot


def LS_raw(f, k, s, tau0, pts=None):
    r"""L^*_S = e^{-i pi s/2}\int_{\gamma^S} f(tau) tau^{s-1} dtau, from -1/tau_0 to tau_0.

    Default path is the straight line, which stays in \HH since \HH is convex; pts overrides
    it with an explicit polyline, which is how the homotopy classes of lem:wall are realised."""
    pts = pts or [-1 / tau0, tau0]
    return mp.e**(-I * pi * s / 2) * path_int(lambda t: f(t) * t**(s - 1), pts)


def LT_raw_path(f, k, s, tau0, pts=None):
    r"""L^*_T along an explicit polyline from tau_0-1 to tau_0."""
    pts = pts or [tau0 - 1, tau0]
    return mp.e**(-I * pi * s / 2) * path_int(lambda t: f(t) * ktil(t, s, k), pts)


def L_raw(f, k, s, tau0):
    r"""The full L-integral of def:Lint at base-point tau_0, both segments by quadrature."""
    return LS_raw(f, k, s, tau0) + LT_raw_path(f, k, s, tau0)


def res_circle(g, p, r=mp.mpf('0.2'), M=None):
    r"""Res_{tau=p} g = (1/2 pi i)\oint g dtau on |tau-p|=r, trapezoid rule.

    The trapezoid rule on a circle converges like (r/R)^M, R being the distance to the next
    singularity.  For the Gamma-orbit of 2i the nearest other image is 2i+1, so R = 1 and the
    original defaults (r=0.3, M=16) gave only 0.3^16 ~ 4e-9 -- which is what made the lem:wall
    checks miss tolerance while agreeing to 8 digits.  r=0.2, M=32 gives 0.2^32."""
    M = M or (64 if SLOW else 32)
    tot = mp.mpc(0)
    for jj in range(M):
        z = mp.e**(2 * pi * I * jj / M)
        tot += g(p + r * z) * z
    return tot * r / M


J2I = mp.mpf(287496)                       # j(2i) = 66^3
_C2I = 1 / res_circle(lambda t: Delta(t) / (jay(t) - J2I), 2 * I)


def f2i(t):
    r"""Delta/(j - j(2i)) in F_12: a simple pole at 2i and at its Gamma-orbit.

    Normalised so that Res_{tau=2i} f2i = 1.  Unnormalised the residue is Delta(2i)/j'(2i)
    ~ 2e-12, which makes every wall-crossing jump ~1e-9 and leaves the comparison hostage to
    the absolute noise floor rather than testing the identity."""
    return _C2I * Delta(t) / (jay(t) - J2I)


# ==================================================== LAYER B: tau_0-independence, walls
@check("prop:indep", "B", desc=r"$\tau_0$-independence on the imaginary axis")
def b1_indep_onaxis():
    r"""L^* = L^*_S + L^*_T is independent of tau_0 for f in M_k^!: both segments move, and
    only the sum is invariant.  Reference value taken at tau_0 = i."""
    s = mp.mpc('2.3', '0.4')
    taus = [mp.mpf('1.3') * I, 2 * I, mp.mpf('0.8') * I]
    for name, f, cf, k in WEAK_FORMS[:4]:
        ref = L_raw(f, k, s, I)
        for tau0 in taus:
            yield ("%s (k=%d) tau0=%s" % (name, k, mp.nstr(tau0, 3)), L_raw(f, k, s, tau0), ref)


@check("prop:indep", "B", desc=r"$\tau_0$-independence off-axis, both signs of ${\rm Re}\,\tau_0$")
def b2_indep_offaxis():
    r"""Same, for tau_0 off the imaginary axis.  Re(tau_0)<0 is included deliberately: the
    CLOSED FORM eq:T-segment-closed-Mk fails there (see A3, report tier), so if the RAW
    L^* is nonetheless tau_0-independent, that failure is confined to the closed form and
    is not a defect of the L-integral itself."""
    s = mp.mpc('2.3', '0.4')
    taus = [mp.mpf('0.3') + mp.mpf('1.3') * I, mp.mpf('-0.35') + mp.mpf('1.45') * I]
    if SLOW:
        taus += [mp.mpf('0.5') + mp.mpf('0.9') * I, mp.mpf('-0.2') + 2 * I]
    for name, f, cf, k in WEAK_FORMS[:4]:
        ref = L_raw(f, k, s, I)
        for tau0 in taus:
            yield ("%s (k=%d) tau0=%s" % (name, k, mp.nstr(tau0, 3)), L_raw(f, k, s, tau0), ref)


@check("prop:indepHomotopy", "B", desc=r"$\tau_0$-independence for $f\in F_k$ within a class")
def b3_indep_meromorphic():
    r"""f = Delta/(j - j(2i)) in F_12, simple pole at 2i and Gamma-images (heights 2, 0.5,
    0.4, ...).  Base-points near i keep both segments clear of every image, so all of these
    lie in one homotopy class and L^* must not move."""
    s = mp.mpc('2.3', '0.4')
    ref = L_raw(f2i, 12, s, mp.mpf('1.1') * I)
    taus = [mp.mpf('1.2') * I, mp.mpf('1.05') * I, mp.mpf('0.2') + mp.mpf('1.15') * I]
    for tau0 in taus:
        yield ("tau0=%s" % mp.nstr(tau0, 4), L_raw(f2i, 12, s, tau0), ref)


@check("lem:wall", "B", desc=r"$T$-segment crossing: jump $=2\pi ie^{-i\pi s/2}r_T$")
def b4_wall_T():
    r"""Same end-points, two homotopy classes.  tau_0 = 0.5+2.2i puts the straight T-segment
    ABOVE the pole at 2i; the detour dips to height 1.8, putting the pole above it instead.
    The enclosed rectangle [-0.5,0.5]x[1.8,2.2] contains 2i and no other image (the orbit
    heights are 2, 0.5, 0.4, ...), so eq:wall predicts a single residue.

    Orientation: straight-then-detour-reversed runs right along the top, down the right, left
    along the bottom, up the left -- CLOCKWISE about 2i.  Hence X_T = -1, not +1."""
    for s in ([mp.mpc('2.3', '0.4')] + ([mp.mpc('1.4', '2.0')] if SLOW else [])):
        tau0 = mp.mpf('0.5') + mp.mpf('2.2') * I
        lo = mp.mpf('1.8')
        straight = LT_raw_path(f2i, 12, s, tau0)
        detour = LT_raw_path(f2i, 12, s, tau0,
                             pts=[tau0 - 1, mp.mpf('-0.5') + lo * I, mp.mpf('0.5') + lo * I, tau0])
        rT = res_circle(lambda t: f2i(t) * ktil(t, s, 12), 2 * I)
        yield ("s=%s" % mp.nstr(s, 3), straight - detour,
               -2 * pi * I * mp.e**(-I * pi * s / 2) * rT)


@check("lem:wall", "B", desc=r"$S$-segment crossing: jump $=2\pi ie^{-i\pi s/2}r_S$")
def b5_wall_S():
    r"""Two S-paths from -1/tau_0 to tau_0 that agree around every image except 2i, which one
    passes to the left and the other to the right.  Their difference encloses 2i alone, so
    eq:wall predicts a single residue.  The loop (up the left, across the top, down the right,
    across the bottom) is again CLOCKWISE, so X_S = -1."""
    for s in ([mp.mpc('2.3', '0.4')] + ([mp.mpc('1.4', '2.0')] if SLOW else [])):
        tau0 = mp.mpf('2.5') * I
        a, b = -1 / tau0, tau0
        x = mp.mpf('0.15')
        left = [a, -x + mp.mpf('0.4') * I, -x + mp.mpf('2.5') * I, b]
        right = [a, -x + mp.mpf('0.4') * I, -x + mp.mpf('1.7') * I,
                 x + mp.mpf('1.7') * I, x + mp.mpf('2.5') * I, b]
        rS = res_circle(lambda t: f2i(t) * t**(s - 1), 2 * I)
        yield ("s=%s" % mp.nstr(s, 3),
               LS_raw(f2i, 12, s, tau0, pts=left) - LS_raw(f2i, 12, s, tau0, pts=right),
               -2 * pi * I * mp.e**(-I * pi * s / 2) * rS)


# ==================================================== LAYER C: thm:weakL, functional eq.
# For f in M_k^!, tau_0 = i makes the S-segment degenerate, so L^* = L^*_T alone and
# LT_raw(f,k,s,I) is the full L^*(f,s) -- the ground truth for this layer.
WEAK_FORMS = [
    ("Delta",    Delta,    cf_Delta,    12),   # cusp form,     c_f(0)=0,  i^k=+1
    ("E4",       E4,       cf_E4,        4),   # c_f(0)=1,                 i^k=+1
    ("E6",       E6,       cf_E6,        6),   # c_f(0)=1,                 i^k=-1
    ("1/Delta",  invDelta, cf_invDelta, -12),  # pole at the cusp, n<0,    i^k=+1
    ("j",        jay,      cf_j,         0),   # k=0: the two poles collide
]


@check("eq:weakL", "C", desc=r"$L^*$ for $f\in M_k^!$ vs raw quadrature")
def c1_weakL():
    r"""eq:weakL against LT_raw at tau_0=i.  Covers c_f(0)=0 and !=0, i^k=+1 and -1,
    negative modes, and the weight-zero case where the s=0 and s=k poles coincide."""
    ss = [mp.mpc('2.3', '0.4')] + ([mp.mpc('5', '2'), mp.mpf('1.5')] if SLOW else [])
    for name, f, cf, k in WEAK_FORMS:
        for s in ss:
            yield ("%s (k=%d) s=%s" % (name, k, mp.nstr(s, 3)),
                   LT_closed(cf, k, s, I), LT_raw(f, k, s, I))


@check("eq:weakL", "C", desc=r"$\mathop{\rm Res}_{s=0}L^*=-c_f(0)$")
def c2_residue():
    r"""Residue extracted from the RAW integral by a contour in s, so this tests the note's
    polar claim against quadrature rather than against its own closed form.
    For k=0 the s=0 and s=k poles coincide and cancel: -c(0)(1/s + i^0/(0-s)) = 0, so the
    residue must vanish even though c_j(0)=744 =/= 0."""
    for name, f, cf, k in WEAK_FORMS:
        r = mp.mpf('0.3') if k != 4 else mp.mpf('0.25')
        got = residue_at(lambda s: LT_raw(f, k, s, I), mp.mpf(0), r=r)
        want = mp.mpf(0) if k == 0 else -cf(0)
        yield ("%s (k=%d) c_f(0)=%s" % (name, k, mp.nstr(cf(0), 6)), got, want)


@check("eq:weakL", "C", desc=r"functional equation $L^*(f,s)=i^kL^*(f,k-s)$")
def c3_functional_eq():
    r"""Checked on RAW quadrature both sides: the cheapest global probe of every phase in the
    construction.  E6 (i^k=-1) is the case that catches a lost sign."""
    ss = [mp.mpc('2.3', '0.4'), mp.mpc('1.1', '3.7')]
    if SLOW:
        ss += [mp.mpf('2.5'), mp.mpc('-0.7', '1.2')]
    for name, f, cf, k in WEAK_FORMS:
        for s in ss:
            yield ("%s (k=%d) s=%s  i^k=%s" % (name, k, mp.nstr(s, 3), mp.nstr(I**k, 3)),
                   LT_raw(f, k, s, I), I**k * LT_raw(f, k, k - s, I))


# ==================================================== LAYER E: the meromorphic machinery
def G_sk(s, k, m):
    r"""G_{s,k}(m) = Gamma(s,2 pi m)/(2 pi m)^s + i^k Gamma(k-s,2 pi m)/(2 pi m)^{k-s}.

    m<0 puts the argument on the negative real axis; the branch is the one pinned in A3,
    (2 pi m)^a := (2 pi|m|)^a e^{i pi a}, with Gamma principal (Re tau_0 = 0 here, so no
    continuation is required -- cf. _gamma_cont)."""
    z = 2 * pi * m
    out = mp.mpc(0)
    for a, pref in ((s, mp.mpf(1)), (k - s, I**k)):
        den = (2 * pi * abs(m))**a * (mp.e**(I * pi * a) if m < 0 else 1)
        out += pref * mp.gammainc(a, z, mp.inf) / den
    return out


def IN_raw(N, s, k, z):
    r"""GROUND TRUTH for I_N(s,z): the L-integral of lambda_{N,z} in the minimal class.

    By prop:indepHomotopy the delta->0 representative may be used; there the block is regular
    at i (since Im z > 1), so the S-integral vanishes and the T-segment is [i-1, i]."""
    g = lambda x: mp.polylog(-N, mp.e**(2 * pi * I * (I - 1 + x - z))) * ktil(I - 1 + x, s, k)
    return mp.e**(-I * pi * s / 2) * mp.quad(g, [0, mp.mpf('0.25'), mp.mpf('0.5'),
                                                 mp.mpf('0.75'), 1])


def IN_closed(N, s, k, z, nmax=None):
    r"""lem:polylogcf: delta_{N,0}(1/s + i^k/(k-s)) + (-1)^{N+1} sum_n n^N e^{2 pi i n z}G(-n)."""
    nmax = nmax or NMAX
    tot = (1 / s + I**k / (k - s)) if N == 0 else mp.mpc(0)
    acc = mp.mpc(0)
    for n in range(1, nmax + 1):
        acc += mp.mpf(n)**N * mp.e**(2 * pi * I * n * z) * G_sk(s, k, -n)
    return tot + (-1)**(N + 1) * acc


@check("eq:gammarec", "E", desc=r"incomplete-$\Gamma$ recursion")
def e1_gammarec():
    r"""Gamma(s,X) = e^{-X} sum_{l<J} (s-1)_{(l)} X^{s-1-l} + (s-1)_{(J)} Gamma(s-J,X)."""
    def fall(a, m):
        out = mp.mpf(1)
        for r in range(m):
            out *= (a - r)
        return out
    for s in (mp.mpc('2.3', '0.4'), mp.mpc('-1.2', '2.0')):
        for X in (mp.mpf('6.28'), mp.mpc('3.0', '2.0'), mp.mpf('20')):
            for J in (1, 3, 5):
                want = mp.gammainc(s, X, mp.inf)
                got = (mp.e**(-X) * sum(fall(s - 1, l) * X**(s - 1 - l) for l in range(J))
                       + fall(s - 1, J) * mp.gammainc(s - J, X, mp.inf))
                yield ("s=%s X=%s J=%d" % (mp.nstr(s, 3), mp.nstr(X, 3), J), got, want)


@check("eq:blockexp", "E", desc=r"block $=$ pure pole $+$ entire $T$-periodic tail")
def e2_blockexp():
    r"""(-2 pi i)^m/(m-1)! Li_{1-m}(e^{2 pi i(tau-i)}) - (tau-i)^{-m} must be REGULAR at tau=i:
    its residue there vanishes, which is what makes r*_{f,z}(m) the Laurent coefficients."""
    for m in (1, 2, 3):
        g = lambda t: ((-2 * pi * I)**m / mp.factorial(m - 1)
                       * mp.polylog(1 - m, mp.e**(2 * pi * I * (t - I))) - (t - I)**(-m))
        yield ("m=%d residue of remainder at i" % m, res_circle(g, I, r=mp.mpf('0.2')), mp.mpf(0))


@check("lem:polylogcf", "E", desc=r"$I_N(s,z)$ for ${\rm Im}\,z>1$ vs raw quadrature")
def e3_polylogcf():
    r"""THE load-bearing claim of section 4: everything downstream is built on it.
    N=0 is the case carrying the delta_{N,0}(1/s + i^k/(k-s)) inversion constant."""
    k = 12
    zs = [mp.mpf('0.3') + mp.mpf('1.4') * I, 2 * I]
    if SLOW:
        zs += [mp.mpf('-0.25') + mp.mpf('1.8') * I, mp.mpf('0.45') + mp.mpf('2.3') * I]
    for z in zs:
        for N in (0, 1, 2):
            for s in ([mp.mpc('2.3', '0.4')] + ([mp.mpc('5', '2')] if SLOW else [])):
                yield ("z=%s N=%d s=%s" % (mp.nstr(z, 4), N, mp.nstr(s, 3)),
                       IN_closed(N, s, k, z), IN_raw(N, s, k, z))


def _poch(a, m):
    out = mp.mpf(1)
    for r in range(m):
        out *= (a - r)
    return out


def FP(N, s, k, J=None, nmax=None):
    r"""def:FP at the minimal cutoff J=N+1 of rem:cutoff.  Four pieces: the harmonic term, the
    zeta-sum, the S-segment polynomial, and the reduced-Gamma remainder E^+_N."""
    # rem:cutoff proves FP_N is J-independent for J >= N+1, so any legal J may be used to
    # EVALUATE it.  J=N+1 is the minimal one but leaves a remainder decaying only as n^{-2},
    # i.e. an error ~1/nmax: at nmax=30 that is 1.5%.  J=N+5 decays as n^{-6} and reaches
    # 1e-12 by nmax=30, so it is the sane default for computing the value.
    J = N + 5 if J is None else J
    nmax = nmax or (60 if SLOW else 30)
    ik = I**k
    c = lambda l: _poch(s - 1, l) + ik * _poch(k - s - 1, l)
    HN = sum(mp.mpf(1) / r for r in range(1, N + 1))
    tot = c(N) * HN / (2 * pi)**(N + 1)
    for l in range(J):
        if l != N:
            tot += c(l) * mp.zeta(l + 1 - N) / (2 * pi)**(l + 1)
    tot += ((-1)**N * mp.factorial(N) / (2 * pi)**(N + 1)
            * sum((-1)**t * _poch(s - 1, t) / (mp.factorial(t) * (N - t)) for t in range(N)))
    for n in range(1, nmax + 1):
        x = 2 * pi * n
        tot += n**N * mp.e**x * (_poch(s - 1, J) * mp.gammainc(s - J, x, mp.inf) / x**s
                                 + ik * _poch(k - s - 1, J) * mp.gammainc(k - s - J, x, mp.inf)
                                 / x**(k - s))
    return tot


def laurent_i(f, P, rc=mp.mpf('0.05'), Np=96):
    r"""Laurent coefficients a_{-m} of f at tau=i, m=1..P, by circle sums."""
    a = {}
    for m in range(1, P + 1):
        a[m] = sum(f(I + rc * mp.e**(I * 2 * pi * j / Np)) * (rc * mp.e**(I * 2 * pi * j / Np))**m
                   for j in range(Np)) / Np
    return a


MERO_AT_I = [
    ("Delta/E6",        lambda t: Delta(t) / E6(t),                6, 1),
    ("E4^2 Delta/E6^2", lambda t: E4(t)**2 * Delta(t) / E6(t)**2,  8, 2),
    ("Delta^2/E6^3",    lambda t: Delta(t)**2 / E6(t)**3,          6, 3),
]


def wallvalue_raw(f, k, P, s, d):
    r"""GROUND TRUTH for lem:wallvalue: L^*(\widehat P_i f, s; i(1+delta)) by quadrature.

    T-segment against ktil, plus the S-segment split into a regular part and the principal
    value of the pure pole (def:null_homotopy takes the pole at i by PV)."""
    a = laurent_i(f, P)
    r = {m: (-2 * pi * I)**m * a[m] / mp.factorial(m - 1) for m in range(1, P + 1)}

    def Pf(t):
        e = mp.e**(2 * pi * I * (t - I))
        return sum(r[m] * mp.polylog(1 - m, e) for m in range(1, P + 1))

    def sing(t):
        return sum(a[m] * (t - I)**(-m) for m in range(1, P + 1))

    def Pf_reg(t):
        u = 2 * pi * I * (t - I)
        if abs(u) < mp.mpf('0.05'):
            return sum(r[m] * sum(mp.zeta(-(m - 1) - p) / mp.factorial(p) * u**p for p in range(8))
                       for m in range(1, P + 1))
        return Pf(t) - sing(t)

    t0 = I * (1 + d)
    T = mp.e**(-I * pi * s / 2) * mp.quad(lambda x: Pf(t0 - 1 + x) * ktil(t0 - 1 + x, s, k),
                                          [0, mp.mpf('0.02'), mp.mpf('0.1'), mp.mpf('0.5'),
                                           mp.mpf('0.9'), mp.mpf('0.98'), 1])
    ylo, yhi = 1 / (1 + d), 1 + d
    S = mp.quad(lambda y: y**(s - 1) * Pf_reg(I * y),
                [ylo, mp.mpf('0.999'), 1, mp.mpf('1.001'), yhi])
    for m in range(1, P + 1):
        coef = a[m] * I**(-m)

        def freg(y, m=m):
            h = y - 1
            if abs(h) < mp.mpf('1e-4'):
                return (_poch(s - 1, m) / mp.factorial(m)
                        + _poch(s - 1, m + 1) / mp.factorial(m + 1) * h)
            return (y**(s - 1) - sum(_poch(s - 1, t) / mp.factorial(t) * h**t
                                     for t in range(m))) / h**m
        val = mp.quad(freg, [ylo, mp.mpf('0.999'), 1, mp.mpf('1.001'), yhi])
        for t in range(m):
            bt = _poch(s - 1, t) / mp.factorial(t)
            p = t - m
            val += (bt * mp.log(1 + d) if p == -1
                    else bt * ((yhi - 1)**(p + 1) - (ylo - 1)**(p + 1)) / (p + 1))
        S += coef * val
    return T + S


@check("lem:wallvalue", "E", tol=mp.mpf('1e-5'),
       desc=r"pole at $i$: $\delta\to0$ limit vs $\sum r^*\mathrm{FP}$")
def e4_wallvalue():
    r"""The fanciest object in the note.  eq:wallvalue asserts a LIMIT, so the raw contour
    integral must be extrapolated, not sampled: L = (8f(d/4) - 6f(d/2) + f(d))/3 annihilates
    both an O(d) and an O(d^2) error term.  Both orders occur, split by the parity dichotomy
    of prop:ellord: 4|k with P even converges as d^2, k = 2 mod 4 with P odd only as d.

    Tolerance is 1e-5, not the global 1e-9: what remains after extrapolation is the O(d^3)
    term the three-point formula cannot annihilate.  The substantive evidence is the
    convergence RATE -- ratios of 2.0 and 4.0 per halving of d, straight to pred."""
    s = mp.mpc('2.3', '0.4')
    d = mp.mpf('0.04') if not SLOW else mp.mpf('0.02')
    for name, f, k, P in MERO_AT_I:
        a = laurent_i(f, P)
        pred = sum((-2 * pi * I)**m * a[m] / mp.factorial(m - 1) * FP(m - 1, s, k)
                   for m in range(1, P + 1))
        f1 = wallvalue_raw(f, k, P, s, d)
        f2 = wallvalue_raw(f, k, P, s, d / 2)
        f4 = wallvalue_raw(f, k, P, s, d / 4)
        yield ("%s (k=%d P=%d) extrapolated" % (name, k, P), (8 * f4 - 6 * f2 + f1) / 3, pred)


@check("rem:cutoff", "E", desc=r"$\mathrm{FP}_N$ independent of $J$ for $J\ge N+1$")
def e5_cutoff():
    r"""rem:cutoff: any larger J moves terms between the zeta-sum and E^+ and must leave FP_N
    unchanged.  Compared among the cutoffs that converge geometrically (J >= N+2); J=N+1 is
    legal but converges only as 1/nmax and is checked separately by e5b_cutoff_minimal."""
    s = mp.mpc('2.3', '0.4')
    for k in (6, 8):
        for N in (0, 1, 2):
            base = FP(N, s, k, J=N + 8)
            for J in (N + 6, N + 7):
                yield ("k=%d N=%d J=%d vs %d" % (k, N, J, N + 8), FP(N, s, k, J=J), base)


@check("rem:cutoff", "E", tier="report", desc=r"slower cutoffs $J=N+1,2,3$ converging to the same $\mathrm{FP}_N$")
def e5b_cutoff_minimal():
    r"""REPORTER.  At the minimal cutoff the remainder decays as n^{-2}, so the truncation
    error falls only as 1/nmax.  Tripling nmax should third the error -- printed as evidence
    that J=N+1 tends to the same value, not as a pass/fail."""
    s = mp.mpc('2.3', '0.4')
    for k, N in ((6, 0), (8, 1)):
        ref = FP(N, s, k, J=N + 5)
        for J in (N + 1, N + 2, N + 3):
            for nm in (100, 400):
                yield ("k=%d N=%d J=N+%d nmax=%d" % (k, N, J - N, nm),
                       FP(N, s, k, J=J, nmax=nm), ref)


@check("lem:Laurent", "E", desc=r"modularity relations on the pole weights at $i$")
def e6_laurent():
    r"""eq:Laurent: sum_{m>=r} r*(m) C(m,r) i^m = i^k sum_{m>=r} r*(m) C(k,m-r) i^{-m}."""
    for name, f, k, P in MERO_AT_I:
        a = laurent_i(f, P)
        for r in range(1, P + 1):
            lhs = sum(a[m] * mp.binomial(m, r) * I**m for m in range(r, P + 1))
            rhs = I**k * sum(a[m] * mp.binomial(k, m - r) * I**(-m) for m in range(r, P + 1))
            yield ("%s (k=%d P=%d) r=%d" % (name, k, P, r), lhs, rhs)


@check("lem:polylogedge", "E", tier="report",
       desc=r"OPEN: edge pole ${\rm Im}\,z=1$, from-above vs from-below")
def e7_edge_route():
    r"""REPORTER.  def:null_homotopy puts the T-contour at 1+delta, ABOVE a pole sitting at
    height exactly 1, so the minimal-class value is the from-BELOW one; lem:polylogedge instead
    continues the from-ABOVE formula of lem:polylogcf down to the edge.  The two differ by the
    crossing residue.  Printed, not scored, until the route is settled."""
    k, N = 12, 0
    s = mp.mpc('2.3', '0.4')
    for x0 in (mp.mpf('0.3'), mp.mpf('0.42')):
        z, eps = x0 + I, mp.mpf('0.0125')
        above = IN_raw(N, s, k, x0 + I * (1 + eps))
        below = mp.mpc(0)
        for n in range(1, (400 if SLOW else 200) + 1):
            X = 2 * pi * n * (1 + eps)
            below += (mp.e**(-2 * pi * I * n * z)
                      * (mp.gammainc(s, X, mp.inf) / (2 * pi * n)**s
                         + I**k * mp.gammainc(k - s, X, mp.inf) / (2 * pi * n)**(k - s)))
        yield ("z=%s+i : below vs above" % mp.nstr(x0, 3), below, above)


# ------------------------------------------------------------------------- driver
def main():
    print("validate.py  tier=%s  dps=%d  tol=%s" % ("SLOW" if SLOW else "fast",
                                                    mp.mp.dps, mp.nstr(TOL, 3)), **FL)
    print("=" * 96, **FL)
    npass = nfail = 0
    failures = []
    for c in CHECKS:
        if ONLY_LAYER and c["layer"] != ONLY_LAYER:
            continue
        if c["tier"] == "slow" and not SLOW:
            continue
        print("\n[%s] %-28s %s" % (c["layer"], c["label"], c["desc"]), **FL)
        for name, got, want in c["fn"]():
            r = rel(got, want)
            if c["tier"] == "report":
                # documents a known-open discrepancy; never counted as pass or fail
                print("   %-4s %-42s rel=%s" % ("note", name, mp.nstr(r, 3)), **FL)
                continue
            ok = r < (c.get("tol") or TOL)
            npass, nfail = npass + ok, nfail + (not ok)
            if not ok:
                failures.append((c["label"], name, got, want, r))
            print("   %-4s %-42s rel=%s" % ("ok" if ok else "FAIL", name, mp.nstr(r, 3)), **FL)
    print("\n" + "=" * 96, **FL)
    print("%d passed, %d failed" % (npass, nfail), **FL)
    for lab, name, got, want, r in failures:
        print("  FAIL %s :: %s\n       got  = %s\n       want = %s\n       rel  = %s"
              % (lab, name, mp.nstr(got, 12), mp.nstr(want, 12), mp.nstr(r, 3)), **FL)
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
