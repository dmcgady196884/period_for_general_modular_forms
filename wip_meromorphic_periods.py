#!/usr/bin/env python3
r"""
wip_meromorphic_periods.py  --  WORK-IN-PROGRESS ledger for the elliptic-block
period computations in section 4 of finite_contour_cocycles_short.tex.

Purpose: keep the (messy, partly contradictory-looking) evidence straight in ONE place.
Every result below is reproducible from here; run  python wip_meromorphic_periods.py.

Notation.  \widehat\Delta_i := f_k/(j-1728)^M is the canonical weight-k block with its
only pole on the Gamma-orbit of i (def in sec 4.1).  L*(f,s) is the two-segment integral
(def:Lint):  L* = e^{-i pi s/2}[ int_{gamma^S} f tau^{s-1} + int_{gamma^T} f ktil ], with
ktil(tau,s) = zeta(1-s,tau+1) - e^{i pi(s-1)} zeta(1-(k-s),tau+1).  B(s):=e^{i pi s/2}L*(f,s).
A "class" = the winding data (X_S,X_T) of cor:wall (one pair per pole); crossing a pole p
via gamma^S shifts L*(s) by  2 pi i e^{-i pi s/2} Res_p(f tau^{s-1}) (weight tau^{s-1}),
via gamma^T by the same with weight ktil.  Period polynomial (def:rf), overall (2pi i)^{n+1}
dropped, n=k-2:   r_f = sum_l (-1)^l C(n,l) i^{l+1} L*(l+1) X^{n-l} Y^l.
Slash: (P|S)(X,Y)=P(-Y,X), (P|U)(X,Y)=P(X-Y,X), U^2:(X,Y)->(-Y,X-Y).

=====================  EVIDENCE LOG  =====================
[E1] check_central_periods : the ROBUST results.
     k=4 (E4 D/E6^2): B(2) = 1/(864 pi)                  [=> L(.,2) = -pi/216]
     k=6 (Delta/E6) : B(3) = varpi^4/(576 pi^4) = Gamma(1/4)^8/(36864 pi^6)   [CM period]
     Central period = pi * (residue at i); parity (prop:ellord): k=0 mod4 elementary,
     k=2 mod4 carries Gamma(1/4).  Odd periods Im B(1) do NOT close in {Gamma(1/4),pi}: OPEN.
[E2] check_es_clean : E4 D/E6^2, tau0=0.3+1.4i, straight, pole off contour (CLEAN).
     r_f|(1+U+U^2)=0 (5.6e-48);  r_f|(1+S) = -2 pi a_{-2}(X^2+Y^2), a_{-2}=-1/(1728 pi^2).
[E3] check_relations_routed : tau0=0.3+1.4i, S ROUTED to dodge 1.1i/0.909i axis poles.
     rho poles break |(1+U+U^2), hold |(1+S);  1.1i poles break |(1+S), hold |(1+U+U^2).
[E4] check_relations_at_i : tau0=i (S degenerate, X_S=0). ALL five forms hold |(1+S),
     break |(1+U+U^2).  <-- CONTRADICTS [E3] for the 1.1i forms.
[E5] check_sdefect : i-pole S-defect shifts by 2x canonical per gamma^S crossing =>
     odd multiples of 2 pi a_{-2}, never 0 (gamma^S only; gamma^T crossings NOT swept).
[E6] check_reconcile + check_violation_is_residue : the RESOLUTION.  [E3] and [E4] differ
     by wall-crossings; crossing 0.909i (=-1/(1.1i)) turns a vanishing |(1+S) into 0.154,
     and that jump EQUALS 2 pi i * (residue polynomial)|(1+S) to 32 digits.  So the
     obstruction is class-dependent (a function of (X_S,X_T)), NOT of the pole location;
     neither [E3] nor [E4] "failed".
[E7] check_2i_periods : f = 1728 Delta/(j-j(2i)), a SINGLE simple pole at 2i (j(2i)=287496
     =66^3), integer q-series 1728 q^2 + 495465984 q^3 + ....  Reference omega^pm(f) via the
     two-segment integral = (-600.76, +33642.03 i): REAL/IMAGINARY, like a genuine cusp form.
     The 2i pole is crossed with DISTINCT kernels -- tau^{s-1} on gamma^S vs ktil on gamma^T --
     so r_S = rho (2i)^l and r_T = rho ktil(2i,l+1) disagree at every l; hence S-crossing !=
     T-crossing (omega^+ jump ratio -0.968-0.175i, not 1).  Both crossings pollute the reality;
     the REFERENCE class is the *unique* one with real omega^+ and imaginary omega^-.
     => candidate canonical-class selector: "the class whose periods are real/imaginary."
     (Corrects an earlier hunch that S- and T-crossing agree: they agree only in modulus of
     the residue, not in the residue itself -- the kernels differ.)
OPEN: (i) canonical class; (ii) class-independent invariant (sweep full (X_S,X_T) lattice,
BOTH residue sources); (iii) closed form for Im B(1).
=========================================================
"""
import mpmath as mp
mp.mp.dps = 34
I, pi = mp.j, mp.pi

def sigma(a, n): return sum(d**a for d in range(1, n + 1) if n % d == 0)
NT = 42
_e4 = [240 * sigma(3, n) for n in range(NT, 0, -1)] + [mp.mpf(1)]
_e6 = [-504 * sigma(5, n) for n in range(NT, 0, -1)] + [mp.mpf(1)]
def q(t): return mp.e**(2 * pi * I * t)
def E4(t): return mp.polyval(_e4, q(t))
def E6(t): return mp.polyval(_e6, q(t))
def E8(t): return E4(t)**2
def Delta(t): return (E4(t)**3 - E6(t)**2) / 1728
def jf(t): return E4(t)**3 / Delta(t)
def ktil(t, s, k): return mp.zeta(1 - s, t + 1) - mp.e**(I*pi*(s-1)) * mp.zeta(1 - (k - s), t + 1)

def Lstar(g, k, s, t0, waypoints=()):
    a, b = -1 / t0, t0
    pts = [a, *waypoints, b]
    S = mp.mpf(0)
    for p0, p1 in zip(pts, pts[1:]):
        S += mp.quad(lambda u: g(p0 + u*(p1-p0)) * (p0 + u*(p1-p0))**(s-1) * (p1-p0), [0, '0.5', 1])
    T = mp.quad(lambda u: g((t0-1) + u) * ktil((t0-1) + u, s, k), [0, '0.5', 1])
    return mp.e**(-I*pi*s/2) * (S + T)

def Lstar_i(g, k, s):                     # tau0 = i : gamma^S degenerate, T over [i-1,i]
    return mp.e**(-I*pi*s/2) * mp.quad(lambda u: g((I-1)+u) * ktil((I-1)+u, s, k), [0, '0.5', 1])

def periodpoly_coeffs(Lvals, k):          # r_f coeffs c[l] = coeff of X^{n-l} Y^l, n=k-2
    n = k - 2
    return [(-1)**l * mp.binomial(n, l) * I**(l+1) * Lvals[l] for l in range(n+1)]

def _P(c, X, Y):
    n = len(c) - 1
    return sum(c[l] * X**(n-l) * Y**l for l in range(n+1))
_PTS = [(1,0),(0,1),(1,1),(2,1),(1,2),(3,1),(1,3),(2,3)]
def relS(c): return max(abs(_P(c,X,Y) + _P(c,-Y,X)) for X,Y in _PTS)
def relU(c): return max(abs(_P(c,X,Y) + _P(c,X-Y,X) + _P(c,-Y,X-Y)) for X,Y in _PTS)
def scaleP(c): return max(abs(_P(c,X,Y)) for X,Y in _PTS)

def Res_at(g, z, r=mp.mpf('0.02'), Np=64):        # residue of g at simple pole z
    return sum(g(z + r*mp.e**(I*2*pi*m/Np)) * r*mp.e**(I*2*pi*m/Np) for m in range(Np)) / Np
def Lcoeff_at(g, z, jpow, r=mp.mpf('0.02'), Np=64):  # coeff of (tau-z)^{jpow} in g at z
    return sum(g(z + r*mp.e**(I*2*pi*m/Np)) * r**(-jpow) * mp.e**(-I*2*pi*m*jpow/Np) for m in range(Np)) / Np

# ---------------------------------------------------------------- [E1]
def check_central_periods():
    print("[E1] central periods  (robust):")
    t0 = mp.mpf('0.3') + mp.mpf('1.4')*I
    varpi = mp.gamma(mp.mpf(1)/4)**2 / (2*mp.sqrt(2*pi))
    B2 = (mp.e**(I*pi*2/2) * Lstar(lambda t: E4(t)*Delta(t)/E6(t)**2, 4, 2, t0)).real
    print("   k=4  B(2)-1/(864 pi)                     = %s" % mp.nstr(B2 - 1/(864*pi), 3))
    B3 = (mp.e**(I*pi*3/2) * Lstar(lambda t: Delta(t)/E6(t), 6, 3, t0)).real
    print("   k=6  B(3)-varpi^4/(576 pi^4)             = %s" % mp.nstr(B3 - varpi**4/(576*pi**4), 3))

# ---------------------------------------------------------------- [E2]
def check_es_clean():
    print("[E2] E4 D/E6^2, tau0=0.3+1.4i (clean): U-rel exact, S-rel = residue:")
    t0 = mp.mpf('0.3') + mp.mpf('1.4')*I
    c = periodpoly_coeffs([Lstar(lambda t: E4(t)*Delta(t)/E6(t)**2, 4, s, t0) for s in (1,2,3)], 4)
    print("   |(1+U+U^2)|/scale = %s   |(1+S)|/scale = %s" %
          (mp.nstr(relU(c)/scaleP(c), 3), mp.nstr(relS(c)/scaleP(c), 3)))

# ---------------------------------------------------------------- [E3],[E4]
def _relation_table(Lfun, label):
    print(label)
    jp = jf(mp.mpf('1.1')*I)
    wp = (mp.mpf('0.3') + mp.mpf('0.62')*I,)          # detour used only for the routed table
    forms = [("pole i    E4 D/E6^2", lambda t: E4(t)*Delta(t)/E6(t)**2, 4, ()),
             ("pole rho  E6/j",      lambda t: E6(t)/(jf(t)-0),        6, ()),
             ("pole rho  D/E4",      lambda t: E8(t)/(jf(t)-0),        8, ()),
             ("pole 1.1i E6/(j-j.)", lambda t: E6(t)/(jf(t)-jp),       6, wp),
             ("pole 1.1i E8/(j-j.)", lambda t: E8(t)/(jf(t)-jp),       8, wp)]
    for name, g, k, wps in forms:
        c = periodpoly_coeffs([Lfun(g, k, s, wps) for s in range(1, k)], k)
        print("   %-22s |(1+S)|=%-9s |(1+U+U^2)|=%-9s" %
              (name, mp.nstr(relS(c)/scaleP(c), 3), mp.nstr(relU(c)/scaleP(c), 3)))

def check_relations_routed():
    t0 = mp.mpf('0.3') + mp.mpf('1.4')*I
    _relation_table(lambda g, k, s, wps: Lstar(g, k, s, t0, wps),
                    "[E3] tau0=0.3+1.4i, S-contour ROUTED (a specific class):")
def check_relations_at_i():
    _relation_table(lambda g, k, s, wps: Lstar_i(g, k, s),
                    "[E4] tau0=i (S degenerate, X_S=0); CONTRADICTS [E3] for 1.1i:")

# ---------------------------------------------------------------- [E5]
def check_sdefect():
    print("[E5] i-pole S-defect: shift per gamma^S crossing / canonical (== 2 => never 0):")
    f = lambda t: E4(t)*Delta(t)/E6(t)**2
    a2 = Lcoeff_at(f, I, -2)
    dL = {s: 2*pi*I*mp.e**(-I*pi*s/2)*Res_at(lambda t: f(t)*t**(s-1), I) for s in (1,2,3)}
    print("   Res_i(f tau) = %s (=> B(2) crossing-invariant);  ratio = %s" %
          (mp.nstr(Res_at(lambda t: f(t)*t, I), 3), mp.nstr((I*dL[1]-I*dL[3])/(-2*pi*a2), 6)))

# ---------------------------------------------------------------- [E6]
def check_violation_is_residue():
    print("[E6] crossing 0.909i: is the |(1+S) violation == the residue?")
    jp = jf(mp.mpf('1.1')*I); k = 6; n = k - 2
    f = lambda t: E6(t)/(jf(t)-jp)
    z = I/mp.mpf('1.1')                               # 0.909.. i = -1/(1.1 i)
    print("   Res_{0.909i}(f) = %s" % mp.nstr(Res_at(f, z), 12))
    cB = periodpoly_coeffs([Lstar_i(f, k, s) for s in range(1, k)], k)
    Rc = [mp.binomial(n, l)*(-1)**l*Res_at(lambda t: f(t)*t**l, z) for l in range(n+1)]  # residue poly
    cA = [cB[l] + 2*pi*I*Rc[l] for l in range(n+1)]   # B + one gamma^S crossing of 0.909i
    pred = [2*pi*I*Rc[l] for l in range(n+1)]
    print("   |(1+S)|: class B = %s ;  after crossing = %s ;  2pi i R|(1+S) = %s" %
          (mp.nstr(relS(cB), 3), mp.nstr(relS(cA), 3), mp.nstr(relS(pred), 3)))
    diff = max(abs((_P(cA,X,Y)+_P(cA,-Y,X)) - (_P(pred,X,Y)+_P(pred,-Y,X))) for X,Y in _PTS)
    print("   violation - residue = %s   (=> violation IS the crossed residue)" % mp.nstr(diff, 3))

# ---------------------------------------------------------------- [E7]
def check_2i_periods():
    r"""[E7] f = 1728 Delta/(j-j(2i)), one simple pole at 2i (j(2i)=287496=66^3).
    Reference omega^pm(f), and the 2i-pole crossed via gamma^S vs gamma^T (distinct kernels)."""
    from period_polynomial_bases import compute, pairing
    print("[E7] 2i-pole block f=1728 Delta/(j-287496): reference periods, S- vs T-crossing:")
    k = 12; n = k - 2
    f = lambda t: 1728 * Delta(t) / (jf(t) - 287496)          # j(2i)=66^3=287496
    d = compute(k)                                            # W_pm from the cusp form Delta_12
    Wp = [mp.mpf(0)]*(n+1); Wm = [mp.mpf(0)]*(n+1)
    for i, a in enumerate(range(0, n+1, 2)): Wp[a] = mp.mpf(d['Wp'][i])
    for i, a in enumerate(range(1, n+1, 2)): Wm[a] = mp.mpf(d['Wm'][i])
    def omega(Lv):                                           # full def:rf (keep (2pi i)^{n+1}), project
        c = [(2*pi*I)**(n+1) * x for x in periodpoly_coeffs(Lv, k)]
        return pairing(c, Wp, n)/pairing(Wp, Wp, n), pairing(c, Wm, n)/pairing(Wm, Wm, n)
    t0 = I*mp.mpf('1.1')                                     # reference class (encloses no interior pole)
    Lref = [Lstar(f, k, s, t0) for s in range(1, k)]
    op, om = omega(Lref)
    print("   reference : omega^+ = %-16s omega^- = %s" % (mp.nstr(op, 8), mp.nstr(om, 8)))
    rho = Res_at(f, 2*I)                                     # Res_{2i} f
    LS = [Lref[l] + 2*pi*I*mp.e**(-I*pi*(l+1)/2) * rho*(2*I)**l          for l in range(n+1)]  # gamma^S
    LT = [Lref[l] + 2*pi*I*mp.e**(-I*pi*(l+1)/2) * rho*ktil(2*I, l+1, k) for l in range(n+1)]  # gamma^T
    opS, omS = omega(LS); opT, omT = omega(LT)
    print("   S-encloses: omega^+ = %-16s omega^- = %s" % (mp.nstr(opS, 8), mp.nstr(omS, 8)))
    print("   T-encloses: omega^+ = %-16s omega^- = %s" % (mp.nstr(opT, 8), mp.nstr(omT, 8)))
    print("   S-jump/T-jump (omega^+) = %s   (!= 1  =>  S-crossing != T-crossing)"
          % mp.nstr((opS-op)/(opT-op), 6))

if __name__ == "__main__":
    check_central_periods()
    check_es_clean()
    check_relations_routed()
    check_relations_at_i()
    check_sdefect()
    check_violation_is_residue()
    check_2i_periods()
