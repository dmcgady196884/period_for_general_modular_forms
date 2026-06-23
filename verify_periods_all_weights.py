"""
Quasi-period extraction at every weight with dim S_k = 1:  k in {12,16,18,20,22,26}.

Generalizes verify_periods.py (which did k=12 only).  For each such k:
  * cusp form     f_k   = E_{k-12} * Delta            = q + O(q^2)  in S_k
  * weak cusp form Dhat_k = f_k * (j^2 + b1 j + b0)   = 1/q + O(q^2) in S_k^!
    (Duke-Jenkins: choose b1,b0 to kill the q^0,q^1 modes)
  * rational period-polynomial basis P^+,P^- in W_{k-2}^pm, from the period
    relations p|(1+S)=0, p|(1+U+U^2)=0, U=ST, split by the Y->-Y parity
  * two-segment cocycle extraction of omega^pm(f_k), eta^pm(Dhat_k)

k=12 is cross-checked against Brown 1710.07912; k=16..26 are computed.
"""

import mpmath as mp
from sympy import Rational, symbols, expand, Poly, Matrix, nsimplify
from math import comb

mp.mp.dps = 50
PREC_Q = 80
X, Y = symbols('X Y')

# ----------------------------------------------------------------------------
# q-series machinery (rationals)
# ----------------------------------------------------------------------------

def lmul(a, b, lo, hi):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = ka + kb
            if lo <= k <= hi:
                out[k] = out.get(k, Rational(0)) + va * vb
    return {k: v for k, v in out.items() if v != 0}

def ldiv(num, den, lo_exp, hi_exp):
    den_min = min(den.keys()); den_lead = den[den_min]; quot = {}
    for target in range(lo_exp + den_min, hi_exp + den_min + 1):
        s = num.get(target, Rational(0))
        for kq in range(lo_exp, target - den_min):
            s -= quot.get(kq, Rational(0)) * den.get(target - kq, Rational(0))
        if s != 0:
            quot[target - den_min] = s / den_lead
    return quot

def sigma(nn, kk):
    s, d = 0, 1
    while d * d <= nn:
        if nn % d == 0:
            s += d ** kk
            if d != nn // d:
                s += (nn // d) ** kk
        d += 1
    return s

def eta24(prec):
    series = {0: Rational(1)}
    for nn in range(1, prec + 1):
        factor = {}
        for kk in range(25):
            te = nn * kk
            if te > prec:
                break
            factor[te] = Rational((-1) ** kk * comb(24, kk))
        series = lmul(series, factor, 0, prec)
    return series

# Delta = q * prod(1-q^n)^24
_e24 = eta24(PREC_Q + 4)
Delta = {k + 1: v for k, v in _e24.items() if 1 <= k + 1 <= PREC_Q + 4}

E4 = {0: Rational(1)}
for nn in range(1, PREC_Q + 5):
    E4[nn] = Rational(240 * sigma(nn, 3))
E6 = {0: Rational(1)}
for nn in range(1, PREC_Q + 5):
    E6[nn] = Rational(-504 * sigma(nn, 5))

E4c = lmul(lmul(E4, E4, 0, PREC_Q + 4), E4, 0, PREC_Q + 4)
j = ldiv(E4c, Delta, -1, PREC_Q + 2)          # j = 1/q + 744 + ...

# Eisenstein E_{k-12} as ring products
E8  = lmul(E4, E4, 0, PREC_Q + 4)
E10 = lmul(E4, E6, 0, PREC_Q + 4)
E14 = lmul(E8, E6, 0, PREC_Q + 4)
EIS = {12: {0: Rational(1)}, 16: E4, 18: E6, 20: E8, 22: E10, 26: E14}

def jpow(m):
    """j^m truncated."""
    out = {0: Rational(1)}
    for _ in range(m):
        out = lmul(out, j, -m, PREC_Q + 2)
    return out

def cusp_and_weak(k):
    """Return (f_k, Dhat_k) q-series dicts: f_k=q+O(q^2), Dhat_k=1/q+O(q^2)."""
    fk = lmul(EIS[k], Delta, 1, PREC_Q)
    # g = j^2 + b1 j + b0 ;  fk*g should be 1/q + O(q^2).
    j2 = jpow(2)
    A = lmul(fk, j2, -1, PREC_Q)   # leading q^{-1}
    Bj = lmul(fk, j, -1, PREC_Q)
    # solve: [A + b1 Bj + b0 fk]_{q^0} = 0, _{q^1} = 0
    # unknowns b1,b0
    M = Matrix([[Bj.get(0, Rational(0)), fk.get(0, Rational(0))],
                [Bj.get(1, Rational(0)), fk.get(1, Rational(0))]])
    rhs = Matrix([-A.get(0, Rational(0)), -A.get(1, Rational(0))])
    b1, b0 = M.solve(rhs)
    g = {}
    for d_, v in j2.items(): g[d_] = g.get(d_, Rational(0)) + v
    for d_, v in j.items():  g[d_] = g.get(d_, Rational(0)) + b1 * v
    g[0] = g.get(0, Rational(0)) + b0
    Dhat = lmul(fk, {d_: v for d_, v in g.items() if v != 0}, -1, PREC_Q)
    return fk, Dhat

# ----------------------------------------------------------------------------
# rational period-polynomial basis  W_{k-2}^pm  from the period relations
# ----------------------------------------------------------------------------

def slash(p_expr, a, b, c, d):
    return expand(p_expr.subs({X: a * X + b * Y, Y: c * X + d * Y}, simultaneous=True))

def period_basis(n):
    """Interior coefficient lists P_plus[j], P_minus[j], j=0..n, in the
    Y->-Y parity eigenbasis of the cuspidal period polynomials."""
    coeffs = symbols(f'c0:{n+1}')
    p = sum(coeffs[jj] * X**(n - jj) * Y**jj for jj in range(n + 1))
    pS = slash(p, 0, -1, 1, 0)                 # S: X->-Y, Y->X
    pU = slash(p, 0, -1, 1, 1)                 # U=ST: X->-Y, Y->X+Y
    pU2 = expand(pU.subs({X: -Y, Y: X + Y}, simultaneous=True))
    rels = expand(p + pS), expand(p + pU + pU2)
    eqs = []
    for r in rels:
        pol = Poly(r, X, Y)
        for mono, co in pol.terms():
            eqs.append(co)
    A = Matrix([[eq.coeff(c) for c in coeffs] for eq in eqs])
    ns = A.nullspace()
    # parity operator eps: Y -> -Y
    def parity_split(vec):
        pv = sum(vec[jj] * X**(n - jj) * Y**jj for jj in range(n + 1))
        ev = expand((pv + pv.subs(Y, -Y)) / 2)
        od = expand((pv - pv.subs(Y, -Y)) / 2)
        return ev, od
    plus_polys, minus_polys = [], []
    for vec in ns:
        ev, od = parity_split(list(vec))
        if ev != 0: plus_polys.append(ev)
        if od != 0: minus_polys.append(od)
    def to_interior_list(poly_expr):
        if poly_expr == 0:
            return None
        pol = Poly(poly_expr, X, Y)
        lst = [Rational(0)] * (n + 1)
        for (ex, ey), co in pol.terms():
            lst[ey] = Rational(co)
        return lst
    # cuspidal pieces: take a basis poly that has nonzero INTERIOR support
    def pick(polys, idxs):
        for pe in polys:
            lst = to_interior_list(pe)
            if lst and any(lst[jj] != 0 for jj in idxs):
                # primitive-integer normalize on interior
                from sympy import lcm, igcd
                dens = [lst[jj].q for jj in range(n + 1) if lst[jj] != 0]
                L = 1
                for dd in dens: L = L * dd // __import__('math').gcd(L, dd)
                ints = [int(lst[jj] * L) for jj in range(n + 1)]
                g = 0
                for jj in idxs:
                    g = __import__('math').gcd(g, abs(ints[jj]))
                if g == 0: g = 1
                res = [Rational(ints[jj], g) for jj in range(n + 1)]
                for jj in idxs:          # normalize: first nonzero interior coeff > 0
                    if res[jj] != 0:
                        if res[jj] < 0:
                            res = [-x for x in res]
                        break
                return res
        return None
    interior_plus = [jj for jj in range(1, n) if jj % 2 == 0]
    interior_minus = [jj for jj in range(1, n) if jj % 2 == 1]
    Pp = pick(plus_polys, interior_plus)
    Pm = pick(minus_polys, interior_minus)
    return Pp, Pm, interior_plus, interior_minus

# ----------------------------------------------------------------------------
# cocycle extraction (mpmath)
# ----------------------------------------------------------------------------

def make_eval(series):
    terms = sorted([(int(nn), mp.mpf(str(c))) for nn, c in series.items() if c != 0])
    TWO = mp.mpc(0, 2) * mp.pi
    def ev(tau):
        tau = mp.mpc(tau); val = mp.mpc(0)
        for nn, c in terms:
            val += c * mp.exp(nn * TWO * tau)
        return val
    return ev

def cocycle(f, ta, tb, n, pref):
    ta = mp.mpc(ta); tb = mp.mpc(tb); dtau = tb - ta
    out = []
    for jj in range(n + 1):
        out.append(mp.quad(lambda s, jj=jj: f(ta + s * dtau) * (ta + s * dtau)**jj * dtau, [0, 1]))
    return [pref * mp.mpf(comb(n, jj)) * mp.mpf((-1)**jj) * out[jj] for jj in range(n + 1)]

def slash_S_num(p, n):
    return [p[n - m] * mp.mpf((-1)**m) for m in range(n + 1)]

def solve_dT(C_T, n):
    A = mp.matrix(n, n); b = mp.matrix(n, 1)
    for m in range(1, n + 1):
        for jj in range(m):
            A[m - 1, jj] = mp.mpf(comb(n - jj, m - jj))
        b[m - 1, 0] = C_T[m]
    ps = mp.lu_solve(A, b)
    return [ps[jj, 0] for jj in range(n)] + [mp.mpc(0)]

def extract(f, tau0, n, pref, Pp, Pm, ip, im):
    C_S = cocycle(f, -1 / tau0, tau0, n, pref)
    C_T = cocycle(f, tau0 - 1, tau0, n, pref)
    P = solve_dT(C_T, n)
    PS = slash_S_num(P, n)
    r = [C_S[jj] - (PS[jj] - P[jj]) for jj in range(n + 1)]
    plus = [r[jj] / mp.mpf(str(Pp[jj])) for jj in ip if Pp[jj] != 0]
    minus = [r[jj] / mp.mpf(str(Pm[jj])) for jj in im if Pm[jj] != 0]
    # cross-check: where the rational basis coeff vanishes -- the central
    # critical L-value forced to 0 by i^k=-1 (k=2 mod 4) -- the numerical
    # period coefficient must vanish too.
    central = ([abs(r[jj]) for jj in ip if Pp[jj] == 0]
               + [abs(r[jj]) for jj in im if Pm[jj] == 0])
    central_res = max(central) if central else mp.mpf(0)
    lock_p = max(abs(v - plus[0]) for v in plus)
    lock_m = max(abs(v - minus[0]) for v in minus)
    return plus[0], minus[0], abs(C_T[0]), lock_p, lock_m, central_res

BROWN = {12: ('-68916772.80959519475431012465533103043907',
              '-5585015.379310401866877139263796275129635',
              '127202100647.1770947773171612986108774951',
              '10276732343.64913275081719307240092090893')}

if __name__ == '__main__':
    tau0 = mp.mpc('0.3', '1.2')
    tau0b = mp.mpc('-0.2', '0.9')
    print(f"dps={mp.mp.dps}, PREC_Q={PREC_Q}, tau0={tau0}\n")
    for k in [12, 16, 18, 20, 22, 26]:
        n = k - 2
        pref = (mp.mpc(0, 2) * mp.pi) ** (n + 1)
        fk, Dhat = cusp_and_weak(k)
        Pp, Pm, ip, im = period_basis(n)
        fev, hev = make_eval(fk), make_eval(Dhat)
        op, om, ct, lp, lm, cv = extract(fev, tau0, n, pref, Pp, Pm, ip, im)
        ep, em, hct, hlp, hlm, hcv = extract(hev, tau0, n, pref, Pp, Pm, ip, im)
        # basepoint independence on Dhat
        ep2, em2, *_ = extract(hev, tau0b, n, pref, Pp, Pm, ip, im)
        bp = max(abs(ep - ep2), abs(em - em2))
        print(f"=== k={k}  (n={n}, dim W^+ interior {ip}, W^- interior {im}) ===")
        print(f"  P^+ interior: {[str(Pp[jj]) for jj in ip]}")
        print(f"  P^- interior: {[str(Pm[jj]) for jj in im]}")
        print(f"  cuspidality |C_T[X^n]|:  f_k {float(ct):.2e}   Dhat {float(hct):.2e}")
        print(f"  lockstep (worst):        f_k {float(max(lp,lm)):.2e}   Dhat {float(max(hlp,hlm)):.2e}")
        if k % 4 == 2:
            print(f"  central L-value vanishes (i^k=-1): |[r_f]_central|  f_k {float(cv):.2e}   Dhat {float(hcv):.2e}")
        print(f"  basepoint indep (Dhat):  {float(bp):.2e}")
        print(f"  omega^+(f_k) = {op}")
        print(f"  omega^-(f_k) = {om}")
        print(f"  eta^+(Dhat)  = {ep}")
        print(f"  eta^-(Dhat)  = {em}")
        # Legendre-type cross-check: det of [[w+,w-],[e+,e-]] should be a clean
        # rational multiple of (2 pi i)^power (catches a global normalization slip).
        det = op * em - om * ep
        tpi = mp.mpc(0, 2) * mp.pi
        for pw in (k - 1, k, k - 2):
            ratio = det / tpi**pw
            try:
                rid = nsimplify(mp.nstr(ratio.real, 30), rational=True)
            except Exception:
                rid = '?'
            print(f"  det/(2pi i)^{pw} = {mp.nstr(ratio, 26)}   (re ~ {rid})")
        if k in BROWN:
            bop, bom, bep, bem = BROWN[k]
            tgt = [mp.mpf(bop), mp.mpc(0, 1) * mp.mpf(bom), mp.mpf(bep), mp.mpc(0, 1) * mp.mpf(bem)]
            for name, got, t in [("omega+", op, tgt[0]), ("omega-", om, tgt[1]),
                                 ("eta+", ep, tgt[2]), ("eta-", em, tgt[3])]:
                print(f"    {name} vs Brown: {float(abs(got - t) / abs(t)):.2e}")
        print()
