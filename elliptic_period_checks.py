#!/usr/bin/env python3
r"""
elliptic_period_checks.py  --  companion numerics for section 4 of
    finite_contour_cocycles_short.tex   (\cite{VerifyScripts})

PUBLISHED companion code: the ROBUST, checkable results only.  The messy / partly
contradictory working evidence (routing-dependent relation tables, etc.) lives in the
separate working ledger  wip_meromorphic_periods.py  and is NOT part of this file.

Verifies:
  [A] B(2) = 1/(864 pi)                     for  E4 Delta/E6^2  (k=4)   [eq:(ell-k4)]
  [B] B(3) = Gamma(1/4)^8/(36864 pi^6)      for  Delta/E6       (k=6)   [eq:(ell-k6)]
  [C] period polynomial in ONE fixed class: r_f|(1+U+U^2)=0, r_f|(1+S)=-2 pi a_{-2}(X^2+Y^2)
  [D] the class-dependence LAW (rem:ellrelations-open): crossing a pole p shifts
      r_f|(1+S) by EXACTLY 2 pi i * (residue polynomial of f at p)|(1+S).  This is
      cor:wall specialized to the period polynomial; it is why "which relation breaks"
      is a property of the homotopy class, not of the pole.

Conventions.  L*(f,s) = e^{-i pi s/2}[ int_{gamma^S} f tau^{s-1} + int_{gamma^T} f ktil ]
(def:Lint), ktil(tau,s)=zeta(1-s,tau+1)-e^{i pi(s-1)} zeta(1-(k-s),tau+1) (def:kernel);
B(s):=e^{i pi s/2} L*(f,s).  Period polynomial (def:rf, overall (2 pi i)^{n+1} dropped,
n=k-2): r_f = sum_l (-1)^l C(n,l) i^{l+1} L*(l+1) X^{n-l} Y^l.  Slash (paper convention):
(P|S)(X,Y)=P(-Y,X), (P|U)(X,Y)=P(X-Y,X), U^2:(X,Y)->(-Y,X-Y).
Run:  /path/to/mmf_venv/bin/python elliptic_period_checks.py   (needs mpmath).
"""
import mpmath as mp
mp.mp.dps = 34
I, pi = mp.j, mp.pi

def sigma(a, n): return sum(d**a for d in range(1, n + 1) if n % d == 0)
NT = 42
_e4 = [240 * sigma(3, n) for n in range(NT, 0, -1)] + [mp.mpf(1)]
_e6 = [-504 * sigma(5, n) for n in range(NT, 0, -1)] + [mp.mpf(1)]
def q(t): return mp.e**(2*pi*I*t)
def E4(t): return mp.polyval(_e4, q(t))
def E6(t): return mp.polyval(_e6, q(t))
def Delta(t): return (E4(t)**3 - E6(t)**2) / 1728
def jf(t): return E4(t)**3 / Delta(t)
def ktil(t, s, k): return mp.zeta(1-s, t+1) - mp.e**(I*pi*(s-1)) * mp.zeta(1-(k-s), t+1)

def Lstar(g, k, s, t0):                          # two-segment L*, straight segments at base t0
    a, b = -1/t0, t0
    S = mp.quad(lambda u: g(a + u*(b-a)) * (a + u*(b-a))**(s-1) * (b-a), [0, '0.25', '0.5', '0.75', 1])
    T = mp.quad(lambda u: g((t0-1) + u) * ktil((t0-1) + u, s, k), [0, '0.5', 1])
    return mp.e**(-I*pi*s/2) * (S + T)
def Lstar_i(g, k, s):                            # base tau0 = i: gamma^S degenerate, T over [i-1,i]
    return mp.e**(-I*pi*s/2) * mp.quad(lambda u: g((I-1)+u) * ktil((I-1)+u, s, k), [0, '0.5', 1])

def rf_coeffs(Lvals, k): return [(-1)**l*mp.binomial(k-2, l)*I**(l+1)*Lvals[l] for l in range(k-1)]
def _P(c, X, Y):
    n = len(c) - 1
    return sum(c[l] * X**(n-l) * Y**l for l in range(n+1))
_PTS = [(1,0),(0,1),(1,1),(2,1),(1,2),(3,1)]
def relS(c): return max(abs(_P(c,X,Y) + _P(c,-Y,X)) for X,Y in _PTS)
def relU(c): return max(abs(_P(c,X,Y) + _P(c,X-Y,X) + _P(c,-Y,X-Y)) for X,Y in _PTS)
def Res_at(g, z, r=mp.mpf('0.02'), Np=64):
    return sum(g(z + r*mp.e**(I*2*pi*m/Np)) * r*mp.e**(I*2*pi*m/Np) for m in range(Np)) / Np

def A_B_central_periods():
    t0 = mp.mpf('0.3') + mp.mpf('1.4')*I
    varpi = mp.gamma(mp.mpf(1)/4)**2 / (2*mp.sqrt(2*pi))
    B2 = (mp.e**(I*pi) * Lstar(lambda t: E4(t)*Delta(t)/E6(t)**2, 4, 2, t0)).real
    B3 = (mp.e**(I*pi*3/2) * Lstar(lambda t: Delta(t)/E6(t), 6, 3, t0)).real
    print("[A] k=4  B(2) - 1/(864 pi)                = %s" % mp.nstr(B2 - 1/(864*pi), 3))
    print("[B] k=6  B(3) - Gamma(1/4)^8/(36864 pi^6) = %s" % mp.nstr(B3 - varpi**4/(576*pi**4), 3))

def C_es_in_one_class():
    t0 = mp.mpf('0.3') + mp.mpf('1.4')*I
    f = lambda t: E4(t)*Delta(t)/E6(t)**2
    c = rf_coeffs([Lstar(f, 4, s, t0) for s in (1,2,3)], 4)
    a2 = sum(f(I+mp.mpf('0.02')*mp.e**(I*2*pi*m/64))*(mp.mpf('0.02')*mp.e**(I*2*pi*m/64))**2 for m in range(64))/64
    sdef = (_P(c,1,0) + _P(c,0,1))            # X^2+Y^2 coefficient of r_f|(1+S)
    print("[C] E4 D/E6^2, one class:  |(1+U+U^2)| = %s ;  S-defect coeff - (-2 pi a_{-2}) = %s" %
          (mp.nstr(relU(c), 3), mp.nstr(sdef - (-2*pi*a2), 3)))

def D_violation_is_residue():
    # crossing p=0.909i=-1/(1.1i) for f=E6/(j-j(1.1i)), k=6: r_f|(1+S) jumps by 2 pi i R|(1+S).
    jp = jf(mp.mpf('1.1')*I); k = 6; n = k-2
    f = lambda t: E6(t)/(jf(t)-jp); z = I/mp.mpf('1.1')
    cB = rf_coeffs([Lstar_i(f, k, s) for s in range(1, k)], k)
    R  = [mp.binomial(n,l)*(-1)**l*Res_at(lambda t: f(t)*t**l, z) for l in range(n+1)]
    cA = [cB[l] + 2*pi*I*R[l] for l in range(n+1)]
    pred = [2*pi*I*R[l] for l in range(n+1)]
    diff = max(abs((_P(cA,X,Y)+_P(cA,-Y,X)) - (_P(pred,X,Y)+_P(pred,-Y,X))) for X,Y in _PTS)
    print("[D] cross 0.909i:  |(1+S)| before = %s, after = %s ;  (violation - 2pi i R)|(1+S) = %s" %
          (mp.nstr(relS(cB), 3), mp.nstr(relS(cA), 3), mp.nstr(diff, 3)))

if __name__ == "__main__":
    A_B_central_periods()
    C_es_in_one_class()
    D_violation_is_residue()
