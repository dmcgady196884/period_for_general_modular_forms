"""
Numerical verification that the explicit period formulas
\eqref{eq:omega-minus-explicit} and \eqref{eq:omega-plus-explicit}
in analytic_periods_weight12.tex hold to working precision.

Specifically:
  ω−(Δ) = -(5/2) i (2π)^9 L(Δ,2) = -(144/5) i (2π)^7 L(Δ,4) = -720 i (2π)^5 L(Δ,6)
                                = -24192 i (2π)^3 L(Δ,8) = -1814400 i π L(Δ,10)
  ω+(Δ) = -90 (2π)^8 L(Δ,3)  = -1680 (2π)^6 L(Δ,5)  = -50400 (2π)^4 L(Δ,7)
                                = -1814400 (2π)^2 L(Δ,9)
"""
from mpmath import mp, mpc, mpf, pi, exp, gamma, gammainc
from sympy import Rational
from math import comb

mp.dps = 45

# Build Δ q-series rationals (as in main notebook)
def laurent_mul(a, b, lo, hi):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = ka + kb
            if lo <= k <= hi:
                out[k] = out.get(k, Rational(0)) + va * vb
    return {k: v for k, v in out.items() if v != 0}

def eta24_qexp(prec):
    series = {0: Rational(1)}
    for n in range(1, prec + 1):
        factor = {}
        for k in range(25):
            te = n * k
            if te > prec: break
            factor[te] = Rational((-1)**k * comb(24, k))
        series = laurent_mul(series, factor, 0, prec)
    return series

PREC = 100
eta24 = eta24_qexp(PREC + 2)
Delta_series = {k+1: v for k, v in eta24.items() if 1 <= k+1 <= PREC + 2}

def L_delta(s, N=100):
    """L(Δ, s) via BFK incomplete-gamma formula at split-point t0=1."""
    two_pi = 2 * pi
    L_star = mpc(0)
    for n in range(1, N + 1):
        c = Delta_series.get(n)
        if c is None: continue
        c = mpf(str(c))
        x = two_pi * n
        L_star += c * (gammainc(s, x) / x**s + gammainc(12 - s, x) / x**(12 - s))
    return two_pi**s / gamma(s) * L_star

# Brown periods (Brown 1710.07912 §8.3, full precision from notebook)
omega_plus_target  = mpf('-68916772.80959519475431012465533103043907')
omega_minus_target = mpc(0, 1) * mpf('-5585015.379310401866877139263796275129635')

two_pi = 2 * pi
i      = mpc(0, 1)

print("=" * 70)
print("VERIFICATION OF eq:omega-minus-explicit (5 equivalent forms)")
print("=" * 70)
# All five forms for ω⁻
formulas_minus = [
    (2,  -mpf(5)/2          * i * two_pi**9, "j=1: -(5/2) i (2π)^9 · L(Δ,2)"),
    (4,  -mpf(144)/5        * i * two_pi**7, "j=3: -(144/5) i (2π)^7 · L(Δ,4)"),
    (6,  -mpf(720)          * i * two_pi**5, "j=5: -720 i (2π)^5 · L(Δ,6)"),
    (8,  -mpf(24192)        * i * two_pi**3, "j=7: -24192 i (2π)^3 · L(Δ,8)"),
    (10, -mpf(1814400)      * i * pi,        "j=9: -1814400 i π · L(Δ,10)"),
]
print()
for s, coeff, label in formulas_minus:
    Lval     = L_delta(s)
    computed = coeff * Lval
    rel_err  = abs(computed - omega_minus_target) / abs(omega_minus_target)
    print(f"  {label}")
    print(f"      L(Δ,{s})   = {mp.nstr(Lval.real, 22)}")
    print(f"      computed  = {mp.nstr(computed.imag, 22)} i")
    print(f"      rel error = {float(rel_err):.3e}")
    print()

print("=" * 70)
print("VERIFICATION OF eq:omega-plus-explicit (4 equivalent forms)")
print("=" * 70)
# All four forms for ω⁺
formulas_plus = [
    (3,  -mpf(90)        * two_pi**8, "j=2: -90 (2π)^8 · L(Δ,3)"),
    (5,  -mpf(1680)      * two_pi**6, "j=4: -1680 (2π)^6 · L(Δ,5)"),
    (7,  -mpf(50400)     * two_pi**4, "j=6: -50400 (2π)^4 · L(Δ,7)"),
    (9,  -mpf(1814400)   * two_pi**2, "j=8: -1814400 (2π)^2 · L(Δ,9)"),
]
print()
for s, coeff, label in formulas_plus:
    Lval     = L_delta(s)
    computed = coeff * Lval
    rel_err  = abs(computed - omega_plus_target) / abs(omega_plus_target)
    print(f"  {label}")
    print(f"      L(Δ,{s})   = {mp.nstr(Lval.real, 22)}")
    print(f"      computed  = {mp.nstr(computed.real, 22)}")
    print(f"      rel error = {float(rel_err):.3e}")
    print()
