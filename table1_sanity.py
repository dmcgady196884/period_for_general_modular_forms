#!/usr/bin/env python3
r"""table1_sanity.py -- resolve DAM's e^{4pi} worry about Table 1.

For Delta_k = q+O(q^2) and Dhat_k = q^{-1}+O(q^2), print:
 - raw L*(f,s) at the two leading-in-X critical points s=1 (coeff of X^n) and s=2 (X^{n-1}Y),
   and (per DAM's words) at s=k-2, k-3;
 - the magnitude ratio |L*(Dhat,s)| / |L*(Delta,s)| -- is it ~e^{4pi}=1e6, or not?
 - the coefficient map: [X^n]r_f = Wp[0]*(omega^+ or eta^+); [X^{n-1}Y]r_f = Wm[0]*(omega^- or eta^-),
   confirming the Appendix numbers match Table 1.
"""
import mpmath as mp
sys_path = "/Users/dmcgady/Documents/math/meromorphic_modular_forms/quasi_period_success"
import sys; sys.path.insert(0, sys_path)
from period_polynomial_bases import Delta_k, Dhat_k, Lstar, compute
mp.mp.dps = 55; I = mp.j; pi = mp.pi
FL = dict(flush=True)

print("e^{2pi} = %.1f ,  e^{4pi} = %.3e\n" % (float(mp.e**(2*pi)), float(mp.e**(4*pi))), **FL)

for k in [12, 16, 18, 20, 22, 26]:
    n = k - 2
    Dk, Dh = Delta_k(k), Dhat_k(k)
    d = compute(k)
    Wp0, Wm0 = d["Wp"][0], d["Wm"][0]          # leading (highest-X) coeffs of W_+, W_-
    print("=== k=%d  (n=%d) ===" % (k, n), **FL)
    # raw L-values
    for s in [1, 2, 3, k-3, k-2, k-1]:
        LD, LH = Lstar(Dk, k, s), Lstar(Dh, k, s)
        ratio = abs(LH)/abs(LD) if abs(LD) > 0 else mp.inf
        print("   s=%-3d |L*(D)|=%.3e  |L*(Dh)|=%.3e   ratio=%.4g"
              % (s, float(abs(LD)), float(abs(LH)), float(ratio)), **FL)
    # coefficient map: leading real coeffs p_0 (X^n) and p_1 (X^{n-1}Y)
    def p01(F):
        c0 = I**1 * Lstar(F, k, 1)                       # a=0
        c1 = -mp.binomial(n, 1) * I**2 * Lstar(F, k, 2)  # a=1
        return c0.imag, c1.real
    p0D, p1D = p01(Dk); p0H, p1H = p01(Dh)
    print("   [X^n]   r_D=%.6e = Wp0*om+ =%.6e  |  r_Dh=%.6e = Wp0*eta+ =%.6e"
          % (float(p0D), float(Wp0*d["omp"]), float(p0H), float(Wp0*d["etp"])), **FL)
    print("   [X^{n-1}Y] r_D=%.6e = Wm0*om- =%.6e | r_Dh=%.6e = Wm0*eta- =%.6e"
          % (float(p1D), float(Wm0*d["omm"]), float(p1H), float(Wm0*d["etm"])), **FL)
    print("   Wp0=%d Wm0=%d ; om+-=(%.4e,%.4e) eta+-=(%.4e,%.4e) ; eta+/om+=%.4g"
          % (Wp0, Wm0, float(d["omp"]), float(d["omm"]), float(d["etp"]), float(d["etm"]),
             float(d["etp"]/d["omp"])), **FL)
    print("", **FL)
