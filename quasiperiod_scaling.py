#!/usr/bin/env python3
r"""quasiperiod_scaling.py -- test eta^pm -> (mp)omega^pm and the Legendre determinant,
for the six weights with dim S_k = 1.  Emits readable + LaTeX tables."""
import mpmath as mp, sys
sys.path.insert(0, "/Users/dmcgady/Documents/math/meromorphic_modular_forms/quasi_period_success")
from period_polynomial_bases import compute
mp.mp.dps = 55; I = mp.j; pi = mp.pi
latex = len(sys.argv) > 1 and sys.argv[1] == "latex"

def frac(x, maxden=10**8):
    f = mp.pslq([x, mp.mpf(1)], maxcoeff=maxden, tol=mp.mpf(10)**(-40))
    if f and f[0] != 0:
        from fractions import Fraction
        fr = Fraction(int(-f[1]), int(f[0]))
        if abs(float(fr) - float(x)) < 1e-30: return fr
    return None

rows = []
for k in [12, 16, 18, 20, 22, 26]:
    d = compute(k)
    op, om, ep, em = d["omp"], d["omm"], d["etp"], d["etm"]
    rp, rm = ep/op, em/om                       # eta^+/omega^+, eta^-/omega^-
    det = op*em - om*ep                          # Legendre determinant (W_pm primitive basis)
    # is det rational?  is det/(2pi)^{k-1} rational?  is det*(2pi)^{k-1}?
    dr = frac(det); dr2 = frac(det*(2*pi)**(k-1)); dr3 = frac(det/(2*pi)**(k-1))
    rows.append((k, float(rp), float(rm), float(det), dr, dr2, dr3))

if not latex:
    print("%3s %14s %14s | %14s  rational?  det*(2pi)^{k-1}?  det/(2pi)^{k-1}?" %
          ("k", "eta+/om+", "eta-/om-", "det=om+et- - om-et+"))
    for (k, rp, rm, det, dr, dr2, dr3) in rows:
        print("%3d %14.6g %14.6g | %14.6g   %-10s  %-14s  %-14s" %
              (k, rp, rm, det, str(dr), str(dr2), str(dr3)))
    print("\npredicted limits: eta+/om+ -> -1 , eta-/om- -> +1")
else:
    print(r"\begin{tabular}{r|rr|rr}")
    print(r"\hline")
    print(r"$k$ & $\eta^+/\omega^+$ & $\eta^-/\omega^-$ & $|\eta^+/\omega^++1|$ & $|\eta^-/\omega^--1|$\\")
    print(r"\hline")
    for (k, rp, rm, det, dr, dr2, dr3) in rows:
        print(r"%d & $%s$ & $%s$ & $%s$ & $%s$\\" %
              (k, mp.nstr(rp, 6), mp.nstr(rm, 6), mp.nstr(abs(rp+1), 3), mp.nstr(abs(rm-1), 3)))
    print(r"\hline")
    print(r"\end{tabular}")
