"""Where does zeta(17) sit in r_{E_18}, and is the odd part rational?

theta_k18_DeltaE4.py found lambda_4/lambda_2 = -25/8, lambda_6/lambda_4 = -26/25 exactly rational,
but lambda_2/lambda_0 irrational.  Theta(z) = (2 pi i)^{n+2} B(z) - 2 pi i (1-z^n) r_E with B(z)
Laurent-rational, so r_E enters ONLY at m = 0 and m = n; hence the irrationality must live in r_E.

REDUCTION.  A single ODD coordinate l0 of c_m equals lambda_m times a fixed scalar, so ratios are
untouched and no W-basis normalisation is needed.  For m not in {0,n}, c_m = (2 pi i)^{n+2} b_m with
b_m rational (checked at k=12 in theta_laurent.py PART 4).  So the whole question is which
coefficients of r_{E_18} are rational multiples of (2 pi i)^{n+1}, and which carry zeta(17).

THE CONFOUNDER, removed here.  Every run so far used E_4^3E_6 as the reference form.  That is
admissible (Phi = 1) but it is NOT the Eisenstein series: M_18 is 2-dimensional and
E_4^3E_6 = E_18 + c Delta E_6 with c rational.  A cusp-form period is transcendental on its own, so
the contamination alone could produce the irrational ratio.  This run uses the TRUE E_18,
    E_18 = E_4^3E_6 + c Delta E_6 ,   c = -216 - 28728/43867 ,
from E_4^3E_6 = 1 + 216q + O(q^2) and E_18 = 1 - (2k/B_k) sum sigma_17(n)q^n with B_18 = 43867/798,
so 2k/B_k = 36*798/43867 = 28728/43867.  PART 1 verifies that against the q-expansion rather than
trusting the arithmetic.

PREDICTION.  L(E_k,s) = zeta(s) zeta(s-k+1), so at s = l+1 the factor is zeta(l+1) zeta(l-16):
  - l = 16 : zeta(17) zeta(0) -- the ONLY odd-zeta appearance, and l = 16 is EVEN;
  - l even, l < 16 : zeta(l-16) at a negative even integer = 0;
  - l odd : zeta(l+1) at a positive even integer is rational * pi^{l+1}, and zeta(l-16) at a
    negative odd integer is rational -- so rational * pi^{even}.
Hence the ODD coordinates of r_{E_18} should be rational multiples of (2 pi i)^{n+1}, and zeta(17)
should appear only at l = 16.  If that holds, lambda_2/lambda_0 is rational for the true Eisenstein
reference and the earlier irrationality was the Delta E_6 contamination.

Caveat kept in view: L* here is the finite-contour L-function on the arc, not the classical
completed L; the prediction above is what the classical factorisation suggests, and the point of
the run is to find out whether the finite-contour version inherits it.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, E6, Delta, rvec, arcint                  # noqa: E402

mp.mp.dps = 40
K, N = 18, 16
cE = -mp.mpf(216) - mp.mpf(28728) / 43867
E18 = lambda t: E4(t)**3 * E6(t) + cE * Delta(t) * E6(t)
E43E6 = lambda t: E4(t)**3 * E6(t)

print("PART 1  is E18 the Eisenstein series?  c = %s" % mp.nstr(cE, 16), flush=True)
for y in (mp.mpf(3), mp.mpf(4)):
    t = I * y
    q = mp.e**(2 * I * pi * t)
    print("   (E18 - 1)/q at tau=%si : %-24s  target -28728/43867 = %s"
          % (mp.nstr(y, 3), mp.nstr((E18(t) - 1) / q, 14),
             mp.nstr(-mp.mpf(28728) / 43867, 14)), flush=True)
sig17_2 = 1 + mp.mpf(2)**17
t = I * mp.mpf(3)
q = mp.e**(2 * I * pi * t)
print("   q^2 coefficient: %-24s target %s"
      % (mp.nstr((E18(t) - 1 + (mp.mpf(28728) / 43867) * q) / q**2, 12),
         mp.nstr(-(mp.mpf(28728) / 43867) * sig17_2, 12)), flush=True)
print("   Phi(E18) = %s   (must be 1)" % mp.nstr(arcint(E18, 11), 14), flush=True)

print("\nPART 2  r_{E_18}[l] / (2 pi i)^{n+1}: rational?  and where is zeta(17)?", flush=True)
rE = rvec(E18, K, 11)
z17 = mp.zeta(17)
norm = (2 * pi * I)**(N + 1)
print("   l    parity  value/(2pi i)^{n+1}            PSLQ vs 1            PSLQ vs {1, zeta(17)}",
      flush=True)
for l in range(N + 1):
    x = rE[l] / norm
    if abs(x) < mp.mpf('1e-28'):
        print("   %-4d %-7s %-30s zero" % (l, "even" if l % 2 == 0 else "ODD", "0"), flush=True)
        continue
    xr = mp.re(x) if abs(mp.im(x)) < mp.mpf('1e-25') * max(abs(mp.re(x)), mp.mpf(1)) else None
    val = xr if xr is not None else x
    s1 = s2 = "-"
    if xr is not None:
        q1 = mp.pslq([xr, mp.mpf(1)], tol=mp.mpf('1e-25'), maxcoeff=10**15, maxsteps=10**5)
        s1 = ("%s/%s" % (-q1[1], q1[0])) if q1 else "no"
        q2 = mp.pslq([xr, mp.mpf(1), z17], tol=mp.mpf('1e-25'), maxcoeff=10**12,
                     maxsteps=10**5)
        s2 = ("%s + %s z17 (/%s)" % (-q2[1], -q2[2], q2[0])) if q2 else "no"
    print("   %-4d %-7s %-30s %-20s %s"
          % (l, "even" if l % 2 == 0 else "ODD", mp.nstr(val, 16), s1, s2), flush=True)

print("\nPART 3  same for the E_4^3E_6 reference, to isolate the Delta E_6 contamination",
      flush=True)
rC = rvec(E43E6, K, 11)
for l in range(1, N, 2):                                   # odd coordinates only
    x = mp.re(rC[l] / norm)
    q1 = mp.pslq([x, mp.mpf(1)], tol=mp.mpf('1e-25'), maxcoeff=10**15, maxsteps=10**5)
    y = mp.re(rE[l] / norm)
    q2 = mp.pslq([y, mp.mpf(1)], tol=mp.mpf('1e-25'), maxcoeff=10**15, maxsteps=10**5)
    print("   l=%-3d  E_4^3E_6: %-14s   E_18: %s"
          % (l, ("%s/%s" % (-q1[1], q1[0])) if q1 else "no small rational",
             ("%s/%s" % (-q2[1], q2[0])) if q2 else "no small rational"), flush=True)
print("\n   difference r_{E_4^3E_6} - r_{E_18} = -c r_{Delta E_6}: its odd part is a cusp-form",
      flush=True)
print("   period, so rationality of the E_18 column and NOT the E_4^3E_6 column is the signal.",
      flush=True)
