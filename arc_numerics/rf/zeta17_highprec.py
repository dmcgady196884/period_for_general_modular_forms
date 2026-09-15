"""The zeta(17) question at enough digits to actually decide it.

zeta17_check.py showed: every ODD coordinate of r_{E_18}/(2 pi i)^{n+1} is rational with 43867 =
numerator(B_18) in the denominator; the even coordinates l = 2..14 vanish; and l = 0, 16 are
purely imaginary and opposite, so the even part is a multiple of p_0 = X^n - Y^n.  The only open
question is what that multiple is.

At 16 digits PSLQ returned [234241, -4914100, 4870775] against {1, zeta(17)}.  That is NOT credible:
a 3-term relation with coefficients ~5e6 needs roughly 3 * 6.7 + margin ~ 25 digits to confirm, and
16 were supplied.  This run supplies 50.

CHEAP BECAUSE ONLY TWO COEFFICIENTS ARE NEEDED.  common.rvec's convention is
    r_f[l] = (2 pi i)^{n+1} (-1)^l C(n,l) L_raw(l+1),
    L_raw(s) = int_arc f tau^{s-1} + int_arc f ktil(.,s,k),
so r_E[0]/(2 pi i)^{n+1} = L_raw(1) and r_E[n]/(2 pi i)^{n+1} = L_raw(n+1), each two quadratures
rather than the seventeen rvec does.  At s = 1 the first integrand is f itself, so that piece is
Phi(E_18) = 1 exactly -- a free internal check.

TESTED, in order of decreasing prior: x rational; x rational multiple of zeta(17); x in
Q + Q zeta(17); and the same for x/(2 pi i) and x*(2 pi i) in case the natural normalisation carries
one more factor.  Any relation whose coefficient height is too large for 50 digits is reported as
untrustworthy rather than as a hit.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, E6, Delta, ktil, arcint                  # noqa: E402

mp.mp.dps = 50
K, N = 18, 16
cE = -mp.mpf(216) - mp.mpf(28728) / 43867
E18 = lambda t: E4(t)**3 * E6(t) + cE * Delta(t) * E6(t)
DEPTH = 12


def Lraw(s):
    return (arcint(lambda t: E18(t) * t**(s - 1), DEPTH)
            + arcint(lambda t: E18(t) * ktil(t, s, K), DEPTH))


print("Phi(E_18) = %s   (internal check: the s=1 first integral)"
      % mp.nstr(arcint(E18, DEPTH), 20), flush=True)
x0 = Lraw(mp.mpf(1))
xn = Lraw(mp.mpf(N + 1))
print("r_E[0]/(2pi i)^{n+1}  = %s" % mp.nstr(x0, 30), flush=True)
print("r_E[16]/(2pi i)^{n+1} = %s" % mp.nstr(xn, 30), flush=True)
print("sum (must vanish if the even part is along p_0) = %s"
      % mp.nstr(abs(x0 + xn) / abs(x0), 6), flush=True)

z17 = mp.zeta(17)
print("\nzeta(17) = %s" % mp.nstr(z17, 30), flush=True)


def report(nm, v):
    print("\n   %s = %s" % (nm, mp.nstr(v, 30)), flush=True)
    q = mp.pslq([v, mp.mpf(1)], tol=mp.mpf('1e-40'), maxcoeff=10**15, maxsteps=10**6)
    print("      rational?            %s"
          % (("%s/%s" % (-q[1], q[0])) if q else "no"), flush=True)
    q = mp.pslq([v, z17], tol=mp.mpf('1e-40'), maxcoeff=10**15, maxsteps=10**6)
    print("      rational * zeta(17)? %s"
          % (("%s/%s" % (q[1] and -q[1], q[0])) if q else "no"), flush=True)
    q = mp.pslq([v, mp.mpf(1), z17], tol=mp.mpf('1e-38'), maxcoeff=10**10, maxsteps=10**6)
    if q:
        h = max(abs(t) for t in q)
        ok = "trustworthy" if h < 10**6 else "height %s too large for 50 digits -- DISCARD" % h
        print("      in Q + Q zeta(17)?   %s   [%s]" % (q, ok), flush=True)
    else:
        print("      in Q + Q zeta(17)?   no", flush=True)


v = -mp.im(x0) if abs(mp.re(x0)) < mp.mpf('1e-30') else x0
report("x := -Im r_E[0]/(2pi i)^{n+1}", v)
report("x / (2 pi)", v / (2 * pi))
report("x * (2 pi)", v * 2 * pi)
report("x * 43867", v * 43867)
