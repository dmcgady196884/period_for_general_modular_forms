"""Closed form for the model integral of eq:elllimit, and the parity factor at i.

RAY (tau_e in SL2(Z).rho, an endpoint of the arc):  C is a ray from 0 to infinity in
direction omega.  With alpha = (a+1)/n, substituting x = omega^n u^n gives

    int_C u^a du / (u^n - 1)^m  =  (-1)^{m+1}/n * e^{i pi alpha} * B(alpha, m - alpha),

independent of omega except through the branch of e^{i pi alpha} -- which is exactly the
n-th root of unity ambiguity of the lemma.  Converges iff Re(m - alpha) > 0, i.e.
m n > a + 1, i.e. P = m n - a > 1.

LINE (tau_e = i, interior of the arc):  C is a full line through 0.  Sending u -> -v on
the second ray gives, at n = 2,

    int_line = ( 1 + (-1)^a ) * int_ray,

so the leading term VANISHES for odd a, hence for odd P = 2m - a.  If that is right, then
at such poles eta^{(P-1)/2} L*(f_eta,s) -> 0 and (P-1)/n is NOT the least exponent of
def:arcsplit.

Both claims are pure calculus; checked here against direct quadrature.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi                                            # noqa: E402

mp.mp.dps = 30
RMAX = mp.mpf('1e7')          # ray truncation; integrand decays like u^{-P}, P >= 2


def ray_quad(n, a, m, omega):
    """int_0^{infty omega} u^a du / (u^n - 1)^m by direct quadrature"""
    def g(t):
        u = omega * t
        return u**a / (u**n - 1)**m * omega
    pts = [mp.mpf(0)] + [mp.mpf(10)**e for e in range(-3, 8)]
    return sum(mp.quad(g, [pts[i], pts[i + 1]]) for i in range(len(pts) - 1))


def ray_closed(n, a, m):
    al = mp.mpf(a + 1) / n
    return (-1)**(m + 1) / mp.mpf(n) * mp.e**(I * pi * al) * mp.beta(al, m - al)


def line_quad(n, a, m, omega):
    """int over the full line through 0 in direction omega, from -inf.omega to +inf.omega"""
    return ray_quad(n, a, m, omega) - ray_quad(n, a, m, -omega)


print("dps = %d\n" % mp.mp.dps, flush=True)
print("RAY: closed form vs quadrature, several directions omega = e^{i psi}", flush=True)
print("  n  a  m   P    psi      |quad|          |closed|        rel.diff", flush=True)
for (n, a, m) in [(3, 0, 1), (3, 1, 1), (3, 2, 2), (3, 0, 2), (2, 0, 1), (2, 1, 2)]:
    P = m * n - a
    if P <= 1:
        continue
    for psi in ('0.7', '1.9', '2.9'):
        om = mp.e**(I * mp.mpf(psi))
        q = ray_quad(n, a, m, om)
        c = ray_closed(n, a, m)
        print("  %d  %d  %d  %2d   %-7s  %-14s  %-14s  %s"
              % (n, a, m, P, psi, mp.nstr(abs(q), 10), mp.nstr(abs(c), 10),
                 mp.nstr(abs(q - c) / abs(c), 6)), flush=True)

print("\nLINE at n = 2: is int_line = (1 + (-1)^a) int_ray?", flush=True)
print("  a  m   P    psi      |line|          predicted       note", flush=True)
for (a, m) in [(0, 1), (1, 2), (0, 2), (1, 3)]:
    n = 2
    P = m * n - a
    if P <= 1:
        continue
    for psi in ('0.7', '1.9'):
        om = mp.e**(I * mp.mpf(psi))
        lq = line_quad(n, a, m, om)
        pred = (1 + (-1)**a) * ray_quad(n, a, m, om)
        note = "VANISHES (a odd)" if a % 2 else ""
        print("  %d  %d  %2d   %-7s  %-14s  %-14s  %s"
              % (a, m, P, psi, mp.nstr(abs(lq), 10), mp.nstr(abs(pred), 10), note),
              flush=True)
