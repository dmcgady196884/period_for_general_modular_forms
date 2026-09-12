"""Does the connector ambiguity at rho reach the PERIOD POLYNOMIAL, or only L*?

basepoint_shift_rho.py: the two connector routings give different (individually flat) L*, and
short - long = -2 pi i e^{-i pi s/2} Res_rho(f ktil).  But eq:arcondefect builds the
(1+U+U^2) defect out of a_p = Res_{tau=p} f -- the residue of f ITSELF, not of f ktil -- and the
same a_p controls the T-winding shift of r_f.  So both "does W select a routing" and "does the
ambiguity reach hat r_f" turn on Res_rho f.

CLAIM.  At an elliptic point of order n the stabiliser V (V'(rho) = rho = omega) with f|_k V = f
forces, for f = sum_j a_j u^{-P+j} in u = tau - rho, the relation a_j omega^{-P+j} =
(c rho + d)^k a_j; taking j = 0 (a_0 != 0) fixes the character, leaving omega^j = 1.  So the
Laurent support is j = 0 mod 3, the residue is a_{P-1}, and

    Res_rho f != 0   <=>   P = 1 mod 3.

PREDICTION.  E_4/j (P=2) and Delta/j (P=3): ZERO.  E_4^2/j^2 (P=4): NONZERO.
Both read off circles of three radii, so a genuine residue is radius-independent and a
quadrature artifact is not.  Res_rho(f ktil) is printed alongside as the control: it is nonzero
in every case (it is what short - long measured), so a zero in the f column is a statement about
f, not about the contour or the quadrature.

If the prediction holds, the W-membership test is informative only at P = 1 mod 3, and elsewhere
the ambiguity plausibly does not reach hat r_f at all.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, Delta, jay, ktil                      # noqa: E402

mp.mp.dps = 30
RHO = mp.e**(2 * I * pi / 3)


def residue(g, c, r, n=60):
    ps = [2 * pi * mp.mpf(j) / n for j in range(n + 1)]
    tot = sum(mp.quad(lambda p: g(c + r * mp.e**(I * p)) * I * r * mp.e**(I * p), [u, v])
              for u, v in zip(ps[:-1], ps[1:]))
    return tot / (2 * I * pi)


CASES = [
    ("P=2  E_4/j      k=4   (P = 2 mod 3 -> expect 0)", lambda t: E4(t) / jay(t), 4),
    ("P=3  Delta/j    k=12  (P = 0 mod 3 -> expect 0)", lambda t: Delta(t) / jay(t), 12),
    ("P=4  E_4^2/j^2  k=8   (P = 1 mod 3 -> expect nonzero)",
     lambda t: E4(t)**2 / jay(t)**2, 8),
]

for lbl, f, k in CASES:
    print("=" * 76, flush=True)
    print(lbl, flush=True)
    for r in (mp.mpf('0.06'), mp.mpf('0.03')):
        a = residue(f, RHO, r)
        b = residue(lambda t: f(t) * ktil(t, mp.mpf(7), k), RHO, r)
        print("   r=%-6s Res_rho f       = %s" % (mp.nstr(r, 3), mp.nstr(a, 12)), flush=True)
        print("           Res_rho(f ktil) = %s   (control, s=7)" % mp.nstr(b, 12), flush=True)
