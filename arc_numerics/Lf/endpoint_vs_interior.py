"""A pole AT rho is not the same case as a pole at i, and the difference is the endpoints.

detour_vs_retract.py settled i: the pole is INTERIOR to gamma^arc, both detours are finite and
r-independent, they differ by 2 pi i Res, and the mean is real.  rho is different.  def:georef
pins tau_0 = rho+1 because S(rho+1) = rho = (rho+1)-1, and lem:arcon shows W-membership forces
that, so tau_0 cannot be slid off.  A pole at SL2(Z).rho is therefore a pole at BOTH ENDPOINTS
of int_{tau_0}^{S tau_0}, and a path cannot go around a point where it terminates.

Cutting the arc at angular distance alpha from each end, the endpoint expansion of f K is
sum_j A_j (tau - rho)^{-P+j} at rho and sum_j B_j (tau - tau_0)^{-P+j} at tau_0, with the same
f-Laurent data (T-periodicity) but different K-data.  Traversal directions differ by the
stabiliser rotation omega = e^{2 pi i/3} and orientations are opposite, so order j contributes
    alpha^{1-P+j} [A_j - omega^{1-P+j} B_j]        (j <= P-2)
    log alpha     [Res_rho(fK) - Res_{tau_0}(fK)]  (j = P-1).
The leading term cancels iff P = 1 mod 3 and K(rho) = K(tau_0); by lem:ellLaurent (2P = k mod 6)
the first forces k = 2 mod 6, which by eq:P1cancel is the second.  Nothing kills the rest.

PREDICTIONS.
 (1) rho, f = E_4/j, k = 4, P = 3m - a = 2:  cutting both ends by alpha gives a 1/alpha
     divergence AND a surviving log alpha.  Both coefficients grid-stable.
 (2) i, f = E_4^3/(j-1728), k = 12, P = 2:  the SYMMETRIC cut about pi/2 gives 1/alpha but NO
     log -- the symmetric cut kills the odd (residue) term, int_{-R}^{-eps} + int_{eps}^{R} of
     du/u being zero.
 (3) i: the finite part of that symmetric cut EQUALS the detour mean already measured,
     0.00086752559522 at s = 7 (and the s = 11 value), to the fit's own precision.  This is what
     ties "Hadamard finite part" to "mean of the two detours" -- for P >= 2 they agree while the
     naive symmetric-cut principal value diverges, so the paper must say detour-mean, not PV.

Log columns are judged by the fit_sanity.py discriminant: c_log(grid/2)/c_log(grid) -> 1 for a
genuine log, -> 2^-(first omitted power) when it is only absorbing truncation.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, E6, Delta, jay, ktil                  # noqa: E402

mp.mp.dps = 30
EXPS = [-1, 0, 1, 2]
DETOUR_MEAN = {7: mp.mpf('0.00086752559521718496505'),      # (in + out)/2 from detour_vs_retract
               11: mp.mpf('0.0026014584280577241819')}


def seg(g, th0, th1, n=20):
    xs = [th0 + (th1 - th0) * mp.mpf(j) / n for j in range(n + 1)]
    return sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
               for u, v in zip(xs[:-1], xs[1:]))


def Lcut(gf, m, x, s, k, pieces):
    f = lambda t: gf(t) / (jay(t) - x)**m
    tot = 0
    for th0, th1 in pieces:
        tot += seg(lambda t: f(t) * (t**(s - 1) + ktil(t, s, k)), th0, th1)
    return mp.e**(-I * pi * s / 2) * tot


def rho_pieces(al):
    """Arc rho -> rho+1 (theta decreasing), both ENDS cut back by al."""
    return [(2 * pi / 3 - al, pi / 3 + al)]


def i_pieces(al):
    """Same arc, symmetric cut of half-width al about theta = pi/2."""
    return [(2 * pi / 3, pi / 2 + al), (pi / 2 - al, pi / 3)]


def solve(ds, vals, withlog):
    cols = list(EXPS) + (['log'] if withlog else [])
    A = mp.matrix(len(cols), len(cols))
    b = mp.matrix(len(cols), 1)
    for i in range(len(cols)):
        for jx, c in enumerate(cols):
            A[i, jx] = mp.log(ds[i]) if c == 'log' else ds[i]**c
        b[i] = vals[i]
    x = mp.lu_solve(A, b)
    return {c: x[jx] for jx, c in enumerate(cols)}


CASES = [
    ("rho  f = E_4/j,          k=4,  P=2, pole at BOTH ENDPOINTS",
     lambda t: E4(t), 1, mp.mpf(0), 4, rho_pieces, False),
    ("i    f = E_4^3/(j-1728), k=12, P=2, pole INTERIOR",
     lambda t: E4(t)**3, 1, mp.mpf(1728), 12, i_pieces, True),
]

for lbl, gf, m, x, k, pieces, cmp_detour in CASES:
    print("=" * 78, flush=True)
    print(lbl, flush=True)
    for s in (mp.mpf(7), mp.mpf(11)):
        need = len(EXPS) + 1
        g1 = [mp.mpf(2)**(-e) for e in range(5, 5 + need)]
        g2 = [mp.mpf(2)**(-e) for e in range(6, 6 + need)]
        out = []
        for ds in (g1, g2):
            v = [Lcut(gf, m, x, s, k, pieces(d)) for d in ds]
            out.append((solve(ds, v, False), solve(ds, v, True)))
        (n1, w1), (n2, w2) = out
        rat = w2['log'] / w1['log'] if w1['log'] != 0 else mp.mpf(0)
        print("  s = %s" % mp.nstr(s, 3), flush=True)
        print("     c_-1        = %s" % mp.nstr(w1[-1], 14), flush=True)
        print("     c_log       = %-24s -> %-24s ratio %s"
              % (mp.nstr(w1['log'], 10), mp.nstr(w2['log'], 10), mp.nstr(rat, 8)), flush=True)
        print("     c_0 (w/log) = %-26s grid2 = %s"
              % (mp.nstr(w1[0], 16), mp.nstr(w2[0], 16)), flush=True)
        print("     c_0 (no log)= %-26s grid2 = %s"
              % (mp.nstr(n1[0], 16), mp.nstr(n2[0], 16)), flush=True)
        if cmp_detour:
            d = DETOUR_MEAN[int(s)]
            print("     detour mean = %-26s  |c_0 - mean| = %s"
                  % (mp.nstr(d, 16), mp.nstr(abs(w1[0] - d), 6)), flush=True)
