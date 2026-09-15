"""Project the Laurent coefficients c_m onto the def:Wpm triple {W_-, W_+, p_0} at k = 18.

common.py records that at k=12 the three basis vectors of def:Wpm have DISJOINT supports --
W_- on odd l, W_+ on l = 2,4,6,8 (the c p_0 subtraction clears l = 0,10), p_0 on l = 0,10 -- so a
decomposition in this basis is a read-off rather than a solve.  PART 1 rebuilds def:Wpm from the
Bernoulli formula at k = 12 (checked against common.W12_minus/W12_plus/P0_12) and at k = 18.

PREDICTION, from zeta17_check.py.  r_{E_18}/(2 pi i)^{n+1} has rational ODD coordinates, VANISHING
even coordinates l = 2..14, and support at l = 0, 16 only in the even extremes.  In the def:Wpm
basis that reads

    r_{E_18} = (rational) W_-^{(18)}  +  0 * W_+^{(18)}  +  x p_0 ,

with x = 0.1848002737722058... the constant z17_highprec.py could not identify.  Since
Theta(z) = (2 pi i)^{n+2} B(z) - 2 pi i (1-z^n) r_E with B(z) Laurent-RATIONAL, it follows that

    the W_- and W_+ coordinates of every c_m are rational multiples of (2 pi i)^{n+2},
    and x can enter ONLY the p_0 coordinate, and only at m = 0 and m = n.

This run uses the TRUE Eisenstein reference E_18 = E_4^3E_6 + c Delta E_6, not E_4^3E_6.  That
matters: theta_k18_DeltaE4.py used E_4^3E_6, whose odd coordinates are irrational
(zeta17_check.py PART 3) because it differs from E_18 by a cusp form whose periods are
transcendental.  So lambda_0/lambda_2 was irrational there for a reason that has nothing to do with
Theta; with E_18 it is predicted to be rational.

Cost: the 19-point Laurent fit plus one rvec at k = 18.  No result here is a read-off of a previous
run -- the c_m are recomputed against the new reference.
"""
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, E4, E6, Delta, ktil, slash, rvec, inW,      # noqa: E402
                    W12_minus, W12_plus, P0_12, arcint)

mp.mp.dps = 45
Sm = (0, -1, 1, 0)


def bern(i):
    p, q = mp.bernfrac(i)
    return F(int(p), int(q))


def B0(m):
    """B^0_m(X) = sum_{i != 1} C(m,i) B_i X^{m-i}, as {degree: Fraction}."""
    out = {}
    for i in range(m + 1):
        if i == 1:
            continue
        b = bern(i)
        if b:
            out[m - i] = out.get(m - i, F(0)) + F(int(mp.binomial(m, i))) * b
    return out


def rev(poly, n):
    """X^n P(1/X)."""
    return {n - d: c for d, c in poly.items()}


def lin(*terms):
    out = {}
    for co, p in terms:
        for d, c in p.items():
            out[d] = out.get(d, F(0)) + co * c
    return {d: c for d, c in out.items() if c}


def pp(poly, n):
    """def:Wpm's pp, returned as an l-indexed integer vector (l = coeff of X^{n-l}Y^l)."""
    from math import gcd, lcm
    L = 1
    for c in poly.values():
        L = lcm(L, c.denominator)
    ints = {d: int(c * L) for d, c in poly.items()}
    g = 0
    for v in ints.values():
        g = gcd(g, abs(v))
    if g:
        ints = {d: v // g for d, v in ints.items()}
    return [ints.get(n - l, 0) for l in range(n + 1)]


def Wpm(k):
    n = k - 2
    p0 = [0] * (n + 1)
    p0[0], p0[n] = 1, -1
    Wm = pp(lin((F(1, n - 1), B0(n - 1)), (F(1, n - 1), rev(B0(n - 1), n)),
               (F(-1, 3), B0(3)), (F(-1, 3), rev(B0(3), n))), n)
    braw = lin((F(1, 2), B0(2)), (F(-1, 2), rev(B0(2), n)),
               (F(1, n), B0(n)), (F(-1, n), rev(B0(n), n)))
    c = braw.get(n, F(0))                      # the X^n coefficient, def:Wpm's c
    Wp = pp(lin((F(1), braw), (-c, {n: F(1), 0: F(-1)})), n)
    return Wm, Wp, p0


print("PART 1  rebuild def:Wpm and check k=12 against common.py", flush=True)
Wm12, Wp12, p012 = Wpm(12)
for nm, got, want in (("W_-", Wm12, W12_minus), ("W_+", Wp12, W12_plus), ("p_0", p012, P0_12)):
    s = "match" if (got == list(want) or got == [-x for x in want]) else "MISMATCH"
    print("   %-4s %-46s %s  (common.py: %s)" % (nm, got, s, list(want)), flush=True)

Wm, Wp, p0 = Wpm(18)
K, N = 18, 16
print("\nPART 2  def:Wpm at k=18", flush=True)
for nm, v in (("W_-", Wm), ("W_+", Wp), ("p_0", p0)):
    vv = [mp.mpf(x) for x in v]
    s_, u_ = inW(vv, K)
    sup = [l for l, x in enumerate(v) if x]
    print("   %-4s in W: %-10s %-10s  support l = %s" % (nm, mp.nstr(s_, 3), mp.nstr(u_, 3), sup),
          flush=True)
print("   supports disjoint? %s"
      % ("yes" if not (set(l for l, x in enumerate(Wm) if x)
                       & set(l for l, x in enumerate(Wp) if x)
                       | set(l for l, x in enumerate(Wp) if x)
                       & set(l for l, x in enumerate(p0) if x)) else "NO"), flush=True)


def proj(v):
    """read-off coordinates in {W_-, W_+, p_0}, using the disjoint supports."""
    out = []
    for b in (Wm, Wp, p0):
        l = max(range(N + 1), key=lambda j: abs(b[j]))
        out.append(v[l] / b[l])
    return out


cE = -mp.mpf(216) - mp.mpf(28728) / 43867
E18 = lambda t: E4(t)**3 * E6(t) + cE * Delta(t) * E6(t)
rE = rvec(E18, K, 11)
print("\nPART 3  r_{E_18} in the def:Wpm basis (predicted: rational, 0, x)", flush=True)
nrm = (2 * pi * I)**(N + 1)
for nm, val in zip(("W_-", "W_+", "p_0"), proj([x / nrm for x in rE])):
    q = mp.pslq([mp.re(val), mp.mpf(1)], tol=mp.mpf('1e-30'), maxcoeff=10**15, maxsteps=10**6) \
        if abs(mp.re(val)) > mp.mpf('1e-30') else None
    print("   %-4s %-40s rational? %s" % (nm, mp.nstr(val, 22),
                                          ("%s/%s" % (-q[1], q[0])) if q else
                                          ("(zero)" if abs(val) < mp.mpf('1e-28') else "no")),
          flush=True)


def Theta(z):
    P = [mp.binomial(N, l) * (-z)**l for l in range(N + 1)]
    PS = slash(P, N, Sm)
    out = []
    for l in range(N + 1):
        br = (P[l] - PS[l]) + mp.binomial(N, l) * (-1)**l * (
            ktil(z, mp.mpf(l + 1), K) - z**N * ktil(-1 / z, mp.mpf(l + 1), K))
        out.append((2 * pi * I)**(N + 1) * br)
    return [2 * pi * I * (a - (1 - z**N) * b) for a, b in zip(out, rE)]


MS = list(range(-1, N + 2))
FIT = [mp.mpf(a) + I * mp.mpf(b) for a, b in
       [('0.31', '1.27'), ('-0.4', '1.1'), ('0.7', '1.9'), ('0.0', '1.3'), ('0.5', '1.4'),
        ('-0.25', '1.6'), ('0.15', '2.1'), ('0.9', '1.05'), ('-0.6', '1.35'), ('0.05', '1.15'),
        ('0.44', '1.72'), ('-0.8', '2.4'), ('0.33', '1.05'), ('0.22', '2.8'), ('-0.15', '1.45'),
        ('0.66', '1.25'), ('-0.55', '1.95'), ('0.12', '1.62'), ('0.85', '2.2')]]
V = mp.matrix(len(MS), len(MS))
for a, z in enumerate(FIT):
    for b, m in enumerate(MS):
        V[a, b] = z**m
TH = [Theta(z) for z in FIT]
cols = [mp.lu_solve(V, mp.matrix([TH[a][l] for a in range(len(MS))])) for l in range(N + 1)]
cm = [[cols[l][b] for l in range(N + 1)] for b in range(len(MS))]

z = mp.mpf('0.62') + I * mp.mpf('1.11')
pred = [sum(cm[b][l] * z**MS[b] for b in range(len(MS))) for l in range(N + 1)]
print("\n   out-of-sample check: %s"
      % mp.nstr(max(abs(a - b) for a, b in zip(Theta(z), pred))
                / max(abs(x) for x in Theta(z)), 6), flush=True)

print("\nPART 4  c_m in the def:Wpm basis, divided by (2 pi i)^{n+2}", flush=True)
nrm2 = (2 * pi * I)**(N + 2)
print("   %-5s %-26s %-26s %s" % ("m", "W_- coord", "W_+ coord", "p_0 coord"), flush=True)
for b, m in enumerate(MS):
    vals = proj([x / nrm2 for x in cm[b]])
    cells = []
    for val in vals:
        if abs(val) < mp.mpf('1e-26'):
            cells.append("0")
            continue
        r = mp.re(val)
        q = mp.pslq([r, mp.mpf(1)], tol=mp.mpf('1e-26'), maxcoeff=10**15, maxsteps=10**6)
        cells.append(("%s/%s" % (-q[1], q[0])) if q else ("irr " + mp.nstr(r, 12)))
    print("   %-5d %-26s %-26s %s" % (m, cells[0], cells[1], cells[2]), flush=True)
