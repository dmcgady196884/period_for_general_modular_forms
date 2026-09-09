"""eq:elllimit's constant, to many digits, by extrapolation.

The loose check in elliptic_constant.py agreed to 9.7e-9 at i and only 6e-4 at rho.  Both
numbers were limited by the MEASUREMENT, not the formula:

  - the compared limits were hardcoded from an 8-significant-figure printout, so 1e-9 was
    already agreement to the last printed digit;
  - the approach to the limit is only O(eta^{1/n}) (delta propto eta^{1/n}), so at
    eta = 1e-5 the residual is ~1e-4 at n = 3.  Brute force would need eta ~ 1e-30.

Fixed here by (i) taking c = j^{(n)}(tau_e)/n! and g_a = g^{(a)}(tau_e)/a! from derivatives
rather than one-sided differences, and (ii) Neville-extrapolating eta^kappa L*(f_eta,s) in
the variable x = eta^{1/n} to x = 0.

Predicted, from eq:ellmodel + eq:raymodel + eta = c delta^n:
    lim = g_a c^{(P-1)/n - m} e^{-i pi s/2} * (-1)^m/n * e^{i pi alpha}
          * Gamma(alpha)Gamma(m-alpha)/Gamma(m) * S,     alpha = (a+1)/n,
with S = (1 + (-1)^a) [w^{s-1} + ktil(w,s)] at w = i (one interior LINE, no free relative
root), and S = [rho-term] - z [tau_0-term] at rho, z one of the three cube roots of unity
(only the RELATIVE root is unconstrained; the overall one is the lemma's stated ambiguity).
Moduli are compared, since the overall n-th root of unity is free by construction.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, Delta, jay, ktil, neville             # noqa: E402

mp.mp.dps = 40
RHO = mp.e**(2 * I * pi / 3)
TAU0 = mp.e**(I * pi / 3)
ETAS = ['1e-2', '1e-3', '1e-4', '1e-5', '1e-6', '1e-7']


def arc_int(g, depth=26):
    a, b = pi / 3, 2 * pi / 3
    mid = (a + b) / 2
    ns = [a, b, mid]
    for base in (a, b, mid):
        d = (b - a) / 4
        for _ in range(depth):
            d /= 2
            ns += [base - d, base + d]
    ns = sorted(set(x for x in ns if a <= x <= b))
    return -sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
                for u, v in zip(ns[:-1], ns[1:]))


def Lstar(f, s, k):
    return mp.e**(-I * pi * s / 2) * (arc_int(lambda t: f(t) * t**(s - 1))
                                      + arc_int(lambda t: f(t) * ktil(t, s, k)))


def Ksum(w, s, k):
    return w**(s - 1) + ktil(w, s, k)


CASES = [
    ("i   a=0 m=1 P=2 k=12", lambda t: E4(t)**3, 1, 0, 12, 2, I, mp.mpf(1728)),
    ("i   a=0 m=2 P=4 k=12", lambda t: E4(t)**3, 2, 0, 12, 2, I, mp.mpf(1728)),
    ("rho a=0 m=1 P=3 k=12", lambda t: Delta(t), 1, 0, 12, 3, RHO, mp.mpf(0)),
    ("rho a=1 m=1 P=2 k=4 ", lambda t: E4(t), 1, 1, 4, 3, RHO, mp.mpf(0)),
    # a = 0 selected w^1 and a = 1 selected w^2, consistent with the relative cube root
    # being w^{a+1}.  Two points do not pin a rule with three possible values, so this
    # third case discriminates: a = 2 predicts w^3 = w^0, unlike either observed value.
    ("rho a=2 m=2 P=4 k=8 ", lambda t: E4(t)**2, 2, 2, 8, 3, RHO, mp.mpf(0)),
]
if len(sys.argv) > 1:
    CASES = [c for c in CASES if c[0].startswith(sys.argv[1])]

print("dps = %d, Neville-extrapolated in x = eta^{1/n}\n" % mp.mp.dps, flush=True)
for lbl, g, m, a, k, n, te, je in CASES:
    P = m * n - a
    kap = mp.mpf(P - 1) / n
    c = mp.diff(jay, te, n) / mp.factorial(n)
    ga = mp.diff(g, te, a) / mp.factorial(a) if a else g(te)
    al = mp.mpf(a + 1) / n
    base = (ga * c**(mp.mpf(P - 1) / n - m) * mp.e**(-I * pi * mp.mpf(3) / 2)
            * (-1)**m / mp.mpf(n) * mp.e**(I * pi * al)
            * mp.gamma(al) * mp.gamma(m - al) / mp.gamma(m))
    for sv in (3, 7):
        s = mp.mpf(sv)
        base_s = (ga * c**(mp.mpf(P - 1) / n - m) * mp.e**(-I * pi * s / 2)
                  * (-1)**m / mp.mpf(n) * mp.e**(I * pi * al)
                  * mp.gamma(al) * mp.gamma(m - al) / mp.gamma(m))
        xs, ys = [], []
        for es in ETAS:
            eta = I * mp.mpf(es)
            f = lambda t, eta=eta, g=g, m=m, je=je: g(t) / (jay(t) - je - eta)**m
            xs.append(mp.mpf(es)**(mp.mpf(1) / n))
            ys.append(eta**kap * Lstar(f, s, k))
        lim, ebar = neville(xs, ys)
        if n == 2:
            preds = [abs(base_s * (1 + (-1)**a) * Ksum(te, s, k))]
            tags = ["line"]
        else:
            z3 = [mp.e**(2 * I * pi * r / 3) for r in range(3)]
            preds = [abs(base_s * (Ksum(RHO, s, k) - z * Ksum(TAU0, s, k))) for z in z3]
            tags = ["rel=w^%d" % r for r in range(3)]
        b = min(range(len(preds)),
                key=lambda idx: abs(preds[idx] - abs(lim)) / abs(lim))
        print("  %s s=%d" % (lbl, sv), flush=True)
        print("      extrapolated |lim| = %-26s (last-order shift %s)"
              % (mp.nstr(abs(lim), 22), mp.nstr(ebar, 4)), flush=True)
        print("      predicted    |lim| = %-26s [%s]  rel %s"
              % (mp.nstr(preds[b], 22), tags[b],
                 mp.nstr(abs(preds[b] - abs(lim)) / abs(lim), 6)), flush=True)
