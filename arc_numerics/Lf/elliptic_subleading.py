"""The subleading elliptic constant, for P = 1 mod n where the leading term dies.

With eps = tau - w = delta u, j - j(w) = c eps^n (1 + C1 eps + ...),
g = g_a eps^a (1 + G1 eps + ...), K(tau) = K(w) + K'(w) eps + ..., and eta = c delta^n,

  f_eta K dtau = g_a/c^m delta^{-(P-1)} u^a du/(u^n-1)^m
                 [ K(w) + delta u ( K'(w) + G1 K(w) - m C1 u^n K(w)/(u^n-1) ) + O(delta^2) ],

  C1 = j^{(n+1)}(w) / ((n+1) j^{(n)}(w)),   G1 = g^{(a+1)}(w) / ((a+1) g^{(a)}(w)).

When P = 1 mod n the leading bracket cancels (parity at i, equal ray weights at rho) and
the survivor is at delta^{-(P-2)}:

  lim eta^{(P-2)/n} L* = g_a c^{(P-2)/n - m} e^{-i pi s/2}
        sum_w { [K'(w) + G1 K(w)] M(a+1, m) - m C1 K(w) M(a+n+1, m+1) },

  M(b,mu) = (-1)^mu / n * e^{i pi beta} Gamma(beta) Gamma(mu-beta)/Gamma(mu),  beta=(b+1)/n,

each M carrying its own n-th root of unity, and at i (a LINE) an extra factor 1+(-1)^b.
The degeneration does not recur: at i, b = a+1 is even so 1+(-1)^b = 2; at rho the relative
root moves from w^{a+1} = 1 to w^{a+2} != 1.

TESTED AT i FIRST -- there Z is one point, so no relative root is free and the prediction is
fully determined.  Only if that passes is rho used to fix the roots.

K = sum of both kernels, so K(w) = w^{s-1} + ktil(w,s) and
K'(w) = (s-1) w^{s-2} - (1-s) zeta(2-s, w+1) + e^{i pi (s-1)} (1-k+s) zeta(2-k+s, w+1).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, E6, jay, ktil, neville                # noqa: E402

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


def Kval(w, s, k):
    return w**(s - 1) + ktil(w, s, k)


def Kder(w, s, k):
    """d/dtau [ tau^{s-1} + ktil(tau,s) ] at w.

    Taken numerically rather than from -w zeta(w+1,a): that identity is correct but
    evaluates to NaN at w = 0, i.e. at s = k-1, where the Hurwitz zeta has its pole and the
    0 * infinity is only removable in the limit (zeta(0,a) = 1/2 - a, derivative -1)."""
    return mp.diff(lambda t: Kval(t, s, k), w)


def M(b, mu, n):
    beta = mp.mpf(b + 1) / n
    return ((-1)**mu / mp.mpf(n) * mp.e**(I * pi * beta)
            * mp.gamma(beta) * mp.gamma(mu - beta) / mp.gamma(mu))


# (label, g, m, a, k, n, tau_e, j(tau_e))  -- both have P = 1 mod n
CASES = [("i   g=E_6   m=2 a=1 P=3 k=6", lambda t: E6(t), 2, 1, 6, 2, I, mp.mpf(1728)),
         ("rho g=E_4^2 m=2 a=2 P=4 k=8", lambda t: E4(t)**2, 2, 2, 8, 3, RHO, mp.mpf(0))]

print("dps = %d.  Subleading constant at P = 1 mod n.\n" % mp.mp.dps, flush=True)
for lbl, g, m, a, k, n, te, je in CASES:
    P = m * n - a
    kap = mp.mpf(P - 2) / n
    c = mp.diff(jay, te, n) / mp.factorial(n)
    c1 = mp.diff(jay, te, n + 1) / ((n + 1) * mp.diff(jay, te, n))
    ga = mp.diff(g, te, a) / mp.factorial(a) if a else g(te)
    G1 = mp.diff(g, te, a + 1) / ((a + 1) * mp.diff(g, te, a)) if a else \
        mp.diff(g, te, 1) / g(te)
    pref = ga * c**(kap - m)
    for sv in (7, 11):
        s = mp.mpf(sv)
        xs, ys = [], []
        for es in ETAS:
            eta = I * mp.mpf(es)
            f = lambda t, eta=eta, g=g, m=m, je=je: g(t) / (jay(t) - je - eta)**m
            xs.append(mp.mpf(es)**(mp.mpf(1) / n))
            ys.append(eta**kap * Lstar(f, s, k))
        lim, ebar = neville(xs, ys)
        M1, M2 = M(a + 1, m, n), M(a + n + 1, m + 1, n)
        if n == 2:
            par1 = 1 + (-1)**(a + 1)
            par2 = 1 + (-1)**(a + n + 1)
            term = ((Kder(te, s, k) + G1 * Kval(te, s, k)) * par1 * M1
                    - m * c1 * Kval(te, s, k) * par2 * M2)
            preds, tags = [abs(pref * mp.e**(-I * pi * s / 2) * term)], ["line"]
        else:
            z3 = [mp.e**(2 * I * pi * r / 3) for r in range(3)]
            preds, tags = [], []
            for r1, z1 in enumerate(z3):
                for r2, z2 in enumerate(z3):
                    tot = 0
                    for w, sgn in ((RHO, mp.mpf(1)), (TAU0, mp.mpf(-1))):
                        zz1 = z1 if w is TAU0 else mp.mpf(1)
                        zz2 = z2 if w is TAU0 else mp.mpf(1)
                        tot += sgn * ((Kder(w, s, k) + G1 * Kval(w, s, k)) * zz1 * M1
                                      - m * c1 * Kval(w, s, k) * zz2 * M2)
                    preds.append(abs(pref * mp.e**(-I * pi * s / 2) * tot))
                    tags.append("w^%d,w^%d" % (r1, r2))
        b = min(range(len(preds)),
                key=lambda idx: abs(preds[idx] - abs(lim)) / abs(lim))
        print("  %s  s=%d" % (lbl, sv), flush=True)
        print("      extrapolated |lim| = %-26s (shift %s)"
              % (mp.nstr(abs(lim), 20), mp.nstr(ebar, 4)), flush=True)
        print("      predicted    |lim| = %-26s [%s]  rel %s"
              % (mp.nstr(preds[b], 20), tags[b],
                 mp.nstr(abs(preds[b] - abs(lim)) / abs(lim), 6)), flush=True)
