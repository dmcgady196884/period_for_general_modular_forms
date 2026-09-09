"""Pin the constant in eq:elllimit, and probe the two P = 1 cases.

PART A -- the constant.  Collecting eq:ellmodel, eq:raymodel and eta = c delta^n,

    lim eta^{(P-1)/n} L*(f_eta,s)
        = g_a c^{(P-1)/n - m} e^{-i pi s/2} sum_{w in Z} sig_w [ w^{s-1} + ktil(w,s) ]
          * int_{C_w} u^a du/(u^n-1)^m,

with int_C = (-1)^m/n e^{i pi alpha} Gamma(alpha)Gamma(m-alpha)/Gamma(m) per ray,
alpha = (a+1)/n, each ray carrying its own n-th root of unity, and a full line at w = i
carrying the factor (1 + (-1)^a) instead.

At tau_e = i, Z is the single interior point, so NO relative root of unity is free and the
modulus is fully determined -- that is a clean test.  At tau_e = rho, Z = {rho, tau_0} and
only the RELATIVE cube root between the two endpoints is unknown; we scan the three and
report which reproduces the measured modulus, thereby pinning it.

Measured limits (moduli) from ell_exponent_scan.py, eta = 1e-5 row:
    i,   a=0,m=1,P=2,k=12  s=3: 0.068342846   s=7: 0.043930493(1e-4 row)
    i,   a=0,m=2,P=4,k=12  s=3: 0.034172932   s=7: 0.021968313
    rho, a=0,m=1,P=3,k=12  s=3: 0.0002273218  s=7: 0.0001461355
    rho, a=1,m=1,P=2,k=4   s=3: 0.01140417    s=7: 0.006946569(1e-4 row)

PART B -- P = 1 at rho (m=1, a=2; g = E_4^2, so f_0 = E_4^2/j = Delta/E_4, k = 8).
eq:raymodel needs P > 1, so this case is outside the derivation.  kappa = (P-1)/n = 0
would mean L*(f_eta,s) simply converges.  But the contour passes through a SIMPLE pole at
both endpoints and the two log divergences cancel only if the kernels agree there, which
they do not: rho^{s-1} != tau_0^{s-1}, and ktil(rho,s) != ktil(tau_0,s) since
zeta(w,rho+2) = zeta(w,rho+1) - (rho+1)^{-w}.  So the honest question is whether
L*(f_eta,s) converges or grows like log eta.  A constant increment per DECADE of eta means
log divergence -- and then no power eta^kappa can regulate it, i.e. def:arcsplit's ansatz
is inadequate at P = 1 on the arc.

PART C -- P = 1 at i with a odd (m=1, a=1; g = E_6, k = 6).  Here eq:ellkappa gives
kappa = (P-2)/n = -1/2 < 0, which predicts L*(f_eta,s) -> 0 like eta^{1/2}: the "simple
pole at i is invisible" dichotomy.  Checked by watching |L*| and |eta^{-1/2} L*|.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, E6, Delta, jay, ktil                  # noqa: E402

mp.mp.dps = 30
RHO = mp.e**(2 * I * pi / 3)
TAU0 = mp.e**(I * pi / 3)


def arc_int(g, depth=22):
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


def local_c(te, n, je):
    """c = lim (j - j(tau_e)) / (tau - tau_e)^n"""
    h = mp.mpf('1e-6')
    return (jay(te + h) - je) / h**n


def local_ga(g, te, a):
    """g_a = lim (tau - tau_e)^{-a} g(tau)"""
    h = mp.mpf('1e-6')
    return g(te + h) / h**a


def Ksum(w, s, k):
    return w**(s - 1) + ktil(w, s, k)


def ray_model(n, a, m):
    al = mp.mpf(a + 1) / n
    return (-1)**m / mp.mpf(n) * mp.e**(I * pi * al) * mp.beta(al, m - al)


print("dps = %d\n" % mp.mp.dps, flush=True)
print("=" * 78, flush=True)
print("PART A: eq:elllimit constant vs the measured limits (moduli)", flush=True)
print("=" * 78, flush=True)

# (label, g, m, a, k, n, tau_e, j(tau_e), {s: measured |limit|})
ACASES = [
    ("i   a=0 m=1 P=2 k=12", lambda t: E4(t)**3, 1, 0, 12, 2, I, mp.mpf(1728),
     {3: '0.068342846', 7: '0.043930493'}),
    ("i   a=0 m=2 P=4 k=12", lambda t: E4(t)**3, 2, 0, 12, 2, I, mp.mpf(1728),
     {3: '0.034172932', 7: '0.021968313'}),
    ("rho a=0 m=1 P=3 k=12", lambda t: Delta(t), 1, 0, 12, 3, RHO, mp.mpf(0),
     {3: '0.0002273218', 7: '0.0001461355'}),
    ("rho a=1 m=1 P=2 k=4 ", lambda t: E4(t), 1, 1, 4, 3, RHO, mp.mpf(0),
     {3: '0.01140417', 7: '0.006946569'}),
]

for lbl, g, m, a, k, n, te, je, meas in ACASES:
    P = m * n - a
    c = local_c(te, n, je)
    ga = local_ga(g, te, a)
    pref_pow = mp.mpf(P - 1) / n - m
    for sv, mv in meas.items():
        s = mp.mpf(sv)
        base = ga * c**pref_pow * mp.e**(-I * pi * s / 2) * ray_model(n, a, m)
        if n == 2:
            pred = [abs(base * (1 + (-1)**a) * Ksum(te, s, k))]
            tag = ["line"]
        else:
            # two rays; only the RELATIVE cube root is unknown, so scan it
            z3 = [mp.e**(2 * I * pi * r / 3) for r in range(3)]
            pred = [abs(base * (Ksum(RHO, s, k) - z * Ksum(TAU0, s, k))) for z in z3]
            tag = ["rel=w^%d" % r for r in range(3)]
        best = min(range(len(pred)),
                   key=lambda idx: abs(pred[idx] - mp.mpf(mv)) / mp.mpf(mv))
        rel = abs(pred[best] - mp.mpf(mv)) / mp.mpf(mv)
        print("  %s s=%d  measured %-14s predicted %-14s  [%s]  rel %s"
              % (lbl, sv, mv, mp.nstr(pred[best], 8), tag[best], mp.nstr(rel, 5)),
              flush=True)

print("\n" + "=" * 78, flush=True)
print("PART B: P = 1 at rho.  f = E_4^2/j = Delta/E_4, k = 8.  Converge, or log eta?",
      flush=True)
print("=" * 78, flush=True)
for s in (mp.mpf(3), mp.mpf(7)):
    print("  s = %s" % mp.nstr(s, 4), flush=True)
    prev = None
    for es in ('1e-2', '1e-3', '1e-4', '1e-5', '1e-6'):
        eta = I * mp.mpf(es)
        f = lambda t, eta=eta: E4(t)**2 / (jay(t) - eta)
        L = Lstar(f, s, 8)
        d = "" if prev is None else "   change %s" % mp.nstr(abs(L - prev), 8)
        print("    eta=%-6s |L*| = %-18s%s" % (es, mp.nstr(abs(L), 12), d), flush=True)
        prev = L

print("\n" + "=" * 78, flush=True)
print("PART C: P = 1 at i, a odd.  f = E_6/(j-1728-eta), k = 6.  Does L* -> 0?",
      flush=True)
print("=" * 78, flush=True)
for s in (mp.mpf(7), mp.mpf(11)):
    print("  s = %s" % mp.nstr(s, 4), flush=True)
    for es in ('1e-2', '1e-3', '1e-4', '1e-5'):
        eta = I * mp.mpf(es)
        f = lambda t, eta=eta: E6(t) / (jay(t) - 1728 - eta)
        L = Lstar(f, s, 6)
        print("    eta=%-6s |L*| = %-18s |eta^{-1/2} L*| = %s"
              % (es, mp.nstr(abs(L), 12), mp.nstr(abs(L / mp.sqrt(eta)), 10)), flush=True)
