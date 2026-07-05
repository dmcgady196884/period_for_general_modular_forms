#!/usr/bin/env python3
r"""t2_generic_s.py -- T2 at GENERIC complex s (k=6, N=0):
    I_0(s,z) = 1/s + i^k/(k-s) - [ Pi_0(s;w) + Lambda_0(s;Q) + E_0(s,z) ]
with w = e^{2 pi i (z-i)}, Q = e^{2 pi i z},
  Pi_0     = sum_{j=0,1} (-1)^{j+1} [(s-1)_(j) + i^k (k-s-1)_(j)]/(2 pi)^{j+1} Li_{j+1}(w)
  Lambda_0 = G(s) (-2pi)^{-s} Li_s(Q) + i^k G(k-s) (-2pi)^{-(k-s)} Li_{k-s}(Q)
  E_0      = - sum_n Q^n [ (s-1)_(2) g(s-2,-2pi n)/(-2pi n)^s
                          + i^k (k-s-1)_(2) g(k-s-2,-2pi n)/(-2pi n)^{k-s} ]
(g = lower incomplete gamma).  Branch: principal powers throughout; sanity-check the
Gamma(s,x) recursion + Gamma/gamma split against mpmath first.
"""
import mpmath as mp

mp.mp.dps = 30
I, pi = mp.j, mp.pi
FL = dict(flush=True)
k = 6

def Gup(sig, x):   return mp.gammainc(sig, x, mp.inf)   # upper incomplete
def glow(sig, x):  return mp.gammainc(sig, 0, x)        # lower incomplete

# ---------------- branch sanity ----------------
s0, x0 = mp.mpc('2.6', '0.7'), mp.mpf(-2) * pi
r1 = Gup(s0, x0) - (mp.e**(-x0) * x0**(s0 - 1) + (s0 - 1) * Gup(s0 - 1, x0))
r2 = Gup(s0, x0) - (mp.gamma(s0) - glow(s0, x0))
print("[sanity] recursion residual = %s ;  Gamma/gamma split residual = %s"
      % (mp.nstr(abs(r1), 3), mp.nstr(abs(r2), 3)), **FL)

# ---------------- kernel + quadrature I_0 ----------------
def ktil(t, s): return mp.zeta(1 - s, t + 1) - mp.e**(I * pi * (s - 1)) * mp.zeta(1 - (k - s), t + 1)
def Li0(u): return u / (1 - u)
def I0_quad(s, z):
    g = lambda u: Li0(mp.e**(2 * pi * I * ((I - 1 + u) - z))) * ktil(I - 1 + u, s)
    return mp.e**(-I * pi * s / 2) * mp.quad(g, [0, '0.25', '0.5', '0.75', 1])

# ---------------- three-layer formula ----------------
def poch(a, j):
    out = mp.mpc(1)
    for r in range(j): out *= (a - r)
    return out

def I0_formula(s, z, nmax=400):
    ik = I**k
    w = mp.e**(2 * pi * I * (z - I))
    Q = mp.e**(2 * pi * I * z)
    Pi = mp.mpc(0)
    for j in (0, 1):
        coef = (-1)**(j + 1) * (poch(s - 1, j) + ik * poch(k - s - 1, j)) / (2 * pi)**(j + 1)
        if abs(coef) > mp.mpf(10)**(-40):
            Pi += coef * mp.polylog(j + 1, w)
    m2pi = mp.mpc(-2 * pi)
    Lam = (mp.gamma(s) * m2pi**(-s) * mp.polylog(s, Q)
           + ik * mp.gamma(k - s) * m2pi**(-(k - s)) * mp.polylog(k - s, Q))
    c1, c2 = poch(s - 1, 2), poch(k - s - 1, 2)
    E = mp.mpc(0)
    for n in range(1, nmax + 1):
        xn = mp.mpf(-2 * pi) * n
        t = c1 * glow(s - 2, xn) / xn**s + ik * c2 * glow(k - s - 2, xn) / xn**(k - s)
        E -= Q**n * t
        if abs(Q**n * t) < mp.mpf(10)**(-32) and n > 20:
            break
    return 1/s + ik/(k - s) - (Pi + Lam + E)

# ---------------- compare ----------------
print("\n[compare] k=6, N=0, generic s:", **FL)
for s, z in [(mp.mpc('2.6', '0.7'), mp.mpc(0, '1.15')),
             (mp.mpc('2.6', '0.7'), mp.mpc('0.23', '1.02')),
             (mp.mpc('0.8', '-1.1'), mp.mpc('0.1', '1.1')),
             (mp.mpc(3, 0),          mp.mpc(0, '1.15'))]:
    a = I0_quad(s, z)
    b = I0_formula(s, z)
    print("   s=%-12s z=%-12s |quad - formula| = %s"
          % (mp.nstr(s, 6), mp.nstr(z, 6), mp.nstr(abs(a - b), 3)), **FL)

# limit at z=i for a generic s (k = 2 mod 4: continuous by parity kill)
sgen = mp.mpc('2.6', '0.7')
print("\n[edge] I_0(s, i) at s=2.6+0.7i via formula (w=1):", **FL)
val = I0_formula(sgen, I, nmax=4000)
print("   I_0(s,i) = %s" % mp.nstr(val, 20), **FL)
for d in ('0.05', '0.02', '0.01'):
    zq = I * (1 + mp.mpf(d))
    print("   quad at z=i(1+%s): diff from edge value = %s"
          % (d, mp.nstr(abs(I0_quad(sgen, zq) - val), 3)), **FL)
