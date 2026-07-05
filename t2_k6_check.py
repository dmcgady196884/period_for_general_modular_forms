#!/usr/bin/env python3
r"""t2_k6_check.py -- T2 first bite: the k = 2 (mod 4), N = 0 edge case (Delta/E6, k=6).

Derived by hand, to be confirmed here:
 [0] ktil(tau, k/2) == 0 identically at k=6   (reflection identity at central s).
 [1] I_0(3, z) == 0 for Im z > 1              (G_{3,6}(-n) == 0 termwise).
 [2] EXACT closed form at integer s=2, all Im z > 1:
       I_0(2,z) = 1/4 + Li2(w)/(2 pi^2) - 3 Li3(w)/(4 pi^3) + 3 Li4(w)/(8 pi^4),
       w = e^{2 pi i (z - i)};   at z = i :  27/80 - 3 zeta(3)/(4 pi^3).
 [3] chamber map at s=3: two-segment L*(Delta/E6, 3) at several tau_0, vs the
     prediction {0, +-2 pi a_{-1} phases} and PV = pi a_{-1} vs known
     B(3) = varpi^4/(576 pi^4).
"""
import mpmath as mp

mp.mp.dps = 30
I, pi = mp.j, mp.pi
FL = dict(flush=True)

# ---------------- Eisenstein machinery ----------------
def sigma(a, m): return sum(d**a for d in range(1, m + 1) if m % d == 0)
NT = 60
_e4 = [240 * sigma(3, m) for m in range(NT, 0, -1)] + [mp.mpf(1)]
_e6 = [-504 * sigma(5, m) for m in range(NT, 0, -1)] + [mp.mpf(1)]
def q(t): return mp.e**(2 * pi * I * t)
def E4(t): return mp.polyval(_e4, q(t))
def E6(t): return mp.polyval(_e6, q(t))
def Delta(t): return (E4(t)**3 - E6(t)**2) / 1728
f = lambda t: Delta(t) / E6(t)          # k = 6, simple pole at i (E6(i) = 0)
k = 6

def ktil(t, s): return mp.zeta(1 - s, t + 1) - mp.e**(I * pi * (s - 1)) * mp.zeta(1 - (k - s), t + 1)

# ---------------- [0] kernel vanishes at central s ----------------
print("[0] ktil(tau, 3) at three points (k=6; predict 0):", **FL)
for t in (mp.mpc('0.3', '1.2'), mp.mpc('-0.4', '1.05'), mp.mpc('0.1', '2.3')):
    print("    ", mp.nstr(abs(ktil(t, 3)), 3), **FL)

# ---------------- residue of f at i ----------------
r_c, Np = mp.mpf('0.05'), 96
a1 = sum(f(I + r_c * mp.e**(I * 2 * pi * m / Np)) * r_c * mp.e**(I * 2 * pi * m / Np)
         for m in range(Np)) / Np
print("\na_{-1} = Res_i(Delta/E6) = %s" % mp.nstr(a1, 20), **FL)

# ---------------- I_0 by direct quadrature ----------------
def Li0(w): return w / (1 - w)
def I0_quad(s, z):
    g = lambda u: Li0(mp.e**(2 * pi * I * ((I - 1 + u) - z))) * ktil(I - 1 + u, s)
    return mp.e**(-I * pi * s / 2) * mp.quad(g, [0, '0.25', '0.5', '0.75', 1])

print("\n[1] I_0(3, z), Im z > 1 (predict 0):", **FL)
for z in (mp.mpc(0, '1.3'), mp.mpc('0.2', '1.15')):
    print("    z=%s :  |I_0| = %s" % (mp.nstr(z, 5), mp.nstr(abs(I0_quad(3, z)), 3)), **FL)

print("\n[2] I_0(2, z) vs exact polylog formula:", **FL)
def I0_exact(z):
    w = mp.e**(2 * pi * I * (z - I))
    return (mp.mpf(1)/4 + mp.polylog(2, w)/(2*pi**2)
            - 3*mp.polylog(3, w)/(4*pi**3) + 3*mp.polylog(4, w)/(8*pi**4))
for z in (mp.mpc(0, '1.3'), mp.mpc(0, '1.1'), mp.mpc('0.15', '1.05')):
    d = I0_quad(2, z) - I0_exact(z)
    print("    z=%s :  quad-exact = %s" % (mp.nstr(z, 5), mp.nstr(abs(d), 3)), **FL)
val_i = mp.mpf(27)/80 - 3*mp.zeta(3)/(4*pi**3)
print("    predicted  I_0(2, i) = 27/80 - 3 zeta(3)/(4 pi^3) = %s" % mp.nstr(val_i, 20), **FL)
print("    formula limit w->1                                = %s" % mp.nstr(I0_exact(I).real, 20), **FL)

# ---------------- [3] chamber map at s = 3 ----------------
print("\n[3] two-segment L*(f,3) at several tau_0  (T-part = 0 since ktil(.,3)=0):", **FL)
def Lstar(s, t0, extra=()):
    a, b = -1/t0, t0
    pts = [a, *extra, b]
    S = mp.mpf(0)
    for p0, p1 in zip(pts, pts[1:]):
        S += mp.quad(lambda u: f(p0 + u*(p1-p0)) * (p0 + u*(p1-p0))**(s-1) * (p1-p0),
                     [0, '0.2', '0.4', '0.6', '0.8', 1])
    T = mp.quad(lambda u: f((t0-1) + u) * ktil((t0-1) + u, s), [0, '0.5', 1])
    return mp.e**(-I*pi*s/2) * (S + T)

cases = [("tau0 = 0.3+0.95i (S-path below i)", mp.mpc('0.3', '0.95'), ()),
         ("tau0 = 0.3+1.4i  (E1 class)",       mp.mpc('0.3', '1.4'),  ()),
         ("tau0 = 0.1+1.2i",                   mp.mpc('0.1', '1.2'),  ()),
         ("tau0 = 0.3+1.4i, S routed ABOVE i", mp.mpc('0.3', '1.4'),
          (mp.mpc('-0.05', '1.15'),))]
vals = []
for nm, t0, wp in cases:
    v = Lstar(3, t0, wp)
    vals.append(v)
    print("    %-36s L* = %s" % (nm, mp.nstr(v, 12)), **FL)

print("\n    reference quantities:", **FL)
print("      2 pi a_{-1}          = %s" % mp.nstr(2*pi*a1, 15), **FL)
print("      pi a_{-1}            = %s" % mp.nstr(pi*a1, 15), **FL)
varpi = mp.gamma(mp.mpf(1)/4)**2 / (2*mp.sqrt(2*pi))
Bknown = varpi**4 / (576 * pi**4)
print("      B(3) known           = %s" % mp.nstr(Bknown, 15), **FL)
print("      B-values (-i L*): %s" % ", ".join(mp.nstr(-I*v, 12) for v in vals), **FL)
