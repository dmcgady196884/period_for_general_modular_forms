#!/usr/bin/env python3
r"""t2_k4_pv.py -- execute the 4|k ledger at k=4, s=2 (central), option-c geometry.

Claim (derived by hand):
  lim_{delta->0} L*(P_i f, 2) = a_{-2} (13 pi/6 - 1),
with P_i f = r1 Li_0(e^{2 pi i (tau - i)}) + r2 Li_{-1}(...),
r1 = -2 pi i a_{-1} = 2 pi a_{-2} (parity a_{-1} = i a_{-2}), r2 = -4 pi^2 a_{-2}.
Contour: tau_0 = i(1+delta); gamma^S = imaginary axis from i/(1+delta) to i(1+delta),
pole at i handled by (Hadamard FP for u^{-2}) + (PV for u^{-1}); gamma^T = [tau_0-1, tau_0].
Also verifies the closed-form T-side: T = r1 T_0 + r2 T_1,
T_N = (1+d)/pi Li_{1-N}(v) + 1/(2 pi^2) Li_{2-N}(v), v = e^{-2 pi delta}.
"""
import mpmath as mp

mp.mp.dps = 30
I, pi = mp.j, mp.pi
FL = dict(flush=True)
k, s = 4, 2

am2 = mp.mpf(-1) / (1728 * pi**2)          # a_{-2} of E4 Delta/E6^2 (paper value)
am1 = I * am2                               # parity
r1 = -2 * pi * I * am1                      # = 2 pi a_{-2}
r2 = -4 * pi**2 * am2

def Li0(u):  return u / (1 - u)
def Lim1(u): return u / (1 - u)**2
def Pf(t):
    e = mp.e**(2 * pi * I * (t - I))
    return r1 * Li0(e) + r2 * Lim1(e)
def sing(t):
    return am2 * (t - I)**(-2) + am1 * (t - I)**(-1)
def Pf_minus_sing(t):
    u = 2 * pi * I * (t - I)
    if abs(u) < mp.mpf('0.05'):     # Bernoulli series of the regular part
        li0reg  = -mp.mpf(1)/2 - u/12 + u**3/720 - u**5/30240
        lim1reg = -mp.mpf(1)/12 + u**2/240 - u**4/6048
        return r1 * li0reg + r2 * lim1reg
    return Pf(t) - sing(t)
def ktil(t, ss): return mp.zeta(1 - ss, t + 1) - mp.e**(I * pi * (ss - 1)) * mp.zeta(1 - (k - ss), t + 1)

target = am2 * (13 * pi / 6 - 1)
print("target  a_{-2}(13 pi/6 - 1) = %s" % mp.nstr(target, 15), **FL)

for d in (mp.mpf('0.1'), mp.mpf('0.05'), mp.mpf('0.025')):
    t0 = I * (1 + d)
    # ---- T-part: quad of Pf * ktil over [t0-1, t0]; pole row i+Z at distance d below
    Tq = mp.e**(-I*pi*s/2) * mp.quad(lambda x: Pf(t0 - 1 + x) * ktil(t0 - 1 + x, s),
                                     [0, '0.02', '0.1', '0.5', '0.9', '0.98', 1])
    # closed form
    v = mp.e**(-2 * pi * d)
    T0 = (1 + d)/pi * mp.polylog(1, v) + mp.polylog(2, v)/(2 * pi**2)
    T1 = (1 + d)/pi * Li0(v)          + mp.polylog(1, v)/(2 * pi**2)
    Tcf = r1 * T0 + r2 * T1
    # ---- S-part: regular quad + closed-form FP/PV of the singular part
    ylo, yhi = 1/(1 + d), 1 + d
    reg = mp.quad(lambda y: Pf_minus_sing(I*y) * (I*y)**(s-1) * I,
                  [ylo, '0.999', 1, '1.001', yhi])
    # int sing * tau dtau over the axis  (derived):
    #   a_{-2} [ -(2+d)/d + ln(1+d) ]  +  i a_{-1} [ ln(1+d) + (2d+d^2)/(1+d) ]
    sing_int = am2 * (-(2 + d)/d + mp.log(1 + d)) \
             + I * am1 * (mp.log(1 + d) + (2*d + d**2)/(1 + d))
    Sq = mp.e**(-I*pi*s/2) * reg + mp.e**(-I*pi*s/2) * sing_int
    tot = Tq + Sq
    print("d=%-6s  T:|quad-cf|=%s   L*(Pf,2) = %s   (dev from target %s)"
          % (mp.nstr(d, 3), mp.nstr(abs(Tq - Tcf), 2), mp.nstr(tot, 10),
             mp.nstr(abs(tot - target), 2)), **FL)
