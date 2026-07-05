#!/usr/bin/env python3
r"""t2_k4_generic.py -- generic-s 4|k ledger, k=4 block, option-c.
Claim: lim_{d->0} L*(P_i f, s) = a_{-2}[ 13 pi/6 - 1 - pi (s-2)^2/6 + 2 pi E0(s) - 4 pi^2 E1(s) ],
E_N(s) = sum_n n^N e^{2 pi n} [ (s-1)_(N+2) G(s-N-2, 2 pi n)/(2 pi n)^s
                              + (3-s)_(N+2) G(2-s-N, 2 pi n)/(2 pi n)^{4-s} ]   (upper G).
Numerics: E_N via deeper tower J=12 (terms O(n^-13)); L*(P_i f, s; d) by quadrature
(T: direct; S: singular-part subtraction with exact FP/PV), d = 0.1, 0.05, 0.025.
"""
import mpmath as mp

mp.mp.dps = 25
I, pi = mp.j, mp.pi
FL = dict(flush=True)
k = 4

am2 = mp.mpf(-1) / (1728 * pi**2)
am1 = I * am2
r1, r2 = 2 * pi * am2, -4 * pi**2 * am2

def Li0(u):  return u / (1 - u)
def Lim1(u): return u / (1 - u)**2
def Pf(t):
    e = mp.e**(2 * pi * I * (t - I))
    return r1 * Li0(e) + r2 * Lim1(e)
def Pf_reg(t):                                  # Pf minus its pole part, series near i
    u = 2 * pi * I * (t - I)
    if abs(u) < mp.mpf('0.05'):
        return r1 * (-mp.mpf(1)/2 - u/12 + u**3/720) + r2 * (-mp.mpf(1)/12 + u**2/240)
    return Pf(t) - (am2 * (t - I)**(-2) + am1 * (t - I)**(-1))
def ktil(t, s): return mp.zeta(1 - s, t + 1) - mp.e**(I * pi * (s - 1)) * mp.zeta(1 - (k - s), t + 1)

def poch(a, j):
    out = mp.mpc(1)
    for r in range(j): out *= (a - r)
    return out

def Eplus(N, s, J=12, nmax=50):
    """E_N via tower: explicit zeta terms j=N+2..J-1 plus Gamma remainder at J."""
    tot = mp.mpc(0)
    for j in range(N + 2, J):
        cj = poch(s - 1, j) + poch(k - s - 1, j)
        tot += cj / (2 * pi)**(j + 1) * mp.zeta(j + 1 - N)
    for n in range(1, nmax + 1):
        x = 2 * pi * n
        t = (poch(s - 1, J) * mp.gammainc(s - J, x, mp.inf) / x**s
             + poch(k - s - 1, J) * mp.gammainc(k - s - J, x, mp.inf) / x**(k - s))
        tot += n**N * mp.e**(x) * t
    return tot

def predicted(s):
    return am2 * (13*pi/6 - 1 - pi*(s - 2)**2/6) + r1*2*pi/(2*pi)*0 \
           + 2*pi*am2 * 0  # assembled below properly
def pred(s):
    return am2 * (13*pi/6 - 1 - pi*(s-2)**2/6 + 2*pi*Eplus(0, s) - 4*pi**2*Eplus(1, s))

def Lquad(s, d):
    t0 = I * (1 + d)
    Tq = mp.e**(-I*pi*s/2) * mp.quad(lambda x: Pf(t0 - 1 + x) * ktil(t0 - 1 + x, s),
                                     [0, '0.02', '0.1', '0.5', '0.9', '0.98', 1])
    ylo, yhi = 1/(1 + d), 1 + d
    # S = int y^{s-1} Pf(iy) dy  (prefactor absorbed; see derivation)
    reg = mp.quad(lambda y: y**(s-1) * Pf_reg(I*y), [ylo, '0.999', 1, '1.001', yhi])
    # -am2 * FP int y^{s-1}(y-1)^{-2}  - i am1 * PV int y^{s-1}(y-1)^{-1}
    def f2(y):
        h = y - 1
        if abs(h) < mp.mpf('1e-4'):
            return (s-1)*(s-2)/2 + (s-1)*(s-2)*(s-3)/6 * h
        return (y**(s-1) - 1 - (s-1)*h) / h**2
    def f1(y):
        h = y - 1
        if abs(h) < mp.mpf('1e-4'):
            return (s-1) + (s-1)*(s-2)/2 * h
        return (y**(s-1) - 1) / h
    fp2 = mp.quad(f2, [ylo, '0.999', 1, '1.001', yhi]) \
          + (-(2 + d)/d) + (s - 1) * mp.log(1 + d)
    pv1 = mp.quad(f1, [ylo, '0.999', 1, '1.001', yhi]) + mp.log(1 + d)
    Sq = reg - am2 * fp2 - I * am1 * pv1
    return Tq + Sq

for s in (mp.mpc('2.6', '0.7'), mp.mpf('1.4')):
    p = pred(s)
    print("s = %s   predicted limit = %s" % (mp.nstr(s, 6), mp.nstr(p, 12)), **FL)
    for d in (mp.mpf('0.1'), mp.mpf('0.05'), mp.mpf('0.025')):
        v = Lquad(s, d)
        print("   d=%-6s  L*(Pf,s) = %s   dev = %s"
              % (mp.nstr(d, 3), mp.nstr(v, 10), mp.nstr(abs(v - p), 2)), **FL)
