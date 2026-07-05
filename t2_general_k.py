#!/usr/bin/env python3
r"""t2_general_k.py -- GENERAL (k,s) wall value of the projection, pole at i, order P.
  lim_{d->0} L*(P_i f, s; tau0=i(1+d)) = sum_{m=1}^P r_m FP_{m-1}(s,k),
  r_m = (-2 pi i)^m a_{-m}/(m-1)!,
  FP_N = pi_N H_N/(2pi)^{N+1} + sum_{j != N} c_j^+ zeta(j+1-N)/(2pi)^{j+1} + E_N^+
         + (-1)^N N!/(2pi)^{N+1} sum_{t<N} (-1)^t (s-1)_(t)/(t! (N-t)),
  c_j^+ = (s-1)_(j) + i^k (k-s-1)_(j),  pi_N = c_N^+.
Forms: Delta/E6 (k=6,P=1), E4^2 Delta/E6^2 (k=8,P=2), Delta^2/E6^3 (k=6,P=3).
Also checks the generalized parity a_{-1} = i(k-2)/2 a_{-2} for P=2.
"""
import mpmath as mp

mp.mp.dps = 25
I, pi = mp.j, mp.pi
FL = dict(flush=True)

def sigma(a, m): return sum(d**a for d in range(1, m + 1) if m % d == 0)
NT = 60
_e4 = [240 * sigma(3, m) for m in range(NT, 0, -1)] + [mp.mpf(1)]
_e6 = [-504 * sigma(5, m) for m in range(NT, 0, -1)] + [mp.mpf(1)]
def q(t): return mp.e**(2 * pi * I * t)
def E4(t): return mp.polyval(_e4, q(t))
def E6(t): return mp.polyval(_e6, q(t))
def Delta(t): return (E4(t)**3 - E6(t)**2) / 1728

def poch(a, j):
    out = mp.mpc(1)
    for r in range(j): out *= (a - r)
    return out

def laurent(f, P, r_c=mp.mpf('0.05'), Np=96):
    """a_{-1..-P} at tau=i by circle sums."""
    a = {}
    for m in range(1, P + 1):
        a[m] = sum(f(I + r_c*mp.e**(I*2*pi*j/Np)) * (r_c*mp.e**(I*2*pi*j/Np))**m
                   for j in range(Np)) / Np
    return a

def FP(N, s, k, J=12, nmax=50):
    ik = I**k
    cj = lambda j: poch(s-1, j) + ik * poch(k-s-1, j)
    HN = sum(mp.mpf(1)/r for r in range(1, N+1))
    tot = cj(N) * HN / (2*pi)**(N+1)
    for j in range(J):
        if j != N:
            tot += cj(j) * mp.zeta(j+1-N) / (2*pi)**(j+1)
    for n in range(1, nmax+1):
        x = 2*pi*n
        tot += n**N * mp.e**x * (poch(s-1, J) * mp.gammainc(s-J, x, mp.inf) / x**s
                                 + ik * poch(k-s-1, J) * mp.gammainc(k-s-J, x, mp.inf) / x**(k-s))
    tot += (-1)**N * mp.factorial(N) / (2*pi)**(N+1) * \
           sum((-1)**t * poch(s-1, t) / (mp.factorial(t) * (N-t)) for t in range(N))
    return tot

def ktil(t, s, k): return mp.zeta(1-s, t+1) - mp.e**(I*pi*(s-1)) * mp.zeta(1-(k-s), t+1)

def run(name, f, k, P, s):
    print("\n== %s  (k=%d, P=%d, s=%s) ==" % (name, k, P, mp.nstr(s, 6)), **FL)
    a = laurent(f, P)
    if P >= 2:
        pr = a[1] / a[2]
        print("   parity a_{-1}/a_{-2} = %s   (predict i(k-2)/2 = %s)"
              % (mp.nstr(pr, 8), mp.nstr(I*(k-2)/mp.mpf(2), 8)), **FL)
    r = {m: (-2*pi*I)**m * a[m] / mp.factorial(m-1) for m in range(1, P+1)}
    pred = sum(r[m] * FP(m-1, s, k) for m in range(1, P+1))
    print("   predicted wall value = %s" % mp.nstr(pred, 12), **FL)

    # tails and their singular parts
    def Pf(t):
        u = 2*pi*I*(t - I)
        e = mp.e**u
        tot = mp.mpc(0)
        li = e/(1-e)                       # Li_0
        tot += r[1] * li
        if P >= 2: tot += r[2] * e/(1-e)**2
        if P >= 3: tot += r[3] * (e + e**2)/(1-e)**3
        return tot
    def sing(t):
        return sum(a[m] * (t - I)**(-m) for m in range(1, P+1))
    def Pf_reg(t):
        u = 2*pi*I*(t - I)
        if abs(u) < mp.mpf('0.05'):        # zeta-series of regular parts
            tot = mp.mpc(0)
            for m in range(1, P+1):
                N = m-1
                tot += r[m] * sum(mp.zeta(-N-p)/mp.factorial(p) * u**p for p in range(8))
            return tot
        return Pf(t) - sing(t)

    for d in (mp.mpf('0.1'), mp.mpf('0.05'), mp.mpf('0.025')):
        t0 = I*(1+d)
        Tq = mp.e**(-I*pi*s/2) * mp.quad(lambda x: Pf(t0-1+x)*ktil(t0-1+x, s, k),
                                         [0, '0.02', '0.1', '0.5', '0.9', '0.98', 1])
        ylo, yhi = 1/(1+d), 1+d
        reg = mp.quad(lambda y: y**(s-1) * Pf_reg(I*y), [ylo, '0.999', 1, '1.001', yhi])
        # FP/PV of  int y^{s-1} sing(iy) dy ;  sing(iy) = sum_m a_m i^{-m} (y-1)^{-m}
        sing_val = mp.mpc(0)
        for m in range(1, P+1):
            coef = a[m] * I**(-m)
            # int y^{s-1}(y-1)^{-m}: subtract Taylor of y^{s-1} to order m-1
            def freg(y, m=m):
                h = y - 1
                if abs(h) < mp.mpf('1e-4'):
                    return poch(s-1, m)/mp.factorial(m) + poch(s-1, m+1)/mp.factorial(m+1)*h
                return (y**(s-1) - sum(poch(s-1, t)/mp.factorial(t)*h**t
                                       for t in range(m))) / h**m
            val = mp.quad(freg, [ylo, '0.999', 1, '1.001', yhi])
            for t in range(m):
                bt = poch(s-1, t)/mp.factorial(t)
                p = t - m
                if p == -1:
                    val += bt * mp.log(1+d)                     # PV
                else:
                    val += bt * ((yhi-1)**(p+1) - (ylo-1)**(p+1)) / (p+1)   # FP endpoint
            sing_val += coef * val
        tot = Tq + reg + sing_val
        print("   d=%-6s  L*(Pf) = %s   dev = %s"
              % (mp.nstr(d, 3), mp.nstr(tot, 10), mp.nstr(abs(tot - pred), 2)), **FL)

s0 = mp.mpc('2.3', '0.4')
run("Delta/E6",        lambda t: Delta(t)/E6(t),            6, 1, s0)
run("E4^2 Delta/E6^2", lambda t: E4(t)**2*Delta(t)/E6(t)**2, 8, 2, s0)
run("Delta^2/E6^3",    lambda t: Delta(t)**2/E6(t)**3,       6, 3, s0)
