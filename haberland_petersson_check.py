#!/usr/bin/env python3
r"""
haberland_petersson_check.py  --  companion numerics for section 3 of
    finite_contour_cocycles_short.tex   (\cite{VerifyScripts}, eq:haberland)

NORMALIZATION-INDEPENDENT external check on the finite-contour periods.

The periods read off the two-segment L-integral,
    r_m(f) = int_0^{i oo} tau^m f(tau) dtau = i^{m+1} L*(f, m+1),
are the classical critical periods, so they must satisfy Haberland's formula
[H. Cohen, "Haberland's formula and numerical computation of Petersson scalar
products", ANTS X, Open Book Ser. 1 (2013), Thm 5.2(2); full modular group, f|T=f]:

  -6 (-2i)^{k-2} <f,f>
        = sum_{m+n<=k-2} C(k-2,m+n) C(m+n,m) (-1)^m Im( r_m(f) conj(r_n(f)) ).

This recovers the Petersson norm <f,f> from the period polynomial and does NOT depend
on any period-polynomial basis normalization (unlike the determinant ratio alpha_k of
Table 1).  We check it two ways:

  [A] k=12:  our L*-values reproduce Zagier's <Delta,Delta> (Cohen p.2, 47 digits)
             to ~14 digits (q-series-truncation limited).
  [B] k=16:  our L*-Haberland value matches the direct fundamental-domain integral
             <f,f> = int_F |f|^2 y^{k-2} dx dy   to ~16 digits (independent method).

Conventions.  L*(f,s) here is the completed period Lambda(f,s) = int_0^oo t^{s-1} f(it) dt,
computed by the upper-incomplete-gamma split that is the s-plane form of def:Lint for a
cusp form (c_f(0)=0); r_m = i^{m+1} L*(m+1).
Run:  /path/to/mmf_venv/bin/python haberland_petersson_check.py   (needs mpmath).
"""
import mpmath as mp
mp.mp.dps = 45
I, pi = mp.j, mp.pi

def sigma(a, n): return sum(d**a for d in range(1, n + 1) if n % d == 0)
NT = 90
def qmul(A, B):
    oa, a = A; ob, b = B; c = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0: continue
        for j, bj in enumerate(b): c[i + j] += ai * bj
    return (oa + ob, c)
def qadd(A, B):
    oa, a = A; ob, b = B; off = min(oa, ob); N = max(oa + len(a), ob + len(b)) - off
    c = [mp.mpf(0)] * N
    for i, ai in enumerate(a): c[oa - off + i] += ai
    for j, bj in enumerate(b): c[ob - off + j] += bj
    return (off, c)
def qscale(A, s): o, a = A; return (o, [s * x for x in a])
def qnorm(A):
    o, a = A; i = 0
    while i < len(a) and a[i] == 0: i += 1
    return (o + i, a[i:])
E4 = (0, [mp.mpf(1)] + [240 * sigma(3, n) for n in range(1, NT + 1)])
E6 = (0, [mp.mpf(1)] + [-504 * sigma(5, n) for n in range(1, NT + 1)])
Delta = qnorm(qscale(qadd(qmul(qmul(E4, E4), E4), qscale(qmul(E6, E6), -1)), mp.mpf(1) / 1728))
AB = {12: (0, 0), 16: (1, 0), 18: (0, 1), 20: (2, 0), 22: (1, 1), 26: (2, 1)}
def Delta_k(k):                                   # normalized weight-k cusp form q + O(q^2)
    a, b = AB[k]; f = Delta
    for _ in range(a): f = qmul(f, E4)
    for _ in range(b): f = qmul(f, E6)
    return qnorm(f)

def Lstar(F, k, s):                               # = completed period Lambda(f,s)
    o, a = F; tot = mp.mpf(0)
    for i, cn in enumerate(a):
        n = o + i
        if n <= 0 or cn == 0: continue
        tot += cn * (mp.gammainc(s, 2*pi*n)/(2*pi*n)**s + I**k * mp.gammainc(k-s, 2*pi*n)/(2*pi*n)**(k-s))
    return tot

def haberland(k):                                 # <f,f> from the period polynomial
    n = k - 2; f = Delta_k(k)
    r = [I**(m+1) * Lstar(f, k, m+1) for m in range(n+1)]
    S = mp.mpf(0)
    for m in range(n+1):
        for nn in range(n+1):
            if m + nn > n: continue
            S += mp.binomial(n, m+nn) * mp.binomial(m+nn, m) * (-1)**m * (r[m]*mp.conj(r[nn])).imag
    return (-S / (6 * (-2*I)**(k-2))).real

def petersson_direct(k, dps=25, nt=32):           # <f,f> = int_F |f|^2 y^{k-2} dx dy
    with mp.workdps(dps):
        f = Delta_k(k); o, a = f; a = a[:nt]
        def fval(x, y):
            q = mp.e**(2*pi*I*x) * mp.e**(-2*pi*y); acc = mp.mpf(0)
            for c in reversed(a): acc = acc*q + c
            return acc * q**o
        def inner(x):
            y0 = mp.sqrt(1 - x*x)
            return mp.quad(lambda y: abs(fval(x, y))**2 * y**(k-2), [y0, y0+mp.mpf('0.4'), y0+1, y0+3, mp.inf])
        return (2 * mp.quad(inner, [0, mp.mpf('0.25'), mp.mpf('0.5')])).real   # |f| even in x

if __name__ == "__main__":
    zag = mp.mpf("1.03536205680430209223478168122251645932249") * mp.mpf(10)**-6   # Cohen p.2
    hb12 = haberland(12)
    print("[A] k=12  Haberland(from L*) = %s" % mp.nstr(hb12, 20))
    print("          Zagier <Delta,Delta> = %s" % mp.nstr(zag, 20))
    print("          rel. error           = %s" % mp.nstr(abs(hb12 - zag)/zag, 3))
    hb16 = haberland(16); di16 = petersson_direct(16)
    print("[B] k=16  Haberland(from L*) = %s" % mp.nstr(hb16, 18))
    print("          direct integral     = %s" % mp.nstr(di16, 18))
    print("          rel. error          = %s" % mp.nstr(abs(hb16 - di16)/hb16, 3))
