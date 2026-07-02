#!/usr/bin/env python3
r"""
period_polynomial_bases.py  --  companion code for Appendix (explicit period-polynomial
bases) of finite_contour_cocycles_short.tex   (\cite{VerifyScripts}).

For each weight k with dim S_k = 1 (k = 12,16,18,20,22,26) it computes, purely from our
finite-contour L-integral (Definition def:Lint / Theorem thm:weakL):
  - the even/odd Kohnen--Zagier rational period polynomials W_+(X,Y), W_-(X,Y) (primitive
    integer coefficients)   [Kohnen--Zagier, \cite{KZ1984}];
  - the periods omega^+-(Delta_k) of the cusp form Delta_k = f_{k,-1} = q + O(q^2), and the
    quasi-periods eta^+-(\widehat\Delta_k) of the weak cusp form \widehat\Delta_k = f_{k,1}
    = q^{-1} + O(q^2), via  r_{Delta_k}   = sum_pm omega^pm W_pm,
                            r_{\widehat\Delta_k} = sum_pm eta^pm  W_pm.
It also emits the Appendix LaTeX (run with argument 'latex').

Period polynomial (def:rf, overall (2 pi i)^{n+1} dropped, n=k-2):
   r_f = sum_l (-1)^l C(n,l) i^{l+1} L*(f,l+1) X^{n-l} Y^l.
The even part (l even) is pure imaginary = omega^+/eta^+ side; the odd part (l odd) is real
= omega^-/eta^- side.  Run: /path/to/mmf_venv/bin/python period_polynomial_bases.py [latex]
"""
import mpmath as mp, sys
from fractions import Fraction
from math import gcd
mp.mp.dps = 55; I = mp.j; pi = mp.pi

def sig(a, n): return sum(d**a for d in range(1, n+1) if n % d == 0)
NT = 60
# q-series as (offset, [coeffs]): sum coeffs[i] q^{offset+i}
def qmul(A, B):
    oa, a = A; ob, b = B; c = [mp.mpf(0)]*(len(a)+len(b)-1)
    for i, ai in enumerate(a):
        if ai == 0: continue
        for j, bj in enumerate(b): c[i+j] += ai*bj
    return (oa+ob, c)
def qadd(A, B):
    oa, a = A; ob, b = B; off = min(oa, ob); N = max(oa+len(a), ob+len(b))-off
    c = [mp.mpf(0)]*N
    for i, ai in enumerate(a): c[oa-off+i] += ai
    for j, bj in enumerate(b): c[ob-off+j] += bj
    return (off, c)
def qscale(A, s): o, a = A; return (o, [s*x for x in a])
def coeff(A, p): o, a = A; i = p-o; return a[i] if 0 <= i < len(a) else mp.mpf(0)
def qnorm(A):
    o, a = A; i = 0
    while i < len(a) and a[i] == 0: i += 1
    return (o+i, a[i:])
def qdiv(A, B):                                   # A/B, series long division
    oa, a = A; ob, b = B; lead = b[0]
    aa = [x/lead for x in a]; bb = [x/lead for x in b]; res = []
    for i in range(NT+1):
        while len(aa) <= i: aa.append(mp.mpf(0))
        r = aa[i]; res.append(r)
        for j in range(1, len(bb)):
            while len(aa) <= i+j: aa.append(mp.mpf(0))
            aa[i+j] -= r*bb[j]
    return (oa-ob, res)

E4 = (0, [mp.mpf(1)]+[240*sig(3, n) for n in range(1, NT+1)])
E6 = (0, [mp.mpf(1)]+[-504*sig(5, n) for n in range(1, NT+1)])
Delta = qnorm(qscale(qadd(qmul(qmul(E4, E4), E4), qscale(qmul(E6, E6), -1)), mp.mpf(1)/1728))
jser = qdiv(qmul(qmul(E4, E4), E4), Delta)        # j = q^{-1}+744+...
AB = {12:(0,0), 16:(1,0), 18:(0,1), 20:(2,0), 22:(1,1), 26:(2,1)}   # Delta_k = E4^a E6^b Delta

def Delta_k(k):
    a, b = AB[k]; f = Delta
    for _ in range(a): f = qmul(f, E4)
    for _ in range(b): f = qmul(f, E6)
    return f
def Dhat_k(k):                                    # q^{-1}+O(q^2)
    Dk = Delta_k(k); Dkj = qmul(Dk, jser); Dkj2 = qmul(Dkj, jser)
    M = mp.matrix([[coeff(Dkj, 0), coeff(Dk, 0)], [coeff(Dkj, 1), coeff(Dk, 1)]])
    al, be = mp.lu_solve(M, mp.matrix([-coeff(Dkj2, 0), -coeff(Dkj2, 1)]))
    return qadd(qadd(Dkj2, qscale(Dkj, al)), qscale(Dk, be))
def Lstar(F, k, s):                               # thm:weakL sum (c_0 = 0 for both forms)
    o, a = F; tot = mp.mpf(0)
    for i, cn in enumerate(a):
        nn = o+i
        if nn == 0 or cn == 0: continue
        tot += cn*(mp.gammainc(s, 2*pi*nn)/(2*pi*nn)**s + I**k*mp.gammainc(k-s, 2*pi*nn)/(2*pi*nn)**(k-s))
    return tot
def to_frac(x):
    if abs(x) < mp.mpf(10)**(-30): return Fraction(0)
    f = mp.pslq([x, mp.mpf(1)], maxcoeff=10**16, tol=mp.mpf(10)**(-44))
    return Fraction(int(-f[1]), int(f[0])) if f and f[0] != 0 else None
def primitive(fracs):
    L = 1
    for f in fracs: L = L*f.denominator//gcd(L, f.denominator)
    ints = [int(f*L) for f in fracs]; g = 0
    for x in ints: g = gcd(g, abs(x))
    return [x//g for x in ints]

def pairing(p, q, n):        # SL2-invariant (Haberland) pairing on V_n: sum_a (-1)^a p_a q_{n-a}/C(n,a)
    return sum((-1)**a * p[a] * q[n-a] / mp.binomial(n, a) for a in range(n+1))

def compute(k):
    n = k-2; Dk = Delta_k(k); Dh = Dhat_k(k)
    # r_f as a real coefficient vector p_a (coeff of X^{n-a}Y^a): even-a is Im (i^{l+1} phase), odd-a is Re
    def rvec(F):
        c = [(-1)**a*mp.binomial(n, a)*I**(a+1)*Lstar(F, k, a+1) for a in range(n+1)]
        return [c[a].imag if a % 2 == 0 else c[a].real for a in range(n+1)]
    rD, rH = rvec(Dk), rvec(Dh)
    # rational period polynomials W_pm (primitive integers) from the cusp form's even/odd coeffs
    Wp = [0]*(n+1); Wm = [0]*(n+1)
    ep = primitive([to_frac(rD[a]/rD[0]) for a in range(0, n+1, 2)])
    om = primitive([to_frac(rD[a]/rD[1]) for a in range(1, n+1, 2)])
    for i, a in enumerate(range(0, n+1, 2)): Wp[a] = ep[i]
    for i, a in enumerate(range(1, n+1, 2)): Wm[a] = om[i]
    # periods/quasi-periods by Haberland-pairing projection (kills the weak form's C_T coboundary)
    omp = pairing(rD, Wp, n)/pairing(Wp, Wp, n); omm = pairing(rD, Wm, n)/pairing(Wm, Wm, n)
    etp = pairing(rH, Wp, n)/pairing(Wp, Wp, n); etm = pairing(rH, Wm, n)/pairing(Wm, Wm, n)
    residD = max(abs(rD[a]-omp*Wp[a]-omm*Wm[a]) for a in range(n+1))     # cusp form: ~0
    residH = max(abs(rH[a]-etp*Wp[a]-etm*Wm[a]) for a in range(n+1))     # weak form: the C_T coboundary
    Wp = [Wp[a] for a in range(0, n+1, 2)]; Wm = [Wm[a] for a in range(1, n+1, 2)]
    return dict(n=n, Wp=Wp, Wm=Wm, omp=omp, omm=omm, etp=etp, etm=etm, residD=residD, residH=residH)

def poly_tex(coeffs, n, parity):                  # compact (anti)palindromic form
    cmap = {2*idx+parity: c for idx, c in enumerate(coeffs)}   # coeff of X^{n-l}Y^l
    def mono(xp, yp):
        return ("X^{%d}" % xp if xp > 1 else ("X" if xp == 1 else "")) + \
               ("Y^{%d}" % yp if yp > 1 else ("Y" if yp == 1 else ""))
    out = []; mid = None
    for l in sorted(cmap):
        c = cmap[l]
        if c == 0: continue
        if 2*l == n: mid = c; continue                          # middle self-paired term
        if l > n-l: continue                                    # take each pair once
        pm = "+" if cmap.get(n-l, 0) == c else "-"              # W_-: '+', W_+: '-'
        out.append(("+" if c > 0 else "-", abs(c), "(%s %s %s)" % (mono(n-l, l), pm, mono(l, n-l))))
    s = ""
    for i, (sg, a, inner) in enumerate(out):
        s += (("" if i == 0 and sg == "+" else sg+" ") + ("" if a == 1 else str(a)) + inner + " ")
    if mid:
        sg = "+" if mid > 0 else "-"
        s += (sg+" " if s else ("" if mid > 0 else "-")) + ("" if abs(mid) == 1 else str(abs(mid))) + mono(n//2, n//2)
    return s.strip()

if __name__ == "__main__":
    latex = len(sys.argv) > 1 and sys.argv[1] == "latex"
    for k in [12, 16, 18, 20, 22, 26]:
        d = compute(k); n = d["n"]
        if latex:
            print(r"\begin{align}")
            print(r"W_+^{(%d)}(X,Y)&= %s,\\" % (k, poly_tex(d["Wp"], n, 0)))
            print(r"W_-^{(%d)}(X,Y)&= %s,\\" % (k, poly_tex(d["Wm"], n, 1)))
            print(r"r_{\Delta_{%d}}&= %s\,W_+^{(%d)} %s\,W_-^{(%d)},\\" %
                  (k, mp.nstr(d["omp"], 10), k, ("+" if d["omm"] >= 0 else "-")+mp.nstr(abs(d["omm"]), 10), k))
            print(r"r_{\widehat\Delta_{%d}}&= %s\,W_+^{(%d)} %s\,W_-^{(%d)}." %
                  (k, mp.nstr(d["etp"], 10), k, ("+" if d["etm"] >= 0 else "-")+mp.nstr(abs(d["etm"]), 10), k))
            print(r"\end{align}")
        else:
            print("k=%2d  residD(cusp,~0)=%s  residH(=C_T coboundary)=%s" %
                  (k, mp.nstr(d["residD"], 2), mp.nstr(d["residH"], 4)))
            print("   W_+ =", d["Wp"], "\n   W_- =", d["Wm"])
            print("   omega^+- = %s, %s ;  eta^+- = %s, %s ;  det = %s" %
                  (mp.nstr(d["omp"], 10), mp.nstr(d["omm"], 10), mp.nstr(d["etp"], 10), mp.nstr(d["etm"], 10),
                   mp.nstr(d["omp"]*d["etm"]-d["omm"]*d["etp"], 10)))
