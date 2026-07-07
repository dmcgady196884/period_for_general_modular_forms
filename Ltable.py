#!/usr/bin/env python3
r"""Ltable.py -- DAM's explicit mode-contribution table, s=1.
Columns: k | f | c_f(n)G(n) for n=-1,1,2,3,4,5 | sum_displayed | L*(f,1) | ratio
G_{s,k}(n) = Gamma(s,2 pi n)/(2 pi n)^s + i^k Gamma(k-s,2 pi n)/(2 pi n)^{k-s}   (eq:gc).
Two rows per k (Delta_k, Dhat_k).  'latex' arg emits the LaTeX tabular."""
import mpmath as mp, sys
sys.path.insert(0, "/Users/dmcgady/Documents/math/meromorphic_modular_forms/quasi_period_success")
from period_polynomial_bases import Delta_k, Dhat_k, Lstar, coeff
mp.mp.dps = 30; I = mp.j; pi = mp.pi
s = 1
MODES = [-1, 1, 2, 3, 4, 5]
latex = len(sys.argv) > 1 and sys.argv[1] == "latex"

def G(n, k):
    return mp.gammainc(s, 2*pi*n)/(2*pi*n)**s + I**k*mp.gammainc(k-s, 2*pi*n)/(2*pi*n)**(k-s)
def r(x, d=5):
    v = mp.re(x)
    return "0" if v == 0 else mp.nstr(v, d)
def tex(x, d=5):
    v = mp.re(x)
    if v == 0: return "0"
    s = mp.nstr(v, d)
    if "e" in s:
        m, e = s.split("e")
        return r"%s{\times}10^{%d}" % (m, int(e))
    return s

rows = []
for k in [12, 16, 18, 20, 22, 26]:
    Dk, Dh = Delta_k(k), Dhat_k(k)
    LD, LH = Lstar(Dk, k, s), Lstar(Dh, k, s)
    ratio = mp.re(LH/LD)
    for nm, F, L in [(r"\Delta_{%d}" % k, Dk, LD), (r"\widehat\Delta_{%d}" % k, Dh, LH)]:
        contribs = [coeff(F, n)*G(n, k) for n in MODES]
        rows.append((k, nm, contribs, sum(contribs), L, ratio))

if not latex:
    hdr = "%3s %-9s " % ("k", "f") + "".join("%12s" % ("n=%+d" % n) for n in MODES) + \
          "%12s%12s%10s" % ("SUM", "L*(f,1)", "ratio")
    print(hdr)
    for (k, nm, contribs, S, L, ratio) in rows:
        fnm = nm.replace(r"\widehat", "^").replace(r"\Delta", "D").replace("_{%d}" % k, "")
        print("%3d %-9s " % (k, fnm) + "".join("%12s" % r(c) for c in contribs) +
              "%12s%12s%10s" % (r(S), r(L), r(ratio, 6)))
else:
    ncol = "c" * (2 + len(MODES) + 3)
    print(r"\begin{tabular}{%s}" % ncol)
    print(r"\toprule")
    head = r"$k$ & $f$ & " + " & ".join(r"$c_f(%d)G(%d)$" % (n, n) for n in MODES) + \
           r" & $\sum_{\rm disp}$ & $L^*(f,1)$ & ratio\\"
    print(head)
    print(r"\midrule")
    for i, (k, nm, contribs, S, L, ratio) in enumerate(rows):
        kcol = str(k) if i % 2 == 0 else ""
        cells = " & ".join("$%s$" % tex(c) for c in contribs)
        print(r"%s & $%s$ & %s & $%s$ & $%s$ & $%s$\\" % (kcol, nm, cells, tex(S), tex(L), r(ratio)))
        if i % 2 == 1 and i < len(rows)-1: print(r"\midrule")
    print(r"\bottomrule")
    print(r"\end{tabular}")
