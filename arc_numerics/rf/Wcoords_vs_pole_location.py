"""Item D, redone with the right diagnostic: the coordinates of hat r_{f_z} inside W.

The previous attempt (parity_vs_pole_location.py) reported max|even l| and max|odd l| of the raw
coefficient vector and got odd/even = 1.13454 at EVERY z.  That was the diagnostic, not the
mathematics: r_f[l] ~ (-1)^l C(n,l) i^{l+1} L*(l+1), so when L* varies slowly in l the profile is
dominated by C(10,l), which does not depend on z at all (C(10,5)/C(10,4) = 1.2 against the
measured 1.135).  It also asked the wrong question: tau -> -conj(tau) is ANTIholomorphic, so its
fixed locus forces hat r to be REAL in a suitable basis, not to lie in a parity eigenspace.

Here instead: decompose in the basis of W.  At k = 12, dim W = 3 = one odd direction W_- plus two
even, W_+ and p_0, and common.decompose12 reads the three coordinates (a, alpha, beta) straight
off.  Then the two questions of item D are:
  (a) REALITY.  For z on the self-conjugate loci -- the imaginary axis, where -conj z = z, and
      Re z = +-1/2, where -conj z = z -+ 1 and so is T-equivalent -- are the coordinates real
      once the overall scale is divided out?
  (b) DIRECTION.  hat r_{f_z} scales with the residue E_10(z), so the scale carries no
      information; the content is the point [a : alpha : beta] in P(W).  Does it move with z?
      If it is constant, hat r_{f_z} sweeps a single LINE in W and the "motion of the W^+-
      coefficients" is trivial.  If it moves, that motion is the answer to D.

Normalisation: coordinates are divided by beta (the p_0 component), giving (a/beta, alpha/beta).
CALIBRATION rows decompose r_Delta and r_{E_12} in the same basis.  This bears on the reference-
form question too: the ambiguity is Phi(f) r(S_k^!), and if r of a cusp form has NO p_0 component
then r(S_k^!) misses that direction, which is the codimension-one claim.

f_z = E_10 j'/(j - j(z)) = -2 pi i j E_6^2/(j - x),  x = j(z),  k = 12, residue E_10(z) at z.
All z sit well inside the fundamental domain so the poles are far from the arc (cf. TRAP 4).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, E4, E6, Delta, jay, ktil, E12,            # noqa: E402
                    arcint, rvec, inW, decompose12, check_orientation)

mp.mp.dps = 30
K, N = 12, 10
rE12 = rvec(E12, K)


def rhat(f, depth=11):
    L = lambda s: mp.e**(-I * pi * s / 2) * arcint(
        lambda t: f(t) * (t**(s - 1) + ktil(t, s, K)), depth)
    rf = [(2 * pi * I)**(N + 1) * (-1)**l * mp.binomial(N, l) * I**(l + 1)
          * L(mp.mpf(l + 1)) for l in range(N + 1)]
    Ph = arcint(f, depth)
    return [a - Ph * b for a, b in zip(rf, rE12)], Ph


def show_coords(tag, v, extra=""):
    a, al, be = decompose12(v)
    s_, u_ = inW(v, K)
    print("   %-30s W_- = %-24s W_+ = %-24s p_0 = %s"
          % (tag, mp.nstr(a, 10), mp.nstr(al, 10), mp.nstr(be, 10)), flush=True)
    if abs(be) > mp.mpf('1e-40'):
        ra, rb = a / be, al / be
        print("       normalised by p_0:  a/beta = %-26s alpha/beta = %-26s"
              % (mp.nstr(ra, 12), mp.nstr(rb, 12)), flush=True)
        print("       imaginary parts rel: %-14s %-14s   in W: %s / %s  %s"
              % (mp.nstr(abs(mp.im(ra)) / max(abs(ra), mp.mpf('1e-300')), 6),
                 mp.nstr(abs(mp.im(rb)) / max(abs(rb), mp.mpf('1e-300')), 6),
                 mp.nstr(s_, 4), mp.nstr(u_, 4), extra), flush=True)
    return (a, al, be)


print("guard: Phi(E_12) = %s" % mp.nstr(check_orientation(), 12), flush=True)
print("=" * 80, flush=True)
print("CALIBRATION", flush=True)
show_coords("r_Delta", rvec(Delta, K))
show_coords("r_{E_12}", rE12)

print("=" * 80, flush=True)
print("hat r_{f_z},  f_z = -2 pi i j E_6^2/(j - j(z))", flush=True)
ZS = [("z = 1.2i           (imag axis)", mp.mpf('1.2') * I),
      ("z = 1.5i           (imag axis)", mp.mpf('1.5') * I),
      ("z = 2.0i           (imag axis)", 2 * I),
      ("z = 0.5+1.2i       (vert edge)", mp.mpf('0.5') + mp.mpf('1.2') * I),
      ("z = 0.5+1.6i       (vert edge)", mp.mpf('0.5') + mp.mpf('1.6') * I),
      ("z = 0.2+1.3i       (generic)", mp.mpf('0.2') + mp.mpf('1.3') * I),
      ("z = 0.35+1.5i      (generic)", mp.mpf('0.35') + mp.mpf('1.5') * I)]

dirs = {}
for lbl, z in ZS:
    f = lambda t, x=jay(z): -2 * I * pi * jay(t) * E6(t)**2 / (jay(t) - x)
    v, Ph = rhat(f)
    a, al, be = show_coords(lbl, v)
    if abs(be) > mp.mpf('1e-40'):
        dirs[lbl] = (a / be, al / be)

print("=" * 80, flush=True)
print("DIRECTION in P(W): does [a:alpha:beta] move with z?", flush=True)
ks = list(dirs)
base = dirs[ks[0]]
for lbl in ks:
    d = dirs[lbl]
    drift = max(abs(d[0] - base[0]) / max(abs(base[0]), mp.mpf('1e-300')),
                abs(d[1] - base[1]) / max(abs(base[1]), mp.mpf('1e-300')))
    print("   %-32s drift from first row = %s" % (lbl, mp.nstr(drift, 8)), flush=True)
