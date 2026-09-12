"""Why does the shifted harness give |hat r| = 312510.344 where arc_eight_classes.py recorded
312512.435, for f = Delta/(j-500), k=12?

x is the j-value (repo shorthand: f = Delta/(j-x)).  x = 500 lies in (0,1728) = j(gamma^arc), so
f has an S-pair of simple poles ON the arc at arg z = 1.8484057 and arg Sz = 1.293187.

DIAGNOSIS.  The two harnesses pass on OPPOSITE sides at BOTH poles.
  arc_eight_classes: log r = sigma * 0.15 sin(6(theta - pi/2)), so at theta = 1.8484 log r = +0.149
      (r ~ 1.161, OUTSIDE) and at theta = 1.2932 log r = -0.149 (r ~ 0.861, INSIDE).
  this harness:      r = (1+eps)^{(pi/2 - theta)/(pi/6)}, INSIDE at 1.8484 and OUTSIDE at 1.2932.
Both are S-symmetric (log r odd about pi/2 either way), so both satisfy (1+S); they differ by one
unit of winding at each pole, on BOTH segments since hat gamma^S = hat gamma^T in each.

So the difference should be an ordinary lem:wall_arc crossing, with X_S = X_T = X at each of z, Sz:
    Delta L*_raw(l+1) = 2 pi i X sum_p [ Res_p(f tau^{s-1}) + Res_p(f ktil(.,l+1)) ],
    Delta Phi         = 2 pi i X sum_p Res_p f,
    Delta r_f[l]      = (2 pi i)^{n+1} (-1)^l C(n,l) Delta L*_raw(l+1),
    Delta hat r       = Delta r_f - Delta Phi * r_{E_k}.
Every residue is a small circle about z or Sz -- nothing from the contour machinery -- so this is a
prediction, not a consistency check.

PREDICTIONS.
 (1) Reproducing their wiggle contour here gives |hat r| = 312512.435, matching the recorded value.
     That is the load-bearing check: it says the two scripts share conventions and the gap is
     geometric, not normalisation (a convention mismatch would not agree to six figures).
 (2) The spiral gives 312510.344, as before.
 (3) hat r(wiggle) - hat r(spiral) equals the wall-crossing above at X = +1 or -1 (sign read off
     the data), to harness precision.
 (4) Both are in W: the flip is S-symmetric on both segments, so neither relation cares.
eps = 1e-1, per the conditioning result (the limit is exact; small eps only costs digits).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, Delta, jay, ktil, E12,                   # noqa: E402
                    rvec, inW, check_orientation)

mp.mp.dps = 30
K, N = 12, 10
X = mp.mpf(500)
EPS = mp.mpf('1e-1')
f = lambda t: Delta(t) / (jay(t) - X)


def theta_of_x():
    """bisect in theta on [pi/3, pi/2] for j(e^{i theta}) = X; j is real and monotone there."""
    lo, hi = pi / 3, pi / 2
    for _ in range(200):
        mid = (lo + hi) / 2
        if mp.re(jay(mp.e**(I * mid))) < X:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def nodes(a, b, n, depth, centres):
    xs = [a + (b - a) * mp.mpf(j) / n for j in range(n + 1)]
    for c in centres:
        d = abs(b - a) / n
        for _ in range(depth):
            d /= 2
            xs += [c - d, c + d]
        xs.append(c)
    return sorted(set(x for x in xs if min(a, b) <= x <= max(a, b)))


def wiggle(theta):
    return mp.e**(mp.mpf('0.15') * mp.sin(6 * (theta - pi / 2)) + I * theta)


def dwiggle(theta):
    return wiggle(theta) * (mp.mpf('0.15') * 6 * mp.cos(6 * (theta - pi / 2)) + I)


def spir(theta):
    L = mp.log(1 + EPS)
    return mp.e**(L * (pi / 2 - theta) / (pi / 6) + I * theta)


def dspir(theta):
    L = mp.log(1 + EPS)
    return spir(theta) * (-L / (pi / 6) + I)


def pathint(g, which, th):
    """theta DECREASING from 2pi/3 to pi/3, clustered at the two pole angles."""
    P, dP = (wiggle, dwiggle) if which == 'wiggle' else (spir, dspir)
    ns = nodes(2 * pi / 3, pi / 3, 24, 12, th)
    ns = sorted(ns, reverse=True)
    return sum(mp.quad(lambda x: g(P(x)) * dP(x), [u, v]) for u, v in zip(ns[:-1], ns[1:]))


def hatr(which, th, rE):
    out = []
    for l in range(N + 1):
        s = mp.mpf(l + 1)
        L = pathint(lambda t: f(t) * t**(s - 1), which, th) + \
            pathint(lambda t: f(t) * ktil(t, s, K), which, th)
        out.append((2 * pi * I)**(N + 1) * (-1)**l * mp.binomial(N, l) * L)
    Ph = pathint(f, which, th)
    return [a - Ph * b for a, b in zip(out, rE)], Ph


def residue(g, c, r, n=48):
    ps = [2 * pi * mp.mpf(j) / n for j in range(n + 1)]
    tot = sum(mp.quad(lambda p: g(c + r * mp.e**(I * p)) * I * r * mp.e**(I * p), [u, v])
              for u, v in zip(ps[:-1], ps[1:]))
    return tot / (2 * I * pi)


print("guard: Phi(E_12) plain arc = %s" % mp.nstr(check_orientation(), 12), flush=True)
t1 = theta_of_x()
z, Sz = mp.e**(I * t1), mp.e**(I * (pi - t1))
th = sorted([t1, pi - t1])
print("poles: arg z = %s, arg Sz = %s   (recorded 1.293187, 1.8484057)"
      % (mp.nstr(t1, 10), mp.nstr(pi - t1, 10)), flush=True)
print("  j(z) = %s" % mp.nstr(jay(z), 12), flush=True)
print("  wiggle r at those angles: %s, %s"
      % (mp.nstr(abs(wiggle(th[0])), 6), mp.nstr(abs(wiggle(th[1])), 6)), flush=True)
print("  spiral r at those angles: %s, %s"
      % (mp.nstr(abs(spir(th[0])), 6), mp.nstr(abs(spir(th[1])), 6)), flush=True)

rE = rvec(E12, K)
rw, Phw = hatr('wiggle', th, rE)
rs, Phs = hatr('spiral', th, rE)
for tag, v, Ph in (("wiggle", rw, Phw), ("spiral", rs, Phs)):
    s_, u_ = inW(v, K)
    print("%-7s |hat r| = %-18s Phi = %-24s |(1+S)| = %-11s |(1+U+U^2)| = %s"
          % (tag, mp.nstr(max(abs(x) for x in v), 12), mp.nstr(Ph, 10),
             mp.nstr(s_, 5), mp.nstr(u_, 5)), flush=True)
print("   recorded by arc_eight_classes.py: 312512.435", flush=True)

rr = mp.mpf('0.05')
a_tot = residue(f, z, rr) + residue(f, Sz, rr)
print("Phi(wiggle) - Phi(spiral) = %-26s  2 pi i sum_p Res_p f = %s"
      % (mp.nstr(Phw - Phs, 12), mp.nstr(2 * I * pi * a_tot, 12)), flush=True)

for Xw in (1, -1):
    dPh = 2 * I * pi * Xw * a_tot
    pred = []
    for l in range(N + 1):
        s = mp.mpf(l + 1)
        rsum = sum(residue(lambda t, s=s: f(t) * t**(s - 1), c, rr)
                   + residue(lambda t, s=s: f(t) * ktil(t, s, K), c, rr) for c in (z, Sz))
        pred.append((2 * pi * I)**(N + 1) * (-1)**l * mp.binomial(N, l) * 2 * I * pi * Xw * rsum)
    dhat = [a - dPh * b for a, b in zip(pred, rE)]
    meas = [a - b for a, b in zip(rw, rs)]
    rel = max(abs(a - b) for a, b in zip(meas, dhat)) / max(abs(a) for a in meas)
    print("X = %+d : predicted vs measured hat r difference, relative = %s"
          % (Xw, mp.nstr(rel, 8)), flush=True)
