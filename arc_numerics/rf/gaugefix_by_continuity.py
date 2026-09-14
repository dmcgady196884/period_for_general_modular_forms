"""Can continuity across the j-cut gauge-fix the reference form in def:rfhat?

DAM's idea.  hat r_f = r_f - Phi(f) r_g depends on the reference form g only through
Phi(f) r_{g-g'}, g-g' in S_k.  Take f_x = Delta/(j-x) at k=12 (simple pole orbit, c_f(0)=0).
For x off the real segment [0,1728] = j(gamma^arc) the poles are off the contour and hat r is
defined; crossing the cut is a wall-crossing and BOTH r_f and Phi jump.  Continuity of hat r
would require

    [hat r] = [r_f] - [Phi] r_g = 0      i.e.      r_g = [r_f] / [Phi],

which is ONE vector equation per x and ONE unknown g.  So the question is not whether some g
works at a single x -- generically one will -- but whether the SAME g works at every x.

WHAT THIS TESTS, precisely:
  (1) is [Phi] non-zero?  If it vanishes the jump carries no g-dependence and the idea is dead
      on arrival, with no computation of r_g needed.
  (2) is v(x) := [r_f](x) / [Phi](x) the same vector at different x?  If it drifts, the system
      is over-determined and continuity CANNOT gauge-fix -- decisive against.
  (3) if v is x-independent, is it an admissible reference?  At k=12 the admissible set is
      r_{E_12} + C r_Delta (any g in M_12 with c_g(0)=1), so the check is whether
      v - r_{E_12} is proportional to r_Delta.

PREDICTION (stated before running).  v(x) DRIFTS with x, and the idea fails.  Reason: the jump
in r_f comes from lem:wall_arc as
    Delta L*(s) = 2 pi i e^{-i pi s/2} a_p [ X_S p^{s-1} + X_T ktil(p,s) ],
which assembles to something carrying the POLE LOCATION p, while [Phi] = 2 pi i X_T sum_p a_p
carries only the residue.  Dividing one by the other leaves p-dependence, and p moves with x.
So I expect the ratio to vary, and the informative number is HOW MUCH -- if it drifts only
slightly, some weakened version (a best-fit g, or a condition at one distinguished x) may still
be worth having.

Caveats.  x is taken off the cut by +- i delta with delta large enough that the poles stay clear
of the arc -- the jump is the delta -> 0 limit, so two deltas are used and the drift in v is
reported against the drift in delta, to separate a real x-dependence from an unconverged limit.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, E4, Delta, jay, ktil, E12,                # noqa: E402
                    arcint, rvec, check_orientation)

mp.mp.dps = 30
K, N = 12, 10


def Lstar(f, s, depth=11):
    g = lambda t: f(t) * (t**(s - 1) + ktil(t, s, K))
    return mp.e**(-I * pi * s / 2) * arcint(g, depth)


def rvec_x(x, depth=11):
    f = lambda t: Delta(t) / (jay(t) - x)
    return [(2 * pi * I)**(N + 1) * (-1)**l * mp.binomial(N, l) * I**(l + 1)
            * Lstar(f, mp.mpf(l + 1), depth) for l in range(N + 1)], arcint(f, depth)


def align(u, v):
    """relative distance from u to the line C v"""
    nv = sum(abs(a)**2 for a in v)
    c = sum(a * mp.conj(b) for a, b in zip(u, v)) / nv
    return max(abs(a - c * b) for a, b in zip(u, v)) / max(abs(a) for a in u), c


print("guard: Phi(E_12) = %s" % mp.nstr(check_orientation(), 12), flush=True)
rE12 = rvec(E12, K)
rDel = rvec(Delta, K)

vs = {}
for x0 in (mp.mpf(500), mp.mpf(900)):
    for d in (mp.mpf(40), mp.mpf(15)):
        rp, Pp = rvec_x(x0 + I * d)
        rm, Pm = rvec_x(x0 - I * d)
        jr = [a - b for a, b in zip(rp, rm)]
        jP = Pp - Pm
        print("=" * 76, flush=True)
        print("x0=%-6s delta=%-5s  [Phi] = %-28s |[r_f]| = %s"
              % (mp.nstr(x0, 4), mp.nstr(d, 3), mp.nstr(jP, 10),
                 mp.nstr(max(abs(a) for a in jr), 10)), flush=True)
        if abs(jP) < mp.mpf('1e-20'):
            print("   [Phi] vanishes: the jump carries NO g-dependence, idea dead here.",
                  flush=True)
            continue
        v = [a / jP for a in jr]
        vs[(x0, d)] = v
        dE, cE = align(v, rE12)
        print("   v := [r_f]/[Phi]   |v| = %-22s dist to line C r_{E_12} = %s"
              % (mp.nstr(max(abs(a) for a in v), 12), mp.nstr(dE, 6)), flush=True)
        w = [a - b for a, b in zip(v, rE12)]
        dD, _ = align(w, rDel)
        print("   admissibility: dist of (v - r_{E_12}) from the line C r_Delta = %s"
              % mp.nstr(dD, 6), flush=True)

print("=" * 76, flush=True)
print("x-INDEPENDENCE of v (the decisive column):", flush=True)
ks = sorted(vs)
for a in range(len(ks)):
    for b in range(a + 1, len(ks)):
        u, w = vs[ks[a]], vs[ks[b]]
        rel = max(abs(p - q) for p, q in zip(u, w)) / max(abs(p) for p in u)
        tag = "same x0, delta drift" if ks[a][0] == ks[b][0] else "DIFFERENT x0"
        print("   %-22s vs %-22s  rel = %-14s  %s"
              % (str((mp.nstr(ks[a][0], 4), mp.nstr(ks[a][1], 3))),
                 str((mp.nstr(ks[b][0], 4), mp.nstr(ks[b][1], 3))),
                 mp.nstr(rel, 6), tag), flush=True)
