"""Audit of eq:symshift: is the SINGLE-pole winding in W, as Corollary~\\ref{cor:symspan} states?

cor:symspan's hypothesis is sound -- "S-symmetric" there means, per lem:georetrace's proof,
"S ghat^S is ghat^S reversed", i.e. S gamma = -gamma as a chain.  But eq:symshift then writes the
shift v_p using r_S(f,p), r_T(f,p) and a_p at the SINGLE point p, while calling it "one unit of
winding about the S-orbit of p".  If the pair must stay S-antisymmetric, winding +1 about p forces
-1 about Sp, so a term at Sp must appear.

Three windings, one script, same f and same quadrature:

    single   c_z              -> algebraically  frak a|_E
    sym      c_z + S c_z      -> frak a|_{(1+S)E}
    anti     c_z - S c_z      -> frak a|_{(1-S)E}

with E = 1 + calT(1-S).  Since E(1+S) = (1+S) + calT(1-S^2) = (1+S) and S^2 = 1 in PSL_2(Z):

    single  S-defect = frak a|_{(1+S)}        != 0
    sym     S-defect = frak a|_{(1+S)^2} = 2 frak a|_{(1+S)}  != 0
    anti    S-defect = frak a|_{(1-S)(1+S)} = 0

So the prediction is that ONLY anti lands in W, and that single and sym have S-defects in the exact
ratio 1 : 2.  That ratio is the sharp part: it is not a magnitude comparison but an identity, and it
discriminates a real algebraic mechanism from a numerical accident.

No residue theorem is used; every integral is trapezoid quadrature on one circle, with
oint_{S c_z} F dtau = oint_{c_z} F(-1/w) w^{-2} dw.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, E4, E6, jay, ktil, rvec, inW, arcint)       # noqa: E402

mp.mp.dps = 40
K, N = 12, 10
E43 = lambda t: E4(t)**3
rE = rvec(E43, K, 11)
jp = lambda t: -2 * I * pi * jay(t) * E6(t) / E4(t)
print("reference E_4^3  Phi = %s" % mp.nstr(arcint(E43, 11), 10), flush=True)


def Xi(f, z, sign, r=mp.mpf('0.02'), M=256):
    """sign = None -> c_z alone;  +1 -> c_z + S c_z;  -1 -> c_z - S c_z."""
    acc = [mp.mpc(0)] * (N + 1)
    tot = mp.mpc(0)
    for m in range(M):
        th = 2 * pi * mp.mpf(m) / M
        w = z + r * mp.e**(I * th)
        dw = 2 * I * pi * r * mp.e**(I * th) / M
        for pt, wt in ((w, mp.mpc(1)),) + (() if sign is None else ((-1 / w, sign / w**2),)):
            fv = f(pt) * wt * dw
            tot += fv
            for l in range(N + 1):
                acc[l] += fv * mp.binomial(N, l) * ((-pt)**l + (-1)**l
                                                    * ktil(pt, mp.mpf(l + 1), K))
    return [(2 * pi * I)**(N + 1) * a - tot * b for a, b in zip(acc, rE)]


for lbl, z in (("z = 0.23+1.41i", mp.mpf('0.23') + mp.mpf('1.41') * I),
               ("z = (1+i sqrt7)/2", (1 + I * mp.sqrt(7)) / 2)):
    f = (lambda t, z=z: E4(t) * E6(t) * jp(t) / (jay(t) - jay(z)))
    print("\n%s" % lbl, flush=True)
    keep = {}
    for tag, sg in (("single  c_z", None), ("sym     c_z + S c_z", +1),
                    ("anti    c_z - S c_z", -1)):
        v = Xi(f, z, sg)
        s_, u_ = inW(v, K)
        keep[tag[:4].strip()] = v
        print("   %-22s |(1+S)| = %-12s |(1+U+U^2)| = %-12s |Xi| = %s"
              % (tag, mp.nstr(s_, 4), mp.nstr(u_, 4),
                 mp.nstr(max(abs(x) for x in v), 8)), flush=True)
    dsym = inW(keep['sym'], K)[0] * max(abs(x) for x in keep['sym'])
    dsing = inW(keep['sing'], K)[0] * max(abs(x) for x in keep['sing'])
    print("   ratio of absolute S-defects  sym/single = %s   (prediction: exactly 2)"
          % mp.nstr(dsym / dsing, 12), flush=True)
