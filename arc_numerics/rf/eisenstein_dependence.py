"""The rho limit is NOT intrinsic: it moves with the reference Eisenstein series.

thm:geoperiod subtracts Phi(f) r_{E_k}.  Both E_12 and E_4^3 have Phi = 1, and they
differ by a cusp form, E_4^3 = E_12 + (432000/691) Delta.  Crucially Phi(f_{j0})
diverges at the SAME j0^{-2/3} rate as tilde r, so the subtraction does not wash out
in the limit -- it shifts it.  With u = j0^{1/3},

    L  := lim u^2 tilde r ,     Pi := lim u^2 Phi(f_{j0}) ,
    L_43 = L_12 - (432000/691) Pi r_Delta .

Valid reference forms are exactly E_k + S_k (since Phi(E_k) = 1 and Phi|S_k = 0), and
h -> r_h is injective by Eichler-Shimura, so the limit determines hat r_f only modulo
r(S_k).  Projectively dim P(W) = 2 dim S_k and the reachable set has dimension
dim S_k: continuity HALVES the ambiguity, it does not kill it.  At k = 12, 2 -> 1.

Do NOT "check" this at dim S_k = 0: there dim W = 1, every nonzero element spans the
same line, and canonicality is automatic and vacuous.  The informative direction is
UPWARD, k = 24 with dim S_k = 2 and dim W = 5, where the prediction is 4 -> 2.

Also here: tilde r_Delta decomposed in the same basis.  Its own beta/alpha must come
out at exactly 36/691 -- the classical r^+_Delta = W_+ + (36/691) p_0 -- which
validates both the def:Wpm basis and the read-off.  v's ratio is NOT 36/691, so v's
even part is not proportional to Delta's even period polynomial.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common import (mp, I, pi, Delta, jay, check_orientation, arcint, rvec,
                    E12, E43, C691, inW, neville, decompose12)

mp.mp.dps = 40
KK, NN, DEPTH = 12, 10, 10
PSI = pi

print("Phi(E_12) =", mp.nstr(check_orientation(), 12), " (guard passed)", flush=True)

rE12, rE43, rD = rvec(E12, KK, DEPTH), rvec(E43, KK, DEPTH), rvec(Delta, KK, DEPTH)
chk = max(abs(a - b - C691 * c) for a, b, c in zip(rE43, rE12, rD))
print("identity r_E4^3 = r_E12 + (432000/691) r_Delta :  %s  (scale %s)"
      % (mp.nstr(chk, 5), mp.nstr(max(abs(x) for x in rE43), 8)))
print("Phi(E_4^3) =", mp.nstr(arcint(E43, DEPTH), 10), "  Phi(Delta) =",
      mp.nstr(arcint(Delta, DEPTH), 6), flush=True)

print()
print("=== Delta's own decomposition (validates the basis) ===")
aD, alD, beD = decompose12(rD)
print("  in W:", " / ".join(mp.nstr(x, 5) for x in inW(rD, KK)))
print("  alpha/a =", mp.nstr(alD / aD, 18))
print("  beta/a  =", mp.nstr(beD / aD, 18))
print("  beta/alpha =", mp.nstr(beD / alD, 18), "   vs 36/691 =",
      mp.nstr(mp.mpf(36) / 691, 18))

print()
print("=== the two limits ===")
us, L12s, L43s, Phis = [], [], [], []
for e in ('1e-2', '1e-3', '1e-4', '1e-5', '1e-6'):
    j0 = mp.mpf(e) * mp.e**(I * PSI)
    f = lambda t, j0=j0: Delta(t) / (jay(t) - j0)
    Phi = arcint(f, DEPTH)
    rf = rvec(f, KK, DEPTH)
    u = j0**(mp.mpf(1) / 3)
    us.append(u)
    L12s.append([u**2 * (a - Phi * b) for a, b in zip(rf, rE12)])
    L43s.append([u**2 * (a - Phi * b) for a, b in zip(rf, rE43)])
    Phis.append(u**2 * Phi)
    print("  |j0|=%-5s |Phi|=%-14s" % (e, mp.nstr(abs(Phi), 8)), flush=True)

L12 = [neville(us, [L12s[m][l] for m in range(len(us))])[0] for l in range(NN + 1)]
L43 = [neville(us, [L43s[m][l] for m in range(len(us))])[0] for l in range(NN + 1)]
Pi, shPi = neville(us, Phis)
print("  Pi =", mp.nstr(Pi, 14), " [shift", mp.nstr(shPi, 3), "]")

n12 = [x / L12[5] for x in L12]
n43 = [x / L43[5] for x in L43]
print()
print("=== DOES THE LINE MOVE? ===")
print("  max|dir(E4^3) - dir(E12)| =",
      mp.nstr(max(abs(a - b) for a, b in zip(n12, n43)), 10))
print("  both in W:  E12", " / ".join(mp.nstr(x, 4) for x in inW(L12, KK)),
      "   E4^3", " / ".join(mp.nstr(x, 4) for x in inW(L43, KK)))

pred = [a - C691 * Pi * b for a, b in zip(L12, rD)]
print()
print("=== quantitative:  L43 = L12 - (432000/691) Pi r_Delta ? ===")
print("  relative residual =",
      mp.nstr(max(abs(a - b) for a, b in zip(L43, pred))
              / max(abs(x) for x in L43), 8))

print()
for nm, L in (("E12 ", L12), ("E4^3", L43)):
    a_, al_, be_ = decompose12(L)
    print("  %s  alpha/a = %-30s beta/a = %s"
          % (nm, mp.nstr(al_ / a_, 18), mp.nstr(be_ / a_, 18)))
