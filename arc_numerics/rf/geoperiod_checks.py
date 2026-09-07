"""thm:geoperiod -- tilde r_f = r_f - Phi(f) r_{E_k} is in W for every f with no pole
on SL2(Z).gamma^arc, with NO condition on where the poles are or how deep.

Two independent checks.

(1) The defect factors through Phi.  If the (1+U+U^2) loop encloses nothing, the map
    f -> r_f|(1+U+U^2) kills a codimension-1 subspace and hence factors through the
    single functional Phi(f) = int_arc f dtau:

        r_f|(1+U+U^2) = Phi(f) . D ,   D a SINGLE universal vector in V_n.

    That is far stronger than "one tuned example gave zero" and fails loudly if the
    cancellation were a coincidence.  D must equal r_{E_k}|(1+U+U^2), which is what
    makes the subtraction in thm:geoperiod work.  Test forms are deliberately
    unrelated: poles near and far, a double pole, two holomorphic forms, and Delta
    (Phi = 0, so its defect must vanish on its own).

(2) tilde r_f lands in W for each of them.

Trap this guards against: condition 2 of eq:Fcirc is int_{gamma^T} f dtau = 0, which
equals c_f(0) = 0 ONLY when no pole sits above the contour.  The arc dips to sqrt3/2,
so it fails generically even for forms with c_f(0) = 0, and testing an untuned
Delta/(j - j0) on the arc gives defects ~0.4 that look like a refutation and are not.
Delta is useless as a control here: being holomorphic it satisfies condition 2 on
every contour, so it validates the harness while being blind to the bug.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common import (mp, I, pi, E4, E6, Delta, jay, check_orientation,
                    arcint, rvec, NU, inW, E12)

mp.mp.dps = 25
KK = 12

print("Phi(E_12) =", mp.nstr(check_orientation(), 12), " (guard passed)\n", flush=True)

PHI = 7 * pi / 12
z1 = mp.mpf('1.02') * mp.e**(I * PHI)      # pole orbit close to the arc
z2 = mp.mpf('1.30') * mp.e**(I * PHI)      # further out
z3 = mp.mpc('0.1', '1.6')                  # high up
j1, j2, j3 = jay(z1), jay(z2), jay(z3)

FORMS = [
    ("A  Delta/(j-j(z1))   near arc ", lambda t: Delta(t) / (jay(t) - j1)),
    ("B  Delta/(j-j(z2))   further  ", lambda t: Delta(t) / (jay(t) - j2)),
    ("C  Delta/(j-j(z3))   far up   ", lambda t: Delta(t) / (jay(t) - j3)),
    ("D  Delta/(j-j(z3))^2 DOUBLE   ", lambda t: Delta(t) / (jay(t) - j3)**2),
    ("E  E4^3              holo     ", lambda t: E4(t)**3),
    ("F  E6^2              holo     ", lambda t: E6(t)**2),
    ("G  Delta             Phi = 0  ", Delta),
]

rE = rvec(E12, KK)
res = []
for name, g in FORMS:
    ph = arcint(g)
    v = rvec(g, KK)
    d = NU(v, KK)
    tr = [a - ph * b for a, b in zip(v, rE)]
    res.append((name, ph, v, d, tr))
    s_, u_ = inW(tr, KK)
    print("%s Phi = %-22s |defect|/|r| = %-12s tilde r in W: %s / %s"
          % (name, mp.nstr(ph, 8),
             mp.nstr(max(abs(x) for x in d) / max(abs(x) for x in v), 5),
             mp.nstr(s_, 4), mp.nstr(u_, 4)), flush=True)

print()
print("=== is defect/Phi the SAME vector D for every f? ===")
ref = None
for name, ph, v, d, tr in res:
    if abs(ph) < mp.mpf('1e-15'):
        print("   %s Phi ~ 0, defect must vanish alone: %s"
              % (name, mp.nstr(max(abs(x) for x in d) / max(abs(x) for x in v), 5)))
        continue
    D = [x / ph for x in d]
    if ref is None:
        ref = D
        print("   %s <- reference D" % name)
    else:
        print("   %s max|D_f - D_ref| / |D_ref| = %s"
              % (name, mp.nstr(max(abs(a - b) for a, b in zip(D, ref))
                               / max(abs(x) for x in ref), 5)))

print()
print("=== and D must be r_{E_12}|(1+U+U^2) ===")
dE = NU(rE, KK)
print("   max|D_ref - r_E12|(1+U+U^2)| / |.| =",
      mp.nstr(max(abs(a - b) for a, b in zip(ref, dE)) / max(abs(x) for x in dE), 6))
