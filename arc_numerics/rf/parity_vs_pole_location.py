"""Item D: how does hat r_f sit inside W^+ (+) W^- as the pole location moves?

THE FORM.  DAM's choice, and it is the right one.  Since j has weight 0,
j' := dj/dtau = -2 pi i j E_6/E_4 has weight 2, and j'/(j - j(z)) = d/dtau log(j - j(z)) has a
SIMPLE pole at each point of SL2(Z).z with residue EXACTLY 1, independent of z: j - j(z) ~
j'(z)(tau - z) and the j'(z) cancels.  So f_z = g j'/(j - j(z)) with g in M_{k-2} has residue
g(z), smooth in z -- unlike Delta/(j-x), whose residue Delta(z)/j'(z) blows up at the elliptic
points where j' vanishes.  The z-dependence of the residue is thereby decoupled from the
z-dependence of the period polynomial, which is what one wants in order to watch the latter.

At k = 12 take g = E_10 = E_4 E_6 (M_10 is one-dimensional).  Then
    f_z = E_4 E_6 * (-2 pi i) (j E_6/E_4)/(j - x) = -2 pi i j E_6^2/(j - x),     x = j(z),
of weight 12, with residue E_10(z) at z.  NB c_{f}(0) = -2 pi i, NOT zero, since j/(j-x) -> 1 at
the cusp; so Phi(f) != 0 and hat r_f = r_f - Phi(f) r_{E_12} genuinely depends on the reference
form (ambiguity A.2, one-dimensional at k = 12).  Both r_f and hat r_f are therefore reported.

THE PREDICTION, and why it is a prediction and not a scan.  On the unit circle
-conj(tau) = e^{i(pi - theta)} = S tau, so the anti-holomorphic reflection tau -> -conj(tau) acts
on the arc exactly as S does.  For g with rational q-coefficients, f_{-conj z}(-conj tau) =
conj(f_z(tau)).  Hence the SELF-CONJUGATE pole locations -- z on the imaginary axis, where
-conj z = z, and z on Re z = +-1/2, where -conj z = z -+ 1 and so is T-equivalent -- are exactly
where the extra symmetry is present, and where hat r_f should be forced into a single parity
eigenspace.  Generic z should switch the other component on.

At k = 12, dim W = 3 = one odd direction plus two even, so the parity content is readable
directly off the coefficients: with n = 10, the coefficient of X^{n-l} Y^l is even in X iff l is
even.  Reported below are the largest even-l and odd-l coefficients and their ratio.

Also checked: the derivative identity.  f_z depends on z only through x = j(z), so
    d/dz hat r_{f_z} = j'(z) * hat r_{ g j'/(j - j(z))^2 },
the period polynomial of the DOUBLE-pole form.  A finite difference in z against the double-pole
computation tests that the z-motion stays inside the theory rather than needing new input.

All z used lie well inside the fundamental domain (Im z >= 1.2), so the poles are far above the
arc and the endpoint-clustered mesh of common.arcint is appropriate (cf. TRAP 4, which bites only
for poles ON or near the contour).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, E4, E6, jay, ktil, E12,                   # noqa: E402
                    arcint, rvec, inW, check_orientation)

mp.mp.dps = 30
K, N = 12, 10
rE12 = rvec(E12, K)


def jprime(t):
    return -2 * I * pi * jay(t) * E6(t) / E4(t)


def f_simple(x):
    return lambda t: -2 * I * pi * jay(t) * E6(t)**2 / (jay(t) - x)


def f_double(x):
    return lambda t: -2 * I * pi * jay(t) * E6(t)**2 / (jay(t) - x)**2


def rhat(f, depth=11):
    L = lambda s: mp.e**(-I * pi * s / 2) * arcint(
        lambda t: f(t) * (t**(s - 1) + ktil(t, s, K)), depth)
    rf = [(2 * pi * I)**(N + 1) * (-1)**l * mp.binomial(N, l) * I**(l + 1)
          * L(mp.mpf(l + 1)) for l in range(N + 1)]
    Ph = arcint(f, depth)
    return [a - Ph * b for a, b in zip(rf, rE12)], rf, Ph


def parity(v):
    ev = max(abs(v[l]) for l in range(0, N + 1, 2))
    od = max(abs(v[l]) for l in range(1, N + 1, 2))
    return ev, od


print("guard: Phi(E_12) = %s" % mp.nstr(check_orientation(), 12), flush=True)
print("f_z = -2 pi i j E_6^2/(j - j(z)),  k=12,  residue E_10(z) at z", flush=True)

ZS = [("imaginary axis  z = 1.2i", mp.mpf('1.2') * I),
      ("imaginary axis  z = 1.5i", mp.mpf('1.5') * I),
      ("imaginary axis  z = 2.0i", 2 * I),
      ("vertical edge   z = 0.5+1.2i", mp.mpf('0.5') + mp.mpf('1.2') * I),
      ("vertical edge   z = 0.5+1.6i", mp.mpf('0.5') + mp.mpf('1.6') * I),
      ("generic         z = 0.2+1.3i", mp.mpf('0.2') + mp.mpf('1.3') * I),
      ("generic         z = 0.35+1.5i", mp.mpf('0.35') + mp.mpf('1.5') * I)]

store = {}
for lbl, z in ZS:
    x = jay(z)
    v, rf, Ph = rhat(f_simple(x))
    store[lbl] = (z, x, v)
    s_, u_ = inW(v, K)
    eh, oh = parity(v)
    er, orr = parity(rf)
    print("=" * 78, flush=True)
    print("%s     x = j(z) = %s" % (lbl, mp.nstr(x, 10)), flush=True)
    print("   Phi = %-28s  |(1+S)| = %-11s |(1+U+U^2)| = %s"
          % (mp.nstr(Ph, 10), mp.nstr(s_, 5), mp.nstr(u_, 5)), flush=True)
    print("   hat r : max|even l| = %-22s max|odd l| = %-22s odd/even = %s"
          % (mp.nstr(eh, 10), mp.nstr(oh, 10), mp.nstr(oh / eh, 6)), flush=True)
    print("   r_f   : max|even l| = %-22s max|odd l| = %-22s odd/even = %s"
          % (mp.nstr(er, 10), mp.nstr(orr, 10), mp.nstr(orr / er, 6)), flush=True)

print("=" * 78, flush=True)
print("derivative identity  d/dz hat r_{f_z} = j'(z) hat r_{double}", flush=True)
z0 = mp.mpf('1.5') * I
h = mp.mpf('0.01')
vp, _, _ = rhat(f_simple(jay(z0 + h)))
vm, _, _ = rhat(f_simple(jay(z0 - h)))
fd = [(a - b) / (2 * h) for a, b in zip(vp, vm)]
vd, _, _ = rhat(f_double(jay(z0)))
pred = [jprime(z0) * a for a in vd]
rel = max(abs(a - b) for a, b in zip(fd, pred)) / max(abs(a) for a in fd)
print("   z0 = 1.5i, h = %s :  |finite difference - j'(z) hat r_double| / |fd| = %s"
      % (mp.nstr(h, 3), mp.nstr(rel, 8)), flush=True)
