"""Does the h+/h- kernel ambiguity see INTERIOR poles?  (P5 of hurwitz_kernel_handoff_20260923.md)

The two one-sided inverses of (1-T) differ by the two-sided Lipschitz sum,
    h+[w] - h-[w] = (-2 pi i)^{-w} Gamma(-w)^{-1} Li_{w+1}(q) ,
so with C(a) := (-2 pi i)^{1-a}/Gamma(1-a),
    ktil_+ - ktil_- = C(s) Li_s(q) - e^{i pi (s-1)} C(k-s) Li_{k-s}(q)   =: B(tau,s),
1-periodic and holomorphic on H.  The S-segment integrand is identical for both kernels, so
    L*_+ - L*_- = e^{-i pi s/2} int_{gamma^T} f B dtau .

PREDICTION, derived before running.  Push the T-segment (the arc, left to right) up past every
pole.  Above all poles f = sum c(n) q^n and the integral is the constant term of f*B, which is
    CT(f Li_s(q)) = sum_{m>=1} c(-m) m^{-s}     (Li_s has no constant term).
Crossing a simple pole z in the strip adds +2 pi i Res_z(f B) = 2 pi i a_z B(z,s): arc minus high
segment bounds the region counterclockwise -- the same orientation that gives the measured
Phi(f_z) = +2 pi i g(z).  Hence

  (1) f_z = Delta E_4 J'/(J - J(z)), k = 18: no principal part at the cusp, so the difference is
      e^{-i pi s/2} 2 pi i g(z) B(z,s).  NONZERO at non-integer s.
  (2) E_4^3 E_6 J, k = 18: cusp principal part q^{-1}, no interior pole; difference
      e^{-i pi s/2} [C(s) - e^{i pi (s-1)} C(k-s)].
  (3) their sum: the sum.
  At s = 6 every prediction is 0, since 1/Gamma(1-s) = 1/Gamma(1-(k-s)) = 0.

If (1) holds, the ambiguity is carried by the principal part ANYWHERE -- cusp or interior -- and
Appendix A must say so.  The pole at z sits well above the arc, so plain arc quadrature applies
(TRAP 4 does not arise).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, E4, E6, Delta, jay, arcint                   # noqa: E402

mp.mp.dps = 30
K = 18
DEPTH = 11

C = lambda a: (-2 * pi * I)**(1 - a) * mp.rgamma(1 - a)
B = lambda t, s: (C(s) * mp.polylog(s, mp.e**(2 * I * pi * t))
                  - mp.e**(I * pi * (s - 1)) * C(K - s) * mp.polylog(K - s, mp.e**(2 * I * pi * t)))


def hplus(w, t):
    return mp.zeta(-w, t + 1)


def hminus(w, t):
    return -mp.e**(I * pi * w) * mp.zeta(-w, -t)


def kern(t, s, h):
    return h(s - 1, t) - mp.e**(I * pi * (s - 1)) * h(K - 1 - s, t)


print("check: h+ - h- = B on the arc, at one point", flush=True)
t0 = mp.e**(I * mp.mpf('1.9'))
for s in (mp.mpf('3.3'), mp.mpc(5, 1.5)):
    d = kern(t0, s, hplus) - kern(t0, s, hminus)
    print("   s=%-10s |(k+ - k-) - B| / |B| = %s"
          % (mp.nstr(s, 4), mp.nstr(abs(d - B(t0, s)) / abs(B(t0, s)), 4)), flush=True)

jp = lambda t: -2 * I * pi * jay(t) * E6(t) / E4(t)
g = lambda t: Delta(t) * E4(t)
cusp = lambda t: E4(t)**3 * E6(t) * (jay(t) - 744)          # E_4^3 E_6 J = q^{-1} + 216 + ...

for lbl, z in (("z=(1+i sqrt7)/2", (1 + I * mp.sqrt(7)) / 2),
               ("z=0.2+1.3i", mp.mpf('0.2') + mp.mpf('1.3') * I)):
    fz = (lambda t, z=z: g(t) * jp(t) / (jay(t) - jay(z)))
    az = g(z)
    print("\n%s   a_z = g(z) = %s" % (lbl, mp.nstr(az, 10)), flush=True)
    for s in (mp.mpf('3.3'), mp.mpc(5, 1.5), mp.mpf(6)):
        pre = mp.e**(-I * pi * s / 2)
        for tag, F, pred in (
                ("(1) interior pole only", fz, pre * 2 * pi * I * az * B(z, s)),
                ("(2) cusp principal only", cusp, pre * (C(s) - mp.e**(I * pi * (s - 1)) * C(K - s))),
                ("(3) sum", lambda t, fz=fz: fz(t) + cusp(t),
                 pre * (2 * pi * I * az * B(z, s) + C(s) - mp.e**(I * pi * (s - 1)) * C(K - s)))):
            got = pre * arcint(lambda t, F=F, s=s: F(t) * B(t, s), DEPTH)
            if abs(pred) < mp.mpf('1e-25'):
                err = "abs %s (prediction is exactly 0)" % mp.nstr(abs(got), 4)
            else:
                err = "rel.err %s" % mp.nstr(abs(got - pred) / abs(pred), 4)
            print("   s=%-10s %-26s L+ - L- = %-34s %s"
                  % (mp.nstr(s, 4), tag, mp.nstr(got, 12), err), flush=True)
