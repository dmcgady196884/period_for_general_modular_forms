"""Super-check of the residue polynomial: honest quadrature, no residue theorem.

residue_polynomial.py collapsed both segment integrals with the residue theorem and then did
exact algebra on Bernoulli polynomials.  That DERIVES the claim; it does not measure it.  Here
every integral is quadrature on an actual circle, so Lemma A (a_{gamma p} = (cp+d)^n a_p) and
Lemma B (K_T = regularised T-orbit sum) are tested in the form they are used, not assumed.

THE CONTOUR.  gamma_o(z) = c_z - S c_z with c_z a small positively oriented circle about z.
Since S is holomorphic, tau = S(w) = -1/w carries c_z to a positively oriented loop about Sz
with dtau = dw/w^2, so for any integrand F

    oint_{gamma_o} F(tau) dtau = oint_{c_z} [ F(w) - F(-1/w)/w^2 ] dw ,

one circle, no second quadrature.  The trapezoid rule is spectrally accurate on an analytic
periodic integrand, so M points on the circle is the right quadrature here.

WHAT IS TESTED
  1. simple pole: quadrature Xi vs the closed form
         Xi = 2 pi i a_z { (2 pi i)^{n+1}[ P_z|_{(1-S)} + K_T(z) - z^n K_T(-1/z) ] - (1-z^n) r_Ek }
     and W-membership of the QUADRATURE value.  Radius-independence is the control: a genuine
     contour integral must not move when the circle is resized.
  2. DOUBLE pole, f = E_4^2 (j')^2/(j - j(z))^2 in F_12.  The simple-pole algebra above does not
     apply -- residues now involve derivatives of both kernels -- but Definition~\\ref{def:respoly}
     is a contour integral and should not care.  If Xi leaves W here, the general claim is false.
  3. the ANTIsymmetric vs SYMmetric loop, by quadrature, at both pole orders.
  4. linearity: Xi_{f+g;z} = Xi_{f;z} + Xi_{g;z}, and two poles at once.

TRAP 4 does not arise (the circles are parametrised uniformly and the poles are at their centres),
but the radius must stay well inside the nearest other pole of the orbit; r = 0.02 with M = 256.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, E4, E6, Delta, jay, ktil, slash, rvec,      # noqa: E402
                    NS, NU, inW, arcint)

mp.mp.dps = 40
K, N = 12, 10
Sm = (0, -1, 1, 0)
E43 = lambda t: E4(t)**3
rE = rvec(E43, K, 11)
print("reference E_4^3   Phi = %s" % mp.nstr(arcint(E43, 11), 12), flush=True)

jp = lambda t: -2 * I * pi * jay(t) * E6(t) / E4(t)          # j' = -2 pi i j E_6/E_4


def loop(F, z, r, M, sign):
    """oint_{c_z + sign * S c_z} F(tau) dtau, by trapezoid on one circle."""
    tot = None
    for m in range(M):
        th = 2 * pi * mp.mpf(m) / M
        w = z + r * mp.e**(I * th)
        dw = 2 * I * pi * r * mp.e**(I * th) / M
        val = F(w)
        oth = F(-1 / w)
        if isinstance(val, list):
            term = [(a + sign * b / w**2) * dw for a, b in zip(val, oth)]
            tot = term if tot is None else [x + y for x, y in zip(tot, term)]
        else:
            term = (val + sign * oth / w**2) * dw
            tot = term if tot is None else tot + term
    return tot


def Pv(t):
    return [mp.binomial(N, l) * (-t)**l for l in range(N + 1)]


def Kv(t):
    return [mp.binomial(N, l) * (-1)**l * ktil(t, mp.mpf(l + 1), K) for l in range(N + 1)]


def Xi_quad(f, z, r=mp.mpf('0.02'), M=256, sign=-1):
    IP = loop(lambda t: [f(t) * c for c in Pv(t)], z, r, M, sign)
    IK = loop(lambda t: [f(t) * c for c in Kv(t)], z, r, M, sign)
    I1 = loop(f, z, r, M, sign)
    pre = (2 * pi * I)**(N + 1)
    return [pre * (a + b) - I1 * c for a, b, c in zip(IP, IK, rE)]


def Xi_alg(z, az, sign=-1):
    P = Pv(z)
    PS = slash(P, N, Sm)
    KA, KB = Kv(z), Kv(-1 / z)
    br = [a + sign * b + c + sign * z**N * d for a, b, c, d in zip(P, PS, KA, KB)]
    return [2 * pi * I * az * ((2 * pi * I)**(N + 1) * x - (1 + sign * z**N) * y)
            for x, y in zip(br, rE)]


def rep(tag, v):
    s_, u_ = inW(v, K) if max(abs(x) for x in v) > 0 else (mp.mpf(0), mp.mpf(0))
    print("   %-34s |(1+S)| = %-12s |(1+U+U^2)| = %-12s  |Xi| = %s"
          % (tag, mp.nstr(s_, 4), mp.nstr(u_, 4), mp.nstr(max(abs(x) for x in v), 8)),
          flush=True)
    return v


ZS = [("z = 0.23+1.41i", mp.mpf('0.23') + mp.mpf('1.41') * I),
      ("z = (1+i sqrt7)/2", (1 + I * mp.sqrt(7)) / 2),
      ("z = i sqrt2", I * mp.sqrt(2))]

print("\n" + "=" * 92, flush=True)
print("PART 1  SIMPLE pole,  f = E_4E_6 j'/(j - j(z))  in F_12", flush=True)
for lbl, z in ZS:
    f = (lambda t, z=z: E4(t) * E6(t) * jp(t) / (jay(t) - jay(z)))
    az = E4(z) * E6(z)
    print("\n%s    a_z = E_4E_6(z) = %s" % (lbl, mp.nstr(az, 10)), flush=True)
    q = rep("ANTI  quadrature", Xi_quad(f, z))
    a = rep("ANTI  closed form", Xi_alg(z, az))
    d = max(abs(x - y) for x, y in zip(q, a)) / max(abs(x) for x in a)
    print("   quadrature vs closed form: %s" % mp.nstr(d, 6), flush=True)
    q2 = Xi_quad(f, z, r=mp.mpf('0.045'), M=384)
    d2 = max(abs(x - y) for x, y in zip(q, q2)) / max(abs(x) for x in q)
    print("   radius 0.02 vs 0.045:      %s" % mp.nstr(d2, 6), flush=True)
    rep("SYM   quadrature", Xi_quad(f, z, sign=+1))

print("\n" + "=" * 92, flush=True)
print("PART 2  DOUBLE pole,  f = E_4^2 (j')^2/(j - j(z))^2  in F_12", flush=True)
print("        (the simple-pole algebra does NOT apply; the contour definition still does)",
      flush=True)
for lbl, z in ZS[:2]:
    f = (lambda t, z=z: E4(t)**2 * jp(t)**2 / (jay(t) - jay(z))**2)
    print("\n%s" % lbl, flush=True)
    q = rep("ANTI  quadrature", Xi_quad(f, z))
    q2 = Xi_quad(f, z, r=mp.mpf('0.045'), M=384)
    print("   radius 0.02 vs 0.045:      %s"
          % mp.nstr(max(abs(x - y) for x, y in zip(q, q2)) / max(abs(x) for x in q), 6),
          flush=True)
    rep("SYM   quadrature", Xi_quad(f, z, sign=+1))

print("\n" + "=" * 92, flush=True)
print("PART 3  linearity, simple + double at the same z", flush=True)
z = ZS[0][1]
f1 = lambda t: E4(t) * E6(t) * jp(t) / (jay(t) - jay(z))
f2 = lambda t: E4(t)**2 * jp(t)**2 / (jay(t) - jay(z))**2
s1, s2 = Xi_quad(f1, z), Xi_quad(f2, z)
s12 = Xi_quad(lambda t: f1(t) + 3 * f2(t), z)
tgt = [a + 3 * b for a, b in zip(s1, s2)]
print("   Xi_{f1 + 3 f2} vs Xi_{f1} + 3 Xi_{f2}: %s"
      % mp.nstr(max(abs(x - y) for x, y in zip(s12, tgt)) / max(abs(x) for x in tgt), 6),
      flush=True)
rep("sum in W?", s12)
