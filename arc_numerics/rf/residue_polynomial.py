"""Is the winding about ONE pole enough to land in W?  Residue form, no quadrature.

For a simple pole at z with residue a_z, the residue theorem collapses both segment integrals
of Proposition~\\ref{prop:periodint}:

    oint_{c_z} f P_tau dtau = 2 pi i a_z P_z,      P_z := (X - zY)^n ,
    oint_{c_z} f K_T dtau   = 2 pi i a_z K_T(z) .

Modularity gives a_{gamma p} = (cp+d)^n a_p, so R_p := a_p P_p obeys R_{gamma p} = R_p|_{gamma^-1}:
the residue-weighted S-kernel is EQUIVARIANT.  In particular a_{Sz} = z^n a_z and
z^n P_{Sz} = P_z|_S.  Hence, with a_z = 1, everything below is a universal function of z alone --
no f is needed, and the claim "one pole's winding lands in W" is a claim about these polynomials.

TWO CANDIDATE LOOPS.  Lemma~\\ref{lem:georetrace} needs S gamma^S = -gamma^S for the S^2- and
(ST)^3-paths to retrace, so a winding added to both segments must be S-ANTIsymmetric:

    sym  (the draft's "positively oriented loop enclosing {z,Sz}")   gamma_o = c_z + S c_z
    anti (what lem:georetrace actually demands)                      gamma_o = c_z - S c_z

PREDICTION, from K_T(tau) = sum_{m>=1} P_{tau+m}|_{(1-S)} (regularised, the only input being the
Hurwitz recursion T calT = 1 + calT) and E := 1 + calT(1-S):

    sym  -> R_z|_{(1+S)E},   E(1+S) = (1+S) + calT(1-S^2) = (1+S),  so the (1+S) defect is
            R_z|_{(1+S)^2} = 2 R_z|_{(1+S)} != 0                     ==> NOT in W
    anti -> R_z|_{(1-S)E},   (1-S)E(1+S) = (1-S)(1+S) = 1 - S^2 = 0  ==> (1+S) defect vanishes

S^2 = 1 holds in PSL_2(Z) and n is even, so this is exact, not asymptotic.  The U-relation is the
open half: it should need the r_{E_k} counterterm, which is what Phi measures.  PART 3 tests that
separately by reporting the U-defect of the raw bracket against 2 pi i (1 -+ z^n) r_{E_k}|_{(1+U+U^2)}.

At z = i the antisymmetric loop DEGENERATES: Si = i, so c_i - S c_i = 0 and Xi_{f;i} = 0
identically.  Consistently, a_{Si} = i^n a_i forces a_i = 0 unless 4 | n, i.e. unless k = 2 mod 4.
At z = rho the orbit point and its S-image are rho and S rho = rho+1, the two ARC ENDPOINTS.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, E4, E6, Delta, jay, ktil, slash, rvec,      # noqa: E402
                    NS, NU, inW, arcint)

mp.mp.dps = 40
Sm = (0, -1, 1, 0)


def Pz(z, n):
    return [mp.binomial(n, l) * (-z)**l for l in range(n + 1)]


def KTz(z, n, k):
    return [mp.binomial(n, l) * (-1)**l * ktil(z, mp.mpf(l + 1), k) for l in range(n + 1)]


def add(*vs):
    return [sum(t) for t in zip(*vs)]


def sc(c, v):
    return [c * x for x in v]


def rel(v, k):
    return inW(v, k) if max(abs(x) for x in v) > 0 else (mp.mpf(0), mp.mpf(0))


# ---------------------------------------------------------------- PART 1: conventions
print("PART 1  kernel identity  K_T(tau) - K_T(tau-1) = -(P_tau - P_tau|_S)", flush=True)
for k in (12, 18):
    n = k - 2
    t = mp.mpf('0.31') + mp.mpf('1.27') * I
    lhs = [a - b for a, b in zip(KTz(t, n, k), KTz(t - 1, n, k))]
    P = Pz(t, n)
    rhs = [-(a - b) for a, b in zip(P, slash(P, n, Sm))]
    err = max(abs(a - b) for a, b in zip(lhs, rhs)) / max(abs(x) for x in rhs)
    print("   k=%-3d relative error %s" % (k, mp.nstr(err, 6)), flush=True)

print("\nPART 1b  residue equivariance  a_{Sz} = z^n a_z, on f = E_4E_6 j'/(j - j(z)), k=12",
      flush=True)
z0 = mp.mpf('0.23') + mp.mpf('1.41') * I
f0 = lambda t: (E4(t) * E6(t) * (-2 * I * pi * jay(t) * E6(t) / E4(t)) /
                (jay(t) - jay(z0)))


def res_at(f, p, r, M=64):
    return sum(f(p + r * mp.e**(2 * I * pi * mp.mpf(m) / M)) * r * mp.e**(2 * I * pi * mp.mpf(m) / M)
               for m in range(M)) / M


az = res_at(f0, z0, mp.mpf('0.05'))
aSz = res_at(f0, -1 / z0, mp.mpf('0.05') / abs(z0)**2)
print("   a_z          = %s" % mp.nstr(az, 12), flush=True)
print("   a_Sz         = %s" % mp.nstr(aSz, 12), flush=True)
print("   z^n a_z      = %s" % mp.nstr(z0**10 * az, 12), flush=True)
print("   g(z)=E4E6(z) = %s   (predicted a_z)" % mp.nstr(E4(z0) * E6(z0), 12), flush=True)
print("   relative error a_Sz vs z^n a_z: %s"
      % mp.nstr(abs(aSz - z0**10 * az) / abs(aSz), 6), flush=True)

# ---------------------------------------------------------------- PART 2 & 3
REF = {12: ('E_4^3', lambda t: E4(t)**3),
       14: ('E_4^2 E_6', lambda t: E4(t)**2 * E6(t)),
       18: ('E_4^3 E_6', lambda t: E4(t)**3 * E6(t))}

s3 = mp.sqrt(3)
ZS = [("d=3   rho      (elliptic, n_z=3)", -mp.mpf(1) / 2 + s3 / 2 * I),
      ("d=4   i        (elliptic, n_z=2)", I),
      ("d=7   (1+isq7)/2 ", (1 + I * mp.sqrt(7)) / 2),
      ("d=8   i sqrt2    ", I * mp.sqrt(2)),
      ("d=11  (1+isq11)/2", (1 + I * mp.sqrt(11)) / 2),
      ("d=19  (1+isq19)/2", (1 + I * mp.sqrt(19)) / 2),
      ("d=43  (1+isq43)/2", (1 + I * mp.sqrt(43)) / 2),
      ("d=67  (1+isq67)/2", (1 + I * mp.sqrt(67)) / 2),
      ("d=163 (1+isq163)/2", (1 + I * mp.sqrt(163)) / 2),
      ("generic 0.2+1.3i ", mp.mpf('0.2') + mp.mpf('1.3') * I),
      ("generic 0.35+1.5i", mp.mpf('0.35') + mp.mpf('1.5') * I)]

for k in (12, 18):
    n = k - 2
    nm, g = REF[k]
    rE = rvec(g, k, 11)
    PhE = arcint(g, 11)
    print("\n" + "=" * 96, flush=True)
    print("k = %d   n = %d   reference %s   Phi = %s   r_E in W? %s / %s"
          % (k, n, nm, mp.nstr(PhE, 10),
             mp.nstr(rel(rE, k)[0], 4), mp.nstr(rel(rE, k)[1], 4)), flush=True)
    print("%-22s %-24s %-24s %s"
          % ("z", "SYM   (1+S) / (1+U+U^2)", "ANTI  (1+S) / (1+U+U^2)",
             "ANTI raw U-defect vs 2pi i(1-z^n) r_E|_U"), flush=True)
    for lbl, z in ZS:
        P, PS = Pz(z, n), slash(Pz(z, n), n, Sm)
        KA, KB = KTz(z, n, k), KTz(-1 / z, n, k)
        pre = (2 * pi * I)**(n + 2)
        out = []
        for sg in (+1, -1):
            raw = sc(pre, add(P, sc(sg, PS), KA, sc(sg * z**n, KB)))
            cor = sc(2 * pi * I * (1 + sg * z**n), rE)
            th = [a - b for a, b in zip(raw, cor)]
            out.append((th, raw, cor))
        (ths, _, _), (tha, rawa, cora) = out
        ss, su = rel(ths, k)
        as_, au = rel(tha, k)
        ru = NU(rawa, k)
        cu = NU(cora, k)
        mr = max(abs(x) for x in ru)
        match = (max(abs(a - b) for a, b in zip(ru, cu)) / mr) if mr > 0 else mp.mpf(0)
        print("%-22s %-11s %-11s %-11s %-11s %s"
              % (lbl, mp.nstr(ss, 4), mp.nstr(su, 4), mp.nstr(as_, 4), mp.nstr(au, 4),
                 mp.nstr(match, 6)), flush=True)
