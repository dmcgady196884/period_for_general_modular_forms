#!/usr/bin/env python3
"""
verify_kernel.py -- reproduces every numerical claim in
  principal_part_ambiguity.tex  and  hurwitz_kernel_handoff.md
(project "Hecke and MFs with poles").

Usage:  pip install mpmath ;  python3 verify_kernel.py [--slow]
  Runtime: ~5 min; --slow adds the horizontal-ray unfolding check (Lemma 3), ~1-2 min more.
  Expected: "47/47 checks passed" with --slow (45/45 without).

Conventions (match the project draft): weight k even, tau0 = i (the S-segment degenerates),
principal branches for tau^w on C \\ (-inf, 0].
  h+[w](t) =  sum_{n>=1} (t+n)^w = zeta(-w, t+1)             (the draft's Hurwitz kernel)
  h-[w](t) = -sum_{n<=0} (t+n)^w = -e^{i pi w} zeta(-w, -t)
  K(t,s)   = hA[s-1](t) - e^{i pi (s-1)} hB[k-1-s](t)          (hA = hB = h+ : draft's k~_T)
  L*(f,s)  = e^{-pi i s/2} int_{i-1}^{i} f(t) K(t,s) dt
Each check prints PASS/FAIL against a stated tolerance.
"""
import sys
from fractions import Fraction as Fr
from mpmath import (mp, mpf, mpc, zeta, exp, pi, gamma, rgamma, polylog, quad,
                    gammainc, bernoulli, ei, conj, im, nsum, factorial, inf, cos)

mp.dps = 30
I = mpc(0, 1)
RESULTS = []

def check(name, err, tol):
    ok = err <= tol
    RESULTS.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {mp.nstr(err, 3)} (tol {mp.nstr(tol, 1)})")

# ---------------------------------------------------------------- kernels
def h_plus(w, t):  return zeta(-w, t + 1)
def h_minus(w, t): return -exp(I*pi*w) * zeta(-w, -t)
def Lstar(c, s, k, hA=h_plus, hB=h_plus):
    f = lambda t: sum(v * exp(2*pi*I*n*t) for n, v in c.items())
    K = lambda t: hA(s - 1, t) - exp(I*pi*(s - 1)) * hB(k - 1 - s, t)
    return exp(-pi*I*s/2) * quad(lambda x: f(x + I) * K(x + I), [-1, -0.5, 0])

def correction(c, s, k):  # Theorem: L*_+ - L*_- = -2 pi i [ D(s)/((2pi)^s G(1-s)) + i^k D(k-s)/(...) ]
    D = lambda a: sum(v * mpf(-n)**(-a) for n, v in c.items() if n < 0)
    return -2*pi*I * ((2*pi)**(-s) * D(s) * rgamma(1 - s)
                      + I**k * (2*pi)**(-(k - s)) * D(k - s) * rgamma(1 - (k - s)))

# ---------------------------------------------------------------- q-series (weight 12 test forms)
N = 40
def sig(n, e): return sum(d**e for d in range(1, n + 1) if n % d == 0)
def mul(a, b):
    out = [0] * N
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b[:N - i]): out[i + j] += x * y
    return out
E4 = [1] + [240 * sig(n, 3) for n in range(1, N)]
E6 = [1] + [-504 * sig(n, 5) for n in range(1, N)]
P = [1] + [0] * (N - 1)                                   # Delta / q
for n in range(1, N):
    fac = [0] * N; fac[0] = 1; fac[n] = -1
    for _ in range(24): P = mul(P, fac)
inv = [0] * N; inv[0] = 1
for n in range(1, N): inv[n] = -sum(P[j] * inv[n - j] for j in range(1, n + 1))
E43, E62 = mul(E4, mul(E4, E4)), mul(E6, E6)
A = mul(mul(E62, E43), inv)                               # q * E6^2 E4^3 / Delta
J = mul(E43, inv)                                         # q * j
# Brown's weak cusp form: Delta' = E6^2E4^3/Delta + (5541/16)E4^3 - (1317/16)E6^2   (email text had 5514: typo)
DPq = {n: Fr(A[n + 1]) + (Fr(5541, 16) * E43[n] - Fr(1317, 16) * E62[n] if n >= 0 else 0) for n in range(-1, N - 2)}
assert DPq[-1] == 1 and DPq[0] == 0 and DPq[1] == 0 and DPq[2] == 47709536
# g = j*Delta' - 744*Delta' = q^-2 + 196884 + 69203296 q + ...   (tests m-dependence: D(g,s) = 2^-s)
gq = {}
for n in range(-2, N - 3):
    gq[n] = sum(Fr(J[a + 1]) * DPq[n - a] for a in range(-1, N - 2) if (n - a) in DPq and a + 1 < N) - 744 * DPq.get(n, 0)
assert gq[-2] == 1 and gq[-1] == 0 and gq[0] == 196884 and gq[1] == 69203296
tomp = lambda c: {n: mpf(v.numerator) / v.denominator for n, v in c.items()}
DP, G = tomp(DPq), tomp(gq)
k = 12
DPNAME = "Delta'"

if __name__ == "__main__":
    print("== 1. Lipschitz identity  h+ - h- = (-2 pi i)^{1-s}/Gamma(1-s) Li_s(q)")
    t = mpc(0.3, 0.8)
    for s in [mpc(-1.3, 0.7), mpc(0.4, 2), mpf(3.7)]:
        d = h_plus(s - 1, t) - h_minus(s - 1, t)
        check(f"s={s}", abs(d - (-2*pi*I)**(1 - s) * rgamma(1 - s) * polylog(s, exp(2*pi*I*t))), mpf(10)**-24)

    print("== 2. On M_k the kernel gives the classical L-function (2pi)^-s Gamma(s) zeta(s) zeta(s-k+1)")
    for kk in [4, 12]:
        cG = {0: -bernoulli(kk) / (2*kk)}; cG.update({n: mpf(sig(n, kk - 1)) for n in range(1, N)})
        for s in [mpf('2.5'), mpc(1.7, 0.9)]:
            check(f"G_{kk}, s={s}", abs(Lstar(cG, s, kk) - (2*pi)**(-s) * gamma(s) * zeta(s) * zeta(s - kk + 1)), mpf(10)**-24)

    print("== 3. Theorem: L*_+ - L*_- = principal-part Dirichlet polynomial (residue at q=0)")
    vals = {}
    for name, c in [("Delta'", DP), ("g", G)]:
        for s in [mpf('3.3'), mpc(5, 1.5), mpf('0.37')]:
            Lp, Lm = Lstar(c, s, k, h_plus, h_plus), Lstar(c, s, k, h_minus, h_minus)
            vals[(name, s)] = (Lp, Lm)
            check(f"{name}, s={s}", abs(Lp - Lm - correction(c, s, k)), mpf(10)**-20)
    print("   sample: L*_+(Delta',3.3) =", mp.nstr(vals[(DPNAME, mpf('3.3'))][0], 16))

    print("== 4. Correction vanishes at critical integers s = 1..k-1 (so it never touches period polynomials)")
    for name, c in [("Delta'", DP), ("g", G)]:
        worst = max(abs(correction(c, mpf(s), k)) for s in range(1, k))
        check(f"{name}: max |closed-form correction| over s=1..11", worst, mpf(10)**-25)
        for s in [3, 6, 9]:
            check(f"{name}: |L*_+ - L*_-| at s={s} (direct)", abs(Lstar(c, mpf(s), k) - Lstar(c, mpf(s), k, h_minus, h_minus)), mpf(10)**-20)
    print(f"   (but nonzero at s=0: correction for Delta' ~ {mp.nstr(correction(DP, mpf('1e-12'), k), 6)})")

    print("== 5. Proposition 9: global difference = sum of per-mode Gamma(s,z) branch jumps")
    eps = mpf(10)**-30
    lat = lambda a, x, side: gammainc(a, mpc(-x, side*eps)) / mpc(-x, side*eps)**a
    for name, c in [("Delta'", DP), ("g", G)]:
        for s in [mpf('3.3'), mpc(5, 1.5), mpf('0.37')]:
            jump = sum(v * ((lat(s, 2*pi*(-n), 1) + I**k * lat(k - s, 2*pi*(-n), 1))
                          - (lat(s, 2*pi*(-n), -1) + I**k * lat(k - s, 2*pi*(-n), -1)))
                       for n, v in c.items() if n < 0)
            Lp, Lm = vals[(name, s)]
            check(f"{name}, s={s}: (L+ - L-) vs sum of jumps (upper - lower)", abs(Lp - Lm - jump), mpf(10)**-20)
            print(f"      control, wrong side: off by {mp.nstr(abs(Lp - Lm + jump), 3)}")

    print("== 6. h+ = BFK series with principal branch (upper side of the cut)")
    for s in [mpf('3.3'), mpc(5, 1.5)]:
        bfk = sum(v * (gammainc(s, 2*pi*n) / (2*pi*n)**s + I**k * gammainc(k - s, 2*pi*n) / (2*pi*n)**(k - s))
                  for n, v in DP.items() if n != 0)
        check(f"Delta', s={s}", abs(vals[(DPNAME, s)][0] - bfk), mpf(10)**-20)

    print("== 7. Reality: at real s, L*_- = conj(L*_+); mean is the principal value")
    for name in ["Delta'", "g"]:
        for s in [mpf('3.3'), mpf('0.37')]:
            Lp, Lm = vals[(name, s)]
            check(f"{name}, s={s}: |L- - conj(L+)|", abs(Lm - conj(Lp)), mpf(10)**-20)
            print(f"      Im L*_+ = {mp.nstr(im(Lp), 6)}  (nonzero: principal-branch value is not real)")
    x = 2*pi
    check("mean of sides of Gamma(0,-2pi) vs -Ei(2pi)",
          abs((gammainc(0, mpc(-x, eps)) + gammainc(0, mpc(-x, -eps))) / 2 + ei(x)), mpf(10)**-25)

    print("== 8. Branch-explicit symmetric formula (L*_sym = mean of L*_+ and L*_-)")
    E = lambda a, m: nsum(lambda j: (2*pi*m)**j / (factorial(j) * (a + j)), [0, inf])
    for s in [mpf('3.3'), mpc(5, 1.5), mpf('-0.7')]:
        pos = sum(v * (gammainc(s, 2*pi*n) / (2*pi*n)**s + I**k * gammainc(k - s, 2*pi*n) / (2*pi*n)**(k - s))
                  for n, v in DP.items() if n > 0)
        neg = sum(v * (cos(pi*s) * (gamma(s) * (2*pi*(-n))**(-s) + I**k * gamma(k - s) * (2*pi*(-n))**(-(k - s)))
                       - E(s, -n) - I**k * E(k - s, -n)) for n, v in DP.items() if n < 0)
        Lp = Lstar(DP, s, k); Lm = Lstar(DP, s, k, h_minus, h_minus)
        check(f"Delta', s={s}", abs(pos + neg - (Lp + Lm) / 2), mpf(10)**-20)

    print("== 9. Functional equation L*(s) = i^k L*(k-s): holds for K+, K-; fails for mixed sides")
    for name, c in [("Delta'", DP), ("g", G)]:
        for s in [mpf('3.3'), mpc(5, 1.5)]:
            for lab, hA, hB in [("K+", h_plus, h_plus), ("K-", h_minus, h_minus)]:
                check(f"{name}, s={s}, {lab}", abs(Lstar(c, s, k, hA, hB) - I**k * Lstar(c, k - s, k, hA, hB)), mpf(10)**-18)
            mix = abs(Lstar(c, s, k, h_plus, h_minus) - I**k * Lstar(c, k - s, k, h_plus, h_minus))
            print(f"      mixed kernel (h+ | h-): FE violated by {mp.nstr(mix, 3)}")

    if "--slow" in sys.argv:
        print("== 10. Lemma 3: segment pairing = horizontal Mellin ray (s=-1.5, 400 periods; truncation ~1e-5)")
        mp.dps = 15
        s = mpf('-1.5'); w = s - 1; X = 400
        f = lambda t: sum(v * exp(2*pi*I*n*t) for n, v in DP.items())
        seg_p = quad(lambda x: f(x + I) * h_plus(w, x + I), [-1, 0])
        seg_m = quad(lambda x: f(x + I) * h_minus(w, x + I), [-1, 0])
        ray_r = sum(quad(lambda x: f(x + I) * (x + I)**w, [j, j + 1]) for j in range(X))
        ray_l = -sum(quad(lambda x: f(x + I) * (x + I)**w, [-j - 1, -j]) for j in range(X))
        check("right ray vs h+ segment", abs(seg_p - ray_r), mpf(10)**-4)
        check("left ray vs h- segment", abs(seg_m - ray_l), mpf(10)**-4)

    print(f"\n{sum(RESULTS)}/{len(RESULTS)} checks passed")
