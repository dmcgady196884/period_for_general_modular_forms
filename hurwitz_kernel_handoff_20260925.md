# Hurwitz T-kernel: the branch choice in L-functions of weak forms, read as a choice of (1−T)⁻¹

Handoff notes for a Claude Code session. Project: L-functions of SL₂(ℤ) modular forms (trivial multiplier) with poles at the cusp and/or in the interior, defined via finite S/T-segment contours. Last updated 2026-09-25. The companion write-up is `claude/principal_part_ambiguity.tex` in the project.

Tags:
- **[verified]**: reproduced by `verify_kernel.py` (Appendix). The current run gives 47/47 PASS with `--slow`.
- **[argued]**: proof sketch here, no machine check.
- **[reasoned, unchecked]**: a plausible argument that has not been tested.
- **[open]**: a research task.

---

## 0. Framing (read this first)

For f ∈ M^!_k with a pole at the cusp, the Bringmann–Fricke–Kent (BFK) L-series contains Γ(s, z) at the negative real argument z = −2πm, once for each polar Fourier mode. That point lies on the branch cut. The series therefore needs a choice of side, which BFK leave implicit. Diamantis–Lee–Raji–Rolen (DLRR) fix the principal branch.

**We do not claim a new term or invariant.** The "principal-part correction" is exactly the summed per-mode jump across that cut (§2.4). Its real-valued average is the principal value (§2.6).

What the work adds:
1. An explicit account of this choice.
2. A dictionary: side of the cut of Γ(s,·) = one-sided inversion of 1−T = direction of a horizontal Mellin ray = choice of T-coboundary on the cocycle side. This is the only item with real structural content.
3. Its consequences (§2.5–2.8).

Discipline on analytic continuation: f lives on ℍ, and ℝ is a natural boundary (dense cusps, exponential blow-up along geodesics into rationals when f has a pole at i∞). Contours of **f** are never deformed off ℍ. Only individual Fourier modes, which are entire, and the elementary kernels τ^w, h± are continued. The only "cusp residue" used is at q = 0, which is an isolated point in the q-disc.

## 1. Conventions (match `finite_contour_cocycles_short.tex`)

- Γ = SL₂(ℤ), S = (0 −1; 1 0), T = (1 1; 0 1). Weight k even, f = Σ_{n≥−M} c(n) qⁿ ∈ M^!_k (interior poles: §2.9).
- Segments: γ^S: −1/τ₀ → τ₀, γ^T: τ₀−1 → τ₀.
- **Kernel** (sign as in the .tex; the 2026-09-23 notebook page has "+", which is a typo):
  k̃_T(τ,s) = ζ(1−s, τ+1) − e^{iπ(s−1)} ζ(1−(k−s), τ+1).
- L*(f,s) = e^{−πis/2}[∫_{γ^S} f τ^{s−1} dτ + ∫_{γ^T} f k̃_T dτ]. At τ₀ = i the S-segment degenerates; all numerics use τ₀ = i.
- Principal branches throughout (arg ∈ (−π, π]). For τ ∈ ℍ, φ_s|_{2−k}S = e^{iπ(s−1)} τ^{k−1−s} exactly, where φ_s = τ^{s−1}.

## 2. Established facts

### 2.1 Operator form [argued]
k̃_T = Σ_{n≥1} (φ_s|_{2−k}(1−S))|Tⁿ. This gives k̃_T(τ) − k̃_T(τ−1) = −(φ_s|(1−S))(τ), i.e. a one-sided (1−T)⁻¹ of the S-coboundary of the Mellin kernel. At integer s this becomes the Bernoulli-polynomial coboundary inversion, which is Brown's recipe of killing C_T by a coboundary.

### 2.2 Non-uniqueness, and the two one-sided inversions [verified]
Basepoint independence uses only the shift equation. The solutions are therefore K₊ + {1-periodic functions}, so the draft's "unique" is too strong. The two natural solutions are:
- h₊[w](τ) = Σ_{n≥1}(τ+n)^w = ζ(−w, τ+1). This is the draft's kernel.
- h₋[w](τ) = −Σ_{n≤0}(τ+n)^w = −e^{iπw} ζ(−w, −τ). Only the Hurwitz function is evaluated at −τ; f never is.

Their difference is given by the Lipschitz formula:
h₊[s−1] − h₋[s−1] = (−2πi)^{1−s} Γ(1−s)^{−1} Li_s(q).

**Unfolding (Lemma 3 of the .tex) [verified, slow check].** For Re w < −1,
- ∫_{τ₀−1}^{τ₀} f h₊[w] dτ = ∫₀^∞ f(τ₀+t)(τ₀+t)^w dt, the horizontal ray to the right, inside ℍ;
- ∫_{τ₀−1}^{τ₀} f h₋[w] dτ = −∫_{−∞}^0 f(τ₀+t)(τ₀+t)^w dt, the ray to the left.

Via S, the second half of the kernel becomes a horocycle arc entering the cusp 0 tangentially.

**Admissible directions [argued].** If c(−M) ≠ 0 with M ≥ 1, a straight ray from τ₀ at angle θ ∈ (0, π) diverges, because of growth at i∞. Rays with θ ∈ (π, 2π) hit ℝ. So θ ∈ {0, π}, the two horizontal rays, are the only choices. They cannot be deformed into each other inside ℍ: the cusp blocks the route above, and the natural boundary blocks the route below.

### 2.3 Theorem: the difference is a residue at q = 0 [verified]
L*₊ − L*₋ = −2πi [ D(f,s)/((2π)^s Γ(1−s)) + i^k D(f,k−s)/((2π)^{k−s} Γ(1−k+s)) ], where D(f,s) = Σ_{m≥1} c(−m) m^{−s}.

Proof: right ray minus left ray is the full horizontal line, which equals ∫_seg f·Λ_w. That equals the residue at q = 0 of F·Λ_w/q, i.e. the coefficient extraction against the principal part.

Verified on Δ′ (principal part q⁻¹) and on g = jΔ′ − 744Δ′ (principal part q⁻², so D = 2^{−s}), at s = 3.3, 5+1.5i and 0.37, with errors ≤ 1e−26.

### 2.4 Proposition 9: global difference = summed per-mode branch jumps [verified]
Evaluate each polar-mode BFK summand at z → −2πm ± i0 using mpmath's own continuation, independently of the closed-form jump. The sum over modes of (upper − lower) reproduces L*₊ − L*₋ to ≤ 1e−25. Taking the wrong side is off by O(10⁻³–1).

So **h₊ is the principal (upper) side**, and it matches the BFK/DLRR convention to ~1e−28; h₋ is the lower side.

The Dirichlet-polynomial formula and the per-mode branch picture are the same ambiguity seen at two scales. The link is the Lipschitz formula, which is a computed identity and not visible from the shape of the formula.

### 2.5 Periodic freedom is larger than the branch choice [argued]
Adding any periodic p = Σ_{n≥0} p_n(s) qⁿ to the kernel shifts L* by e^{−πis/2} Σ_m c(−m) p_m(s). For example, p = q shifts L* by e^{−πis/2}c(−1), which is not a branch effect.

The functional equation only requires e^{−πis/2} p_m(s) to be symmetric under s ↔ k−s, so it leaves infinitely many choices. The branch choice is the one-parameter family λK₊ + (1−λ)K₋ built from one-sided T-orbit sums. Singling out one-sided sums is extra input.

### 2.6 Reality = principal value, and what it costs [verified]
At real s the two sides are exact complex conjugates, so L*₋ = conj(L*₊). Their mean is real and is the principal value. For s = 0 this is literal: the mean of the two sides of Γ(0, −2π) equals −Ei(2π) to 30+ digits, and Ei is a Cauchy PV integral. This is generic branch-cut behaviour, not a discovery.

The cost worth recording: the standard principal-branch values are **not real** at real s when D ≢ 0. Examples: Im L*₊(Δ′, 3.3) = 0.00308…, Im L*₊(g, 3.3) = 0.000507…. Intuition imported from M_k would say otherwise.

### 2.7 Where the branch choice does and doesn't matter [verified]
- It is invisible on M_k: D = 0, and the kernel reproduces (2π)^{−s}Γ(s)ζ(s)ζ(s−k+1) for G₄ and G₁₂ to 1e−33.
- It vanishes at every critical integer s = 1…k−1, closed form and direct: at s = 3, 6, 9, |L*₊ − L*₋| ≤ 1e−30. It is nonzero at s = 0 and s = k; for Δ′ at s → 0 the correction is −2πi.
- The functional equation holds for K₊, for K₋ and for affine combinations. It **fails** for the mixed kernel (h₊ in one half, h₋ in the other), by 0.014 and 0.38 for Δ′ and by 0.001 and 0.006 for g. So the functional equation forces the same side at both cusps.
- Hecke does not pick a side [argued]: D(f|T_p, s) = Σ_m c(−pm) m^{−s} + p^{k−1−s} D(f,s) has the same shape as the positive-index series.

### 2.8 Branch-explicit formula at τ₀ = i [verified]
Write E(s,m) = Σ_{j≥0} (2πm)^j/(j!(s+j)). This is single-valued and equals ∫₀¹ t^{s−1}e^{2πmt}dt for Re s > 0; use the series for Re s ≤ 0. Then

L*_sym(f,s) = Σ_{n>0} c(n)[Γ(s,2πn)/(2πn)^s + i^k Γ(k−s,2πn)/(2πn)^{k−s}] − c(0)(1/s + i^k/(k−s))
  + Σ_{m≥1} c(−m)[cos(πs)(Γ(s)(2πm)^{−s} + i^k Γ(k−s)(2πm)^{−(k−s)}) − E(s,m) − i^k E(k−s,m)],

and L*_± = L*_sym ∓ iπ[D(f,s)/((2π)^sΓ(1−s)) + i^k D(f,k−s)/((2π)^{k−s}Γ(1−k+s))].

### 2.9 Interior poles [argued, unchecked]
In the residue proof, the circle |q| = e^{−2πy₀} also encloses the images of poles lying above height y₀. The difference therefore acquires their residues and becomes height-dependent, consistent with the homotopy classes in the draft.

### 2.10 Data fixes
- Brown's Δ′ = E₆²E₄³/Δ + **(5541/16)**E₄³ − (1317/16)E₆² = q⁻¹ + 47709536q² + …. The 2019 email text reads 5514, which is a typo.
- The notebook's kernel sign is also a typo (§1).

## 3. Relation to the period/quasi-period theory (the long draft)

As described by the user in voice discussion. The long draft was not re-read in this session, so details are unverified here.
- Two ambiguities for periods of meromorphic forms:
  - (i) **Wall-crossing**: a contour crossing an interior pole shifts every special value by a residue.
  - (ii) **Subtraction**: the regularising subtraction is defined only modulo S^!_k, which has dimension 2·dim S_k.
- For dim S_k = 1 (k = 12, 16, 18, 20, 22, 26), the periods of Δ_k = ΔE_{k−12} and the quasi-periods of Δ̂_k = q⁻¹ + O(q) span the W⁺⊕W⁻ plane, because det pd_k ≠ 0.
- **No lattice.** The subtraction coefficients are arbitrary complex numbers, so the period vector modulo subtraction is a smear, i.e. a single point. The meaningful question is whether a **canonical section** exists, meaning a canonical rule for the subtraction. The h₊/h₋ choice is a small instance of the same question.
- **[reasoned, unchecked] The branch choice is transverse to the subtraction smear.** It vanishes at s = 1…k−1 (§2.7), so it never touches the V_{k−2} period polynomial or the W± coordinates. This assumes the long draft builds periods only from s = 1…k−1.
- **[reasoned, unchecked] Wall-crossing jumps are subtraction-invariant.** Subtracting a weakly holomorphic form adds no interior residues, so differences of periods across walls are well-defined even though absolute periods are smeared. Relative periods may therefore be where canonical content survives.

## 4. Open problems (ranked)

**Q1. Canonical section [open].** Is there a principled rule fixing the subtraction in the long draft, and the side in λK₊ + (1−λ)K₋? Candidates: Hecke (formally does not distinguish sides), a converse theorem in the style of DLRR, the cocycle side, reality (which picks λ = ½, the principal value).

**Q2. Interior poles, numerically [open].** Take a form with a pole in ℍ, preferably at a CM point (i or ρ) where the Laurent data is algebraic. Compute L*₊ − L*₋ at heights y₀ on either side of the pole and confirm the jump is the residue of f·Λ_w at the pole's q-image (§2.9). In the same computation, test the two unchecked claims of §3.

**Q3. Novelty check [open].** Do BFK, DLRR or later papers state (a) that the negative-index terms need a branch, and (b) that principal-branch values are non-real at real s? Also relate this to Bringmann's paper on L-functions for forms with poles, and to the missing citation of arXiv:1806.09874 flagged in the dim-S_k-one note.

**Q4. Hecke via the distribution relation [open].** Prove a segment-level Hecke formula using Σ_{j<N} ζ(w,(a+j)/N) = N^w ζ(w,a), and check it on Δ′ with T₂.

**Q5. Cohomological home of τ^{s−1} [open].** See Bruggeman–Choie–Diamantis, Mem. AMS 253 (2018) no. 1212 (arXiv:1404.6718). Express k̃_T as a (1−T)-coboundary inversion in their modules, and identify which modules the h₊ and h₋ choices live in (their cuts differ).

**Q6. Zagier/Kronecker identity on S^!_k with quasi-periods [open, low priority].** The Kronecker lead itself deflated: the Bernoulli factor cancels on pairing, leaving only consistency with Zagier '91 (the kernel reproduces his extended r_f on M_k). The existing literature (Choie–Park–Zagier; Choie; Blakestad–Choie) covers holomorphic forms only.

## 5. Skepticism notes
- All ingredients are classical: Lipschitz, Hankel, the branch jump of Γ(s,z), conjugate lateral values, principal values, natural boundaries. The contribution is explicitness plus the (1−T)⁻¹ ↔ branch-side dictionary.
- §3's two reasoned claims have not been tested. §2.9 has not been computed.

## References
- K. Bringmann, K.-H. Fricke, Z. Kent, Special L-values and periods of weakly holomorphic modular forms, Proc. AMS (2014).
- N. Diamantis, M. Lee, W. Raji, L. Rolen, L-series of harmonic Maass forms and a summation formula for harmonic lifts, IMRN 2023(18); arXiv:2107.12366.
- R. Bruggeman, Y. Choie, N. Diamantis, Holomorphic automorphic forms and cohomology, Mem. AMS 253 (2018) no. 1212; arXiv:1404.6718.
- F. Brown, A class of non-holomorphic modular forms III, arXiv:1710.07912.
- D. Zagier, Periods of modular forms and Jacobi theta functions, Invent. Math. 104 (1991).
- J. Lewis, D. Zagier, Period functions for Maass wave forms I, Ann. Math. 153 (2001).
- Y. Choie, Y. K. Park, D. Zagier, JEMS (2019); C. Blakestad, Y. Choie, arXiv:2404.06016.

## Appendix: `verify_kernel.py`
Requires `pip install mpmath`. Runtime is ~5 min, plus ~1–2 min with `--slow`. Expected output: `47/47 checks passed` with `--slow` (45/45 without). Each check prints PASS/FAIL against a stated tolerance.

```python
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
```
