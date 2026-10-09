# Hurwitz T-kernel: verified structure, a branch ambiguity, and open problems

Handoff notes for a Claude Code session. Project: L-functions of SL₂(ℤ) modular forms (trivial multiplier) with poles at the cusp and/or in the interior, defined via finite S/T-segment contours. Date: 2026-09-23.

Everything marked **[verified]** is reproduced by `verify_kernel.py` (Appendix). Everything marked **[argued]** has a short proof sketch here but no machine check. **[open]** items are the research tasks.

---

## 0. Conventions (match `finite_contour_cocycles_short.tex`)

- Γ = SL₂(ℤ), S = (0 −1; 1 0), T = (1 1; 0 1). Weight k even, f ∈ M^!_k (or F_k with interior poles).
- Segments: γ^S_{τ₀}: −1/τ₀ → τ₀, γ^T_{τ₀}: τ₀−1 → τ₀.
- **Kernel (sign as in the .tex; the notebook page of 2026-09-23 has "+", which is a typo):**

  k̃_T(τ,s) = ζ(1−s, τ+1) − e^{iπ(s−1)} ζ(1−(k−s), τ+1).

- L*(f,s) = e^{−πis/2} [ ∫_{γ^S} f(τ) τ^{s−1} dτ + ∫_{γ^T} f(τ) k̃_T(τ,s) dτ ],  L(f,s) = (2π)^s Γ(s)^{-1} L*(f,s).
- At τ₀ = i the S-segment drops out, so L* = e^{−πis/2} ∫_{i−1}^{i} f k̃_T dτ. All numerics use this.
- Principal branches throughout (arg ∈ (−π, π]).

## 1. Structural facts

### 1.1 Operator form [argued; consistent with Lemma 4.1 of the project draft]
Let φ_s(τ) = τ^{s−1} and use the weight 2−k action (g|_{2−k}γ)(τ) = (cτ+d)^{k−2} g(γτ). For τ ∈ ℍ, −1/τ ∈ ℍ, and arg(−1/τ) = π − arg τ. Therefore

  φ_s |_{2−k} S = e^{iπ(s−1)} τ^{k−1−s}  (exactly, with no branch fudge),

and

  **k̃_T = Σ_{n≥1} (φ_s |_{2−k}(1−S)) | T^n**,  so  k̃_T(τ) − k̃_T(τ−1) = −(φ_s|(1−S))(τ).

In words, the kernel is the one-sided formal inverse (1−T)^{−1} applied to the S-coboundary of the Mellin kernel. At integer s it becomes the Bernoulli-polynomial coboundary inversion. This is exactly Brown's recipe of killing C_T with a coboundary P, as in his 2019 email and in arXiv:1710.07912.

The first half, Σ_{n≥1}(τ+n)^{s−1} = ζ(1−s, τ+1), is Mayer's transfer operator for the Gauss map applied to the constant function 1. Lewis–Zagier period functions are built from the same Σ_{n≥1} T^n.

### 1.2 The kernel is NOT unique: the Lipschitz ambiguity [verified]
Any 1-periodic p(τ) can be added to the kernel without breaking the shift identity or basepoint independence, because d/dτ₀ ∫_{τ₀−1}^{τ₀} f p = 0. The draft's phrase "unique V_n-analog" is too strong.

The natural competitor is the other one-sided inverse:
- h₊[w](τ) = Σ_{n≥1}(τ+n)^w = ζ(−w, τ+1). This is the draft's choice.
- h₋[w](τ) = −Σ_{n≤0}(τ+n)^w = −e^{iπw} ζ(−w, −τ).

Both satisfy h(τ) − h(τ−1) = −τ^w. Their difference is given by the Lipschitz formula (checked to 1e−26):

  **h₊[s−1] − h₋[s−1] = (−2πi)^{1−s} Γ(1−s)^{−1} Li_s(e^{2πiτ}).**

Let L*₊ and L*₋ be L* built with the h₊ and h₋ kernels. For f ∈ M^!_k with f = Σ c(n) qⁿ, the pairing over the T-segment extracts the constant term of f·Li_s(q), which gives

  L*₊(f,s) − L*₋(f,s) = e^{−πis/2} [ C(s) D₋(f,s) − e^{iπ(s−1)} C(k−s) D₋(f,k−s) ],
  where C(a) = (−2πi)^{1−a}/Γ(1−a) and D₋(f,s) = Σ_{m≥1} c(−m) m^{−s}.

This is a finite Dirichlet polynomial in the principal part. Checked on Brown's Δ' at s = 3.3, 6, and 5+1.5i to about 1e−26. With interior poles there are extra residue terms from f·Li_s(q) [argued, not checked].

Consequences:
- **Cusp forms and M_k don't see the ambiguity** (D₋ = 0). [verified: on G₄ and G₁₂ the kernel reproduces (2π)^{−s}Γ(s)ζ(s)ζ(s−k+1) to about 1e−28.]
- **The ambiguity vanishes at every integer s with 1 ≤ s ≤ k−1**, because 1/Γ(1−s) = 0 and 1/Γ(1−(k−s)) = 0. So critical values and periods are canonical. [verified at s = 6 for k = 12]
- **Asymptotics cannot distinguish h₊ from h₋.** Their difference decays like q, so both have the same power-law expansion as Im τ → ∞. What differs is the direction of analytic continuation across ℝ: h₊ is cut along (−∞, −1] and h₋ along [0, ∞).

### 1.3 The branch in Bringmann–Fricke–Kent / DLRR is this same ambiguity [verified]
BFK's per-mode formula contains Γ(s, 2πn t₀)/(2πn)^s for n < 0, which needs a branch. The jump across the negative axis of z^{−s}Γ(s,z) is −2πi a^{−s}/Γ(1−s) at z = −a. With a = 2π|n| this is exactly the Lipschitz term of §1.2. Numerically:

- **The Hurwitz kernel h₊ equals the BFK formula with principal branches (arg = +π on the negative axis)**, to about 1e−23 on Δ'. This is the convention Diamantis–Lee–Raji–Rolen state explicitly (arXiv:2107.12366, eq. (3.3), −π < arg ≤ π). Their L_f(I_s 1_{t₀}) recovers the unsymmetrised BFK series.
- **h₋ corresponds to the arg = −π branch.**

### 1.4 Reality [verified numerically; easy proof]
For f with real Fourier coefficients and real s, **L*₋ = conj(L*₊)**. For example, for Δ' at s = 3.3:

  L*₊ = −4.269149131848… + 0.003080441651…·i.

So the Hurwitz-kernel L-function (equivalently principal-branch BFK/DLRR) is **not real on the real axis** once f has a principal part. Its imaginary part is exactly half the Lipschitz term. The reflection-symmetric kernel ½(h₊ + h₋) = h₊ − ½·Lipschitz gives a real L*, equal to Re L*₊ for real s. Proof idea: τ ↦ −τ̄ swaps the two one-sided sums up to conjugation, and f(−τ̄) = conj f(τ).

### 1.5 Hecke does not fix the branch [argued]
Under T_p, c_{f|T_p}(−m) = c(−pm) + p^{k−1} c(−m/p). So D₋(f|T_p, s) = Σ_m c(−pm) m^{−s} + p^{k−1−s} D₋(f,s). This has the same formal shape as the positive-index Dirichlet series. Any kernel h₊ − λ(s)·Lipschitz therefore behaves the same way under Hecke, and Hecke equivariance alone cannot single out λ = 0 or λ = ½. The functional equation s ↔ k−s is also respected by the whole family: the two halves of the ambiguity swap.

### 1.6 Data fix
Brown's Δ' has constants **(5541/16) E₄³ − (1317/16) E₆²**. The 2019 email text reads 5514/16, which is a typo. Only 5541 gives Δ' = q^{−1} + 0 + 0·q + 47709536 q² + ….

## 2. The Kronecker-function lead: checked, mostly deflated

Hypothesis from the earlier chat: the kernel's integer-s Bernoulli structure is the q → 0 degeneration of Zagier's Kronecker function F_τ(u,v) = θ'(0)θ(u+v)/(θ(u)θ(v)) (Invent. Math. 104, 1991).

What actually holds:
1. Generating function over integer s. Since ζ(1−s,a) = −B_s(a)/s,
   Σ_{s≥1} ζ(1−s, τ+1) z^{s−1}/(s−1)! = 1/z − e^{τz}/(1−e^{−z}).
   Pairing over the T-segment, the factor 1/(1−e^{−z}) **cancels exactly** against ∫_{τ₀−1}^{τ₀} e^{(z+2πin)τ} dτ:
   ∫_{τ₀−1}^{τ₀} f(τ)[1/z − e^{τz}/(1−e^{−z})] dτ = c(0)/z − e^{zτ₀} Σ_n c(n) q₀ⁿ/(z+2πin).
   The result is a partial-fraction (Laplace/Borel) transform with simple poles at z = −2πin and residue −c(n), independent of τ₀. Principal-part coefficients give poles on the positive imaginary axis. No theta function appears.
2. The q⁰ term of Zagier's Main Theorem is P(X,T) = ¼ coth(XT/2) coth(T/2), a product of two Bernoulli generating functions. It comes from L(G_k,s) = ζ(s)ζ(s−k+1). One factor is "Hurwitz at a = 1" and is kernel-like; the other comes from the divisor sums of G_k, not from the kernel.
3. What survives: on M_k the Hurwitz-kernel L* equals the classical continuation [verified]. Its poles at s = 0, k have residues −c(0) and (−1)^{k/2}c(0), which match Zagier's eq. (9)–(10). So **the kernel reproduces Zagier's extended period map r: M_k → Ŵ_k**, including the X^{−1} and X^{k−1} terms. That is consistency, not a new structure.

Verdict: there is no evidence that the kernel "is" a Kronecker degeneration. The interesting question the lead points to is §3 P4.

## 3. Open problems (ranked)

**P1. Canonical normalisation on M^!_k.** Classify the 1-periodic corrections p(τ, s) that preserve (a) basepoint independence, (b) the functional equation L*(s) = i^k L*(k−s), (c) the classical values on S_k, and (d) Hecke compatibility in Guerzhoy's weak-Hecke sense (f|T_p ≡ λ_p f mod D^{k−1}M^!_{2−k}).
- Conjecture to test: (a)–(d) leave exactly the one-parameter family h₊ − λ(s)·Lipschitz, plus constant-term corrections c(s)·c_f(0). Reality then picks λ = ½.
- Concrete test: compute L* on D^{k−1}g for g ∈ M^!_{2−k}. For example, k = 12 with g = E₄E₆²/Δ² or similar, or the Bol image of a weight −10 form. See whether any λ makes L* vanish, or behave well, on the Bol image at non-critical s. Bol's identity implies L*(D^{k−1}g, s) is a shift of L*(g, s−k+1), so first work out what "should" happen.

**P2. Hecke via the distribution relation.** Prove a segment-level Hecke formula for the kernel pairing using Σ_{j<N} ζ(w, (a+j)/N) = N^w ζ(w, a). Check it numerically on Δ' with T₂.

**P3. Cohomological home of τ^{s−1}.** The draft's Outlook conjectures a principal-series local system. First check Bruggeman–Choie–Diamantis, *Holomorphic Automorphic Forms and Cohomology*, Mem. AMS 253 (2018), no. 1212 (arXiv:1404.6718). It handles holomorphic functions with the |_{2−k} action for arbitrary real weight, and mixed parabolic cohomology. Task: express k̃_T as a (1−T)-coboundary inversion in their modules, and identify which module (boundary-germ / analytic-vector type) the h₊ vs h₋ choice lives in. The different cuts in §1.2 suggest it is exactly their distinction between modules of functions extending across different parts of ℙ¹(ℝ).

**P4. A Zagier/Kronecker identity on S^!_k.** S^!_k / D^{k−1}M^!_{2−k} has dimension 2·dim S_k (Bringmann–Guerzhoy–Kent–Ono) and carries periods and quasi-periods (η±; for Δ', η₊ = 127202100647.177…, η₋ = 10276732343.649…). Question: does Zagier's Main Theorem admit an extension whose coefficients involve quasi-periods, with the Hurwitz-kernel L-values as input? Known literature (Choie–Park–Zagier 2019 for Γ₀(N); Choie 2021 for Hilbert forms; Blakestad–Choie 2024 for twisted versions) covers **holomorphic forms only**. I found no weakly holomorphic version in a quick search, but this is not a thorough novelty check. Search Pasol–Popa, Kent, Guerzhoy, and Brown's "multiple modular values" before investing.

**P5. Interior poles.** Extend §1.2 to f ∈ F_k. The ambiguity should gain residues of f·Li_s(q) at poles inside the T-strip. Check against Proposition 5.7 of the draft (the polylog projection P̂).

## 4. Skepticism notes
- §1.2–1.4 are elementary once seen. The Lipschitz formula and the incomplete-gamma branch jump are classical. The contribution is identifying that the branch implicit in BFK's L-series equals the choice of one-sided (1−T)^{−1}, and that the resulting L* is non-real. I have not checked whether BFK, DLRR, or later papers remark on this non-reality. **Check before claiming.**
- P1's conjecture is a guess. P3's identification with BCD is plausible but unverified.

## References
- D. Zagier, Periods of modular forms and Jacobi theta functions, Invent. Math. 104 (1991) 449–465.
- K. Bringmann, K.-H. Fricke, Z. Kent, Special L-values and periods of weakly holomorphic modular forms, PAMS (2014).
- N. Diamantis, M. Lee, W. Raji, L. Rolen, L-series of harmonic Maass forms and a summation formula for harmonic lifts, IMRN 2023(18); arXiv:2107.12366.
- R. Bruggeman, Y. Choie, N. Diamantis, Holomorphic automorphic forms and cohomology, Mem. AMS 253 (2018) no. 1212; arXiv:1404.6718.
- J. Lewis, D. Zagier, Period functions for Maass wave forms I, Ann. Math. 153 (2001).
- F. Brown, A class of non-holomorphic modular forms III, arXiv:1710.07912.
- Y. Choie, Y. K. Park, D. Zagier, Periods of modular forms on Γ₀(N) and products of Jacobi theta functions, JEMS (2019).
- C. Blakestad, Y. Choie, Twisted Kronecker series and period polynomials on Γ₀(N), arXiv:2404.06016.

## Appendix: `verify_kernel.py`
Requires `pip install mpmath`. Runtime is about 1 minute. Expected output: all errors ≤ 1e−22.

```python
# Reproduces every numerical claim in hurwitz_kernel_handoff.md.  pip install mpmath
from mpmath import mp, mpf, mpc, zeta, exp, pi, gamma, rgamma, polylog, quad, gammainc, bernoulli
from fractions import Fraction as Fr
mp.dps = 25; I = mpc(0, 1)

# --- one-sided (1-T)^{-1} inverses of tau^w, both solve h(t)-h(t-1) = -t^w ---
def h_plus(w, t):  return zeta(-w, t + 1)                    #  sum_{n>=1} (t+n)^w
def h_minus(w, t): return -exp(I*pi*w) * zeta(-w, -t)        # -sum_{n<=0} (t+n)^w

def kernel(t, s, k, h):  # k~_T with the .tex sign:  h[s-1] - e^{i pi (s-1)} h[k-1-s]
    return h(s - 1, t) - exp(I*pi*(s - 1)) * h(k - 1 - s, t)

def Lstar(c, s, k, h):   # tau0 = i: S-segment vanishes, T-segment = [i-1, i]
    f = lambda t: sum(mpf(v.numerator)/v.denominator * exp(2*pi*I*n*t) for n, v in c.items())
    return exp(-pi*I*s/2) * quad(lambda x: f(x + I) * kernel(x + I, s, k, h), [-1, -0.5, 0])

# --- q-series ---
N = 40
def sig(n, e): return sum(d**e for d in range(1, n + 1) if n % d == 0)
def mul(a, b):
    out = [0]*N
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b[:N - i]): out[i + j] += x*y
    return out
E4 = [1] + [240*sig(n, 3) for n in range(1, N)]
E6 = [1] + [-504*sig(n, 5) for n in range(1, N)]
P = [1] + [0]*(N - 1)                      # Delta/q
for n in range(1, N):
    fac = [0]*N; fac[0] = 1; fac[n] = -1
    for _ in range(24): P = mul(P, fac)
inv = [0]*N; inv[0] = 1
for n in range(1, N): inv[n] = -sum(P[j]*inv[n - j] for j in range(1, n + 1))
E43, E62 = mul(E4, mul(E4, E4)), mul(E6, E6)
A = mul(mul(E62, E43), inv)                # q * E6^2 E4^3 / Delta
# Brown's Delta' = E6^2E4^3/Delta + (5541/16)E4^3 - (1317/16)E6^2 = q^-1 + 47709536 q^2 + ...
DP = {n: Fr(A[n + 1]) + (Fr(5541, 16)*E43[n] - Fr(1317, 16)*E62[n] if n >= 0 else 0) for n in range(-1, N - 2)}
assert DP[-1] == 1 and DP[0] == 0 and DP[1] == 0 and DP[2] == 47709536

if __name__ == "__main__":
    # (1) Lipschitz: h+ - h- = (-2 pi i)^{1-s}/Gamma(1-s) Li_s(q)
    t = mpc(0.3, 0.8)
    for s in [mpc(-1.3, 0.7), mpc(0.4, 2), mpf(3.7)]:
        d = h_plus(s - 1, t) - h_minus(s - 1, t)
        print("Lipschitz err", abs(d - (-2*pi*I)**(1 - s)*rgamma(1 - s)*polylog(s, exp(2*pi*I*t))))
    # (2) kernel ambiguity on Delta' and match with BFK/DLRR principal branch
    k = 12
    C = lambda a: (-2*pi*I)**(1 - a) * rgamma(1 - a)
    for s in [mpf('3.3'), mpf(6), mpc(5, 1.5)]:
        Lp, Lm = Lstar(DP, s, k, h_plus), Lstar(DP, s, k, h_minus)
        pred = exp(-pi*I*s/2) * (C(s) - exp(I*pi*(s - 1))*C(k - s))     # times sum_m c(-m) m^{-s} = 1
        bfk = sum(mpf(v.numerator)/v.denominator*(gammainc(s, 2*pi*n)/(2*pi*n)**s
                  + gammainc(k - s, 2*pi*n)/(2*pi*n)**(k - s)) for n, v in DP.items() if n != 0)
        print(f"s={s}\n  L+={Lp}\n  L-={Lm}\n  (L+ - L-) - pred = {abs(Lp - Lm - pred)}\n  L+ - BFK_principal = {abs(Lp - bfk)}")
    # (3) on M_k the kernel gives the classical continuation (2pi)^-s Gamma(s) zeta(s) zeta(s-k+1) for G_k
    for k in [4, 12]:
        c0 = -bernoulli(k)/(2*k)
        cG = {n: Fr(sig(n, k - 1)) for n in range(1, N)}
        f0 = lambda s: exp(-pi*I*s/2)*c0*quad(lambda x: kernel(x + I, s, k, h_plus), [-1, 0])
        for s in [mpf('2.5'), mpc(1.7, 0.9)]:
            L = Lstar(cG, s, k, h_plus) + f0(s)
            print(f"G_{k} s={s} err", abs(L - (2*pi)**(-s)*gamma(s)*zeta(s)*zeta(s - k + 1)))
```
