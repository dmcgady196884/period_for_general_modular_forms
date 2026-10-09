# Is there a modular origin for the Hurwitz T-kernel?

Kickoff notes for a fresh Claude Code session. Project: finite-contour L-functions and period
polynomials for SL₂(ℤ) forms, `finite_contour_cocycles_short.tex`. Date: 2026-09-23.

Tags: **[verified]** = reproduced numerically this session; **[proved in draft]** = has a label
in the .tex; **[argued]** = short argument here, no machine check; **[lead]** = plausible,
unverified, possibly already in the literature; **[open]** = the work.

---

## 0. The question

DAM's framing, lightly edited:

> The Riemann zeta function is the L-function of the Jacobi theta function. The Hurwitz zeta
> function is a deformation of it, and it shows up in exactly the kernel needed to make
> finite-contour integrals of generic modular forms land in W. Is there some hidden modular
> origin story for this kernel — a duality from the space of modular forms to itself?

A previous session (see `hurwitz_kernel_handoff_20260923.md`) asked a neighbouring question —
*is the kernel unique?* — and found a one-parameter ambiguity. That is useful input, but it is not
this question. This session should answer the original one.

## 1. Working rules (DAM's standing preferences)

- Reason analytically first. Numerics only to test predictions that were derived and stated
  before the run.
- No edits to the `.tex` without showing the diff. Never compile it; DAM compiles.
- Never `git commit`. Edit tracked files with Edit/Write, not `sed -i` or heredoc Python.
- Python: `/Users/dmcgady/Documents/math/mmf_venv/bin/python` (mpmath). Run numerics in the
  foreground with flushed prints and bounded loops. Say what a slow run computes before firing it.
- Literature-standard terms only. No invented jargon.
- Before believing a PSLQ hit, check coefficient height against digits supplied.

## 2. What is established

**Kernel.** [proved in draft, `def:kernel`, line ~190]

  k̃_T(τ,s) = ζ(1−s, τ+1) − e^{iπ(s−1)} ζ(1−(k−s), τ+1),

and at integer s = ℓ+1, 0 ≤ ℓ ≤ n = k−2, it assembles into the V_n-valued kernel
**K**_T(τ;X,Y) = Σ_ℓ C(n,ℓ)(−1)^ℓ k̃_T(τ,ℓ+1) X^{n−ℓ}Y^ℓ of `prop:periodint`.

**Shift identity.** [proved, `lem:resorbit`; verified to 1e−30 today]
k̃_T(τ) − k̃_T(τ−1) = −(φ_s|_{2−k}(1−S))(τ) with φ_s(τ) = τ^{s−1}. Equivalently
**K**_T(τ) = Σ_{m≥1} P_{τ+m}|_{(1−S)} (Hurwitz-regularised), with P_τ = (X−τY)^n.
So the kernel is the **one-sided** (1−T)^{−1} applied to the S-coboundary of the Mellin kernel.

**The S-relation survives at complex s.** [proved, `eq:ktilflip`, line 1474]
k̃_T(τ,k−s) = −e^{iπ(k−s−1)} k̃_T(τ,s). This is the complex-s form of **K**_T|_{(1+S)} = 0. Note
it is *pointwise in τ* and relates s to k−s.

**Lipschitz ambiguity.** [verified; see the other handoff §1.2]
The two one-sided inverses h₊ (draft) and h₋ differ by the two-sided sum

  Σ_{n∈ℤ}(τ+n)^{s−1} = (−2πi)^{1−s} Γ(1−s)^{−1} Li_s(q),

which is the c = 1 row of an Eisenstein series at complex weight. It vanishes identically at
integer 1 ≤ s ≤ k−1, so periods and residue polynomials do not see the choice.

**Growth.** [verified today] |k̃_T(iy+x,s)| ~ C·y^{max(Re s, k−Re s)} with C = 1/|s| or 1/|k−s|
as the Hurwitz asymptotic predicts (s = 3.3: 0.11494 ≈ 1/8.7; s = 6: 0.33333 = 1/6+1/6). Polynomial,
never exponential. A general 1-periodic correction with negative q-powers *does* blow up, so
"polynomial growth at the cusp" is a natural first constraint on the ambiguity.

**Holomorphic case.** [verified] On M_k the kernel gives the classical continuation
(2π)^{−s}Γ(s)ζ(s)ζ(s−k+1) for G_k to 1e−28, including the 1/s, 1/(k−s) poles.

## 3. What was tried today, and why it failed

**Conjecture tested.** The scalar kernel g(τ) = k̃_T(τ,s), with the weight 2−k action
(g|γ)(τ) = (cτ+d)^{k−2} g(γτ), is itself a cocycle value at complex s.

**Result.** [verified] False, at O(1), even at s = 6 where g is a Bernoulli polynomial:

- s = 3.3: |g + g|S|/|g| ≈ 0.10 and 0.03; |g|(1+U+U²)|/|g| ≈ 5 and 40.
- s = 5+1.5i: ≈ 1.2 and 1.0; ≈ 15 and 89.
- s = 6: ≈ 0.80 and 0.80; ≈ 5 and 41.
- Shift identity at the same points: ≤ 1e−30.

**Why — the trap to avoid.** The integer-s relation **K**_T|_{(1+S)} = 0 is a slash in (X,Y). That
permutes *coefficients* (ℓ ↔ n−ℓ) at fixed τ; it is not a transformation of τ. On P_τ the two
actions are linked — P_τ|_{(X,Y)}γ = (−cτ+a)^n P_{γ^{−1}τ} — so under the τ-action S sends the
s-coefficient function to the (k−s) one, which is φ_s|_{2−k}S = e^{iπ(s−1)}φ_{k−s}. But **K**_T is a
cocycle piece, not an equivariant object, so no single-coefficient τ-relation should be expected.
The correct complex-s shadow of the S-relation is `eq:ktilflip`, pointwise in τ.

## 4. The correct setting for the U-relation [argued]

Identify V_n with polynomials in x of degree ≤ n under (p|γ)(x) = (cx+d)^n p((ax+b)/(cx+d)).
The ℓ-th coefficient sits on the monomial x^{n−ℓ}. At complex s the natural extension is to
monomials x^α of complex degree α = n−ℓ = k−1−s — the weight-(−n) principal series.

On these, (x^α|γ)(x) = (ax+b)^α (cx+d)^{n−α}:

- S: x^α ↦ e^{±iπα} x^{n−α}. This is s ↔ k−s with a phase; branch choice matters.
- T: x^α ↦ (x+1)^α = Σ_{j≥0} C(α,j) x^{α−j} for |x| > 1 (or Σ C(α,j) x^j for |x| < 1).
- U = TS: x^α ↦ (x−1)^α x^{n−α} = Σ_{j≥0} C(α,j)(−1)^j x^{n−j} for |x| > 1.

At integer α these series are finite and reproduce the binomial matrices of `common.opmat`. At
complex α they are **infinite**. U sends a complex-degree monomial into *integer* degrees n−j,
including negative ones, and U² sends it to the coset n−α−ℤ≥0. So the complex-s U-relation is an
infinite-series identity linking k̃_T at s, at k−s+ℤ, and at integers.

Two consequences to keep in view:

- The choice of expansion region (|x| > 1 or |x| < 1) is a choice of module. It is a natural
  suspect for the h₊/h₋ ambiguity, and for Bruggeman–Choie–Diamantis's distinction between
  modules of functions extending across different parts of ℙ¹(ℝ).
- **Complex s cannot be moved into the exponent of the period kernel.** For f of weight k,
  f(τ)(x−τ)^w dτ is Γ-invariant under the joint action only for w = n. Under the full PSL₂(ℝ),
  which acts transitively on pairs (x, τ) ∈ ℙ¹(ℝ) × ℍ, a jointly equivariant holomorphic kernel is
  forced to be c·(x−τ)^n. Under Γ alone it may carry a Γ-invariant factor, which only changes f.
  So complex s lives in the Mellin/coefficient decomposition, never in a deformed period kernel.

## 5. The strongest lead: Kohnen–Zagier kernels and the Cohen kernel [lead]

This is the most promising answer to §0, and it may already be in the literature.

- **Integer s.** Kohnen–Zagier (1984, already cited in the draft as `KZ1984` for `def:Wpm`)
  define cusp forms
  R_n(τ) ∝ Σ_{(a b; c d)∈Γ} (aτ+b)^{−n−1}(cτ+d)^{−(k−1−n)}, 0 < n < k−2,
  with ⟨f, R_n⟩ = r_n(f), the n-th period of f. Their "rational periods" theorem computes the
  periods of R_n in terms of **Bernoulli numbers**. The draft's kernel at integer s *is* a
  Bernoulli-polynomial object, and `def:Wpm` already uses KZ's Bernoulli bases. **Check whether
  **K**_T at integer s reproduces KZ's period formula for R_n.** If so, the integer-s origin
  has been known since 1984.
- **Complex s.** The same sum with exponents −s and −(k−s),
  C_k(τ; s) = Σ_{γ∈Γ} (aτ+b)^{−s}(cτ+d)^{−(k−s)} = Σ_γ (τ^{−s})|_k γ,
  is a weight-k Poincaré series of the Mellin seed, attributed to H. Cohen. Diamantis–O'Sullivan,
  "Kernels of L-functions of cusp forms", Math. Ann. 346 (2010) [**verify citation**], show
  ⟨f, C_k(·;s)⟩ ∝ L*(f,s) for cusp forms, in 1 < Re s < k−1.
- **The conjectural link.** The (k−1)-fold antiderivative of τ^{−s} is ∝ τ^{k−1−s} = φ_{k−s},
  so the Eichler integral of C_k(·;s) is, modulo polynomials, the full Γ-orbit sum
  Σ_γ φ|_{2−k}γ. The Hurwitz kernel is its one-sided Γ_∞-partial sum of the S-coboundary.
  And Stokes applied to ⟨f, C_k⟩ over the standard fundamental domain gives boundary integrals
  over its walls: the bottom arc (the S-wall) and the vertical sides (the T-walls). That is
  exactly the draft's geodesic reference contour, γ^{arc} and its T-translates.

So the concrete conjecture is:

> The finite-contour L*(f,s), with S-kernel τ^{s−1} and the Hurwitz T-kernel, is the boundary
> form of the Petersson pairing ⟨f, C_k(·;s)⟩. For f ∈ M^!_k or F_k it is a *regularised*
> Petersson pairing, and the one-sided choice h₊ corresponds to a choice of regularisation.

If true, the "duality from modular forms to itself" DAM asked about is f ↦ ⟨f, C_k(·;s)⟩, the
Hurwitz sum is the part of the Cohen kernel's Γ-orbit that the finite contour sees, and the
Lipschitz ambiguity is the missing half of the c = 1 Eisenstein row.

## 6. Tasks, in order

**T1 — literature gate (do first; may end the project).**
Read KZ1984 §1 and Diamantis–O'Sullivan. Extract: the definition and convergence range of R_n and
C_k(·;s); the exact constant in ⟨f, C_k(·;s)⟩ ∝ L*(f,s); any computation of the period polynomial
or Eichler integral of C_k(·;s), especially one written with Hurwitz zeta. Also search for
extensions to weakly holomorphic or meromorphic f (Bringmann–Kane–von Pippich on regularised
inner products are the natural place). **If D–O'S or a successor already express the Eichler
integral of the Cohen kernel through Hurwitz zeta, the modular origin is known: report that, with
section/equation references, and stop.**

**T2 — integer-s check against KZ.** Compare **K**_T at s = ℓ+1 (k = 12, and one weight with
dim S_k = 2, e.g. k = 24) with KZ's formula for the periods of R_n. Prediction: they agree up to
normalisation and possibly a coboundary. This is cheap: exact Bernoulli arithmetic.

**T3 — normalisation at complex s.** For f = Δ (cusp form, no Lipschitz ambiguity), compare the
draft's L*(Δ, s) with ⟨Δ, C_12(·;s)⟩ at s = 3.3 and 5+1.5i. Prediction: equal up to the D–O'S
constant. C_12(·;s) converges slowly; either sum it with care or use the D–O'S constant from T1
and treat this as a normalisation check only.

**T4 — the new content: derive the finite contour from ⟨f, C_k⟩.** On paper first. Write
⟨f, C_k(·;s)⟩ as a boundary integral over ∂𝓕 via the Eichler integral of C_k and Stokes, and
identify the T-wall contribution with ∫_{γ^T} f k̃_T. The prediction to test: the one-sided
Hurwitz sum comes from organising the Γ_∞-cosets along the T-walls, and the one-sidedness is a
choice of fundamental domain for Γ_∞ there. Only then check numerically.

**T5 — weakly holomorphic and meromorphic f.** For Brown's Δ′ (other handoff §1.6), ask whether
L*₊(Δ′,s), L*₋(Δ′,s) or their average equals a regularised ⟨Δ′, C_12(·;s)⟩ (Borcherds, or BKvP).
This is where the other handoff's P1 (canonical normalisation) and the review's §8 (BKvP) meet.
Caution: memory records that eq:inner ≠ Petersson at k = 24 and that Haberland fails for
meromorphic f by ~24%. "L* = regularised Petersson" for meromorphic f is genuinely uncertain.

**T6 — optional: the U-relation on the s-line (§4).** Formalise the module of §4, apply the
infinite binomial action to {k̃_T(τ,s)}, and test (1+U+U²) up to an Eisenstein counterterm. Do
this only if T4 fails to give a clean statement; T4 would imply it.

## 7. Pitfalls

- **Wrong action** (§3): (X,Y)-slash ≠ τ-action on a single coefficient.
- **Branches**: (aτ+b)^{−s}, x^α and e^{iπ(s−1)} all need consistent principal-branch
  conventions. The h₊/h₋ split *is* a branch choice (other handoff §1.3).
- **Sign bookkeeping** between (−1)^ℓ at integer ℓ and e^{iπ(s−1)} at complex s: redo it
  explicitly; don't trust memory of it.
- **Normalising a near-zero vector by its own magnitude** manufactures spurious rank. Use a global
  scale.
- **Mesh clustering** (common.py TRAP 4): cluster where the poles are, not at parameter endpoints.
- **PSLQ**: a 3-term relation of height 5·10⁶ needs ~25 digits. A 16-digit "hit" at that height was
  spurious on 2026-09-15.

## 8. Pointers

- Draft: `finite_contour_cocycles_short.tex` — `def:kernel` (~190), `prop:periodint` (508),
  `lem:cocycle` (577), `lem:rfW` (628), `def:Wpm` (654), `lem:tau0indep` (807), `lem:arcfunceq`
  (1455), `eq:ktilflip` (1474), `lem:resorbit` (1715). Line numbers as of 2026-09-23.
- `hurwitz_kernel_handoff_20260923.md` — Lipschitz ambiguity, BFK/DLRR branch, reality, P1–P5.
  Its pointers into the draft are to an older version.
- `closed_contour_review_20260922.md` — §1's elementary cocycle (c_S, c_U) = (R|(1−S), 0); §8 BKvP.
- `arc_numerics/common.py` — `ktil`, `slash`, `opmat`, `rvec`, `arcint`, `E4`, `E6`, `Delta`.
- `ai_collab/dr_b_llm_episodes.md` — Episode 14 and its addendum.
- Memory: `haberland_petersson_match`, `haberland_fails_for_meromorphic`,
  `lit_scan_2026-08-17_brown_fonseca`, `residue_polynomials`.

## Appendix: reproduce §2's growth check and §3's failed test

Run from `arc_numerics/`.

```python
from common import mp, I, ktil
mp.mp.dps = 30; k = 12
sl = lambda F, a, b, c, d: (lambda t: (c*t + d)**(k - 2) * F((a*t + b)/(c*t + d)))
for s in (mp.mpf('3.3'), mp.mpc(5, 1.5), mp.mpf(6)):
    G = lambda t, s=s: ktil(t, s, k)
    GS, GU = sl(G, 0, -1, 1, 0), sl(G, 1, -1, 1, 0)
    GU2 = sl(GU, 1, -1, 1, 0)
    for t in (mp.mpc('0.3', '1.2'), mp.mpc('-0.2', '0.9')):
        g = G(t)
        shift = G(t) - G(t - 1) + t**(s - 1) - t**(k - 2) * (-1/t)**(s - 1)
        print(mp.nstr(s, 3), mp.nstr(t, 3),
              'S:', mp.nstr(abs(g + GS(t))/abs(g), 4),
              'U:', mp.nstr(abs(g + GU(t) + GU2(t))/abs(g), 4),
              'shift:', mp.nstr(abs(shift)/abs(g), 4))
    e = max(mp.re(s), k - mp.re(s))
    for y in (10, 100, 1000):
        print('   growth y=%d:' % y, mp.nstr(abs(ktil(mp.mpf('0.3') + I*y, s, k))/mp.mpf(y)**e, 8))
```

Expected: S and U columns O(1); shift column ≤ 1e−29; growth ratios converging to 0.1149 (s=3.3),
0.3333 (s=6), 0.0132 (s=5+1.5i).
