# Closed-contour object $O(z)$: objections, refinements, cross-checks, and open tasks

**Purpose.** This is a hand-off for a Claude Code session. It records a critical assessment of one claim from the meromorphic period-polynomial project, and gives a self-contained script to reproduce the checks locally.

**Scope.** The question is whether the claim is *surprising or new*. Correctness of your normalisations is only checked where the tests below touch it.

**Repo context the session should look for:**
- `finite_contour_cocycles_short.tex`, especially `lem:cocycle`, `lem:rfW` and `prop:periodint`.
- The `dr_b_periods_and_Lfunctions` draft (Cor. 5.6, the wall-crossing residue).
- Companion scripts `wip_meromorphic_periods.py` and `period_polynomial_bases.py`.

---

## 0. The claim under review

Let $f \in F_k$ have a simple (or higher-order) pole at a generic $z \in \mathbb H$.

- **Cycle.** Take a small loop $c_z$ around $z$ and form $\gamma_\circ(z) = c_z - S c_z$. This satisfies $S\gamma_\circ = -\gamma_\circ$.
- **Kernel.** Use the single kernel
$$\mathcal K = (X-\tau Y)^n + \mathbf K_T(\tau) - r_{E_k}, \qquad n = k-2,\quad \mathbf K_T = \chi - \chi|_S .$$
- **Object.** Define $O(z) := \int_{\gamma_\circ(z)} f\,\mathcal K\,d\tau$.

**Claim:** $O(z) \in W = \ker(1+S)\cap\ker(1+U+U^2)$.

This was reported as a shock.

## 1. Objection: W-membership is a formal consequence of `lem:cocycle` + `lem:rfW`

Notation for this section:
- Deform **both** segments of an admissible pair $(\gamma^S_{\tau_0},\gamma^T_{\tau_0})$ by the same cycle $\gamma_\circ$.
- Drop the $(2\pi i)^{n+1}$ prefactors.
- Let $\omega = f(X-\tau Y)^n d\tau$, and set
$$R=\int_{c_z}\omega,\qquad D=\int_{\gamma_\circ}\omega,\qquad Q_\circ=\int_{\gamma_\circ} f\chi,\qquad \lambda_\circ=\int_{\gamma_\circ} f\,d\tau .$$

The argument has four steps.

1. **The raw jump is a cocycle.** Equivariance gives $D = R - R|_S$, so $D|(1+S)=0$. The jumps in the raw cocycle are $\Delta C_S = \Delta C_T = D$, which forces $\Delta C_U = \Delta C_T|_S + \Delta C_S = D|(1+S) = 0$.
   - Geometrically, $S\gamma_\circ=-\gamma_\circ$ means the composite path for $U=TS$ is undeformed in homology.
   - Hence the relations $S^2 = U^3 = -I$ survive, and $\Delta C$ is a genuine 1-cocycle.
   - By contrast, the single-segment wall-crossings of Cor. 5.6 break the $S$-relation. The antisymmetrised pair is exactly the relation-preserving combination.
2. **The gauge-fixed jump is a cocycle.** Since $\Delta E = Q_\circ$, we have $\Delta\hat C = \Delta C - \delta Q_\circ$, which is again a cocycle. Its $S$-value is
$$\Delta\hat C_S = D + Q_\circ|(1-S) = \int_{\gamma_\circ} f\,[(X-\tau Y)^n + \mathbf K_T],$$
   i.e. $O(z)$ before the Eisenstein term.
3. **Its $T$-value is the Eisenstein obstruction.** Use $\chi(\tau-1)-\chi(\tau) = (X-\tau Y)^n$ and $\chi|_T = \chi(\tau-1) + A_0$, with $A_0$ a constant polynomial. Then
$$\Delta\hat C_T = \int_{\gamma_\circ} f\,[\chi(\tau-1) - \chi|_T] = -A_0\,\lambda_\circ .$$
   This is the same obstruction as in `lem:cocycle`, with $c_f(0)\to\lambda_\circ = 2\pi i(\mathrm{Res}_z f - \mathrm{Res}_{Sz} f)$.
4. **Subtracting Eisenstein gives a parabolic cocycle.** Subtract $(\lambda_\circ / c_{E_k}(0))\,\hat C^{E_k}$. The $T$-value of the result vanishes, and `lem:rfW` then puts its $S$-value in $W$. The "$-r_{E_k}$ inside the kernel" is exactly this subtraction, since $\lambda_\circ = \int_{\gamma_\circ} f\cdot 1\,d\tau$.

**Most elementary form.** On $\mathrm{PSL}_2(\mathbb Z) = \mathbb Z/2 * \mathbb Z/3$, the pair $(c_S, c_U) = (R|(1-S),\,0)$ is a cocycle for every $R\in V_n$. $O(z)$ is the parabolic normal form of that cocycle, corrected by Eisenstein.

**Assessment.** W-membership should not be presented as new or deep.

## 2. Refinements

- **R1 (localisation).** $O(z) = \sum \mathrm{Res}_z\big(f\cdot\tilde H\big) - \lambda_\circ\,\tilde r$, with $\tilde H(\tau) = G(\tau) - \tau^n G(-1/\tau)$ and $G = (X-\tau Y)^n + \mathbf K_T$.
  - Only the principal part of $f$ at $z$ matters.
  - The only modularity used is $\mathrm{Res}_{Sz} f = z^n\,\mathrm{Res}_z f$.
  - The claim is therefore equivalent to a Bernoulli-polynomial identity: $F(\tau) := \tilde H(\tau) - (1-\tau^n)\tilde r \in W$ for all $\tau$.
  - Higher-order poles follow automatically, because they only involve $\tau$-derivatives of $F$.
- **R2 (shape of the U-defect).** $\tilde H|(1+U+U^2) = (1-\tau^n)\,D_k$ with $D_k$ independent of $\tau$. This is the signature of step 3 above.
- **R3 (not degenerate).** For $k = 12, 16, 18$, the span of $\{F(\tau)\}$ is all of $W$ (dimension 3). So $O$ has genuine cuspidal components, not only the coboundary line $X^n - Y^n$.
- **R4 (depends on the lift).** $O(z+1) \ne O(z)$, and similarly for other $\Gamma$-translates, with the difference **not** in $\mathbb C(X^n-Y^n)$ when $\dim S_k \ge 1$. So $z \mapsto O(z)$ is a *non-canonical* map from pole data to $W$ that depends on which lift of the orbit you pick and on the choice of $S$ as antisymmetriser.
- **R5 (natural home).** Brown–Fonseca, arXiv:2508.04844, §5 (in your bibliography) give the residue exact sequence
$$0\to H^1(\mathcal Y;\mathcal V_n)\to H^1(\mathcal Y\setminus\{w\};\mathcal V_n)\xrightarrow{\mathrm{Res}}\mathcal V_n^{\Gamma_w}(-1)\to 0.$$
  - The $\Psi^{p,q}_\Gamma$ give a canonical de Rham splitting (their Cor. 5.16).
  - $O(z)$ looks like a choice-dependent Betti-side splitting/extension datum in your gauge.
  - The explicit level-1 Hurwitz/Bernoulli formula for $z\mapsto O(z)$ was not found in the literature, but the search was not exhaustive. Kohnen–Zagier Bernoulli bases of $W$ are close relatives. Neutral skepticism on novelty.

## 3. Tests

### Setup

```bash
pip install sympy mpmath numpy
python closed_contour_tests.py                 # all k in {4,...,26}; numerics for k>=12, k!=14
python closed_contour_tests.py 12 16 --nonum   # symbolic only
```

Numerics use 40-digit mpmath and take a few minutes per weight. The test form is
$$f = \frac{E_4^6\,E_{k-12}}{E_4^3 - j(z_0)\,\Delta},\qquad z_0 = 0.2+1.3i,$$
(note: $z_0\in\mathbb Q(i)$ is a CM point, see §7; irrelevant for T1–T7),
which has simple poles on $\Gamma z_0$ and $c_f(0) = 1 \neq 0$, so the Eisenstein term is exercised.

### Test list

| ID | What | Expected | Observed (k = 12, 16 run; 4, 10 symbolic) |
|---|---|---|---|
| T1a | $\tilde H\vert(1+S) \equiv 0$ | pass | pass |
| T1b | $\tilde H\vert(1+U+U^2) = (1-\tau^n) D_k$, $D_k$ constant | pass | pass |
| T2a | a particular $\tilde r$ exists and $F(\tau)\in W$ for all $\tau$ | pass | pass |
| T2b | rank of $\mathrm{span}_\tau F(\tau)$ vs $\dim W$ | equal | 3 = 3 (k = 12, 16, 18) |
| T2c | weights with no cusp forms: $F(\tau)\in\mathbb C(X^n-Y^n)$ | pass | pass (k = 4, 10) |
| T3a–c | ingredients of §1: shift identity, $A_0$ constant, $[D]_{X^n}=\lambda_\circ$ | pass | pass |
| T4a | $z^n O(Sz) = -O(z)$ (sanity) | pass | pass |
| T4b,c | lift dependence: $O(z+1)-O(z)$, $O(gz)(cz+d)^n - O(z)$ | nonzero, outside coboundary line iff $S_k\ne0$ | as expected |
| T6a | end-to-end numeric $O$ **without** $r_E$ term | $S$-relation ok, $U$-defect $\ne 0$ | $\sim 10^{3}$ |
| T6b | end-to-end numeric $O - \lambda_\circ\tilde r$ | both relations $\approx 0$ | $\lesssim 10^{-33}$ |
| T6c | numeric $O - \lambda_\circ\tilde r$ equals $2\pi i\,\mathrm{Res}_z(f)\,F(z_0)$ (localisation) | $\approx 0$ | $\lesssim 10^{-34}$ |
| T7 | $\tilde r - \hat C^{E_k}_S(\tau_0=i) \in W$, i.e. the counterterm **is** the Eisenstein $S$-value from the same recipe (per unit $c_{E_k}(0)$) | pass | pass (double precision) |

### Tasks for the Claude Code session (not yet done)

- **T8 (reconcile with repo).**
  - Replace the unit normalisation with the repo's $(2\pi i)^{n+1}$ prefactors and its $r_{E_k}$.
  - Confirm that T6b and T7 still pass with *your* $r_{E_k}$.
  - Cross-check T6 against `wip_meromorphic_periods.py` on the same form and pole.
- **T9 (antisymmetriser choice).**
  - Replace $S$ by a conjugate $gSg^{-1}$, i.e. cycle $c_z - gSg^{-1}c_z$, and repeat T6.
  - Question: which choices preserve W-membership, and how does the output change?
  - §1 predicts that any order-2 element works, with the appropriate gauge.
- **T10 (cuspidal content).**
  - For $k=12$, project $F(z)$ onto $W_\pm$ (basis in `period_polynomial_bases.py`), or pair it with $r_\Delta$ via Haberland.
  - Result: explicit Laurent polynomials in $z$ with coefficients in $\mathbb Q\,\omega^\pm(\Delta)$.
  - Question: do these match a Betti period of Brown–Fonseca's $\Psi^{p,q}_\Gamma(\cdot, z)$ (their eq. (3.2), Cor. 3.11) for some path choice?
  - If yes, $O$ is known extension data in your gauge. If nothing natural matches, that is the more interesting outcome.
- **T11 (lift dependence, quantified).**
  - Compute $O(\gamma z)(cz+d)^n - O(z)$ for generators $\gamma$ and project the result to $W/\mathbb C(X^n-Y^n)$.
  - Check whether $\gamma \mapsto$ (that difference) is itself a cocycle or coboundary in $\gamma$. This would pin down exactly what non-canonical choice $O$ encodes.

---

## 4. Overview of §§5–11

- **§5:** checks of the k = 18 draft section.
- **§6:** the Mayer–Vietoris identification.
- **§§7–9:** open research tasks (jump lattice, Haberland/BKvP, Brown–Fonseca overlap).
- **§10:** the scripts.
- **§11:** a canonical W-valued object via Brown–Fonseca's residue splitting, tested on Bengoechea's forms.

Suggested order for the session: run §10 scripts → fix the §5 pushbacks in the draft → T13 (it gates §7) → T18 (it gates everything) → T16–T17 → T19–T21.

## 5. Cross-checks of the draft section `sec:resex` (k = 18, $f_z=\Delta E_4\,J'/(J-J(z))$)

Run `python k18_crosscheck.py`. It takes about 2–3 minutes and needs `closed_contour_tests.py` alongside. Results from my run are below. The session should re-run it and then fix the draft text where marked **PUSHBACK**.

### Confirmed

- **K0.** $W^{(18)}_\pm$ as in `def:Wpm` lie in $W$, with supports $\ell$ odd and $\ell\in\{2,4,6,10,12,14\}$.
- **K3.** $\mu_-/P_- = 2/131601$ and $\mu_+/P_+ = 2/3$, exactly as in the draft. This holds even in unit normalisation.
- **K5.** $R = 43867\,P_+/P_-$, and all four special values in `eq:Rheeg`/`eq:Rheeg7` match exactly.
- **K6.** $R(-1/z)=R(z)$ and $R(-z)=-R(z)$.
- **K2 (Eisenstein structure).** Take $\tilde r$ to be the rational particular solution chosen by the script. Then $r_E-\tilde r\in W$ has:
  - a rational $W_-$ coordinate;
  - a $W_+$ coordinate equal to 0 to about 45 digits;
  - a purely imaginary transcendental $p_0$ coordinate ($\approx -0.18480027377220582\,i$ in unit normalisation, not yet identified in closed form).

### PUSHBACKS to fix in the draft

- **P1 (growth).** "$|\Xi_{f_z}|\asymp e^{-4\pi\,\mathrm{Im}\,z}$" is wrong: $\Delta\sim q$, so the decay is $e^{-2\pi\,\mathrm{Im}\,z}$ (K9). There is also polynomial growth $\sim|z|^{16}$ from $\Xi^\circ$.
- **P2 (normalisation).** "the normalisations of $W_\pm$ both cancel in the ratio" is false. Rescaling $W_+$ against $W_-$ rescales $R$.
  - The primitive-part convention fixes $R$, so it is well-defined.
  - The specific numbers ($43867/8129$, …) are convention-dependent.
  - Symmetries, the "purely imaginary" loci and the zero/pole locations are convention-independent.
- **P3 (lift dependence).** "depends only on the $\mathrm{SL}_2(\mathbb Z)$-orbit of the pole through its $S$-class" is misleading. $R(z+1)\ne R(z)$ (K6), so $R$ depends on the chosen lift and is invariant only under $z\mapsto-1/z$. Reword.
- **P4 (elliptic point).** $(1+\sqrt{-3})/2$ is elliptic (equivalent to $\rho$).
  - There $E_4=0$ (K8), and $J'/(J-J(z))$ has a simple pole, so $f_z$ is holomorphic at that point.
  - More generally, no weight-18 form has a nonzero residue polynomial at the $\rho$-orbit, since $16\not\equiv0\pmod 3$.
  - So $R$ there is a formal value of the rational function. Lead with $\sqrt{-2}$ instead, or add a remark.
- **P5 ($z=i$).** "$R$ is undefined" is imprecise. $P_\pm$ share the factor $z^2+1$ (K7), and after cancelling it $R(i)=833473\,i/89174$. $\Xi$ vanishes at $i$, but its direction has a limit.
- **P6 ($r_{E_{18}}$ is not in $W$).** $r_{E_{18}}$ has a nonzero $U$-defect (K1); that is why it works as a counterterm.
  - "the $W_+$ coordinate of $r_{E_{18}}$ vanishes…" therefore needs a convention.
  - Phrase it via the coefficients $c_0,c_n$ of $\Xi^\circ$, which do lie in $W$, or via $r_E$ minus a stated rational complement.
- **P7 (minor).** The $p_0$ coordinate of $\Xi^\circ$ is a Laurent polynomial with odd powers from $z^{-1}$ to $z^{17}$ (K4). These come from the degree-$(n+1)$ Bernoulli polynomial in $\chi$. This only matters if $\sum_m c_m z^m$ is meant to run over $0\le m\le n$.
- **P8 (observation, not an error).** $P_+$ vanishes at $z=0,\pm\tfrac12,\pm1,\pm2$, with double zeros at $\pm1$ (K7). This may be worth a sentence, or a check of whether it persists at other weights.

## 6. The object is the Mayer–Vietoris image of the residue

`mv_check.py` checks the following for k = 12, 16, 18. Let
$$\mathrm{NF}(v) := \text{normal form of the cocycle } (S\mapsto v|(1-S),\ U\mapsto 0),$$
where "normal form" means: $T$-gauge-fixed, with the Eisenstein term subtracted using the same $\tilde r$.

Then:

- **T12a.** $O(z)\equiv \mathrm{NF}(R_z)$ modulo $\mathbb C(X^n-Y^n)$, where $R_z=(X-zY)^n$ is the unit-residue polynomial.
- **T12b.** $V^S$ and $V^U$ map to the coboundary line with $\lambda=0$. Dimensions: $5+3$ at $k=12$, $7+5$ at $k=16$, $9+5$ at $k=18$. Each matches $\operatorname{codim}=\dim M_k+\dim S_k$.

This is the classical isomorphism $V/(V^S+V^U)\cong H^1(\mathrm{PSL}_2(\mathbb Z),V)$ for $\mathbb Z/2*\mathbb Z/3$ (Serre, *Trees*).

**Task.** Re-run, and check the sign and normalisation against the repo conventions.

## 7. Item 1: are the jumps discrete?

**Question.** Moving contours around different copies $\gamma z$ of the pole gives jumps. Do all the jumps together form a lattice (so "period polynomial mod jumps" is meaningful, like $\log$ mod $2\pi i$) or a dense set?

**Mechanism.** The jump from lift $\gamma z$ has $W_\pm$ coordinates $g(z)\,(cz+d)^{16}\mu_\pm(\gamma z)$. That is a polynomial of degree $\le16$ in $z$ with rational coefficients of bounded denominator. Hence:

- **CM case.** If $z$ is imaginary quadratic, all values lie in $\mathbb Q(z)^2$, which has $\mathbb Q$-dimension 4. So the jumps form a lattice.
- **Transcendental case.** If $z$ is transcendental and the jump polynomials have $\mathbb Q$-rank $>4$, the jumps cannot be discrete in $\mathbb C^2\cong\mathbb R^4$.
- The $p_0$ (coboundary) direction is set aside.

**Results (`lattice_test.py`).**

- **(A)** 36 lifts (words of length ≤ 4 in $S,T,T^{-1}$) give jump polynomials of $\mathbb Q$-rank **15** > 4. So the jumps are non-discrete for transcendental $z$.
- **(B)** At $z=\sqrt{-2}$ the $(a,b)$-vectors have $\mathbb Q$-rank 4 and common denominator 43867, so they form a lattice.
- **(B′)** At $z=(2+13i)/10$ they also form a lattice (rank 4), but a very fine one: common denominator $\approx 8.2\times10^{18}$.
  - **Note:** $0.2+1.3i$, the "generic" test point used in §§1–3, is actually a CM point in $\mathbb Q(i)$. The W-membership tests don't care, but the label "generic" was wrong.

**Caveat (must be tested).** All of this assumes that each lift's jump is a genuine ambiguity of a W-valued meromorphic period polynomial. That assumption comes from the "sections of $\pi_1^{\rm orb}(Y\setminus\{z\})\to\Gamma$" argument, which is unverified.

**Tasks.**

- **T13.** For $f_z$ at k = 18:
  1. Construct two admissible contour systems that differ by passing on opposite sides of a *translated* pole copy $\gamma z$, e.g. $\gamma=T$ or $ST$.
  2. Apply the same $T$-gauge fix and Eisenstein subtraction to both.
  3. Check that both results lie in $W$, and that their difference equals the predicted jump $g(z)(cz+d)^{16}\,(\mu_-,\mu_+)(\gamma z)$ modulo $p_0$.
  4. If the difference is not in $W$, the section argument needs repair before anything else in this section is used.
- **T14.** Extend (A) to longer words and to k = 22, 26 ($\dim S_k=1$, so $W^\pm$ are again lines), and record the $\mathbb Q$-rank.
  - Try a genuinely non-quadratic point such as $z=0.2+\tfrac{\pi}{2}i$.
  - As a numerical illustration, run an integer-relation / LLL search showing ever-smaller nonzero combinations as the precision increases. My quick LLL attempt was dominated by rounding noise, so do this carefully (mpmath `pslq`, or `fpylll` at high precision).
- **T15.** At CM points, compute the lattice (covolume, discriminant) for $\sqrt{-2}$, $(1+\sqrt{-7})/2$ and $(1+\sqrt{-11})/2$, and see whether anything arithmetic appears beyond the $B_{18}$ numerator.
  - Note that the overall scalar $g(z)=\Delta(z)E_4(z)$ multiplies the lattice.

## 8. Item 2: pair with $\Delta$ (Haberland) and compare with regularized inner products

**Idea.** Pairing a W-valued period polynomial of meromorphic $f$ with the period polynomial of a cusp form $g$ should give something related to a regularized Petersson product $\langle f,g\rangle$, well-defined only modulo the pairing of $g$ against the jump set of §7.

**Literature to locate.** Bringmann–Kane–von Pippich, *Regularized inner products of meromorphic modular forms and higher Green's functions*. Brown–Fonseca cite it as [BKvP19]; verify the bibliographic details. Extract:

- their regularisation;
- their normalisation;
- whether their product is defined for forms like $f_z$ (a cusp form at $\infty$ with simple poles at a non-elliptic orbit).

**Tasks.**

- **T16.** Implement the Haberland pairing exactly as in the draft (`finite_contour_cocycles_short.tex`).
  1. Verify it on holomorphic cusp forms against a direct numerical Petersson product; this also fixes the normalisation. Use k = 12 ($\Delta$) or k = 18 ($\Delta E_6$).
  2. Pair the meromorphic W-valued period polynomial of $f_z$ (from T13) with $r_{\Delta E_6}$ at k = 18.
  3. Compare with the BKvP regularized product, numerically if their formula is computable.
- **T17.** Compute the pairing of $r_{\Delta E_6}$ with the jump lattice at $z=\sqrt{-2}$. This gives a lattice in $\mathbb C$, and any match with BKvP should hold modulo it.
- **Stopping rule.** If BKvP's product is not defined for these $f$, or needs a different class of $g$, say so and stop.

## 9. Item 3: check the Brown–Fonseca overlap before investing further

**Source.** Brown–Fonseca, arXiv:2508.04844. Only §§1–5.3 have been read so far.

**What they do (per those sections).**

- A residue exact sequence for $H^1$ of the punctured modular curve.
- A canonical *de Rham* splitting via the meromorphic forms $\Psi^{p,q}_\Gamma$.
- *Single-valued* periods (Green's functions) of the relative biextension $H^1_{\rm cusp}(Y\setminus\{w\},\{z\})$.

The present project looks like the *multivalued Betti* side: jumps are the analogue of $2\pi i$ for $\log x$, and BF's objects are the analogue of $\log|x|^2$.

**Tasks.**

- **T18.** Skim BF §§6–8 (Betti realisation, biextension, periods) and answer:
  1. Do they give an explicit Betti/group-cohomology description of $H^1(Y\setminus\{w\};V)$ or its splittings? Any Mayer–Vietoris, Serre-tree or $\mathbb Z/2*\mathbb Z/3$ presentation?
  2. Do they compute the (non-single-valued) Betti periods of $\Psi^{p,q}_\Gamma(\cdot,w)$, or the lattice of loop-periods around punctures? If so, compare directly with $\Xi_{f}$ and §7.
  3. Remark 7.6 mentions meromorphic modular forms built as sums of negative powers of binary quadratic forms of negative discriminant, which have poles exactly at CM points (compare Zagier's $f_{k,D}$ for $D<0$ and Bengoechea's work on them; verify references). These are natural test forms for §7. Does anything there bear on discreteness at CM points?
  4. Anything on Hecke equivariance of these constructions?
- **Report.** A short verdict: same object / Betti counterpart of their object / unrelated. Include section and equation references.

## 10. Scripts

Put all four scripts in one directory and run them in this order:

1. `closed_contour_tests.py`: §§1–3 tests T1–T7.
2. `mv_check.py`: §6.
3. `k18_crosscheck.py`: §5.
4. `lattice_test.py`: §7.

Dependencies: `pip install sympy mpmath numpy`.

### `closed_contour_tests.py`

```python
"""
Tests for the S-antisymmetrised closed-contour object O(z) and the claim O(z) in W.
Conventions (match finite_contour_cocycles_short.tex, Lemma lem:cocycle):
  V_n = C[X,Y]_n, right action (X,Y)|g = (aX+bY, cX+dY)
  S=[[0,-1],[1,0]], T=[[1,1],[0,1]], U=TS=[[1,-1],[1,0]]
  chi(tau) = sum_l C(n,l)(-1)^l zeta(-l,tau+1) X^{n-l} Y^l,  zeta(-l,a) = -B_{l+1}(a)/(l+1)
  K_T = chi - chi|S,   G(tau) = (X - tau Y)^n + K_T(tau)
All (2 pi i)^{n+1} prefactors are DROPPED (unit normalisation). Reconcile before comparing
with repo numbers.
"""
import sys
import sympy as sp
import mpmath as mp

X, Y, t = sp.symbols('X Y tau')
S = (0, -1, 1, 0); T = (1, 1, 0, 1); U = (1, -1, 1, 0)
_A, _B = sp.symbols('A_ B_')

def slash(P, M):
    a, b, c, d = M
    return sp.expand(P.subs({X: _A, Y: _B}, simultaneous=True)
                      .subs({_A: a*X + b*Y, _B: c*X + d*Y}, simultaneous=True))

def coeffs(P, n):
    P = sp.Poly(sp.expand(P), X, Y)
    return [P.coeff_monomial(X**(n-i)*Y**i) for i in range(n+1)]

def from_coeffs(v, n):
    return sum(v[i]*X**(n-i)*Y**i for i in range(n+1))

def chi(n, tt):
    return sum(sp.binomial(n, l)*(-1)**l*(-sp.bernoulli(l+1, tt+1)/(l+1))*X**(n-l)*Y**l
               for l in range(n+1))

def G(n, tt):
    c = chi(n, tt)
    return sp.expand((X - tt*Y)**n + c - slash(c, S))

def Htilde(n, tt):
    # S-antisymmetrised kernel: integrating f*G over c_z - S c_z  ==  Res_z of f*Htilde
    return sp.expand(sp.cancel(sp.together(G(n, tt) - tt**n*G(n, -1/tt))))

def relS(P):  return sp.expand(P + slash(P, S))
def relU(P):  return sp.expand(P + slash(P, U) + slash(slash(P, U), U))

def particular_r(n, Dk):
    cs = sp.symbols('c0:%d' % (n+1)); r = from_coeffs(cs, n)
    eqs = coeffs(relS(r), n) + coeffs(relU(r) - Dk, n)
    sol = sp.solve(eqs, cs, dict=True)
    if not sol: return None
    return sp.expand(r.subs(sol[0]).subs({c: 0 for c in cs}))

def W_basis(n):
    cs = sp.symbols('c0:%d' % (n+1)); r = from_coeffs(cs, n)
    eqs = coeffs(relS(r), n) + coeffs(relU(r), n)
    M = sp.Matrix([[sp.diff(e, c) for c in cs] for e in eqs])
    return [from_coeffs(list(v), n) for v in M.nullspace()]

def in_span(P, basis, n):
    if sp.expand(P) == 0: return True
    M = sp.Matrix([coeffs(b, n) for b in basis]).T
    return M.rank() == M.row_join(sp.Matrix(coeffs(P, n))).rank()

PASS = {True: 'PASS', False: 'FAIL'}

def symbolic_tests(k):
    n = k - 2
    print(f'\n===== k = {k} (n = {n}) =====')
    H = Htilde(n, t)
    # T1: S-relation identically; U-defect has shape (1 - tau^n) * D_k
    t1a = relS(H) == 0
    uq = sp.cancel(relU(H)/(1 - t**n))
    t1b = sp.simplify(sp.diff(uq, t)) == 0
    print(f'[T1a] Htilde|(1+S) == 0                         {PASS[t1a]}')
    print(f'[T1b] Htilde|(1+U+U^2) = (1-tau^n) D_k, D_k const {PASS[t1b]}')
    Dk = sp.expand(uq)
    # T3: cocycle-derivation ingredients
    c0 = chi(n, t)
    t3a = sp.expand(chi(n, t-1) - c0 - (X - t*Y)**n) == 0
    A0 = sp.expand(slash(c0, T) - chi(n, t-1))
    t3b = sp.diff(A0, t) == 0
    print(f'[T3a] chi(tau-1)-chi(tau) == (X-tau Y)^n            {PASS[t3a]}')
    print(f'[T3b] A_0 := chi|T - chi(tau-1) is tau-independent  {PASS[t3b]}   A_0 = {sp.factor(A0)}')
    # lambda_circ per unit residue: X^n-coefficient of the pure (X-tau Y)^n|(1-S) part
    Dpart = sp.expand((X - t*Y)**n - slash((X - t*Y)**n, S))
    t3c = sp.simplify(coeffs(Dpart, n)[0] - (1 - t**n)) == 0
    print(f'[T3c] [ (X-tauY)^n|(1-S) ]_{{X^n}} == 1 - tau^n        {PASS[t3c]}')
    # T2: particular r~, W-membership of F(tau), span
    r0 = particular_r(n, Dk)
    if r0 is None:
        print('[T2]  no r~ solves r|(1+S)=0, r|(1+U+U^2)=D_k       FAIL'); return None
    F = sp.expand(sp.cancel(H - (1 - t**n)*r0))
    t2a = relS(F) == 0 and relU(F) == 0
    print(f'[T2a] F(tau) := Htilde - (1-tau^n) r~ lies in W      {PASS[t2a]}')
    Wb = W_basis(n)
    Fm = sp.Poly(sp.expand(F*t**(n+2)), t)
    M = sp.Matrix([coeffs(c, n) for c in Fm.all_coeffs()])
    print(f'[T2b] dim W = {len(Wb)};  rank span_tau F(tau) = {M.rank()}')
    cob = X**n - Y**n
    if len(Wb) == 1:
        t2c = all(in_span(sp.expand(c), [cob], n) for c in Fm.all_coeffs())
        print(f'[T2c] (no cusp forms) F(tau) in C*(X^n - Y^n)         {PASS[t2c]}')
    # T4: lift dependence. Res_{gz} f = (cz+d)^n Res_z f  =>  O(gz)*(cz+d)^n vs O(z)
    O = lambda z: sp.expand(sp.cancel(F.subs(t, z)))
    t4s = sp.simplify(sp.cancel(t**n*O(-1/t) + O(t))) == 0
    print(f'[T4a] antisymmetry: z^n O(Sz) == -O(z)                {PASS[t4s]}')
    dT = sp.expand(sp.cancel(O(t+1) - O(t)))
    print(f'[T4b] O(z+1) - O(z) == 0 ?  {dT == 0};   in C*(X^n-Y^n)? {in_span(dT,[cob],n) if dT!=0 else True}')
    g = (1, 0, 1, 1)  # tau -> tau/(tau+1)
    dg = sp.expand(sp.cancel((t+1)**n*O(t/(t+1)) - O(t)))
    print(f'[T4c] O(gz)(z+1)^n - O(z), g=[[1,0],[1,1]]: zero? {dg == 0};  in C*(X^n-Y^n)? {in_span(dg,[cob],n) if dg!=0 else True}')
    return dict(n=n, r0=r0, F=F, A0=A0, Dk=Dk, Wb=Wb)

# ---------------- numerics ----------------
def Ek_num(k, tau, N=60):
    q = mp.e**(2j*mp.pi*tau)
    Bk = mp.bernoulli(k)
    s = mp.mpf(0)
    for m in range(1, N):
        s += mp.mpf(sum(d**(k-1) for d in range(1, m+1) if m % d == 0))*q**m
    return 1 - (2*k/Bk)*s

def Delta_num(tau, N=60):
    q = mp.e**(2j*mp.pi*tau)
    p = mp.mpc(1)
    for m in range(1, N): p *= (1 - q**m)**24
    return q*p

def lambdify_V(P, n):
    fs = [sp.lambdify(t, c, 'mpmath') for c in coeffs(P, n)]
    return lambda z: [f(z) if callable(f) else f for f in fs]

def to_mp(c):
    c = sp.nsimplify(c)
    re_, im_ = sp.re(c), sp.im(c)
    q = lambda r: mp.mpf(sp.Rational(r).p)/mp.mpf(sp.Rational(r).q)
    return mp.mpc(q(re_), q(im_))

def slash_matrix(M, n):
    # column j = coefficients of (basis_j)|M
    cols = [coeffs(slash(X**(n-j)*Y**j, M), n) for j in range(n+1)]
    return mp.matrix([[to_mp(cols[j][i]) for j in range(n+1)] for i in range(n+1)])

def rel_residual(v, n):
    Sm, Um = slash_matrix(S, n), slash_matrix(U, n)
    I = mp.eye(n+1); x = mp.matrix(v)
    a = (I + Sm)*x; b = (I + Um + Um*Um)*x
    return float(max(abs(e) for e in a)), float(max(abs(e) for e in b))

def numeric_tests(k, data, z0=mp.mpc(0.2, 1.3), rad=mp.mpf('0.05'), M=256):
    mp.mp.dps = 40
    n = data['n']
    print(f'\n----- numerics, k = {k} -----')
    Gvec = lambdify_V(G(n, t), n)
    # T7: Eisenstein normalisation.  C^E_S at tau0 = i:  int_{i-1}^{i} E_k(tau) K_T(tau) dtau
    KT = lambdify_V(sp.expand(chi(n, t) - slash(chi(n, t), S)), n)
    CE = [mp.quad(lambda s, i=i: Ek_num(k, mp.mpc(s, 1))*KT(mp.mpc(s, 1))[i], [-1, 0]) for i in range(n+1)]
    CEp = from_coeffs([sp.Float(mp.re(x), 35) + sp.I*sp.Float(mp.im(x), 35) for x in CE], n)
    diff = sp.expand(data['r0'] - CEp)
    # is r~ - C^E_S in W (numerically)?
    Wb = data['Wb']
    Mw = sp.Matrix([[complex(c) for c in coeffs(b, n)] for b in Wb]).T
    import numpy as np
    A = np.array(Mw.tolist(), dtype=complex); bvec = np.array([complex(c) for c in coeffs(diff, n)])
    sol, res, *_ = np.linalg.lstsq(A, bvec, rcond=None)
    resid = np.linalg.norm(A@sol - bvec)/max(1e-300, np.linalg.norm(bvec) + 1)
    print(f'[T7]  r~ - C^{{E_k}}_S(tau0=i) in W   relative residual = {resid:.2e}  ({PASS[resid < 1e-12]})')
    # T6: end-to-end with a genuine meromorphic form of weight k, simple poles on Gamma*z0:
    #   f = E_4^6 * E_{k-12} / (E_4^3 - j(z0) Delta)      (E_0 := 1; k-12 must not be 2)
    j0 = Ek_num(4, z0)**3/Delta_num(z0)
    def f(tau):
        num = Ek_num(4, tau)**6 * (Ek_num(k-12, tau) if k > 12 else 1)
        return num/(Ek_num(4, tau)**3 - j0*Delta_num(tau))
    def loop(center, h):
        out = [mp.mpc(0)]*(n+1); scal = mp.mpc(0)
        for m in range(M):
            th = 2*mp.pi*m/M; w = center + rad*mp.e**(1j*th)
            dw = 1j*rad*mp.e**(1j*th)*(2*mp.pi/M)
            fv = f(w); gv = h(w)
            out = [out[i] + fv*gv[i]*dw for i in range(n+1)]
            scal += fv*dw
        return out, scal
    Sz0 = -1/z0
    I1, l1 = loop(z0, Gvec); I2, l2 = loop(Sz0, Gvec)
    lam = l1 - l2
    Ovec = [I1[i] - I2[i] for i in range(n+1)]
    r0v = [to_mp(c) for c in coeffs(data['r0'], n)]
    Ocorr = [Ovec[i] - lam*r0v[i] for i in range(n+1)]
    s1, u1 = rel_residual(Ovec, n); s2, u2 = rel_residual(Ocorr, n)
    scale = max(abs(x) for x in Ocorr)
    print(f'[T6a] raw O (no Eisenstein term):  |.(1+S)|={s1:.2e}  |.(1+U+U^2)|={u1:.2e}  (U-defect expected nonzero)')
    print(f'[T6b] O - lambda r~            :  |.(1+S)|={s2:.2e}  |.(1+U+U^2)|={u2:.2e}  scale={float(scale):.2e}')
    # localisation: compare with 2 pi i Res_z f * F(z0)
    res = l1/(2j*mp.pi)
    Fz = lambdify_V(data['F'], n)(z0)
    err = max(abs(Ocorr[i] - 2j*mp.pi*res*Fz[i]) for i in range(n+1))
    print(f'[T6c] O - lambda r~  ==  2 pi i Res_z(f) F(z0)   max err = {float(err):.2e}')

if __name__ == '__main__':
    ks = [int(a) for a in sys.argv[1:] if a.isdigit()] or [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 26]
    for k in ks:
        d = symbolic_tests(k)
        if d is not None and k >= 12 and k != 14 and '--nonum' not in sys.argv:
            numeric_tests(k, d)
```

### `mv_check.py`

```python
import sympy as sp
from closed_contour_tests import *
def nullspace_polys(M, n):
    return M.nullspace()
def check(k):
    n=k-2
    H=Htilde(n,t); Dk=sp.expand(sp.cancel(relU(H)/(1-t**n))); r0=particular_r(n,Dk)
    F=sp.expand(sp.cancel(H-(1-t**n)*r0))
    A0=sp.expand(slash(chi(n,t),T)-chi(n,t-1))
    cob=X**n-Y**n
    def NF(v):
        cS=sp.expand(v-slash(v,S)); cT=cS
        lam=coeffs(cT,n)[0]
        rhs=sp.expand(cT+lam*A0)
        cs=sp.symbols('p0:%d'%(n+1)); P=from_coeffs(cs,n)
        sol=sp.solve(coeffs(sp.expand(slash(P,T)-P-rhs),n),cs,dict=True)[0]
        P=sp.expand(P.subs(sol).subs({c:0 for c in cs}))
        out=sp.expand(cS-lam*r0-(slash(P,S)-P))
        return out, lam
    z=sp.Rational(3,10)+sp.Rational(7,5)*sp.I
    R=sp.expand((X-z*Y)**n)
    nf,lam=NF(R)
    Oz=sp.expand(F.subs(t,z))
    print(k,'NF(R) in W:', relS(nf)==0 and relU(nf)==0,
          '| O(z)-NF(R) in line:', in_span(sp.expand(Oz-nf),[cob],n),
          '| O(z)+NF(R) in line:', in_span(sp.expand(Oz+nf),[cob],n))
    # kernel check: V^S + V^U  -> line, lam=0
    for g,name in [(S,'S'),(U,'U')]:
        cs=sp.symbols('c0:%d'%(n+1)); v=from_coeffs(cs,n)
        M=sp.Matrix([[sp.diff(e,c) for c in cs] for e in coeffs(sp.expand(slash(v,g)-v),n)])
        ok=True
        for b in M.nullspace():
            w,l=NF(from_coeffs(list(b),n)); ok &= in_span(w,[cob],n) and l==0
        print('   V^%s maps to line with lambda=0:'%name, ok, ' dim V^%s ='%name, len(M.nullspace()))
for k in [12,16,18]: check(k)
```

### `k18_crosscheck.py`

```python
"""
Cross-checks of the draft section 'An example localised in W at weight k=18' (sec:resex).
Unit normalisation (no (2 pi i) factors); reference form E_18 (constant term 1).
Needs closed_contour_tests.py in the same directory.
"""
import sympy as sp, mpmath as mp
from closed_contour_tests import (X, Y, t, S, U, slash, coeffs, from_coeffs, chi, Htilde,
                                  relS, relU, particular_r, lambdify_V, Ek_num, Delta_num, to_mp)
k, n = 18, 16
mp.mp.dps = 40
H = Htilde(n, t)
Dk = sp.expand(sp.cancel(relU(H)/(1 - t**n))); r0 = particular_r(n, Dk)
F = sp.expand(sp.cancel(H - (1 - t**n)*r0))

# ---- draft's basis (Definition def:Wpm) ----
x = sp.symbols('x')
def B0(m): return sum(sp.binomial(m, i)*sp.bernoulli(i)*x**(m-i) for i in range(m+1) if i != 1)
def refl(p, s): return sp.expand(p + s*x**n*p.subs(x, 1/x))
def pp(p):
    cs = sp.Poly(sp.expand(p), x).all_coeffs()
    den = sp.ilcm(*[sp.fraction(sp.Rational(c))[1] for c in cs]); q = [sp.Rational(c)*den for c in cs]
    g = sp.igcd(*[int(c) for c in q if c != 0]); return sp.Poly([c/g for c in q], x).as_expr()
p0 = x**n - 1
Wm = pp(sp.Rational(1, n-1)*refl(B0(n-1), 1) - sp.Rational(1, 3)*refl(B0(3), 1))
tmp = sp.expand(sp.Rational(1, 2)*refl(B0(2), -1) + sp.Rational(1, n)*refl(B0(n), -1))
Wp = pp(tmp - sp.Poly(tmp, x).coeff_monomial(x**n)*p0)
hom = lambda p: sp.expand(Y**n*p.subs(x, X/Y))
cm, cp, c0 = coeffs(hom(Wm), n), coeffs(hom(Wp), n), coeffs(hom(p0), n)
print('[K0] W-, W+ in W:', relS(hom(Wm)) == 0 and relU(hom(Wm)) == 0, relS(hom(Wp)) == 0 and relU(hom(Wp)) == 0)
print('     supports (Y-power l): W-', [i for i in range(n+1) if cm[i]], ' W+', [i for i in range(n+1) if cp[i]])

# ---- Eisenstein S-value from the finite-contour recipe at tau0 = i ----
KT = lambdify_V(sp.expand(chi(n, t) - slash(chi(n, t), S)), n)
CE = [mp.quad(lambda s, i=i: Ek_num(k, mp.mpc(s, 1))*KT(mp.mpc(s, 1))[i], [-1, 0]) for i in range(n+1)]
CEp = from_coeffs([sp.Float(mp.re(v), 35) + sp.I*sp.Float(mp.im(v), 35) for v in CE], n)
uE = max(abs(complex(c)) for c in coeffs(relU(CEp), n))
print(f'[K1] r_E itself is NOT in W: max|r_E|(1+U+U^2)| = {uE:.3e}  -> "coordinates of r_E" need a convention')
w = [CE[i] - to_mp(coeffs(r0, n)[i]) for i in range(n+1)]
Bm = mp.matrix([[to_mp(cm[i]), to_mp(cp[i]), to_mp(c0[i])] for i in range(n+1)])
sol = mp.lu_solve(Bm.T*Bm, Bm.T*mp.matrix(w))
print(f'[K2] (r_E - r~) coords: W- = {mp.nstr(sol[0], 20)}, W+ = {mp.nstr(sol[1], 5)}, p0 = {mp.nstr(sol[2], 20)}')
wm = sp.nsimplify(sp.Float(mp.re(sol[0]), 35), rational=True, tolerance=1e-28)

# ---- mu_± and P_± ----
Fc = coeffs(F, n)
fm, fp, f0 = sp.cancel(Fc[1]/cm[1]), sp.cancel(Fc[2]/cp[2]), sp.cancel(Fc[0]/c0[0])
assert all(sp.expand(sp.cancel(Fc[i] - fm*cm[i] - fp*cp[i] - f0*c0[i])) == 0 for i in range(n+1))
mum, mup = sp.expand(fm - (1 - t**n)*wm), sp.expand(fp)
Pm = -18000+350936*t**2-1096675*t**4+1140542*t**6-1140542*t**10+1096675*t**12-350936*t**14+18000*t**16
Pp = -24*t+154*t**3-273*t**5+143*t**7+143*t**9-273*t**11+154*t**13-24*t**15
print('[K3] mu_-/P_- =', sp.simplify(mum/Pm), '   mu_+/P_+ =', sp.simplify(mup/Pp), '   (draft: 2/131601, 2/3)')
print('[K4] p0-coordinate of Xi° has Laurent terms: powers', sorted(m[0] for m in sp.Poly(sp.expand(f0*t), t).monoms()), '(shifted by -1)')
R = sp.cancel(mup/mum)
print('[K5] R == 43867 P+/P- :', sp.simplify(R - 43867*Pp/Pm) == 0)
for zz, want in [((1+sp.sqrt(-3))/2, sp.Rational(43867, 8129)*sp.sqrt(-3)),
                 (sp.sqrt(-2), sp.Rational(43867*18, 119513)*sp.sqrt(-2)),
                 ((1+sp.sqrt(-7))/2, 43867*(-107+316769*sp.sqrt(-7))/sp.Integer(3931333248)),
                 ((1+sp.sqrt(-11))/2, 43867*(-30607+4873558*sp.sqrt(-11))/sp.Integer(75795816699))]:
    print(f'     R({zz}) matches draft:', sp.simplify(sp.radsimp(R.subs(t, zz)) - want) == 0)
print('[K6] R(-1/z)=R(z):', sp.simplify(R.subs(t, -1/t) - R) == 0,
      ' R(-z)=-R(z):', sp.simplify(R.subs(t, -t) + R) == 0,
      ' R(z+1)=R(z):', sp.simplify(R.subs(t, t+1) - R) == 0, '  <- lift dependence')
print('[K7] factors: P- =', sp.factor(Pm), '\n               P+ =', sp.factor(Pp))
print('     R(i) after cancelling (z^2+1):', sp.simplify(R.subs(t, sp.I)), '  (finite limit, not undefined)')
zr = (1 + mp.sqrt(-3))/2
print(f'[K8] (1+sqrt-3)/2 is elliptic: |E_4| there = {mp.nstr(abs(Ek_num(4, zr)), 3)} -> f_z has no pole there')
print('[K9] |Delta(x+iy)| e^{2 pi y} for y=1,2,3:',
      [mp.nstr(abs(Delta_num(mp.mpc(0.1, y)))*mp.e**(2*mp.pi*y), 6) for y in (1, 2, 3)],
      ' -> decay is e^{-2 pi Im z}, not e^{-4 pi Im z}')
```

### `lattice_test.py`

```python
"""
Item 1: is the set of W-jumps from different copies (lifts) of the pole discrete?
k = 18, using the draft's mu_-(z) ∝ P_-(z), mu_+(z) ∝ P_+(z) (overall scalars irrelevant here).
Jump from the lift gamma*z (gamma = [[a,b],[c,d]]) has W± coordinates
    g(z) * (c z + d)^16 * mu_±(gamma z)   -- a polynomial in z of degree <= 16 with Q-coefficients.
(A) exact: Q-rank of these polynomial pairs over many gamma.  Rank > 4  =>  for transcendental z
    the Z-span inside C^2 = R^4 cannot be discrete.  For imaginary-quadratic z all values lie in
    Q(z)^2 (Q-dim 4), so the span is automatically a lattice (bounded denominators).
(B) numeric illustration via LLL: shortest nonzero integer combination, CM vs generic z.
CAVEAT: assumes each lift's jump is a genuine ambiguity of the meromorphic period polynomial.
"""
import itertools, sympy as sp
z = sp.symbols('z')
Pm = -18000+350936*z**2-1096675*z**4+1140542*z**6-1140542*z**10+1096675*z**12-350936*z**14+18000*z**16
Pp = -24*z+154*z**3-273*z**5+143*z**7+143*z**9-273*z**11+154*z**13-24*z**15
mum, mup = sp.Rational(2,131601)*Pm, sp.Rational(2,3)*Pp
Sm, Tm, Ti = sp.Matrix([[0,-1],[1,0]]), sp.Matrix([[1,1],[0,1]]), sp.Matrix([[1,-1],[0,1]])
gens = {'S': Sm, 'T': Tm, 't': Ti}
def words(L):
    out = {}
    for l in range(L+1):
        for w in itertools.product('STt', repeat=l):
            if any(w[i]+w[i+1] in ('SS','Tt','tT') for i in range(len(w)-1)): continue
            M = sp.eye(2)
            for ch in w: M = M*gens[ch]
            key = tuple(M) if M[1,0] > 0 or (M[1,0] == 0 and M[1,1] > 0) else tuple(-M)
            out.setdefault(key, ''.join(w) or 'id')
    return out
CM = [sp.Poly(mum, z).coeff_monomial(z**m) for m in range(17)]
CP = [sp.Poly(mup, z).coeff_monomial(z**m) for m in range(17)]
def jump(M):
    # (cz+d)^16 * mu(gamma z) = sum_m c_m (az+b)^m (cz+d)^(16-m)
    a, b, c, d = M
    u, v = sp.Poly(a*z+b, z), sp.Poly(c*z+d, z)
    pw_u = [sp.Poly(1, z)]; pw_v = [sp.Poly(1, z)]
    for _ in range(16): pw_u.append(pw_u[-1]*u); pw_v.append(pw_v[-1]*v)
    return [sum((cc*pw_u[m]*pw_v[16-m] for m, cc in enumerate(C) if cc != 0), sp.Poly(0, z)).as_expr()
            for C in (CM, CP)]
W = words(4)
vecs = []
for key in W:
    jm, jp = jump(key)
    vecs.append([sp.Poly(jm, z).coeff_monomial(z**i) for i in range(17)] +
                [sp.Poly(jp, z).coeff_monomial(z**i) for i in range(17)])
rk = sp.Matrix(vecs).rank()
print(f'(A) {len(vecs)} lifts; Q-rank of jump polynomials = {rk}  (> 4 => non-discrete for transcendental z)')

def cm_check(minpoly, name):
    # reduce each jump polynomial mod the minimal polynomial of z: value = a + b z with a, b rational
    mp_ = sp.Poly(minpoly, z)
    rows = []
    for key in W:
        row = []
        for j in jump(key):
            r = sp.Poly(j, z).rem(mp_)
            row += [r.coeff_monomial(1), r.coeff_monomial(z)]
        rows.append(row)
    M = sp.Matrix(rows)
    den = sp.ilcm(*[sp.fraction(x)[1] for x in rows for x in x])
    print(f'(B) {name:>22}: Q-rank of (a,b)-vectors = {M.rank()} (<= 4 => lattice); common denominator = {den}')
cm_check(z**2 + 2, 'z = sqrt(-2)')
cm_check(100*z**2 - 40*z + 173, 'z = (2+13i)/10 [in Q(i)]')
```

---

## 11. A canonical W-valued object for meromorphic forms: remove the residues first

### Idea (plain language)

The open-contour period polynomial of a meromorphic form is ambiguous because every pole crossing adds an $O$-type jump, and those jumps are built from the residue polynomials. So: subtract from $f$ a *canonical* meromorphic form with the same residue polynomials, then take the period polynomial of what's left.

- **Where the canonical subtraction comes from.** Brown–Fonseca (arXiv:2508.04844) give canonical meromorphic forms $\Psi^{p,q}_\Gamma(\tau,w)$ that split off the residue part (de Rham side; their Cor. 5.16, as read earlier; **verify**).
- **What we hope to get.** If this works, you get a canonical W-valued period polynomial for meromorphic forms, and canonical critical L-values $L^*(f^\circ,m)$, $m=1,\dots,k-1$. The L-function at non-integer $s$ stays homotopy-dependent, because $\tau^{s-1}f^\circ$ still has residues.

### Key simplification (checked, `residue_injectivity.py`)

At a pole of order $m$, the residue polynomial determines the principal part **iff $m \le n+1 = k-1$**. Checked for $n = 10, 16$: rank $m$ for $m\le n+1$, and rank $n+1$ at $m=n+2$.

So if all poles of $f$ have order $\le k-1$ and $f^\circ := f - (\Psi\text{-part})$ has zero residue polynomials, then **$f^\circ$ has no poles in $\mathbb H$ at all**. Given suitable behaviour at the cusp, $f^\circ \in M_k$, and its canonical W-element is just the classical period polynomial of a holomorphic form.

So in this regime, the "canonical meromorphic period polynomial" is really a **canonical projection $F_k \to M_k$**, $f \mapsto f^\circ$. The interesting object is that projection: which one is it? A natural guess is the regularized Petersson projection onto cusp forms (Bringmann–Kane–von Pippich; cf. Bringmann–Kane on Petersson's problem for weight-0 meromorphic forms). **Verify.**

For pole orders $> k-1$, $f^\circ$ can keep poles, and it is only weakly-holomorphic-plus-exact in cohomology. Treat that case separately.

### Prediction to test with Bengoechea's forms

For $f_{6,-7}$ (weight 12, order-6 poles at the orbit of $z_0=(-1+\sqrt{-7})/2$), a direct computation in the originating chat found (exact residue formula, confirmed by contour integration of her series truncated at $a\le400$ to $8\times10^{-11}$; closed-contour $O$ has $W_\mp$ coordinates $\propto\sqrt{-7}$ with rational ratio $31095/14$)
$$R = -\tfrac{36\sqrt{-7}}{7^5}\,Q_0(X,Y)^5,\qquad Q_0=(X-z_0Y)(X-\bar z_0Y)=X^2+XY+2Y^2 .$$

- **Which $\Psi$.** If BF's $\Psi^{p,q}(\tau,w)$ has residue polynomial $\propto (X-\bar wY)^p(X-wY)^q$, as recalled from their Cor. 3.11 (**verify** convention and normalisation), then a single $\Psi^{5,5}(\tau,z_0)$ matches $R$.
- **Prediction.** $f^\circ = f_{6,-7} - c\,\Psi^{5,5}(\cdot,z_0) = \alpha E_{12} + \beta\Delta$. The conjecture is:
  - $\beta = c_D$, the cusp-form coefficient in Bengoechea's decomposition (meromorphic part with algebraic Fourier coefficients $+\,c_D\Delta$);
  - $\alpha$ is either $0$ or something explicit.

  Equivalently: BF's $\Psi$ at a CM point is Bengoechea's algebraic part. BF's Remark 7.6 mentions exactly such sums over quadratic forms.
- **Caveats.**
  - BF's $\Psi$ may depend real-analytically on $w$, since their residues involve $\bar w$. Canonicity then costs holomorphic dependence on the pole position (the $\log|x|^2$ side of the analogy).
  - If $\Psi^{5,5}(\cdot,z_0)$ turns out to be exactly proportional to $f_{6,-7}$ (both are elliptic-Poincaré-type sums), then $f^\circ=0$ and the prediction fails as stated. Report this. It would mean BF's canonical choice differs from Bengoechea's split.

### Tasks

- **T19 (read).** From BF, extract:
  1. the definition of $\Psi^{p,q}_\Gamma(\tau,w)$ (series or construction);
  2. the weight convention (their $k+2$ versus our $k$);
  3. the residue formula and which of $w,\bar w$ pairs with $p$;
  4. the behaviour at the cusp (constant term?);
  5. the exact splitting statement and what makes it canonical (orthogonality? Hodge filtration?);
  6. Remark 7.6 and anything relating $\Psi$ at CM points to sums $\sum Q(\tau)^{-k}$.

  From Bengoechea (Math. Res. Lett. 22 (2015) 337–352; arXiv:1304.5653), extract the decomposition theorem and an explicit formula for the cusp-form part $c_D$ at weight 12, if given.
- **T20 (done here).** Run `residue_injectivity.py` to confirm the injectivity statement.
- **T21 (compute).** In `psi_projection_scaffold.py`:
  1. Implement `psi()`.
  2. Fix $c$ (and $p,q$) by matching residue polynomials exactly.
  3. Evaluate $f_{6,-7}-c\Psi$ at the sample points and fit to $\alpha E_{12}+\beta\Delta$. The fitter is self-tested: a synthetic form gives residual $3\times10^{-16}$, while $f_{6,-7}$ alone gives residual $0.77$, as it should.

  A small residual confirms the projection picture. Then compare $\beta$ with $c_D$ from T19.
- **T22 (cross-check).** If computable, compare $\beta$ with the regularized Petersson ratio $\langle f_{6,-7},\Delta\rangle^{\rm reg}/\langle\Delta,\Delta\rangle$ (Bringmann–Kane–von Pippich). This ties in with §8/T16.
- **T23 (extend).**
  - Repeat at $D=-8$ ($z_0=\sqrt{-2}$, form $[1,0,2]$) and $D=-11$.
  - Repeat at a non-CM pole, using the elliptic Poincaré series or the draft's $f_z$ at $k=18$ (simple pole: $p=0$ or $q=0$ per BF's convention).
  - Question: is $\beta$ algebraic×period-like only at CM points?
- **T24 (optional, Hecke).** Using the Hecke relations for $f_{k,D}$ (Zagier / Kohnen–Zagier for $D>0$; check Bengoechea for $D<0$), test whether $f\mapsto f^\circ$ commutes with $T_2$.
- **Stopping rule.** If BF's $\Psi$ are not explicitly computable, or their canonical splitting is only defined up to something that reintroduces the ambiguity, stop and report what the splitting is canonical *for*.

### `residue_injectivity.py`

```python
"""
T20: does the residue polynomial determine the principal part?
For f with principal part sum_{j=1}^{m} a_j (tau-z)^{-j} at z, the residue polynomial is
  R = sum_j a_j * (1/(j-1)!) d^{j-1}/dtau^{j-1} (X - tau Y)^n |_{tau=z}.
Claim: the map (a_1..a_m) -> R is injective iff m <= n+1.  So a form with zero residue
polynomials and poles of order <= n+1 = k-1 has no poles at all.
"""
import sympy as sp
X, Y, tau, z = sp.symbols('X Y tau z')
for n in (10, 16):
    for m in (n, n+1, n+2):
        vecs = []
        for j in range(1, m+1):
            d = sp.diff((X - tau*Y)**n, tau, j-1).subs(tau, z)/sp.factorial(j-1)
            P = sp.Poly(sp.expand(d), X, Y)
            vecs.append([P.coeff_monomial(X**(n-i)*Y**i) for i in range(n+1)])
        r = sp.Matrix(vecs).rank()
        print(f'n={n:2d}, pole order m={m:2d}: rank {r} of {m}  ->  injective: {r == m}')
```

### `psi_projection_scaffold.py`

```python
"""
T21 scaffold: canonical 'holomorphic projection' of Bengoechea's f_{6,-7} via Brown-Fonseca Psi.
  f°(tau) := f_{6,-7}(tau) - c * Psi^{p,q}(tau, z0)      (c, p, q fixed by matching residue polynomials)
If residue polynomials match and pole order <= k-1 (see residue_injectivity.py), f° is holomorphic in H,
so (given suitable cusp behaviour) f° = alpha*E_12 + beta*Delta.  Fit alpha, beta and check the residual.
TODO (Claude Code): implement psi() from Brown-Fonseca arXiv:2508.04844 (definition of Psi^{p,q}_Gamma,
their eq. (3.2)?; residue formula Cor. 3.11?) -- VERIFY equation numbers, conventions (which of w, w-bar
goes with p vs q), weight convention (their weight k+2 <-> our V_k), and normalisation.
"""
import numpy as np

def forms_disc(D, amax=400, bmin=120):
    out = []
    for a in range(1, amax+1):
        for b in range(-max(2*a, bmin), max(2*a, bmin)+1):
            if (b*b - D) % (4*a) == 0: out.append((a, b, (b*b - D)//(4*a)))
    return np.array(out, dtype=float)

FORMS = forms_disc(-7)
def f_beng(tau, k=6, F=FORMS):
    """Bengoechea normalisation: pi^{-k} sum_{b^2-4ac=D, a>0} (a tau^2 + b tau + c)^{-k}; weight 2k."""
    return np.sum((F[:, 0]*tau**2 + F[:, 1]*tau + F[:, 2])**(-float(k)))/np.pi**k

def sigma(m, r): return sum(d**r for d in range(1, m+1) if m % d == 0)
def E12(tau, N=40):
    q = np.exp(2j*np.pi*tau); return 1 + 65520/691*sum(sigma(m, 11)*q**m for m in range(1, N))
def Delta(tau, N=40):
    q = np.exp(2j*np.pi*tau); p = 1.0+0j
    for m in range(1, N): p *= (1 - q**m)**24
    return q*p

def psi(tau, w, p, q):
    raise NotImplementedError('implement Brown-Fonseca Psi^{p,q}(tau, w) here')

def fit_M12(g, pts):
    A = np.array([[E12(t), Delta(t)] for t in pts]); b = np.array([g(t) for t in pts])
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    return sol, np.max(np.abs(A @ sol - b))/np.max(np.abs(b))

# sample points away from the pole orbit of z0 = (-1+sqrt(-7))/2
PTS = [0.1+1.1j, -0.3+1.6j, 0.37+0.95j, 0.05+2.2j, -0.45+1.05j, 0.2+1.4j]

if __name__ == '__main__':
    # self-test of the fitter on a synthetic holomorphic form
    sol, res = fit_M12(lambda t: 0.3*E12(t) + 0.7*Delta(t), PTS)
    print('fitter self-test: alpha, beta =', np.round(sol, 12), ' residual', f'{res:.1e}')
    # sanity: f_beng itself is NOT in M_12 (it has poles) -> large residual expected
    sol, res = fit_M12(f_beng, PTS)
    print('f_{6,-7} alone: residual', f'{res:.2e}', '(should be O(1))')
    z0 = (-1 + np.sqrt(7)*1j)/2
    try:
        c, p_, q_ = None, 5, 5      # R(f_{6,-7}) = -(36 sqrt(-7)/7^5) * Q0^5  with Q0 = (X - z0 Y)(X - z0bar Y)
        raise NotImplementedError
    except NotImplementedError:
        print('TODO: fix c from the residue polynomial of Psi^{5,5}(., z0) and fit f_beng - c*Psi.')
```
