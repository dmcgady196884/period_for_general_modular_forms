# Handoff: Holomorphic Regularization, BFK, and Period-Polynomial Failure

## TL;DR

BFK's regularized period polynomial for weakly holomorphic forms with cusp poles is **not in the cocycle space $W_k$** — despite their stated claim. It satisfies $(1+S) = 0$ but fails $(1+U+U^2) = 0$. The structural reason is that BFK's single-parameter regulator $e^{iu\tau}$ has $S$-symmetry (under $\tau \to -1/\tau$ it becomes $e^{-iu/\tau}$, regulating the other cusp) but no analogous $U$-symmetry. A two-parameter holomorphic regulator $e^{\alpha\tau + \beta/\tau}$ with **independent** analytic continuations $\alpha \to 0$, $\beta \to 0$ is proposed as a unifying scheme that:

1. Reproduces BFK's incomplete-Gamma formulas in their domain (cusp poles only)
2. Reproduces McGady 1806's incomplete-Gamma formulas in his domain (interior poles, regular at cusps)
3. Potentially recovers the Brown / Diamantis-Rolen quasi-periods $\eta^\pm$ that BFK's $\alpha = \beta$ constraint biases away from

## Setup and notation

- $\Delta(\tau) = q\prod(1-q^n)^{24} \in S_{12}$, the unique weight-12 cusp form
- $\hat\Delta = \Delta'$ (Brown's notation): the unique weakly holomorphic form in $S_{12}^!$ with principal part $q^{-1}$, satisfying $\hat\Delta = q^{-1} + 0 + 0 + 47709536 q^2 + \cdots$. Construction: $\hat\Delta = \Delta \cdot (j^2 - 1464 j + 142236)$ or equivalently $\Delta(J^2 + 24 J - 393444)$ with $J = j - 744$.
- Brown's periods $\omega^\pm$ of $\Delta$ and quasi-periods $\eta^\pm$ of $\hat\Delta$, with $|\omega^+| \approx 6.89 \times 10^7$, $|\eta^+| \approx 1.27 \times 10^{11}$, satisfying $\det\begin{pmatrix}\eta^+ & \omega^+ \\ i\eta^- & i\omega^-\end{pmatrix} = 10!\,(2\pi i)^{11}$.
- $W_{12} = \{P \in V_{10} : P|(1+S) = P|(1+U+U^2) = 0\}$, with $S = \begin{pmatrix}0&-1\\1&0\end{pmatrix}$, $U = \begin{pmatrix}1&-1\\1&0\end{pmatrix}$. Slash action of weight $-10$.
- $\dim W_{12} = 3$: $W_{12}^+ = \langle P_0, P_1\rangle$ (dim 2), $W_{12}^- = \langle P_2\rangle$ (dim 1), where:
  - $P_0(X) = X^{10} - 1$ (Eisenstein direction)
  - $P_1(X) = X^8 - 3X^6 + 3X^4 - X^2$ (cuspidal even)
  - $P_2(X) = 4X^9 - 25X^7 + 42X^5 - 25X^3 + 4X$ (cuspidal odd)

## Three regularization schemes

### BFK (Bringmann-Fricke-Kent, 2011)

Non-holomorphic regulator $e^{2\pi \nu \mathrm{Im}(\tau)}$, equivalently $e^{iu\tau}$ on the imaginary axis, taken in the limit $\nu \to 0$ ($u \to 0$). Defines for $f \in S_k^!$:

$$L^*_f(s) = \sum_{m \geq m_0} a_f(m)\frac{\Gamma(s, 2\pi m t_0)}{(2\pi m)^s} + i^k \sum_{m \geq m_0} a_f(m) \frac{\Gamma(k-s, 2\pi m / t_0)}{(2\pi m)^{k-s}}$$

Key BFK results:
- Theorem 2.2: $L^*_f(s) = R.\int_0^\infty f(iy) y^{s-1} dy$, functional equation $L^*_f(k-s) = i^k L^*_f(s)$ ✓
- Theorem 2.4 *identity*: the polynomial $r(f;z) := c_k(E_f - E_f|_{2-k}S)$ equals $\sum i^{1-n} \binom{k-2}{n} L^*_f(n+1) z^{k-2-n}$ ✓
- Claim (asserted, not proved in paper): $r(f;z) \in W$ — cited from BGKO and Kohnen-Zagier. **This is the claim that fails.**

### McGady 1806.09874 (2018)

Holomorphic contour deformation: integrate $\int_{i/\Lambda}^{i\Lambda} f(\tau) \tau^{s-1} d\tau$, deforming around interior poles, take $\Lambda \to \infty$. Works for $f \in F_k$ with poles in $\mathcal{H}$, regular at cusps. Functional equation $L^*_f(s) = i^k L^*_f(k-s)$ proved (Cor 2.15). Says nothing about cusp poles (Lemma 2.2: divergent for those).

### Two-parameter holomorphic regulator (proposed)

$$e^{\alpha \tau + \beta/\tau} \cdot t^{-s-1}$$

inserted in the integrand, with the integral split at a base point $c \in \mathbb{H}$:

$$L^*_f(s) = I_1(\alpha, \beta) + I_2(\beta, -\alpha)$$

where $I_2$ is the S-image of the $0 \to c$ piece. Initial domain of convergence:
- $I_1$: $\mathrm{Re}(\alpha) > 2\pi N$ where $f$ has $q^{-N}$ leading term at $i\infty$
- $I_2$: $\mathrm{Re}(\beta) < -2\pi N$ (after S, this becomes $\mathrm{Re}(-\beta) > 2\pi N$)

**Key claim:** Each $I_j$ is separately analytic in $(\alpha, \beta)$ in its domain, and admits independent analytic continuation:
- $I_1$: continue $\alpha \to 0$ with $\beta$ fixed at some $\beta_0$ outside the convergence boundary
- $I_2$: continue $\beta \to 0$ with $\alpha$ fixed at some $\alpha_0$

Each piece reduces to **incomplete Gamma functions** in this asymmetric limit (cleaner than the symmetric $\alpha = -\beta$ case which gives K-Bessel functions). The asymmetric limit is "S-symmetric in the appropriate sense": both limits are "remove the cusp regulator at the appropriate cusp" in the $i\infty$-coordinate of each piece.

## The numerical experiment (`period_synthesis_bfk.pdf`, §8)

Implementation: BFK formula with corrected sum over $n \neq 0$ (including $n = -1$ for the principal part), naked $n$ values (not $|n|$), at high precision via mpmath. Computed $L^*_{\hat\Delta}(s)$ for $s = 1, \ldots, 11$, then assembled the period polynomial coefficients via Theorem 2.4.

### What the data shows

**Odd-part coefficient ratios** (from $X^9, X^7, X^5, X^3, X$):
- $\Delta$: $(1, -6.25, 10.5, -6.25, 1) = (4, -25, 42, -25, 4)/4$ — exact match to $P_2$
- $\hat\Delta$ (BFK): $(1, -11.029, 22.544, -11.029, 1)$ — symmetric but **not** $P_2$

**Even-part coefficient ratios** (from $X^{10}, X^8, X^6, X^4, X^2, X^0$):
- $\Delta$: $(1, -19.194, 57.583, -57.583, 19.194, -1)$ — combination of $P_0$ and $P_1$ with the $691$
- $\hat\Delta$ (BFK): $(1, -40.051, 176.931, -176.931, 40.051, -1)$ — antisymmetric but not in $\langle P_0, P_1\rangle$

Internal accuracy: $< 10^{-27}$ for $\Delta$; numerics for $\hat\Delta$ are stable, showing factors of 1.5–3.0 deviation from $P_j$ ratios — far outside numerical error.

### Diagnosis

$\dim W_{12}^- = 1$. Anything in $W_{12}^-$ is a scalar multiple of $P_2$. BFK's odd-part polynomial for $\hat\Delta$ has different ratios → **not in $W_{12}^-$**. Symmetric coefficient pattern → $(1+S) = 0$ does hold. By exclusion, $(1+U+U^2) = 0$ **fails**.

(This should be confirmed by directly computing $(1+U+U^2)$ on the polynomial. Mechanical: $P|U(X) = X^{10} P((X-1)/X)$, $P|U^2(X) = (X-1)^{10} P(-1/(X-1))$.)

### Why this happens

BFK's regulator $e^{iu\tau}$ has the property that under $\tau \to -1/\tau$, it becomes $e^{-iu/\tau}$, regulating the *other* cusp. So the analytic continuation in $u$ commutes with $S$, and $(1+S)$ relations survive. But under $\tau \to (\tau-1)/\tau$ (the $U$-action), the regulator does not transform into something that regulates a cusp. The three legs of the $U$-cycle ($\int_0^{i\infty} + \int_{i\infty}^1 + \int_1^0$) each require regularization, and a single parameter $u$ cannot regulate all three correctly. The analytic continuation in $u$ does not commute with $U$, breaking $(1+U+U^2) = 0$.

This is the structural obstruction. It cannot be fixed within a single-parameter holomorphic regularization.

### What works: Brown / Diamantis-Rolen

Finite-endpoint cocycle: $C^f_\gamma(X, Y) = \int_{\tau_0}^{\gamma\tau_0} f(\tau)(X - \tau Y)^{k-2} d\tau$ for $\tau_0 \in \mathcal{H}$. The contour stays in $\mathbb{H}$, never approaching cusps; modular symmetries are preserved exactly because the regularization mechanism is "don't go to cusps" rather than "regulate at cusps." Numerical reproduction of Brown's $\eta^\pm$ to ~28 digits is in `periods.pdf`.

## What the proposed two-parameter scheme could deliver

**Plausible:**
- Reproduces BFK's L-values where BFK applies (forms with cusp poles, regular at interior)
- Reproduces 1806's L-values where 1806 applies (interior poles, regular at cusps)
- For forms with both pole types, gives a unified formula that's just BFK + 1806 incomplete-Gamma terms glued

**Open / requires verification:**
- Whether the resulting period polynomial for $\hat\Delta$ lands in $W_{12}$
- Whether it matches Brown's quasi-period polynomial (i.e., reproduces $\eta^\pm$ from a regularized cusp-to-cusp integral)
- Whether $(1+U+U^2) = 0$ is recovered with two independent parameters where it failed with one

If the two-parameter scheme *does* recover $\eta^\pm$ from a regularized cusp-to-cusp integral, this is a stronger result than "regularization unification" — it would mean quasi-periods *are* reachable by direct integration when you use the right regulator, contrary to the naive view that they're inherently cocycle objects.

## Concrete next steps

1. **Verify $(1+U+U^2) = 0$ failure on BFK's polynomial.** Take the polynomial coefficients from `period_synthesis_bfk.pdf` Tables 2 and 3 ($\hat\Delta$ data), apply $(1+U+U^2)$ via the slash formulas, confirm nonzero. This pins down the diagnosis.

2. **Implement two-parameter regularized L-values for $\hat\Delta$.** The integrand on the imaginary axis is $\hat\Delta(iy) y^{s-1} e^{-\alpha y - \beta/y}$ with appropriate signs. Compute:
   - $I_1(\alpha, \beta_0) = \int_{c}^\infty \hat\Delta(iy) y^{s-1} e^{-\alpha y - \beta_0/y} dy$ for $\mathrm{Re}(\alpha) > 2\pi$, $\beta_0 < -2\pi$
   - Continue $\alpha \to 0$ holding $\beta_0$ fixed
   - Symmetrically for $I_2$
   - Sum, evaluate at $s = 1, \ldots, 11$, assemble polynomial coefficients

3. **Compare resulting polynomial to:**
   - BFK's polynomial (should differ by a polar correction term)
   - Brown's quasi-period polynomial $\eta^+(\alpha_0 P_0 + \alpha_1 P_1) + i\eta^- \alpha_2 P_2$ (this is the test of whether the new scheme lands in $W$ and matches cohomology)

4. **Test $(1+S)$ and $(1+U+U^2)$ on the new polynomial.** If both hold, the unification claim has teeth. If $(1+U+U^2)$ still fails, the regulator design needs revision.

## References (in `/mnt/project/`)

- `1806_09874.pdf` — McGady, "L-functions of meromorphic modular forms" (regular at cusps, interior poles)
- `BFK.pdf` — Bringmann-Fricke-Kent, "Special L-values and periods of weakly holomorphic modular forms"
- `1710_07912.pdf` — Brown, "A class of non-holomorphic modular forms III" (defines $\eta^\pm$, gives numerical values)
- `MF_Periods_and_Jacobi_theta_Zagier_91.pdf` — Zagier, classical periods and $W_k$
- `period_rationality_20251128.pdf` — earlier exposition of rationality patterns and BFK setup
- `period_synthesis_bfk.pdf` — the actual numerical experiment showing BFK polynomial fails template match
- `periods.pdf` — finite-endpoint cocycle reproduction of Brown's $\omega^\pm, \eta^\pm$ (reference implementation that *does* land in $W$)
- `francis_brown_quasiperiod_email.pdf` — Brown's email explaining the cocycle construction

## Convention note

In Brown's period matrix the basis is asymmetric: $P = \begin{pmatrix}\eta^+ & \omega^+ \\ i\eta^- & i\omega^-\end{pmatrix}$. Numerical extractions yield $\omega^-, \eta^-$ as $i$ times Brown's reported real values; magnitudes match exactly. See `periods.pdf` §6 convention note.
