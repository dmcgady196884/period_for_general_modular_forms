# Research Summary: Periods, Cocycles, and L-functions for Modular Forms
## A working note on the mathematical program and its context

---

## 1. What the current note (`period_polynomials_dim_Sk_one`) actually contains

### Proven (rigorous):

- **Lemma 2.2** (cocycle property): The map $C^f_\gamma = (2\pi i)^{11}\int_{\gamma^{-1}\tau_0}^{\tau_0} f(\tau)(X-\tau Y)^{10}d\tau$ is a genuine 1-cocycle of $\Gamma$. The proof is path-additivity plus modularity of the integrand — the weight-12 transformation of $f(\tau)d\tau$ cancels the weight-$(-10)$ slash on $(X-\tau Y)^{10}$, making the integrand $\Gamma$-equivariant. This is rigorous.

- **Theorem 1.1** (period extraction): The universal kernel formula $(2)$ extracts $\omega^\pm(f)$ from two finite line segments, for all $f\in S^!_{12}$ (weakly holomorphic, vanishing constant term). The Bernoulli-operator derivation of the kernels is rigorous; basepoint-independence follows from an explicit coboundary identity.

- **Gap identity (20)**: For holomorphic $\Delta$, splitting the Mellin integral at the DR-B endpoints and applying S-modularity gives an exact identity relating $\Lambda^*(\Delta,s)$ to the BFK incomplete-gamma sum split at $\tau_0$. Rigorous, modulo standard analytic continuation.

### Sketch-level (proof is incomplete):

- **Theorem 8.1** (complex-$s$ extension): The Hurwitz-zeta formula $(22)$ for $\Lambda^*(\Delta, s)$ is proved by a $\tau_0 \to i\infty$ limiting argument that requires: (a) $\tilde{k}_T(\tau,s) = \zeta(1-s,\tau+1) - e^{i\pi(s-1)}\zeta(s-11,\tau+1)$ grows at most polynomially as $\mathrm{Im}(\tau)\to\infty$. This is true (standard Hurwitz asymptotics give $|\zeta(\sigma, \tau)| = O(\mathrm{Im}(\tau)^{1/2-\mathrm{Re}(\sigma)})$ by Euler-Maclaurin), but the explicit bound needs to be written to make it a proof rather than a sketch.

### Numerical conjecture (not proven analytically):

- **$\Lambda^*_\mathrm{DR}(\hat\Delta, s) = \Lambda^*_\mathrm{BFK}(\hat\Delta, s)$ for all $s\in\mathbb{C}$**: verified numerically at integer $s\in\{1,\ldots,11\}$ and several complex test points. There is no current analytical proof extending the gap identity (20) to weakly holomorphic forms.

---

## 2. The DR-B contours as cycles: an expository account for a string theorist

### The setup

You want to integrate the differential form
$$\omega_f = f(\tau)\,(X - \tau Y)^{k-2}\,d\tau$$
over paths in the upper half-plane $\mathbb{H}$.

For a string theorist: think of $\mathbb{H}$ as the worldsheet, $f(\tau)$ as the "matter" insertion, and $(X-\tau Y)^{k-2}$ as a polynomial "polarisation tensor" in auxiliary variables $(X,Y)$. The form $\omega_f$ is a holomorphic 1-form on $\mathbb{H}$ — closed, because it has no $d\bar\tau$ component. It's essentially the same structure as integrating a holomorphic 3-form over a 3-cycle in a Calabi-Yau, just one dimension down.

### The modular group acts, and makes a mess

The group $\Gamma = \mathrm{SL}_2(\mathbb{Z})$ acts on $\mathbb{H}$ by $\gamma:\tau\mapsto\frac{a\tau+b}{c\tau+d}$, and it acts on $(X,Y)$ by $(X,Y)\mapsto(aX+bY,cX+dY)$. The key fact is that $\omega_f$ is **equivariant**: the pullback of $\omega_f$ under $\tau\mapsto\gamma\tau$ equals $\omega_f|_\gamma$ in the coefficient system. So $\omega_f$ descends to a well-defined object on the modular curve $Y(\Gamma) = \Gamma\backslash\mathbb{H}$ — but it's a section of a nontrivial **local system** (the Sym$^{k-2}$ bundle over the modular curve), not a plain differential form.

### Periods = integrals over relative cycles

In classical Eichler-Shimura, you integrate $\omega_f$ along the geodesic from $0$ to $i\infty$. These are the two cusps of $\mathrm{SL}_2(\mathbb{Z})$, i.e. the two boundary points of the fundamental domain on $\partial\mathbb{H}$. The result is the period polynomial $r_f(X)$, whose coefficients are the critical $L$-values $L(f, 1), L(f,2),\ldots, L(f, k-1)$.

Why this path? Because in homology, the path $[0\to i\infty]$ is a **relative 1-cycle** in $H_1(\overline{X(\Gamma)}, \text{cusps}; V_{k-2})$ — a path with both endpoints at cusps. The cusps are the "marked points" on the compactified modular curve $X(\Gamma)$. Integrating a section of the local system over this relative cycle gives a period, in exactly the same sense as integrating $\Omega^{3,0}$ over a 3-cycle in a CY.

### What goes wrong for weakly holomorphic forms

If $f$ has a pole at the cusp $i\infty$ (i.e. $f\in S^!_k$), the integral $\int_0^{i\infty}\omega_f$ **diverges** — the integrand blows up at the upper limit. The form $\omega_f$ is still holomorphic on $\mathbb{H}$ (the pole is at the boundary, not in the interior), but the relative cycle $[0\to i\infty]$ hits the singularity.

The fix: don't push the endpoint to the cusp. Keep a **finite basepoint** $\tau_0\in\mathbb{H}$, strictly in the interior, and integrate over finite segments.

### The DR-B contour pair: cycles for a group cocycle

Here's the key geometric point. A **1-cocycle** of $\Gamma$ with coefficients in $V_n$ is a function $C:\Gamma\to V_n$ satisfying
$$C_{gh} = C_g|_h + C_h \quad\text{for all }g,h\in\Gamma.$$

Such a cocycle is **completely determined** by its values on the two generators $S = \bigl(\begin{smallmatrix}0&-1\\1&0\end{smallmatrix}\bigr)$ and $T = \bigl(\begin{smallmatrix}1&1\\0&1\end{smallmatrix}\bigr)$, subject to the relations $S^2 = (ST)^3 = 1$.

The DR-B construction: fix $\tau_0\in\mathbb{H}$ and define
$$C^f_\gamma = (2\pi i)^{k-1}\int_{\gamma^{-1}\tau_0}^{\tau_0} \omega_f.$$

This is a 1-cocycle because paths concatenate:
$$\int_{(gh)^{-1}\tau_0}^{\tau_0} = \int_{h^{-1}g^{-1}\tau_0}^{h^{-1}\tau_0} + \int_{h^{-1}\tau_0}^{\tau_0},$$
and the first piece, after the substitution $\tau = h\sigma$, gives $C^f_g|_h$ by modularity of $\omega_f$.

So the two contour segments are:

| Generator | Path | Physical meaning |
|-----------|------|-----------------|
| $S$ | $[-1/\tau_0,\, \tau_0]$ | the S-segment: from $S^{-1}\tau_0 = -1/\tau_0$ to $\tau_0$ |
| $T$ | $[\tau_0-1,\, \tau_0]$ | the T-segment: from $T^{-1}\tau_0 = \tau_0-1$ to $\tau_0$ |

These are not a clever trick. They are the **unique canonical paths** that read off the cocycle on the generators of $\Gamma$. The two segments are the generators of $H_1(\mathbb{H}; \tau_0\text{ and }\gamma^{-1}\tau_0)$ for $\gamma = S, T$ — relative cycles with prescribed endpoints related by $\Gamma$.

*See the companion conceptual document for the topological refinement of this picture: the contours project to the two generators of $\pi_1$ of the modular curve, and the period integrals are holonomies of the modular-form local system — the direct mathematical analog of Wilson loops in gauge theory.*

### The "gauge freedom" in $\tau_0$

Changing $\tau_0\to\tau_0'$ changes $C^f_\gamma$ by
$$C^f_\gamma(\tau_0') - C^f_\gamma(\tau_0) = (2\pi i)^{k-1}\left(\int_{\gamma^{-1}\tau_0'}^{\gamma^{-1}\tau_0} + \int_{\tau_0}^{\tau_0'}\right)\omega_f = (F(\tau_0') - F(\tau_0))|_\gamma - (F(\tau_0') - F(\tau_0)),$$
where $F(\tau) = (2\pi i)^{k-1}\int^{\tau}\omega_f$ (any primitive). This is a **coboundary**: the difference $C^f(\tau_0') - C^f(\tau_0) = \delta P$ for $P = F(\tau_0') - F(\tau_0)\in V_n$.

So $\tau_0$ is literally a gauge parameter. The cohomology class $[C^f]\in H^1(\Gamma; V_n)$ is gauge-invariant (basepoint-independent). The periods $\omega^\pm(f)$, which are coordinates of $[C^f]$ in the cuspidal subspace, are the "physical observables."

### The Mellin contour as a degenerate limit

The classical Mellin integral $\int_0^{i\infty}f(\tau)\tau^j d\tau$ arises by taking $\tau_0\to i\infty$ in the S-segment integral, with the polynomial kernel $k^j_S(\tau)$ collapsing to the bare monomial $\tau^j$. Simultaneously the T-segment integral vanishes exponentially, because $f$ is cuspidal and the segment $[\tau_0-1,\tau_0]$ contracts to a point at the cusp.

So: **the Mellin contour is the limit of the DR-B S-segment when the basepoint is pushed to the cusp.** This is why it works for holomorphic cusp forms (the limit exists) and fails for weakly holomorphic ones (the limit diverges).

### What changes when $f$ has interior poles

For meromorphic $f$ with a pole at $\tau_p\in\mathbb{H}$:

- $\omega_f$ is a **meromorphic** 1-form on $\mathbb{H}$, with a pole at $\tau_p$
- It is **no longer closed** everywhere; it has a residue $\mathrm{Res}(\omega_f, \tau_p) = \mathrm{Res}(f,\tau_p)\cdot(X-\tau_p Y)^{k-2}$
- The integral over a path crossing $\tau_p$ picks up $2\pi i$ times this residue

Concretely: moving $\tau_0$ across a pole $\tau_p$ changes $C^f_S$ by
$$\delta C^f_S = 2\pi i(2\pi i)^{k-1}\,\mathrm{Res}(f,\tau_p)\cdot(X-\tau_p Y)^{k-2},$$
a "wall-crossing" term. This breaks $\tau_0$-independence — the "gauge invariance" is anomalous, with the anomaly controlled by the residues of $f$.

In physics language: the cohomological "gauge invariance" has a **monodromy defect** at each interior pole. The correct framework is group cohomology of $\Gamma$ acting on $\mathbb{H}\setminus\Gamma\cdot\{\tau_p\}$, or equivalently relative cohomology of the punctured modular curve.

The residue formula (24) admits a clean OPE interpretation worth flagging: $\mathrm{Res}(f,\tau_p)$ plays the role of an OPE coefficient, and $(X-\tau_p Y)^{k-2}$ is its tensor structure on the auxiliary $(X,Y)$ polarization variables. The wall-crossing identity *is* an operator product expansion: it measures how the contour integral changes when the path brushes past the singular insertion. This is the technically-meaningful content of the §9 OPE analogy in the draft, not just rhetorical decoration.

---

## 3. What a paper would need to contain

### Tier 1: Completing the current note (relatively tractable)

1. **Explicit Hurwitz growth bound**: prove $|\zeta(1-s, \tau+1)| = O(\mathrm{Im}(\tau)^{\mathrm{Re}(s)-1})$ as $\mathrm{Im}(\tau)\to\infty$ (via Euler-Maclaurin), making Theorem 8.1 fully rigorous for holomorphic $f$.

2. **Analytical BFK equivalence for $\hat\Delta$**: extend the gap identity (20) to weakly holomorphic forms. The key step is showing the BFK regulator at the cusps agrees with the DR-B finite-contour construction analytically, not just numerically at integer $s$.

   **Concrete proof strategy** (developed in chat — needs verification by execution):
   
   **Step 1.** Fourier-decompose $\hat\Delta = q^{-1} + \sum_{n\geq 2} a(n)q^n$ and use linearity of the contour integrals in the $q$-modes. This is mode-by-mode linearity inside an integral over a fixed contour — not a decomposition of the modular form into non-modular pieces.
   
   **Step 2.** For each positive Fourier mode $n\geq 1$: apply the S-fold $\tau\mapsto -1/\tau$ to convert $\int_{-1/\tau_0}^{\tau_0}$ into $\int_i^{\tau_0}$ with the symmetrized kernel $[\tau^{s-1} - e^{i\pi(s-1)}\tau^{11-s}]$ (this is the BFK-style symmetric kernel; the second term is the S-image of the first). The T-segment dies exponentially against the cuspidal modes. The $\tau_0\to i\infty$ limit gives the BFK incomplete-gamma representation directly. This direction is straightforward.
   
   **Step 3.** For the polar mode $n=-1$: both S-segment and T-segment blow up like $e^{2\pi T} T^{s-1}$ as $\tau_0 = iT \to i\infty$, but $\tau_0$-independence of the full DR-B sum (which holds because $\hat\Delta(\tau)(X-\tau Y)^{10}d\tau$ is holomorphic on $\mathbb{H}$ — no interior poles) forces the divergent pieces to cancel. The Hurwitz-zeta T-kernel asymptotic
   $$\tilde k_T(\tau, s) \sim \frac{\tau^{s-1}}{s-1} - e^{i\pi(s-1)}\frac{\tau^{11-s}}{11-s} \quad\text{as }\mathrm{Im}(\tau)\to\infty$$
   matches exactly the divergent pieces of the symmetrized S-kernel — this is what makes Theorem 8.1 work, and it's the same fact that makes BFK = DR-B. Subtract these asymptotic divergent parts and show the remainder matches BFK's regulated polar value.
   
   This last step is the entire content of the proof; the rest is bookkeeping. It is tedious but appears tractable. **If it closes, this alone is a publishable paper**: "the L-value of a weakly holomorphic form is a finite-cycle holonomy, not a regulated divergent integral that happens to converge."

3. **Functional equation for $\Lambda^*_\mathrm{DR}(\hat\Delta, s)$**: the classical $\Lambda^*(\Delta, s) = i^k\Lambda^*(\Delta, k-s)$ should lift to a statement about the Hurwitz-zeta kernel using $\zeta(s,a) \leftrightarrow \zeta(1-s, \cdot)$ under the functional equation of Hurwitz zeta. This is likely tractable and currently unwritten.

### Tier 2: The genuinely new content

4. **Interior pole case, systematic theory**: 
   - The wall-crossing formula (24) is in the note; a theorem about what $\Lambda^*_\mathrm{DR}(f, s)$ actually computes for meromorphic $f$ with interior poles is not.
   - Analogy with your 1806 paper (regular at cusps, poles interior): what is the "right" L-function for such $f$? How does it relate to residues at poles?
   - Closed-form expressions analogous to your 1806 theorem, using Hurwitz zeta for the T-kernel and a modified S-kernel incorporating the pole residues.

5. **Functional equation for the meromorphic case**: likely involves monodromy data at the poles in an essential way. This is where things could get large.

### Tier 3: The bigger picture (long-term hobby territory)

6. **L-functions on non-compact modular curves**: the modular curve $Y(\Gamma) = \Gamma\backslash\mathbb{H}$ is non-compact; adding interior poles to $f$ is equivalent to working with sections of a bundle on a modular curve with additional punctures. The wall-crossing data is cohomological data for the punctured curve, which connects to:
   - Relative cohomology of $X(\Gamma)$ minus the $\Gamma$-orbit of the poles
   - Possibly Aomoto-Gelfand hypergeometric integrals, twisted de Rham theory
   - Connection to work of Hain-Brown on iterated Eichler integrals and $\mathcal{M}_{1,1}$

---

## 4. Literature connections

### Direct ancestors
- **Eichler (1957), Shimura (1959)**: the classical isomorphism $S_k\oplus\overline{S_k}\cong H^1_\mathrm{par}(\Gamma; V_{k-2})$.
- **Manin (1972)**: period polynomials and modular symbols; rationality of periods.
- **Zagier (1991)**: Periods of modular forms and Jacobi theta functions — the generating function perspective. The basis $P_0, P_1, P_2$ for $W_{12}$ is here.
- **Diamantis-Rolen**: the finite-endpoint cocycle for weakly holomorphic forms. The DR-B construction is the explicit unpacking of their (abstract) framework.
- **Brown (arXiv)**: quasi-periods $\eta^\pm(\hat\Delta)$ via the DR cocycle; the values the note verifies against.
- **BFK (Bringmann-Fricke-Kent, or similar)**: the incomplete-gamma regularised L-functions for weakly holomorphic forms.
- **McGady-arXiv:1806**: your own paper — L-functions for meromorphic forms regular at cusps, via Mellin-type contours. The interior-pole extension of the current note is the natural sequel.

### Structural context
- **BRST/BV formalism** (Becchi-Rouet-Stora-Tyutin, 1970s; Henneaux-Teitelboim): the general principle that cohomology classes = gauge-invariant observables, coboundaries = gauge transformations. The analogy in the note is not novel at this level of generality.
- **Lewis-Zagier period functions for Maass forms**: the note itself points to this as the Maass-form analog of Theorem 8.1. The Hurwitz-zeta kernel is the weight-$s$ analog of the Bernoulli-polynomial kernel, in the same way Lewis-Zagier period functions extend period polynomials to non-holomorphic settings.
- **Brown-Hain, Keilthy-Raum**: iterated Eichler integrals on $\mathcal{M}_{1,1}$; connections to Feynman amplitudes at genus 1. The "manifestly cohomological" construction of periods is part of this program.
- **Amplituhedron/on-shell methods** (Arkani-Hamed et al.): the analogy in the note's §9 — cocycle = off-shell data, periods = on-shell observables — is a rhetorical framing, not a technical connection. But the open question (construct $\omega^\pm(f)$ directly from $\{a_f(n)\}$ without the cocycle intermediate) is a real mathematical question.

---

## 5. The one-line version

**What you have**: A finite-contour formula for periods of holomorphic and weakly holomorphic weight-12 forms, with a numerical extension to all complex $s$ via Hurwitz zeta, and a clear picture of what breaks and how for meromorphic forms with interior poles.

**What the minimal small paper needs**: Just close the BFK = DR-B equivalence for $\hat\Delta$ via the Step 1–2–3 strategy above (chiefly: the Step 3 polar-term asymptotic cancellation). That alone is publishable as "regulator-free L-values for weakly holomorphic forms via finite-cycle holonomies." Add Tier 1 items 1 and 3 for completeness.

**What a richer paper needs**: The above plus a clean treatment of the interior-pole case (Tier 2) with closed-form L-function expressions, analogous to your 1806 paper but extended to forms with both cusp and interior singularities.

**What a larger project is**: An L-function theory for meromorphic modular forms on non-compact modular curves, where the "periods" are controlled partly by standard Eichler-Shimura data and partly by residue data at interior poles — a theory that doesn't currently exist in clean form.
