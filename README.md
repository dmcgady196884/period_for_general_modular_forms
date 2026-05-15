# Period polynomials for weight-12 modular forms

Numerical and analytic period extraction for weight-12 modular forms on
$\mathrm{SL}_2(\mathbb{Z})$ via the **Diamantis–Rolen finite-endpoint
cocycle**, in the $\dim S_k = 1$ setting.

The main result is a single linear functional that gives both the
classical periods $\omega^\pm$ of the discriminant
$\Delta \in S_{12}$ and Brown's quasi-periods $\eta^\pm$ of the weakly
holomorphic cusp form $\widehat\Delta \in S_{12}^!$ — same formula, same
kernels, different modular form under the integral.

## Main result

For any $f \in S_{12}^!$ (weakly holomorphic weight 12 with vanishing
constant Fourier coefficient — including the holomorphic
$\Delta \in S_{12} \subset S_{12}^!$), any basepoint $\tau_0 \in \mathbb{H}$,
and any interior monomial $j \in \{1, 2, \ldots, 9\}$:

$$
  \omega^{\varepsilon_j}(f) \;=\; \frac{(2\pi i)^{11}}{c_S^j}\,
  \Bigl[\,\int_{-1/\tau_0}^{\tau_0}\!f(\tau)\,k_S^j(\tau)\,d\tau
     \;+\; \int_{\tau_0-1}^{\tau_0}\!f(\tau)\,k_T^j(\tau)\,d\tau\,\Bigr]
$$

where $\varepsilon_j \in \{+, -\}$ is the parity sign of $j$, $c_S^j \in \mathbb{Q}$
is Brown's basis coefficient at the monomial $X^{10-j} Y^j$, and
$(k_S^j, k_T^j)$ are universal $\mathbb{Q}$-polynomial kernels of degree
$\le 10$ in $\tau$ — independent of $f$ and of $\tau_0$, tabulated in
[`period_polynomials_dim_Sk_one.pdf`](period_polynomials_dim_Sk_one.pdf).

Specialisations:

- $f = \Delta \;\Rightarrow\; \omega^\pm(\Delta)$, the Eichler–Shimura periods.
- $f = \widehat\Delta = \Delta(J^2 + 24J - 393444) \;\Rightarrow\; \eta^\pm(\widehat\Delta)$, Brown's quasi-periods.

## Numerical headline

All four values match Brown ([arXiv:1710.07912](https://arxiv.org/abs/1710.07912)
§ 8.3) to relative error $< 10^{-31}$ at 45-digit mpmath precision:

|                            | extracted                                | Brown § 8.3                   |
| -------------------------- | ---------------------------------------- | ----------------------------- |
| $\omega^+(\Delta)$         | $-68\,916\,772.80959519475431\ldots$     | $-68\,916\,772.809595194754\ldots$ |
| $\omega^-(\Delta)/i$       | $-5\,585\,015.379310401866877\ldots$     | $-5\,585\,015.3793104018668\ldots$ |
| $\eta^+(\widehat\Delta)$   | $127\,202\,100\,647.1770947773\ldots$    | $127\,202\,100\,647.17709477\ldots$ |
| $\eta^-(\widehat\Delta)/i$ | $10\,276\,732\,343.6491327508\ldots$     | $10\,276\,732\,343.649132750\ldots$ |

## Documents

1. **[`periods.pdf`](periods.pdf)** — bare numerics. Implementation of the
   DR cocycle in mpmath, verification against Brown's reported values.
2. **[`analytic_periods_weight12.pdf`](analytic_periods_weight12.pdf)** —
   long-form companion. Pedagogical derivation of the linear-functional
   formula via the Bernoulli-operator inverse of the $T$-coboundary
   $\delta_T$, worked numerical polynomials, full cohomological framework,
   appendix walkthrough.
3. **[`period_polynomials_dim_Sk_one.pdf`](period_polynomials_dim_Sk_one.pdf)** —
   math-paper restatement. Theorem 1 (the linear functional) front-and-center,
   supporting lemmas and propositions, outlook on (i) extending to
   meromorphic modular forms with interior poles and (ii) a gauge-invariance
   analogy with the on-shell S-matrix programme.

Each is self-contained but they share results; read in order or jump
directly to (3) for the punchline.

## Computational toolchain

- **`period_extraction_DR_general.ipynb`** — the main pipeline. A single
  $q$-series-agnostic function `extract_periods(coeffs_dict, tau0, k)`
  that takes any weight-$k$ form's Laurent-series coefficients and returns
  the periods. Blind to whether $f$ is holomorphic or weakly holomorphic
  — same code, same precision, $\omega^\pm$ or $\eta^\pm$ as appropriate.
- **`delta_and_hatdelta_periods_executed.ipynb`** — the original
  verification notebook (the codebase the main pipeline grew out of).
- **`compute_period_kernels.py`** — exact sympy-rational construction
  of the 9 universal kernel pairs $(k_S^j, k_T^j)$.
- **`compute_numerical_polynomials.py`** — explicit numerical
  $(P|_S - P)$ and $C'^f_S$ polynomials at $\tau_0 = 0.3 + 1.2i$ for
  both $\Delta$ and $\widehat\Delta$.
- **`verify_bernoulli_formula.py`** — confirms the Bernoulli closed-form
  for $P$ agrees with the upper-triangular recursion at $10^{-42}$
  relative residual.
- **`verify_omega_plus_formula.py`** — confirms the 9 explicit period
  formulas reproduce $\omega^\pm(\Delta)$ via $L(\Delta, s)$ at $10^{-42}$
  relative residual.

## Handoff notes (BFK gap diagnosis)

- **`bfk_regularization_handoff.md`** — why BFK's regularization gives a
  polynomial outside $W_{12}$: single-parameter regulator has $S$-symmetry
  but no $U$-symmetry; $(1 + S) = 0$ survives, $(1 + U + U^2) = 0$ fails.
- **`bgko_bfk_W_membership_gap.md`** — structural proof that
  BGKO Prop. 2.3's "$r(F; z) \in W$" silently discards exponentially-growing
  pieces of the Eichler integral. Numerical evidence: BFK polynomial for
  $\widehat\Delta$ fails $(1 + U + U^2) = 0$ by $\sim 22\times$ the typical
  coefficient magnitude.
- **`BFK_numerics_and_flaws/`** — the underlying numerical experiments,
  Sage scripts, and writeup that established the BFK $\notin W$ diagnosis.

## Reproducing

```bash
git clone https://github.com/dmcgady196884/period_for_general_modular_forms.git
cd period_for_general_modular_forms
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Register the venv as a Jupyter kernel named "mmf_venv" (matches what
# the notebooks reference)
python -m ipykernel install --user --name mmf_venv --display-name mmf_venv

# Execute the main pipeline (~2-3 minutes at 45-digit precision)
jupyter nbconvert --to notebook --execute --inplace \
        --ExecutePreprocessor.kernel_name=mmf_venv \
        period_extraction_DR_general.ipynb

# Compile the math-paper document (requires a TeX Live distribution)
pdflatex period_polynomials_dim_Sk_one.tex
pdflatex period_polynomials_dim_Sk_one.tex   # second pass for cross-refs
```

Tested with Python 3.14, mpmath 1.3.0, sympy 1.14.0, TeX Live 2025 on macOS.

## References

- F. Brown, *A class of non-holomorphic modular forms III*,
  [arXiv:1710.07912](https://arxiv.org/abs/1710.07912) (2017) — definition
  and numerical values of $\eta^\pm(\widehat\Delta)$, the period-matrix
  convention, Brown's basis at $S$ for $H^1_{\rm cusp}(\Gamma; V_{10})$.
- N. Diamantis & L. Rolen, *Period polynomials, derivatives of $L$-functions,
  and zeros of polynomials*,
  [arXiv:1707.04814](https://arxiv.org/abs/1707.04814) (2017) — the
  finite-endpoint cocycle on which everything here rests.
- W. Kohnen & D. Zagier, *Modular forms with rational periods*,
  in *Modular Forms* (R. A. Rankin, ed.), Ellis Horwood, 1984 — Zagier's
  basis $\\{P_0, P_1, P_2\\}$ for $W_{12}$.
- K. Bringmann, K.-H. Fricke, Z. Kent, *Special $L$-values and periods
  of weakly holomorphic modular forms*, *Ramanujan J.* **24** (2011) —
  the BFK regularisation we diagnose.

## Status and successor problems

The weight-12 / $\dim S_k = 1$ case is fully worked out. Open threads
flagged for future work:

- **Higher weights with $\dim S_k > 1$** (weights 24, 28, 30, …) —
  cuspidal cohomology of higher dimension; multiple Hecke eigenforms;
  what does the linear functional look like there?
- **Fully meromorphic modular forms** with interior poles
  (e.g. $1/E_6$) — residue contributions weighted by $(X - \tau_p Y)^j$
  polynomials; relative-cohomology formulation; the natural successor
  to the present repo, foreshadowed in
  [`period_polynomials_dim_Sk_one.pdf`](period_polynomials_dim_Sk_one.pdf)
  §8 Outlook.
- **A "manifestly on-shell" formulation of periods** — bypassing the
  cocycle/coboundary scaffolding entirely. Modular-form analog of
  BCFW / on-shell unitarity; see same §8 Outlook for the gauge-invariance
  framing.
