# Periods and $L$-functions of meromorphic modular forms via finite contours

A single, manifestly finite contour-integral functional — the **$L$-integral**
$L^*(f, s; \gamma^S_{\tau_0}, \gamma^T_{\tau_0})$ — that unifies three previously
separate constructions for $\mathrm{SL}_2(\mathbb{Z})$ modular forms:

1. the periods $\omega^\pm$ and quasi-periods $\eta^\pm$ of (weakly holomorphic)
   cusp forms — Brown, [arXiv:1710.07912](https://arxiv.org/abs/1710.07912);
2. the $L$-function of weakly holomorphic forms $f \in M^!_k$ —
   Bringmann–Fricke–Kent (BFK);
3. the $L$-function of meromorphic forms with poles inside the fundamental
   domain — McGady, [arXiv:1806.09874](https://arxiv.org/abs/1806.09874).

All three fall out of the *same* functional — at integer $s$ it is the
period extractor, at general $s$ it is the completed $L$-function — finite at
every basepoint $\tau_0 \in \mathbb{H}$, with no regulator and no
analytic-continuation step.

**Active paper: [`dr_b_periods_and_Lfunctions.tex`](dr_b_periods_and_Lfunctions.tex)** (~31 pp).

## The functional

For $f \in F_k$ and a basepoint $\tau_0 \in \mathbb{H}$ whose two contour
segments avoid the $\Gamma$-orbit of $f$'s interior poles,

$$
  L^*(f, s; \gamma^S_{\tau_0}, \gamma^T_{\tau_0})
  \;=\; i^{-s}\Bigl[\,
    \int_{\gamma^S_{\tau_0}}\! f(\tau)\,\tau^{s-1}\,d\tau
    \;+\;
    \int_{\gamma^T_{\tau_0}}\! f(\tau)\,\widetilde{k}_T(\tau, s)\,d\tau
  \,\Bigr],
$$

where $\gamma^S_{\tau_0}$ runs from $-1/\tau_0$ to $\tau_0$,
$\gamma^T_{\tau_0}$ from $\tau_0 - 1$ to $\tau_0$, and the $T$-kernel is the
two-term Hurwitz zeta

$$
  \widetilde{k}_T(\tau, s)
  \;=\; \zeta(1 - s,\, \tau + 1)
  \;-\; e^{i\pi(s-1)}\,\zeta(s - (k-1),\, \tau + 1).
$$

Here $F_k$ is the space of **meromorphic modular forms** of weight $k$ — the
weight-$k$ part of the field of fractions of $M_* = \bigoplus_k M_k$ (quotients
$f_1/f_2$ with $f_i \in M_{k_i}$, $k = k_1 - k_2$) — so
$M_k \subset M^!_k \subset F_k$.

## Three theorems

- **Theorem 1.1 (periods, $\dim S_k = 1$).** For $f \in S^!_k$, the integer-$s = \ell$
  specialisation reduces to an explicit $\mathbb{Q}$-rational linear functional in
  the Fourier coefficients of $f$ that produces $\omega^\pm(f)$ / $\eta^\pm(f)$,
  independent of $\tau_0$ and of the interior index $\ell$. At $k = 12$ it
  reproduces Brown's $\omega^\pm(\Delta)$ and $\eta^\pm(\widehat\Delta)$.
- **Theorem 1.2 ($L$-function on $M^!_k$, all even $k$).** For $f \in M^!_k$,
  $L^*(f, s)$ is $\tau_0$-independent for every $s \in \mathbb{C}$; at the
  elliptic fixed point $\tau_0 = i$ the $S$-segment collapses and $L^*$ becomes
  a single integral on $[i-1, i]$ equal to the BFK incomplete-gamma sum — the
  completed $L$-function, recovering Hecke–Weil for holomorphic forms and BFK
  for weak forms.
- **Theorem 1.3 (meromorphic forms).** For general meromorphic $f$, $L^*$
  depends on the contour pair only through its homotopy class in
  $\mathbb{H} \setminus \Gamma\!\cdot\! f^{-1}(\infty)$, with explicit
  wall-crossing residue jumps across interior poles; on forms regular at the
  cusp it coincides with the deformed-Mellin $L$-function of arXiv:1806.09874.

The §5 material further gives a closed-form per-mode $L^*$ for forms with poles
**both** at interior points and at the cusp (via a Hauptmodul split
$f = g + h$, valid whenever $S_{2-k} = \{0\}$ — i.e. all $k \ge 0$ and
$k \in \{-2,-4,-6,-8,-12\}$).

## Numerical verification (50-digit `mpmath`)

| script | checks | residual |
| --- | --- | --- |
| [`verify_periods.py`](verify_periods.py) | Thm 1.1 at $k=12$: Brown's $\omega^\pm(\Delta)$, $\eta^\pm(\widehat\Delta)$; basepoint- and $\ell$-independence; Bernoulli $\delta_T^{-1}$; the kernel identity (sympy) | $\lesssim 10^{-40}$ |
| [`dr_b_Lfunction_complex_s.py`](dr_b_Lfunction_complex_s.py) | Thm 1.2 at non-integer $s$ (half-integer, off-axis, outside-strip) | $\sim 10^{-53}$ ($\Delta$), $\sim 10^{-47}$ ($\widehat\Delta$) |
| [`verify_s6.py`](verify_s6.py), [`verify_multi_s.py`](verify_multi_s.py) | Thm 1.2 per-mode at $\tau_0 = i$ = BFK, $s \in \{1,3,6,9,11\}$ incl. boundaries and the polar mode $n=-1$ | $\lesssim 10^{-40}$ |

Brown's reference values (1710.07912 §8.3):
$\omega^+(\Delta) = -68\,916\,772.8096\ldots$,
$\omega^-(\Delta)/i = -5\,585\,015.3793\ldots$,
$\eta^+(\widehat\Delta) = 1.27202\times10^{11}$,
$\eta^-(\widehat\Delta)/i = 1.02767\times10^{10}$.

## File map

- **[`dr_b_periods_and_Lfunctions.tex`](dr_b_periods_and_Lfunctions.tex)** — the paper (three theorems + full proofs).
- **`wall_crossing_and_poincare_sums.tex`** — companion working notes: cumulative
  wall-crossing as a weight-$(2-k)$ Poincaré-style orbit sum.
- Verification scripts: the four listed above, plus older helpers
  `compute_period_kernels.py`, `compute_numerical_polynomials.py`,
  `verify_bernoulli_formula.py`, `verify_omega_plus_formula.py`,
  `dr_vs_bfk_modewise_nogo.py`.
- Superseded / historical (kept for the record): `period_polynomials_dim_Sk_one.tex`
  (the older $\dim S_k=1$ note), `analytic_periods_weight12.tex`, `periods.tex`,
  `tau0_eq_i_proof.tex` (now folded into Theorem 1.2).
- `me_1806.pdf` — McGady, arXiv:1806.09874.

## Reproducing

```bash
git clone https://github.com/dmcgady196884/period_for_general_modular_forms.git
cd period_for_general_modular_forms
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# numerics (~30 s each)
python verify_periods.py
python dr_b_Lfunction_complex_s.py

# compile the paper (TeX Live; second pass for cross-refs)
pdflatex dr_b_periods_and_Lfunctions.tex
pdflatex dr_b_periods_and_Lfunctions.tex
```

Tested with mpmath 1.3.0, sympy 1.14.0, TeX Live 2025.

## References

- F. Brown, *A class of non-holomorphic modular forms III*,
  [arXiv:1710.07912](https://arxiv.org/abs/1710.07912) (2017).
- K. Bringmann, K.-H. Fricke, Z. Kent, *Special $L$-values and periods of weakly
  holomorphic modular forms*, *Ramanujan J.* **24** (2011).
- D. A. McGady, *$L$-functions for meromorphic modular forms and sum rules in
  CFT*, [arXiv:1806.09874](https://arxiv.org/abs/1806.09874) (2019).
- N. Diamantis & L. Rolen, *Period polynomials, derivatives of $L$-functions,
  and zeros of polynomials*, [arXiv:1707.04814](https://arxiv.org/abs/1707.04814) (2017).
- W. Kohnen & D. Zagier, *Modular forms with rational periods*, in *Modular Forms*
  (R. A. Rankin, ed.), Ellis Horwood, 1984.

## Status

The $\dim S_k = 1$ periods (Thm 1.1) and the $M^!_k$ / meromorphic $L$-functions
(Thms 1.2–1.3) are worked out with full proofs and numerics. Open threads:
higher-dimensional $S_k$ (matrix-valued kernels); a systematic per-mode closed
form on interior-pole forms parallel to the BFK sum; and the geometric content
of the bare Mellin kernel $\tau^{s-1}\,f\,d\tau$ (conjecturally a section of a
Lewis–Zagier-type local system on $Y(1)$).
