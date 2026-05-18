# RETRACTED 2026-05-18

This document originally argued that BFK's regularized period polynomial
for $\widehat\Delta$ is **not** in the cocycle space $W_{12}$ — that it
satisfies $(1+S)=0$ but fails $(1+U+U^2)=0$ — and proposed a
two-parameter holomorphic regulator $e^{\alpha\tau + \beta/\tau}$ to
"fix" this.

**The motivating numerical failure was a bug, not a real failure.**
The polynomial in question was produced by
`BFK_numerics_and_flaws/minimal_delta_comparison.sage`, which silently
dropped the $n=-1$ Fourier mode of $\widehat\Delta$.  With the $q^{-1}$
mode properly included, $r_{\mathrm{BFK}}(\widehat\Delta)$ satisfies
both period relations to machine precision (~$10^{-48}$ relative); see
`dr_b_Lfunction_complex_s.py` §(F–H).

Consequences:

1. BFK's regulator does **not** have an $S$-symmetry-without-$U$-symmetry
   problem at weight 12 with a $q^{-1}$ principal part.  The
   single-parameter formula in BFK Theorem 2.4 already lands in $W_{12}$.
2. The "two-parameter regulator $e^{\alpha\tau + \beta/\tau}$" proposal
   was a solution to a problem that does not exist.  It may or may not
   be useful for genuinely meromorphic modular forms (with interior
   poles), but it is not needed to land BFK in $W_{12}$ for cusp-pole
   weakly holomorphic forms.
3. The bgko_bfk_W_membership_gap.md companion is also retracted.

The canonical comparison between BFK and DR-B for $\Delta$ and
$\widehat\Delta$ is now in `dr_b_Lfunction_complex_s.py` §(F–H):
the two methods agree at all nine interior period-polynomial
coefficients (Brown's $\eta^\pm$-encoding monomials); they differ
only at the boundary $X^0$ and $X^{10}$ by an Eisenstein-direction
correction, which reflects the $\tau_0$-dependent boundary ambiguity
of the polynomial-$V_{10}$ DR-B representative (pinned down by the
Hurwitz-zeta $T$-kernel of the complex-$s$ theorem in
`period_polynomials_dim_Sk_one.tex`).
