# RETRACTED 2026-05-18

This document originally argued that BGKO Prop 2.3 and BFK p.6's claim
"$r(F;z) \in W$ for $F \in S_k^!$" is wrong, on the basis that
$\varepsilon(S)(ti) \sim t^{10} e^{2\pi t}$ grows exponentially and so
the polynomial-only projection can't be in $W$.  The §6 verification
protocol — "$r_{\mathrm{BFK}}(\widehat\Delta)$ passes $(1+S)=0$ exactly,
fails $(1+U+U^2)=0$ by ~22× typical coefficient magnitude" — was the
numerical anchor for that conclusion.

**That verification was wrong.**  It was performed on the polynomial
produced by `BFK_numerics_and_flaws/minimal_delta_comparison.sage`,
which silently dropped the $n=-1$ Fourier mode of $\widehat\Delta$
(bare `except: pass` swallowing the negative-argument incomplete-Gamma
coercion to RealField).

With the $q^{-1}$ mode properly included, $r_{\mathrm{BFK}}(\widehat\Delta)$
satisfies **both** $(1+S)=0$ and $(1+U+U^2)=0$ to machine precision
(~$10^{-48}$ relative).  See `dr_b_Lfunction_complex_s.py` §(F–H) for
the canonical computation.

The asymptotic statement $\varepsilon(S)(ti) \sim t^{10} e^{2\pi t}$ as a
*function* on $\mathbb{H}$ may still be defensible, but the conclusion
"therefore $r_{\mathrm{BFK}} \notin W$" does not follow — the
polar-mode contribution to $L^*(j+1)$ in BFK's formula evidently
restores the $U$-symmetry of the polynomial projection.  Working out
*why* structurally is an open thread.

BGKO Prop 2.3, BFK p.6, and Brown's framework are all consistent with
the corrected numerics; the "gap" this document claimed to identify
does not exist.
