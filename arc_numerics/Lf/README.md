# `Lf/` — the explicit formula for $L^*(f,s)$ in the geodesic class

**Empty on purpose.** There is no closed form for $L^*$ at $\tau_0=\rho+1$ yet, so
every $L^*$ value used anywhere in `arc_numerics/` is raw quadrature (see
`../common.py`). This directory is where the `thm:mero` analogue and its checks go
once it exists.

The authoritative inventory is the `%TODO(3)` block in
`finite_contour_cocycles_short.tex`, which mirrors the numbered results of §4.1:

- **(i) DONE** — `def:geoproj`, the projector thresholded at ${\rm Im}\,\tau_p\ge\sqrt3/2$
  (that is ${\rm Im}\,\tau_0$, because the suppression is $e^{-2\pi n\,{\rm Im}\,\tau_0}$).
- **(ii) open** — the `lem:proj` / `def:cRf` / `lem:projspace` analogue: the split
  $f=\widehat P f+Rf$ and the growth bound.
- **(iii) DONE** — `def:geoIN` and `lem:geopolylog`: `eq:geoG`, `eq:geoG0`, `eq:geoIN`,
  with the sheet identity `eq:sheet` $\Gamma^{\rm cont}(s,x)=e^{2\pi is}\Gamma(s,x)+(1-e^{2\pi is})\Gamma(s)$
  and the alternative form `eq:geoSupper`. Note the block integrals need **not** be done
  on the arc: both kernels are holomorphic on all of $\HH$, so evaluate on the
  horizontal chord $\rho\to\rho+1$ instead; crossing a block pole costs exactly
  $(-1)^N(2\pi i)^{-N}g^{(N)}(z)$.
- **(iv) open** — the `lem:polylogcf` analogue on the chord. Easier than at $\tau_0=i$:
  the chord is horizontal at $\sqrt3/2$, so the inversion criterion is a single global
  height comparison rather than a varying one.
- **(v) open** — `lem:polylogedge` / `def:FP`, the on-contour case at height $\sqrt3/2$.
- **(viii) open** — the `thm:mero` analogue itself.

Two coverage gaps worth flagging, because they are not obvious from the item list:
`lem:geopolylog` assumes ${\rm Im}\,z>\sqrt3/2$ **and $|z|\ge1$**, so it says nothing
about poles in the lens below the arc but above $\sqrt3/2$; and the inequality is
strict, so it says nothing on the threshold line itself.

Practical note carried over from the $\tau_0=i$ work: the $\gamma^*$ series loses about
$2.7n$ digits to cancellation, so `eq:geoSupper` is the better form if large $n$ is
ever needed.
