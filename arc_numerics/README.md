# arc_numerics — the geodesic reference class, $\tau_0=\rho+1$

Numerics for §4.4–4.5 of `finite_contour_cocycles_short.tex`: the reference class in
which $\gamma^S=\gamma^T=\gamma^{\rm arc}$, the unit arc from $\rho$ to $\rho+1$.

Everything in `end_to_end_numerical_testing/` and the tracked top-level `*.py` belongs
to the **other** class, $\tau_0=i$ with the horocyclic triangle (§4.1–4.3). The two
have different contours, different projector thresholds and different game-boards; do
not mix results between them.

Run with `/Users/dmcgady/Documents/math/mmf_venv/bin/python` (mpmath 1.3).

## How $L^*$ is computed here

By **raw quadrature** of the two contour integrals of `def:Lint` — `mp.quad` over the
arc. There is no closed form for $L^*$ in this class: the `thm:mero` analogue is
unwritten (item (viii) of the `%TODO(3)` block in the .tex). The only special function
anywhere is the kernel `ktil` $=\zeta(1-s,\tau+1)-e^{i\pi(s-1)}\zeta(1-(k-s),\tau+1)$,
inside the integrand.

Consequence worth stating: **the $r_f\in W$ results are logically independent of the
unfinished $L^*$ machinery** and cannot inherit an error from it.

## Layout

- `common.py` — the only definition of the arc contour, plus the `def:Wpm` basis,
  the $V_n$ operators, and Neville extrapolation. **Import it; do not re-derive.**
- `rf/` — period polynomials, $\tilde r_f\in W$.
- `Lf/` — the explicit formula for $L^*$. Currently empty; see its README.

## Two traps, both of which cost real time

1. **Orientation.** $\gamma^{\rm arc}$ runs $\rho\to\rho+1$, i.e. $\theta$ *decreasing*
   from $2\pi/3$ to $\pi/3$. Reversing it flips $r_f$, $\Phi$ and $r_{E_k}$ together,
   turning $r_f-\Phi r_{E_k}$ into $-(r_f+\Phi r_{E_k})$ — a wrong combination that
   still looks plausible and surfaces only as a spurious $(1+U+U^2)$ defect of order 1.
   `common.check_orientation()` asserts $\Phi(E_{12})=+1$; every script calls it.
2. **Which $\Phi$.** `def:Phi` is over the **full arc**. The `eq:geodelta` chord
   integral is a different number — using it makes $\Phi(f)=0$ for a form with a pole
   at $\rho$ and nothing lands in $W$.

A third, from the earlier round: condition 2 of `eq:Fcirc` is
$\int_{\gamma^T}f\,d\tau=0$, which equals $c_f(0)=0$ **only** when no pole sits above
the contour. The arc dips to $\sqrt3/2$, so it fails generically. $\Delta$ is a useless
control for this — being holomorphic it satisfies condition 2 on every contour.

## Status board, by pole location

Relative to $\gamma^{\rm arc}$ (height $\sqrt3/2$ at the endpoints, $1$ at $i$):

| where the pole is | $r_f\in W$ | $L^*$ closed form |
|---|---|---|
| above the arc, $\lvert\tau\rvert>1$ | `thm:geoperiod`, no condition on location or order | `lem:geopolylog` (needs ${\rm Im}\,z>\sqrt3/2$ **and** $\lvert z\rvert\ge1$) |
| below the arc, ${\rm Im}>\sqrt3/2$ | `thm:geoperiod` | **open** — `lem:geopolylog` assumes $\lvert z\rvert\ge1$ |
| on ${\rm Im}\,\tau=\sqrt3/2$ | `thm:geoperiod` | **open** — item (v), and the inequality is strict |
| on the arc, $\tau\ne i,\rho$ | $S$-symmetric indentation ($\log r$ odd about $\theta=\pi/2$); legitimate here because the pole and its $S$-image are **distinct**, so nothing pinches. **Numerics + prose only, not a lemma** | **open** |
| at $i$ or $\rho$ (elliptic, order $n$) | deform and rescale: $\hat r_f=\lim\eta^{(n-1)/n}\tilde r_{f_\eta}\in W$. Verified $n=3$ at $\rho$ ($k=12$) and $n=2$ at $i$ ($k=20$, to $8\!\times\!10^{-25}$). **Not in the .tex** | **open** |

One-liner: **$r_f\in W$ is solved everywhere; $L^*$ is solved almost nowhere.**

### The elliptic points are not a special case, they are a degeneration

At an elliptic point of order $n$ a pole is never stably *on* the contour: it is the merger
of the **stabiliser orbit**, $n$ poles pinching the arc. $j-j(e)$ has an $n$-fold zero, so
the poles sit at distance $\delta\sim|\eta|^{1/n}$ with $\eta=j_0-j(e)$, and the residues
$1/j'$ blow up like $\delta^{-(n-1)}$. Hence

$$\tilde r\;\sim\;\eta^{-(n-1)/n}\;\sim\;\delta^{-(n-1)},\qquad
\hat r_f:=\lim_{\eta\to0}\eta^{(n-1)/n}\,\tilde r_{f_\eta}\in W .$$

Nothing can go wrong here: for $\eta\ne0$ the poles are off the arc, so `thm:geoperiod`
puts **every** member of the family in $W$ exactly, and $W$ is closed. The branch
ambiguity is an $n$-th root of unity — a sign at $i$, a cube root at $\rho$ — hence an
overall scalar, so the line is untouched.

**Rule: never evaluate at an elliptic point on the contour.** Doing so produced a
spurious $4\mid k$ "obstruction" (`pole_at_i_indented.py`): a fixed indentation $h>\delta$
encircles *both* members of the pinching pair, which is a different homotopy class and not
the continuation of the arc.

Consequence for the draft: `lem:georho`, `lem:Qinv` and `lem:kersum` are all true but
describe the degenerate on-contour object at $\rho$. The deformation limit supersedes
them.

## What each script establishes

- `rf/cohomology_ranks.py` — `lem:kersum` (rank $=\dim V_n^U$ at $k=4..28$;
  $V_n^S\cap V_n^U=0$) and `lem:Qinv` ($\Res_p$ is $U$-fixed; ~1e-28 for modular
  forms vs $O(1)$ for non-modular control weights). Pure linear algebra plus four
  circle integrals.
- `rf/geoperiod_checks.py` — `thm:geoperiod`. The defect factors through $\Phi$ with a
  single universal $D=r_{E_k}|_{(1+U+U^2)}$ across seven unrelated forms.
- `rf/rho_limit.py` — the continuity limit at $\rho$: quadrature floor, Neville
  extrapolation in $u=j_0^{1/3}$, family-independence, and the $E_4^3/(j-j_0)$
  control that must stay bounded.
- `rf/eisenstein_dependence.py` — the limit is **not** intrinsic; it moves with the
  reference form, and the shift is exactly $-\tfrac{432000}{691}\Pi\,r_\Delta$.

## Known open items

- $4\mid k$ double pole at $i$ — the PV of a double pole is not automatically finite.
- The arc-interior indentation result is numerics only; it wants a lemma.
- The continuity limit is in memory and in this directory but **not in the .tex**.
- The two numbers characterising $v$, $\alpha/a=-12.3386299797\,i$ and
  $\beta/a=0.6419871475\,i$, are unidentified. PSLQ finds nothing rational, and they
  are close to but demonstrably not $\Delta$'s own $-12.3395851451\,i$,
  $0.6428727427\,i$.
- All of `Lf/`.
