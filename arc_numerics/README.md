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
| on the arc, $\tau\ne i$ | $S$-symmetric indentation ($\log r$ odd about $\theta=\pi/2$); **numerics + prose only, not a lemma** | **open** |
| at $i$ | PV, $\tilde r_f\in W$ with no $Q_f$. **Simple poles only** ($k=18$, $E_4^6/E_6$); $4\mid k$ **never run** | **open** |
| at $\rho$ (endpoints) | `lem:georho` + `lem:Qinv` + `lem:kersum`; value pinned by the continuity limit, **which is not yet in the .tex** | **open** |

One-liner: **$r_f\in W$ is solved almost everywhere; $L^*$ is solved almost nowhere.**

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
