# Handoff: moving the reference base point to $\rho+1$ and the bottom arc of $\mathcal{F}$

**Session date:** 2026-09-03
**Branch:** `filling_in_proofs` @ `c0cd715` ("fixing the figure and ancillary language")
**File under discussion:** `finite_contour_cocycles_short.tex` (this is the newest copy of it on any branch; `main` and `claude/modular-forms-repo-0FQpy` are both back at `eca0c90`, 2026-08-18)

Nothing in the repo was modified this session except the addition of this file. No numerics were run beyond
the verification described below.

---

## 0. Numbering key (verified against the compiled PDF)

`remark` shares the theorem counter — easy to miscount without it.

| label | number | what |
|---|---|---|
| `def:null_homotopy` | 4.3 | reference homotopy class, contour along $\mathrm{Im}\,\tau=1+\delta$ |
| `lem:wall` | 4.4 | wall-crossing |
| `lem:polylogcf` | 4.12 | $I_{m-1}(s,\tau_p)$, generic |
| `lem:polylogedge` | 4.13 | $I_{m-1}$, edge case |
| `lem:Laurent` | 4.15 | Laurent relations at the $S$-fixed point $i$ |
| `lem:wallvalue` | 4.16 | $L^*(\widehat P_i f,s)$, pole at $i$ |
| `thm:mero` | 4.18 | main explicit formula for $f\in F_k$ |
| `def:Tri` | 4.21 | the region $\mathcal{T}$ |
| `def:Fcirc` | 4.22 | $F_k^\circ$, eq. (81) |
| `lem:rfWmero` | **4.24** | $r_f\in W$ for $f\in F_k^\circ$ |
| `def:Qf` | 4.25 | polar defect polynomial $Q_f$ |
| `lem:rfWFk` | 4.26 | $r_f|_{(1+U+U^2)}=-Q_f$ |

The hard-coded numbers in `figures/T_region.tex:1` ("Definition 4.21") and `figures/verify_regions.py:3,5`
("Definition 4.21", "Lemma 4.26") are **correct**. Do not "fix" them.

---

## 1. Figure 1 (`fig:regions`), for orientation

- **(a)** `figures/T_region.pdf` — $\mathcal{T}$ drawn in $\HH$ against the tessellation. Curvilinear triangle,
  corners $i$, $1+i$, $\tfrac{1+i}{2}$; flat top is $T\hat\gamma^T_f$ (*not* $\hat\gamma^T_f$, which runs
  $i-1\to i$); lower sides the tangent circles $|\tau-\tfrac i2|=\tfrac12$, $|\tau-1-\tfrac i2|=\tfrac12$.
  Straddles $\mathcal{F}$, $T\mathcal{F}$, $S\mathcal{F}$, $TS\mathcal{F}$.
- **(b)** `figures/TF_region.pdf` — $\mathcal{T}_\mathcal{F}=\{\tau\in\mathcal{F}:\mathrm{Im}\,\tau\le1\}$,
  the sliver between the unit arc and $\mathrm{Im}\,\tau=1$, pinching at $i$. Area $\pi/3-1$, i.e.
  $1-3/\pi\approx4.5\%$ of $\mathcal{F}$.

`figures/verify_regions.py` was run and passes clean: `area(T)=0.1415926534` vs $\pi-3$; `area(T_F)` vs
$\pi/3-1$; fold multiplicity histogram `{3: 300}`; `eq:foldT` 0 mismatches on 4000 points.

**Key geometric fact, verified:** $U=TS=\left(\begin{smallmatrix}1&-1\\1&0\end{smallmatrix}\right)$ has fixed
point $\tau^2-\tau+1=0\Rightarrow\tau=\rho+1=e^{i\pi/3}$, and $\rho+1$ is a **strict interior point of
$\mathcal{T}$** (distance to both excluded disks $=0.6197>\tfrac12$). So $U$ is an elliptic rotation of order
3 *about a point inside $\mathcal{T}$* — that is the real source of the 3:1 fold, sharper than the main
caption's "$U$ cycles the corners". The figure source comments (`TF_region.tex:6-7`, `verify_regions.py:25`)
say this correctly.

---

## 2. The proposal

**Take $\tau_0=\rho+1=\tfrac12+i\tfrac{\sqrt3}{2}$ and $\gamma^S=\gamma^T=$ the bottom arc of $\mathcal{F}$
(the unit-circle geodesic $\rho\to\rho+1$ through $i$).**

Legal because `def:segments` (line ~176) asks only for *smooth paths* from $\sigma^{-1}\tau_0$ to $\tau_0$ in
$\HH\setminus\Gamma\!\cdot\!f^{-1}(\infty)$ — homotopy class only, no shape constraint. At $\tau_0=\rho+1$:
$-1/(\rho+1)=\rho$ and $(\rho+1)-1=\rho$, so **both segments have the same endpoints** and may be taken to be
the same path.

### Why both relator loops degenerate

Two facts do everything:

1. $S$ maps the arc **onto itself, reversing orientation** ($e^{i\theta}\mapsto e^{i(\pi-\theta)}$ on
   $\theta\in[\pi/3,2\pi/3]$, fixing $i$).
2. $U=TS$ **fixes** $\rho+1$.

Then, in the loop language of `lem:rfWFk` (which is explicit that *"the loops, not the group elements, carry
the content"*):

- **$S^2$-loop:** $\gamma^S$ ($\rho\to\rho+1$) followed by $S\gamma^S$ ($\rho+1\to\rho$, same arc reversed) —
  retraces itself $\Rightarrow C_{-I}=0\Rightarrow r_f|_{(1+S)}=0$. Same mechanism as `lem:rfWFk`'s $\delta=0$
  argument on the imaginary axis, with the arc substituted.
- **$U^3$-loop:** from $C_U=C_T|_S+C_S$, the $U$-path is $S^{-1}\gamma^T$ ($\rho+1\to\rho$) then $\gamma^S$
  ($\rho\to\rho+1$) — again traversed and retraced. $U$ fixes $\rho+1$, so every $U^j$-transport is another
  retracing loop at the same point. The $(1+U+U^2)$ loop bounds nothing $\Rightarrow$ **$Q_f\equiv0$**.

No contradiction with $\hat C_U=r_f$: the geometric loop $C_U$ vanishes while the coboundary $E|_U-E=-r_f$
carries $r_f$. This is exactly the bookkeeping `lem:cocycle` describes at $\tau_0=i$ (lines 607–608, "$C_S=0$
and $E|_S-E=-r_f$").

So $\mathcal{T}$ is not merely *small* at this base point — $\partial\mathcal{T}$ is degenerate and bounds no
area for any $f$ whatsoever.

### What this buys, precisely

`lem:rfWmero` (4.24) needs **both** conditions of (81). Only the first dies:

1. $f^{-1}(\infty)\cap\mathcal{T}=\emptyset$ — **gone.** Replaced by "no pole on the arc / at $\rho+1$":
   measure zero, principal-value-able as $i$ already is.
2. $\int_{\hat\gamma^T_f}f\,d\tau=0$ — **survives.** It is what forces $\hat C_T=0$ (line ~624); without it
   `lem:rfW` eq. (646) gives $\hat C_U\ne r_f$ and the transfer fails. Becomes $\int_{\text{arc}}f\,d\tau=0$,
   still one linear condition on $F_k$.
   ($S$-stability of the arc only yields $\int_{\text{arc}}f(\tau)(1+\tau^{k-2})d\tau=0$, which is *not* the
   same condition.)

So the claim is **not** "trivially applies to all $f\in F_k$" — it is "applies to all $f$ in a codimension-1
subspace, with no constraint on pole location". But condition 1 is precisely what $\mathcal{T}$, $Q_f$ and the
Brown–Fonseca comparison exist to service, so killing it is the substantive win.

---

## 3. Scope of work

The $q$-expansion / strip structure is load-bearing **only in §4.1** (lines 878–1356). Its own framing says so
(line 883: *"the $q$-series of $f$ diverges on the reference contour as soon as a pole sits at or above
$\mathrm{Im}\,\tau=1$"*), and the projector is literally indexed by contour height,
$\widehat P:=\sum_{\mathrm{Im}\,\tau_p\ge1}\widehat P_{\tau_p}$.

Unaffected by a base-point move:

- **§3's explicit formulas** — $\tau_0$-independent by `prop:indep`; `cor:Sclosed` already gives arbitrary
  $\tau_0$ in closed form. The $q$-expansion there is a derivation device used once at a convenient point, not
  a structural dependence. **For $f\in M^!_k$, moving to $\rho+1$ is substitution, not new analysis.**
  ($(\tau_0/i)^s\to e^{-i\pi s/6}$ at $\tau_0=e^{i\pi/3}$.)
- **§3's cocycle / period-polynomial machinery** (`prop:periodint`, `lem:cocycle`, `lem:rfW`, $W$-basis,
  periods/quasi-periods) — pure contour topology + cocycle law.
- **§4.2, §4.3** — Cauchy on loops.

So the edit is: **one definition** (`def:null_homotopy`, 4.3 — the reference class) **plus roughly four
computational lemmas in §4.1** (`lem:projspace`, `lem:polylogcf` 4.12, `lem:polylogedge` 4.13,
`lem:wallvalue` 4.16), and probably `lem:Laurent` (4.15) as a fifth.

The needed tool is already in the paper. `lem:wallvalue` lines 1242–1244:

> *"As the contour $\mathrm{Im}(\hat\gamma_f^T)$ has height $1+\delta>1$, the analogous step of $w$-inversion,
> $\Li_m(w)=-\delta_{m,0}-(-1)^m\Li_m(1/w)$ from the proof of Lemma 4.12, does not need to be invoked. Thus
> there is no $\tfrac1s+\tfrac{i^k}{k-s}$ in the final result."*

Moving the contour to the floor mostly flips which branch is generic.

### The obstruction trade-off (worth stating in the paper)

The two obstructions run in **opposite** directions with contour height, and are the same geometry read from
two sides — line 1518's $3/\pi\approx95.5\%$ serves both:

- **§4.3 (cohomological, $Q_f$)** wants poles *above* the contour → wants it **low**.
- **§4.1 (analytic, $q$-convergence)** wants poles *below* the contour → wants it **high**.

$\mathrm{Im}\,\tau=1$ is the compromise, and a good one because that horocycle is the unique one tangent to
$\partial\mathcal{F}$ (at $i$) — which is exactly why $\mathcal{T}_\mathcal{F}$ is a pinched two-lobe sliver
rather than a band.

The arc is the extreme point in the §4.3 direction. Price: the arc bottoms out at $\sqrt3/2$, so the
projection threshold drops there, and since every pole orbit has an $\mathcal{F}$-representative with
$\mathrm{Im}\ge\sqrt3/2$, **every orbit now needs projecting** — the clean remainder $(I-\widehat P)f$ shrinks
toward nothing. The projector itself is unaffected (built from $f$'s Laurent data and
$\Li_{1-m}(e^{2\pi i(\tau-\tau_p)})$, nothing to do with the contour); it is the $I_{m-1}$ *evaluations* that
become arc integrals.

---

## 4. DO THIS FIRST

**Mod-3 Laurent relations + wall value at $\rho+1$.**

Right now the reference contour meets **one** elliptic orbit ($i$, disc $-4$) and `lem:wallvalue` exists
solely to handle it. The arc meets **two**: $\rho,\rho+1$ (order 3, disc $-3$) as endpoints and $i$ (order 2,
disc $-4$) at the midpoint. So `lem:wallvalue` acquires a sibling at $\rho+1$.

That sibling is load-bearing in a way the $i$ one is not. `lem:wallvalue` lines 1223–1225:

> *"The limit exists precisely because the pole weights obey the relations of Lemma 4.15 at the $S$-fixed
> point $i$: they cancel the individually divergent contour pieces."*

Finiteness at an elliptic point on the contour is a theorem, not a formality. The mod-3 analogue's shape is
already stated in the paper (line ~1599): the admissible pole orders are $\alpha\equiv k/2\bmod2$ at $i$ and
$\alpha\equiv-k\bmod3$ at $e^{i\pi/3}$.

$\rho+1$ is simultaneously (i) on the contour, (ii) the base point, and (iii) the point whose $U$-loop
degeneration delivers $Q_f\equiv0$. Those three roles coincide there, so the principal-value convention at
$\rho+1$ carries the whole scheme. **If the divergent pieces do not cancel at the order-3 point the way they
do at the order-2 point, that changes the plan rather than lengthening it.** Settle it before doing the
algebra in §4.1.

---

## 5. The BF2025 rationality conjecture (speculative — flagged as such)

Hope: this base point may explain Brown–Fonseca 2025's apparently out-of-nowhere *rational* normalization.

Two threads:

**(a) The prefactor is a $\tau_0=i$ artifact.** `def:Lint` carries $e^{-\pi is/2}$, and
$e^{-\pi is/2}\tau_0^{\,s}=1$ *precisely* when $\tau_0=i$ — chosen to trivialize at the $S$-fixed point. The
$\rho+1$-adapted analogue is $e^{-i\pi s/3}$, leaving $e^{-i\pi s/6}$: twelfth roots of unity where there are
currently fourth.

*Cheap diagnostic, no numerics:* check whether BF's normalization carries disc $-3$ fingerprints
($\sqrt{-3}$, cube/sixth roots of unity, $\Gamma(1/3)$) versus disc $-4$ ($i^k$, $\sqrt{-1}$, $\Gamma(1/4)$).
Disc $-3$ would be strong circumstantial evidence they have implicitly parked at $\rho$.

**(b) The arc is a mirror, and sits inside a rational modular symbol.** $\mathrm{PSL}_2(\ZZ)$ is index 2 in
the $(2,3,\infty)$ reflection group and $|\tau|=1$ is the fixed locus of the anti-holomorphic reflection
$\tau\mapsto1/\bar\tau$ — defined over $\QQ$. Better: the *complete* geodesic containing the arc runs from
$-1$ to $+1$, both cusps, both $\Gamma$-equivalent to $\infty$. So the arc sits inside the modular symbol
$\{-1,+1\}$, and $\Gamma$-translates of it tile that geodesic. Modular symbols with $\mathbb{P}^1(\QQ)$
endpoints are where Manin rationality lives and are the classical substrate of Eichler–Shimura.

**Honest gap:** the contour is a *segment between elliptic points*, and it has not been checked that the
rational structure of the ambient cusp-to-cusp symbol descends to it.

**Note in favor:** if §4 goes through, the contour meets both CM discriminants $-3$ and $-4$. A normalization
compatible with *both* elliptic orbits is far more constrained than one that need only respect $i$ — a better
candidate for something forced rather than merely clean.

### Framing that seems worth keeping

Two canonical contours, each adapted to a different generator, canonical in different senses:

- Horocycle $\mathrm{Im}\,\tau=1$ ↔ $T$ ↔ Hurwitz zeta ↔ Fourier/Mellin. Canonical **metrically** (unique
  horocycle tangent to $\partial\mathcal{F}$). Buys the $q$-expansion, hence §4.1.
- Geodesic $\rho\to\rho+1$ ↔ $S$ ↔ $\tau^{s-1}$ ↔ cycle integrals. Canonical **arithmetically** (reflection
  wall through all three elliptic points, defined over $\QQ$). Buys $Q_f\equiv0$.

Rational normalizations come from arithmetic canonicity, not metric canonicity — and BF's constant is
rational, not merely clean.

---

## 6. Dead ends / corrections from this session (do not re-derive)

- **A height-threshold analysis is the wrong frame.** An earlier line of reasoning generalized
  $\mathcal{T}_\mathcal{F}(h)=\{\tau\in\mathcal{F}:\mathrm{Im}\,\tau\le h\}$ for a *horizontal* contour at
  height $h$, and found it collapses at $h=\sqrt3/2=\min_\mathcal{F}\mathrm{Im}$ (attained only at
  $\rho,\rho+1$), with area
  $2\big(\tfrac\pi6-\tfrac1{2h}-\arcsin\sqrt{1-h^2}+\tfrac{\sqrt{1-h^2}}{h}\big)$ for
  $\sqrt3/2\le h\le1$ — numerically $0$, $0.00045$, $0.0096$, $0.0264$, $0.0472$ at
  $h=\sqrt3/2,\,0.88,\,0.93,\,0.97,\,1$. That is all *correct* but beside the point: paths need not be
  horizontal, and the real mechanism is loop degeneracy (§2 above), which is cleaner and stronger.
- **"$\gamma^T$ is the horizontal chord at height $\sqrt3/2$" is false** — over-specification; `def:segments`
  constrains endpoints and homotopy class only.
- **"The horizontal contour buys every explicit formula in the paper" is false** — only §4.1.
- The main caption's "$U$ cycles $i\to1+i\to\tfrac{1+i}{2}$" is true but weaker than the figure sources'
  "$U$ rotates $\mathcal{T}$ about its interior fixed point $\rho+1$"; the latter is what matters here.

---

## 7. Open questions for the laptop session

1. Do the divergent pieces cancel at $\rho+1$ (order 3) as they do at $i$ (order 2)? — **§4 above, do first.**
2. Are the arc integrals of $\Li_{1-m}$ tractable in closed form (polylog argument inversion, etc.)?
3. Does the $r_f$ produced at $\tau_0=\rho+1$ agree with the reference-class $r_f$ where both are defined?
   On $S^!_k$ it must, by `prop:indep`; on $F_k$ they differ by wall-crossing (`lem:wall`), so the difference
   needs characterizing. This is where Brown–Fonseca remains necessary — not as machinery, but as the
   comparison point telling you whether you have extended the classical $r_f$ or landed on a cousin.
4. Numerical check (deferred by choice this session): compute $r_f$ in both classes for
   $E_{12}(j-j_1)/(j-j(2i))$ (the worked example at line ~1526) with a pole placed inside $\mathcal{T}$,
   confirm $Q_f$ absent at $\tau_0=\rho+1$, and measure the difference. `figures/verify_regions.py` is a
   reasonable harness to build on.
