# Raw run outputs, 2026-09-11/12 (Episode 13)

Verbatim stdout, no commentary. Each run is ~20–70 min, so these are kept rather than
regenerated. Reasoning and predictions live in the generating script's docstring.

| file | script | what it establishes |
|---|---|---|
| `fpscan.txt` | `Lf/finite_part_scan.py` | $\delta$-retraction finite parts at both elliptic points, $P=2,3,4$; leading coefficient vanishes iff $P\equiv1\pmod n$ |
| `bpshift.txt` | `Lf/basepoint_shift_rho.py` | $L^*$ finite and flat to **20 digits** over $\epsilon=10^{-2}..10^{-5}$ at a pole on $\SLZ\cdot\rho$, $P=2,3,4$, $s=7,11$; pole-free control reproduces `plain_arc`; $\Phi(E_{12})=+1$ |
| `connres.txt` | `Lf/connector_residue.py` | short $-$ long $=-2\pi i e^{-i\pi s/2}\mathrm{Res}_\rho(f\tilde k)$, 18 digits, radius-independent |
| `bpW.txt` | `rf/basepoint_W_rho.py` | first (pre-validation) $W$ test of the two connector routings; $O(1)$ defects at $k=12,8$ |
| `harness.txt` | `rf/shifted_harness_validation.py` | shifted harness vs `common.tilde_r` on off-contour poles ($7i$, danger zone, below-arc): agreement $5$–$7\times10^{-16}$ |
| `onarc75.txt` | `rf/shifted_onarc_75.py` | on-arc $S$-pairs (75°, $x=500$) in $W$ at $2\times10^{-20}$; spiral supplies the $S$-symmetric indentation for free |
| `encl.txt` | `rf/connector_enclosing_check.py` | enclosing regime: 11-component difference vector predicted from independent residues to $1.02\times10^{-25}$; residues obey $R(n-\ell)=(-1)^{\ell+1}R(\ell)$ |
| `wind.txt` | `rf/winding_plus_one.py` | four-winding sweep giving $D=-A\,m$ exactly, $A=501060.84529$ |
| `imzero.txt` | `rf/i_pole_m_zero.py` | at $i$: single sides fail (identical raw defects $2744182.69$), the mean ($m=0$) is in $W$ at $10^{-28}$, radius-independent |
| `pksweep.txt` | `rf/m_zero_Pk_sweep.py` | $m=0$ at $(\rho,P{=}2,k{=}4)$, $(\rho,P{=}4,k{=}8)$, $(i,P{=}4,k{=}12)$; $D=-A\,m$ at each; $i$ clean at $10^{-26}$ so the prescription is not $P$-dependent |
| `epscond.txt` | `rf/m_zero_eps_conditioning.py` | the $\rho$ $m=0$ residuals are $\epsilon$-conditioning, not real: clean power laws ($10\times$/decade at $P{=}2$, $10^4\times$ at $P{=}4$) while $\vert D(-1)\vert$, $\vert D(2)\vert$ stay $\epsilon$-flat to 10–12 digits with ratio exactly 2. **Run the shifted geometry at $\epsilon=10^{-1}$, not $10^{-3}$.** |

| `x500.txt` | `rf/x500_discrepancy.py` | diagnoses the $x=500$ gap: wiggle reproduces `arc_eight_classes.py`'s $312512.434941$ exactly (conventions shared), spiral gives $212068.947897$ — **not** the $312510.344$ from the endpoint-clustered mesh. See TRAP 4 in `common.py` |
| `onarcres.txt` | `rf/onarc_resolved.py` | on-arc poles with **pole-centred** clustering. $75°$: wiggle $286185.562742$ and spiral $198819.70067$, **both in $W$** at $1.4\times10^{-28}$ / $4.4\times10^{-29}$ and differing — `lem:arcon`'s free side choice. Also $a_{Sz}=-\overline{a_z}$ at both $x$ values |

`res_at_rho.py`, `fit_sanity.py`, `log_ratio_test.py`, `detour_vs_retract.py` and
`endpoint_vs_interior.py` ran in the foreground; their numbers are quoted in Episode 13.

## 2026-09-15 — residue polynomials $\Xi_{f;z}$, and a bug in `eq:symshift`

| file | script | what it establishes |
|---|---|---|
| `respoly.txt` | `rf/residue_polynomial.py` | residue form, no quadrature. Kernel identity $\mathbf K_T(\tau)-\mathbf K_T(\tau-1)=-P_\tau\vert_{(1-S)}$ to $10^{-41}$; $a_{Sz}=z^na_z$ to $8\times10^{-41}$. **ANTI loop in $W$ at $10^{-41}/10^{-38}$; SYM loop is not** — $k=12,18$, nine CM points ($d=3,4,7,8,11,19,43,67,163$) plus controls. Raw $U$-defect $=2\pi i(1-z^n)r_{E_k}\vert_{(1+U+U^2)}$, so the counterterm IS the $U$-anomaly |
| `respolyq.txt` | `rf/residue_polynomial_quad.py` | same by honest quadrature, no residue theorem. Matches the closed form to $10^{-40}$, radius-independent; **double pole** $E_4^2(j')^2/(j-j(z))^2$ also in $W$; linearity to $10^{-40}$ |
| `symaudit.txt` | `rf/symshift_audit.py` | **`eq:symshift` is wrong.** Winding about a single pole $c_z$ gives $(1+S)$ defect $1.49$, $U$-defect $66$. sym/single $S$-defect ratio is exactly $2.0$ at both $z$, the predicted $\mathfrak a\vert_{(1+S)^2}=2\mathfrak a\vert_{(1+S)}$ |

| `thetalaurent.txt` | `rf/theta_laurent.py` | $\Xi_{f;z}=a_z\Theta(z)$ and **$\Theta(z)=\sum_{m=-1}^{n+1}c_mz^m$ with every $c_m\in W$** — verified out of sample to $10^{-48}$. The 13 coefficients have **rank 3 $=\dim W$**. Functional equation $\Theta(-1/z)=-z^{-n}\Theta(z)$ to $10^{-51}$, i.e. $c_{n-m}=(-1)^{m+1}c_m$. $L(z):=[\Theta+2\pi i(1-z^n)r_{E_k}]/(2\pi i)^{n+2}$ has **rational** coefficients ($\pm1/11$, $\pm2$, $\pm4$, $\pm3/2$, $\pm12$, $\pm42$, $\pm126$), so $r_{E_k}$ touches only $m=0,n$. Directions in $\mathbb P(W)$ move by $10^{-9}$–$10^{-7}$: real, but nearly flat |

| `thetak18.txt` | `rf/theta_k18_DeltaE4.py` | DAM's example $f_z=\Delta E_4\,J'/(J-J(z))\in F_{18}$: $\Xi_{f_z}=\Delta E_4(z)\cdot\Theta_{18}(z)$, one orbit so local $=$ global. Laurent form out of sample to $10^{-48}$; **$c_8=0$** (predicted from the functional-equation midpoint at $k\equiv2\bmod4$), 36 orders below its neighbours; **$\Xi$ surjective onto $W$** — seven $\Theta(z)$ of rank $3=\dim W$ down to tol $10^{-10}$. Odd coordinate $\lambda_-$ has only even powers of $z$, with $\lambda_4/\lambda_2=-25/8$, $\lambda_6/\lambda_4=-26/25$, $\lambda_6/\lambda_2=13/4$ exactly, while $\lambda_2/\lambda_0$ is irrational — $r_{E_k}$ enters only at $m=0,n$ |
| `paritydims.txt` | `rf/parity_dims.py` | $\dim W^{\rm even}=\dim S_k+1$, $\dim W^{\rm odd}=\dim S_k$ at $k=12,16,18,20,24,26$; $\varepsilon$-stability of $W$ to $10^{-29}$–$10^{-36}$; $p_0=X^n-Y^n$ in $W$ exactly and even |

| `zeta17.txt` | `rf/zeta17_check.py` | $r_{E_{18}}/(2\pi i)^{n+1}$: every **odd** coordinate is rational with $43867=\mathrm{numer}(B_{18})$ in the denominator ($-\tfrac{1443183}{7457390}$, $-\tfrac{5586}{43867}$, $-\tfrac{26258}{219335}$, $-\tfrac{5187}{43867}$, palindromic); even coordinates $\ell=2..14$ vanish; $\ell=0,16$ purely imaginary and opposite, so the even part is a single multiple of $p_0$. The $\zeta(17)$ coefficient is $0$ in every odd slot. PART 3: the same coordinates are **irrational** for the $E_4^3E_6$ reference — that is the $\Delta E_6$ contamination |
| `z17hp.txt` | `rf/zeta17_highprec.py` | the $p_0$ coefficient $x=0.1848002737722058226\ldots$ at 30 digits: **not** rational, **not** a rational multiple of $\zeta(17)$, **not** in $\mathbb Q+\mathbb Q\zeta(17)$ — nor for $\zeta(18),\zeta(19),\pi,1/\pi,\zeta(3)$, nor rational over $\pi^{e}$, $|e|\le3$. So $\zeta(17)$ does NOT appear, and $x$ is unidentified. Internal checks: $\Phi(E_{18})=1$ to $10^{-51}$, $r_E[0]+r_E[16]=0$ to $2.9\times10^{-50}$ |

A 16-digit PSLQ in `zeta17.txt`'s predecessor returned `[234241, -4914100, 4870775]` against
$\{1,\zeta(17)\}$. That was **spurious** — a 3-term relation of height $5\times10^6$ needs ~25 digits
and had 16. `z17hp.txt` supplies 50 and finds nothing. Height-vs-precision is the check to apply
before believing any PSLQ hit.

**Open, and predicted:** since $r_{E_{18}}$'s odd part is rational, $\lambda_-$ computed against the
TRUE $E_{18}$ reference should have all coefficients rational over $(2\pi i)^{n+2}$, hence
$\lambda_0/\lambda_2\in\QQ$ — the irrational value in `thetak18.txt` being the $E_4^3E_6$
contamination. Not yet run; needs the 19-point fit redone with the $E_{18}$ reference.

**Two numbers in `thetak18.txt` are WRONG — do not reuse them.** Its PART 6 prints
`dim W^even = 2, dim W^odd = 2` (sum 4 > $\dim W=3$) and its PART 4 prints
`rank of the 19 coefficients c_m = 4`. Both are the same artifact: each row was normalised by its
OWN maximum, so a row that is essentially zero — $c_8$, and the near-zero odd part of a basis
vector — has its roundoff amplified to $O(1)$ and contributes a spurious direction. The true values
are $\dim W^{\rm odd}=1$ and rank $3$, established in `paritydims.txt` with a global normalisation.
The $\lambda_-$ extraction in that same run is NOT affected: it self-validated at $3.7\times10^{-48}$
against $\Theta^{\rm odd}(z)=\lambda_-(z)w_{\rm odd}$, which only holds if $W^{\rm odd}$ is a line.
Same family as TRAP 4: never normalise a vector by its own magnitude when it may be zero.

**Impact on §4 and its numerics.** `cor:symspan`'s *hypothesis* is sound: "$S$-symmetric" there means,
per `lem:georetrace`'s proof, "$S\hat\gamma^S$ is $\hat\gamma^S$ reversed", i.e. $S\gamma=-\gamma$ as a
chain. What is wrong is `eq:symshift`, which writes $v_p$ from $r_S(f,p)$, $r_T(f,p)$, $a_p$ at the
single point $p$ while calling it winding about the $S$-orbit; the $Sp$ term, with the opposite sign,
is missing. Wording to fix too: state the hypothesis as $S\hat\gamma^S=-\hat\gamma^S$, not
"$S$-symmetric".

Not affected, because they never used `eq:symshift`: every direct $W$-membership result
(`onarcres.txt`, `imzero.txt`, `pksweep.txt`, `wind.txt`, validate.py layer H) computes $\hat r_f$ on
an actual contour and tests it. `encl.txt` is also clean — it runs $X_S=0$, $X_T=-1$ (winding on the
$T$-segment alone, so $\hat\gamma^S\neq\hat\gamma^T$), which is outside the `cor:symspan` family and
makes no $W$-claim; it tests `lem:wall_arc`, which is correct.

Still to re-check: whether any §4 statement other than `eq:symshift` silently assumes a one-pole
winding stays admissible — `prop:arcside` and `lem:arcon` change the indentation side at a single
$p$, and $S$ maps outside-at-$p$ to inside-at-$Sp$, so the admissible move is a side flip at $p$
AND $Sp$ together. `onarcres.txt` verified both realised contours land in $W$, so the conclusion
holds; it is the description that needs auditing.

## 2026-09-14 — item D, $\hat r_{f_z}$ against pole location

| file | script | what it establishes |
|---|---|---|
| `wcoords.txt` | `rf/Wcoords_vs_pole_location.py` | $k=12$: $\hat r_{f_z}\in W$ at seven $z$; direction in $\mathbb P(W)$ drifts only $6\times10^{-5}$ |
| `rank.txt` | `rf/rank_of_rSkbang.py` | $\Delta j^m-c_0E_4^3$ at $k=12$: rank 3 at $10^{-18}$, 2 at $10^{-12}$, 1 at $10^{-8}$ — inconclusive as a $\dim r(S_k^!)$ measurement, rows span 12 orders |
| `cm18.txt` | `rf/cm_vs_generic_k18.py` | $k=18$, $g=\Delta E_4$: $\dim W=3$; $\hat r_{f_z}\in W$ at $d=7,8,11,19$ and four non-CM controls; **$\Phi(f_z)=2\pi i\,g(z)$** at every point |
| `dirfloor.txt` | `rf/direction_floor.py` | quadrature floor is $1.1\times10^{-28}$ (depth $10\to12$, fixed dps), so the $\sigma_2/\sigma_1=7.7\times10^{-7}$ spread across $z$ is real; the eight $\hat r_{f_z}$ lie on one complex line to that precision |
| `scalarlaw.txt` | `rf/scalar_law.py` | no product of cusp forms is the scalar: $\Delta^2$ is closest at $0.94\%$, $gh$ at $10\%$ — guessing was the wrong move |
| `jexp.txt` | `rf/jexpansion.py` | **$\hat r_{f_z}=\sum_{m\ge1}v_m\,j(z)^{-m}$, $v_m=-\hat r_{gj^{m-1}j'}$.** $v_1=0$ to $7.7\times10^{-34}$ (since $gj'=-2\pi i E_4^3E_6$ = the reference form); truncation residuals over 11 decades at $d=7,11,19$, decay ratio matching $1387/\vert x\vert$ to 3 digits each |

$\Phi(E_4^3E_6\,j^{m-1})$ came back as the exact integers $1,960,1068480,1309056000$ — so
$\Phi=c_f(0)$ holds on $M^!_{18}$ with pole order up to 4, not just on holomorphic forms.

The effective series radius is $1387$, not the naive bound $\max_{\rm arc}\vert j\vert=j(i)=1728$:
$E_6$ vanishes at $i$, killing the integrand exactly where $\vert j\vert$ peaks.

**Do not read `cm18.txt`'s coordinate ratios or its reality column as properties of the forms.**
`common.nullspace` returns SVD right singular vectors, fixed only up to a unitary of the null
space, so coordinates are comparable *within* one run and meaningless across runs (visible in
`dirfloor.txt`, where dps 30 and dps 40 give unrelated ratios for the same $z$). The reality
figures $0.81493/0.57957$ are identical at all eight points, generic ones included — that is the
complex basis, not CM. The $v_m$ are also mutually parallel to $3\times10^{-6}$, so every
$\hat r_{f_z}$ sits nearly on one line whatever $z$ is.

**CORRECTED 2026-09-15 — why the CM test really found nothing.** It compared coordinate RATIOS,
which are scale-free by construction. For $f_z=g\,j'/(j-j(z))$ the residue polynomial factors as
$\Xi_{f;z}=g(z)\,\Xi^\circ_k(z)$, with $\Xi^\circ_k$ universal (depending on $k$ alone) and rational
in the $E_k$ reference — so ALL arithmetic sits in the scalar $g(z)$, and taking ratios divides it
out exactly. The test therefore discarded the only quantity where CM could live. The CM content is
Chowla–Selberg for $g(z)=\Delta E_4(z)$, a CM value of a weight-$(k-2)$ cusp form, lying in
$\overline{\QQ}\cdot\Omega_d^{k-2}$; the position in $W$ is that scalar times a rational point on a
rational curve, so the period polynomial adds nothing beyond $g$. $\Delta(z)E_4(z)$ is not
contamination — it is the whole arithmetic content, and the rational skeleton is the empty half.

**Do not trust `shifted_onarc_75.py`'s two on-arc rows** — its mesh clustered at the parameter
endpoints while the poles sit mid-interval, so its $|\hat r|$ values and its $2\times10^{-20}$
$W$-membership are artifacts. Superseded by `onarc_resolved.py`. Curiously its numbers reproduced
the *wiggle's* class rather than its own ($286185.562742$ at $75°$), unexplained.

## 2026-09-23 — the $h_\pm$ kernel ambiguity at interior poles

| file | script | what it establishes |
|---|---|---|
| `kamb.txt` | `Lf/kernel_ambiguity_interior.py` | $k=18$: $L^*_+-L^*_-=e^{-i\pi s/2}[C(s)\mathcal D(s)-e^{i\pi(s-1)}C(k-s)\mathcal D(k-s)]$ with $\mathcal D(s)=\sum_m c(-m)m^{-s}+2\pi i\sum_{\mathcal F}\Res(f\,{\rm Li}_s(q))$ — interior pole ($f_z$, $d=7$ and $z=0.2+1.3i$), cusp principal part ($E_4^3E_6(j-744)$), and their sum, all to $\le2.3\times10^{-30}$ at $s=3.3$, $5+1.5i$; exactly $0$ at $s=6$. $h_+-h_-={\rm Li}$ term to $4\times10^{-26}$ |
