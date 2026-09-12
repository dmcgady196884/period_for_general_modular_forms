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

**Do not trust `shifted_onarc_75.py`'s two on-arc rows** — its mesh clustered at the parameter
endpoints while the poles sit mid-interval, so its $|\hat r|$ values and its $2\times10^{-20}$
$W$-membership are artifacts. Superseded by `onarc_resolved.py`. Curiously its numbers reproduced
the *wiggle's* class rather than its own ($286185.562742$ at $75°$), unexplained.
