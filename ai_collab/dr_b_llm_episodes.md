# Candidate episodes — LLM-assisted milestones in the DR-B project

Found via `conversation_search` over this project's chat history (snippets only,
not full reads — treat dates/content as leads to verify, not final).

---

## Episode 1 — PRIORITY: Opus 4.7 reproduces Brown's quasi-periods (40 digits vs Brown's ~20)

**Chat:** [Successful reproduction of Brown's quasi-periods](https://claude.ai/chat/984e5092-8155-4f08-b17d-98145d25a311)
**Updated:** 2026-05-06 (~5.5 weeks before this conversation)

**What the snippet shows:**
- Confirms the Diamantis-Rolen "bottom formula" cocycle $\sigma_f(\gamma)$ is
  mathematically equivalent to Francis Brown's email recipe ($C_f^\gamma$):
  both integrate $f(\tau)(\tau-z)^{k-2}d\tau$ between finite points, no
  regularization, no cusps.
- Walks through the coboundary-modification step (kill $C_T$, read off
  $\eta^\pm$ from $C'_S$) applied to $\hat\Delta = \Delta' = 1/q + O(q^2)$.
- This is the numerical-reconstruction episode: matching/exceeding Brown's
  published $\eta^+ \approx 127202100647.18$, $\eta^- \approx 10276732343.65$
  to higher precision than his stated ~20 digits.

**For Claude Code follow-up:** search any Sage/PARI scripts from around this
date for `eta_plus`, `eta_minus`, `delta_periods`, `cocycle`, `coboundary`,
`Diamantis`, `Rolen` — likely candidates: `delta_periods_*` or
`*BROWN_METHOD*.sage`-style filenames (per earlier chats in this project).

---

## Episode 2 — Generalization to s-plane, Hurwitz kernel, T-gauge-fixing, incomplete gammas (BFK match, finite)

Two candidates, both ~1.5–2 weeks after Episode 1. **Cross-check against your
dated notebook entry for the $\tau_0 = i$ simplification** — whichever of
these postdates that note is the better match for "after I asked Claude to
look at DLMF."

### 2a. [Linear functional and period polynomial coefficients](https://claude.ai/chat/74dee996-f9a6-4262-b7a4-1bca50e52d8f)
**Updated:** 2026-05-17

- Discusses matching the linear-functional coefficients (finite-segment
  contour integration against polynomial kernels) to BFK's incomplete-gamma
  representation of regularized $L$-values.
- Notes that for the holomorphic case (Δ) this matches at all $s$; for the
  weakly-holomorphic case ($\hat\Delta$) the naive matching *fails* — i.e.
  this looks like the point right before the fix.

**ACTION NEEDED: open in browser and read manually.**
`https://claude.ai/chat/74dee996-f9a6-4262-b7a4-1bca50e52d8f`
Claude Code cannot access claude.ai web chats; this one was never verified
beyond the original snippet. The following is the case for why it matters:

**Why 2a is probably the missing episode.** Commit `c30552f` (2026-05-18
12:00 UTC) is the first commit containing the entire-in-$s$ Hurwitz extension
(commit message: "generalize DR-B period polynomial extraction from special
L-values in critical strip... to an entire candidate L-function... Extension
involves Hurwitz-zetas..."). The oldest local Claude Code session (`841ef028`)
starts at 12:25 on the same day and bootstraps by reading `c30552f` as already
committed — so the session that produced the integer-$s$ → entire-in-$s$
generalization is **outside the local log record**. Chat 2a (May 17) sits
in exactly the right window and its described content ("matching... to BFK's
incomplete-gamma... for the holomorphic case... fails for $\hat\Delta$")
is consistent with working through the Hurwitz ansatz. **Check this chat
directly to confirm.**

### 2b. τ₀ = i elegance + Hurwitz unfolding closes the proof (CONFIRMED from local logs)

**Session:** `e50a7aea` (local Claude Code log,
`~/.claude/projects/-Users-dmcgady-Documents-math-meromorphic-modular-forms-quasi-period-success/`)
**Date:** 2026-05-21, ~20:52–21:34 UTC

**What happened, in order:**

1. **You arrived with new untracked files** (~20:52) and said "big news: I think
   we have the kernel of a solid proof that the DR-B L-function for the whole
   complex s-plane should equal BFK."

2. **Claude identified the τ₀ = i observation as the key insight** (~20:55):
   > *"This is a real step forward — the τ₀ = i idea is genuinely elegant."*
   The argument: $S \cdot i = i$, so the S-segment $[-1/\tau_0, \tau_0]$
   collapses to a point. You keep only the T-segment. The dual choice
   (τ₀ → cusp) kills T exponentially. Both are fixed points of $Y(1)$; the
   two choices are the two generators swapped. Claude also flagged a gap: the
   Carlson / polynomial argument only covered $s \in \{1, \ldots, 11\}$, not
   all $s \in \mathbb{C}$.

3. **You asked it integrated into the main .tex** (~20:57); Claude wrote §8.1
   as a sketch, with the mode-by-mode formula and an honest "Conjecture 8.x."

4. **You intuited the Hurwitz → incomplete-gamma identity** (~21:05):
   > *"I bet that there's a standard hurwitz zeta function argument that the
   > integral of one of those fuckers, of the sort defined in Eq. (25) of the
   > note, is just an incomplete gamma — or diff of incomplete gammas."*
   (Eq. 25 being $\widetilde{k}_T(\tau,s) := \zeta(1-s,\tau+1) -
   e^{i\pi(s-1)}\zeta(s-11,\tau+1)$.)

5. **Claude derived the identity in two moves** (~21:07), no DLMF lookup needed:
   - *Move 1 (telescope):* integer-$n$ periodicity of $e^{2\pi i n m}$ collapses
     $\sum_m \int_{i+m}^{i+m+1} e^{2\pi i n u} u^{-\nu} du$ to a single ray
     integral $\int_i^{i+\infty}$.
   - *Move 2 (contour deformation):* substitution $t = -2\pi i n u$ + Cauchy
     deformation to the real axis gives $\Gamma(1-\nu, 2\pi n)$.
   - Result: $\displaystyle\int_i^{i+1} e^{2\pi i n \sigma}\,\zeta(\nu,\sigma)\,
     d\sigma = i^{1-\nu}(2\pi n)^{\nu-1}\,\Gamma(1-\nu, 2\pi n)$.

6. **Phases collapse perfectly** (~21:07): applying the lemma at $\nu = 1-s$
   and $\nu = s-11$ and combining with the $-e^{i\pi(s-1)}$ relative phase,
   $i^{12} = 1$ absorbs all complex factors, and per-mode DR-B = per-mode BFK
   as an exact analytic identity in $s$. Both sides entire; apparent poles at
   $s = 0, 12$ removable because $\int_{i-1}^i e^{2\pi i n \tau} d\tau = 0$.
   **Carlson gap closed; Conjecture became Theorem.**

7. **DLMF cited in the .tex** (~21:34) for the analytic-continuation steps
   (Hurwitz zeta entirety §25.11; incomplete gamma entirety §8.2 / §8.2.7;
   Hurwitz growth bound on vertical strips §25.11.36). DLMF was the citation,
   not the discovery source.

**Disclosure note:** τ₀ = i was your observation (arrived with the new files);
the two-move proof was Claude's derivation in response to your stated intuition
about the Hurwitz → Γ link. Both contributions are verified by the final theorem
and by `verify_multi_s.py`.

---

## Episode 3 — Meromorphic §4: $L^*$ finite + basepoint/path-independent within a homotopy class; wall-crossing (CONFIRMED from local logs)

**Session:** `f3bb5d3c` (local Claude Code log, this project), 2026-06-24 → 07-01
**Paper:** `finite_contour_cocycles_short.tex`, §4 (`sec:mero`)

**What happened, in order:**

1. **You set the frame** that carries the whole meromorphic section:
   > *"we have $L^*(f,s)$ that is finite for all $\tau_0$ and path independent,
   > within a fixed cohomology class."*
   Claude formalized this as **`prop:indepHomotopy`** (the $\partial_{\tau_0}L^*=0$
   argument uses only $T$-periodicity and modularity at the moving endpoints, so it
   holds for all $f\in F_k$ within one chamber of $\HH\setminus\Gamma\!\cdot\!f^{-1}(\infty)$),
   with cross-class change governed by residues.

2. **Your $S$-contour insight** (drove `cor:wall` + `rem:Snonvanish`):
   > *"take heed of the $S$-contour. Even if its end-points are identified/equal,
   > because we're working with a meromorphic case, this does not mean that the
   > $S$-contour vanishes!"*
   Claude wrote `cor:wall` with winding data $X_S,X_T:\Gamma\!\cdot\!f^{-1}(\infty)\to\ZZ$
   and `rem:Snonvanish` (a loop enclosing interior poles $=2\pi i\sum$residues, so
   $\gamma^S$ carries genuine residue data weighted by $\tau^{s-1}$, on the same
   footing as $\gamma^T$ weighted by $\tk$).

3. **You asked for a "skeleton, point-by-point outline for how this hangs
   together"** when the elliptic-pole handling felt like hand-waving. Claude's
   unlock: separate **the object** ($L^*$ = the contour integral, finite by
   `lem:finite`, chamber-constant by `prop:indepHomotopy`) from **the recipe**
   (the incomplete-Γ/polylog series, anchored at $\tau_0=i$). The elliptic point
   is special only because $\tau_0=i$ is where the *recipe* is written; the *object*
   never breaks.

**Disclosure note:** the meromorphic framing (finiteness + path-independence within
a homotopy class; the non-vanishing $S$-contour) is yours; Claude formalized it into
`prop:indepHomotopy` / `cor:wall` / `rem:Snonvanish` / `prop:ellfull` with proofs,
and verified basepoint-independence numerically (`verify_fullL_elliptic.py`, agreement
$\le$ 1e-28).

---

## Episode 4 — The elliptic pole $\tau=i$: parity rule, the polylog-projection flip-flop (your skepticism corrects repeated Claude errors), and the CM-period target (CONFIRMED from local logs)

**Session:** `f3bb5d3c`, 2026-06-30 → 07-01. **Primary-source record (direct transcript).**

**What happened, in order:**

1. **Parity of the pole order at $i$.** You conjectured a pole at $i$ "must be even
   order." Claude proved the correct refinement, now **`prop:ellord`**:
   ${\rm ord}_{\tau=i}(f)\equiv k/2 \pmod 2$ (even iff $4\mid k$, odd iff $k\equiv2\pmod4$),
   via the elliptic coordinate $w=(\tau-i)/(\tau+i)$ conjugating $S$ to $w\mapsto-w$.
   You then offered a **disproof**: *"Consider $E_6/(j-j(i))$. That has a 2nd order
   pole at $i$, and weight $k=6\equiv2\pmod4$."* Claude showed this **confirms** the
   rule: $j-1728=E_6^2/\Delta$ has a *double* zero at $i$, so $E_6/(j-1728)=\Delta/E_6$
   is a **simple** pole (the $E_6$ in the numerator vanishes at $i$). Verified
   `parity_check.py`, `disproof_check.py`.

2. **The flip-flop — you were right, Claude was wrong (repeatedly).** Claude claimed,
   several times and from *flawed* numerics, that the polylog projection *fails* at a
   pole at $i$. You pushed back:
   > *"I REALLY cannot believe that the polylog projection doesn't capture the polar
   > structure at $i$. tf is it doing if it does not?"*
   and
   > *"we don't much care if we do set $\tau_0=i$... We only care if our polylog
   > projections, evaluated along the $S$- and $T$-kernel, behave 'nicely' when there
   > is a pole at $i$."*
   You were correct. Claude's "divergences" were **wrong objects**: (a) the degenerate
   $\tau_0=i$ anchor ($S$ collapses, $T$-alone), and (b) a *single* bare tail $\Li_{-1}$,
   which mismatches the order-1 Laurent coefficient and so carries a spurious leftover
   pole. The **proper full-principal-part projection** $\widehat P_i f = c_2\Li_{-1}+c_1\Li_0$
   gives ${\rm ord}_i(f-\widehat P_i f)=0$ and $L^*(\widehat P_i f)$ **converges** on the
   full $S{+}T$ contour as $\tau_0\to i$ (`proj_clean.py`; the bare-tail divergences
   cancel in the right combination, as $S$ cancels $T$ for the modular form). The
   modular $L^*(f)$ is itself dead-constant to 1e-28 in $\tau_0$ (`modconv.py`).

3. **"What did you integrate??"** Claude floated $c_2=1/432$ as evidence of "pretty
   functions." You caught that this specifies only the RHS:
   > *"when you say 'the matching result came out rational', you're only specifying
   > the RHS. What did you integrate??"*
   Correct — $c_2=a_{-2}(f)/a_{-2}(\Li_{-1})$ is a ratio of **local Laurent coefficients**
   (a small circle around $i$), not the $S{+}T$ $L$-integral. Rational *residue*
   ($a_{-2}(E_4\Delta/E_6^2)=-1/(1728\pi^2)$), no bearing on the $L$-value. Retracted.

4. **CM-period target (the "known thing to match").** You flagged you didn't follow
   the CM/Chowla–Selberg claim: *"I am not at all sure what most of what is going on in
   'the CM period at $i$'... even IS, let alone why it might be (vibe-)true."* Claude
   explained: $\tau=i$ is a CM point (the Gaussian lattice $\mathbb Z+\mathbb Z i$
   admits multiplication by $i$), which rigidifies its periods to a single transcendental
   — the lemniscate constant $\varpi=\Gamma(1/4)^2/(2\sqrt{2\pi})$ (Chowla–Selberg),
   the $\tau=i$ analog of $2\pi$. Proposed hunt (mirroring how the Hurwitz kernel was
   found — *demand the numerics match a known quantity*): PSLQ $L^*(\widehat\Delta_i)$
   against $\{\Gamma(1/4)^a\pi^b\}$.

5. **PSLQ hunt — open, first pass negative.** $L^*(E_4\Delta/E_6^2)$ computed to ~54
   digits at integer $s$ (`pslq_hunt.py`, basepoint check 3.3e-55). It does **not**
   close in the pure period basis $\{\varpi^a\pi^b\}$ (PSLQ returns large junk
   coefficients), so $\Gamma(1/4)$ alone is insufficient; the incomplete-Γ / $e^{-2\pi}$
   transcendentals from the Hurwitz kernel are almost certainly in the ring. No closed
   form yet.

**Disclosure note:** across this episode the *mathematical intuition and error-correction*
were yours — you caught the parity over-claim's true form, drove the projection
correction against Claude's repeated wrong numerics, and caught the local-vs-integral
conflation. Claude supplied the `prop:ellord` proof, the numerics (several initially
wrong, corrected under your pushback), and the CM exposition. The pretty-function hunt
is unresolved as of 2026-07-01. This episode is a candid record of the AI being wrong
and the human being right; it should be represented as such in any disclosure.

---

## Episode 5 — Period polynomials of the elliptic blocks: the $|(1+S)$ / $|(1+U+U^2)$ anomalies, transcendental periods ($\pi$ and $\Gamma(1/4)$), and a canonical-class criterion (CONFIRMED from local logs)

**Session:** `f3bb5d3c`, 2026-07-01. **Primary-source record (direct transcript).**

**What happened, in order:**

1. **DAM proposed the Eichler--Shimura cross-check** --- *"see if any of these results when combined into
   some genuine period polynomial [is] annihilated by $|(1+S)$ and $|(1+U+U^2)$... might point towards an
   actual-canonical definition of which homotopy class is 'the' canonical one."* Claude assembled $r_f$ from
   $L^*(1),L^*(2),L^*(3)$ for $E_4\Delta/E_6^2$ ($k=4$) and found (`period_poly.py`):
   $r_f|(1+U+U^2)=0$ to 5.6e-48 (exact), and **$r_f|(1+S)=\tfrac{1}{864\pi}(X^2+Y^2)$** to 4e-48. This is the
   first appearance of the $X^2+Y^2$ violation of $|(1+S)=0$. Geometry: the pole at $i$ (the $S$-fixed point)
   breaks the $S$-relation by exactly the residue $-2\pi a_{-2}$ ($a_{-2}=-1/(1728\pi^2)$); $f$ is regular at
   $\rho$ ($E_4(\rho)=0$) so the $U$-relation is clean.

2. **"A period $\propto\pi$?"** DAM: *"Whoa, did I read you right? Is it possible that even just one of the
   periods... is propto $\pi$?? genuine transcendental, with no need for fields-medal-level musings."*
   Confirmed: $B(2)=1/(864\pi)$, i.e. $L(E_4\Delta/E_6^2,2)=-\pi/216$. It is elementary precisely because it
   *is* the residue (rational$/\pi^2$ at the CM point, $\times2\pi$). Claude had prematurely written off the
   $\Gamma(1/4)$/Chowla--Selberg guess after seeing only $k=4$.

3. **$k=6$ (DAM's example $\Delta/E_6$) and the parity dichotomy.** PSLQ (`k6_id.py`, 44 digits, clean
   coefficients) gave $B(3)=\varpi^4/(576\pi^4)=\Gamma(1/4)^8/(36864\pi^6)$ --- the CM/lemniscate period. So:
   $k\equiv0\pmod4$ (even-order pole) $\Rightarrow$ residue rational$/\pi^2$ $\Rightarrow$ elementary period;
   $k\equiv2\pmod4$ (odd-order pole) $\Rightarrow$ residue carries $\Gamma(1/4)$ $\Rightarrow$ CM period. DAM's
   Chowla--Selberg instinct was right for $k\equiv2$ (Claude's write-off was premature). The **odd** period
   $\mathrm{Im}\,B(1)$ closed in nothing tried ($\Gamma(1/4),\pi$, Catalan, $\zeta(3),\dots$) at either weight
   --- a genuine open transcendental.

4. **Paper.** New subsection \S4.1 ``Transcendental periods of the elliptic blocks'' (both closed forms, the
   residue/parity mechanism, the ES cross-check). At DAM's request Claude also added the canonical-class
   sentence (below) and the two relation checks explicitly.

5. **"Can we route the residue away?"** DAM: *"Can we find a contour where there is no residue that pollutes
   the $|(1+S)$-vanishing?"* Claude verified (`sdefect.py`): $\Res_i(f\tau)=0$ (the parity relation
   $a_{-1}=i\,a_{-2}$), so the central value is crossing-invariant; the $S$-defect shifts by $2\times$ its
   canonical value per $i$-pole crossing, so it lives on **odd** multiples of $2\pi a_{-2}$ and is **never
   zero**. No contour zeroes $|(1+S)$. But $|(1+U+U^2)=0$ holds in exactly one chamber, so the $U$-relation is
   the clean canonical-class selector; and the $\pi$-period $1/(864\pi)$ is itself crossing-invariant.

6. **The proposition + four examples.** DAM: *"this deserves a Theorem or Proposition: Whenever $f\in F_k$ has
   a pole at either of the elliptic points, the period polynomial will violate the relevant relation
   ($|(1+S)=0$ for poles at $i$, $|(1+U+U^2)=0$ for poles at $\rho$)"* --- and requested
   $E_6/(j-j(\rho))$, $E_8/(j-j(\rho))$, $E_6/(j-j(1.1i))$, $E_8/(j-j(1.1i))$. Claude computed
   (`relations.py`): $\rho$-poles ($E_6/j$, $\Delta/E_4$) $\Rightarrow$ $S$ holds ($\sim$1e-28), $U$ violated
   ($\sim$1.1); $1.1i$-poles $\Rightarrow$ $S$ violated ($3$--$6$), $U$ holds ($\sim$1e-28).
   **Refinement:** $1.1i$ is *not* an elliptic point but lies on the imaginary axis --- the $S$-invariant
   geodesic through $i$ --- so it too breaks $S$. The proposition is really about the $S$-/$U$-invariant loci
   (imaginary axis, $|\tau|=1$ for $S$), with the elliptic fixed points as special cases; broader than
   ``elliptic points only.''

**Disclosure note:** the entire period-polynomial program in this episode --- the ES cross-check, the
``period $\propto\pi$'' observation, the canonical-class question, the ``can we zero the $|(1+S)$ defect''
question, the proposition, and the four-example battery --- is DAM's. Claude executed the numerics
(`period_poly.py`, `k6_id.py`, `sdefect.py`, `relations.py`), supplied the CM/Chowla--Selberg exposition
(with one premature write-off, corrected by the $k=6$ result), and drafted \S4.1. Findings verified to
28--48 digits. The odd period and the precise general form of the proposition (invariant loci vs.\ fixed
points) remain open.

---

## Episode 6 — Explicit $W_\pm$ bases (appendix), the literature cross-check, and the Petersson-norm smash: the finite-contour periods ARE the canonical periods (CONFIRMED from local logs)

**Session:** `f3bb5d3c`, 2026-07-02. **Primary-source record (direct transcript).**

**What happened, in order:**

1. **DAM chose T5** (of a five-item slate) — *"put the explicit bases into an appendix... one equation per basis
   element... $r_{\Delta_k}=\sum_\pm\omega^\pm W_\pm$ (cusp) and $r_{\widehat\Delta_k}=\sum_\pm\eta^\pm W_\pm$
   (weak)."* Claude computed the Kohnen--Zagier rational period polynomials $W_\pm$ for $k=12,16,18,20,22,26$ via
   Haberland-pairing projection (`period_polynomial_bases.py`); the Bernoulli numerators $691,3617,43867,174611$
   are the literature fingerprint. DAM supplied the aesthetic fix: write them in **(anti)palindromic** form
   $W_+=\sum c_\ell(X^{n-\ell}Y^\ell - X^\ell Y^{n-\ell})$, $W_-=\sum c_\ell(X^{n-\ell}Y^\ell + X^\ell Y^{n-\ell})$
   --- halving the length. (DAM: *"stop with the formatting, I can do that."*)

2. **DAM's consistency question** --- *"the period-poly basis decomp at ALL weights looks almost identical, weak
   vs cusp. Consistent with Table 1?"* Claude's comparison (`cmp_table1.py`): the **odd** ratios
   $\eta^-/\omega^-$ match Table 1 to **5--7 digits** at every weight; the near-equality $\omega(\Delta_k)\approx\pm
   \eta(\widehat\Delta_k)$ at high weight (the $k=26$ ratio $\to1$) is REAL and already in Table 1 --- DAM's read
   was right. The **even** ratio $\eta^+/\omega^+$ is off $0.2$--$2\%$ = the $C_T$ coboundary. So the normalization
   ambiguity is localized to the even quasi-period $\eta^+$.

3. **Literature (DAM: "move towards looking up things").** Web + Cohen's ANTS X paper (`cohen_period_paper.pdf`):
   the odd side is Haberland/Petersson-canonical; the even-$\eta^+$ ambiguity is the known BGKO / Paşol--Popa
   *"extra relation on even periods of weakly holomorphic cusp forms"*; Brown--Hain fix it via de Rham.

4. **The Petersson smash (DAM: "yes, please check this!").** The finite-contour periods $r_m(f)=i^{m+1}L^*(f,m+1)$
   must satisfy Haberland's formula [Cohen Thm 5.2(2)], which recovers $\langle f,f\rangle$ from the period
   polynomial. `haberland_petersson_check.py`: at $k=12$ our $L^*$ reproduces **Zagier's
   $\langle\Delta,\Delta\rangle$ to 14 digits**; at $k=16$ it matches the **direct fundamental-domain integral to
   25 digits**. This is a *normalization-independent* validation --- the finite-contour $L^*$ IS the classical
   completed period. Now in the paper as a remark + `eq:haberland` after Table 1, with bibitem `Cohen2013`.

5. **DAM's sharp follow-up** --- *"does this mean $\alpha_k$ is noise?"* Answer: no, but not an invariant either.
   $\det(\mathrm{pd})\in\QQ^\times(2\pi i)^{k-1}$ is Brown--Hain's theorem, so $\alpha_k$ is a definite rational;
   but its value carries basis scaling ($\alpha_k\to\alpha_k/(cd)$) and the $\eta^+$ gauge, neither canonically
   pinned. The $k\ge16$ rationals are gauge-relative, not new intrinsic numbers --- which is why no clean pattern.
   The Petersson check (which sees only the cusp-form periods) is the trustworthy invariant.

**Disclosure note:** the directions in this episode are DAM's --- the appendix scope, the Table-1 consistency
probe, the palindromic form, the "move to literature" push, the Petersson-check idea, and the "$\alpha_k$ = noise?"
question. Claude executed the numerics (`period_polynomial_bases.py`, `cmp_table1.py`, `haberland_petersson_check.py`),
read Cohen's paper, and drafted the remark. Verified: $k=12$ to 14 digits vs Zagier, $k=16$ to 25 digits vs the
direct integral. Open: the canonical $\eta^+$ gauge (whether it makes $\alpha_k$ clean) and $\mathrm{Im}\,B(1)$.

---

## Episode 7 — The winding-class lattice: [E5] falsified, two exact kernel identities, the Bernoulli
realization, and the joint-exactness proposition (PRIMARY-SOURCE: this session)

**Session:** `0fe6428f` (Fable 5), 2026-07-03.  **Deliverable:** `defect_lattice_note.tex/.pdf`.

**What happened, in order:**

1. **Claude's lattice sweep** (T1, DAM's flagship): defects are linear in the winding data
   (cor:wall at integer $s$), so "which class kills which relation" is integer linear algebra.
   Three scratchpad passes (`t1_lattice.py`, `t1_structure.py`, `t1_exact.py`; LS -> PSLQ ->
   Smith NF) found: **[E5] is FALSE** --- the $S$-defect weights of the $i$-orbit windings are
   the integers $(2,3,3)$, $\gcd=1$, so crossing $i$ once and $i{-}1$ once kills the $S$-defect
   (verified 7e-33).  [E5]'s "odd multiples, never zero" swept only $S@i$.

2. **DAM's pushback redirects the method** --- *"There has to be some THINKING through this
   before just churning-out brazillions of lines of code."* And the sharp question: *"why
   doesn't just the straight line from $i{-}1$ to $i$ for the T-kernel and the trivial loop
   at $i$ for the S-kernel land in the trivial space where both relations are satisfied?"*
   Claude's answer, forced by the question, produced the session's two exact identities:
   (i) the kernel reflection $\tk(\tau,k{-}1{-}m)=-(-1)^m\tk(\tau,m{+}1)$, hence
   $\mathbf K_T|(1+S)=0$ identically --- $T$-windings never move the $S$-defect, and the
   $\tau_0=i$ class holds $|(1+S)$ as a THEOREM (upgrading the empirical [E4] row);
   (ii) at integer $s$ the Hurwitz kernel is a **Bernoulli polynomial**
   ($\zeta(-m,x)=-B_{m+1}(x)/(m+1)$), so every winding weight is a residue of $f$ against an
   explicit polynomial --- closed form, no sweeps.  DAM: *"THAT I buy... groks with sums of
   contours I've seen in e.g. KZ1984."*

3. **DAM catches a Claude overstatement.**  Claude had written "there is no class in which
   [the U-defect] vanishes"; DAM cited the intact-U evidence ([E2]: $5.6$e-48).  Corrected
   picture: the defect pair is a point on a lattice coset; [E2] sits on the $U{=}0$ axis,
   $\tau_0{=}i$ on the $S{=}0$ axis.  In the clean class DAM's ORIGINAL proposition (residue
   breaks the resonant relation) is restored and now hand-derived:
   $R_i = ia_{-2}(X^2{+}Y^2)$, so $r_f|(1+S)=2\pi iR_i$ = the measured $-2\pi a_{-2}(X^2{+}Y^2)$.

4. **A DAM misreading, corrected with receipts:** he briefly read Claude as having computed
   period polynomials of polylog-subtracted forms (*"I want the fucking period polynomial of
   the actual fucking MF"*).  All session numerics are direct quadrature of the actual
   meromorphic forms; the subtraction discussion was Claude's answer to DAM's own "why not
   subtract the pole?" (answer: breaks modularity; the modular subtraction localizes the
   defect on the blocks, it does not remove it).

5. **DAM demands a proposition** --- *"I see a lot of examples, but no grand prop... with a
   _claimed_ proof, which one can inspect and interrogate... try and prove this with
   diophantine equation set-ups/solutions."*  Result: Claude hand-derived the remaining
   $T$-columns ($Q_i$ via Bernoulli values at $i$: $7/12$; $Q_{i-1}=Q_i$ a parity identity;
   $Q_{i+1}-Q_i$ contributes exactly $+1$, giving $19/12$) and proved **prop:joint** in the
   note: for $E_4\Delta/E_6^2$ the jointly-exact classes are exactly the coset
   $(1,-1,0;-1,0,1)+\ZZ\langle 4\text{ explicit kernel vectors}\rangle$ --- nonempty, rank 4,
   so ES-exactness alone does not select a contour.  Every weight in the defining equations
   is derived, not fitted.

6. **Honestly open** (note \S7): the $E_6/j$ weight table $(-11,-3,-8/5,-48/5)$ is
   PSLQ-measured, not derived; its joint obstruction (mod-8) is depth-one-window-limited;
   $\Delta/E_4$ and the generic pole are undecided; a free-product $H^2$ heuristic predicting
   "generic poles unobstructed" is flagged speculative; `t1_depth2.py` written, never run.

**Disclosure note:** the falsification method (linear lattice + PSLQ + Smith NF) and the two
kernel identities are Claude's; every course correction is DAM's --- the reason-first
redirect, the $\tau_0=i$ challenge that produced the identities, the intact-U catch that
fixed Claude's overstatement, the KZ1984 anchoring, and the demand for an inspectable
Diophantine proposition.  Claude also wrongly framed [E5]'s failure before DAM's questions
sharpened it into the clean-class/coset picture.  All derived weights were independently
confirmed numerically at 34 digits.

---

## Episode 8 — Rebuilding §4 ($L^*$ for meromorphic $F_k$) statement-by-statement, and the general-$P$ Laurent relations at the elliptic point (PRIMARY-SOURCE: this session)

**Session:** `0fe6428f` (Fable 5 → Opus 4.8), 2026-07-03 → 07-07.
**Paper:** `finite_contour_cocycles_short.tex`, §4 (`sec:mero`).

**What happened, in order:**

1. **The injection protocol.** DAM cleared all of §4 to below `\end{document}` and rebuilt it
   fresh, statement-by-statement: Claude proposes each def/lemma + proof, DAM injects after
   détente ("just propose wording, I'll paste it in --- saves tokens"). Rebuilt in order:
   `prop:indepHomotopy`, `def:class`, `def:null_homotopy` (the minimal homotopy class; basepoint
   lifted to $i(1+\delta)$ for an on-line pole so no principal value is needed --- DAM's design
   call), `lem:wall`, the projection defs (`def:proj`/`cRf`/`Finfty`), `lem:growth`,
   `eq:gc`+`def:IN`, `lem:polylogcf`, `def:tailpieces`, `lem:polylogedge`, `lem:Laurent`,
   `def:FP`+`lem:wallvalue`. Still open: `thm:mero` (the assembly) and a commented `rem:edgevalues`.

2. **DAM's corrections along the way.** (a) *Jargon:* Claude named the three tail pieces "polar /
   Lerch / remainder **layers**"; DAM --- *"what is a Lerch layer? a polar layer?"* --- and the
   christenings were dropped (and "polar" was outright wrong: $\Pi_0$ carries a **log**
   $\Li_1(w)$, not a pole). (b) *A hanging sentence:* `lem:polylogedge` first advertised a
   $w=1$ divergence "vanishing iff $k\equiv2$" and walked away; DAM --- *"why is this left
   hanging. Don't we deal with this?"* --- trimmed to continuity away from $i+\ZZ$, with $z=i$
   handed to the next lemma. (c) *1806 regime:* Claude claimed the new $\mathrm{Res}_{s=0}$ term
   "trivially reduces to 1806's cusp-pole case"; DAM corrected --- *"1806's whole shtick was for
   MFs with poles inside $\mathcal F$ but away from cusps"* --- so it is the **same** regime and
   a genuine consistency check (tracker #10).

3. **A Claude self-catch (def:IN basepoint-dependence).** Writing `lem:wallvalue`, Claude
   realized `def:IN` had over-reached: the tail $\Li_{-N}(e^{2\pi i(\cdot-z)})$ is **not
   modular**, so `prop:indepHomotopy` does not apply and a single $I_N(s,z)$ is
   basepoint-dependent; at $z=i$ the individual tail integral **diverges** as the contour
   tightens onto the pole. Only the whole pole's contribution
   $L^*(\widehat P_i f)=\sum_m r^*(m)I_{m-1}(s,i)$ is finite, its divergences cancelling
   *across* the tails --- which is exactly `lem:Laurent`. Fix: `def:IN` narrowed to
   ${\rm Im}\,z>1$; `lem:wallvalue` restated as the **combined** contribution
   $\sum_m\frac{(-2\pi i)^m a_{-m}}{(m-1)!}\mathrm{FP}_{m-1}$, not a per-tail value. DAM: *"yeah,
   that looks fine."*

4. **THE TECHNICAL CONTRIBUTION --- general-$P$ Laurent relations at $i$.** DAM flagged that
   `lem:Laurent` (first written explicit only for $P=1,2$) *"doesn't deal properly with generic
   pole orders"* and asked *"is it hard to pursue?"* It is not. With $\varepsilon=\tau-i$, the
   key simplification is $-1/\tau-i=-\varepsilon/(1-i\varepsilon)$, hence
   $(-1/\tau-i)^{-1}=i-\varepsilon^{-1}$, and matching the coefficient of $\varepsilon^{-r}$ in
   $f(-1/\tau)=\tau^k f(\tau)$ gives, for $r=1,\dots,P$,
   $$(-1)^r\sum_{m=r}^P a_{-m}\binom{m}{r}i^{m-r}=i^k\sum_{m=r}^P a_{-m}\binom{k}{m-r}(-i)^{m-r}.$$
   Checks: $r=P\Rightarrow\big((-1)^P-i^k\big)a_{-P}=0$ (pole order matches $k\bmod4$);
   $r=1,P=1\Rightarrow(1+i^k)a_{-1}=0$; $r=1,P=2\Rightarrow a_{-1}=\tfrac{i(k-2)}2 a_{-2}$ ---
   all matching the hand-derived low-$P$ cases. **Structural reading:** in the coordinate
   $u=1/(\tau-i)$ the map $S:\tau\mapsto-1/\tau$ acts as $u\mapsto i-u$, so the relation is
   precisely weight-$k$ $S$-covariance of the principal-part polynomial $g(u)=\sum a_{-m}u^m$:
   $g(i-u)=\big[i^k u^{-k}(u-i)^k g(u)\big]_+$. `lem:Laurent` was made general and the
   `rem:higherP` punt deleted --- nothing left to defer.

**Disclosure note:** the §4 rebuild is DAM's protocol, and every structural course-correction is
his --- the "layer" jargon call, the dangling-sentence catch, the 1806-regime correction, and
the *"is it hard to pursue?"* that turned a hedge into a theorem. Claude proposed the
statements and proofs, self-caught the `def:IN` basepoint-dependence mid-derivation (before it
shipped), and derived the general-$P$ Laurent closed form together with its $S$-covariance
reading. The general-$P$ relation is a Claude derivation prompted by DAM's question, verified
against the independently hand-derived $P=1,2$ cases; DAM flagged it as *"the sorta thing that's
hard to vibe, and a technical contribution."* §4 rebuild ongoing (`thm:mero` + `rem:edgevalues`
remain).

---

## Episode 9 — Period polynomials for MEROMORPHIC forms: $r_f\in W$ on a characterised subclass $F_k^{\circ}\subset F_k$ (plus a defect in 1806 and the end-to-end validation suite) (PRIMARY-SOURCE: this session)

**Session:** `0fe6428f` (Opus 5), 2026-08-09 → 08-13.
**Paper:** `finite_contour_cocycles_short.tex`; new `fix_previous_1806_file/`, `end_to_end_numerical_testing/`.

**HEADLINE RESULT (item 6 below).** Eichler--Shimura period polynomials are pushed past the
weakly-holomorphic ceiling. `lem:rfW` gave $r_f\in W$ for $f\in S^!_k$; this session isolates
*exactly* what that proof needs, and the two conditions are checkable geometry rather than
holomorphy. The outcome is a subclass $F_k^{\circ}\subset F_k$ of genuinely **meromorphic**
forms --- poles in $\HH$, not merely at the cusp --- whose period polynomials satisfy the full
Eichler--Shimura relations $r_f|_{(1+S)}=r_f|_{(1+U+U^2)}=0$, together with an explicit
Eisenstein correction covering the constant-mode failure. Injected as `def:Fcirc` and
`thm:rfWmero`. (Adjacent literature treats periods of meromorphic forms via polar harmonic
Maass forms --- Bringmann--Kane, BKV --- so this is the finite-contour framework reaching the
same frontier by different means, not a claim of priority over that line.)

Items 1--5 are the infrastructure and course-corrections that made item 6 reachable and
trustworthy.

**What happened, in order:**

1. **A defect in published work (arXiv:1806.09874, `LemRegX3`).** Chasing DAM's long-standing
   worry about "the naive tension between the 1806 sum-rules and the $1/s$ poles", Claude
   calibrated the lemma's two branches against quadrature. The $y<B$ branch is exact to
   $10^{-32}$; the $y>B$ branch is wrong for $N=0,1,2$, on- and off-axis, by a factor $\sim12$.
   Cause: for $y>B$ the ray crosses $|e(it-\tau_p)|=1$ at $t=y$, where neither the $q$-series nor
   its inversion is valid, so the integral must be split there; the published proof evaluates a
   single antiderivative at $t=B$ and $t=\infty$, dropping the $t=y$ end-point, the $\delta_{N,0}$
   inversion constant (its inversion identity is stated "for $N$ a *positive* integer" while
   `ThmLfn` sums $m>0$, i.e. needs $N=0$), and the entire $\int_y^\infty$ tail. **Nothing DAM
   ever claimed is affected:** Thm 1, the Hurwitz class-number sum rule and every special value
   are *residue* statements, and $1/\Gamma(s)$ annihilates regular terms at $s=0$; the block
   integrals are entire regardless. DAM: *"BAH! Don't care... as it affects NONE of the claims I
   ever gave a shit about."* Corrigendum written to `fix_previous_1806_file/`.

2. **DAM overrules Claude on `lem:Tclosed-Mk` — and is right.** Claude reported the lemma "wrong
   as stated" for ${\rm Re}\,\tau_0<0$. DAM: *"I actually highly doubt this conclusion of yours...
   Absent a clear source for where it might be wrong, I would suspect the numerics."* He was
   correct. For $n<0$, $z=-2\pi in\tau_0$ has ${\rm Im}\,z$ carrying the sign of
   ${\rm Re}\,\tau_0$, so the cut of $\Gamma$ is crossed exactly as ${\rm Re}\,\tau_0$ changes
   sign; the integral is continuous, the *principal* branch jumps. With the standard monodromy
   $\Gamma_{\rm cont}=\Gamma(a)(1-e^{2\pi ia})+e^{2\pi ia}\Gamma_{\rm prin}$ the lemma reproduces
   quadrature to $10^{-31}$. `prop:HGG`'s "the identity continues to $n<0$" was doing real work.
   Claim retracted.

3. **The end-to-end validation suite** (`validate.py` + `validation_notes.tex`, seven layers
   A–G, every check tagged to the `\label` it validates, fast/slow tiers): kernel and unfolding;
   $\tau_0$-independence and both wall-crossings; `thm:weakL`, its residues and the functional
   equation; periods (Haberland–Petersson against direct fundamental-domain integration,
   $6\times10^{-16}$); all of §4; external anchors; and hygiene tests on the suite itself.
   **193 fast / 311 slow, zero failures.** Design lesson recorded in the write-up: a tolerance
   must reflect a test's own numerics, not the working precision — a check limited by
   quadrature or truncation does not improve when mpmath is handed more digits.

4. **Task #10 discharged, twice.** The $1/s$ residue is a property of the *strip the contour
   lies in*: $L^*_S$ is entire, so the pole comes from $\tk\sim-1/s$ times the constant Fourier
   mode, giving $\mathop{\rm Res}_{s=0}L^*=-\int_{\hat\gamma^T_f}f\,d\tau$. Confirmed on a form
   with a pole above the contour, and then on the Hurwitz pair: for $\Lambda_3=E_6/3E_4$ every
   pole lies below and the residue is $-H(3)=-1/3$ (the published sum rule); for $\Lambda_7$ the
   CM point $\alpha_7$ lies *above*, contributes $2\pi ia_{-1}=-1$, and the residue **vanishes**.
   $\Lambda_d$'s $q$-expansion reproduces $H(d)$ and the traces $t_1(3)=-248$, $t_1(7)=-4119$
   exactly. Same machinery, residue $-H(d)$ or $0$ purely by which side of the contour the CM
   point falls on: the whole 1806-vs-note question in one number.

5. **DAM's raised-contour proposal.** *"Could we re-define the reference contour... $\epsilon$-wide
   rectangles which go above all poles... which should essentially match 1806?"* It matches more
   than essentially: lifting $\hat\gamma^T_f$ above every pole makes $L^*$ **exactly** the
   classical $\Lambda(f,s)$ to $10^{-25}$ (above every pole $f$ has a convergent $q$-series in its
   *cusp* coefficients, so `thm:weakL` applies verbatim). But that is also its cost — with no
   poles above the contour there is nothing to project, no inversion, and §4 collapses into
   `thm:weakL`. **The flat class is not an awkward constraint; it is the paper's content.**
   Conclusion: do not switch.

6. **THE RESULT — $F_k^{\circ}$, and period polynomials beyond $S^!_k$.** DAM asked
   whether the raised class gives honest period polynomials. Measurement: $S$ is free in both
   classes (it is the functional equation), $U$ fails in both. Localising the failure to the two
   places `lem:rfW` invokes $f\in S^!_k$:
   - $\hat C_T=-(2\pi i)^{n+1}A_0\int_{\hat\gamma^T_f}f\,d\tau$, which vanishes iff the constant
     Fourier mode *of the strip* vanishes — **not** $c_f(0)$: $f_7=\Delta/(j+3375)$ vanishes at
     the cusp yet has $c=2\pi ia_{-1}(\alpha_7)\ne0$;
   - $\hat C_{\pm I}$, the integral of $f(\tau)(X-\tau Y)^n$ around a closed loop.

   DAM's fix for the first — *"just subtract-off $E_k$ times $c_f(0)$, eh?"* — is exactly right:
   $r_f|_{(1+U+U^2)}=c\,r_{E_k}|_{(1+U+U^2)}$ to $10^{-23}$, so $r_f-c\,r_{E_k}\in W$, and
   $r_{E_k}$ is classical Bernoulli data already imported in `def:Wpm`. It does **not** rescue the
   second: a form with a pole *inside* the loop has $c=0$ and an unrepairable defect (winding
   $-1$, residual unchanged by the subtraction). The loop is the curvilinear triangle bounded by
   the **height-one horocycles at the cusps $\infty$, $0$, $1$** — a $U$-orbit of cusps — with
   sides the $T$-, $UT$- and $S^{-1}$-translates of $\hat\gamma^T_f$ and vertices $i$, $1+i$,
   $\tfrac{1+i}{2}$ at their mutual tangencies. Tangency requires base-point height exactly $1$,
   which is what makes $\tau_0=i$ special. Injected as `def:Fcirc` ($f$ with no poles in
   $\mathcal{T}$ and vanishing strip-constant) and `thm:rfWmero` ($r_f\in W$ for such $f$),
   proved by naming the two hypotheses inside the existing `lem:rfW` argument. **What this buys:
   the period-polynomial theory of §3, previously available only for $f\in S^!_k$, now covers
   meromorphic $f\in F_k$ with poles in $\HH$ --- the objects §4 was built for --- at the price
   of one geometric condition and one Eisenstein subtraction, and with no new machinery
   ($\mathcal{T}$ is built from translates of the existing contour; $r_{E_k}$ is the Bernoulli
   data already in `def:Wpm`).**

**Disclosure note:** the 1806 defect, the residue–strip rule, the suite, the two-obstruction
localisation and the horocycle identification are Claude's; the decisive corrections and the two
structural ideas are DAM's — the `lem:Tclosed-Mk` pushback (where Claude was simply wrong), the
raised-contour proposal, the Eisenstein subtraction, and the standing calls on conventions and
scope. Claude's counterexample bounding DAM's "cusp forms with poles anywhere" reading came after
Claude had itself overstated the result one message earlier. **Transparency point worth recording:
across the suite's construction the apparatus failed far more often than the mathematics — one
genuine defect in published work, none in the note, and nine bugs in Claude's own test code
(ill-conditioned coefficient extraction, a branch folded into a base, a relative-error floor, a
winding orientation, an under-resolved residue, a sampled limit, a $q$-series evaluated outside
its disc, a silently non-applying patch, and a tolerance tracking precision instead of method).
Every red the suite ever showed against the paper was the suite's fault.**

---

## Episode 9 addendum — period polynomials for ALL simple-pole $f\in F_k$ (`thm:rtildeW`), same session, Fable 5 + Opus 5

**The result.** Every $f\in F_k$ with simple poles (none on $\partial\mathcal T$) now has a
canonical period polynomial $\tilde r_f\in W$, built from residue data alone. This is the
paper's closing theorem and the substantive mathematical gain of the session: `lem:rfW` covered
$S^!_k$ only; `lem:rfWmero` covered the $f$ with nothing inside the triangle; `thm:rtildeW`
covers everything, with the correction written out.

After Opus 5's residue formula missed by 15x, Fable 5 located the missing pole via the
$U$-symmetry of the loop (enclosed sets are unions of $U$-orbits; the third member sat outside
Opus's translate window), confirmed the defect identity to $10^{-24}$ three independent ways,
and injected `def:Qf` + the defect identity. DAM then pushed twice, decisively: *"So we cannot
add a polynomial to $r_f$ to 'fix' it?"* and *"WHAT is that unique $w_f$? Should be fairly
straightforward in terms of the residues."* The first push produced the abstract existence
proof (three lines: $\mathrm{im}(1{+}S)=V^S$, self-adjointness under the Haberland pairing,
$V^S\cap V^U=V^\Gamma=0$); the second produced the closed form — Fable's wall-crossing
construction: winding $\hat\gamma^T_f$ about $q=Sp$ (one pole per enclosed orbit) trivialises
the $(TS)^3$-loop at the explicit price of `lem:wall`, and $r_{E_k}$ absorbs the shifted
strip-constant. Result, now `thm:rtildeW` in the note and verified to $10^{-22}$ in both
regimes (suite 203/0): $\tilde r_f = r_f - c_f r_{E_k} + \sum_{\mathcal O} 2\pi i\,
r^*_{f,q}(1)[(2\pi i)^{n+1}\mathbf K_T(q;X,Y) - r_{E_k}] \in W$ for arbitrary simple-pole
$f\in F_k$. Two failed pretty candidates ($\tfrac13 Q|_\Lambda$, $\tfrac12 Q|_{(1-S)}$) are
recorded with their numbers. Attribution: DAM's two questions set both targets; the loop-
trivialisation mechanism and the proof are Claude's; every ingredient in the formula
(`lem:wall`, $\mathbf K_T$, $r_{E_k}$) was already DAM's machinery.

---

### Endgame: DAM catches a vacuous uniqueness claim

The first `thm:rtildeW` said $\tilde r_f\in W$, ``unique modulo $W$''. DAM, reading it cold:
*``tf?? The period polynomial should be unique.''* He was right, and the failure was total, not
cosmetic: the ambiguity lived in exactly the space the theorem claimed to land in, so the
statement carried no information about **which** element of $W$ one gets. Claude had verified
that each choice of representative lands in $W$ and had never compared two choices against each
other. Measured after the challenge: the three poles of a single triple give $\tilde r_f$
differing by $83$--$93\times$ the scale of $r_f$ itself.

The repair is the pairing DAM already uses in `def:periods`: demanding the corrector be
orthogonal to $W$ picks a unique representative, and all three choices then agree to $3\times
10^{-23}$. He also rejected the first phrasing outright --- *``i cannot follow tf ... even
MEANS''* --- forcing the $U$-orbit/representative jargon down to ``the poles inside $\mathcal T$
come in triples; pick one from each triple and set $q=-1/p$'', plus a worked example
(`rem:rtildeexample`) with the actual three poles, the actual residue, and the observation that
the sum has exactly one term. Suite check `e9_explicit_corrector` was rebuilt around
choice-independence rows --- the previous version tested one representative and by construction
could not have caught the bug. Full fast suite 208/0.

**Attribution.** The construction (wall-crossing about $q=-1/p$ trivialising the $(TS)^3$-loop,
$r_{E_k}$ absorbing the shifted strip-constant), the existence proof, and the
Haberland-orthogonal normalisation are Claude's; every ingredient was DAM's pre-existing
machinery (`lem:wall`, $\mathbf K_T$, the pairing, $r_{E_k}$). DAM set both targets with two
questions, then caught the one defect that would have shipped a vacuous theorem. Worth
recording plainly: on this result the model supplied speed and the explicit formula, and the
author supplied the two corrections that made it true and legible.

---

## Episode 10 — The $\Psi$-form subtraction: a canonical period polynomial with no projection, and a literature collision found by asking (PRIMARY-SOURCE: this session, Opus 5 + Fable 5)

**The result.** For $f\in F_k$ with simple poles on one orbit $\mathrm{SL}_2(\ZZ)\,p$, let
$\Psi_p:=E_k/(j-j(p))$ and let $c$ match residues at $p$. Then $f-c\Psi_p$ is holomorphic, so
$r_f-c\,r_{\Psi_p}\in W$ — no projection, no choice of orbit representative, no ambiguity. At
$k=12$, using $E_4^3=j\Delta$ and $E_6^2=(j-1728)\Delta$, one gets
$E_{12}=\Delta(691j-432000)/691$, hence $c=691/(691j_0-432000)$: **rational whenever $j_0$ is, with
no complex multiplication assumed**, and $r_f-c\,r_{\Psi_p}=-c\,r_\Delta$ identically. Rational
$j_0\in(0,1728)$ puts the poles on the unit arc, inside $\mathcal T$ — exact test cases in the
ambiguous locus. Verified at $j_0=200,1000$: naive $U$-defect $7.5$–$7.9$, corrected residual
$1.6\times10^{-23}$.

This supersedes `thm:rtildeW`'s winding corrector and its Haberland-orthogonal normalisation from
Episode 9's endgame. The earlier construction worked but required choosing a representative and
then projecting; this one requires neither.

**Both halves came from DAM's questions, and neither was Claude's idea.** He asked (i) to dig into
Brown–Fonseca for anything usable, and (ii) for an $f$ with poles inside $\mathcal T$ and rational
residues so the subtraction could be reverse-engineered exactly. Those turned out to be one
question: the reference form is Brown–Fonseca's Poincaré series $\Psi^{0,n}$ (their Def. 3.1,
identified with $E_k/(j-j(p))$ in their Ex. 3.13), and the rationality falls out of the $E_{12}$
identity for free.

**The literature check DAM ordered, and what it cost.** Asked to verify novelty before anything
went in a letter, Claude found that arXiv:2508.04844 (Brown–Fonseca, Aug 2025) — already in the
bibliography, cited nowhere in the text — contains `def:Qf` verbatim as their Lemma 5.10, the
residue exact sequence as their (5.7), and a strict generalisation of the session's elliptic
computation as their Remark 3.12 (character orthogonality on $\mathbb{C}(X-wY)^p(X-\bar wY)^q$, both
elliptic points, any $\Gamma$). Their Example 3.13 is the test-form family. What survives is the
$L$-function side: the finite contour, $r_f$ from special values, and the chamber structure —
their residue runs over every pole orbit, $Q_f$ only over poles enclosed by $\mathcal T$. Also
found: their Remark 5.17 concedes the splitting is Hodge-only and **not** of the $\QQ$-MHS, so the
rational normalisation is open in the state of the art, not a hole in this note. Matthes (2101.11491)
and Löbrich–Schwagenscheidt (1907.04024) were checked and do not collide.

**DAM's corrections, five of them, all load-bearing.**
1. *"Numerics have no role, here, in formal proofs and remarks"* — Claude had put a grid-scan
   minimum ($0.49$, unnamed test function) inside a remark. Chasing the objection produced the
   actual proof: at $p=e^{i\pi/3}$, $A^{n-j}B^j$ is a $U$-eigenvector, $1+U+U^2$ is the norm of
   $\ZZ/3$, and $T_i$ dies iff $i\le d=(k/2-1)\bmod 3$, which is exactly the orders the valence
   formula forbids stopping at. The complaint upgraded a measurement to `prop:noQzero`.
2. *"Don't we have an unambiguous $r_f$ for any $f$ which doesn't have poles in the horotope?"* —
   yes, and Claude had not noticed: the corrector sum is then **empty**. This became `cor:unamb`
   ($\ge95.5\%$ of pole positions, area$(\mathcal T)=\pi-3$ exactly). It also killed Claude's
   definition $\tilde r_f:=\Pi_W(r_f)$, since $|\Pi_W(r_{E_k})|=1.28|r_{E_k}|\neq0$ would have
   *moved* the answer in the case DAM had just shown was unambiguous.
3. *"wouldn't the framing be … independent of the reference homotopy class?"* — half right, and
   the measurement (class shift $=1.26\times|\tilde r_f|$) showed the second half false. The
   correct statement is stronger: single windings **leave** $W$, so the class is selected rather
   than stipulated.
4. *"Did you check that there is NO extra residue whose contribution might be in $W$?"* — Claude
   had not; the scan then found one, at $q=\rho=S\,e^{i\pi/3}$, contradicting a claim made one
   message earlier.
5. *"doesn't make it vanish — only in the context of the period polynomial, not the general
   $L$-function at generic $s$"* — Claude had just written that the $\Psi$-subtraction dissolves
   the homotopy-class question. It does so only for the composite object built from $k-1$ special
   values; $L^*(f,s)$ at generic $s$ stays class-dependent, and that dependence is §3's content.

**Claude errors worth recording.** A projection was asserted canonical on the grounds that it
preserves Petersson observables; three independent computations then showed `eq:inner` is *not*
the Petersson pairing (correlation $+0.123$ vs $-0.335$ at $k=24$, confirmed against direct 2D
fundamental-domain integrals to $10^{-10}$) and that Haberland's formula does **not** extend to
meromorphic $f$ (24% gap, converged in every control). A counting argument then closed the route
entirely: $\dim S_k$ constraints cannot pin $2\dim S_k+1$ unknowns. Separately, a guessed character
rule ($T_i$ survives iff $i\equiv k/2$) was falsified by its own table before being used, and a
string-slice bug silently duplicated 536 lines of `validate.py`.

**Attribution.** The $\Psi$-subtraction, its proof, the $E_{12}$ rationality identity, the norm
computation, and the literature findings are Claude's; the two questions that produced the result,
and all five corrections above, are DAM's. The reference form and the cohomological frame are
Brown–Fonseca's, and are now cited as such in `rem:literature`.

---

## Episode 11 — DAM asks for a picture, and the picture picks the base point $\tau_0=\rho+1$ (PRIMARY-SOURCE: this session, Opus 5)

**What happened.** After a stretch of splitting definitions and arguing about the horocyclic
triangle $\mathcal{T}$ in prose, DAM stopped and asked for a drawing: *"You know what I really want?
A friggin drawing/image of $\mathcal{T}$, superimposed on the fundamental domain and its $S$-image
… Don't bother with making it within the narrow range of packages we have in `mmf_venv`. Just make
it."* That produced standalone TikZ figures (`figures/T_region.tex`, `figures/TF_region.tex`,
compiled and rasterised to check label collisions) and a dependency-free checker
(`figures/verify_regions.py`), landing as Figure 1a/1b.

**The picture then did the mathematics.** Drawn rather than described, $\mathcal{T}$ makes one thing
immediate that the prose had buried: its size tracks *how low the reference contour sits*. The
region whose boundary is the $(1+U+U^2)$ loop — the region whose enclosed poles obstruct
$r_f\in W$ — shrinks as the contour is lowered. Push the contour all the way down to the bottom arc
of $\mathcal{F}$ and it collapses: the cohomological obstruction vanishes outright. Hence the base
point $\tau_0=\rho+1$ with both segments the unit arc, which is admissible because
$-1/(\rho+1)=\rho=(\rho+1)-1$ — the two segments have the same endpoints and may be taken to be the
same path. §4.5 now opens on exactly this trade-off: lowering the contour shrinks the cohomological
obstruction and grows the analytic one, and $\mathrm{Im}\,\tau=1$ is the compromise while the arc is
the opposite extreme.

The chronology is in the git log: `c0cd715` ("fixing the figure and ancillary language") precedes
`f116aeb` ("Handoff: base point at rho+1 with the bottom arc of F as reference contour").

**The figure also caught a fact worth checking.** Looking at it, DAM asked: *"With the T-contour of
the reference class being from $i-1$ to $i$, it is slightly surprising that the top of the horo-tope
is from $i$ to $i+1$…"* It is surprising, and it is correct — the top edge is the $T$-translate, and
saying so out loud is only possible once the thing is drawn. Verification, all
dependency-free: $\mathrm{area}(\mathcal{T})=\pi-3$ exactly (0.1415926534),
$\mathrm{area}(\mathcal{T}_{\mathcal{F}})=\pi/3-1$, the forbidden fraction $1-3/\pi\approx4.5\%$, a
Möbius-orbit multiplicity histogram of `{3: 300}` confirming the sides are one path's $U$-orbit, and
`eq:foldT` holding with 0 mismatches on 4000 points.

**Why this is DAM's and not Claude's.** Claude had been reasoning about $\mathcal{T}$ symbolically
for many exchanges without noticing that its area was a *function of the contour height* and could
be driven to zero. The request was for an illustration; what came back was the next base point.
Everything downstream — `thm:geoperiod` ($\tilde r_f\in W$ with no condition on pole location or
order), and then the elliptic degeneration of Episode 12 — sits on that choice.

---

## Episode 12 — Elliptic points are degenerations, not obstructions: deform-and-rescale at $\rho$ and $i$ (PRIMARY-SOURCE: this session, Opus 5)

**The result.** In the geodesic reference class ($\tau_0=\rho+1$, both segments the unit arc), a
pole at an elliptic point is never stably *on* the contour. It is the merger of the **stabiliser
orbit**: $j-j(e)$ has an $n$-fold zero at an order-$n$ elliptic point, so nearby the form has $n$
simple poles at distance $\delta\sim|\eta|^{1/n}$, $\eta:=j_0-j(e)$, pinching the arc from all
sides. The residues $1/j'$ blow up like $\delta^{-(n-1)}$, and

$$\tilde r_f\sim\eta^{-(n-1)/n}\sim\delta^{-(n-1)},\qquad
\hat r_f:=\lim_{\eta\to0}\eta^{(n-1)/n}\,\tilde r_{f_\eta}\in W .$$

Nothing can go wrong: for $\eta\neq0$ the poles are off the arc, so `thm:geoperiod` puts **every**
member of the family in $W$ exactly, and $W$ is closed. The branch ambiguity is an $n$-th root of
unity — a sign at $i$, a cube root at $\rho$ — hence an overall scalar, so the line is untouched.
Verified at $n=3$ ($\rho$, $k=12$, $f_{j_0}=\Delta/(j-j_0)$: growth $4.665,4.653,4.647$ per decade
against $10^{2/3}=4.6416$) and at $n=2$ ($i$, $k=20$, $f_{j_0}=E_4^5/(j-j_0)$, $\eta>0$: growth
$3.179,3.168,3.164,3.163$ against $\sqrt{10}=3.16228$; $\delta=6.35$e-4/-5/-6; $Sz=(1-\delta)i$ to
the digit; $\tilde r\in W$ at 3–7e-25 for every $\eta$, with the $(1+S)$ residual **exactly** $0$ on
the real ray; limit in $W$ at 3.0e-25 and the two rays agreeing to 8.3e-25).

This subsumes `lem:georho`'s $Q_f$, `lem:Qinv`, `lem:kersum`, the principal-value prescription at
$i$ and the apparent $4\mid k$ obstruction, all at once. The board for $r_f$ collapses to: off the
contour, `thm:geoperiod`; on the contour away from $i,\rho$, an $S$-symmetric indentation (which is
legitimate there precisely because the pole and its $S$-image are *distinct*, so nothing pinches);
at an elliptic point, deform and rescale.

**Both moves were DAM's, and Claude resisted the second one with numerics.** (i) Confronted with the
$\dim W$ ambiguity in correcting a pole at $\rho$, DAM asked whether a continuity argument could
remove it, "like what I tried to build when confronted with the ambiguity that Brown–Fonseca found."
That produced the $\rho$ limit. (ii) Shown a measured $(1+U+U^2)$ defect of $2.67$ for an order-2
pole at $i$, DAM refused it — *"Unbelievable that a simple fucking order-2 pole sinks the ship and
pulls $r_f$ outta $W$. I do not believe it for a second"* — and asked the question that resolved it:
*"what about a pole at $(1+\delta)i$ or a pole at $(1-\delta)i$? What happens when $\delta\to0$?"*
Those are not two cases. $z=(1+\delta)i$ has $Sz=(1+\delta)^{-1}i$, so the pair is one orbit and
cannot be separated. Claude had been treating "indent above $i$" and "indent below $i$" as a genuine
choice for several rounds.

**Claude's errors, in order, all caught by DAM or by a check DAM demanded.**
1. Wrote `prop:georhoW` with a reference form $g$ matching $f$'s polar part — but the hypotheses
   admit $g=f$, forcing $\hat r_f=0$. DAM: *"the period polynomial … is just fucking zero?"*
2. Called the extrapolated odd coefficients $(2/21,-25/42,1)$ an independent confirmation. They are
   **forced**: $\dim W=3$ splits as one odd direction and two even, so every element of $W$ has odd
   part proportional to $W_-$. A tautology reported as evidence.
3. Attached "as predicted" to an $O(|j_0|^{1/3})$ convergence rate never checked; the data did not
   support that exponent. DAM caught the tell — *"your numerical confirmations were $10^{-7}$ not
   your usual $10^{-22/23/24}$. We OK dude?"* — and separately caught a claim of statistical
   significance computed from ray-to-ray spread, which is blind to a truncation offset common to
   all rays.
4. Called the $\rho$ limit line canonical before checking that it moves with the reference
   Eisenstein series. It does: $L_{43}=L_{12}-\tfrac{432000}{691}\Pi\,r_\Delta$, verified to
   relative 4.2e-41, so the ambiguity is halved ($2\dim S_k\to\dim S_k$), not killed. Also proposed
   testing this at $\dim S_k=0$, where $\dim W=1$ and the test is vacuous.
5. Predicted the Hadamard finite part would rescue $4\mid k$. It misses $W$ by $O(1)$.
6. Blamed the $k=20$ numbers on under-resolved quadrature. A resolution check showed depth 8 and
   depth 13 agreeing to 1.0e-31 — the numbers were converged all along, and the non-monotonicity was
   component-wise cancellation inside a max-norm.
7. Built indented contours in the wrong homotopy class (a fixed $h>\delta$ encircles *both* members
   of the pinching pair) and reported the resulting defect as a property of $f$.
8. Drafted a `.tex` comment endorsing DAM's deformation idea before testing it. DAM: *"do some
   numerics and try'n see before drafting any shit."*

The one measurement that survived is the crude one: symmetric excision, distrusted at the time,
gives $A=\lim_{\epsilon\to0}\epsilon\,\tilde r_\epsilon\in W$ to 5.8e-5 — the same limit through a
cruder regulator, with $\epsilon$ playing the role of $\delta$.

**Process, also DAM's.** He ordered the arc-class numerics out of ephemeral job scratch and into the
repo (`arc_numerics/`, with `common.py` the single definition of the contour and a
`check_orientation()` assertion, since a silent orientation flip — $\Phi(E_{12})=-1$ — had converted
$r_f-\Phi r_{E_k}$ into $-(r_f+\Phi r_{E_k})$ and read as `thm:geoperiod` failing). He also imposed
a standing rule after a lemma landed in the paper on the back of an approval to *delete* something
else: **no `.tex` edit without an explicit diff shown.**

**Not settled.** The limit is not intrinsic (it moves with the reference form within $E_k+S_k$);
$\rho$ is tested at $k=12$ and $i$ at $k=20$ only; the two numbers characterising $v$ at $\rho$,
$\alpha/a=-12.3386299797\,i$ and $\beta/a=0.6419871475\,i$, resist PSLQ and are close to but
demonstrably not $\Delta$'s own $-12.3395851451\,i$, $0.6428727427\,i$; and none of this is in the
draft. The $L^*$ column in the geodesic class remains essentially empty — items (ii), (iv), (v) and
(viii) of the `%TODO(3)` block.

**Code.** `arc_numerics/rf/{i_limit,rho_limit,eisenstein_dependence,cohomology_ranks,geoperiod_checks}.py`;
`pole_at_i.py` and `pole_at_i_indented.py` are retained and marked SUPERSEDED because their failure
is the instructive part.

---

## Episode 13 — Two manufactured divergences, and the orbifold weight: DAM's skepticism kills $\delta^\kappa$, the base-point shift makes $L^*$ finite at $\rho$, and the prescription is $m=0$ — $1/n$ of the elliptic residue (PRIMARY-SOURCE: this session, Opus 5)

**The result.** $L^*(f,s)$ is finite at a pole on *either* elliptic point of $\gamma^{\rm arc}$, and
neither case needs a regulator. The whole $\delta^\kappa$ apparatus — `lem:ellLstar`, `eq:ellkappa`,
the ray/line models `eq:raymodel`/`eq:linemodel`, `eq:ellvalue`, `eq:ellsub` — was measuring the
divergence rate of a pinch Claude had created, and is not repaired but **deleted**.

*At $i$* the pole is interior to the arc: detour inside or outside, both finite and **exactly
radius-independent** over $r=0.12,0.06,0.02$ (spread 1.4e-31), differing by $2\pi i\,{\rm Res}_i(fK)$
to 4.6e-33, with the mean exactly real. The $\delta^{1-P}$ blow-up came from retracting the **pole**:
$f$ is modular, so a pole cannot move alone, and at an order-$n$ elliptic point the stabiliser acts
as rotation by $2\pi/n$ about a point the contour runs *through*, so every displacement direction
sends one orbit member to each side. §4.1 never hit this because it moved the **path**.

*At $\rho$ the pole sits at both pinned endpoints* of $\int_{\tau_0}^{S\tau_0}$ — a genuinely
different case, and the plain arc really does diverge ($c_{-1}=1.0{\rm e}{-}4+1.2{\rm e}{-}4\,i$ at
$k=4$, $P=2$). DAM's move: displace the base point, as §4.1 did at $i$. It works, and the reason is
a four-line identity — moving $\tau_0\to\tau_0'$, the $S$-substitution ($f|_kS=f$) contributes
$(-1)^{s-1}f\tau^{k-s-1}$, $T$-periodicity plus
$\tk(\tau)-\tk(\tau-1)=-\tau^{s-1}+e^{i\pi(s-1)}\tau^{k-s-1}$ contributes the rest, and the bracket
vanishes identically. **$L^*$ is exactly $\tau_0$-independent**, so $L^*(\epsilon)$ is *constant*,
not merely convergent; $\epsilon=0$ is the one degenerate base point. Verified flat to **20 digits**
across $\epsilon=10^{-2}..10^{-5}$ at $P=2,3,4$ and $s=7,11$, with the noise floor tracking
$\epsilon^{1-P}$ (2e-35, 1.2e-33, 8e-29) exactly as predicted — the guard against quadrature that
steps over the pole and returns identical junk at every $\epsilon$.

Geometrically, displacing $\tau_0$ separates the two segments' far endpoints onto **different points
of the stabiliser orbit** ($V:=T^{-1}S$ fixes $\rho$, $V'(\rho)=\rho$, $VS\tau_0=T^{-1}\tau_0$
exactly), so $K_S$ and $K_T$ pick up different powers of $\omega$. The leading cancellation condition
reduces to $k\equiv2P\pmod6$ — which is `lem:ellLaurent`, automatic — and at $i$ to
$k+2P\equiv0\pmod4$, automatic since $k$ is even. The old $P=1$ argument is the special case
$\omega^{P-1}=1$. The mechanism is uniform in $n$; DAM's worry that the order-3 point might be worse
controlled than the order-2 one was unfounded.

**Both dissolutions were DAM's, from the same instinct, stated twice.** (i) *"there should be no
scaling of the L-function via delta to any power other than fucking zero… Imagine I have a modular
form with two poles. One at $7i$, the other at $i+\epsilon$… WHY THE FUCK would the L-function be
scaled like $\epsilon^{\rm positive}$?? This would OBLITERATE the linearity of the functional."*
(ii) *"was this because even the contour integral blew up in the $\delta\to0$ limit of the old
formulation? Because I am extremely skeptical that this happens…"* and then, when Claude conceded a
real divergence at $\rho$: *"we never did any of that bullshit with the reference contour which ran
from $i(1+\delta)-1\to i(1+\delta)$… (Section 4.1)."* Each time the object was fine and Claude was
sampling it at the one place it isn't defined.

**Claude's errors, in order.**
1. Defined $L^*:=\lim\delta^\kappa L^*(f_\delta,s)$ with $\kappa>0$ and wrote it into the draft.
   **This destroys linearity of $L^*$, and with it everything downstream.** $\kappa$ depends on the
   form — on $P$ and on which elliptic point carries the pole — so $L^*(f+g)\ne L^*(f)+L^*(g)$ for
   two forms with different pole orders. $L^*$ is a *linear functional* on $F_k$; $r_f$ and
   $\hat r_f$ are assembled linearly from its special values; $W$ is a linear subspace and the
   $(1+S)$, $(1+U+U^2)$ relations are linear conditions; `lem:wall_arc` is linear in the residues.
   A form-dependent prefactor invalidates the whole chain, not merely the definition. DAM's
   counterexample was immediate and needed no computation: take $f$ with poles at $7i$ and
   $i+\epsilon$ — the $7i$ contribution varies smoothly in $\epsilon$ and has a finite limit, so
   multiplying the functional by $\epsilon^{\kappa}$ annihilates it. *"This would OBLITERATE the
   linearity of the functional, and would OBLITERATE the contribution from the pole at $7i$… and
   TOTALLY VIOLATES the whole fact that contour integrals of meromorphic functions over finite
   contours are finite, with perhaps residue terms."*
   And the failure is **discontinuous at $\epsilon=0$**, which is worse than non-additive. For every
   $\epsilon\ne0$ the pole at $i+\epsilon$ is not elliptic, so $\kappa=0$, $L^*$ is the ordinary
   finite value, both poles contribute and additivity holds. At $\epsilon=0$ exactly, $\kappa$ jumps
   to $P-1$ and the functional throws away everything but the leading pinch coefficient. So
   $\lim_{\epsilon\to0}L^*(f_\epsilon)\ne L^*(f_0)$: the definition is not even the limit of itself,
   and linearity holds on a punctured neighbourhood of $\epsilon=0$ and fails at the single point —
   precisely the point the definition was invented to cover.
   Claude had instead measured $\kappa$ to eight digits at four configurations. Compute in place of
   thought: measuring the exponent of a divergence is not the same as knowing which object you want.
2. Offered "the $c_0$ drift equals $-c_{\log}\log2$ under a grid halving" as evidence of a genuine
   $\log\delta$. It is an **exact identity of the interpolation scheme** (the two no-log fits are one
   Neville step from the with-log fit), ratio 1.0 in every synthetic case with or without a log. The
   real discriminant is $c_{\log}({\rm grid}/2)/c_{\log}({\rm grid})\to1$ vs $\to2^{-3}$.
   (`fit_sanity.py`.) The conclusion it supported — a scale-dependent finite part at $\rho$ — was
   withdrawn.
3. Recommended "PV at both elliptic points" in one line without checking that $\rho$ is an endpoint
   case. Also called the object a principal value: for $P\ge2$ the symmetric-cut PV **diverges**
   while the detour mean is finite, so the right name is the mean of the two detours (they agree —
   the cut's finite part converges to it at the $\alpha^3$ truncation rate).
4. Floated reality as the canonical selector at $\rho$ when the short connector came out real and the
   long one complex. Dead at $P=3$, where both are real. What survives is one-directional and not a
   proof: short is real in all 6 rows, long in 3 of 6.
5. Reported $\Phi=0$ as possibly structural. $\Phi(f)=c_f(0)$ — every $n\ne0$ mode integrates to zero
   across a period — and all three test forms had $j$ downstairs. A one-glance check of the
   $q$-expansion would have caught it. **Then stated that identity twice too broadly**: it needs *no
   pole of $f$ above the contour anywhere in the strip $|{\rm Re}\,\tau|\le\tfrac12$*, since the
   argument deforms the $T$-segment to a horizontal line at large height and crosses every such
   pole. Two validation rows caught it, both with $c_f(0)=0$: a pole at $0.2+0.99i$ (DAM's "danger
   zone", $|\tau|>1$, ${\rm Im}<1$) gives $\Phi=1.739{\rm e}{-}6-1.085{\rm e}{-}6\,i$, and
   $\Delta/(j-j(e^{i\pi/4}))$, whose orbit sits at $\pm\tfrac12+1.21i$, gives
   $\Phi=-2.740{\rm e}{-}7$.
6. Predicted a surviving $\log\alpha$ at $\rho$ from ${\rm Res}_\rho(fK)\ne{\rm Res}_{\tau_0}(fK)$.
   The discriminant says no log; the mechanism is unexplained.
7. Wrote the wall-crossing sign backwards: ${\rm short}-{\rm long}=-2\pi ie^{-i\pi s/2}
   {\rm Res}_\rho(f\tk)$, i.e. $X_T(\rho)=-1$.

**What the base-point shift costs.** At $\epsilon\ne0$ the segments are distinct paths, so
`def:georef`'s $\hat\gamma^S=\hat\gamma^T$ is gone and the Figure-1(b) region reopens at size
$O(\epsilon)$. $\hat\gamma^T$ needs a connector around $\rho$, and the two routings ($2\pi/3$,
$4\pi/3$) are each flat but differ by **exactly one `lem:wall_arc` crossing** ($X_S=0$, $X_T=-1$),
18 digits at both $s$, against a residue computed independently at three radii — so it is ordinary
wall-crossing, not new machinery.

### Resolution: the orbifold weight. DAM again, and the prescription is $m=0$

**The harness first.** The $O(1)$ defects below were reported before anything had checked that the
spiral+connector machinery reproduces a known answer — DAM: *"I am unconvinced that what you did
establishes what you say it establishes. Second, we have not checked for pedestrian cases."* He was
right to suspect it and the check came back the other way: six cases pass. Off-contour far ($7i$),
off-contour near (the danger zone), below-arc ($e^{i\pi/4}$ — DAM proposed it as on-arc; it is at
$45^\circ$ and the arc is $60$–$120^\circ$), and two genuinely on-arc ($75^\circ$ per DAM's
correction, and $x=500$). All give $\hat r_f\in W$ at $10^{-17}..10^{-20}$ and reproduce
`common.tilde_r` to $5$–$7\times10^{-16}$. The enclosing regime was then checked against
independently-computed residues: the whole 11-component difference vector predicted from small
circles about $\rho$, matching to $1.02\times10^{-25}$. **Bonus structure**: those residues satisfy
$R(n-\ell)=(-1)^{\ell+1}R(\ell)$ exactly — the $s\leftrightarrow k-s$ functional equation — which
predicted the last three rows before they printed.

**The four-point law.** Sweeping the connector winding for $\Delta/j$, $k=12$, $P=3$:
$|D|=501060.84529,\ 1002121.69058,\ 2505304.22645,\ 2004243.38116$ at turns $0,+1,+2,-1$, each to
nine digits, with $\Phi$ exactly linear at $-9.89994{\rm e}{-}6$ per turn. So $D$ is affine in the
winding with zero at turns $=1/3$.

**DAM saw what the $1/3$ was.** *"wouldn't the winding of 1/3 be exactly the right winding that we
should pick up in the contour-integral for the L-function too? I mean, $\rho$ is only 1/3 in the
fundamental domain, which is why e.g. the valence formula has that funny factor of 3 … Not quite
sure what the tension is, yet."* Locally at an order-$n$ elliptic point $\HH\to\SLZ\backslash\HH$ is
$z\mapsto z^n$, so $2\pi/n$ upstairs is ONE loop downstairs. Writing $m$ for the **quotient**
winding, $m=-1+3\,$turns, and the law collapses to
$$D=-A\,m,\qquad A=501060.84529,$$
exact at all four points. The endpoints $S\tau_0$ and $T^{-1}\tau_0$ are one stabiliser step apart
($V=T^{-1}S$ fixes $\rho$, $VS\tau_0=T^{-1}\tau_0$ exactly in $\rm PSL_2$), so paths realise only
$m\equiv2\pmod 3$ and $m=0$ is unreachable by any single contour.

**But $m=0$ is reachable as a weighted average, and $L^*$ is linear.** The prescription:
$$L^*(f,s):=L^*(m{=}0),\qquad
\text{at }i:\ \tfrac12 L^*(-1)+\tfrac12 L^*(+1),\qquad
\text{at }\rho:\ \tfrac23 L^*(-1)+\tfrac13 L^*(2),$$
weights fixed uniquely by $\sum w_j=1$, $\sum w_jm_j=0$. Equivalently: add $1/n$ of the elliptic
residue to the winding-free routing,
$L^*_{\rm can}=L^*_{m=-1}-\tfrac{2\pi i}{n}e^{-i\pi s/2}{\rm Res}(f\tk)$. This is LINEAR in $f$ —
the $1/n$ depends only on the elliptic point's order, never on $P$ or $f$ — which is precisely what
$\delta^\kappa$ was not. It is $\epsilon$-independent, representative-independent, and respects
$s\leftrightarrow k-s$.

**Verified directly at $i$** ($\Delta/(j-1728)$, $k=12$, $P=2$, unshifted base point, both segments
the arc, detour radius $0.1$ and $0.05$ giving identical numbers): the two single sides FAIL —
$|(1+S)|=1.62594$ and $0.534116$, $|(1+U+U^2)|=68.2894$ and $22.4329$, with **identical** raw
defects $2744182.69$ (i.e. $D=-Am$ at $m=\mp1$) — while the mean is in $W$ at $|(1+S)|=8.4$e$-31$,
$|(1+U+U^2)|=2.4$e$-28$. At $i$ the $(1+S)$ relation is the discriminator, because $S$ carries an
inside detour to an outside one so neither single side is $S$-symmetric; at $\rho$ it never was,
the spiral being $S$-symmetric by construction (DAM: *"It only ever was so, dear silconaceuous
comrade"*). And this retroactively explains the earlier finding that the detour mean at $i$ is
exactly real and exactly radius-independent: it is $m=0$.

So the L-function and the period polynomial are both defined at $m=0$ and cannot disagree. The
"architectural problem" reported an hour earlier was an artifact of insisting on a single contour.

**Two more errors, both mine, both caught by arithmetic DAM prompted.**
8. Reported "both routings fail and no winding can work", built a story about finiteness trading
   against $W$-membership on top of it, and credited that story to DAM's own Figure-1(b) framing —
   from comparing max-NORMS rather than vectors, having already noted the norm arithmetic was
   unreliable.
9. Then flipped to "it's a lock, the zero is at an integer" on an exact factor of 2 — which was a
   mislabelling: the old `longway` branch ADDED a turn (the reduced angular travel $d_0=-2\pi/3$ is
   negative), so old-long was turns $=+1$, not $-1$. The $\Phi$ sign was sitting there the whole
   time. Refitting all four points gave $-3$, not $-1$, i.e. the orbifold factor.

**$D(0)=0$ IS GENERAL — $(P,k)$ sweep, `rf/m_zero_Pk_sweep.py`.** Both elliptic points, three pole
orders, four weights. $D=-A\,m$ reconfirmed at each ($|D(2)|/|D(-1)|$ exactly $2$), and $m=0$ kills
the defect: $(i,P{=}2)$ $2.4$e$-28$, $(i,P{=}4)$ $9.3$e$-26$, $(\rho,P{=}2,k{=}4)$ $9.1$e$-14$,
$(\rho,P{=}4,k{=}8)$ $2.4$e$-10$ (both $\rho$ at $\epsilon=10^{-1}$), $(\rho,P{=}3)$ by the
four-point law. At $i$ both single sides fail — including on $(1+S)$ — with **identical** raw
defects, $2744182.69$ at $P=2$ and $171222.5141$ at $P=4$, i.e. $D=-Am$ at $m=\mp1$.
**$P\bmod n$ is irrelevant to the prescription**: $P=4$ at $\rho$ is the class where the old
$\delta^\kappa$ analysis degenerated, and the $i$ control at the same $P$ came back clean.

The one loose number, $(\rho,P{=}4,k{=}8)$ reading $1.9$e$-2$ at $\epsilon=10^{-3}$, is
**$\epsilon$-conditioning, not a defect** (`rf/m_zero_eps_conditioning.py`): clean power laws in
$\epsilon$ — exactly $10\times$ per decade at $P=2$ ($9.1$e$-14$, $8.8$e$-13$, $8.7$e$-12$) and
$10^4\times$ at $P=4$ ($2.4$e$-10$, $2.0$e$-6$, $1.9$e$-2$) — while $|D(-1)|$ and $|D(2)|$ stay
$\epsilon$-flat to 10–12 digits with ratio exactly $2$. A genuine nonzero $D(0)$ could not scale
with $\epsilon$ at all, $\epsilon$ being no parameter of the problem. (Claude proposed
$\epsilon^{-P}$ as the mechanism; it fits $P=4$ and NOT $P=2$, where the rate is one order per
decade, so the exponent's $P$-dependence is unexplained and two points do not fix it.)
**Practical consequence: run the shifted geometry at $\epsilon=10^{-1}$.** Small $\epsilon$ is not
"closer to the limit" — the limit is exact — it only costs precision, and every earlier shifted
number was left several orders on the table for nothing.

**Two parameters, do not conflate them** (DAM caught Claude listing the second as an open tension):
$\epsilon$ is the BASE POINT, $\tau_0=(1+\epsilon)e^{i\pi/3}$ — a contour deformation with the
poles untouched, and exactly flat by the $\tau_0$-independence identity. $\eta=j_0-j(\rho)$ is the
POLE POSITION, Episode 12's family — the abandoned route, since moving a pole drags its stabiliser
orbit across the contour. Episode 12's $\tilde r_{f_\eta}\sim\eta^{-2/3}$ divergence is therefore
no tension with the finite $\hat r_{f_0}(m{=}0)$: one describes a family of DIFFERENT forms as
$j_0\to0$, the other describes the single form at $j_0=0$. $\hat r$ is simply not continuous in
$j_0$ there (the pole order jumps from three simple to one triple) and nothing requires it to be.
Episode 12's $v$ is **superseded**, not competing: it was a workaround for a question that now has
a direct answer.

**ON-ARC NON-ELLIPTIC POLES, and a mesh trap that voided a claim.** `rf/shifted_onarc_75.py`
clustered its quadrature mesh at the PARAMETER endpoints $t=0,1$ — the arc's ends — while an
on-arc pole at angle $\phi$ sits at $t=(2\pi/3-\phi)/(\pi/3)$, mid-interval. At $\epsilon=10^{-3}$
the spiral passes $5\times10^{-4}$ from such a pole against a panel width $0.05$, so the
quadrature never resolved it and returned a smooth wrong number that still looked like
$W$-membership at 2e-20. **Verbatim the trap already recorded for `lem:arcdeform`**, second
occurrence, now TRAP 4 in `common.py`. NOT affected: the off-contour validations, and — because
$\rho$ and $\rho+1$ ARE the parameter endpoints — none of the elliptic work.

Redone with pole-centred clustering (`rf/x500_discrepancy.py`, `rf/onarc_resolved.py`): the
hand-indented wiggle contour $\log r=0.15\sin(6(\theta-\pi/2))$ reproduces
`arc_eight_classes.py`'s recorded $|\hat r|=312512.435$ **to all digits**, so conventions are
shared across scripts and the $x=500$ gap was never normalisation. At $75^\circ$ (DAM's angle)
wiggle gives $286185.562742$ and spiral $198819.70067$, **both in $W$** at $1.4$e$-28$ and
$4.4$e$-29$, differing — i.e. `lem:arcon`'s free side choice, NOT the stronger "the shifted spiral
picks the side for free" that Claude claimed off the broken mesh. At a non-elliptic pole the
stabiliser is trivial, so there is no $2\pi/n$ to weight and no endpoint pinning to define $m$
against; the side is genuinely free and the draft already says so in case 2 of `thm:arcperiod`.

Two incidental facts: $a_{Sz}=-\overline{a_z}$ for the $S$-pair at both $x=500$ and $75^\circ$
(hence $a_z+a_{Sz}$ purely imaginary, $a_z-a_{Sz}$ purely real), and $\Phi({\rm wiggle})=
2\pi i\,a_z$ exactly at $x=500$. The wiggle/spiral difference is NOT a minimal winding — the four
$(X_z,X_{Sz})\in\{\pm1\}^2$ combinations give relative $1.49$, $0.257$, $2.26$, $1.47$, with
$(+1,-1)$ closest but not matching — unsurprising, since the wiggle swings 16% off the circle and
likely crosses other points of the pole orbit. DAM's call: drop it, two arbitrary classes have no
reason to differ minimally. Also unexplained: the broken-mesh run reproduced the WIGGLE's value
($286185.562742$) rather than its own class.

**Still open.** §5.2 (comparison to [1806]) is still a placeholder, and **none of this is in the
draft**. The deletion list is settled (`lem:ellLstar` and its apparatus, `def:retract`,
`def:deltafam`, `def:arcsplit`'s $\delta^\kappa$, the $\kappa$ language at 1953–1966, 2033–2034,
2473–2479, `lem:arcell`'s dead references to `def:ellfam`/`eq:ellkappa`, case 3 of
`thm:arcperiod`); the additions are two — the $\tau_0$-independence identity promoted from the
assertion at line 1885 to a proved lemma, and the $m=0$ prescription. Per DAM's standing rules the
§5 skeleton gets agreed before drafting and every diff shown. No `.tex` was edited this session.

**Code.** `arc_numerics/Lf/{detour_vs_retract,endpoint_vs_interior,basepoint_shift_rho,connector_residue,res_at_rho,fit_sanity,log_ratio_test,finite_part_scan}.py`;
`arc_numerics/rf/{basepoint_W_rho,shifted_harness_validation,shifted_onarc_75,connector_enclosing_check,winding_plus_one,i_pole_m_zero}.py`.
Raw outputs in `arc_numerics/logs/`.

---

## Episode 14 — Residue polynomials: DAM's skepticism finds a real bug, and one pole's winding turns out to be enough (PRIMARY-SOURCE: this session, Opus 5)

**The result.** For $f\in F_k$ and a pole $z\in\HH$ of **any order**, the closed-contour integral

$$\Xi_{f;z}:=\oint_{c_z-Sc_z}f(\tau)\,\mathcal K(\tau;X,Y)\,d\tau,\qquad
\mathcal K:=(2\pi i)^{n+1}\big[(X-\tau Y)^n+\mathbf K_T\big]-r_{E_k}$$

lies in $W$. No base point, no homotopy class, no period polynomial in the construction — so the two
ambiguities that §4.3 could not remove, both of which are statements about *open* contours, do not
touch it. It vanishes identically on $M^!_k$ and at $z=i$ (where $Si=i$ kills the cycle). New §5.

**How it started.** DAM, on the §4.4 draft: *"I really am pretty skeptical of this core claim that the
proper winding of the S- and T-contours of a single fucking pole ... is enough to give an element in
$W$."* Asked for (1) a direct proof via the Hurwitz kernel on the $T$-contour and the standard kernel
on the $S$-contour, and (2) an explicit computation at the class-number-one CM points.

**The proof, and what the Hurwitz kernel actually is.** Two inputs. *(A)* $\int_{\gamma c}fP_\tau
d\tau=(\int_cfP_\tau d\tau)|_{\gamma^{-1}}$, from $f|_k\gamma=f$ and $k-n-2=0$. *(B)* At $s=\ell+1$,
$\tk(\tau,\ell+1)=\zeta(-\ell,\tau+1)-(-1)^\ell\zeta(\ell-n,\tau+1)$, and $\zeta(-m,x)$ regularises
$\sum_{j\ge0}(j+x)^m$, so $\mathbf K_T(\tau)=\sum_{m\ge1}P_{\tau+m}|_{(1-S)}$ — **the $T$-kernel is
the regularised $T$-orbit sum of the $S$-kernel**, its only algebraic content being the Hurwitz
recursion $T\mathcal T=1+\mathcal T$. With $\mathfrak a:=\oint_{c_z}fP_\tau d\tau$ (which is $2\pi i$
times Brown–Fonseca Lem. 5.10's residue polynomial, so *any* Laurent data) the bracket collapses to
$\mathfrak a|_{(1-S)E}$, $E=1+\mathcal T(1-S)$, and the $S$-relation is one line:
$(1-S)E(1+S)=(1-S^2)+(1-S)\mathcal T(1-S^2)=0$, since $S^2=\pm I$ acts trivially for $n$ even. The
$U$-relation is where $r_{E_k}$ earns its place — the loop shifts $\Phi$ by $\oint f\,d\tau$, breaking
the $\hat C_T=0$ hypothesis of `lem:rfW`'s third display — which answers DAM's earlier question
*"what does the closed contour integral of $f$ around a pole at $z$ need to know about $\Phi,r_{E_k}$?"*
Measured: raw $U$-defect $=2\pi i(1-z^n)\,r_{E_k}|_{(1+U+U^2)}$. **$r_{E_k}$ is the $U$-anomaly, not a
bolt-on.**

**The bug the skepticism found.** The draft said "a simple positively oriented loop enclosing
$\{z,Sz\}$" — the *symmetric* cycle. But `lem:georetrace` makes the relator paths retrace via
$S\hat\gamma^S=-\hat\gamma^S$, so a winding added to both segments must **reverse** under $S$. And
$E(1+S)=(1+S)$, so the symmetric loop's $S$-defect is $\mathfrak a|_{(1+S)^2}=2\mathfrak a|_{(1+S)}$
— exactly twice the single loop's. Measured at two $z$: ratio **2.0**, not approximately. Single loop
$(1+S)$ defect 1.49, $U$-defect 66; symmetric 0.36 / 16.0; antisymmetric $2$e$-41$ / $3$e$-38$.

DAM then asked the right follow-up: *"if the contour suggested by section 4 corollary doesn't match
the actual contour, this might point towards a bug in section 4.3, eh?"* It does. `cor:symspan`'s
*hypothesis* is sound (its "$S$-symmetric" means, per `lem:georetrace`'s proof, "$S\hat\gamma^S$ is
$\hat\gamma^S$ reversed"), but **`eq:symshift` is wrong**: it writes $v_p$ from $r_S,r_T,a_p$ at the
single point $p$ while calling it winding about the $S$-orbit. The $Sp$ term, opposite in sign, was
missing; a spurious $i^{\ell+1}$ was also present, having already cancelled against `def:Lint`'s
$e^{-i\pi s/2}$. Both fixed. Nothing verified by direct contour computation is affected —
`encl.txt` runs $X_S=0,X_T=-1$, outside the `cor:symspan` family entirely, and tests `lem:wall_arc`,
which is correct. Flagged but not fixed: `prop:arcside`/`lem:arcon` flip the indentation side at a
single $p$, and $S$ carries outside-at-$p$ to *inside*-at-$Sp$.

**Numerics.** Two independent routes agreeing to $10^{-40}$. *Residue form* (`respoly.txt`): exact
Bernoulli algebra, $k=12,18$, all nine class-number-one CM points $d=3,4,7,8,11,19,43,67,163$ plus
generic controls; kernel identity to $10^{-41}$, $a_{Sz}=z^na_z$ to $8$e$-41$. *Honest quadrature*
(`respolyq.txt`): trapezoid on one circle using $\oint_{Sc_z}F=\oint_{c_z}F(-1/w)w^{-2}dw$, no
residue theorem — matches the closed form to $10^{-40}$, radius-independent, **double pole also in
$W$**, linearity to $10^{-40}$. Audit in `symaudit.txt`. New end-to-end layer I, 19/19, carrying the
negative control: $W$-membership of the antisymmetric loop alone would not have caught the bug.

**Explicit form.** At a simple pole everything is Bernoulli polynomials: $\Xi_{f;z}$ is
$a_z(2\pi i)^{n+2}\times$(polynomial in $z$ over $\QQ$) minus $2\pi i\,a_z(1-z^n)r_{E_k}$, so its only
transcendental content beyond $a_z$ is $r_{E_k}$, with coefficient $(1-z^n)$ — `eq:respolyexp`.

**Novelty.** Brown–Fonseca arXiv:2508.04844 Lem. 5.10 computes the same local residue polynomial, and
their residue sequence (5.7) runs class $\to$ local data with a splitting (Cor. 5.16) landing in
$H^1_{dR}(U_\Gamma)$ and satisfying ${\rm Res}\circ s={\rm id}$. $\Xi$ runs the other way and lands in
$W\cong\ker({\rm Res})$, so it is not their splitting; it looks like the Betti counterpart of their de
Rham sequence. Not found elsewhere in the scan; cannot be ruled out as a known map in disguise.

**Also this session.** `prop:Jexp`: $\hat r_{f_z}=-\sum_{m\ge1}\hat r_{hJ^{m-1}J'}J(z)^{-m}$ for
$|J(z)|>984$, with $v_1=0$ exactly when the reference is $hJ'/c_0(hJ')$ — which killed the CM-arithmetic
question (the $v_m$ are parallel to $3$e$-6$, so every $\hat r_{f_z}$ sits on one line regardless of
$z$). And the uniqueness of $\tau_0=\rho+1$ written into `def:georef`: the two segments share endpoints
iff $-1/\tau_0=\tau_0-1$, whose only root in $\HH$ is $\rho+1$, which is also $U$'s fixed point; $T$
admits no counterpart, being parabolic. DAM: *"the ONLY reference contour which is $\SLZ$ invariant, in
the sense that it exactly divides the fundamental domain from its $S$-image."*

**Code.** `arc_numerics/rf/{residue_polynomial,residue_polynomial_quad,symshift_audit,jexpansion,direction_floor,scalar_law,cm_vs_generic_k18}.py`;
`end_to_end_numerical_testing/validate.py` layer I. Raw outputs in `arc_numerics/logs/`.

---

## Episode 14 addendum — the $E$-operator: period polynomial and residue polynomial are one functional on two cycles (PRIMARY-SOURCE: same session, Opus 5)

**The unification.** Put $A(c):=\int_c f(\tau)(X-\tau Y)^n\,d\tau$ for a $1$-chain $c$, let
$\mathcal T:=\sum_{m\ge1}T^{-m}$ (Hurwitz-regularised, the only relation being $T\mathcal T=1+\mathcal T$),
and set
$$E:=1+\mathcal T(1-S).$$
Then, with $\gamma^S=\gamma^T=\gamma^{\rm arc}$,
$$r_f=(2\pi i)^{n+1}A(\gamma^{\rm arc})\big|_E ,\qquad
\Xi_{f;z}+\text{(counterterm)}=(2\pi i)^{n+1}A(c_z-Sc_z)\big|_E .$$
**Same functional, different cycle.** And the $S$-relation needs only $Sc=-c$: that gives
$A(c)|_{(1+S)}=0$, hence $A(c)|_{E(1+S)}=A(c)|_{(1+S)}+A(c)|_{\mathcal T(1-S^2)}=0$. Both the arc
(`lem:georetrace`) and the antisymmetrised pole loop satisfy $Sc=-c$, so one argument covers both.

**What the Hurwitz kernel IS, geometrically.** Lemma~B read backwards through
$A(\gamma c)=A(c)|_{\gamma^{-1}}$ gives
$\int_{\gamma^T}f\mathbf K_T\,d\tau=\big[\sum_{m\ge1}\int_{T^m\gamma^T}fP_\tau\,d\tau\big]\big|_{(1-S)}$,
and $\sum_{m\ge1}T^m\gamma^T$ is the chain $\tau_0\to T\tau_0\to T^2\tau_0\to\cdots$ --- in
$\Gamma\backslash\HH$, the loop about the **cusp** traversed infinitely often. So the Hurwitz
$T$-kernel is the $\zeta$-regularised winding about the cusp, and $\Xi$ is to an interior pole what
$r_f$ is to the cusp: the cusp winding is infinite and needs regularisation, the pole winding is
finite and needs none. This is the framing to put to Brown; it closes the loop on the history DAM
describes --- Brown's cohomological quasi-periods for $\Delta'\in S^!_{12}$, then the Hurwitz
polynomials reverse-engineered at integer $s$, then the generic-$s$ kernel.

**The $W^\pm$ localisation, §5.1.** For $f_z=\Delta E_4\,J'/(J-J(z))\in F_{18}$ (one pole orbit, so
local $=$ global; $\dim S_{18}=1$, so $W$ has a cuspidal part), $\Xi_{f_z}=\Delta(z)E_4(z)\,
\Xi^\circ_{18}(z)$ and, in the `def:Wpm` basis whose supports are disjoint,
$$R(z)=\frac{\mu_+(z)}{\mu_-(z)}=43867\,\frac{P_+(z)}{P_-(z)},$$
$P_\pm$ coprime-integer (`eq:mum`/`eq:mup`); the $43867=\mathrm{numer}(B_{18})$ is a universal
prefactor since $\mu_-$ carries $2/131601$ and $\mu_+$ carries $2/3$. `lem:Wpmloc`: $R(-1/z)=R(z)$
--- the localisation sees only the $S$-orbit --- $R(-z)=-R(z)$, and $R$ is purely imaginary exactly
on the imaginary axis and $|z|=1$. Exactly, at class number one:
$R(\tfrac{1+\sqrt{-3}}2)=\tfrac{43867}{8129}\sqrt{-3}$,
$R(\sqrt{-2})=\tfrac{43867\cdot18}{119513}\sqrt{-2}$ (both on reflection loci, hence pure), and e.g.
$R(\tfrac{1+\sqrt{-7}}2)=\tfrac{43867(-107+316769\sqrt{-7})}{3931333248}$. At $z=i$ both $\mu_\pm$
vanish --- the algebraic shadow of $\Xi_{f;i}=0$.

**DAM's correction, and it reframes the whole CM story.** Claude called the irrational
$\lambda_2/\lambda_0$ "$\Delta E_6$ contamination" and reported the rationality of $\Xi^\circ_{18}$
as the prize. DAM: *"don't lose the forest for the trees. $\Delta(z)E_4(z)$ is not contamination."*
Right, and deeper than the word: $\Xi^\circ_k$ depends on $k$ alone, so proving it rational is
proving the geometric skeleton carries no arithmetic. ALL information about $f$ sits in the scalar
$g(z)=\Delta(z)E_4(z)$. That also finally explains why `cm_vs_generic_k18.py` found nothing: it
compared coordinate RATIOS, which are scale-free and divide $g(z)$ out exactly --- it discarded the
only quantity where CM could live. The CM content is Chowla--Selberg for $g(z)$, i.e.
$g(z)\in\overline\QQ\cdot\Omega_d^{16}$; the residue polynomial adds nothing on top. Corrected in
`memory/j_expansion_of_rhat.md` and `arc_numerics/logs/README.md`, both of which had recorded the
wrong reason.

**$\zeta(17)$: a clean negative.** Every odd coordinate of $r_{E_{18}}/(2\pi i)^{17}$ is rational
with $43867$ in the denominator; even coordinates $\ell=2..14$ vanish; the even part is a single
multiple of $p_0$. That much is standard $\zeta(s)\zeta(s-k+1)$ bookkeeping, as DAM noted, and is a
consistency check rather than a find. The $p_0$ coefficient $x=0.1848002737722058226\ldots$ is at 30
digits NOT rational, NOT $\QQ\cdot\zeta(17)$, NOT in $\QQ+\QQ\zeta(17)$, and likewise for
$\zeta(18),\zeta(19),\pi,1/\pi,\zeta(3)$ and $\pi^{e}$, $|e|\le3$. A 16-digit PSLQ had returned
$[234241,-4914100,4870775]$; that was spurious (height $5\times10^6$ needs ~25 digits, had 16) and
is the cautionary case for height-vs-precision.

**Open questions, ranked.** (1) The two-variable kernel
$\langle\Xi^\circ_k(z),\Xi^\circ_k(w)\rangle$ --- universal, rational, exactly computable from the
$c_m$; Brown--Fonseca's paper is about matrix-valued higher Green's functions $\Psi^{m,n}(z,w)$, so
compare. (2) Hecke equivariance, $\Xi_{T_pf}=T_p\Xi_f$? (3) Higher-order poles --- $\Xi^\circ_k$
exists only because a simple pole gives $\mathfrak a=2\pi i\,a_zP_z$; at order $P$ the principal
part has $P$ independent numbers and no factorisation, so that is the only place the $W$-position
can carry more than the residue. (4) $\ker\Xi$ beyond $M^!_k$. (5) Identify $x$.

**For Brown.** Their (5.7) runs $H^1_{dR}(U_\Gamma)\to\bigoplus(\mathrm{Sym}^kH^1)^{\Gamma_w}$ with a
splitting that is Hodge-only (Cor.~5.16, Rem.~5.17); their Lem.~5.10 is exactly $\mathfrak a$. $\Xi$
runs the other way into $\ker(\mathrm{Res})\cong W$, so it looks like the Betti counterpart. Precise
question: *is "local residue data $\to W$, by winding an $S$-antisymmetric cycle with an Eisenstein
counterterm", the Betti splitting of (5.7), and is it known?*

**Search terms.** Start with Bruggeman--Choie--Diamantis, *Holomorphic automorphic forms and
cohomology* (Mem.~AMS 253, arXiv:1404.6718) --- already in the bibliography, and the standard
reference for Eichler cohomology WITH singularities. Then: "Eichler cohomology meromorphic modular
forms"; "rational period functions" (Knopp); "Sczech Eisenstein cocycle"; "residue exact sequence
local system modular curve"; "periods of meromorphic Poincaré series" (Bringmann--Kane--von
Pippich); "higher Green's functions period polynomials"; "Manin symbols coefficients residues".

**Code.** `arc_numerics/rf/{zeta17_check,zeta17_highprec,Wpm18_projection,Wpm_direction,Wpm_heegner_exact,theta_laurent,theta_k18_DeltaE4,parity_dims}.py`;
logs `zeta17.txt`, `z17hp.txt`, `wpm18.txt`, `thetalaurent.txt`, `thetak18.txt`, `paritydims.txt`.

---

## Episode 15 — "Why is $\Phi$ not just a lens-pole correction?", and the proof of `thm:arcperiod` gets audited (PRIMARY-SOURCE: session `1a93c382`, 2026-09-23, Opus 5.5)

**DAM's question.** He was "99.999% sure" the $\Phi(f)\,r_{E_k}$ subtraction was needed only
for a pole in the lens, the region between $\gamma^{\rm arc}$ and the horizontal segment
$\rho\to\rho+1$.

**Answer.**

- $\Phi(f)=\int_{\rm arc}f\,d\tau$. Pushing the arc up to $i\infty$ sweeps exactly $\mathcal F$,
  so $\Phi(f)=c_f(0)+2\pi i\sum_{z\in\mathcal F}\operatorname{Res}_zf$. Every pole orbit has a
  representative in $\mathcal F$, so **every** interior pole moves $\Phi$ off $c_f(0)$.
- Lens poles lie *below* the arc and never enter $\Phi$ directly. The lens is what separates the
  arc from the flat segment (`eq:geocross`), and that is the only place it matters.
- Evidence: $f_z=\Delta E_4J'/(J-J(z))$ at $z=(1+\sqrt{-7})/2$, height 1.32, far from the lens
  (whose top is height 1). Here $c_f(0)=0$ but $\Phi=2\pi i\,g(z)$, measured in `cm18.txt`.

**Where the belief came from:** the commented-out block at `.tex` 1629–1632, which states the
lens rule. The live `.tex` at 1170 also has the sign wrong: it says "less $2\pi i$", and it
should be plus.

**Kernel check (requested first).** The $h_\pm$ ambiguity is carried by the principal part
*anywhere*, cusp or interior:
$L^*_+-L^*_-=e^{-i\pi s/2}[C(s)\mathcal D(s)-e^{i\pi(s-1)}C(k-s)\mathcal D(k-s)]$, with
$\mathcal D(s)=\sum c(-m)m^{-s}+2\pi i\sum_{\mathcal F}\operatorname{Res}(f\,{\rm Li}_s(q))$.

- Verified to $2.3\times10^{-30}$; exactly zero at critical $s$.
- Script: `arc_numerics/Lf/kernel_ambiguity_interior.py`; log `logs/kamb.txt`.

**Audit of the $W$-membership proof (DAM: it must be "a (Claude-claimed) proof, rather than some
vibe-sketched piece of shit").** Five gaps:

- **G1.** `lem:arcon` asserts $\hat C_T=0$, which contradicts `lem:CTPhi`.
- **G2.** `eq:arcondefect` has the wrong shape. The true defect carries $|_{S(1+U+U^2)}$ and has
  no $D_E$ term.
- **G3.** `def:arczero` handles $i$ with the spiral, but the spiral passes through $i$. The
  symbol $W$ is also used for two different windings.
- **G4.** `eq:arczero` at $i$ omits $r_S$.
- **G5.** `cor:arcell`'s $c\neq0$ is circular. In fact $c\propto\operatorname{Res}_{\rho+1}\omega_f$,
  which is non-zero for pole order $P\le k-1$ and may vanish above that.

**Fix, derived by hand: `lem:relator`.** For any admissible pair at any $\tau_0$,

$$\hat r_f|_{(1+S)}=(2\pi i)^{n+1}\oint_{\gamma^S+S\gamma^S}\omega_f,\qquad
\hat r_f|_{(1+U+U^2)}=(2\pi i)^{n+1}\oint_{(1+U^{-1}+U^{-2})(S\gamma^T+\gamma^S)}\omega_f.$$

- With $\gamma^S=\gamma^T=c$, a chain with $Sc=-c$ (at $i$, the mean of the two indentations),
  both chains are zero. That is $W$-membership in two lines.
- At $\rho$ the $U$-chain is a small loop about $\rho+1$, and the value is read off at winding 0.
- By-product: $\hat r_f=\int_c f\,\mathcal K\,d\tau$, where $\mathcal K$ is the residue kernel.
  So the side choice at an $S$-pair costs exactly $\pm\Xi_{f;p}$.

**Status: PROPOSED, not applied.** DAM: "propose a detailed rewrite, rather than do the rewrite
yourself. This is delicate." The full proposal, with the structure (main §4 + Appendices A, B, C),
label moves, draft text and seven decisions for DAM, is in
`section4_rewrite_proposal_20260923.md`.

---

## Related but probably-downstream chats (for context, not primary episodes)

- [Thematic summary of DR-B](https://claude.ai/chat/bd99e11d-2039-478e-b04b-54bdfd0e2b43) — 2026-06-05
- [Retry request](https://claude.ai/chat/91f8be0f-ce3b-4717-aca8-6dc0ce08ea17) — 2026-06-06
- [Understanding cohomology classes and period polynomials](https://claude.ai/chat/d8bf2d6e-7a24-4e29-b83e-bbd11c1bb11a) — 2026-06-13

These appear to be later consolidation/exposition passes on the same results
(the $\tau_0=i$ / T-gauge-fixing framing in `dr_b_periods_and_Lfunctions.tex`),
not new discoveries.

---

## Caveat

Episodes 1 and the claude.ai chat candidates (2a, old 2b) are based on search
snippets, not full transcripts — treat their content characterizations as
leads to verify against the actual chats.

**Local log floor:** the oldest local Claude Code session is `841ef028`,
starting 2026-05-18 at 12:25 UTC. Anything before that — including the
integer-$s$ → entire-in-$s$ Hurwitz step — is not in local logs and must
be recovered from claude.ai chat history. Episode 2b (the τ₀ = i and
Hurwitz-unfolding account) is from direct reads of local `.jsonl` logs and
should be treated as a primary-source record.
