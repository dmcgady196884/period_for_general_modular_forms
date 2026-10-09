# §4 rewrite proposal — 2026-09-23

> **RESUME HERE.** Session `1a93c382-1900-4a2b-a605-960f23537aa4`; resume with
> `claude --resume` from the repo root.
>
> State at the close of that session:
>
> - DAM asked for a *proposal*, not edits. **Nothing in the `.tex` has been changed.**
> - The Φ question is answered. Φ sees every pole in 𝓕; lens poles don't enter Φ. The wrong
>   lens rule lives in the commented block at `.tex` 1629–1632. See also Episode 15 in
>   `ai_collab/dr_b_llm_episodes.md`.
> - The kernel interior-pole check is done (`arc_numerics/logs/kamb.txt`).
> - **Waiting on DAM's answers to §8 below** before any `.tex` edit. When editing, show a diff
>   for each block first.
>
> Files from other sessions dated 2026-09-25 (`principal_part_ambiguity_20260925.tex`,
> `appendix_branch_choice_20260925.tex`, `hurwitz_kernel_handoff_20260925.md`,
> `verify_kernel_20260925.py`) postdate this proposal. Reconcile Appendix A (§6) with them before
> drafting it.

A proposal only. No `.tex` has been edited. Line numbers refer to
`finite_contour_cocycles_short.tex` as of this session: §4 runs 750–2331, and the bibliography
starts at 2332.

Contents:

1. Audit of `thm:arcperiod` as it stands: five gaps.
2. A replacement proof that closes all five and is shorter.
3. Target structure: main text §4, then Appendices A, B and C.
4. Where each label goes.
5. Draft text for the new statements, the choices subsection and the Brown–Fonseca paragraph.
6. Appendix A: statements, each with its provability status.
7. Errors found in passing.
8. Decisions for you.

---

## 1. Audit: the current proof of `thm:arcperiod` has five gaps

The theorem is true: every measurement agrees with it. Five steps of the written proof do not
hold up.

**G1. `lem:arcon`'s proof asserts $\hat C_T=0$ (2149–2150).**
The proof says: "$\Phi$ is by construction the quantity $\hat C_T$ measures, so
Definition~\ref{def:rfhat} gives $\hat C_T=0$".

This is false. $\hat C_T$ is a cocycle component, and Definition `def:rfhat` does not touch it.
Lemma `lem:CTPhi` says $\hat C_T=-(2\pi i)^{n+1}A_0\Phi(f)$, which is non-zero whenever
$\Phi(f)\neq0$.

The next displayed formula, $\hat r_f|_{(1+U+U^2)}=\hat C_{U^3}-\Phi(f)D_E$, inherits the error.
The correct statement is $\hat r_f|_{(1+U+U^2)}=\hat C_{U^3}$, with no $D_E$ term (see §2 below).

**G2. `eq:arcondefect` has the wrong shape.**
The left side, $\hat r_f|_{(1+U+U^2)}$, is in the image of $|_{(1+U+U^2)}$. The right side is a
bare sum $\sum\nu_ja_{p_j}(X-p_jY)^n$ plus a $D_E$ term. The correct defect, from §2, is

$$(2\pi i)^{n+1}\Big(\oint_{\hat\gamma^T-\hat\gamma^S}\omega_f\Big)\Big|_{S(1+U+U^2)},$$

which carries $|_{S(1+U+U^2)}$ and has no $D_E$ term. The Vandermonde non-vanishing argument
rests on the wrong formula. None of this is on the theorem's critical path, since case 2 uses
only $m=0$.

**G3. `def:arczero` puts the pole at $i$ in the shifted class, where the spiral passes through
$i$.**
`eq:arcspiral` has $\log r=0$ at $\theta=\pi/2$, so the spiral runs through $i$ exactly, and a
pole at $i$ sits on the contour.

Two different windings are also being used under one name:

- $W(f,\epsilon)$ of `def:arcwind` is the $U^3$-loop's winding about $\rho+1$.
- The "$W$ odd at $i$" in the paragraph after `def:arczero` is the $S^2$-loop's winding about $i$.

**G4. `eq:arczero` at $i$ is missing $r_S$.**
At $i$ both segments pass through the pole ($\hat\gamma^S=\hat\gamma^T$), so switching
indentation moves both. The correction should carry
$\operatorname{Res}_i\big(f\,(\tau^{s-1}+\tilde k)\big)$, as in `eq:sidediff`. As written it
carries only $\operatorname{Res}(f\tilde k)$.

At $\rho$ the formula is right: only the connector $\beta$, which is part of $\hat\gamma^T$
alone, winds.

**G5. `cor:arcell`'s uniqueness claim is circular, and may be false.**

- $c\neq0$ is justified by "the defect is non-zero at $W=-1$", which is the claim itself.
- Correct statement: $c$ is a non-zero multiple of $R:=\operatorname{Res}_{\rho+1}\omega_f$. $U$
  fixes $\rho+1$, so $R|_U=R$ and $R|_{(1+U+U^2)}=3R$.
- $R\neq0$ holds when the pole order satisfies $P\le k-1$: the leading Laurent coefficient
  appears in $\operatorname{Res}\big(f\,(\tau-\rho-1)^{P-1}\big)$, and $P-1\le n$.
- For $P>k-1$ with vanishing low Laurent coefficients, $R=0$ is possible. Then every class
  satisfies both relations, and uniqueness fails.
- The at-$i$ half of the corollary is one sentence ("the same argument ... gives $W$ odd").

Uniqueness is not needed for $W$-membership.

**Structural problems:**

- `lem:CTPhi` is a lemma environment nested inside the proof of `lem:arcoff` (2084–2100 sits
  inside 2054–2122).
- `lem:arcoff` writes $\hat C_{S^2}=\hat C_{-I}=\hat C_{U^3}$. For meromorphic $f$ these are
  integrals over different closed paths, so equating them as "$\hat C_{-I}$" assumes what is
  being proved. The loop integrals must be kept separate, which §2 does.
- `lem:georetrace` and `def:arcwind` say $(ST)^3$; everything else says $U=TS$, $U^3$.

---

## 2. Replacement proof: one identity, valid in every class

The new ingredient is a single identity. It needs no retracing, no fixed point, and no
particular class.

**Proposed Lemma `lem:relator`.**
Let $f\in F_k$ and $\tau_0\in\HH$, and let $(\hat\gamma^S,\hat\gamma^T)$ be admissible segments
(a pair of 1-chains is allowed) avoiding $f^{-1}(\infty)$. Put

$$\ell_S:=\hat\gamma^S+S\hat\gamma^S,\qquad
\ell_U:=(1+U^{-1}+U^{-2})\big(S\hat\gamma^T+\hat\gamma^S\big).$$

Both are closed 1-chains in $\HH\setminus f^{-1}(\infty)$. Then

$$\hat r_f\big|_{(1+S)}=(2\pi i)^{n+1}\oint_{\ell_S}\omega_f,\qquad
\hat r_f\big|_{(1+U+U^2)}=(2\pi i)^{n+1}\oint_{\ell_U}\omega_f .$$

**Proof.** I derived this by hand. It uses only facts already in the draft.

1. By `eq:rfsplit`, $r_f=C_S-(E|_S-E)$ with $C_S=(2\pi i)^{n+1}\int_{\hat\gamma^S}\omega_f$.
   By `lem:CTPhi`, $C_T-(E|_T-E)=-(2\pi i)^{n+1}A_0\Phi(f)$. Neither proof uses $c_f(0)=0$ or
   anything about the paths beyond their endpoints.
2. The substitution in `lem:resequiv` gives $(\int_c\omega_f)|_g=\int_{g^{-1}c}\omega_f$.
   Hence $C_S|_S=(2\pi i)^{n+1}\int_{S\hat\gamma^S}\omega_f$ and
   $C_T|_S=(2\pi i)^{n+1}\int_{S\hat\gamma^T}\omega_f$, since $-I$ acts trivially on $\HH$.
3. **S-relation.** $r_f|_{(1+S)}=C_S|_{(1+S)}-(E|_{S^2}-E)=C_S|_{(1+S)}$, because $S^2=-I$ acts
   trivially on $V_n$ ($n$ even). This equals $(2\pi i)^{n+1}\oint_{\ell_S}\omega_f$.
4. **U-relation.** Put $Q:=(C_T-(E|_T-E))|_S+r_f=(C_T|_S+C_S)-(E|_U-E)$. Then
   $Q|_{(1+U+U^2)}=(C_T|_S+C_S)|_{(1+U+U^2)}-(E|_{U^3}-E)=(2\pi i)^{n+1}\oint_{\ell_U}\omega_f$.
   The other expression for $Q$ gives
   $r_f|_{(1+U+U^2)}=(2\pi i)^{n+1}\oint_{\ell_U}\omega_f+\Phi(f)\,D$, with
   $D:=(2\pi i)^{n+1}A_0|_{S(1+U+U^2)}$.
5. **Apply to $E_k$.** It is holomorphic on $\HH$, so both loop integrals vanish, and
   $\Phi(E_k)=1$. Hence $r_{E_k}|_{(1+S)}=0$ and $r_{E_k}|_{(1+U+U^2)}=D$. Subtracting
   $\Phi(f)$ times these gives the lemma. ∎

This single lemma replaces `lem:geoloops`, `lem:arcoff`, the relation parts of `lem:arcon` and
`lem:arcwind`, and all of `eq:arcoffS`, `eq:arcoffU`, `eq:arcoffprop` and `eq:arcondefect`.
It also makes $D=r_{E_k}|_{(1+U+U^2)}$ come out of the argument rather than being identified
afterwards.

**Proposed Lemma `lem:georetrace`, replacing the current one.** Let $c$ be a 1-chain in
$\HH\setminus f^{-1}(\infty)$ from $\rho$ to $\rho+1$ with $Sc=-c$. With $\tau_0=\rho+1$ and
$\hat\gamma^S=\hat\gamma^T=c$, both $\ell_S$ and $\ell_U$ are the zero chain.

*Proof.* $\ell_S=c+Sc=0$, and $S\hat\gamma^T+\hat\gamma^S=Sc+c=0$. ∎

That is the whole proof: two chain identities, with no pictures and no "retraces".

**Proposed Lemma `lem:symchain`.** If $f\in F_k$ has no pole on $\SLZ\cdot\rho$, such a chain
exists. It equals $\gamma^{\rm arc}$ outside $\epsilon$-discs about the poles of $f$ on the arc.

*Construction.*

- **An $S$-pair $\{p,Sp\}$ on the arc with $p\neq i$.** Replace the arc near $p$ by a
  semicircle $\delta_p$, on either side. Replace it near $Sp$ by $-S\delta_p$.
- **The point $i$.** Replace the arc near $i$ by $\tfrac12(\delta_+-S\delta_+)$, where
  $\delta_+$ is the outside semicircle. Then $-S\delta_+$ is the inside one.
  Check: $S(\delta_+-S\delta_+)=S\delta_+-\delta_+$, using $S^2=-I$ on $\HH$.
- **The rest.** $S$ maps the arc to itself with the orientation reversed, which gives $Sc=-c$
  on the unindented part. ∎

**Proposed proof of `thm:arcperiod`.**

- **No pole on $\SLZ\cdot\rho$.** Take $c$ from `lem:symchain` as both segments at
  $\tau_0=\rho+1$. By `lem:georetrace` both relator chains vanish, and `lem:relator` gives
  $\hat r_f\in W$.
- **A pole on $\SLZ\cdot\rho$.** Work at $\tau_0(\epsilon)$ of `def:arcshift`:
  - $\hat\gamma^S=\sigma$, the spiral. If $f$ also has a pole at $i$, replace $\sigma$ near $i$
    by the mean chain as in `lem:symchain`.
  - $\hat\gamma^T=\beta\cdot\sigma$.
  - $S\sigma=-\sigma$ gives $\ell_S=0$ and $\ell_U=(1+U^{-1}+U^{-2})S\beta$.
  - $S\beta$ lies in a disc about $\rho+1$ that holds no other pole once $\epsilon$ is small.
    So `lem:relator` gives $\hat r_f|_{(1+U+U^2)}=(2\pi i)^{n+2}N\,R$, where $N$ is the winding
    of `def:arcwind` and $R=\operatorname{Res}_{\rho+1}\omega_f$.
  - By `lem:wall_arc` and `lem:windclass`, $\hat r_f$ is affine in $N$. So its value at $N=0$,
    as `def:arczero` prescribes, satisfies both relations. ∎

What the theorem then depends on:

- from §3: `def:segments`, `def:Lint`, `def:kernel`, `def:rf`, `prop:periodint`,
  `eq:rfsplit`, `lem:resequiv`, `def:w`;
- from the §4 preamble: `def:georef`, `def:arcshift`, `def:arcwind`, `lem:windclass`,
  `def:arczero` (restricted to $\rho$), `lem:wall_arc`, `lem:tau0indep` (only so that
  `def:arczero` is well defined);
- from §4.2: `def:Phi`, `def:rfhat`, `lem:CTPhi`, `lem:relator`, `lem:georetrace`,
  `lem:symchain`.

**Two by-products that fall out for free:**

- **(a)** `lem:arcfunceq` holds under the same hypothesis as $W$-membership. Its $S$-segment
  step uses exactly $S\hat\gamma^S=-\hat\gamma^S$. Worth one sentence: the functional equation
  and $W$-membership rest on the same property of the contour.
- **(b)** With $\hat\gamma^S=\hat\gamma^T=c$, `prop:periodint` and `def:Phi` give

  $$\hat r_f=\int_c f(\tau)\,\mathcal K(\tau;X,Y)\,d\tau,$$

  with $\mathcal K$ the residue kernel of `def:reskernel`. So **the period polynomial is the
  residue-polynomial kernel integrated along the arc.**
  - Two admissible chains that differ at one $S$-pair differ by $\pm(c_p-Sc_p)$, so the side
    choice moves $\hat r_f$ by $\pm\Xi_{f;p}$.
  - This replaces `cor:symspan` and `eq:symshift`, and ties §4.1 to §4.2 in one line.

**Before any of this goes in:** recheck `lem:windclass` in the §3 path convention. In that
convention the path for $C_g$ runs from $g^{-1}\tau_0$ to $\tau_0$, whereas `lem:windclass` uses
$\gamma_g$ from $\tau_0$ to $g\tau_0$. The residue class of the winding may come out as $+1$
rather than $-1\pmod3$.

This does not affect the proof, which needs only two distinct realised values. It does affect the
weights $\tfrac23,\tfrac13$ and `eq:arczero`'s explicit form. It is a two-line check.

---

## 3. Target structure

Main text:

- **§4 preamble**: $L^*$ in a class, the base point, the geodesic class, the prescription at
  $\rho$, the wall formula, $\Phi$, and the functional equation.
- **§4.1 Residue polynomials**: current §4.3, essentially unchanged.
- **§4.2 Period polynomials in the geodesic class**: $\hat r_f$, `lem:CTPhi`, `lem:relator`,
  `lem:georetrace`, `lem:symchain`, `thm:arcperiod`, and $\hat r_f=\int_c f\mathcal K$.
- **§4.3 Choices**: short prose covering the three ambiguities, then the Brown–Fonseca
  paragraph.

Appendices:

- **A. The kernel ambiguity** (new).
- **B. Explicit formulae on the arc**: old §4.1, minus `def:Phi` and `lem:arcfunceq`, plus
  `prop:arcside`.
- **C. Relation to the L-functions of [1806]**: old §4.2.

Material off the theorem's critical path: I recommend parking it below `\end{document}`, as you
did in §4 before, rather than creating an Appendix D. `lem:relator` supersedes it and parts of it
are wrong (G1, G2, G5). This covers:

- `lem:geoloops`, `lem:arcoff` (the text; the content is absorbed);
- `lem:arcon`, `eq:arcondefect`;
- `lem:arcdeform`;
- `cor:symspan`, `eq:symshift`, `eq:symspan`, and the prose at 2184–2193 and 2227–2241;
- `cor:arcell`, `lem:arcwind`;
- the two prose paragraphs at 2281–2285 and 2314–2330, whose content moves to §4.3.

Only `sec:arc` is cited from outside §4 (at line 165, in the intro), so none of the moves breaks
an outside reference.

---

## 4. Where each label goes

**Stays in the §4 preamble:**

- `prop:indepHomotopy_arc`, `def:class_arc`, `eq:Lclass_arc`. Cut the vague parenthetical at
  800–802 ("subject to natural constraints ... monodromies"); `lem:tau0indep` says it
  precisely.
- `lem:tau0indep`, `eq:ktilshift`.
- `def:georef`. Drop its bracket title. The two paragraphs after it (848–856) stay.
  - In 858–865, cut the last sentence, "All three cases are treated the same way, by moving
    the poles off the contour rather than moving the contour". It says the opposite of 867:
    "We move the contour off the poles, never the poles off the contour."
- `lem:georetrace`, replaced by the chain version from §2. Moving it to §4.2 is also fine,
  since `def:respoly`'s prose at 1686 cites it.
- `def:arcshift`, `eq:arcspiral`, and the paragraph beginning "Three features" (904–911).
  Delete that paragraph's "...and Lemma~\ref{lem:georetrace} applies to the $S^2$-path"; the
  new `lem:georetrace` doesn't need it.
- `def:arcwind`: rename the winding, e.g. to $N$, to end the clash with the subspace $W$, and
  write $U^3$ for $(ST)^3$.
- `lem:windclass` with `eq:relatorcollapse`. Recheck the convention (§2).
- `def:arczero`, restricted to $\SLZ\cdot\rho$; $i$ is handled by `lem:symchain`'s mean chain.
  - In the paragraph with `eq:arczero`, fix G4 or restrict it to $\rho$.
  - Cut "and by Section~\ref{sec:arcperiod} either separation puts $\hat r_f$ in $W$"
    (980–981). It is a forward reference, and §4.3 now says it properly.
- `eq:Lstararc`, `lem:wall_arc`.
- **Moved up from old §4.1:**
  - `def:Phi`, with its prose fixed (§7, E1 and E2). Optionally add one line: since
    $L^*_S$ is entire and $\operatorname{Res}_{s=0}\tilde k=-1$, we get
    $\Phi(f)=-\operatorname{Res}_{s=0}L^*(f,s)$. That gives $\Phi$ a meaning without reference
    to the contour. This is currently buried in `rem:arcresidue`.
  - `lem:arcfunceq` with `eq:ktilflip`. Note its hypothesis: $S\hat\gamma^S=-\hat\gamma^S$ on
    the chain.

**Moves to Appendix B:**

- The paragraph at 1050–1062, "This provides the basics ... polylogarithmic projections",
  which becomes B's opening.
- `prop:arcside` with `eq:sidediff`, and the exclusion paragraph at 1046–1048.
- `def:Finfty_arc`, `def:proj_arc`, `def:geoproj`, `lem:proj_arc`, `def:cRf_arc`,
  `lem:projspace_arc`, `lem:ellLaurent`, `def:geoIN`, `lem:geopolylog` with `eq:geoG0` and
  `eq:geocross`, `lem:geoedge`, `lem:geowhole`, `lem:geotail`, `thm:geomero`.
- In `lem:ellLaurent`'s lead-in (1173–1175), cut "and the period polynomials constructed later".
  Period polynomials no longer use it.
- `thm:geomero` cites `lem:arcfunceq`, which would now be in the main text. That is a backward
  reference from the appendix, which is fine.

**Moves to Appendix C:** 1548–1618 as they stand, with three changes.

- **`prop:arc1806`'s proof does not prove its claim.** "$L^{*\uparrow}$ is the L-function
  of [1806]" is identified only by matching the residue at $s=0$, and two meromorphic functions
  are not identified by one residue. The numerics say the claim is true (1e-25, see
  `raised_contour_equals_1806`).
  - Either prove it by expanding term by term above all poles, against [1806]'s Thm 1 (which
    survives the LemRegX3 bug);
  - or state it as a numerically verified comparison, without a proposition.
- **The remarks.** `rem:arcresidue` and `rem:arcextremal` are remark blocks; convert them to
  prose. Also, `rem:arcresidue` says "the quantity that must be subtracted ... is the same object
  as the departure from the residue convention". The subtracted quantity is $\Phi(f)$; the
  departure is $\Phi(f)-c_f(0)$. Reword.
- **The paragraph "One might then ask why the arc ..." (1610–1616).** With `lem:relator` it
  becomes a clean statement. In the raised class $\ell_U$ encloses the poles between the arc and
  the raised segment, so the defect is those residues, carried by $|_{S(1+U+U^2)}$.

**§4.1 (residue polynomials), minimal changes:**

- Drop the bracket titles of `def:respoly` and `def:respolyfull`.
- Delete the commented block at 1621–1651. It contains the false lens claim you were
  remembering (1629–1632): "That subtraction is $c_f(0)\,r_{E_k}$ unless a pole lies between
  $\gamma^{\rm arc}$ and the horizontal line from $\rho\to\rho+1$". It is harmless in the PDF
  but it is where the wrong belief lives.
- The prose at 1686–1688 cites `lem:georetrace`. Under the chain version the sentence still
  reads correctly.

---

## 5. Draft text

These drafts are for you to cut down. There is no editorializing in them.

### 5a. After `def:rfhat` (§4.2)

> With $\hat\gamma^S=\hat\gamma^T=c$, Proposition~\ref{prop:periodint} and
> Definition~\ref{def:Phi} give
> \begin{equation}\label{eq:rfhatK}
> \hat r_f=\int_c f(\tau)\,\mathcal K(\tau;X,Y)\,d\tau ,
> \end{equation}
> with $\mathcal K$ as in Definition~\ref{def:reskernel}.

### 5b. §4.3 Choices (prose, no environments except one short lemma if you want the $\Xi$ statement formal)

> Three choices enter the construction.
>
> *The contour.* $L^*(f,s)$ depends on the homotopy class of the segments
> (Lemma~\ref{lem:wall_arc}), and outside the class of Definition~\ref{def:georef} the
> relator chains of Lemma~\ref{lem:relator} enclose poles and $\hat r_f$ leaves $W$. Within it,
> at $i$ and on $\SLZ\cdot\rho$ nothing is free. At each $S$-pair $\{p,Sp\}$ of poles on
> $\gamma^{\rm arc}$ with $p\ne i$ one chooses which member is passed outside. The two choices
> differ by $\pm(c_p-Sc_p)$, so by \eqref{eq:rfhatK} they change $\hat r_f$ by $\pm\Xi_{f;p}$,
> which lies in $W$ by Proposition~\ref{prop:respolyW}.
>
> *The reference form.* $E_k$ in Definition~\ref{def:rfhat} may be replaced by any $g$ with
> $\Phi(g)=1$; two choices shift $\hat r_f$ by $\Phi(f)\,r_{g-g'}$. For $g\in M_k$ the shifts
> run over $r(S_k)$, of dimension $\dim S_k$; for $g\in M^!_k$ over $r(S^!_k)$, of dimension
> $2\dim S_k$~\cite{BGKO}. [family-selected reference: the $hJ'/(-2\pi i\,c_h(1))$ sentence
> from 2326–2330 goes here]
>
> *The kernel.* $\tilde k$ is one of two one-sided inverses of $1-T$; the other changes
> $L^*(f,s)$ at non-integer $s$ by a sum over the principal parts of $f$, and leaves
> $L^*(f,\ell+1)$, $1\le\ell+1\le k-1$, unchanged, hence $\hat r_f$ and $\Xi_{f;z}$ as well
> (Appendix~\ref{app:kernel}).

The $\dim S_k$ claim for $g\in M_k$ is Eichler–Shimura injectivity on $S_k$. For the
$2\dim S_k$ claim, check that BGKO's theorem is literally $r(S^!_k)$ of that dimension before
citing it; the bib entry exists, is not yet cited anywhere, and I have not re-read the paper.

### 5c. Brown–Fonseca, closing §4.3

> Brown and Fonseca~\cite{BrownFonseca2025} also attach period polynomials in $W$ to meromorphic
> modular forms, and resolve the same ambiguity differently. Their reference forms are
> Petersson's Poincaré series $\Psi^{p,q}_\Gamma(z,w)$, meromorphic in $z$ but real-analytic in
> the pole position $w$; subtracting them gives a canonical splitting that does not depend
> holomorphically on the poles. The construction here keeps every quantity holomorphic in the
> positions of the poles of $f$, and pays with the choices above.

Check the $\Psi^{p,q}_\Gamma$ notation and wording against arXiv:2508.04844 before this goes in.
You might also want one sentence saying the two constructions differ even where both apply:
at $\dim S_k=0$, BF gives $\tilde r_f\equiv0$ and the arc gives a non-zero element of $W$
(`geodesic_contour_findings`).

---

## 6. Appendix A — kernel ambiguity: statements and status

Notation: $h_+[w](\tau):=\zeta(-w,\tau+1)$ and $h_-[w](\tau):=-e^{i\pi w}\zeta(-w,-\tau)$, with
$\tilde k_\pm:=h_\pm[s-1]-e^{i\pi(s-1)}h_\pm[k-1-s]$, so that $\tilde k_+=\tilde k$ of
`def:kernel`. Also $C(a):=(-2\pi i)^{1-a}/\Gamma(1-a)$.

- **A.1 Both kernels solve the shift identity.** $h_\pm[w](\tau)-h_\pm[w](\tau-1)=-\tau^w$.
  - Status: **provable in one line.** For $h_-$ use $(-\tau)^w=e^{-i\pi w}\tau^w$ on $\HH$.
- **A.2 Their difference is a Lipschitz sum.** $h_+[w]-h_-[w]=C(w+1)\,{\rm Li}_{w+1}(q)$, hence
  $\tilde k_+-\tilde k_-=C(s){\rm Li}_s(q)-e^{i\pi(s-1)}C(k-s){\rm Li}_{k-s}(q)$.
  - Status: **classical** (Hurwitz/Lipschitz); verified to 4e-26 (`kamb.txt`).
- **A.3 Characterization.** $h_+$ is the unique solution of A.1 that is holomorphic on
  $\CC\setminus(-\infty,-1]$ with polynomial growth as $|{\rm Im}\,\tau|\to\infty$, up to an
  additive constant. $h_-$ is the mirror image, on $\CC\setminus[0,\infty)$.
  - Status: **provable.**
  - Proof: the difference of two solutions is 1-periodic. Its boundary values across the cut
    agree after a shift by an integer, so it is entire. It is then a Laurent series in $q$,
    and growth in both vertical directions leaves only the constant.
  - This is from the other session's review (`kernel_modular_origin_review`); I have re-derived
    it.
- **A.4 The formula for the difference.**
  $L^*_+-L^*_-=e^{-i\pi s/2}\big[C(s)\mathcal D_f(s)-e^{i\pi(s-1)}C(k-s)\mathcal D_f(k-s)\big]$,
  where

  $$\mathcal D_f(s)=\sum_{m\ge1}c_f(-m)m^{-s}+2\pi i\sum_{z\in\mathcal F}\operatorname{Res}_{\tau=z}\big(f\,{\rm Li}_s(q)\big).$$

  - Status: **provable** (push $\gamma^T$ to $i\infty$; the $S$-segment is common to both
    kernels).
  - **Verified today to ≤2.3e-30:** interior pole alone, cusp principal part alone, and their
    sum.
  - The text must say the difference is carried by the principal part **anywhere, cusp or
    interior**.
- **A.5 Critical values agree.** $L^*_+=L^*_-$ at $s=1,\dots,k-1$, because $1/\Gamma(1-s)$ and
  $1/\Gamma(1-(k-s))$ both vanish there. Hence $r_f$, $\hat r_f$ and $\Xi_{f;z}$ do not depend
  on the kernel.
  - Status: **immediate** from A.4; exactly 0 at $s=6$ in `kamb.txt`.
- **A.6 $L^*_+$ is BFK's L-function.** It equals the Bringmann–Fricke–Kent L-function with
  principal branches, the convention stated in DLRR, arXiv:2107.12366, eq. (3.3).
  - Status: **verified numerically** to 1e-23 on $\Delta'$ (`hurwitz_kernel_handoff`), not
    proved. Either state it as a check, or prove it by expanding BFK's series.
- **A.7 Reality.** For $f$ with real Fourier coefficients and real $s$,
  $L^*_-=\overline{L^*_+}$, so $L^*_+$ is not real once $f$ has a principal part, and
  $\tfrac12(\tilde k_++\tilde k_-)$ gives ${\rm Re}\,L^*_+$.
  - Status: **verified numerically**; the proof idea is the reflection $\tau\mapsto-\bar\tau$.
  - The handoff flags: check whether BFK or DLRR already remark on the non-reality before
    claiming it.
- **A.8 Growth.** $|\tilde k_\pm|$ grows polynomially, of degree $\max({\rm Re}\,s,k-{\rm Re}\,s)$
  (times $e^{-\pi{\rm Im}\,s/2}$).
  - Status: **provable** from Hurwitz asymptotics; include it only if something uses it.
- **Open, prose only:** the whole family $\lambda\tilde k_++(1-\lambda)\tilde k_-$ plus constants
  satisfies A.1; only $\lambda\in\{0,1\}$ has a slit-plane domain. Whether Hecke equivariance
  selects $\lambda$ is not known.

---

## 7. Errors found in passing

- **E1 (1170) — sign.** It should read "$c_f(0)$ **plus** $2\pi i$ times the residues". The arc,
  traversed left to right, minus a high segment bounds the region counterclockwise. `cm18.txt`
  measures $\Phi(f_z)=+2\pi i\,g(z)$.
- **E2 (1168–1169) — the condition.** "whenever no pole of $f$ lies above
  $\hat\gamma^T_f$", said of the arc, is exactly "$f$ has no pole in $\HH$ off the contour". The
  region swept is $\mathcal F$, which meets every orbit.
  - Suggested wording: $\Phi(f)=c_f(0)+2\pi i\sum_z\operatorname{Res}_{\tau=z}f$, the sum over
    the poles of $f$ in $\mathcal F$ off $\gamma^{\rm arc}$, one per orbit; in particular
    $\Phi=c_f(0)$ on $M^!_k$.
  - Also, "$\Phi(f)=0$ for $f\in S_k$" can be widened to $S^!_k$.
- **E3 (864–865 vs 867).** The two sentences contradict each other (§4 above).
- **E4 (2039–2042).** "The remaining cases are reached by continuity from it" is not how they
  are reached; they are reached by the prescriptions. Cut it.
- **E5 (2078–2079)**, "This is where the hypothesis ... it is the only place". This is
  commentary; it goes with `lem:arcoff`.
- **Pending audits, unchanged:** the `onarcres.txt` discrepancy, and T8–T24 of the review. The
  `prop:arcside` side-flip audit is partly overtaken: under `lem:symchain` the only side choice
  left is per $S$-pair, and it costs $\pm\Xi_{f;p}$.

---

## 8. Decisions for you

1. Adopt `lem:relator` + chain `lem:georetrace` + `lem:symchain`, replacing the five-lemma
   route? (Recommended: it closes G1–G5 and is about a page shorter.)
2. Drop `cor:arcell`'s uniqueness from the paper, or keep it with the hypothesis
   $P\le k-1$ at $\rho$?
3. Main-text order: residue polynomials, then period polynomials, as now? The new
   \eqref{eq:rfhatK} reads well in that order.
4. Park the superseded material below `\end{document}`, or delete it outright?
5. Old §4.2 as its own Appendix C, or as the last part of B?
6. `prop:arc1806`: prove it or demote it (§4 above)?
7. Put $\Phi(f)=-\operatorname{Res}_{s=0}L^*(f,s)$ into the main text?
