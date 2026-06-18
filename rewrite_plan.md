# Rewrite plan — `finite_contour_cocycles.tex`

Working TODO for the holistic cleanup. Scope: §§2–4 now; §5 deferred until DAM says.
Calibration = `previous_paper_examples_for_style/{20190115lfn,mfdhp_s2}.tex` (DAM's own).

## Done (committed, branch `intro_fable_vibe`)
- Path #2 intro rewrite (L0–L5): primitive/branch framing, finite τ₀-cocycles explicit, τ₀ = anchor.
- Gauge scrubbed; duplicate `eq:bfk-eq-statement` + `eq:Lambda-DR-general-restated` refs fixed; dead-comment + 790-line scratch + editing-changelog swept.
- `thm:periods` parity VERIFIED vs Brown (ℓ even→ω⁺, ℓ odd→ω⁻).
- Retitle + `git mv` → `finite_contour_cocycles.{tex,pdf}`.
- **Bernoulli nuked from the intro** (3 spots → Hurwitz-primary).
- **Bernoulli de-named across §3–§4** (0f6f635): titles renamed after what they do; machinery (`\widehat{k}_T`, `lem:T-piece-equiv`) KEPT per DAM ("had to be there for the generalization"); Bernoulli now only where it computes. → item 4 DONE.
- **De-paragraph §§2–4 + appendix** (f541601, e4ac3da, 46fe7b0, a825187): all 28 body+proof `\paragraph` headers dissolved into flowing prose (§2 narrative; §3 + §4 + App proofs flowed from `Step N` staccato; H^1 now displayed `eq:H1-def`). Only §5's 2 `\paragraph{Remark}` remain (deferred). → item 1 DONE (ex-§5).
- **Thm 3.9 (`thm:periods3`) realigned to Thm 1.4** (46fe7b0): ℓ-indexing, bloat/unref'd display removed. Overhang scan: bfk pair aligned; §5 restatements deferred.

## Holistic agenda (DAM's directives + found issues)

1. **De-`\paragraph{}` → flowing prose.** Kill `\paragraph{}` fragmentation across §§2–5; sections/subsections fine, but prose should flow one to the next (not staccato). Smooth "Step 1:…/Step 6:…" proof scaffolding where natural. DAM: "Paragraphs are more fragmented, and jarring. Not my style."

2. **Numerics → appendix, passing mention only.** Concrete numerical statements live in App `sec:numerics`; body only points at them. Word as *evidence*, NEVER *verification* (we are proving). Partially done.

3. **Tight statements.** Theorem/lemma STATEMENTS short + tight: setup sentence → display → one-line "where". Proofs CAN be hairy (but pretty when possible). Template = `20190115lfn` Thm 1. Offenders: bloated multi-claim titles, e.g. `lem:T-piece-equiv` = "Bernoulli-form vs. Hurwitz-form T-kernel: polynomial identity and integral equivalence" → split/trim. NB: draft uses `\begin{theorem}/\begin{lemma}`; DAM's style uses short `\begin{thm}/\begin{lem}` (labels informal, e.g. `LemReg2`).

4. **Bernoulli = passing breath.** Only the integer-s bridge to the Hurwitz T-kernel; "vital as an intermediate result… doesn't deserve time in the sun longer than a passing breath." ~39 usages remain in §§3–5 + numerics. Big §3 restructure: collapse the parallel Bernoulli machinery — `def:khat-T` (Bernoulli T-segment kernel), `lem:Bernoulli-inv`, `lem:T-piece-equiv` (Bernoulli↔Hurwitz equivalence) — so the Hurwitz kernel is the primary object and Bernoulli appears once, as ζ(−m,a)=−B_{m+1}(a)/(m+1).

5. **Every name-drop EARNED.** "Brown cocycle"/"DR-B"/"de Rham–Betti" are NOT in the draft (good). Audit remaining "Brown" usages (basis/normalization/values) for earned-ness. Per `context.tex`: also cite Bringmann–Kane–Kohnen 2015 + Löbrich–Schwagenscheidt 2020 (closest, currently absent); position the wall-crossing residue vs the LS local polynomial.

6. **Maybe merge §2 + §3.** (DAM: "Not sure.")

7. **Indexing reconciliation (found).** Intro `thm:periods` is ℓ-indexed; body `sec:periods` still uses old s-value indexing → mismatch (e.g. intro τ^{ℓ−1} at line ~465 vs thm τ^ℓ). Reconcile when reworking §3 proofs, or the restated theorem mismatches.

## §5 — DEFERRED (extensive rewrite, only when DAM says)
- Polylog subtraction of ALL pole-heights at/above Im(τ)=1 — required by (a) incomplete-gamma convergence and (b) q-series NON-commutativity at the endpoints τ=i−1, i when poles have Im ≥ 1. This subtlety contaminates 1806 Lemma 2.1 (innocuously) — note it.
- The decomposition f = (f−g) + g (g carries f's cusp poles, regular in the fundamental domain) — needs a careful/rigorous revisit.
- `cor:1806-coincide` lives here (currently the lone undefined ref).

## Open housekeeping
- `cor:1806-coincide` undefined (resolve in §5 rewrite).
