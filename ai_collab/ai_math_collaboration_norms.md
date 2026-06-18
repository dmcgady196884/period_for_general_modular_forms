# AI/LLM collaboration in math research — current landscape (as of June 2026)

## 1. The general publisher consensus (Elsevier, Springer, Wiley, T&F, SAGE; COPE-aligned)

Three points, repeated almost verbatim across publishers:

- **No AI authorship.** An LLM cannot be a co-author — it cannot take legal/scientific
  responsibility, consent, or approve a final manuscript.
- **Full human responsibility.** Authors are responsible for *everything* in the paper,
  regardless of how it was produced. This is the load-bearing clause — disclosure doesn't
  shift responsibility, it just makes the record honest.
- **Disclosure, not prohibition.** Use of generative AI/LLMs should be declared — typically
  tool name + version + purpose — in acknowledgments or a methods/supplementary statement.
  "AI-assisted copy editing" (grammar/wording polish on human-written text) is usually
  *exempt* from disclosure.

## 2. Pure-math-specific data points

- **Annals of Mathematics** has an explicit AI & LLM policy: won't consider papers
  *generated* by AI products, but separately requires authors to declare, as supplemental
  info, the **extent and purpose** of any AI/LLM use. Authorship is restricted to those who
  can take full responsibility for the content. This is the most direct top-tier pure math
  precedent found.
- **AMS** formed an "Advisory Group on AI and the Mathematical Community" (2023); its 2024
  white paper is mostly education-focused. No AMS-wide research-disclosure policy found as
  of this search — MSP/Springer/Elsevier number theory journals (Alg. & Number Theory,
  Research in Number Theory, J. Number Theory) appear to inherit their publishers' general
  templates rather than having bespoke math policies.
- **arXiv** (Dec 2025 / early 2026): not a ban on AI use. A "one-strike" policy — a
  one-year submission ban if a paper contains "incontrovertible evidence" the authors did
  not check AI output (fabricated refs, etc.). The stated principle: using AI is fine,
  *not verifying it* is the violation.

## 3. Precedent for *substantive* (not just drafting) contributions

- **Tao + DeepMind, "Mathematical exploration and discovery at scale"** (Nov 2025): AlphaEvolve
  used as an active research tool on open problems; reported as a transparent experimental
  study, with prompts/outputs published as part of the paper's *subject matter* (not as a
  generic provenance requirement).
- **Don Knuth** (March 2026): credited Claude (and a multi-agent Claude/GPT setup) with
  solving an open Hamiltonian-decomposition problem for TAOCP; wrote up the contribution
  narratively, coined "Claude-like decompositions," and had results independently/formally
  verified (one variant via Lean 4).
- **VUB "vibe-proving" paper** (March 2026): ChatGPT-5.2 produced an original proof of a
  2024 conjecture; framed explicitly around "the bottleneck is now human verification."

## 4. Bottom line for this project

- The norm is **disclosure of role + author verification**, not prohibition, and not a
  transcript/provenance archive. The standard the field is converging on is closer to "I
  used [model/version] for [specific role], and I have independently verified the results"
  than "here is the full chat log."
- Full prompt-response logs are *not* expected — Tao's repo was published because the
  paper's subject *was* a systematic tool-capability study, which is a different genre from
  organic, months-long ideation.
- **Practical implication**: rather than auditing ~30k exchanges, the tractable task is a
  per-claim audit — for each genuinely novel result/technique in the draft, note (a)
  whether/how an LLM contributed, and (b) that it's been independently verified. This is
  squarely in line with what Annals, COPE-aligned publishers, and arXiv all actually ask for.

## 5. Open questions / things not yet resolved community-wide

- No subfield-specific (analytic number theory / modular forms journal) policy was found
  beyond inheriting general publisher templates — worth re-checking the *specific* target
  journal's current guidelines close to submission, since this is moving fast.
- The boundary between "idea generation" (generally requires disclosure) and "AI-assisted
  copy editing" (generally exempt) is still publisher-dependent and not yet settled for
  cases where an LLM contributed a *technique* (e.g., finding the Hurwitz-kernel
  reformulation) rather than text.
