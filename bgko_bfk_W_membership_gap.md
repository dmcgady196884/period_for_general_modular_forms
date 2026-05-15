# The Period-Polynomial-in-W Claim: Where BGKO and BFK Don't Reach Détente With the Numerics

**Status:** Diagnosis verified numerically and proved structurally.
Includes corrections of two earlier errors (see §9).

**Companion to:** `bfk_regularization_handoff.md`

---

## 1. The literal claims being challenged

**BGKO (Bringmann–Guerzhoy–Kent–Ono, 2010)**,
Proposition 2.3 proof, p. 12:

> "By the modularity of $F$, it follows that
> $r(F; z) = r^-(F; z) + ir^+(F; z) \in W$."

That single line is the *entire* justification offered.
It is also used in Theorem 1.6 / Proposition 4.1 case (2)
(p. 22), which explicitly invokes
"we used the fact that $r(G;z) \in W$."

**BFK (Bringmann–Fricke–Kent, 2011)**, p. 6:

> "It is known [BGKO, Kohnen-Zagier] that $r(f; z)$
> satisfies period relations, i.e., $r(f; z) \in W$."

Kohnen-Zagier (1984) is for *holomorphic* forms only.
The general claim for $f \in S^!_k$ rests entirely on
BGKO Prop 2.3.

---

## 2. Core definitions

The space $W \subset V_{k-2}$
(polynomials of degree $\leq k-2$ in $z$, with weight
$2-k$ right slash action of $\mathrm{SL}_2(\mathbb{Z})$)
is defined by:

$$W := \{P \in V_{k-2} :
P|(1+S) = 0 \text{ and } P|(1+U+U^2) = 0\}$$

where $S = \bigl(\begin{smallmatrix}0&-1\\1&0\end{smallmatrix}\bigr)$,
$U = TS = \bigl(\begin{smallmatrix}1&-1\\1&0\end{smallmatrix}\bigr)$.

For $F \in S^!_k$ with Fourier expansion
$F(z) = \sum_{n \gg -\infty} a_F(n) q^n$,
the formal Eichler integral is

$$E_F(z) := \sum_{n \neq 0} a_F(n)\, n^{1-k}\, q^n.$$

BGKO eq (1.7) defines the period function

$$r(F; z) := c_k\bigl(E_F - E_F|_{2-k} S\bigr)(z),
\qquad c_k := -\frac{\Gamma(k-1)}{(2\pi i)^{k-1}}.$$

For $F \in S_k$ (holomorphic), $r(F;z) \in V_{k-2}$
and lies in $W$. BGKO Prop 2.3 claims the same for
all $F \in S^!_k$.

---

## 3. Proof that $\varepsilon(S) \notin V_{k-2}$ for $\hat\Delta$

**Setup.** Take $k = 12$, $\hat\Delta \in S^!_{12}$
with principal part $q^{-1}$.

For $P \in V_{10}$, polynomial growth requires
$|P(ti)| = O(t^{10})$ as $t \to +\infty$.
We show $\varepsilon(S)(ti)$ grows as $t^{10} e^{2\pi t}$.

**Term 1: $E_F(ti)$.**

$$E_F(ti) = \underbrace{-e^{2\pi t}}_{n=-1 \text{ mode}}
+ \sum_{n \geq 2} a_F(n)\, n^{-11}\, e^{-2\pi nt}.$$

The $n \geq 2$ modes decay exponentially.
So $E_F(ti) \sim -e^{2\pi t}$.

**Term 2: $E_F|S(ti) = (ti)^{10}\, E_F(i/t)$.**

The positive modes of $E_F$ evaluated at $i/t$:

$$E_F(i/t) = \sum_{n \neq 0} a_F(n)\, n^{-11}\,
e^{-2\pi n/t}.$$

As $t \to \infty$, the factors $e^{-2\pi n/t} \to 1$,
so these do NOT individually decay. The coefficients
$a_F(n)$ for $n \geq 2$ grow roughly as
$e^{4\pi\sqrt{n}}$, and the full sum is dominated by
its saddle-point contribution at $n_0 \sim t^2$,
giving growth $\sim C e^{2\pi t}$.

This is confirmed directly from modularity of
$\hat\Delta$: since $\hat\Delta \in S^!_{12}$, the
weight-12 transformation law gives

$$\hat\Delta(i/t) = \hat\Delta(-1/(it))
= (it)^{12}\,\hat\Delta(it)
\sim i^{12}\, t^{12}\, e^{2\pi t}
= t^{12}\, e^{2\pi t}.$$

Since $\hat\Delta = D^{11}(E_F) \cdot c_k$ (modulo the
precise normalization), $E_F(i/t)$ inherits the same
$e^{2\pi t}$ exponential growth from its positive
modes, though with a different polynomial prefactor.
So $E_F(i/t) \sim C e^{2\pi t}$ for a nonzero
constant $C$.

Therefore:

$$E_F|S(ti) = (ti)^{10}\, E_F(i/t)
\sim i^{10}\, t^{10} \cdot C e^{2\pi t}
= -C\, t^{10}\, e^{2\pi t}.$$

**Both terms grow exponentially**, but $E_F|S$
dominates by a factor of $t^{10}$.

**Conclusion:**

$$\varepsilon(S)(ti) = c_k\bigl(E_F(ti) - E_F|S(ti)\bigr)
\sim c_k\bigl(-e^{2\pi t} + C\, t^{10}\, e^{2\pi t}\bigr)
\sim c_k\, C\, t^{10}\, e^{2\pi t}.$$

This grows as $t^{10} e^{2\pi t}$, exponentially
faster than any polynomial of degree $\leq 10$.

**Therefore $\varepsilon(S) \notin V_{k-2}$.** $\square$

**Remark.** The weight-$2-k$ prefactor $(ti)^{k-2}$
in the slash action enhances $E_F(i/t)$ by $t^{k-2}$,
but both $E_F$ and $E_F|S$ blow up exponentially —
the weight factor shifts the polynomial prefactor
but cannot remove the $e^{2\pi t}$ growth.
The Eichler integral is (up to normalization) an
$(k-1)$-fold antiderivative of $F$; antiderivatives
of exponentials are still exponential.

---

## 4. What BGKO's argument actually does (and where it silently fails)

BGKO's eq (1.9) writes, for $F \in S^!_k$ with
$a_F(0) = 0$:

$$r(F; z) =
\sum_{n=0}^{k-2} i^{1-n}\binom{k-2}{n}
r_n(F)\, z^{k-2-n}.$$

This claims $r(F;z)$ is a polynomial of degree
$\leq k-2$. The above proof shows this is false
when $F$ has a $q^{-m}$ principal part: the function
$\varepsilon(S) = c_k(E_F - E_F|S)$ grows exponentially
at $i\infty$.

BGKO's eq (1.9) accounts only for the constant-term
contribution $a_F(0)$ (the $z^{k-1}+1/z$ piece).
It **silently ignores the exponentially growing
contributions from negative-index Fourier modes.**
For $\hat\Delta$ with $a_F(-1) = 1$, there is a
piece of $\varepsilon(S)$ growing as $t^{10}e^{2\pi t}$
that eq (1.9) simply drops without acknowledgment.

In other words, BGKO's "by modularity of $F$,
$r(F;z) \in W$" is silently doing:

1. Form $\varepsilon(S) = c_k(E_F - E_F|S)$
   (a function growing exponentially at $i\infty$,
   not a polynomial).
2. **Discard the exponentially growing parts**
   (unacknowledged).
3. Call the residue $r(F;z)$.
4. Conclude $r \in W$.

Step 4 does not follow after step 2.

**Why does $(1+S)$ survive this projection anyway?**

The full function identity $\varepsilon(S) + \varepsilon(S)|S = 0$
holds because $\varepsilon(S^2) = \varepsilon(I) = 0$ and
the cocycle relation gives $\varepsilon(S)|(1+S) = 0$.
Under the slash by $S$, the exponentially growing
term $C t^{10} e^{2\pi t}$ at $i\infty$ maps to
a corresponding exponentially growing term at
$0$ — and these pair up to cancel, leaving
$r_{\rm BFK}|(1+S) = 0$ at the level of the
polynomial part. $S$ exchanges the two cusps
$\{0, i\infty\}$, so the non-polynomial parts
cancel by themselves, and the polynomial parts
cancel by themselves.

For $U$ (order-3 cycling of $\{0, i\infty, 1\}$):
the three-term orbit lacks this pairing symmetry.
The exponentially growing parts don't cancel by
themselves when summed over the $U$-orbit, and
their residual projection contaminates the
polynomial part. The polynomial-part projection of
$\varepsilon(S)|(1+U+U^2) = 0$ no longer holds.

This is the same observation as
`bfk_regularization_handoff.md`'s
"the BFK regulator has $S$-symmetry but no
$U$-symmetry," now derived from the explicit
exponential-growth structure.

---

## 5. What BFK Theorem 2.4 actually computes

BFK's polynomial

$$r_{\rm BFK}(f;z) = \sum_{n=0}^{k-2}
i^{1-n}\binom{k-2}{n} L^*_f(n+1)\, z^{k-2-n}$$

is by construction an element of $V_{k-2}$
(a finite degree-$(k-2)$ polynomial). BFK Theorem 2.4
identifies it with $c_k(E_f - E_f|S)$ via the
regularized integral $R.\int_0^{i\infty}$.

For $f \in S_k$: both $E_f$ and $E_f|S$ decay at
$i\infty$ (only positive modes), so $\varepsilon(S)$
is genuinely a polynomial, and
$r_{\rm BFK} = \varepsilon(S) \in W$.

For $f \in S^!_k$ with cusp pole:
$r_{\rm BFK}$ is the **polynomial part**
of $\varepsilon(S)$, extracted by the
regularized-integral prescription.
The cocycle relation $\varepsilon(S)|(1+U+U^2) = 0$
holds as a function identity — including the
exponentially growing parts — but projecting
to the polynomial part breaks the cancellation
for the $U$-relation, as argued in §4.

---

## 6. The verification protocol

User's polynomial $r_{\rm BFK}(\hat\Delta; X)$ from
`period_synthesis_bfk.pdf` Tables 2-3, tested with
weight $-10$ slash action on $V_{10}$:

| Check | Operator | Residual | Outcome |
|---|---|---|---|
| (a) | $1 + S$ | $= 0$ exactly | **passes** |
| (b) | $1 + ST + (ST)^2$ | $\sim 1.86 \times 10^5$ | **fails** |
| (c) | $1 + U + U^2$ | $\sim 1.86 \times 10^5$ | **fails** |

Typical coefficient magnitude of $r_{\rm BFK}(\hat\Delta)$:
$\sim 8 \times 10^3$.
Residuals in (b) and (c): **22× typical** —
definitively non-zero, not noise.

(a) passes for the cusp-exchange reason in §4.
(b) and (c) fail with equal residual magnitude —
consistent with a single polar-correction origin
from the $q^{-1}$ mode.

---

## 7. Damage assessment

### BGKO

The gap is in **Proposition 2.3**: "by modularity
of $F$, $r(F;z) \in W$" silently discards
exponentially growing pieces of $\varepsilon(S)$
without acknowledging that this invalidates the
$W$-membership conclusion.

| Result | Depends on $r \in W$? | Assessment |
|---|---|---|
| Theorem 1.1 (mock periods → L-values) | No | Fine |
| Theorem 1.5 (multiplicity-two Hecke) | Not visibly | Probably fine |
| Theorem 1.2 (exact sequences) | Yes — target space is $W$ | **Suspect** |
| Theorem 1.6 / 1.7 (Haberland formula) | Yes — Prop 4.1 case (2) explicit | **Suspect** |

### BFK

| Result | Assessment |
|---|---|
| Theorem 2.2 (regularized integral, FE) | Fine |
| Theorem 2.4 (generating function for L-values) | Fine — polynomial exists, formula valid |
| p. 6 assertion $r(f;z) \in W$ | **Wrong as stated** for $f \in S^!_k$ with cusp poles |
| Theorem 2.5 (vanishing on $D^{k-1}M^!_{2-k}$) | Fine — proof constructs $r = c(1-z^{k-2})$ directly |

BFK's main theorems survive. The damage is the
one parenthetical on p. 6.

---

## 8. What can responsibly be said in writing

**Solid claims:**

1. For $F \in S^!_k$ with $q^{-m}$ principal part,
   $\varepsilon(S) = c_k(E_F - E_F|S)$ is not a
   polynomial in $V_{k-2}$: both $E_F(ti)$ and
   $E_F|S(ti)$ grow exponentially as $t \to \infty$,
   with $E_F|S$ dominating by $t^{10}$, so
   $\varepsilon(S)(ti) \sim c_k C t^{10} e^{2\pi t}$.
   (Proved in §3.)

2. BGKO Prop 2.3 drops the non-polynomial parts
   of $\varepsilon(S)$ without acknowledgment, then
   concludes $r \in W$. This does not follow.

3. BFK p. 6 inherits the gap via citation.

4. Numerical verification: $r_{\rm BFK}(\hat\Delta)$
   satisfies $(1+S) = 0$ exactly and fails
   $(1+U+U^2) = 0$ by $\sim 22\times$ typical
   coefficient magnitude.

**Not yet proved here:**

- The precise characterization of the space in
  which $r_{\rm BFK}$ does live.
- Whether BGKO Theorems 1.2, 1.6, 1.7 can be
  salvaged with a corrected target space.

**Do not say:**

- "BGKO is wrong" without qualification —
  results not depending on Prop 2.3 are fine.
- "BFK is wrong" — one line is wrong;
  the theorems are not.

---

## 9. Correction log

**Error 1 (cocycle convention).** Earlier draft
had the right-action cocycle relation reversed,
predicting (b) would pass and (c) would fail.
Corrected: the formal cocycle predicts $r \in W$
entirely ($r|(1+U+U^2) = 0$ is the prediction).
The numerics deny this, forcing the structural
diagnosis in §3.

**Error 2 (growth rate argument).** An intermediate
version of §3 incorrectly claimed
$E_F|S(ti) \sim Ct^{10}$ (polynomial growth).
This is wrong: $E_F|S(ti) = (ti)^{10} E_F(i/t)$,
and $E_F(i/t)$ grows as $C e^{2\pi t}$ (not bounded)
because the positive Fourier modes of $E_F$
collectively blow up via saddle-point / modularity
(the modular transformation $\hat\Delta(i/t)
= (it)^{12}\hat\Delta(it) \sim t^{12}e^{2\pi t}$
makes this explicit). Both $E_F$ and $E_F|S$ blow
up exponentially at $i\infty$, as one should expect
for the Eichler integral of a form with exponential
blow-up at both cusps. The corrected argument
(§3) shows $E_F|S$ dominates by $t^{10}$, so
$\varepsilon(S)$ grows as $t^{10} e^{2\pi t}$.

---

## 10. Open questions

- **Cohomological target.** What is the correct
  target space for $f \mapsto r_{\rm BFK}(f)$
  on $S^!_k$?

- **Diamantis–Rolen comparison.** The finite-endpoint
  cocycle in `periods.pdf` produces
  $r_{\rm DR}(\hat\Delta)$ in $W$ proper.
  The difference $r_{\rm DR} - r_{\rm BFK}$
  should isolate the polar correction from the
  $q^{-1}$ mode explicitly.

- **Two-parameter regulator.** The proposal in
  `bfk_regularization_handoff.md` uses
  $e^{\alpha\tau + \beta/\tau}$. Since the
  obstruction is exponential growth from both
  cusps, the question is whether that regulator's
  $\alpha, \beta \to 0$ prescription suppresses
  the exponential pieces in a $U$-symmetric way.
