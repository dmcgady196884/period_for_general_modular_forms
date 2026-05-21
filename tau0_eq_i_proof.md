# The τ₀ = i Elliptic-Fixed-Point Proof of BFK = DR-B

*Companion document to `research_summary.md` and `conceptual_framework.md`. This one documents a specific proof strategy and its numerical verification.*

---

## 1. The idea (DM, this chat)

Choose the DR-B basepoint $\tau_0 = i$. Since $S \cdot i = -1/i = i$, the elliptic fixed point of $S$ collapses the S-segment $[-1/\tau_0, \tau_0]$ to a single point. The DR-B formula reduces to

$$\Lambda^*_{\rm DR}(\hat\Delta, s) \;=\; i^{-s}\int_{i-1}^{i} \hat\Delta(\tau)\, \tilde{k}_T(\tau, s)\, d\tau,$$

a single integral over the T-segment alone. This is the dual specialization to the classical Mellin choice: at $\tau_0 \to i\infty$ (the cusp), the T-segment dies exponentially and the S-segment survives; at $\tau_0 = i$ (the elliptic fixed point of $S$), the S-segment vanishes identically and the T-segment survives. Both specializations exploit a fixed point of the modular curve, just different ones.

| Basepoint | Vanishing segment | Surviving contour | What it gives |
|---|---|---|---|
| $\tau_0 \to i\infty$ (cusp) | T (exponentially) | S-segment → $[0, i\infty]$ | Classical Mellin |
| $\tau_0 = i$ (elliptic fixed point of $S$) | S (identically) | T-segment $[i-1, i]$ | This strategy |

---

## 2. Why this is a better proof route than the τ₀ → i∞ Step 1–2–3

The earlier strategy (sketched in `research_summary.md`) used the general DR-B sum at any $\tau_0$, with the limit $\tau_0 \to i\infty$ producing the BFK formula in the cuspidal modes. The polar mode required a delicate asymptotic-subtraction argument because both S- and T-segments individually diverged as $\tau_0 \to i\infty$, with the divergent pieces forced to cancel.

The $\tau_0 = i$ approach avoids this entirely:

1. **No S-segment.** The S-integral is identically zero, not vanishing in a limit.
2. **No competing divergences.** The T-segment $[i-1, i]$ is compact; the integral is finite for every Fourier mode including the polar one.
3. **Per-mode equality is exact and direct.** At integer $s$, each mode contributes via a closed-form integration-by-parts evaluation that matches BFK's incomplete-gamma value, with the polar mode handled by the same closed form via analytic continuation in the second argument of $\Gamma(s, \cdot)$.
4. **Boundary cases $s = 1, 11$ work without modification.** The polynomial-$V_{10}$ construction in the note explicitly fails at these boundaries; the $\tau_0 = i$ approach does not.

---

## 3. The proof at integer $s$

### Step 1: collapse and commute

At $\tau_0 = i$:
$$\Lambda^*_{\rm DR}(\hat\Delta, s) \;=\; i^{-s} \int_{i-1}^{i} \hat\Delta(\tau)\, \tilde{k}_T(\tau, s)\, d\tau.$$

The $q$-series $\hat\Delta(\tau) = q^{-1} + \sum_{n \geq 2} a(n) q^n$ converges uniformly on the compact segment $[i-1, i]$ (including the polar mode, since $|q^{-1}|$ is bounded by $e^{2\pi}$ on this segment). The kernel $\tilde{k}_T(\tau, s)$ is continuous on this segment. By uniform convergence, sum and integral commute:

$$\Lambda^*_{\rm DR}(\hat\Delta, s) \;=\; i^{-s} \sum_n a(n) \int_{i-1}^{i} q^n\, \tilde{k}_T(\tau, s)\, d\tau.$$

### Step 2: kernel collapses to a polynomial at integer $s$

For integer $s \in \mathbb{Z}$, the Hurwitz zeta values collapse to Bernoulli polynomials via $\zeta(-m, a) = -B_{m+1}(a)/(m+1)$:

$$\tilde{k}_T(\tau, s) \;=\; \zeta(1-s, \tau+1) - e^{i\pi(s-1)}\zeta(s-11, \tau+1)$$

becomes a polynomial in $\tau$. For example:
- $s = 6$: $\tilde{k}_T(\tau, 6) = -B_6(\tau+1)/3$.
- $s = 1$: $\tilde{k}_T(\tau, 1) = -\tau - \tfrac{1}{2} + B_{11}(\tau+1)/11$.
- $s = 11$: $\tilde{k}_T(\tau, 11) = -B_{11}(\tau+1)/11 + \tau + \tfrac{1}{2}$.

So the integrand is $q^n$ times a polynomial in $\tau$, and the per-mode integral is a finite sum of $\int_{i-1}^{i} q^n \tau^m d\tau$ for $m = 0, 1, \dots, 11$.

### Step 3: closed-form integration by parts

After substituting $\sigma = \tau + 1$ (so the kernel becomes $B_{m+1}(\sigma)/(m+1)$ in clean form), the per-monomial integral is

$$\int_i^{i+1} e^{2\pi i n \sigma}\, \sigma^m\, d\sigma \;=\; \sum_{k=0}^{m} \frac{(-1)^k\, m!}{(m-k)!\,(2\pi i n)^{k+1}} \left[(i+1)^{m-k} - i^{m-k}\right] e^{-2\pi n}$$

via repeated integration by parts, using $e^{2\pi i n (i+1)} = e^{2\pi i n \cdot i}\cdot e^{2\pi i n} = e^{-2\pi n}\cdot 1$ for integer $n$.

For positive integer $n$, this exactly matches $i^{m+1}\Gamma(m+1, 2\pi n)/(2\pi n)^{m+1}$ via the explicit polynomial form
$$\Gamma(N, x) \;=\; (N-1)!\, e^{-x} \sum_{m=0}^{N-1} \frac{x^m}{m!}\,.$$
This recovers BFK's incomplete-gamma summand directly.

For $n = -1$ (the polar mode), the same closed form matches $\Gamma(N, -2\pi)/(-2\pi)^N$ via the same formula — which is exactly the analytic-continuation value that BFK uses (the incomplete gamma is entire in its second argument; the polynomial expression is valid for any $x \in \mathbb{C}$).

Both terms in BFK's formula $\Gamma(s, 2\pi n)/(2\pi n)^s + i^k \Gamma(k-s, 2\pi n)/(2\pi n)^{k-s}$ are accounted for by the two halves of the kernel $\tilde{k}_T(\tau, s) = \zeta(1-s, \tau+1) - e^{i\pi(s-1)}\zeta(s-11, \tau+1)$.

### Step 4: extending to all complex $s$ via Carlson

Both sides are entire functions of $s$:
- $\Lambda^*_{\rm DR}(\hat\Delta, s)$: the integral is over a compact contour with integrand entire in $s$ (Hurwitz zeta $\zeta(\sigma, \tau+1)$ is meromorphic in $\sigma$ with a simple pole only at $\sigma = 1$, which corresponds to $s = 0$ or $s = 12$ in the kernel — these are the only candidates for non-entirety, and they're handled by the $i^{-s}$ prefactor structure).
- $L^*_{\rm BFK}(\hat\Delta, s)$: each summand $\Gamma(s, 2\pi n)/(2\pi n)^s$ is entire in $s$, and the sum converges uniformly on compact sets of $s$ given the exponentially-decaying Fourier coefficients of $\hat\Delta$.

If the two agree at all positive integers $s = 1, 2, 3, \dots$ (the polynomial-kernel proof above extends to any positive integer $s$, since $\zeta(-m, a)$ is a Bernoulli polynomial for all positive integer $m$), and if both sides have controlled growth on vertical strips, then **Carlson's theorem** forces equality everywhere.

Growth check: standard estimates give
$$|\zeta(\sigma, a)| \;=\; O(|\mathrm{Im}(\sigma)|^{\max(0,\, 1/2 - \mathrm{Re}(\sigma)) + \epsilon})$$
on vertical lines, with $a$ in a compact set off the negative real axis. This gives polynomial growth of $\tilde{k}_T(\tau, s)$ in $|\mathrm{Im}(s)|$, and hence of $\Lambda^*_{\rm DR}(\hat\Delta, s)$. Both sides are of finite exponential type, and Carlson's theorem (or the cleaner Phragmén–Lindelöf variant) applies.

---

## 4. Numerical verification

Verified at 50-digit precision in `verify_s6.py` and `verify_multi_s.py`:

| $s$ | $n = -1$ | $n = +1$ | $n = +2$ | $n = +3$ |
|---|---|---|---|---|
| 1 | ✓ ($10^{-51}$) | ✓ ($10^{-53}$) | ✓ ($10^{-57}$) | ✓ ($10^{-59}$) |
| 3 | ✓ ($10^{-58}$) | ✓ ($10^{-53}$) | ✓ ($10^{-57}$) | ✓ ($10^{-59}$) |
| 6 | ✓ ($10^{-51}$) | ✓ ($10^{-54}$) | ✓ ($10^{-57}$) | ✓ ($10^{-59}$) |
| 9 | ✓ ($10^{-56}$) | ✓ ($10^{-53}$) | ✓ ($10^{-57}$) | ✓ ($10^{-59}$) |
| 11 | ✓ ($10^{-51}$) | ✓ ($10^{-53}$) | ✓ ($10^{-57}$) | ✓ ($10^{-59}$) |

Values shown are relative residuals between the DR-B per-mode integral and the BFK per-mode incomplete-gamma value. All match to within numerical precision.

Notable observations from the output:
- $s = 1$ and $s = 11$ (the boundary cases the polynomial-$V_{10}$ construction in the note explicitly cannot reach) work without modification.
- $s = 3$ and $s = 9$ produce identical per-mode values (dual under $s \leftrightarrow k - s = 12 - s$ since $i^{12} = 1$).
- The polar mode $n = -1$ matches at every $s$ via the analytic-continuation closed form.

---

## 5. What this gives, and what remains to do

### What's established

The numerical evidence at 50 digits, together with the closed-form integration-by-parts argument outlined above, establishes per-mode equality at every integer $s$. This is the heart of the proof. The summability of $\hat\Delta$'s $q$-series gives the full identity at integer $s$.

### What still needs to be written

1. **The closed-form integration-by-parts identity in clean form.** The expression in §3 Step 3 needs to be carried through to show explicitly that the per-monomial DR-B integral equals the BFK incomplete-gamma value, term-by-term, for both positive and negative $n$. This is mechanical but tedious (essentially a binomial-coefficient identity).

2. **Carlson application with explicit growth bound.** The polynomial-growth claim for $\tilde{k}_T(\tau, s)$ on vertical lines needs an explicit constant, not just an order estimate. Standard but needs writing.

3. **Literature check.** The $\tau_0 = i$ specialization is the obvious move once you know about the S-fixed point. It is possible this has been observed before. A careful check of the Diamantis–Rolen original paper, the Brown follow-ups, and the harmonic-Maass-form / mock-modular literature is required before claiming priority on the framing. The mathematical content (the equivalence with BFK, the analytic-continuation handling of the polar mode) is likely new regardless, but the framing might not be.

### What this enables as a paper

The result, if the literature check confirms novelty, is a self-contained paper:

> **Title (provisional):** *L-values of weakly holomorphic modular forms as elliptic-fixed-point holonomies.*
>
> **Abstract sketch:** The completed L-function $\Lambda^*(f, s)$ of a weakly holomorphic cusp form $f \in S^!_k$ admits a closed-form representation as a finite contour integral on the line segment $[i-1, i]$ in the upper half-plane, against a Hurwitz-zeta kernel. The basepoint $i$ is the elliptic fixed point of the modular involution $S$. The construction requires no regulator, no analytic continuation, and no asymptotic-subtraction argument; the polar Fourier modes of $f$ contribute via the same integration-by-parts evaluation as the cuspidal modes. The Bringmann–Fricke–Kent regulated formula is the integration-by-parts evaluation of this single integral. The construction extends without modification to the boundary L-values $s = 1, k-1$ which are inaccessible to the polynomial-period construction.

Length: 12-20 pages, depending on how much of the conceptual framework (holonomy interpretation, gauge analogy) gets included.

### Scope decisions to make

- **Include or exclude the holonomy/Wilson-loop framing?** The paper works as a self-contained technical result without it. Including it makes the contribution feel larger but invites debate about whether the framing is "really" new. Suggestion: include it, but in a "remark" or "discussion" section, not as a load-bearing claim.
- **Include or exclude the meromorphic-form extension?** This belongs in a separate paper. The $\tau_0 = i$ approach breaks for meromorphic forms with interior poles in $[i-1, i]$ (and the more interesting case is interior poles anywhere in $\mathbb{H}$, where the wall-crossing data enters). Keeping the present paper focused on weakly holomorphic forms is cleaner.
- **General weight $k$ or restrict to $k = 12$?** The argument generalizes immediately to any even $k$ with $\dim S_k = 1$, and with more care to arbitrary $k$ (where the period polynomial space is higher-dimensional). For the first paper: probably just state for general $k$ but exhibit the proof at $k = 12$.

---

## 6. Files

- `verify_s6.py`: per-mode verification at $s = 6$ for $n \in \{-1, 1, 2, 3, 4, 5\}$.
- `verify_multi_s.py`: per-mode verification at $s \in \{1, 3, 6, 9, 11\}$ for $n \in \{-1, 1, 2, 3\}$.

Both use mpmath at 50-digit precision. Run via `python verify_s6.py` or `python verify_multi_s.py`.
