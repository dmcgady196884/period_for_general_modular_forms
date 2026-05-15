"""
Numerical verification that the Bernoulli closed-form formulas in
analytic_periods_weight12.tex (eq:P-bernoulli, eq:p-explicit) agree with
the upper-triangular recursion solver used in the notebook.

For Δ at τ₀ = 0.3 + 1.2i, compares:
  - p_0, p_1, p_2 from solve_coboundary (recursion)
  - p_0, p_1, p_2 from explicit closed forms in .tex
  - The full P from the Bernoulli operator formula
"""
from mpmath import mp, mpc, mpf, pi, exp, quad, matrix, lu_solve, bernoulli, factorial
from sympy import Rational
from math import comb

mp.dps = 45
TWOPII = mpc(0, 2) * pi
PREF   = TWOPII ** 11

# ---- Δ q-series ----
def laurent_mul(a, b, lo, hi):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = ka + kb
            if lo <= k <= hi:
                out[k] = out.get(k, Rational(0)) + va * vb
    return {k: v for k, v in out.items() if v != 0}

def eta24_qexp(prec):
    series = {0: Rational(1)}
    for n in range(1, prec + 1):
        factor = {}
        for k in range(25):
            te = n * k
            if te > prec: break
            factor[te] = Rational((-1)**k * comb(24, k))
        series = laurent_mul(series, factor, 0, prec)
    return series

eta24 = eta24_qexp(72)
Delta_series = {k + 1: v for k, v in eta24.items() if 1 <= k + 1 <= 72}
terms = sorted([(int(n), mpf(str(c))) for n, c in Delta_series.items() if c != 0])

def Delta_eval(tau):
    tau = mpc(tau)
    val = mpc(0)
    for n, c in terms:
        val += c * exp(n * TWOPII * tau)
    return val

tau0 = mpc('0.3', '1.2')

# ---- Moments M_j = ∫_{τ₀-1}^{τ₀} Δ(τ) τ^j dτ ----
def Mj(j):
    a, b = tau0 - 1, tau0
    dtau = b - a
    def integrand(s, j=j):
        tau = a + s * dtau
        return Delta_eval(tau) * tau**j * dtau
    return quad(integrand, [0, 1])

print("Moments of Δ at τ₀ = 0.3 + 1.2i:")
M = [Mj(j) for j in range(11)]
for j in range(4):
    print(f"  M_{j} = {M[j]}")

# ---- C_T via the cocycle_poly routine ----
def cocycle_poly(tau_a, tau_b, n=10):
    a, b = mpc(tau_a), mpc(tau_b)
    dtau = b - a
    cs = []
    for j in range(n + 1):
        def integrand(s, j=j):
            tau = a + s * dtau
            return Delta_eval(tau) * tau**j * dtau
        moment = quad(integrand, [0, 1])
        cs.append(PREF * mpf((-1)**j) * mpf(comb(n, j)) * moment)
    return cs

C_T = cocycle_poly(tau0 - 1, tau0)
print(f"\n[C_T]_{{X^10}} = {C_T[0]}  (cuspidality, ≈ 0)")

# ---- p from recursion (notebook's solve_coboundary) ----
def solve_coboundary(C_T, n=10):
    A = matrix(n, n)
    b = matrix(n, 1)
    for m in range(1, n + 1):
        for j in range(m):
            A[m-1, j] = mpf(comb(n - j, m - j))
        b[m-1, 0] = C_T[m]
    psol = lu_solve(A, b)
    return [psol[j, 0] for j in range(n)] + [mpc(0)]

p_recursion = solve_coboundary(C_T)

# ---- p from explicit closed forms (.tex eq:p-explicit) ----
p0_tex = -PREF * M[1]
p1_tex = 5 * PREF * (M[1] + M[2])
p2_tex = -PREF * (mpf('7.5') * M[1] + mpf('22.5') * M[2] + 15 * M[3])

# ---- p_j (j = 0..9) from full Bernoulli operator formula ----
# P = (1/Y) ∫_0^X C_T(t,Y) dt + Σ_{k=1}^{10} (B_k/k!) Y^{k-1} ∂_X^{k-1} C_T(X,Y)
#
# Coefficient of X^{10-m} Y^m in P, for m = 0, ..., 10:
#   antideriv: c_{m+1} / (10-m)   (for m = 0, ..., 9; from c_j X^{11-j} Y^{j-1}/(11-j) with j = m+1)
#   k-th term: (B_k/k!) · c_{m-k+1} · (9-m+k)! / (10-m)!  (for k = 1, ..., m+1; needs c_{m-k+1} defined, i.e., m-k+1 in 0..10)

def bernoulli_pm(m):
    """Coefficient of X^{10-m} Y^m in P via the closed-form Bernoulli operator."""
    val = mpc(0)
    # antiderivative term: only for m = 0,...,9
    if m <= 9:
        val += C_T[m + 1] / mpf(10 - m)
    # Bernoulli sum
    for k in range(1, m + 2):
        idx = m - k + 1
        if 0 <= idx <= 10:
            val += (mpf(bernoulli(k)) / mpf(factorial(k))
                    * C_T[idx]
                    * mpf(factorial(9 - m + k)) / mpf(factorial(10 - m)))
    return val

p_bernoulli = [bernoulli_pm(m) for m in range(10)] + [mpc(0)]

# ---- Compare ----
print("\n" + "─" * 65)
print("Comparison (working precision: 45 digits)")
print("─" * 65)

for m, name in [(0, "p_0"), (1, "p_1"), (2, "p_2")]:
    print(f"\n{name}:")
    print(f"  recursion        = {p_recursion[m]}")
    print(f"  Bernoulli (full) = {p_bernoulli[m]}")
    if m == 0: tex_val = p0_tex
    if m == 1: tex_val = p1_tex
    if m == 2: tex_val = p2_tex
    print(f"  .tex closed-form = {tex_val}")
    diff_b = abs(p_recursion[m] - p_bernoulli[m])
    diff_t = abs(p_recursion[m] - tex_val)
    print(f"  |rec − Bernoulli| = {float(diff_b):.3e}")
    print(f"  |rec − .tex|      = {float(diff_t):.3e}")

# Compare all p_j (m=0..9) for the full Bernoulli formula:
print("\nAll p_m (m=0..9):  |recursion − Bernoulli formula|:")
for m in range(10):
    d = abs(p_recursion[m] - p_bernoulli[m])
    print(f"  m={m}:  {float(d):.3e}")
