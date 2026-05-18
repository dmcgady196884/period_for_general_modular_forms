"""
DR-B extension of period polynomial machinery to extract L(Δ, s) at complex s.

Plan:
 (A) Reference: L*(Δ,s) = (2π)^{-s} Γ(s) L(Δ,s) via the entire BFK incomplete-Γ
     formula split at τ=i.  This is our ground-truth canonical L-function.
 (B) Integer-s sanity: the existing DR-B period-polynomial extraction, with the
     proper normalization, must reproduce Λ*(Δ, j+1) for j ∈ {0,..,10}.
 (C) Complex-s extension: replace τ^j → τ^{s-1} in the S-kernel and Bernoulli
     numbers B_n → Hurwitz-zeta continuation in the T-kernel.  Evaluate the
     DR-B linear functional on Δ at complex s and compare to Λ_class.
"""
from mpmath import (mp, mpc, mpf, pi, exp, quad, gammainc, gamma as mpgamma,
                    zeta as mpzeta, factorial, rgamma)
from sympy import Rational
from math import comb

mp.dps = 50
TWOPII = mpc(0, 2) * pi
TAU0   = mpc('0.3', '1.2')
K      = 12        # weight
N      = K - 2     # = 10  (degree of period polynomial)

# ============================================================================
# Δ q-series and evaluator
# ============================================================================
def laurent_mul(a, b, lo, hi):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            kk = ka + kb
            if lo <= kk <= hi:
                out[kk] = out.get(kk, Rational(0)) + va * vb
    return {kk: v for kk, v in out.items() if v != 0}

def build_Delta(prec):
    series = {0: Rational(1)}
    for n in range(1, prec + 1):
        factor = {}
        for k in range(25):
            te = n * k
            if te > prec: break
            factor[te] = Rational((-1)**k * comb(24, k))
        series = laurent_mul(series, factor, 0, prec)
    return {kk + 1: v for kk, v in series.items() if 1 <= kk + 1 <= prec + 1}

PREC_Q = 60
Delta_series = build_Delta(PREC_Q)
delta_terms = sorted([(int(n), mpf(str(c))) for n, c in Delta_series.items() if c != 0])
print(f"Δ q-series: {len(delta_terms)} nonzero modes up to q^{PREC_Q}")

def Delta_eval(tau):
    tau = mpc(tau)
    q   = exp(TWOPII * tau)
    val = mpc(0)
    qn  = mpc(1)
    n_prev = 0
    for n, c in delta_terms:
        while n_prev < n:
            qn *= q
            n_prev += 1
        val += c * qn
    return val

# ============================================================================
# (A) Reference Λ*(Δ, s) via BFK incomplete-Γ split at τ=i.
#     Entire in s; functional equation Λ*(s) = Λ*(K - s).
# ============================================================================
def Lambda_class(s):
    s = mpc(s)
    total = mpc(0)
    for n, c in delta_terms:
        x = 2 * pi * n
        term = gammainc(s, x) / x**s + gammainc(K - s, x) / x**(K - s)
        total += c * term
    return total

print("\n=== (A) BFK reference Λ*(Δ, s) ===")
print(f"  Functional eq check: |Λ*(2) - Λ*(10)| = {float(abs(Lambda_class(2) - Lambda_class(K-2))):.3e}")
for s_int in [1, 2, 3, 4, 5, 6]:
    print(f"  Λ*(Δ, {s_int}) = {Lambda_class(s_int)}")

# Also try at non-integer s to confirm Λ_class is finite there.
print("\n  At non-integer s (entire function — must be finite):")
for s_test in [mpc('2.5'), mpc('5.5'), mpc('6'), mpc('1.5', '0.7')]:
    print(f"  Λ*(Δ, {s_test}) = {Lambda_class(s_test)}")

# ============================================================================
# (B) Integer-s sanity check: the existing DR-B linear functional reproduces
#     Λ*(Δ, j+1) at all 11 integer-j positions.
# ============================================================================
def cocycle_moment(f_eval, tau_a, tau_b, s):
    """∫_a^b f(τ) τ^{s-1} dτ along the straight line."""
    tau_a, tau_b = mpc(tau_a), mpc(tau_b)
    dtau = tau_b - tau_a
    def integrand(t):
        tau = tau_a + t * dtau
        return f_eval(tau) * tau**(s - 1) * dtau
    return quad(integrand, [0, 1])

def cocycle_polynomial(f_eval, tau_a, tau_b):
    """Returns coefficients c_j = (2πi)^{n+1} (-1)^j C(n,j) ∫ f(τ) τ^j dτ, for j = 0..n."""
    coeffs = []
    PREF = TWOPII ** (N + 1)
    for j in range(N + 1):
        mom = cocycle_moment(f_eval, tau_a, tau_b, s=j+1)
        coeffs.append(PREF * mpf((-1)**j) * mpf(comb(N, j)) * mom)
    return coeffs

def slash_S(p):
    """[p_0,..,p_N] → coefficients of (P|S)(X) = X^N P(-1/X), for N even (here N=10)."""
    n = len(p) - 1
    return [p[n - m] * mpf((-1)**m) for m in range(n + 1)]

def solve_coboundary(C_T):
    """Solve P (polynomial coefficients) so that (P|T - P)_j = C_T[j], j=1..N."""
    from mpmath import matrix as mp_matrix, lu_solve
    A = mp_matrix(N, N); b = mp_matrix(N, 1)
    for m in range(1, N + 1):
        for j in range(m):
            A[m - 1, j] = mpf(comb(N - j, m - j))
        b[m - 1, 0] = C_T[m]
    psol = lu_solve(A, b)
    return [psol[j, 0] for j in range(N)] + [mpc(0)]

# Run the integer extraction.
print("\n=== (B) Integer-s DR-B extraction of Λ*(Δ, j+1) ===")
C_S = cocycle_polynomial(Delta_eval, -1/TAU0, TAU0)
C_T = cocycle_polynomial(Delta_eval,  TAU0 - 1, TAU0)
P   = solve_coboundary(C_T)
PS  = slash_S(P)
Cp  = [C_S[j] - (PS[j] - P[j]) for j in range(N + 1)]

# Period polynomial coefficient of X^{N-j} = (2πi)^{N+1} · C(N,j) · (-1)^{N-j} · i^{j+1} · Λ*(Δ, j+1)
# So Λ*(Δ, j+1) = Cp[j] / [(2πi)^{N+1} · C(N,j) · (-1)^{N-j} · i^{j+1}]
def Cp_to_Lambda(Cp_j, j):
    denom = TWOPII**(N + 1) * mpf(comb(N, j)) * mpf((-1)**(N - j)) * mpc(0,1)**(j + 1)
    return Cp_j / denom

print(f"  τ_0 = {TAU0}")
print(f"  {'j':>3} {'s=j+1':>5}  {'Λ*_DR(Δ,s)':>30}  {'Λ*_BFK(Δ,s)':>30}  {'|diff|':>10}")
print("  " + "-" * 95)
for j in range(N + 1):
    L_DR  = Cp_to_Lambda(Cp[j], j)
    L_ref = Lambda_class(j + 1)
    diff  = abs(L_DR - L_ref)
    # Λ* should be real for real s in critical strip (cusp form ⇒ real L*-values).
    print(f"  {j:>3} {j+1:>5}  {complex(L_DR):>30.18g}  {complex(L_ref):>30.18g}  {float(diff):>10.3e}")

# ============================================================================
# (C-prep) Quantify the "missing piece":
#   how much does the bare S-segment moment alone differ from i^s · Λ*?
#   At integer s, T-coboundary correction kills this gap; want to see it
#   structurally so we can extend to complex s.
# ============================================================================
print("\n=== (C-prep) Bare S-moment vs full Λ* — gap at integer s ===")
print(f"  {'j':>3} {'s':>3}  {'M_S/[norm]':>30}  {'Λ*_BFK':>30}  {'gap (=> needs k_T)':>22}")
print("  " + "-" * 95)
for j in range(N + 1):
    M_S = cocycle_moment(Delta_eval, -1/TAU0, TAU0, s=j+1)
    # Period-polynomial coefficient at X^{N-j} (bare S-segment only):
    #   [r_bare]_{X^{N-j}} = C(N,j)(-1)^{N-j} M_S = norm · L*_bare
    norm_factor = mpf(comb(N, j)) * mpf((-1)**(N-j)) * mpc(0,1)**(j + 1)
    Lstar_bare  = M_S / mpc(0,1)**(j + 1)  # = i^{-(j+1)} ∫ Δτ^j dτ
    # Compare directly: i^{-s} M_S vs Λ*(s).
    iminus_s_MS = M_S * mpc(0,-1)**(j+1)
    L_ref = Lambda_class(j + 1)
    gap   = iminus_s_MS - L_ref
    print(f"  {j:>3} {j+1:>3}  {complex(iminus_s_MS):>30.10g}  {complex(L_ref):>30.10g}  {complex(gap):>22.6g}")

# ============================================================================
# (C-prep-2) Verify the analytical claim:
#    -gap(s, τ_0) = BFK incomplete-Γ sum SPLIT AT τ_0 (not at τ=i).
#    This is the closed-form expression of the contour-tail-to-i∞ + tail-to-0
#    that the T-coboundary correction must absorb to make the DR-B finite-
#    contour functional reproduce Λ*(Δ, s).
# ============================================================================
def Lambda_BFK_split(s, tau_0):
    """BFK incomplete-Γ sum with splitting point at tau_0 (not at i).

    Formula (derived from i^{-s} (tail_{i∞} + tail_0) using Δ(-1/u)=u^12 Δ(u)):
       Λ* - i^{-s} M_S(τ_0)
         = Σ_n τ(n) [ Γ(s, -2πin τ_0)/(2πn)^s + Γ(12-s, -2πin τ_0)/(2πn)^{12-s} ].
    At τ_0=i: -2πin·i = 2πn (real), recovers standard BFK split at τ=i.
    """
    s = mpc(s)
    tau_0 = mpc(tau_0)
    total = mpc(0)
    for n, c in delta_terms:
        z = -TWOPII * n * tau_0       # = -2πin·τ_0
        norm_s = (2 * pi * n)**s
        norm_k = (2 * pi * n)**(K - s)
        term = gammainc(s, z) / norm_s + gammainc(K - s, z) / norm_k
        total += c * term
    return total

print("\n=== (C-prep-2) Gap vs BFK-split-at-τ_0 (closed-form tail sum) ===")
print(f"  {'j':>3} {'s':>3}  {'-gap (numerical)':>34}  {'BFK-split (closed form)':>34}  {'|diff|':>10}")
print("  " + "-" * 100)
for j in range(N + 1):
    M_S = cocycle_moment(Delta_eval, -1/TAU0, TAU0, s=j+1)
    iminus_s_MS = M_S * mpc(0,-1)**(j+1)
    L_ref       = Lambda_class(j + 1)
    neg_gap     = L_ref - iminus_s_MS
    closed      = Lambda_BFK_split(j+1, TAU0)
    diff        = abs(neg_gap - closed)
    print(f"  {j:>3} {j+1:>3}  {complex(neg_gap):>34.18g}  {complex(closed):>34.18g}  {float(diff):>10.3e}")

# ============================================================================
# (C-main) Build the complex-s T-kernel via Hurwitz-zeta extension of β,
#   integrate Δ against it over [τ_0-1, τ_0], and check whether this finite-
#   contour integral reproduces the BFK-split-at-τ_0 sum.
#
# β(m, r) at integer m:
#   β(m, r) = (-1)^r C(10,r) B_{m-r+1} (10-r)! / [(m-r+1)! (10-m)!]
#
# Hurwitz extension at complex m = s-1:
#   B_n at complex n via:  B_n = -n · ζ(1-n)  (Hurwitz/Riemann; convention B_1=+1/2).
#   To match sympy's B_1=-1/2 used in the existing project code, we subtract 1
#   from B_1 in the formula — equivalently, use B_n^- = B_n^+ - δ_{n,1}.
#
#   In smooth-complex form:
#     β(s-1, r) = -(-1)^r C(10,r) (10-r)! · ζ(r-s+1) / [Γ(s-r) Γ(12-s)]
#   with a correction at r = s (i.e., the n=1 Bernoulli case).
# ============================================================================

def bernoulli_minus(n_complex):
    """B_n at complex n, sympy convention (B_1 = -1/2).
    Uses B_n = -n·ζ(1-n) for the "+ convention", then corrects at n=1.
    Returns 1 at n=0 (handled by limit). Otherwise B_n^- = -n·ζ(1-n) - δ_{n,1}."""
    n = mpc(n_complex)
    if abs(n) < mpf('1e-30'):
        return mpc(1)
    val = -n * mpzeta(1 - n)
    # n=1 correction: B_1^- = -1/2 vs B_1^+ = +1/2
    if abs(n - 1) < mpf('1e-30'):
        val = val - 1
    return val

def beta_complex(m, r):
    """β(m, r) at complex m, integer r in {1, ..., 10}.
    Uses rgamma=1/Γ which is entire and vanishes at non-positive integers,
    naturally handling the "r > m+1" and "m > 10" cutoffs of the integer β."""
    if r < 1 or r > 10:
        return mpc(0)
    m   = mpc(m)
    Bn  = bernoulli_minus(m - r + 1)        # B_{m-r+1}
    # (-1)^r · C(10,r) · B_{m-r+1} · (10-r)! · rgamma(m-r+2) · rgamma(11-m)
    sign = mpf((-1)**r)
    binom_10_r = mpf(comb(10, r))
    fact_10_r  = mpf(int(factorial(10 - r)))
    return sign * binom_10_r * Bn * fact_10_r * rgamma(m - r + 2) * rgamma(11 - m)

def kT_complex(tau, s):
    """k_T(τ, s) = Σ_{r=1}^{10} [β(s-1, r) - (-1)^{s-1} β(11-s, r)] τ^r."""
    s = mpc(s)
    tau = mpc(tau)
    sgn = exp(mpc(0, 1) * pi * (s - 1))   # = (-1)^{s-1} at integer s
    total = mpc(0)
    for r in range(1, 11):
        c = beta_complex(s - 1, r) - sgn * beta_complex(11 - s, r)
        total += c * tau**r
    return total

def kS_complex_normalized(tau, s):
    """k_S(τ, s) WITHOUT the (-1)^{s-1} C(10, s-1) prefactor — just τ^{s-1}.
    The prefactor moves into the overall L-value normalization below."""
    return mpc(tau)**(mpc(s) - 1)

# At integer j, the existing extraction uses:
#   ω^± = (2πi)^11 / [P^±]_j · [ ∫_S Δ · (-1)^j C(10,j) τ^j dτ + ∫_T Δ · kT^j(τ) dτ ]
# In L-value form, [r_Δ]_{X^{N-j}} = (2πi)^{11} · C(10,j) · (-1)^{N-j} · i^{j+1} · Λ*(Δ, j+1).
# Working backwards, the linear functional that directly extracts Λ*(Δ, s) is
#
#   Λ*(Δ, s) = i^{-s} · [ ∫_S Δ τ^{s-1} dτ + ∫_T Δ · kT(τ, s) / [(-1)^{s-1} C(10, s-1)] dτ ]
#
# i.e. divide kT by the prefactor (-1)^{s-1} C(10, s-1) so that the S-piece is just τ^{s-1}.
#
# In compute_period_kernels.py kS *includes* (-1)^j C(10,j); kT does not. So to match,
# rescale kT by 1/[(-1)^j C(10,j)] at integer j.

def kT_normalized(tau, s):
    """kT divided by (-1)^{s-1} · C(10, s-1) — so it pairs with the bare τ^{s-1} S-kernel.
    Uses rgamma for the binomial denominator to handle integer-s poles smoothly."""
    s = mpc(s)
    sgn   = exp(mpc(0, 1) * pi * (s - 1))                              # (-1)^{s-1}
    binom = mpgamma(11) * rgamma(s) * rgamma(12 - s)                   # C(10, s-1) = 10! · 1/Γ(s) · 1/Γ(12-s)
    return kT_complex(tau, s) / (sgn * binom)

def kT_hurwitz(tau, s):
    """The 'right' complex-s T-kernel as a Hurwitz-zeta function OF τ:

       k̃_T(τ, s) = ζ(1-s, τ+1)  -  e^{iπ(s-1)} ζ(s-11, τ+1).

    Derived from the τ_0-invariance condition
       k̃_T(τ,s) - k̃_T(τ-1,s) = -τ^{s-1} + e^{iπ(s-1)} τ^{11-s}
    using Hurwitz's raising identity ζ(σ, a) - ζ(σ, a+1) = a^{-σ}.

    At integer s, both ζ-terms reduce to Bernoulli polynomials in τ,
    matching the original V_N-polynomial kernel.  At complex s, the kernel
    is a genuine Hurwitz-zeta function of τ, NOT a polynomial.
    """
    s   = mpc(s)
    tau = mpc(tau)
    sgn = exp(mpc(0, 1) * pi * (s - 1))   # = (-1)^{s-1}
    return mpzeta(1 - s, tau + 1) - sgn * mpzeta(s - 11, tau + 1)

def Lambda_DR_hurwitz(s, tau_0):
    """Λ*(Δ, s) via DR-B linear functional with the Hurwitz-zeta-of-τ T-kernel."""
    s = mpc(s); tau_0 = mpc(tau_0)
    MS = cocycle_moment(Delta_eval, -1/tau_0, tau_0, s)
    a, b = tau_0 - 1, tau_0
    dtau = b - a
    def integrand(t):
        tau = a + t * dtau
        return Delta_eval(tau) * kT_hurwitz(tau, s) * dtau
    MT = quad(integrand, [0, 1])
    return mpc(0, -1)**s * (MS + MT)

def Lambda_DR_via_kernels(s, tau_0):
    """Λ*(Δ, s) via DR-B linear functional with complex-s Hurwitz-extended T-kernel."""
    s = mpc(s); tau_0 = mpc(tau_0)
    # S-piece: ∫_{-1/τ_0}^{τ_0} Δ(τ) τ^{s-1} dτ
    MS = cocycle_moment(Delta_eval, -1/tau_0, tau_0, s)
    # T-piece: ∫_{τ_0-1}^{τ_0} Δ(τ) kT_normalized(τ, s) dτ
    a, b = tau_0 - 1, tau_0
    dtau = b - a
    def integrand(t):
        tau = a + t * dtau
        return Delta_eval(tau) * kT_normalized(tau, s) * dtau
    MT = quad(integrand, [0, 1])
    return mpc(0, -1)**s * (MS + MT)

print("\n=== (C-main) DR-B with Hurwitz-extended kT at integer s — internal consistency ===")
print(f"  {'s':>4}  {'Λ_DR (Hurwitz kT)':>30}  {'Λ_BFK reference':>30}  {'|diff|':>10}")
print("  " + "-" * 85)
for s_test in range(1, K):
    val = Lambda_DR_via_kernels(s_test, TAU0)
    ref = Lambda_class(s_test)
    diff = abs(val - ref)
    print(f"  {s_test:>4}  {complex(val):>30.18g}  {complex(ref):>30.18g}  {float(diff):>10.3e}")

print("\n=== (C-main′) DR-B with Hurwitz-zeta-OF-τ T-kernel at integer s ===")
print(f"  {'s':>4}  {'Λ_DR (Hurwitz-of-τ)':>30}  {'Λ_BFK reference':>30}  {'|diff|':>10}")
print("  " + "-" * 85)
for s_test in range(1, K):
    val = Lambda_DR_hurwitz(s_test, TAU0)
    ref = Lambda_class(s_test)
    diff = abs(val - ref)
    print(f"  {s_test:>4}  {complex(val):>30.18g}  {complex(ref):>30.18g}  {float(diff):>10.3e}")

print("\n=== (D) HEADLINE TEST: DR-B with Hurwitz-zeta-OF-τ kT at NON-INTEGER s ===")
print(f"  {'s':>20}  {'Λ_DR (Hurwitz)':>36}  {'Λ_BFK':>36}  {'|diff|':>10}")
print("  " + "-" * 110)
test_s = [
    mpc('1.5'),     mpc('2.5'),  mpc('3.5'),  mpc('4.5'),  mpc('5.5'),
    mpc('6'),       mpc('6.5'),
    mpc('2.7', '0.4'),   mpc('5', '1'),    mpc('5.5', '2.5'),
    mpc('0.5'),     mpc('-0.3'),   # outside critical strip
]
for s in test_s:
    try:
        val = Lambda_DR_hurwitz(s, TAU0)
        ref = Lambda_class(s)
        diff = abs(val - ref)
        s_str = f"{complex(s):.4g}"
        print(f"  {s_str:>20}  {complex(val):>36.18g}  {complex(ref):>36.18g}  {float(diff):>10.3e}")
    except Exception as e:
        print(f"  s={s}: ERROR {e}")

# Also: verify τ_0-independence at non-integer s.
print("\n=== (E) τ_0-independence at non-integer s ===")
print("  At a fixed s, varying τ_0 should give the same Λ_DR (because Λ* is intrinsic).")
TAU0_VARIATIONS = [mpc('0.3', '1.2'), mpc('0', '1.5'), mpc('-0.4', '0.8'), mpc('0.2', '2.0')]
for s in [mpc('2.5'), mpc('5.5'), mpc('3', '1.5')]:
    print(f"\n  s = {complex(s)}:")
    ref = Lambda_class(s)
    for t0 in TAU0_VARIATIONS:
        val = Lambda_DR_hurwitz(s, t0)
        diff = abs(val - ref)
        print(f"    τ_0 = {complex(t0)}:   Λ_DR = {complex(val)}   |diff vs Λ_BFK| = {float(diff):.3e}")
