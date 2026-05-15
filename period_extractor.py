"""
The goal of this function is to make a function that does the following:

1. Inputs the modular weight k, and q-series, for a specific modular form for SL2(Z), to some order in its q-series expansion
2. Outputs numerical approximations for its two periods, i.e. the two coefficients of P_S^\pm of the modular forms period polynomial

Note: Subtleties might arise in the case where the modular weight is such that the space of cusp forms is greater than unity. Not sure.

Basis, which this code will generalize:

# Periods and Quasi-Periods of $\Delta$ and $\hat\Delta$

Reproduce Brown's numerical values (1710.07912 §8.3) using the **finite-endpoint Diamantis-Rolen cocycle formula**:
$$\sigma_f(\gamma)(z) = \int_{\gamma^{-1}\tau_0}^{\tau_0} f(\tau)(\tau - z)^{k-2}\,d\tau$$

**Targets:**
- $\Delta$: $\omega^+ = -68916772.809595194754\ldots$, $\omega^- = -5585015.3793104018668\ldots$
- $\hat\Delta = \Delta(J^2 + 24 J - 393444)$: $\eta^+ = 127202100647.17709477\ldots$, $\eta^- = 10276732343.649132750\ldots$

**Procedure:**
1. Fix $\tau_0 \in \mathbb{H}$ (we use $\tau_0 = 0.3 + 1.2 i$).
2. Compute $C_\gamma(X,Y) = (2\pi i)^{11}\int_{\gamma^{-1}\tau_0}^{\tau_0} f(\tau)(X - \tau Y)^{10}\,d\tau$ for $\gamma \in \{S, T\}$.
3. Solve $P|_T - P = C_T$ for the coboundary polynomial $P \in V_{10}$.
4. Form $C'_S = C_S - (P|_S - P)$.
5. Decompose $C'_S = \omega^+ P^+_S + \omega^- P^-_S$ (mod $X^{10}-Y^{10}$ junk direction).
6. Read off periods.

All integration is over **finite-length line segments in $\mathbb{H}$**. Nothing approaches the cusp.

from mpmath import mp, mpc, mpf, pi, exp, quad, matrix, lu_solve
from sympy import Rational
from math import comb
import time

mp.dps = 40
n = 10
TWOPII = mpc(0, 2*pi)
PREFACTOR = TWOPII**(n+1)
print(f'Working at {mp.dps} decimal digits.')

## 1. Build $q$-series for $\Delta$, $j$, $J = j-744$, and $\hat\Delta$

$\hat\Delta = \Delta \cdot (J^2 + a_1 J + a_0)$, with $a_1, a_0$ determined by demanding
$\hat\Delta = q^{-1} + 0 + 0\cdot q + O(q^2)$.

Since $\Delta$ has $a_0 = 0$ and $a_1 = 1$, the linear conditions read:
- $q^0$ of $\hat\Delta = 0$: $\;[\Delta J^2]_0 + a_1 [\Delta J]_0 = 0$.
- $q^1$ of $\hat\Delta = 0$: $\;[\Delta J^2]_1 + a_1 [\Delta J]_1 + a_0 = 0$.

We work with exact rationals.

PREC_QSERIES = 60

def laurent_mul(a, b, lo, hi):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            nn = ka + kb
            if lo <= nn <= hi:
                out[nn] = out.get(nn, Rational(0)) + va * vb
    return {k: v for k, v in out.items() if v != 0}

def laurent_div(num, den, lo_exp, hi_exp):
    den_min = min(den.keys())
    den_lead = den[den_min]
    quot = {}
    for target_n in range(lo_exp + den_min, hi_exp + den_min + 1):
        s = num.get(target_n, Rational(0))
        for k_quot in range(lo_exp, target_n - den_min):
            s -= quot.get(k_quot, Rational(0)) * den.get(target_n - k_quot, Rational(0))
        if s != 0:
            quot[target_n - den_min] = s / den_lead
    return quot

def divisor_sigma(nn, k):
    s = 0; d = 1
    while d*d <= nn:
        if nn % d == 0:
            s += d**k
            if d != nn // d:
                s += (nn // d)**k
        d += 1
    return s

# Delta = q * prod (1-q^n)^24
def eta24_qexp(prec):
    series = {0: Rational(1)}
    for nn in range(1, prec + 1):
        factor = {}
        for k in range(25):
            te = nn * k
            if te > prec: break
            factor[te] = Rational((-1)**k * comb(24, k))
        series = laurent_mul(series, factor, 0, prec)
    return series

eta24 = eta24_qexp(PREC_QSERIES + 2)
Delta_series = {k+1: v for k, v in eta24.items() if k+1 <= PREC_QSERIES + 2}

# E4, j, J
E4 = {0: Rational(1)}
for nn in range(1, PREC_QSERIES + 3):
    E4[nn] = Rational(240 * divisor_sigma(nn, 3))
E4_cubed = laurent_mul(laurent_mul(E4, E4, 0, PREC_QSERIES + 2), E4, 0, PREC_QSERIES + 2)
j_series = laurent_div(E4_cubed, Delta_series, -1, PREC_QSERIES)
J_series = dict(j_series); J_series[0] = J_series.get(0, Rational(0)) - 744
J_series = {k: v for k, v in J_series.items() if v != 0}
J_sq = laurent_mul(J_series, J_series, -2, PREC_QSERIES)

print('Delta:', {k: Delta_series[k] for k in sorted(Delta_series)[:5]})
print('J = j - 744:', {k: J_series[k] for k in sorted(J_series)[:5]})
print('J^2:', {k: J_sq[k] for k in sorted(J_sq)[:5]})

# Solve for a_1, a_0 in hat_Delta = Delta * (J^2 + a_1 J + a_0)
D_Jsq = laurent_mul(Delta_series, J_sq, -2, PREC_QSERIES)
D_J   = laurent_mul(Delta_series, J_series, -1, PREC_QSERIES)

# q^0 eqn:  D_Jsq[0] + a_1 * D_J[0] = 0
a_1 = -D_Jsq[0] / D_J[0]
# q^1 eqn:  D_Jsq[1] + a_1 D_J[1] + a_0 * Delta[1] = 0  (Delta[1] = 1)
a_0 = -(D_Jsq[1] + a_1 * D_J[1])
print(f'a_1 = {a_1},  a_0 = {a_0}')
assert a_1 == 24 and a_0 == -393444

# hat_Delta
hat_Delta_series = {}
for k, v in D_Jsq.items():
    hat_Delta_series[k] = hat_Delta_series.get(k, Rational(0)) + v
for k, v in D_J.items():
    hat_Delta_series[k] = hat_Delta_series.get(k, Rational(0)) + a_1 * v
for k, v in Delta_series.items():
    hat_Delta_series[k] = hat_Delta_series.get(k, Rational(0)) + a_0 * v
hat_Delta_series = {k: v for k, v in hat_Delta_series.items() if v != 0}

print('\nhat_Delta first coefficients:')
for k in sorted(hat_Delta_series)[:8]:
    print(f'  q^{k}: {hat_Delta_series[k]}')
print('\nBrown reports: q^{-1}=1, q^2=47709536, q^3=39862705122, q^4=7552626810624')

## 2. Numeric evaluators for $\Delta$ and $\hat\Delta$

$\Delta$ via the $\eta^{24}$ product (fast convergence). $\hat\Delta$ via its truncated $q$-expansion (also fast since $|q| \sim 10^{-3}$ on our paths).

def delta_form(tau, N_eta=300):
    tau = mpc(tau)
    q = exp(TWOPII * tau)
    log_prod = mpc(0)
    qn = q
    for _ in range(N_eta):
        log_prod += mp.log(1 - qn)
        qn *= q
    eta = exp(TWOPII*tau/24 + log_prod)
    return eta**24

# Convert hat_Delta to mpf coefficients
hat_coeffs_mp = sorted([(int(e), mpf(int(c))) for e, c in hat_Delta_series.items()])

def hat_delta(tau):
    tau = mpc(tau)
    q = exp(TWOPII * tau)
    val = mpc(0)
    for e, c in hat_coeffs_mp:
        val += c * (q ** e)
    return val

# Quick checks
tau_test = mpc(0, 1)
print(f'Delta(i) = {delta_form(tau_test)}')
print(f'hat_Delta(i) = {hat_delta(tau_test)}')

## 3. Cocycle integration and slash operators

Generic over the form $f$: pass in any callable `f(tau)` returning the numeric value. We evaluate
$$\int_{\tau_a}^{\tau_b} f(\tau)(X - \tau Y)^{10}\,d\tau$$
as a polynomial in $X, Y$ along a straight-line path in $\mathbb{H}$.

def cocycle_polynomial(f, tau_a, tau_b, n_max=10):
    '''(2 pi i)^{n+1} int_{tau_a}^{tau_b} f(tau) (X - tau Y)^n dtau, as coefficient list for X^{n-j}Y^j.'''
    tau_a = mpc(tau_a); tau_b = mpc(tau_b)
    dtau = tau_b - tau_a
    moments = []
    for n_pow in range(n_max + 1):
        def integrand(s, n_pow=n_pow):
            tau = tau_a + s * dtau
            return f(tau) * tau**n_pow * dtau
        moments.append(quad(integrand, [0, 1]))
    coeffs = []
    for j in range(n_max + 1):
        c = mpf(comb(n_max, j)) * mpf((-1)**j) * moments[j]
        coeffs.append(PREFACTOR * c)
    return coeffs

def slash_T(p):
    n_loc = len(p) - 1
    out = [mpf(0)] * (n_loc + 1)
    for j in range(n_loc + 1):
        for i in range(n_loc - j + 1):
            out[j + i] += p[j] * mpf(comb(n_loc - j, i))
    return out

def slash_S(p):
    n_loc = len(p) - 1
    return [p[n_loc - m] * mpf((-1)**m) for m in range(n_loc + 1)]

## 4. Generic period extraction routine

Wraps the whole pipeline: feed in $f$ and $\tau_0$, get $(\omega^+, \omega^-)$ (or $(\eta^+, \eta^-)$).

# Brown's basis (1710.07912 §8.2)
P_plus = [mpf(0)]*11
P_plus[0]  = -mpf(36)/mpf(691); P_plus[10] = mpf(36)/mpf(691)
P_plus[2]  =  mpf(1); P_plus[4]  = -mpf(3); P_plus[6]  =  mpf(3); P_plus[8]  = -mpf(1)
P_minus = [mpf(0)]*11
P_minus[1] = mpf(4); P_minus[3] = -mpf(25); P_minus[5] = mpf(42)
P_minus[7] = -mpf(25); P_minus[9] = mpf(4)

def extract_periods(f, tau0, label='f', verbose=True):
    if verbose:
        print(f'\n========== Extracting periods of {label} at tau_0 = {tau0} ==========')
    t0 = time.time()
    C_S = cocycle_polynomial(f, -1/tau0, tau0, n)
    if verbose: print(f'  C_S done ({time.time()-t0:.1f} s)')
    t0 = time.time()
    C_T = cocycle_polynomial(f, tau0 - 1, tau0, n)
    if verbose: print(f'  C_T done ({time.time()-t0:.1f} s)')

    if verbose:
        print(f'  |C_T|_max = {max(abs(c) for c in C_T)}')
        print(f'  [C_T]_X^10 = {abs(C_T[0])} (cuspidality, should be ~0)')

    # Solve P|_T - P = C_T for P
    A = matrix(10, 10); b = matrix(10, 1)
    for m in range(1, 11):
        for j in range(0, m):
            A[m-1, j] = mpf(comb(10 - j, m - j))
        b[m-1, 0] = C_T[m]
    psol = lu_solve(A, b)
    p_coeffs = [psol[j, 0] for j in range(10)] + [mpf(0)]

    PS = slash_S(p_coeffs)
    C_S_prime = [C_S[j] - (PS[j] - p_coeffs[j]) for j in range(11)]

    # Internal consistency: read off from each available monomial
    minus_estimates = [C_S_prime[j] / P_minus[j] for j in [1,3,5,7,9]]
    plus_estimates  = [C_S_prime[j] / P_plus[j]  for j in [2,4,6,8]]

    if verbose:
        print('  - estimates (odd):')
        for j, e in zip([1,3,5,7,9], minus_estimates):
            print(f'    X^{10-j}Y^{j}: {e}')
        print('  + estimates (even, interior):')
        for j, e in zip([2,4,6,8], plus_estimates):
            print(f'    X^{10-j}Y^{j}: {e}')

    return plus_estimates[0], minus_estimates[0]

## 5. Run for $\Delta$

tau0 = mpc('0.3', '1.2')
omega_plus, omega_minus = extract_periods(delta_form, tau0, 'Delta')
print('\n=== EXTRACTED ===')
print('omega^+ =', omega_plus)
print('omega^- =', omega_minus)
print('\n=== BROWN ===')
print('omega^+ = -68916772.809595194754...')
print('omega^- = -5585015.3793104018668...')
print('\n=== RATIOS (i factor expected on omega^-) ===')
print('omega^+ / Brown:', omega_plus / mpf('-68916772.809595194754'))
print('omega^- / Brown:', omega_minus / mpf('-5585015.3793104018668'))

## 6. Run for $\hat\Delta$

Same procedure, same $\tau_0$, just a different $f$.

eta_plus, eta_minus = extract_periods(hat_delta, tau0, 'hat_Delta')
print('\n=== EXTRACTED ===')
print('eta^+ =', eta_plus)
print('eta^- =', eta_minus)
print('\n=== BROWN ===')
print('eta^+ = 127202100647.17709477...')
print('eta^- = 10276732343.649132750...')
print('\n=== RATIOS (i factor expected on eta^-) ===')
print('eta^+ / Brown:', eta_plus  / mpf('127202100647.17709477'))
print('eta^- / Brown:', eta_minus / mpf('10276732343.649132750'))

## 7. Basepoint independence sanity check

Cohomology class shouldn't depend on $\tau_0$. Pick a different one and verify.

tau0_alt = mpc('-0.2', '0.9')
eta_plus_alt, eta_minus_alt = extract_periods(hat_delta, tau0_alt, 'hat_Delta (alt tau_0)', verbose=False)
print(f'eta^+ at tau_0 = {tau0}:  {eta_plus}')
print(f'eta^+ at tau_0 = {tau0_alt}: {eta_plus_alt}')
print(f'  difference: {eta_plus - eta_plus_alt}')
print()
print(f'eta^- at tau_0 = {tau0}:  {eta_minus}')
print(f'eta^- at tau_0 = {tau0_alt}: {eta_minus_alt}')
print(f'  difference: {eta_minus - eta_minus_alt}')

## Summary of conventions used

- **Cocycle**: $C_\gamma(X,Y) = (2\pi i)^{11}\int_{\gamma^{-1}\tau_0}^{\tau_0} f(\tau)(X-\tau Y)^{10}\,d\tau$. Direct port of Diamantis-Rolen eq. (2.16) bottom formula, with cusp form normalization.
- **Slash actions**: $P|_T(X,Y) = P(X+Y, Y)$, $\;P|_S(X,Y) = P(Y, -X)$.
- **Coboundary modification**: solve $P|_T - P = C_T$ for $P$, then $C'_S = C_S - (P|_S - P)$.
- **Brown basis** (Brown 1710.07912 §8.2): $P^+_S = \tfrac{36}{691}(Y^{10}-X^{10}) + X^2Y^2(X^2-Y^2)^3$, $\;P^-_S = 4X^9Y - 25X^7Y^3 + 42X^5Y^5 - 25X^3Y^7 + 4XY^9$.

**Convention discrepancy**: extracted $\omega^-, \eta^-$ come out as $i$ times Brown's reported real values. This matches Brown's own period-matrix convention in eq. (2.13) where the entries are $\omega^+, i\omega^-$ — i.e., he absorbs the $i$ into the basis on one side and not the other depending on context. The magnitudes match to the full reported precision.
"""