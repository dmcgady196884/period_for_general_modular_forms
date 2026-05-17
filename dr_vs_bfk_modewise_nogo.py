from mpmath import mp, mpc, mpf, pi, exp, quad, matrix, lu_solve, gammainc, gamma as mpgamma, factorial
from sympy import Rational
from math import comb

mp.dps = 45
TWOPII = mpc(0, 2) * pi

print(f'Precision: {mp.dps} decimal digits')
print(f'2πi = {TWOPII}')

def laurent_mul(a, b, lo, hi):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = ka + kb
            if lo <= k <= hi:
                out[k] = out.get(k, Rational(0)) + va * vb
    return {k: v for k, v in out.items() if v != 0}

def build_Delta(prec):
    series = {0: Rational(1)}
    for n in range(1, prec + 1):
        factor = {}
        for k in range(25):
            te = n * k
            if te > prec: break
            factor[te] = Rational((-1)**k * comb(24, k))
        series = laurent_mul(series, factor, 0, prec)
    return {k + 1: v for k, v in series.items() if 1 <= k + 1 <= prec + 1}

PREC_Q = 20
Delta_series = build_Delta(PREC_Q)
tau_values = {n: int(Delta_series[n]) for n in range(1, 8) if n in Delta_series}
print('τ(n) for n=1..7:')
for n, t in tau_values.items():
    print(f'  τ({n}) = {t}')

def make_evaluator(coeffs_dict):
    terms = sorted([(int(n), mpf(str(c))) for n, c in coeffs_dict.items() if c != 0])
    def f(tau):
        tau = mpc(tau); val = mpc(0)
        for n, c in terms: val += c * exp(n * TWOPII * tau)
        return val
    return f

def cocycle_poly(f_eval, tau_a, tau_b, n):
    tau_a, tau_b = mpc(tau_a), mpc(tau_b); dtau = tau_b - tau_a
    PREF = TWOPII ** (n + 1); coeffs = []
    for j in range(n + 1):
        def integrand(s, j=j):
            tau = tau_a + s * dtau; return f_eval(tau) * tau**j * dtau
        moment_j = quad(integrand, [0, 1])
        coeffs.append(PREF * mpf((-1)**j) * mpf(comb(n, j)) * moment_j)
    return coeffs

def slash_S(p):
    n = len(p) - 1; return [p[n - m] * mpf((-1)**m) for m in range(n + 1)]

def solve_coboundary(C_T, n):
    A = matrix(n, n); b = matrix(n, 1)
    for m in range(1, n + 1):
        for j in range(m): A[m-1, j] = mpf(comb(n - j, m - j))
        b[m-1, 0] = C_T[m]
    psol = lu_solve(A, b)
    return [psol[j, 0] for j in range(n)] + [mpc(0)]

def brown_basis_k12():
    P_plus  = [mpf(0)] * 11; P_minus = [mpf(0)] * 11
    P_plus[0]  = -mpf(36)/mpf(691); P_plus[10] =  mpf(36)/mpf(691)
    P_plus[2]  =  mpf(1);  P_plus[4]  = -mpf(3)
    P_plus[6]  =  mpf(3);  P_plus[8]  = -mpf(1)
    P_minus[1] =  mpf(4);  P_minus[3] = -mpf(25)
    P_minus[5] =  mpf(42); P_minus[7] = -mpf(25); P_minus[9] = mpf(4)
    return P_plus, P_minus

def extract_omega(coeffs_dict, tau0=mpc('0.3', '1.2')):
    """Returns dicts {j: ω estimate} for odd j (→ ω⁻) and even j (→ ω⁺)."""
    n_deg = 10; f_eval = make_evaluator(coeffs_dict)
    C_S   = cocycle_poly(f_eval, -1/tau0, tau0,   n_deg)
    C_T   = cocycle_poly(f_eval,  tau0-1, tau0,   n_deg)
    p     = solve_coboundary(C_T, n_deg)
    PS    = slash_S(p)
    Cp    = [C_S[j] - (PS[j] - p[j]) for j in range(n_deg + 1)]
    P_plus, P_minus = brown_basis_k12()
    om = {j: Cp[j] / P_minus[j] for j in [1,3,5,7,9]}
    op = {j: Cp[j] / P_plus[j]  for j in [2,4,6,8]}
    return op, om

print('DR machinery loaded.')

def bfk_omega_single_mode(n_mode):
    """BFK contribution of mode e^{2πi·n_mode·τ} (unit coefficient) to ω±."""
    two_pi = 2 * pi; n_deg = 10
    P_plus, P_minus = brown_basis_k12()

    def L_star(s):
        x = two_pi * n_mode
        return gammainc(s, x) / x**s + gammainc(12-s, x) / x**(12-s)

    def C_prime_j(j):
        s     = j + 1
        jfact = mpf(int(factorial(j)))
        Lval  = two_pi**s / mpgamma(s) * L_star(s)
        return (two_pi)**(n_deg-j) * mpc(0,-1)**j * mpf(comb(n_deg,j)) * jfact * Lval

    om = {j: C_prime_j(j) / P_minus[j] for j in [1,3,5,7,9]}
    op = {j: C_prime_j(j) / P_plus[j]  for j in [2,4,6,8]}
    return op, om

print('BFK single-mode machinery loaded.')

BROWN_OMEGA_MINUS = mpf('-5585015.3793104018668')
TAU0 = mpc('0.3', '1.2')

# DR full sum
op_dr_full, om_dr_full = extract_omega(dict(Delta_series), tau0=TAU0)
omega_minus_DR = om_dr_full[1]   # j=1 estimate

# BFK full sum (N=150 modes)
two_pi = 2 * pi
L_star_full = mpc(0)
for mm in range(1, 151):
    c = Delta_series.get(mm)
    if c is None: continue
    c = mpf(str(c)); x = two_pi * mm
    L_star_full += c * (gammainc(2, x)/x**2 + gammainc(10, x)/x**10)   # s=2, k-s=10
n_deg = 10; j = 1; jfact = mpf(1)
Lval_full  = two_pi**2 / mpgamma(2) * L_star_full
Cp_j_full  = (two_pi)**(n_deg-j) * mpc(0,-1)**j * mpf(comb(n_deg,j)) * mpf(1) * Lval_full
_, P_minus = brown_basis_k12()
omega_minus_BFK = Cp_j_full / P_minus[j]

print('Sanity check — full sums vs Brown reference')
print(f'  Brown ref ω⁻/i  = {BROWN_OMEGA_MINUS}')
print(f'  DR  full  ω⁻/i  = {(omega_minus_DR / mpc(0,1)).real}')
print(f'  BFK full  ω⁻/i  = {(omega_minus_BFK / mpc(0,1)).real}')
print(f'  |DR  - ref| / |ref| = {float(abs((omega_minus_DR/mpc(0,1)).real - BROWN_OMEGA_MINUS)/abs(BROWN_OMEGA_MINUS)):.2e}')
print(f'  |BFK - ref| / |ref| = {float(abs((omega_minus_BFK/mpc(0,1)).real - BROWN_OMEGA_MINUS)/abs(BROWN_OMEGA_MINUS)):.2e}')

N_MODES = 5

print('Mode-by-mode comparison at j=1 (ω⁻ sector)')
print('='*100)
print(f'{"n":>3}  {"τ(n)":>8}  {"DR_n (Im part)": >28}  {"BFK_n (Im part)":>28}  {"ratio |DR|/|BFK|":>18}  {"arg(DR/BFK)": >14}')
print('-'*100)

dr_modes  = {}
bfk_modes = {}

for n_mode in range(1, N_MODES + 1):
    tau_n = int(Delta_series.get(n_mode, 0))
    op_dr,  om_dr  = extract_omega({n_mode: 1}, tau0=TAU0)
    op_bfk, om_bfk = bfk_omega_single_mode(n_mode)

    dr_modes[n_mode]  = om_dr
    bfk_modes[n_mode] = om_bfk

    dr_val  = om_dr[1]
    bfk_val = om_bfk[1]
    ratio   = dr_val / bfk_val
    from mpmath import arg, log
    print(f'{n_mode:>3}  {tau_n:>8}  {float(dr_val.imag):>28.6f}  {float(bfk_val.imag):>28.6f}  '
          f'{float(abs(ratio)):>18.6f}  {float(arg(ratio)):>14.6f} rad')

from mpmath import arg

print('Ratio DR_n / BFK_n across all interior j (n=1 and n=2)')
print('If conjecture held: all ratios = 1+0i')
print('='*80)

for n_mode in [1, 2]:
    tau_n = int(Delta_series.get(n_mode, 0))
    print(f'\n  n={n_mode}, τ({n_mode})={tau_n}')
    print(f'  {"j":>3}  {"sector":>6}  {"Re(DR/BFK)":>20}  {"Im(DR/BFK)":>20}  {"  |DR/BFK - 1|": >18}')
    print('  ' + '-'*72)
    for j in [1,3,5,7,9]:
        ratio = dr_modes[n_mode][j] / bfk_modes[n_mode][j]
        print(f'  {j:>3}    ω⁻   {float(ratio.real):>20.10f}  {float(ratio.imag):>20.10f}  {float(abs(ratio-1)):>18.6e}')
    for j in [2,4,6,8]:
        ratio = dr_modes[n_mode][j] / bfk_modes[n_mode][j]
        print(f'  {j:>3}    ω⁺   {float(ratio.real):>20.10f}  {float(ratio.imag):>20.10f}  {float(abs(ratio-1)):>18.6e}')

print('Cumulative partial sums Σ_{n=1}^N τ(n)·ω⁻_n  (j=1 estimate, Im part)')
print(f'Reference ω⁻/i = {float(BROWN_OMEGA_MINUS):.6f}')
print()
print(f'{"N":>3}  {"DR partial sum":>22}  {"BFK partial sum":>22}  {"DR-ref":>14}  {"BFK-ref":>14}')
print('-'*80)

cum_dr  = mpc(0)
cum_bfk = mpc(0)

for n_mode in range(1, N_MODES + 1):
    tau_n = int(Delta_series.get(n_mode, 0))
    cum_dr  += tau_n * dr_modes[n_mode][1]
    cum_bfk += tau_n * bfk_modes[n_mode][1]

    dr_val  = float((cum_dr  / mpc(0,1)).real)
    bfk_val = float((cum_bfk / mpc(0,1)).real)
    ref     = float(BROWN_OMEGA_MINUS)
    print(f'{n_mode:>3}  {dr_val:>22.6f}  {bfk_val:>22.6f}  {dr_val-ref:>14.6f}  {bfk_val-ref:>14.6f}')
