"""
Numerical verification of Theorem~\\ref{thm:periods} (period extraction at
dim S_k = 1) for the active paper dr_b_periods_and_Lfunctions.tex.

Six claims verified at 50-digit mpmath precision:

  (i)  Period extraction for f = Delta (weight 12) via the finite-contour
       two-segment cocycle reproduces Brown's published omega^pm(Delta)
       to ~10^-40 relative residual.
  (ii) Same for f = hat-Delta:  extraction reproduces Brown's eta^pm(hat-Delta).
  (iii) Lockstep across all 9 interior monomials j = 1..9 in V_{10}:  each
       gives the same omega^pm (or eta^pm), parity-selected.
  (iv) Basepoint independence between tau_0 = 0.3 + 1.2 i and
       tau_0 = -0.2 + 0.9 i.
  (v)  The Bernoulli closed-form for P (Proposition prop:Bernoulli-inv in
       the active paper; eq:P-bernoulli in analytic_periods_weight12.tex)
       agrees with the LU-decomposition recursion on every p_m (m = 0..9),
       for both Delta and hat-Delta.
  (vi) Clean form of Theorem 1.1 (eq:periods-statement):
            omega^pm(f) = (2 pi i)^{k-1} * (-1)^{j-1} * binom(k-2, j-1) / c_S^{j-1}
                          * [ int_S f(tau) tau^{j-1} dtau
                            + int_T f(tau) ktilde_T(tau, j) dtau ]
       with ktilde_T(tau, s) the Hurwitz-zeta kernel of Theorem 1.2 evaluated at
       integer s = j.  Bare Mellin tau^{j-1} and bare Hurwitz ktilde_T(tau, j)
       inside the bracket; all combinatorial structure outside.  Verified at
       j = 2 (omega^-) and j = 3 (omega^+) for f = Delta.

Bernoulli closed-form (operator identity delta_T = e^{Y d_X} - 1):

    P  =  (1/Y) int_0^X C_T(t, Y) dt
        + sum_{k>=1} (B_k/k!) Y^{k-1} (d/dX)^{k-1} C_T(X, Y),

equivalently in coefficients (p_m = [P]_{X^{10-m} Y^m}, q_r = [C_T]_{X^{10-r} Y^r}):

    p_m  =  q_{m+1}/(10-m)
          + sum_{k=1}^{m+1} (B_k/k!) * q_{m-k+1} * (9-m+k)! / (10-m)!

Brown's published targets (Brown 2017, 1710.07912 sec 8.3):
    Delta:     omega+ = -68916772.80959519475431012465533103043907
               omega- = -5585015.379310401866877139263796275129635 * i
    hat-Delta: eta+   =  127202100647.1770947773171612986108774951
               eta-   =   10276732343.64913275081719307240092090893 * i

Run:  python verify_periods.py
"""

import mpmath as mp
from sympy import Rational
from math import comb

mp.mp.dps = 50

TWOPII = mp.mpc(0, 2) * mp.pi
PREF = TWOPII ** 11
N = 10  # weight k = 12, V_{k-2} = V_10

# -------- 1. Build Delta and hat-Delta q-series (exact rationals) --------

def laurent_mul(a, b, lo, hi):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = ka + kb
            if lo <= k <= hi:
                out[k] = out.get(k, Rational(0)) + va * vb
    return {k: v for k, v in out.items() if v != 0}


def laurent_div(num, den, lo_exp, hi_exp):
    den_min = min(den.keys())
    den_lead = den[den_min]
    quot = {}
    for target in range(lo_exp + den_min, hi_exp + den_min + 1):
        s = num.get(target, Rational(0))
        for k_quot in range(lo_exp, target - den_min):
            s -= quot.get(k_quot, Rational(0)) * den.get(target - k_quot, Rational(0))
        if s != 0:
            quot[target - den_min] = s / den_lead
    return quot


def divisor_sigma(nn, kk):
    s, d = 0, 1
    while d * d <= nn:
        if nn % d == 0:
            s += d ** kk
            if d != nn // d:
                s += (nn // d) ** kk
        d += 1
    return s


def eta24_qexp(prec):
    series = {0: Rational(1)}
    for nn in range(1, prec + 1):
        factor = {}
        for kk in range(25):
            te = nn * kk
            if te > prec:
                break
            factor[te] = Rational((-1) ** kk * comb(24, kk))
        series = laurent_mul(series, factor, 0, prec)
    return series


PREC_Q = 60
eta24 = eta24_qexp(PREC_Q + 2)
Delta_series = {k + 1: v for k, v in eta24.items() if 1 <= k + 1 <= PREC_Q + 2}

E4 = {0: Rational(1)}
for nn in range(1, PREC_Q + 3):
    E4[nn] = Rational(240 * divisor_sigma(nn, 3))
E4_cubed = laurent_mul(laurent_mul(E4, E4, 0, PREC_Q + 2), E4, 0, PREC_Q + 2)
j_series = laurent_div(E4_cubed, Delta_series, -1, PREC_Q)
J_series = dict(j_series)
J_series[0] = J_series.get(0, Rational(0)) - 744
J_series = {k: v for k, v in J_series.items() if v != 0}

J_sq = laurent_mul(J_series, J_series, -2, PREC_Q)
D_Jsq = laurent_mul(Delta_series, J_sq, -2, PREC_Q)
D_J = laurent_mul(Delta_series, J_series, -1, PREC_Q)

hat_Delta_series = {}
for k_, v in D_Jsq.items():
    hat_Delta_series[k_] = hat_Delta_series.get(k_, Rational(0)) + v
for k_, v in D_J.items():
    hat_Delta_series[k_] = hat_Delta_series.get(k_, Rational(0)) + 24 * v
for k_, v in Delta_series.items():
    hat_Delta_series[k_] = hat_Delta_series.get(k_, Rational(0)) + Rational(-393444) * v
hat_Delta_series = {k: v for k, v in hat_Delta_series.items() if v != 0}

delta_terms = sorted([(int(n_), mp.mpf(str(c))) for n_, c in Delta_series.items() if c != 0])
hat_terms = sorted([(int(n_), mp.mpf(str(c))) for n_, c in hat_Delta_series.items() if c != 0])


def Delta_eval(tau):
    tau = mp.mpc(tau)
    val = mp.mpc(0)
    for n_, c in delta_terms:
        val += c * mp.exp(n_ * TWOPII * tau)
    return val


def hat_Delta_eval(tau):
    tau = mp.mpc(tau)
    val = mp.mpc(0)
    for n_, c in hat_terms:
        val += c * mp.exp(n_ * TWOPII * tau)
    return val


# -------- 2. Cocycle integration and slash operators --------

def cocycle_poly(f, tau_a, tau_b):
    """C_gamma(X, Y) = (2 pi i)^11 * int_{tau_a}^{tau_b} f(tau) (X - tau Y)^10 dtau,
    returned as a list of 11 coefficients indexed by j (X^{10-j} Y^j)."""
    tau_a = mp.mpc(tau_a); tau_b = mp.mpc(tau_b)
    dtau = tau_b - tau_a
    moments = []
    for j in range(N + 1):
        def integrand(s, j=j):
            tau = tau_a + s * dtau
            return f(tau) * tau ** j * dtau
        moments.append(mp.quad(integrand, [0, 1]))
    return [PREF * mp.mpf(comb(N, j)) * mp.mpf((-1) ** j) * moments[j]
            for j in range(N + 1)]


def slash_T(p):
    """(P|_T)(X, Y) = P(X + Y, Y).  Returns coefficient list."""
    out = [mp.mpc(0)] * (N + 1)
    for j in range(N + 1):
        for i in range(N - j + 1):
            out[j + i] += p[j] * mp.mpf(comb(N - j, i))
    return out


def slash_S(p):
    """(P|_S)(X, Y) = P(-Y, X).  Returns coefficient list."""
    return [p[N - m] * mp.mpf((-1) ** m) for m in range(N + 1)]


# -------- 3. Two ways of solving delta_T P = C_T --------

def solve_coboundary_LU(C_T):
    """Solve delta_T P = C_T by LU decomposition of the strictly-lower-
    triangular matrix M_{m, j} = binom(10 - j, m - j) of
    Lemma lem:solvability-cuspidality(b)."""
    A = mp.matrix(N, N)
    b = mp.matrix(N, 1)
    for m in range(1, N + 1):
        for j in range(m):
            A[m - 1, j] = mp.mpf(comb(N - j, m - j))
        b[m - 1, 0] = C_T[m]
    psol = mp.lu_solve(A, b)
    return [psol[j, 0] for j in range(N)] + [mp.mpc(0)]


def solve_coboundary_bernoulli(C_T):
    """Solve delta_T P = C_T via the operator-identity Bernoulli closed form
    (Proposition prop:Bernoulli-inv in dr_b_periods_and_Lfunctions.tex;
    eq:P-bernoulli in analytic_periods_weight12.tex):

        p_m  =  q_{m+1}/(N - m)
              + sum_{k=1}^{m+1} (B_k/k!) * q_{m-k+1} * (N-1-m+k)! / (N-m)!,

    with q_r = C_T[r] and B_k the k-th Bernoulli number."""
    p = [mp.mpc(0)] * (N + 1)
    for m in range(N):
        val = mp.mpc(0)
        val += C_T[m + 1] / mp.mpf(N - m)  # k = 0 / antiderivative term
        for k in range(1, m + 2):
            idx = m - k + 1
            if 0 <= idx <= N:
                val += (mp.mpf(mp.bernoulli(k)) / mp.mpf(mp.factorial(k))
                        * C_T[idx]
                        * mp.mpf(mp.factorial(N - 1 - m + k))
                        / mp.mpf(mp.factorial(N - m)))
        p[m] = val
    return p


# -------- 4. Brown's basis (1710.07912 §8.2) --------

P_plus = [mp.mpf(0)] * (N + 1)
P_plus[0] = -mp.mpf(36) / mp.mpf(691); P_plus[10] = mp.mpf(36) / mp.mpf(691)
P_plus[2] = mp.mpf(1); P_plus[4] = -mp.mpf(3); P_plus[6] = mp.mpf(3); P_plus[8] = -mp.mpf(1)
P_minus = [mp.mpf(0)] * (N + 1)
P_minus[1] = mp.mpf(4); P_minus[3] = -mp.mpf(25); P_minus[5] = mp.mpf(42)
P_minus[7] = -mp.mpf(25); P_minus[9] = mp.mpf(4)

INTERIOR_PLUS = [2, 4, 6, 8]
INTERIOR_MINUS = [1, 3, 5, 7, 9]


# -------- 5. Full extraction pipeline --------

def extract_periods(f, tau0):
    C_S = cocycle_poly(f, -1 / tau0, tau0)
    C_T = cocycle_poly(f, tau0 - 1, tau0)
    cuspidality_residual = abs(C_T[0])

    P_LU = solve_coboundary_LU(C_T)
    P_ber = solve_coboundary_bernoulli(C_T)
    bernoulli_max_res = max(abs(P_LU[m] - P_ber[m]) for m in range(N + 1))

    PS = slash_S(P_LU)
    r_f = [C_S[j] - (PS[j] - P_LU[j]) for j in range(N + 1)]

    plus_vals = [r_f[j] / P_plus[j] for j in INTERIOR_PLUS]
    minus_vals = [r_f[j] / P_minus[j] for j in INTERIOR_MINUS]

    return {
        'cuspidality': cuspidality_residual,
        'bernoulli_max_res': bernoulli_max_res,
        'plus_vals': plus_vals,
        'minus_vals': minus_vals,
        'plus': plus_vals[0],
        'minus': minus_vals[0],
        'lockstep_plus': max(abs(plus_vals[i] - plus_vals[0]) for i in range(len(plus_vals))),
        'lockstep_minus': max(abs(minus_vals[i] - minus_vals[0]) for i in range(len(minus_vals))),
    }


# -------- 6. Brown's published targets --------

BROWN_OMEGA_PLUS = mp.mpf('-68916772.80959519475431012465533103043907')
BROWN_OMEGA_MINUS = mp.mpc(0, 1) * mp.mpf('-5585015.379310401866877139263796275129635')
BROWN_ETA_PLUS = mp.mpf('127202100647.1770947773171612986108774951')
BROWN_ETA_MINUS = mp.mpc(0, 1) * mp.mpf('10276732343.64913275081719307240092090893')


# -------- 7. Run and report --------

def report(label, f, tau0, target_plus, target_minus, mark):
    print(f"\n  --- {label} at tau_0 = {tau0} ---")
    res = extract_periods(f, tau0)
    print(f"  cuspidality [C_T]_X^10:     {float(res['cuspidality']):.3e}")
    print(f"  Bernoulli vs LU residual:   {float(res['bernoulli_max_res']):.3e}")
    print(f"  {mark}+ extracted:           {res['plus']}")
    print(f"  {mark}- extracted:           {res['minus']}")
    rel_plus = abs(res['plus'] - target_plus) / abs(target_plus)
    rel_minus = abs(res['minus'] - target_minus) / abs(target_minus)
    print(f"  |{mark}+ - Brown|/|Brown|:   {float(rel_plus):.3e}")
    print(f"  |{mark}- - Brown|/|Brown|:   {float(rel_minus):.3e}")
    print(f"  lockstep spread ({mark}+):  {float(res['lockstep_plus']):.3e}")
    print(f"  lockstep spread ({mark}-):  {float(res['lockstep_minus']):.3e}")
    return res


if __name__ == '__main__':
    print(f"Working at {mp.mp.dps}-digit mpmath precision.")
    print(f"V_{N}-cocycle prefactor: (2 pi i)^{N+1}.")

    print()
    print("=" * 76)
    print(" Theorem 1.1 verification, weight k = 12, dim S_k = 1.")
    print(" Pipeline:  cocycle integration  ->  delta_T P = C_T  ->  r_f = C_S - dP_S")
    print("            ->  read omega^pm (eta^pm) at 9 interior monomials.")
    print("=" * 76)

    tau0_main = mp.mpc('0.3', '1.2')
    r_D_main = report("Delta", Delta_eval, tau0_main, BROWN_OMEGA_PLUS, BROWN_OMEGA_MINUS, "omega")
    r_H_main = report("hat-Delta", hat_Delta_eval, tau0_main, BROWN_ETA_PLUS, BROWN_ETA_MINUS, "eta")

    print()
    print("=" * 76)
    print(" Basepoint independence: re-extract eta^pm(hat-Delta) at second tau_0.")
    print("=" * 76)
    tau0_alt = mp.mpc('-0.2', '0.9')
    r_H_alt = report("hat-Delta", hat_Delta_eval, tau0_alt, BROWN_ETA_PLUS, BROWN_ETA_MINUS, "eta")
    diff_plus = abs(r_H_main['plus'] - r_H_alt['plus'])
    diff_minus = abs(r_H_main['minus'] - r_H_alt['minus'])
    print(f"\n  |eta^+(tau_0_main) - eta^+(tau_0_alt)|:  {float(diff_plus):.3e}")
    print(f"  |eta^-(tau_0_main) - eta^-(tau_0_alt)|:  {float(diff_minus):.3e}")

    # --- (vi) Clean Hurwitz-kernel identity (Theorem 1.1, eq:periods-statement) ---
    print()
    print("=" * 76)
    print(" Theorem 1.1 clean-form identity: bare Mellin + bare Hurwitz inside,")
    print(" all combinatorial structure outside (eq:periods-statement).")
    print("=" * 76)

    def ktilde_T(tau, s, k=12):
        z1 = mp.zeta(1 - s, tau + 1)
        z2 = mp.zeta(s - (k - 1), tau + 1)
        phase = mp.exp(mp.mpc(0, 1) * mp.pi * (s - 1))
        return z1 - phase * z2

    def hurwitz_extract(f, tau0, j, c_S_val):
        """Apply eq:periods-statement clean form at integer s = j."""
        ta_S, tb_S = -1 / tau0, tau0
        ta_T, tb_T = tau0 - 1, tau0
        dS = tb_S - ta_S; dT = tb_T - ta_T
        def f_tau_pow(t):
            tau = ta_S + t * dS
            return f(tau) * tau ** (j - 1) * dS
        def f_hurwitz(t):
            tau = ta_T + t * dT
            return f(tau) * ktilde_T(tau, j) * dT
        I_S = mp.quad(f_tau_pow, [0, 1])
        I_T = mp.quad(f_hurwitz, [0, 1])
        bracket = I_S + I_T
        prefactor = PREF * mp.mpf((-1) ** (j - 1)) * mp.mpf(comb(N, j - 1)) / mp.mpf(c_S_val)
        return prefactor * bracket

    # (j, c_S^{j-1}, parity, Brown target):  j -> V_n monomial index j-1
    hurwitz_cases = [
        (2, 4,  "omega^-", BROWN_OMEGA_MINUS),
        (3, 1,  "omega^+", BROWN_OMEGA_PLUS),
        (4, -25, "omega^-", BROWN_OMEGA_MINUS),
        (5, -3, "omega^+", BROWN_OMEGA_PLUS),
    ]
    for jval, c_S_val, parity, target_use in hurwitz_cases:
        computed = hurwitz_extract(Delta_eval, tau0_main, jval, c_S_val)
        rel = abs(computed - target_use) / abs(target_use)
        print(f"  Delta, j = {jval}, parity = {parity}, c_S^{{j-1}} = {c_S_val}:")
        print(f"    computed: {computed}")
        print(f"    Brown:    {target_use}")
        print(f"    |diff|/|Brown|: {float(rel):.3e}")
    print()
    print("=" * 76)
    print(" Summary")
    print("=" * 76)
    rel_omega_plus = abs(r_D_main['plus'] - BROWN_OMEGA_PLUS) / abs(BROWN_OMEGA_PLUS)
    rel_omega_minus = abs(r_D_main['minus'] - BROWN_OMEGA_MINUS) / abs(BROWN_OMEGA_MINUS)
    rel_eta_plus = abs(r_H_main['plus'] - BROWN_ETA_PLUS) / abs(BROWN_ETA_PLUS)
    rel_eta_minus = abs(r_H_main['minus'] - BROWN_ETA_MINUS) / abs(BROWN_ETA_MINUS)
    print(f" (i)   omega^+(Delta)  matches Brown to:   {float(rel_omega_plus):.3e}")
    print(f"       omega^-(Delta)  matches Brown to:   {float(rel_omega_minus):.3e}")
    print(f" (ii)  eta^+(hat-Delta) matches Brown to:  {float(rel_eta_plus):.3e}")
    print(f"       eta^-(hat-Delta) matches Brown to:  {float(rel_eta_minus):.3e}")
    print(f" (iii) lockstep across 9 interior monomials (worst case):")
    print(f"         Delta:     {float(max(r_D_main['lockstep_plus'], r_D_main['lockstep_minus'])):.3e}")
    print(f"         hat-Delta: {float(max(r_H_main['lockstep_plus'], r_H_main['lockstep_minus'])):.3e}")
    print(f" (iv)  basepoint independence (eta^pm at two tau_0):")
    print(f"         max diff:  {float(max(diff_plus, diff_minus)):.3e}")
    print(f" (v)   Bernoulli closed-form vs LU recursion (worst |p_m diff|):")
    print(f"         Delta:     {float(r_D_main['bernoulli_max_res']):.3e}")
    print(f"         hat-Delta: {float(r_H_main['bernoulli_max_res']):.3e}")
