"""
Compute the explicit numerical polynomials at tau_0 = 0.3 + 1.2i:

  M_j, N_j         (moment integrals — the "inputs")
  p_j              (coboundary polynomial coefficients, via Bernoulli operator)
  (P|_S − P)_m     (the polynomial we SUBTRACT from C_S to get C'_S)
  C'_S coefficient (the result, decomposed in Brown's basis)

For both f = Delta and f = DeltaHat.  Output is LaTeX-ready tabular data.
All arithmetic in mpmath at 30 decimal digits.
"""
from mpmath import mp, mpc, mpf, pi, exp, quad, matrix, lu_solve, nstr
from sympy import Rational
from math import comb

mp.dps = 30
TWOPII = mpc(0, 2) * pi
PREF   = TWOPII ** 11

# ----- q-series builders (exact rationals via sympy) -----
def laurent_mul(a, b, lo, hi):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = ka + kb
            if lo <= k <= hi:
                out[k] = out.get(k, Rational(0)) + va * vb
    return {k: v for k, v in out.items() if v != 0}

def laurent_div(num, den, lo_exp, hi_exp):
    den_min  = min(den.keys())
    den_lead = den[den_min]
    quot = {}
    for target in range(lo_exp + den_min, hi_exp + den_min + 1):
        s = num.get(target, Rational(0))
        for kq in range(lo_exp, target - den_min):
            s -= quot.get(kq, Rational(0)) * den.get(target - kq, Rational(0))
        if s != 0:
            quot[target - den_min] = s / den_lead
    return quot

def divisor_sigma(n, k):
    s = 0; d = 1
    while d*d <= n:
        if n % d == 0:
            s += d**k
            if d != n // d: s += (n // d)**k
        d += 1
    return s

def Delta_over_q_qexp(prec):
    series = {0: Rational(1)}
    for n in range(1, prec + 1):
        factor = {}
        for k in range(25):
            te = n * k
            if te > prec: break
            factor[te] = Rational((-1)**k * comb(24, k))
        series = laurent_mul(series, factor, 0, prec)
    return series

PREC = 70
DoQ = Delta_over_q_qexp(PREC + 2)
Delta_series = {k + 1: v for k, v in DoQ.items() if 1 <= k + 1 <= PREC + 2}

E4 = {0: Rational(1)}
for n in range(1, PREC + 3):
    E4[n] = Rational(240 * divisor_sigma(n, 3))
E4c   = laurent_mul(laurent_mul(E4, E4, 0, PREC + 2), E4, 0, PREC + 2)
j_ser = laurent_div(E4c, Delta_series, -1, PREC)
J_ser = dict(j_ser); J_ser[0] = J_ser.get(0, Rational(0)) - Rational(744)
J_ser = {k: v for k, v in J_ser.items() if v != 0}
J_sq  = laurent_mul(J_ser, J_ser, -2, PREC)

D_Jsq = laurent_mul(Delta_series, J_sq,  -2, PREC)
D_J   = laurent_mul(Delta_series, J_ser, -1, PREC)
a1 = -D_Jsq[0] / D_J[0]; a0 = -(D_Jsq[1] + a1 * D_J[1])
assert a1 == 24 and a0 == -393444

DH = {}
for k, v in D_Jsq.items(): DH[k] = DH.get(k, Rational(0)) + v
for k, v in D_J.items():   DH[k] = DH.get(k, Rational(0)) + a1 * v
for k, v in Delta_series.items(): DH[k] = DH.get(k, Rational(0)) + a0 * v
DeltaHat_series = {k: v for k, v in DH.items() if v != 0}


# ----- Make numerical evaluators -----
def make_evaluator(coeffs):
    terms = sorted([(int(n), mpf(str(c))) for n, c in coeffs.items() if c != 0])
    def f(tau):
        tau = mpc(tau); val = mpc(0)
        for n, c in terms:
            val += c * exp(n * TWOPII * tau)
        return val
    return f

f_Delta    = make_evaluator(dict(Delta_series))
f_DeltaHat = make_evaluator(dict(DeltaHat_series))


# ----- Cocycle / coboundary machinery (notebook-identical) -----
def moments(f_eval, tau_a, tau_b, n=10):
    """Returns [m_j]_{j=0..n} where m_j = int_{tau_a}^{tau_b} f(tau) tau^j dtau."""
    tau_a, tau_b = mpc(tau_a), mpc(tau_b)
    dtau = tau_b - tau_a
    out = []
    for j in range(n + 1):
        def integrand(s, j=j):
            tau = tau_a + s * dtau
            return f_eval(tau) * tau**j * dtau
        out.append(quad(integrand, [0, 1]))
    return out

def cocycle_polynomial_from_moments(moments_list, n=10):
    """C(X,Y) coefficient list [c_j] = (2πi)^{n+1} (-1)^j C(n,j) * moment_j."""
    return [PREF * mpf((-1)**j) * mpf(comb(n, j)) * moments_list[j] for j in range(n+1)]

def slash_T(p):
    n = len(p) - 1
    out = [mpc(0)] * (n + 1)
    for j in range(n + 1):
        for i in range(n - j + 1):
            out[j + i] += p[j] * mpf(comb(n - j, i))
    return out

def slash_S(p):
    n = len(p) - 1
    return [p[n - m] * mpf((-1)**m) for m in range(n + 1)]

def solve_coboundary(C_T, n=10):
    A = matrix(n, n); b = matrix(n, 1)
    for m in range(1, n + 1):
        for j in range(m):
            A[m-1, j] = mpf(comb(n - j, m - j))
        b[m-1, 0] = C_T[m]
    psol = lu_solve(A, b)
    return [psol[j, 0] for j in range(n)] + [mpc(0)]


# ----- Compute everything at tau_0 -----
tau0 = mpc('0.3', '1.2')

def all_polys(f_eval, label):
    """Return dict with M_j, N_j, p_m, (P|_S - P)_m, C'_S coefficients."""
    Mj = moments(f_eval, tau0 - 1, tau0)         # T-moments
    Nj = moments(f_eval, -1/tau0, tau0)          # S-moments
    C_T = cocycle_polynomial_from_moments(Mj)
    C_S = cocycle_polynomial_from_moments(Nj)
    p = solve_coboundary(C_T)
    PS = slash_S(p)
    coboundary = [PS[j] - p[j] for j in range(11)]   # (P|_S - P)_{X^{10-j}Y^j}
    C_prime_S = [C_S[j] - coboundary[j] for j in range(11)]
    return {
        'label': label, 'M': Mj, 'N': Nj,
        'C_T': C_T, 'C_S': C_S, 'p': p,
        'coboundary': coboundary, 'C_prime_S': C_prime_S
    }

print("Computing polynomials for Delta ...")
Delta_data    = all_polys(f_Delta,    'Delta')
print("Computing polynomials for DeltaHat ...")
DeltaHat_data = all_polys(f_DeltaHat, 'DeltaHat')


# ----- Format & print -----
def fmt(z, digits=8):
    """Format mpc as (re ± im i) with `digits` significant figures."""
    r, i = z.real, z.imag
    sr = mp.nstr(r, digits, strip_zeros=False)
    si = mp.nstr(abs(i), digits, strip_zeros=False)
    sign = '-' if i < 0 else '+'
    return f"{sr} {sign} {si}\\,i"


def latex_table(data, name):
    """Tabulate (P|_S - P)_m and C'_S coefficients for one form."""
    print(f"\n% ----- {name} -----")
    print(r"\begin{tabular}{c|l|l}")
    print(r"\toprule")
    print(r"monomial & $(P|_S - P)$ coefficient & $C'_S$ coefficient \\")
    print(r"\midrule")
    monomials = [r'X^{10}', r'X^9 Y', r'X^8 Y^2', r'X^7 Y^3', r'X^6 Y^4',
                 r'X^5 Y^5', r'X^4 Y^6', r'X^3 Y^7', r'X^2 Y^8', r'X Y^9', r'Y^{10}']
    for j, mon in enumerate(monomials):
        cb_str = fmt(data['coboundary'][j])
        cp_str = fmt(data['C_prime_S'][j])
        print(rf"${mon}$ & ${cb_str}$ & ${cp_str}$ \\")
    print(r"\bottomrule")
    print(r"\end{tabular}")


def latex_moments_table(d_data, dh_data):
    """Tabulate M_j and N_j moments for both forms."""
    print()
    print(r"\begin{tabular}{c|ll|ll}")
    print(r"\toprule")
    print(r" & \multicolumn{2}{c|}{$\Delta$} & \multicolumn{2}{c}{$\widehat\Delta$} \\")
    print(r"$j$ & $M_j$ ($T$-moment) & $N_j$ ($S$-moment) & $M_j$ & $N_j$ \\")
    print(r"\midrule")
    for j in range(11):
        print(rf"{j} & ${fmt(d_data['M'][j], 6)}$ & ${fmt(d_data['N'][j], 6)}$ "
              rf"& ${fmt(dh_data['M'][j], 6)}$ & ${fmt(dh_data['N'][j], 6)}$ \\")
    print(r"\bottomrule")
    print(r"\end{tabular}")


# ----- Verification: extract omega/eta via interior monomials -----
print("\n--- Numerical verification ---")
P_minus = [mpf(0)]*11
P_minus[1]=mpf(4); P_minus[3]=-mpf(25); P_minus[5]=mpf(42); P_minus[7]=-mpf(25); P_minus[9]=mpf(4)
P_plus  = [mpf(0)]*11
P_plus[2]=mpf(1); P_plus[4]=-mpf(3); P_plus[6]=mpf(3); P_plus[8]=-mpf(1)

for data, sym in [(Delta_data, ('omega', '\\omega')), (DeltaHat_data, ('eta', '\\eta'))]:
    print(f"\n{data['label']}:")
    print(f"  {sym[0]}^- estimates (odd interior):")
    for j in [1, 3, 5, 7, 9]:
        val = data['C_prime_S'][j] / P_minus[j]
        print(f"    j={j}:  {nstr(val, 16)}")
    print(f"  {sym[0]}^+ estimates (even interior):")
    for j in [2, 4, 6, 8]:
        val = data['C_prime_S'][j] / P_plus[j]
        print(f"    j={j}:  {nstr(val, 16)}")

print()
print(r"% =========================================================")
print(r"% LaTeX tables for appendix")
print(r"% =========================================================")
latex_moments_table(Delta_data, DeltaHat_data)
latex_table(Delta_data,    r"\Delta")
latex_table(DeltaHat_data, r"\widehat\Delta")
