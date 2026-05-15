"""
Compute the 9 polynomial kernel pairs (k_S^{pm, j}, k_T^{pm, j}) for which

    omega^pm(f)  =  (2 pi i)^11 / [P^pm_S]_{X^{10-j} Y^j}
                   * [ int_{-1/tau0}^{tau0}  f(tau) k_S^{pm,j}(tau)  dtau
                       + int_{tau0-1}^{tau0} f(tau) k_T^{pm,j}(tau) dtau ]

holds for any weight-12 modular form f on SL_2(Z), at any tau0 in H.

The S-kernel is a monomial:  k_S^{pm, j}(tau) = (-1)^j binom(10, j) tau^j.
The T-kernel is a Q-polynomial of degree <= 10 obtained from the Bernoulli
operator that inverts T - 1.

We use sympy throughout for exact rational arithmetic.
"""
from sympy import Rational, binomial, factorial, bernoulli, symbols, Poly, latex


def beta(m, r):
    """
    Coefficient of M_r in  p_m / (2 pi i)^11, where p_m is the X^{10-m}Y^m
    coefficient of the Bernoulli closed-form coboundary polynomial P.

    Formula (m = 0, ..., 9; r = 1, ..., m+1):
        beta(m, r) = (-1)^r binom(10, r) B_{m-r+1} (10-r)!
                       / [ (m-r+1)! (10-m)! ]
    Outside that range, beta(m, r) = 0.

    Cuspidality (M_0 = 0) lets us drop r = 0.
    """
    if r < 1 or r > m + 1:
        return Rational(0)
    k = m - r + 1
    return ((-1)**r * binomial(10, r) * bernoulli(k) * factorial(10 - r)
            / (factorial(k) * factorial(10 - m)))


def beta_full(l, r):
    """Extends beta to l = 0..10.  p_10 = 0 by convention, so beta_full(10, r) = 0."""
    if l == 10:
        return Rational(0)
    return beta(l, r)


# Brown's basis coefficients [P^pm_S]_{X^{10-j} Y^j} at the interior monomials.
P_basis = {
    -1: {1: Rational(4),  3: Rational(-25), 5: Rational(42),
         7: Rational(-25), 9: Rational(4)},          # P^-_S
    +1: {2: Rational(1),  4: Rational(-3),  6: Rational(3),
         8: Rational(-1)}                            # P^+_S
}


def kernel_pair(j):
    """
    For interior j (1 <= j <= 9), return (sign, [P^pm_S]_j, kS_dict, kT_dict).

    Kernels are *universal* (NOT divided by [P^pm_S]_j); the formula is

        omega^pm(f) = (2 pi i)^11 / [P^pm_S]_j
                    * [ int_{-1/tau0}^{tau0}  f(tau) kS(tau) dtau
                      + int_{tau0-1}^{tau0}   f(tau) kT(tau) dtau ].
    """
    sign = -1 if j % 2 == 1 else +1
    P_j  = P_basis[sign][j]

    # S-kernel:  k_S(tau) = (-1)^j binom(10, j) tau^j   (universal in j; no P_j)
    kS = {j: Rational((-1)**j * binomial(10, j), 1)}

    # T-kernel:  k_T(tau) = sum_r [ beta(j, r) - (-1)^j beta_full(10-j, r) ] tau^r
    kT = {}
    for r in range(1, 11):
        coef = beta(j, r) - (-1)**j * beta_full(10 - j, r)
        if coef != 0:
            kT[r] = coef

    return sign, P_j, kS, kT


def fmt_poly(d, var='tau'):
    """Render a dict {power: rational coefficient} as a readable polynomial string."""
    if not d:
        return "0"
    terms = []
    for r in sorted(d):
        c = d[r]
        if c == 0:
            continue
        sign = '-' if c < 0 else '+'
        ac = abs(c)
        if r == 0:
            t = f"{ac}"
        elif ac == 1:
            t = f"{var}^{r}" if r > 1 else var
        else:
            t = f"({ac}) {var}^{r}" if r > 1 else f"({ac}) {var}"
        terms.append((sign, t))
    out = []
    for i, (s, t) in enumerate(terms):
        if i == 0:
            out.append(f"-{t}" if s == '-' else t)
        else:
            out.append(f" {s} {t}")
    return ''.join(out)


def fmt_latex(d, var=r'\tau'):
    """Same but for LaTeX output."""
    if not d:
        return "0"
    terms = []
    for r in sorted(d):
        c = d[r]
        if c == 0:
            continue
        sign = '-' if c < 0 else '+'
        ac = abs(c)
        ac_str = latex(ac)
        if r == 0:
            t = ac_str
        elif ac == 1:
            t = f"{var}^{{{r}}}" if r > 1 else var
        else:
            t = f"{ac_str}\\,{var}^{{{r}}}" if r > 1 else f"{ac_str}\\,{var}"
        terms.append((sign, t))
    out = []
    for i, (s, t) in enumerate(terms):
        if i == 0:
            out.append(f"-{t}" if s == '-' else t)
        else:
            out.append(f" {s} {t}")
    return ''.join(out)


if __name__ == '__main__':
    print("=" * 78)
    print("Polynomial kernels for periods omega^pm in weight 12")
    print("=" * 78)
    print()
    print("Formula:  omega^pm(f) = (2 pi i)^11 / [P^pm_S]_j *")
    print("           [ int_{-1/tau0}^{tau0}  f(tau) k_S(tau) dtau")
    print("           + int_{tau0-1}^{tau0}   f(tau) k_T(tau) dtau ]")
    print()
    print("-" * 78)
    print(f"{'j':>3}  {'sgn':>4}  {'[P_S]_j':>8}  kS(tau)                  kT(tau)")
    print("-" * 78)
    for j in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
        sign, Pj, kS, kT = kernel_pair(j)
        sgn_str = '-' if sign == -1 else '+'
        print(f"{j:>3}  {sgn_str:>4}  {str(Pj):>8}  {fmt_poly(kS):<22}  {fmt_poly(kT)}")
    print()

    print("=" * 78)
    print("LaTeX-ready table rows:")
    print("=" * 78)
    print()
    for j in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
        sign, Pj, kS, kT = kernel_pair(j)
        sgn_str = '-' if sign == -1 else '+'
        print(rf"$j={j}$ ({sgn_str}) & ${latex(Pj)}$ & ${fmt_latex(kS)}$ & ${fmt_latex(kT)}$ \\")
