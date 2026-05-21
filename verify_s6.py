"""
Per-mode verification at s = 6 of DR-B = BFK for hat-Delta = q^{-1} + sum_{n>=2} a(n) q^n.

DR-B formula at tau_0 = i (S-segment collapses to a point):
  Lambda*_DR(f, s) = i^{-s} * integral_{i-1}^{i} f(tau) * kT_tilde(tau, s) dtau
  where kT_tilde(tau, 6) = 2*zeta(-5, tau+1) = -B_6(tau+1)/3.

BFK formula:
  L*_BFK(f, s) = sum_{n != 0} a(n) [Gamma(s, 2*pi*n)/(2*pi*n)^s + Gamma(k-s, 2*pi*n)/(2*pi*n)^(k-s)]
  At s = 6, k = 12: both terms equal Gamma(6, 2*pi*n)/(2*pi*n)^6.

Per-mode claim:
  DR-B mode-n integral:  -i^{-s}/3 * integral_{i-1}^{i} q^n * B_6(tau+1) dtau
  BFK mode-n value:       2 * Gamma(6, 2*pi*n)/(2*pi*n)^6

These should be equal for every n in Z\{0}.

We verify this numerically for n = -1 (polar mode), n = 1, 2, ..., 5 (cuspidal modes).
"""

import mpmath as mp

mp.mp.dps = 50  # 50-digit precision

# Bernoulli polynomial B_6(x) = x^6 - 3x^5 + (5/2)x^4 - (1/2)x^2 + 1/42
def B6(x):
    return x**6 - 3*x**5 + mp.mpf(5)/2 * x**4 - mp.mpf(1)/2 * x**2 + mp.mpf(1)/42

# Kernel at s = 6
def kT_tilde_s6(tau):
    # kT_tilde(tau, 6) = -B_6(tau+1)/3
    return -B6(tau + 1) / 3

# DR-B per-mode integral: integral_{i-1}^{i} q^n * kT_tilde(tau, 6) dtau
# tau = i-1 + t for t in [0, 1] is the horizontal segment
# So dtau = dt, tau goes from i-1 (t=0) to i (t=1)
def DRB_mode(n):
    """Compute the DR-B mode-n integral at s = 6, tau_0 = i."""
    def integrand(t):
        tau = mp.mpc(-1, 1) + t   # tau = (i - 1) + t, real parameter t in [0,1]
        q_n = mp.exp(2 * mp.pi * mp.mpc(0, 1) * n * tau)  # q^n = e^{2 pi i n tau}
        return q_n * kT_tilde_s6(tau)
    integral = mp.quad(integrand, [0, 1])
    # Apply the i^{-s} = i^{-6} = -1 prefactor that sits in front of Theorem 8.1
    return mp.power(mp.mpc(0, 1), -6) * integral

# BFK per-mode value at s = 6 for one n
def BFK_mode(n):
    """Compute the BFK mode-n contribution at s = 6, k = 12."""
    # For n > 0: standard incomplete gamma
    # For n < 0: use the explicit polynomial formula Gamma(6, x) = 5! e^{-x} sum_{m=0}^{5} x^m/m!
    # mpmath handles negative arguments of Gamma(s, x) via analytic continuation
    x = 2 * mp.pi * n
    # Use mpmath's gammainc for upper incomplete gamma
    # mpmath.gammainc(s, a, b) = integral_a^b t^{s-1} e^{-t} dt
    # We want Gamma(s, x) = integral_x^infty t^{s-1} e^{-t} dt = gammainc(s, x, mp.inf)
    # But for negative x, we use the explicit closed form
    if n > 0:
        Gamma_6_x = mp.gammainc(6, x, mp.inf)
        return 2 * Gamma_6_x / x**6
    else:
        # Use Gamma(6, x) = 120 * e^{-x} * sum_{m=0}^5 x^m/m!
        # This is valid by analytic continuation in x
        Gamma_6_x = 120 * mp.exp(-x) * sum(x**m / mp.factorial(m) for m in range(6))
        return 2 * Gamma_6_x / x**6

# Test cases
print("=" * 80)
print(f"Per-mode equality check at s = 6, weight k = 12, tau_0 = i")
print(f"Precision: {mp.mp.dps} digits")
print("=" * 80)
print()

modes_to_test = [-1, 1, 2, 3, 4, 5]
results = []

for n in modes_to_test:
    drb = DRB_mode(n)
    bfk = BFK_mode(n)
    diff = drb - bfk
    rel_err = abs(diff) / max(abs(drb), abs(bfk), mp.mpf(1))
    results.append((n, drb, bfk, rel_err))
    print(f"n = {n:+d}")
    print(f"  DR-B  = {mp.nstr(drb, 30)}")
    print(f"  BFK   = {mp.nstr(bfk, 30)}")
    print(f"  diff  = {mp.nstr(diff, 10)}")
    print(f"  |rel error| = {mp.nstr(rel_err, 5)}")
    print()

# Aggregate
print("=" * 80)
all_match = all(rel_err < mp.mpf(10)**(-40) for _, _, _, rel_err in results)
print(f"All modes match to < 1e-40: {all_match}")
print("=" * 80)
