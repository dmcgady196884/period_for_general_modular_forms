#!/usr/bin/env python3
r"""gamma_sign_check.py -- is the magnitude match sign-independent, and why?

For the mode n=+1 (dominant for Delta) and n=-1 (dominant for Dhat), split the summand
  G_{s,k}(n) = A(n) + B(n),
  A(n) = Gamma(s, 2 pi n)/(2 pi n)^s              [first term: small order s]
  B(n) = i^k Gamma(k-s, 2 pi n)/(2 pi n)^{k-s}    [second term: large order k-s]
Test the claim: |B(+1)| ~ |B(-1)| ~ (k-s-1)!/(2pi)^{k-s}, sign-independent, because
Gamma(k-s, +-2pi) ~ (k-s-1)! when k-s >> 2pi (nearly-complete gamma), whereas the
e^{-+2pi} asymmetry lives only in the SUBDOMINANT A(n).
"""
import mpmath as mp
mp.mp.dps = 40; I = mp.j; pi = mp.pi
FL = dict(flush=True)

def A(n, s, k): return mp.gammainc(s, 2*pi*n)/(2*pi*n)**s
def B(n, s, k): return I**k * mp.gammainc(k-s, 2*pi*n)/(2*pi*n)**(k-s)

print("e^{2pi}=%.1f\n" % float(mp.e**(2*pi)), **FL)
for k in [12, 16, 20, 26]:
    for s in [1]:
        fac = mp.factorial(k-s-1)/(2*pi)**(k-s)     # (k-s-1)!/(2pi)^{k-s}
        print("k=%2d s=%d :  (k-s-1)!/(2pi)^{k-s} = %.4e" % (k, s, float(fac)), **FL)
        for n in (1, -1):
            print("   n=%+d :  |A|=%.4e   |B|=%.4e   |B|/fac=%.4f   |A|/|B|=%.3e"
                  % (n, float(abs(A(n, s, k))), float(abs(B(n, s, k))),
                     float(abs(B(n, s, k))/fac), float(abs(A(n, s, k))/abs(B(n, s, k)))), **FL)
        # ratio of full dominant summands |G(-1)|/|G(1)|
        G1 = A(1, s, k) + B(1, s, k); Gm1 = A(-1, s, k) + B(-1, s, k)
        print("   |G(-1)|/|G(+1)| = %.4f   (vs naive e^{4pi}=%.2e)"
              % (float(abs(Gm1)/abs(G1)), float(mp.e**(4*pi))), **FL)
    print("", **FL)

# and: Gamma(a, +2pi) vs Gamma(a,-2pi) magnitude, vs (a-1)!, as a grows
print("a       Gamma(a,+2pi)   |Gamma(a,-2pi)|   (a-1)!", **FL)
for a in [2, 5, 10, 19, 25]:
    print("a=%2d :  %.4e   %.4e   %.4e"
          % (a, float(mp.gammainc(a, 2*pi)), float(abs(mp.gammainc(a, -2*pi))),
             float(mp.factorial(a-1))), **FL)
