"""R(z) = mu_+/mu_- at the class-number-one CM points, EXACTLY, in Q(sqrt(-D)).

mu_pm of eq:mum/eq:mup have rational coefficients and every Heegner point lies in an imaginary
quadratic field, so R(z_d) lies in that same field and is an exact element of Q(sqrt(-D)) -- no
decimals needed.  Arithmetic here is Fraction pairs (p,q) standing for p + q sqrt(-D), with
(p1,q1)(p2,q2) = (p1 p2 - D q1 q2, p1 q2 + p2 q1) and division by the conjugate.

The points, written as A + B sqrt(-D) with D the field generator rather than the discriminant:
    d = 3   (1 + sqrt(-3))/2     d = 4   sqrt(-1)         d = 7   (1 + sqrt(-7))/2
    d = 8   sqrt(-2)             d = 11  (1+sqrt(-11))/2  d = 19  (1+sqrt(-19))/2
    d = 43,67,163  (1 + sqrt(-d))/2
For d = 8 the point is i sqrt2 = sqrt(-2), so D = 2 and not 8; likewise D = 1 at d = 4.

Cross-checked against the mpmath decimals already reported, so a transcription slip in the
coefficient table cannot pass silently.
"""
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I                                                   # noqa: E402

mp.mp.dps = 30
N = 16
MU_M = {0: F(-12000, 43867), 2: F(16, 3), 4: F(-50, 3), 6: F(52, 3),
        10: F(-52, 3), 12: F(50, 3), 14: F(-16, 3), 16: F(12000, 43867)}
MU_P = {1: F(-16), 3: F(308, 3), 5: F(-182), 7: F(286, 3),
        9: F(286, 3), 11: F(-182), 13: F(308, 3), 15: F(-16)}


class Q:
    """p + q sqrt(-D) with p, q rational."""

    def __init__(self, p, q, D):
        self.p, self.q, self.D = F(p), F(q), D

    def __mul__(self, o):
        return Q(self.p * o.p - self.D * self.q * o.q,
                 self.p * o.q + o.p * self.q, self.D)

    def __add__(self, o):
        return Q(self.p + o.p, self.q + o.q, self.D)

    def scale(self, c):
        return Q(self.p * c, self.q * c, self.D)

    def inv(self):
        n = self.p * self.p + self.D * self.q * self.q
        return Q(self.p / n, -self.q / n, self.D)

    def __str__(self):
        if self.q == 0:
            return str(self.p)
        if self.p == 0:
            return "(%s) sqrt(-%d)" % (self.q, self.D)
        return "%s + (%s) sqrt(-%d)" % (self.p, self.q, self.D)

    def num(self):
        return mp.mpf(self.p.numerator) / self.p.denominator + \
            I * mp.sqrt(self.D) * (mp.mpf(self.q.numerator) / self.q.denominator)


def ev(c, z, D):
    tot = Q(0, 0, D)
    pw = {0: Q(1, 0, D)}
    for m in range(1, max(c) + 1):
        pw[m] = pw[m - 1] * z
    for m, v in c.items():
        tot = tot + pw[m].scale(v)
    return tot


PTS = [("d=3   (1+sqrt-3)/2", F(1, 2), F(1, 2), 3),
       ("d=4   sqrt-1", F(0), F(1), 1),
       ("d=7   (1+sqrt-7)/2", F(1, 2), F(1, 2), 7),
       ("d=8   sqrt-2", F(0), F(1), 2),
       ("d=11  (1+sqrt-11)/2", F(1, 2), F(1, 2), 11),
       ("d=19  (1+sqrt-19)/2", F(1, 2), F(1, 2), 19),
       ("d=43  (1+sqrt-43)/2", F(1, 2), F(1, 2), 43),
       ("d=67  (1+sqrt-67)/2", F(1, 2), F(1, 2), 67),
       ("d=163 (1+sqrt-163)/2", F(1, 2), F(1, 2), 163)]

print("R(z) = mu_+(z)/mu_-(z) at the class-number-one points, exactly\n", flush=True)
for lbl, A, B, D in PTS:
    z = Q(A, B, D)
    dn = ev(MU_M, z, D)
    nu = ev(MU_P, z, D)
    if dn.p == 0 and dn.q == 0:
        print("%-22s mu_- = 0   (mu_+ = %s)" % (lbl, nu), flush=True)
        continue
    R = nu * dn.inv()
    chk = abs(R.num() - (nu.num() / dn.num()))
    # mu_- carries 2/131601 and mu_+ carries 2/3, so R = (131601/3) P_+/P_- = 43867 P_+/P_-
    # with P_pm the integer polynomials; pulling that universal factor out shrinks the height.
    S = R.scale(F(1, 43867))
    print("%-22s R = %s" % (lbl, R), flush=True)
    print("%-22s   = 43867 * [ %s ]" % ("", S), flush=True)
    print("%-22s   numerically %s   (field vs direct: %s)"
          % ("", mp.nstr(R.num(), 18), mp.nstr(chk, 4)), flush=True)
