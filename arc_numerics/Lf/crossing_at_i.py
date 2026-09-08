"""Does the crossing term X_N of eq:geocross vanish at z = i?

arc_pole_side_test.py finds that at a generic arc pole the two side-choices differ by
X_N, but at z = i BOTH sides agree with the bare eq:geoIN to the dps floor.  That is only
possible if X_N(s,i) = 0, i.e. a block pole at i costs no residue and there is nothing to
choose.  If so it is the S-vs-T cancellation of prop:ellfull, seen at block level: the
tau^{s-1} and ktil contributions to the residue cancel against each other at i, and at i
only.

Printed per (s,N): |X_N| at z = i and at z = e^{1.4 i}, with the two kernel pieces shown
separately so a cancellation is visible rather than inferred.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi                                            # noqa: E402

mp.mp.dps = 40
KK = 12
SVALS = [mp.mpf(3), mp.mpf(7), mp.mpc('2.4', '0.6')]


def prod(terms):
    p = mp.mpf(1)
    for t in terms:
        p *= t
    return p


def pieces(s, z, N):
    """the tau^{s-1} piece and the ktil piece of eq:geocross, before summing"""
    pre = mp.e**(-I * pi * s / 2) / (2 * I * pi)**N
    pw = (-1)**N * z**(s - 1 - N) * prod(s - 1 - j for j in range(N))
    kt = (mp.zeta(1 - s + N, z + 1) * prod(1 - s + j for j in range(N))
          - mp.e**(I * pi * (s - 1)) * mp.zeta(1 - KK + s + N, z + 1)
          * prod(1 - KK + s + j for j in range(N)))
    return pre * pw, pre * kt


print("dps = %d, k = %d" % (mp.mp.dps, KK), flush=True)
for s in SVALS:
    print("\ns = %s" % mp.nstr(s, 8), flush=True)
    for lbl, z in (("z = i        ", I), ("z = e^{1.4i} ", mp.e**(I * mp.mpf('1.4')))):
        for N in (0, 1, 2):
            a, b = pieces(s, z, N)
            tot = a + b
            sc = max(abs(a), abs(b))
            print("  %s N=%d  |pow| %-12s |ktil| %-12s |X_N| %-12s  ratio %s"
                  % (lbl, N, mp.nstr(abs(a), 6), mp.nstr(abs(b), 6),
                     mp.nstr(abs(tot), 6), mp.nstr(abs(tot) / sc, 4)), flush=True)
