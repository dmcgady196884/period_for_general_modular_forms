"""Is the W-direction of hat r_{f_z} constant in z, or is that a quadrature floor?

cm_vs_generic_k18.py returned coordinate ratios agreeing to ~3e-6 across eight z -- CM,
control and generic alike -- with no monotone dependence on z and with two points at equal
Im z differing by as much as the full spread.  Two readings:

  (a) the direction in P(W) is genuinely z-independent and 3e-6 is the error floor, or
  (b) the direction moves with z and 3e-6 is the signal.

The nullspace basis is orthonormal (common.nullspace returns right singular vectors), so
`coords` is perfectly conditioned and cannot be blamed.  What settles it is the SELF-
consistency of ONE z: recompute its ratios at two mesh depths and at two precisions.  If a
single z moves by the same ~2e-4 absolute as the spread between different z, reading (a)
holds and the CM test measured nothing.  If a single z is stable to 1e-20, the spread is real.

PART 1 is free: the printed coordinates from logs/cm18.txt, tested for rank one and for
Phi = 2 pi i g(z).  PART 2 is the recompute, at z = 1.3i (real, no pole subtleties, above
the arc so plain arcint applies).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (mp, I, pi, E4, E6, Delta, jay, ktil, arcint,           # noqa: E402
                    NS, NU, opmat, nullspace, inW)

K, N = 18, 16

# ---------------------------------------------------------------- PART 1 (printed data)
# (label, Im z, g(z), Phi, c0, c1, c2) transcribed from logs/cm18.txt
ROWS = [
    ("d=7",   1.32288, -0.00023250678,      -0.0014608832j,
     -4.111937873e+13 - 2.924308739e+13j, -2.924356958e+13 + 4.111870072e+13j,
     5.603817921e+11j),
    ("d=8",   1.41421,  0.00014246937,       0.00089516148j,
     -1.297501082e+13 - 9.227538141e+12j, -9.227659648e+12 + 1.297483997e+13j,
     1.768260047e+11j),
    ("d=11",  1.65831, -2.9652784e-5,       -0.00018631394j,
     -6.053634848e+11 - 4.305203744e+11j, -4.305266698e+11 + 6.053546328e+11j,
     8250005107.0j),
    ("d=19",  2.17945, -1.1290558e-6,       -7.0940667e-6j,
     -866410352.9 - 616170964.3j, -616179821.3 + 866397899.0j, 11807601.82j),
    ("1.3i",  1.3,      0.00030083431,       0.0018901977j,
     -5.438713663e+13 - 3.867896642e+13j, -3.867942716e+13 + 5.438648879e+13j,
     7.411992199e+11j),
    (".5+1.4i", 1.4,   -0.00014631361,      -0.00091931554j,
     -1.55788669e+13 - 1.107931471e+13j, -1.107948836e+13 + 1.557862273e+13j,
     2.12311524e+11j),
    (".2+1.3i", 1.3,    7.3632604e-5 + 0.00027991918j, -0.0017587841 + 0.0004626473j,
     6.692557152e+13 - 9.446932464e+11j, -9.437796913e+11 - 6.692543297e+13j,
     -4.3932336e+11 - 5.996388018e+11j),
    (".35+1.5i", 1.5,  -4.7870349e-5 + 6.3948905e-5j, -0.00040180282 - 0.00030077827j,
     -1.621309781e+12 + 5.184916551e+12j, 5.184871246e+12 + 1.621376107e+12j,
     5.736248295e+10 - 1.86993276e+10j),
]

mp.mp.dps = 20
print("PART 1   Phi(f_z) =? 2 pi i g(z)", flush=True)
for lbl, _, gz, Ph, _, _, _ in ROWS:
    pred = 2 * pi * I * mp.mpc(gz)
    print("   %-9s pred %-34s got %-34s rel %s"
          % (lbl, mp.nstr(pred, 9), mp.nstr(mp.mpc(Ph), 9),
             mp.nstr(abs(pred - mp.mpc(Ph)) / abs(pred), 3)), flush=True)

print("\nPART 1   rank of the 8 x 3 coordinate matrix", flush=True)
M = mp.matrix(len(ROWS), 3)
for i, r in enumerate(ROWS):                       # unit-normalise each row first
    v = [mp.mpc(r[4]), mp.mpc(r[5]), mp.mpc(r[6])]
    nv = mp.sqrt(sum(abs(x)**2 for x in v))
    for j in range(3):
        M[i, j] = v[j] / nv
sv = mp.svd_c(M, compute_uv=False)
print("   singular values: %s" % "  ".join(mp.nstr(sv[i], 8) for i in range(3)), flush=True)
print("   sigma2/sigma1 = %s   sigma3/sigma1 = %s"
      % (mp.nstr(sv[1] / sv[0], 6), mp.nstr(sv[2] / sv[0], 6)), flush=True)
print("   (rank 1 at the 1e-6 level would read sigma2/sigma1 ~ 1e-6)", flush=True)

# ---------------------------------------------------------------- PART 2 (recompute)
g = lambda t: Delta(t) * E4(t)
ref = lambda t: E4(t)**3 * E6(t)
z = mp.mpf('1.3') * I


def run(dps, depth):
    mp.mp.dps = dps
    A, B_ = opmat(lambda v: NS(v, K), K), opmat(lambda v: NU(v, K), K)
    MM = mp.matrix(2 * (N + 1), N + 1)
    for r in range(N + 1):
        for c in range(N + 1):
            MM[r, c], MM[r + N + 1, c] = A[r, c], B_[r, c]
    B = nullspace(MM, mp.mpf('1e-18'))
    rE = [x for x in _rvec(ref, depth)]
    x = jay(z)
    f = lambda t: g(t) * (-2 * I * pi * jay(t) * E6(t) / E4(t)) / (jay(t) - x)
    rf = _rvec(f, depth)
    Ph = arcint(f, depth)
    v = [a - Ph * b for a, b in zip(rf, rE)]
    G = mp.matrix(3, 3)
    rhs = mp.matrix(3, 1)
    for a in range(3):
        for b in range(3):
            G[a, b] = sum(mp.conj(B[a][l]) * B[b][l] for l in range(N + 1))
        rhs[a] = sum(mp.conj(B[a][l]) * v[l] for l in range(N + 1))
    c = mp.lu_solve(G, rhs)
    return c[0] / c[2], c[1] / c[2], inW(v, K), Ph


def _rvec(h, depth):
    out = []
    for l in range(N + 1):
        s = mp.mpf(l + 1)
        L = mp.e**(-I * pi * s / 2) * arcint(
            lambda t: h(t) * (t**(s - 1) + ktil(t, s, K)), depth)
        out.append((2 * pi * I)**(N + 1) * (-1)**l * mp.binomial(N, l) * I**(l + 1) * L)
    return out


print("\nPART 2   z = 1.3i recomputed.  Reference row from cm18.txt:", flush=True)
print("   dps 30 depth 10   c0/c2 = -52.184305357752 + 73.377217859608j", flush=True)
print("                     c1/c2 =  73.376343812515 + 52.184926968724j", flush=True)
prev = None
for dps, depth in ((30, 10), (30, 12), (40, 12)):
    r1, r2, w, Ph = run(dps, depth)
    print("   dps %-3d depth %-3d  c0/c2 = %s" % (dps, depth, mp.nstr(r1, 18)), flush=True)
    print("                     c1/c2 = %s" % mp.nstr(r2, 18), flush=True)
    print("                     inW %s / %s    Phi = %s"
          % (mp.nstr(w[0], 4), mp.nstr(w[1], 4), mp.nstr(Ph, 12)), flush=True)
    if prev is not None:
        print("                     MOVED from previous line by %s , %s"
              % (mp.nstr(abs(r1 - prev[0]), 4), mp.nstr(abs(r2 - prev[1]), 4)), flush=True)
    prev = (r1, r2)
print("\n   spread between DIFFERENT z in cm18.txt was 2.2e-4 absolute.", flush=True)
print("   If one z moves that much on its own, the CM comparison measured noise.", flush=True)
