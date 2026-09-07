"""Poles at i by DEFORMATION, the order-2 analogue of rho_limit.py.

DAM's point.  i is an order-2 elliptic point, so j - 1728 has a DOUBLE zero there and a
double pole at i is the merger of an S-orbit PAIR {z, Sz}.  The pair is forced:
z = (1+delta) i has Sz = -1/((1+delta)i) = (1+delta)^{-1} i, just BELOW i.  Poles above
and below i are the SAME orbit, so as delta -> 0 they pinch the arc from both sides and
the arc threads BETWEEN them.  Indenting past i on one side by a fixed h passes around
BOTH members once h > delta -- a different homotopy class, which is what the defect in
pole_at_i_indented.py was measuring.

For j0 != 1728 the poles are off the arc, so thm:geoperiod gives tilde r in W EXACTLY,
and W is closed, so any limit that exists stays in W.  The only question is whether the
limit exists and at what rate.

PREDICTION (analogy with rho: order 3, exponent 2/3):
    residues ~ 1/j'(z), and j' has an (n-1)-fold zero at an order-n elliptic point,
    so at i (n = 2)   tilde r ~ eta^{-1/2},   eta := j0 - 1728,
    invariant v = lim eta^{1/2} tilde r  in W,  branch ambiguity a SIGN (overall scalar).

Family: f_{j0} = E4^5/(j - j0) at k = 20, whose limit form at j0 = 1728 is
E4^5 Delta/E6^2, a genuine double pole at Gamma.i (E6 has a simple zero at i).
Reference form with Phi = 1: E4^5 (holomorphic, constant term 1).

NB eta must avoid the negative reals: j maps the arc onto [0, 1728], so j0 < 1728 has
preimages ON the arc.  eta > 0 puts the poles on the imaginary axis at (1 +- delta) i,
exactly the configuration in question.

Quadrature clusters at theta = pi/2 (the MIDDLE of the arc), not at the endpoints as
common.arcint does -- the poles approach the middle here.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common import mp, I, pi, E4, E6, Delta, jay, ktil, inW, neville

mp.mp.dps = 30
KK, NN = 20, 18
DEPTH = 13

ref = lambda t: E4(t)**5                      # Phi = 1


def nodes_mid(depth):
    """nodes on [pi/3, 2pi/3] clustered geometrically toward pi/2 from both sides"""
    a, b, mid = pi / 3, 2 * pi / 3, pi / 2
    left, right, d = [], [], mid - a
    for _ in range(depth):
        d /= 2
        left.append(mid - d)
        right.append(mid + d)
    return sorted(set([a, mid, b] + left + right))


NODES = nodes_mid(DEPTH)


def arcint(g):
    return -sum(mp.quad(lambda x: g(mp.e**(I * x)) * I * mp.e**(I * x), [u, v])
                for u, v in zip(NODES[:-1], NODES[1:]))


def rvec(g):
    out = []
    for l in range(NN + 1):
        s = mp.mpf(l + 1)
        out.append((2 * pi * I)**(NN + 1) * (-1)**l * mp.binomial(NN, l)
                   * (arcint(lambda t: g(t) * t**(s - 1))
                      + arcint(lambda t: g(t) * ktil(t, s, KK))))
    return out


print("Phi(E4^5) =", mp.nstr(arcint(ref), 12), " (must be +1)", flush=True)
rR = rvec(ref)

print()
print("=" * 74)
print("pole geometry: eta = j0 - 1728 > 0 puts poles at (1 +- delta) i")
print("=" * 74)
for e in ('1e-2', '1e-4', '1e-6'):
    eta = mp.mpf(e)
    zz = mp.findroot(lambda z: jay(z) - (1728 + eta), I * mp.mpf('1.01'))
    print("  eta=%-6s  z = %s   delta = %s   S z = %s"
          % (e, mp.nstr(zz, 10), mp.nstr(abs(zz - I), 6), mp.nstr(-1 / zz, 10)),
          flush=True)

print()
print("=" * 74)
print("does eta^{1/2} tilde r converge, and is tilde r in W for eta != 0?")
print("=" * 74)
store = {}
for psi_lbl, psi in (("0.0", mp.mpf(0)), ("1.1", mp.mpf('1.1'))):
    us, vs = [], []
    print("  ray arg(eta) = %s" % psi_lbl)
    for e in ('1e-2', '1e-3', '1e-4', '1e-5', '1e-6'):
        eta = mp.mpf(e) * mp.e**(I * psi)
        j0 = 1728 + eta
        f = lambda t, j0=j0: E4(t)**5 / (jay(t) - j0)
        Phi = arcint(f)
        tr = [a - Phi * b for a, b in zip(rvec(f), rR)]
        sc = max(abs(x) for x in tr)
        s_, u_ = inW(tr, KK)
        us.append(eta**(mp.mpf(1) / 2))
        vs.append([x / tr[9] for x in tr])
        print("    |eta|=%-5s |Phi|=%-13s |tilde r|=%-15s in W: %s / %s"
              % (e, mp.nstr(abs(Phi), 7), mp.nstr(sc, 8),
                 mp.nstr(s_, 4), mp.nstr(u_, 4)), flush=True)
        store.setdefault(psi_lbl, []).append((e, sc, tr))
    print("    |tilde r| ratios per decade (sqrt(10) = 3.1623 confirms eta^{-1/2}):")
    row = store[psi_lbl]
    print("     ", ", ".join(mp.nstr(row[m + 1][1] / row[m][1], 6)
                             for m in range(len(row) - 1)))
    lim = [neville(us, [vs[m][l] for m in range(len(us))]) for l in range(NN + 1)]
    store[psi_lbl + "_lim"] = [x[0] for x in lim]
    print("    extrapolated direction: max Neville shift = %s"
          % mp.nstr(max(x[1] for x in lim), 4))
    print("    limit in W: %s / %s"
          % tuple(mp.nstr(x, 5) for x in inW(store[psi_lbl + "_lim"], KK)))
    print()

a, b = store["0.0_lim"], store["1.1_lim"]
print("  ray independence of the limiting direction:  max|a - b| = %s"
      % mp.nstr(max(abs(x - y) for x, y in zip(a, b)), 8))
