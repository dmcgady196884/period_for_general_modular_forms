"""Rank check for the converse of lem:arcon.

The defect when gamma^T parts from gamma^S at poles p_1..p_m is

    -(2 pi i)^{n+2} sum_j nu_j a_j (X - p_j Y)^n  -  2 pi i (sum_j X_T(p_j) a_j) D_E ,
    D_E := r_{E_k}|(1+U+U^2)  != 0 .

Vandermonde makes {(X - p_j Y)^n} independent for distinct p_j with m <= n+1.  What is NOT
automatic is that D_E avoids their span: if it did not, the two pieces could cancel and the
defect could vanish with the segments still parted.

So: for x scanning (0,1728), with z, Sz the two arc preimages, compute the relative
distance from D_E to
    span{(X - zY)^n}                     (m = 1, only z flipped)
    span{(X - Sz Y)^n}                   (m = 1, only Sz flipped)
    span{(X - zY)^n, (X - Sz Y)^n}       (m = 2, both flipped)
If these stay bounded away from 0 over the whole range, the converse holds without proviso
for k = 12; if some x drives one to 0, the proviso is real and must be stated.
"""
import os, sys
sys.path.insert(0, "arc_numerics")
from common import mp, I, pi, jay, E12, rvec, NU, inW
mp.mp.dps = 25
KK, NN = 12, 10

rE = rvec(E12, KK)
DE = NU(rE, KK)
print("|r_E12| =", mp.nstr(max(abs(x) for x in rE), 9))
print("|D_E|   =", mp.nstr(max(abs(x) for x in DE), 9))
print("D_E normalised by its largest entry:")
big = max(range(NN+1), key=lambda i: abs(DE[i]))
print("  ", ", ".join(mp.nstr(DE[l]/DE[big], 6) for l in range(NN+1)))
print()

def powvec(p):
    return [mp.binomial(NN, l)*(-p)**l for l in range(NN+1)]

def dist_to_span(v, basis):
    """relative distance from v to span(basis), by SVD projection"""
    m = len(basis)
    A = mp.matrix(NN+1, m)
    for j,b in enumerate(basis):
        for r in range(NN+1): A[r,j] = b[r]
    U_, sv, Vh = mp.svd_c(A, full_matrices=False)
    proj = [mp.mpc(0)]*(NN+1)
    top = max(abs(sv[i]) for i in range(len(sv)))
    for i in range(len(sv)):
        if abs(sv[i])/top < mp.mpf('1e-20'): continue
        c = sum(mp.conj(U_[r,i])*v[r] for r in range(NN+1))
        for r in range(NN+1): proj[r] += c*U_[r,i]
    return max(abs(a-b) for a,b in zip(v,proj))/max(abs(x) for x in v)

print("     x        arg z     dist(D_E, <z>)  dist(D_E, <Sz>)  dist(D_E, <z,Sz>)")
worst = None
for xv in ('20','100','300','500','800','1200','1500','1700'):
    x = mp.mpf(xv)
    # bracket the root: arg in (pi/2, 2pi/3) for the "upper" preimage
    z = mp.findroot(lambda t: jay(t)-x, mp.e**(I*mp.mpf('1.85')))
    if not (pi/3 < mp.arg(z) < 2*pi/3) or abs(abs(z)-1) > mp.mpf('1e-12'):
        z = mp.findroot(lambda t: jay(t)-x, mp.e**(I*mp.mpf('1.7')))
    sz = -1/z
    vz, vs = powvec(z), powvec(sz)
    d1 = dist_to_span(DE, [vz]); d2 = dist_to_span(DE, [vs]); d3 = dist_to_span(DE, [vz,vs])
    if worst is None or d3 < worst[0]: worst = (d3, xv)
    print("  %-8s %-9s %-15s %-16s %s"
          % (xv, mp.nstr(mp.arg(z),7), mp.nstr(d1,6), mp.nstr(d2,6), mp.nstr(d3,6)), flush=True)
print()
print("worst (smallest) two-pole distance over the scan: %s at x = %s" % (mp.nstr(worst[0],6), worst[1]))
print()
print("control: distance from a genuine member of the family to the span (should be ~0)")
z = mp.findroot(lambda t: jay(t)-mp.mpf(500), mp.e**(I*mp.mpf('1.85')))
print("   dist((X-zY)^n, <z,Sz>) =", mp.nstr(dist_to_span(powvec(z), [powvec(z), powvec(-1/z)]), 4))
