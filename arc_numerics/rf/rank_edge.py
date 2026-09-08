"""Push the rank check right up to both ends of the cut [0, 1728].

At x -> 1728 the two arc preimages MERGE at i, so v_z -> v_Sz and the two-dimensional
span collapses to one dimension.  A collapsing span makes the distance go UP, not down,
so this end should be safe -- but the two-pole projection becomes numerically rank
deficient there, so we print the singular values of the span matrix to make the collapse
visible rather than silent.

At x -> 0 the preimages go to rho and rho+1, which are DISTINCT (indeed S rho = rho+1),
so the span stays two-dimensional and nothing degenerates.

Reported per x: arg z, |z| (must be 1: the root must lie on the arc), the two singular
values of [v_z, v_Sz], and the relative distances from D_E to the three spans.
"""
import os, sys
sys.path.insert(0, "arc_numerics")
from common import mp, I, pi, jay, E12, rvec, NU
mp.mp.dps = 40
KK, NN = 12, 10

rE = rvec(E12, KK); DE = NU(rE, KK)
print("|D_E| =", mp.nstr(max(abs(x) for x in DE), 9), "\n")

def powvec(p): return [mp.binomial(NN,l)*(-p)**l for l in range(NN+1)]

def dist_to_span(v, basis):
    A = mp.matrix(NN+1, len(basis))
    for j,b in enumerate(basis):
        for r in range(NN+1): A[r,j] = b[r]
    U_, sv, Vh = mp.svd_c(A, full_matrices=False)
    top = max(abs(sv[i]) for i in range(len(sv)))
    proj = [mp.mpc(0)]*(NN+1)
    for i in range(len(sv)):
        if abs(sv[i])/top < mp.mpf('1e-30'): continue
        c = sum(mp.conj(U_[r,i])*v[r] for r in range(NN+1))
        for r in range(NN+1): proj[r] += c*U_[r,i]
    return max(abs(a-b) for a,b in zip(v,proj))/max(abs(x) for x in v), sv

def svals(basis):
    A = mp.matrix(NN+1, len(basis))
    for j,b in enumerate(basis):
        for r in range(NN+1): A[r,j] = b[r]
    return mp.svd_c(A, compute_uv=False)

print("  x                      arg z       ||z|-1|    s2/s1      d(<z>)    d(<Sz>)   d(<z,Sz>)")
for lbl, xv, guess in [
    ("1e-6",        mp.mpf('1e-6'),        mp.e**(I*mp.mpf('2.0944'))),
    ("1e-3",        mp.mpf('1e-3'),        mp.e**(I*mp.mpf('2.0943'))),
    ("1e-1",        mp.mpf('0.1'),         mp.e**(I*mp.mpf('2.093'))),
    ("1",           mp.mpf(1),             mp.e**(I*mp.mpf('2.09'))),
    ("1728-1",      mp.mpf(1728)-1,        mp.e**(I*mp.mpf('1.60'))),
    ("1728-1e-1",   mp.mpf(1728)-mp.mpf('0.1'), mp.e**(I*mp.mpf('1.58'))),
    ("1728-1e-3",   mp.mpf(1728)-mp.mpf('1e-3'), mp.e**(I*mp.mpf('1.5716'))),
    ("1728-1e-6",   mp.mpf(1728)-mp.mpf('1e-6'), mp.e**(I*mp.mpf('1.57081'))),
]:
    try:
        z = mp.findroot(lambda t: jay(t)-xv, guess)
    except Exception as ex:
        print("  %-22s findroot failed: %s" % (lbl, ex)); continue
    sz = -1/z
    vz, vs = powvec(z), powvec(sz)
    s = svals([vz, vs])
    d1,_ = dist_to_span(DE,[vz]); d2,_ = dist_to_span(DE,[vs]); d3,_ = dist_to_span(DE,[vz,vs])
    print("  %-22s %-11s %-10s %-10s %-9s %-9s %s"
          % (lbl, mp.nstr(mp.arg(z),8), mp.nstr(abs(abs(z)-1),3),
             mp.nstr(abs(s[1]/s[0]),4), mp.nstr(d1,5), mp.nstr(d2,5), mp.nstr(d3,5)), flush=True)

print()
print("limiting spans, evaluated exactly at the two ends:")
rho = mp.e**(2*I*pi/3)
d,_ = dist_to_span(DE, [powvec(rho), powvec(rho+1)])
print("  x = 0     : d(D_E, <(X-rho Y)^n, (X-(rho+1)Y)^n>) =", mp.nstr(d,8))
d,_ = dist_to_span(DE, [powvec(I)])
print("  x = 1728  : d(D_E, <(X-iY)^n>)                    =", mp.nstr(d,8))
