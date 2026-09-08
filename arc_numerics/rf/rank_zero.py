"""The x -> 0 end, with the root found in theta rather than in tau.

j has a TRIPLE zero at rho, so j'(rho) = 0 and Newton on j(tau) = x diverges near there.
But we know |z| = 1, so solve the REAL one-dimensional problem j(e^{i theta}) = x for
theta in (pi/2, 2pi/3] by bisection, which is insensitive to the vanishing derivative.
Then Sz = -1/z as before.
"""
import os, sys
sys.path.insert(0, "arc_numerics")
from common import mp, I, pi, jay, E12, rvec, NU
mp.mp.dps = 40
KK, NN = 12, 10
rE = rvec(E12, KK); DE = NU(rE, KK)

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
    return max(abs(a-b) for a,b in zip(v,proj))/max(abs(x) for x in v)

def theta_of(x):
    lo, hi = pi/2, 2*pi/3           # j real, decreasing from 1728 at i to 0 at rho
    f = lambda th: mp.re(jay(mp.e**(I*th))) - x
    for _ in range(300):
        mid = (lo+hi)/2
        if f(mid) > 0: lo = mid
        else: hi = mid
    return (lo+hi)/2

print("  x            theta        |j(z)-x|/x    ||z|-1|   s2/s1     d(<z>)   d(<Sz>)  d(<z,Sz>)")
for lbl in ('1e-8','1e-6','1e-4','1e-2','1e-1','1'):
    x = mp.mpf(lbl)
    th = theta_of(x)
    z = mp.e**(I*th); sz = -1/z
    err = abs(jay(z)-x)/x
    vz, vs = powvec(z), powvec(sz)
    A = mp.matrix(NN+1,2)
    for j,b in enumerate((vz,vs)):
        for r in range(NN+1): A[r,j] = b[r]
    s = mp.svd_c(A, compute_uv=False)
    print("  %-12s %-12s %-13s %-9s %-9s %-8s %-8s %s"
          % (lbl, mp.nstr(th,8), mp.nstr(err,3), mp.nstr(abs(abs(z)-1),3),
             mp.nstr(abs(s[1]/s[0]),4),
             mp.nstr(dist_to_span(DE,[vz]),5), mp.nstr(dist_to_span(DE,[vs]),5),
             mp.nstr(dist_to_span(DE,[vz,vs]),5)), flush=True)

rho = mp.e**(2*I*pi/3)
print()
print("  exact limit x=0:  d =", mp.nstr(dist_to_span(DE,[powvec(rho),powvec(rho+1)]),12))
print("  1/sqrt(2)      =", mp.nstr(1/mp.sqrt(2),12))
