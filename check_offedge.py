"""
Capstone check for lem:offedge:  the full meromorphic L-function, pole cleanly off the edge.

Test form (weight 12):  f = Delta / (j - j(tau_p)),  tau_p = 1.5 i   (Im = 1.5 > 1, above edge).
 - simple pole at tau_p (and its T-translates), residue Res = Delta(tau_p)/j'(tau_p);
 - polylog data of a simple pole: r*(1) = -2 pi i Res, Phat f = r*(1) Li_0(e(tau-tau_p));
 - other Gamma-translates of the pole sit BELOW the edge (e.g. -1/tau_p = 0.667 i).

Claim (lem:offedge):
   L*_direct(s) = i^{-s} int_{i-1}^{i} f ktilde dtau
                =  sum_{n>=1} c_{Rf}(n) G(s,n)            [regularized cusp sum, thm:weakL on (I-Phat)f]
                +  r*(1) * I_0(s, tau_p)                  [lem:polylogcf, Im tau_p>1]
with c_{Rf}(n) the Fourier modes of (I-Phat)f on the line Im=1, and
   G(s,m) = Gamma(s,2pi m)/(2pi m)^s + i^k Gamma(k-s,2pi m)/(2pi m)^{k-s},
   I_0(s,z) = (1/s + i^k/(k-s)) - sum_{n>=1} e^{2pi i n z} G(s,-n).
"""
import mpmath as mp
mp.mp.dps = 30

k = 12
TWOPII = 2j * mp.pi
IK = mp.power(1j, k)
CONTOUR = [-1 + 1j, 1j]
TAUP = mp.mpc(0, '1.5')

def e(x):
    return mp.e ** (TWOPII * x)

def Delta(tau):
    q = mp.e ** (TWOPII * tau)
    p = mp.mpf(1)
    for n in range(1, 60):
        p *= (1 - q ** n) ** 24
    return q * p

def jj(tau):
    return mp.kleinj(tau)

def ktilde(tau, s):
    return mp.zeta(1 - s, tau + 1) - mp.e ** (1j * mp.pi * (s - 1)) * mp.zeta(s - (k - 1), tau + 1)

def gpair(s, m):
    a = 2 * mp.pi * m
    return mp.gammainc(s, a) / mp.power(a, s) + IK * mp.gammainc(k - s, a) / mp.power(a, k - s)

JP = jj(TAUP)
def f(tau):
    return Delta(tau) / (jj(tau) - JP)

# principal-part data of the simple pole
JPRIME = mp.diff(jj, TAUP)
RES = Delta(TAUP) / JPRIME
RSTAR1 = -TWOPII * RES

def Phat(tau):
    return RSTAR1 * mp.polylog(0, e(tau - TAUP))

def Rf(tau):                       # (I - Phat) f
    return f(tau) - Phat(tau)

def I0_closed(s, z, nmax=80):
    return (1 / s + IK / (k - s)) - mp.fsum(mp.e ** (TWOPII * n * z) * gpair(s, -n) for n in range(1, nmax + 1))

def Lstar_direct(s):
    return mp.power(1j, -s) * mp.quad(lambda tau: f(tau) * ktilde(tau, s), CONTOUR)

# Fourier modes of (I-Phat)f on Im=1 via DFT of samples Rf(x+i), x in [0,1)
M = 128
SAMPLES = [Rf(mp.mpf(jx) / M + 1j) for jx in range(M)]
def cRf(n):
    dft = mp.fsum(SAMPLES[jx] * mp.e ** (-TWOPII * n * (mp.mpf(jx) / M)) for jx in range(M)) / M
    return dft * mp.e ** (2 * mp.pi * n)     # the e^{-2pi i n (i)} = e^{2pi n} factor

def Lstar_RHS(s, nmax=42):
    cusp = mp.fsum(cRf(n) * gpair(s, n) for n in range(1, nmax + 1))
    interior = RSTAR1 * I0_closed(s, TAUP)
    return cusp, interior, cusp + interior

print(f"tau_p={TAUP},  j(tau_p)={mp.nstr(JP,8)},  r*(1)={mp.nstr(RSTAR1,8)}")
print(f"sanity  c_Rf(0)={mp.nstr(cRf(0),3)}  c_Rf(-1)={mp.nstr(cRf(-1),3)}  (should be ~0)")
for s in [mp.mpf('2.5'), mp.mpc('3.0', '1.3')]:
    ld = Lstar_direct(s)
    cusp, interior, rhs = Lstar_RHS(s)
    print(f"\ns={s}:")
    print(f"  cusp sum   = {mp.nstr(cusp,12)}")
    print(f"  interior   = {mp.nstr(interior,12)}")
    print(f"  L*_RHS     = {mp.nstr(rhs,14)}")
    print(f"  L*_direct  = {mp.nstr(ld,14)}")
    print(f"  |diff|     = {mp.nstr(abs(ld-rhs),3)}")
