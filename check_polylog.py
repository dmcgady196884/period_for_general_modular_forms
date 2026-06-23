"""
Check lem:polylogcf  (closed form of the polylog integral I_N over OUR contour [i-1,i]).

I_N(s,z) := i^{-s} * int_{i-1}^{i} Li_{-N}(e(tau-z)) ktilde_T(tau,s) dtau,   e(x)=exp(2pi i x)
ktilde_T(tau,s) = zeta(1-s, tau+1) - exp(i pi (s-1)) zeta(s-(k-1), tau+1).

Claim (Im z < 1):
  I_N(s,z) = sum_{n>=1} n^N e^{-2pi i n z} [ Gamma(s,2pi n)/(2pi n)^s + i^k Gamma(k-s,2pi n)/(2pi n)^{k-s} ].

Building block (n>0):
  i^{-s} int_{i-1}^{i} e^{2pi i n tau} ktilde_T dtau  ==  Gamma(s,2pi n)/(2pi n)^s + i^k Gamma(k-s,2pi n)/(2pi n)^{k-s}.
"""
import mpmath as mp
mp.mp.dps = 30

k = 12
TWOPII = 2j * mp.pi
IK = mp.power(1j, k)          # i^k  (=1 for k=12)
CONTOUR = [-1 + 1j, 1j]       # the T-segment at tau_0=i, i.e. [i-1, i]

def e(x):
    return mp.e ** (TWOPII * x)

def ktilde(tau, s):
    return mp.zeta(1 - s, tau + 1) - mp.e ** (1j * mp.pi * (s - 1)) * mp.zeta(s - (k - 1), tau + 1)

def gpair(s, n):              # the per-mode Gamma pair (thm:weakL building block)
    a = 2 * mp.pi * n
    return mp.gammainc(s, a) / a ** s + IK * mp.gammainc(k - s, a) / a ** (k - s)

def block_direct(s, n):       # i^{-s} int e^{2pi i n tau} ktilde
    f = lambda tau: mp.e ** (TWOPII * n * tau) * ktilde(tau, s)
    return mp.power(1j, -s) * mp.quad(f, CONTOUR)

def I_N_direct(s, z, N):
    f = lambda tau: mp.polylog(-N, e(tau - z)) * ktilde(tau, s)
    return mp.power(1j, -s) * mp.quad(f, CONTOUR)

def I_N_sum(s, z, N, nmax=400):
    return mp.fsum(n ** N * mp.e ** (-TWOPII * n * z) * gpair(s, n) for n in range(1, nmax + 1))

print(f"k={k},  i^k={mp.chop(IK)},  contour={CONTOUR},  dps={mp.mp.dps}")

for s in [mp.mpf('2.5'), mp.mpc('3.0', '1.3')]:
    print(f"\n=== s = {s} ===")
    print(" building block  i^{-s} int e^{2pi i n tau} ktilde   vs   Gamma-pair:")
    for n in [1, 2, 3, 5, 8]:
        bd, gp = block_direct(s, n), gpair(s, n)
        print(f"   n={n}: |block - pair| = {mp.nstr(abs(bd - gp), 3)}")
    for z in [mp.mpc('0.3', '0.5'), mp.mpc('-0.2', '0.8'), mp.mpc('0.41', '0.05')]:
        print(f"  z = {z}  (Im z = {mp.im(z)} < 1):")
        for N in [0, 1, 2, 3]:
            di = I_N_direct(s, z, N)
            su = I_N_sum(s, z, N)
            print(f"    N={N}: direct={mp.nstr(di,10)}  sum={mp.nstr(su,10)}  |diff|={mp.nstr(abs(di-su),3)}")
