"""
Check the CORRECT regime for lem:polylogcf:  Im z >= 1  (the projected poles).

On the contour Im tau = 1, for Im z >= 1 we have |e(tau-z)| >= 1, so the forward series of
Li_{-N}(e(tau-z)) diverges.  Polylog inversion:
    Li_{-N}(w) = (-1)^{N+1} Li_{-N}(1/w)   (N>=1),     Li_0(w) = -1 - Li_0(1/w).
With 1/w = e(z-tau) (|.|<=1) this gives a CONVERGENT closed form:

    G(s,m) := i^{-s} int_{i-1}^{i} e^{2pi i m tau} ktilde(tau,s) dtau
    N>=1:  I_N(s,z) = (-1)^{N+1} sum_{n>=1} n^N e^{2pi i n z} G(s,-n)
    N=0 :  I_0(s,z) = -G(s,0) - sum_{n>=1} e^{2pi i n z} G(s,-n)

Two claims to verify:
  (A) the NEGATIVE-mode building block has the same closed form, G(s,-n) = gpair(s,-n)
      with the incomplete gammas continued to negative argument;
  (B) I_N_direct(s,z) == inversion sum, for Im z > 1.
"""
import mpmath as mp
mp.mp.dps = 30

k = 12
TWOPII = 2j * mp.pi
CONTOUR = [-1 + 1j, 1j]
IK = mp.power(1j, k)

def e(x):
    return mp.e ** (TWOPII * x)

def ktilde(tau, s):
    return mp.zeta(1 - s, tau + 1) - mp.e ** (1j * mp.pi * (s - 1)) * mp.zeta(s - (k - 1), tau + 1)

def G_direct(s, m):
    f = lambda tau: mp.e ** (TWOPII * m * tau) * ktilde(tau, s)
    return mp.power(1j, -s) * mp.quad(f, CONTOUR)

def gpair(s, m):                    # candidate closed form; m may be negative
    a = 2 * mp.pi * m
    return mp.gammainc(s, a) / mp.power(a, s) + IK * mp.gammainc(k - s, a) / mp.power(a, k - s)

def I_N_direct(s, z, N):
    f = lambda tau: mp.polylog(-N, e(tau - z)) * ktilde(tau, s)
    return mp.power(1j, -s) * mp.quad(f, CONTOUR)

def I_N_inv(s, z, N, nmax=200):     # inversion closed form, using gpair for G(s,-n)
    if N == 0:
        return -G_direct(s, 0) - mp.fsum(mp.e ** (TWOPII * n * z) * gpair(s, -n) for n in range(1, nmax + 1))
    return ((-1) ** (N + 1)) * mp.fsum(n ** N * mp.e ** (TWOPII * n * z) * gpair(s, -n) for n in range(1, nmax + 1))

for s in [mp.mpf('2.5'), mp.mpc('3.0', '1.3')]:
    print(f"\n=== s = {s} ===")
    print(" (A) negative-mode block  G_direct(s,-n)  vs  gpair(s,-n):")
    for n in [1, 2, 3, 5]:
        print(f"   n={n}: |G_direct - gpair| = {mp.nstr(abs(G_direct(s, -n) - gpair(s, -n)), 3)}")
    for z in [mp.mpc('-0.3', '1.4'), mp.mpc('0.25', '2.0'), mp.mpc('0.1', '1.15')]:
        print(f"  (B) z = {z}  (Im z = {mp.im(z)} > 1):")
        for N in [0, 1, 2, 3]:
            di, iv = I_N_direct(s, z, N), I_N_inv(s, z, N)
            print(f"     N={N}: direct={mp.nstr(di,10)}  inversion={mp.nstr(iv,10)}  |diff|={mp.nstr(abs(di-iv),3)}")
