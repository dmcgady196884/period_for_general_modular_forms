"""
Check the wall residue (lem:cross / the Im z >= 1 half of lem:polylogcf).

I_N(s,z) = i^{-s} int_{i-1}^{i} Li_{-N}(e(tau-z)) ktilde(tau,s) dtau  is analytic in z off the
contour line Im tau = 1, with a jump across it.  CLAIM:

  lim_{eps->0} [ I_N(s, z0+i eps) - I_N(s, z0-i eps) ]  =  2 pi i * Res_{tau=z0} [ i^{-s} Li_{-N}(e(tau-z0)) ktilde(tau,s) ]

for z0 = x0 + i on the contour (x0 in (-1,0)).  The residue is computed independently by a
small circle integral around z0.  If jump/(2 pi i Res) -> a constant +-1 as eps->0, the wall
residue is confirmed.
"""
import mpmath as mp
mp.mp.dps = 40

k = 12
TWOPII = 2j * mp.pi
CONTOUR = [-1 + 1j, 1j]        # [i-1, i]

def e(x):
    return mp.e ** (TWOPII * x)

def ktilde(tau, s):
    return mp.zeta(1 - s, tau + 1) - mp.e ** (1j * mp.pi * (s - 1)) * mp.zeta(s - (k - 1), tau + 1)

def I_N_direct(s, z, N):
    f = lambda tau: mp.polylog(-N, e(tau - z)) * ktilde(tau, s)
    return mp.power(1j, -s) * mp.quad(f, CONTOUR)

def twopii_Res(s, z0, N, r=mp.mpf('0.05')):
    # circle integral  oint i^{-s} Li_{-N}(e(tau-z0)) ktilde(tau,s) dtau = 2 pi i Res
    def f(th):
        tau = z0 + r * mp.e ** (1j * th)
        return (mp.power(1j, -s) * mp.polylog(-N, e(r * mp.e ** (1j * th))) * ktilde(tau, s)
                * (1j * r * mp.e ** (1j * th)))
    return mp.quad(f, [0, mp.pi, 2 * mp.pi])

for s in [mp.mpf('2.5'), mp.mpc('3.0', '1.3')]:
    z0 = mp.mpc('-0.3', '1.0')      # on the contour [i-1,i]
    print(f"\n=== s={s},  z0={z0} ===")
    for N in [0, 1, 2, 3]:
        tpr = twopii_Res(s, z0, N)
        print(f"  N={N}:  2pi i Res = {mp.nstr(tpr, 10)}")
        for eps in [mp.mpf('0.12'), mp.mpf('0.06'), mp.mpf('0.03')]:
            jump = I_N_direct(s, z0 + 1j * eps, N) - I_N_direct(s, z0 - 1j * eps, N)
            print(f"      eps={float(eps):.3f}: jump={mp.nstr(jump,10)}   jump/(2pi i Res)={mp.nstr(jump/tpr,8)}")
