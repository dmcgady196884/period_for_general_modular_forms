"""
Option (c):  pole AT tau=i handled by the full two-segment contour at tau_0 = i(1+delta)
            (i.e. Lambda = 1/(1+delta), so 1-Lambda ~ delta small), with the S-segment
            PRINCIPAL-VALUED through i.

  f = Delta/E_6  (weight 6, simple pole at i).
  tau_0 = i(1+delta);   -1/tau_0 = i/(1+delta).
  gamma^S : i/(1+delta) -> i(1+delta)   (vertical, straddles i)  -> PV (symmetric exclusion).
  gamma^T : tau_0-1 -> tau_0            (horizontal at height 1+delta; pole i sits BELOW it).
  L* = i^{-s} [ int_S f tau^{s-1} + int_T f ktilde ].

Watch:  does i^{-s} S -> 0  ("S vanishingly small")?   does L* settle to a finite limit?
"""
import mpmath as mp
mp.mp.dps = 30
k = 6
TWOPII = 2j * mp.pi

def Delta(tau):
    q = mp.e ** (TWOPII * tau)
    p = mp.mpf(1)
    for n in range(1, 60):
        p *= (1 - q ** n) ** 24
    return q * p

def E6(tau):
    q = mp.e ** (TWOPII * tau)
    ss = mp.mpf(0)
    for n in range(1, 90):
        sig5 = sum(d ** 5 for d in range(1, n + 1) if n % d == 0)
        ss += sig5 * q ** n
    return 1 - 504 * ss

def f(tau):
    return Delta(tau) / E6(tau)

def ktilde(tau, s):
    return mp.zeta(1 - s, tau + 1) - mp.e ** (1j * mp.pi * (s - 1)) * mp.zeta(s - (k - 1), tau + 1)

RES = Delta(1j) / mp.diff(E6, 1j)            # Res_i[f] = Delta(i)/E6'(i)
print(f"Res_i[f] = Delta(i)/E6'(i) = {mp.nstr(RES, 10)}")

epsS = mp.mpf('1e-8')                         # symmetric-exclusion radius for the S-PV
for s in [mp.mpf('2.5'), mp.mpc('3.0', '1.3')]:
    print(f"\n=== s = {s} ===")
    pref = mp.power(1j, -s)
    prevL = None
    for delta in [mp.mpf('0.1'), mp.mpf('0.03'), mp.mpf('0.01'), mp.mpf('0.003'), mp.mpf('0.001')]:
        tau0 = 1j * (1 + delta)
        mtau0 = -1 / tau0
        # S-segment PV: straddle i, exclude +-i*epsS
        S = (mp.quad(lambda t: f(t) * mp.power(t, s - 1), [mtau0, 1j - 1j * epsS])
             + mp.quad(lambda t: f(t) * mp.power(t, s - 1), [1j + 1j * epsS, tau0]))
        T = mp.quad(lambda t: f(t) * ktilde(t, s), [tau0 - 1, tau0])
        L = pref * (S + T)
        msg = f"  delta={float(delta):.4f}:  i^-s*S={mp.nstr(pref*S,7)}  i^-s*T={mp.nstr(pref*T,8)}  L*={mp.nstr(L,12)}"
        if prevL is not None:
            msg += f"   dL={mp.nstr(L - prevL, 4)}"
        print(msg)
        prevL = L
