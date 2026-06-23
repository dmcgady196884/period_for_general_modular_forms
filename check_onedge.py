"""
Two checks for the edge case:

  PART A  (lem:onedge, generic on-edge pole):  f = Delta/(j - j(tau_p)), tau_p = -0.5 + i
          (weight 12), pole sitting ON the edge but in the segment interior.  Route the T-path
          through a small semicircle ABOVE vs BELOW tau_p and verify
             L*_below - L*_above = 2 pi i i^{-s} Res_{tau_p}[ f ktilde ].

  PART B  (the bad corner):  f = Delta/E_6  (weight 6), simple pole AT tau = i = the basepoint
          and the order-2 elliptic point.  Probe how the straight T-integral behaves as the
          endpoint approaches i along the edge (expect log divergence).
"""
import mpmath as mp
mp.mp.dps = 30
TWOPII = 2j * mp.pi

def e(x):
    return mp.e ** (TWOPII * x)

def Delta(tau):
    q = mp.e ** (TWOPII * tau)
    p = mp.mpf(1)
    for n in range(1, 60):
        p *= (1 - q ** n) ** 24
    return q * p

def E6(tau):
    q = mp.e ** (TWOPII * tau)
    s = mp.mpf(0)
    for n in range(1, 90):
        sig5 = sum(d ** 5 for d in range(1, n + 1) if n % d == 0)
        s += sig5 * q ** n
    return 1 - 504 * s

def ktilde(tau, s, k):
    return mp.zeta(1 - s, tau + 1) - mp.e ** (1j * mp.pi * (s - 1)) * mp.zeta(s - (k - 1), tau + 1)

# ---------- PART A : generic on-edge jump ----------
print("==== PART A : on-edge jump, tau_p = -0.5 + i, weight 12 ====")
k = 12
taup = mp.mpc('-0.5', '1.0')
JP = mp.kleinj(taup)
def fA(tau):
    return Delta(tau) / (mp.kleinj(tau) - JP)
RESA = Delta(taup) / mp.diff(mp.kleinj, taup)            # Res_{tau_p} f

def Lstar_routed(s, above, d=mp.mpf('0.08')):
    g = lambda t: fA(t) * ktilde(t, s, k)
    seg1 = mp.quad(g, [-1 + 1j, taup - d])
    seg2 = mp.quad(g, [taup + d, 1j])
    if above:
        arc = mp.quad(lambda th: g(taup + d * mp.e ** (1j * th)) * (1j * d * mp.e ** (1j * th)),
                      [mp.pi, mp.pi / 2, 0])
    else:
        arc = mp.quad(lambda th: g(taup + d * mp.e ** (1j * th)) * (1j * d * mp.e ** (1j * th)),
                      [mp.pi, 3 * mp.pi / 2, 2 * mp.pi])
    return mp.power(1j, -s) * (seg1 + arc + seg2)

for s in [mp.mpf('2.5'), mp.mpc('3.0', '1.3')]:
    La = Lstar_routed(s, above=True)
    Lb = Lstar_routed(s, above=False)
    res_term = 2j * mp.pi * mp.power(1j, -s) * RESA * ktilde(taup, s, k)   # 2pi i i^-s Res[f ktilde]
    print(f"  s={s}:  L_below-L_above={mp.nstr(Lb-La,10)}")
    print(f"          2pi i i^-s Res[f kt]={mp.nstr(res_term,10)}   ratio={mp.nstr((Lb-La)/res_term,8)}")

# ---------- PART B : pole AT i ----------
print("\n==== PART B : pole AT tau=i, f = Delta/E_6, weight 6 ====")
k6 = 6
print(f"  E6(i) = {mp.nstr(E6(1j), 6)}   (forced zero, i^6=-1)")
def fB(tau):
    return Delta(tau) / E6(tau)
s = mp.mpf('2.5')
prev = None
for eps in [mp.mpf('0.1'), mp.mpf('0.03'), mp.mpf('0.01'), mp.mpf('0.003')]:
    val = mp.quad(lambda t: fB(t) * ktilde(t, s, k6), [-1 + 1j, -eps + 1j])
    msg = f"  int_[i-1, {float(-eps)}+i] f kt = {mp.nstr(val, 8)}"
    if prev is not None:
        msg += f"   delta vs prev = {mp.nstr(val - prev, 6)}"
    print(msg)
    prev = val
