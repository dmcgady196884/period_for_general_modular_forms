"""Does eq:geoG hold for a pole in the LENS between the arc and the chord?

lem:geopolylog assumes Im z > sqrt3/2 AND |z| >= 1.  I claimed the second could always be
arranged by a T-translation, since the block Li_{-N}(e^{2 pi i(tau - z)}) is literally the
same function for z and z+1.  DAM called that cleverness-as-self-trickery.

The point at issue: the block's poles are ALL of z + Z.  If z lies in the lens
    Lambda = { |tau| < 1, Im tau > sqrt3/2, |Re tau| <= 1/2 }
then a pole lies in Lambda no matter which translate we name, so deforming arc -> chord
sweeps across it.  If so eq:geoG returns the CHORD value, and the true ARC value differs
by 2 pi i times the residue there.

Computed here, at k = 12, s = 3, N = 0 (so Li_0(w) = w/(1-w)):
  (1) L* of the block by direct quadrature on the ARC        <- ground truth
  (2) L* of the block by direct quadrature on the CHORD
  (3) eq:geoIN / eq:geoG evaluated at z
  (4) the same evaluated at z+1
  (5) 2 pi i e^{-i pi s/2} (Res_S + Res_T) at z              <- predicted arc-chord gap
CONTROL: z = 0.3 + 1.4i, well outside the lens, where all of these must agree.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import mp, I, pi, ktil

mp.mp.dps = 30
KK = 12
S = mp.mpf(3)
RHO = mp.e**(2*I*pi/3); TAU0 = mp.e**(I*pi/3)

def block(tau, z):                      # Li_0(w) = w/(1-w)
    w = mp.e**(2*I*pi*(tau - z))
    return w/(1-w)

def arc_int(g, panels=60):
    a,b = pi/3, 2*pi/3
    th=[a+(b-a)*m/panels for m in range(panels+1)]
    return -sum(mp.quad(lambda x: g(mp.e**(I*x))*I*mp.e**(I*x), [u,v])
                for u,v in zip(th[:-1],th[1:]))

def chord_int(g, panels=60):
    A,B = RHO, TAU0                     # both at height sqrt3/2, horizontal
    pts=[A+(B-A)*m/panels for m in range(panels+1)]
    return sum(mp.quad(lambda t: g(t)*(B-A), [mp.mpf(m)/panels,(mp.mpf(m)+1)/panels])
               for m in range(panels)) if False else \
           sum(mp.quad(lambda x: g(A+(B-A)*x)*(B-A), [mp.mpf(m)/panels,(mp.mpf(m)+1)/panels])
               for m in range(panels))

def Lstar(g, integrator):
    return mp.e**(-I*pi*S/2)*(integrator(lambda t: g(t)*t**(S-1))
                              + integrator(lambda t: g(t)*ktil(t,S,KK)))

def gstar(s,x,J=200):
    return sum((-x)**j/(mp.factorial(j)*(j+s)) for j in range(J))

def Ggeo(n):                            # eq:geoG for mode -n, n >= 1
    a = 2*pi*I*n
    Spart = mp.e**(-I*pi*S/2)*(TAU0**S*gstar(S,a*TAU0) - RHO**S*gstar(S,a*RHO))
    base = mp.e**(mp.log(2*pi*n)+I*pi)   # (-2 pi n) with principal log
    Tpart = mp.gammainc(S, a*TAU0, mp.inf)/base**S \
          + I**KK*mp.gammainc(KK-S, a*TAU0, mp.inf)/base**(KK-S)
    return Spart + Tpart

def Ggeo0():                            # eq:geoG0
    return mp.e**(-I*pi*S/2)*(TAU0**S - RHO**S)/S \
           - (TAU0/I)**S/S - I**KK*(TAU0/I)**(KK-S)/(KK-S)

def formula(z, NMAX=60):
    return -Ggeo0() - sum(mp.e**(2*I*pi*n*z)*Ggeo(n) for n in range(1,NMAX+1))

def residue_gap(z):
    """2 pi i e^{-i pi s/2} [Res_z(block * tau^{s-1}) + Res_z(block * ktil)]"""
    # Li_0(w) = w/(1-w) has residue -1/(2 pi i) in tau at tau = z
    R = -1/(2*pi*I)
    return 2*pi*I*mp.e**(-I*pi*S/2)*R*(z**(S-1) + ktil(z,S,KK))

for lbl, z in (("CONTROL  z = 0.3+1.4i (outside lens, high)", mp.mpc('0.3','1.4')),
               ("LENS     z = 0.2+0.93i (|z|<1, Im>sqrt3/2)", mp.mpc('0.2','0.93'))):
    print("="*74); print(lbl)
    print("   |z| = %s   Im z = %s   sqrt3/2 = %s"
          % (mp.nstr(abs(z),8), mp.nstr(z.imag,8), mp.nstr(mp.sqrt(3)/2,8)))
    arc  = Lstar(lambda t: block(t,z), arc_int)
    cho  = Lstar(lambda t: block(t,z), chord_int)
    f_z  = formula(z)
    f_z1 = formula(z+1)
    print("   (1) arc quadrature   = %s" % mp.nstr(arc, 14))
    print("   (2) chord quadrature = %s" % mp.nstr(cho, 14))
    print("   (3) eq:geoG at z     = %s" % mp.nstr(f_z, 14))
    print("   (4) eq:geoG at z+1   = %s" % mp.nstr(f_z1, 14))
    sc = abs(arc)
    print("   |(3)-(1)|/|arc| = %-12s   |(3)-(2)|/|arc| = %s"
          % (mp.nstr(abs(f_z-arc)/sc,6), mp.nstr(abs(f_z-cho)/sc,6)))
    print("   |(4)-(1)|/|arc| = %-12s   |(4)-(3)|/|arc| = %s"
          % (mp.nstr(abs(f_z1-arc)/sc,6), mp.nstr(abs(f_z1-f_z)/sc,6)))
    print("   arc - chord      = %s" % mp.nstr(arc-cho, 12))
    print("   predicted gap    = %s" % mp.nstr(residue_gap(z), 12))
    print("   |gap - predicted|/|arc| = %s" % mp.nstr(abs((arc-cho)-residue_gap(z))/sc, 6))
    print()
