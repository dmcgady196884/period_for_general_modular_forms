#!/usr/bin/env python3
"""Checks behind Figure 1 of finite_contour_cocycles_short.tex: the passage from
the region T of Definition 4.21 (panel a) to the locus T_F inside the standard
fundamental domain F (panel b).  Poles in T_F are not excluded; they are the ones
carrying the correction Q_f.

    T = { 0 <= Re t <= 1,  1/2 <= Im t <= 1,  |t - i/2| >= 1/2,  |t - 1 - i/2| >= 1/2 }
    F = { |Re t| <= 1/2,  |t| >= 1 }

The claim proved in the note and checked here is

    (eq:foldT)      SL2(Z) . T   cap  F   =   { t in F : Im t <= 1 },

together with the two areas.  The proof of "no point of T can be lifted above
Im t = 1" runs by cases on the lower row (c,d) of gamma, since
Im(gamma t) = Im t / |c t + d|^2 exceeds 1 only if |c t + d|^2 < Im t <= 1:

    c = 0      needs |d| < 1, impossible;
    |c| >= 2   gives |c t + d|^2 >= c^2 (Im t)^2 >= 1  because Im t >= 1/2;
    |c| = 1    asks |t - m - i/2| < 1/2 for some integer m, which on the strip
               0 <= Re t <= 1 can only happen for m = 0 or 1 -- and those two
               disks are exactly what the last two inequalities of T exclude.

So T's defining inequalities are precisely the conditions making it unliftable,
and the fold is 3:1 because U = TS fixes rho+1 in T and rotates T onto itself.

Run:  python3 verify_regions.py
"""
import math
import random

# --------------------------------------------------------------- the regions --
def in_T(z, eps=0.0):
    x, y = z.real, z.imag
    return (-eps <= x <= 1 + eps and 0.5 - eps <= y <= 1 + eps
            and abs(z - 0.5j) >= 0.5 - eps and abs(z - 1 - 0.5j) >= 0.5 - eps)

def in_F(z, eps=1e-12):
    return abs(z.real) <= 0.5 + eps and abs(z) >= 1 - eps

# ------------------------------------------------------- a ball in SL2(Z) ----
def group_ball(depth=8, cap=600):
    S, T, Ti = (0, -1, 1, 0), (1, 1, 0, 1), (1, -1, 0, 1)
    def mul(a, b):
        return (a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
                a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3])
    seen = {(1, 0, 0, 1)}; frontier = [(1, 0, 0, 1)]
    for _ in range(depth):
        nxt = []
        for g in frontier:
            for h in (S, T, Ti):
                p = mul(h, g)
                if p not in seen and max(map(abs, p)) < cap:
                    seen.add(p); nxt.append(p)
        frontier = nxt
    return sorted(seen)

BALL = group_ball()

def act(g, z):
    a, b, c, d = g
    return (a * z + b) / (c * z + d)

def moebius_orbit_in_T(z, eps=1e-11):
    """Distinct images of z lying in T (deduplicated: g and -g act alike)."""
    out = []
    for g in BALL:
        w = act(g, z)
        if in_T(w, eps) and not any(abs(w - v) < 1e-9 for v in out):
            out.append(w)
    return out

# ------------------------------------------------------------------ checks ---
def check_area_T(n=2_000_000):
    """area(T) = int_0^1 [1/y_lo(x) - 1] dx,  y_lo = 1/2 + sqrt(1/4 - min(x,1-x)^2)."""
    tot = 0.0
    for k in range(n):
        x = (k + 0.5) / n
        u = min(x, 1.0 - x)
        tot += 1.0 / (0.5 + math.sqrt(max(0.0, 0.25 - u * u))) - 1.0
    return tot / n

def check_fold_multiplicity(trials=300, seed=11):
    rng = random.Random(seed)
    hist = {}
    done = 0
    while done < trials:
        z = complex(rng.uniform(0, 1), rng.uniform(0.5, 1))
        if not in_T(z, -1e-3):                 # strictly interior
            continue
        done += 1
        m = len(moebius_orbit_in_T(z))
        hist[m] = hist.get(m, 0) + 1
    return hist

def check_foldT(trials=4000, seed=3):
    """SL2(Z).T cap F  ==  { t in F : Im t <= 1 }."""
    rng = random.Random(seed)
    bad = []
    done = 0
    while done < trials:
        z = complex(rng.uniform(-0.5, 0.5), rng.uniform(0.80, 6.0))
        if not in_F(z) or abs(z.imag - 1.0) < 2e-3:      # skip the boundary
            continue
        done += 1
        if bool(moebius_orbit_in_T(z)) != (z.imag <= 1.0):
            bad.append(z)
    return done, bad

if __name__ == "__main__":
    aT = check_area_T()
    print("area(T)                    = %.10f   (pi - 3     = %.10f)"
          % (aT, math.pi - 3))
    print("area{t in F : Im t <= 1}   = %.10f   (pi/3 - 1   = %.10f)"
          % (2 * math.asin(0.5) - 1, math.pi / 3 - 1))
    print("   ... and (pi-3)/3        = %.10f, so the fold is 3:1"
          % ((math.pi - 3) / 3))
    print("   T_F as a fraction of F  = %.6f   (1 - 3/pi)" % (1 - 3 / math.pi))
    print()
    print("|Moebius orbit of an interior point of T, intersected with T|:")
    print("   histogram over 300 points:", check_fold_multiplicity())
    print("   (3 = the order of U = TS, which fixes rho+1 and rotates T onto itself)")
    print()
    n, bad = check_foldT()
    print("eq:foldT tested on %d points of F: %d mismatches" % (n, len(bad)))
    for z in bad[:5]:
        print("   MISMATCH at %.6f%+.6fi" % (z.real, z.imag))
