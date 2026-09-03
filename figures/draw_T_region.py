#!/usr/bin/env python3
"""Draw the region T of Definition 4.21 of finite_contour_cocycles_short.tex,
superimposed on the standard fundamental domain F, its S-image SF, and a few
T-translates of both.  Emits a standalone SVG; no third-party dependencies.

    T = { 0 <= Re t <= 1,  1/2 <= Im t <= 1,  |t - i/2| >= 1/2,  |t - 1 - i/2| >= 1/2 }

Its three sides:
    top    Im t = 1                 (the T-segment of the reference contour)
    right  |t - 1 - i/2| = 1/2      (circle through 1 and 1+i)
    left   |t - i/2|     = 1/2      (circle through 0 and i)
the last two meeting tangentially at (1+i)/2.

Also prints the hyperbolic area of T as a check on the geometry: it should be
exactly pi - 3.
"""
import math

# ----------------------------------------------------------------- geometry --
R3 = math.sqrt(3.0) / 2.0            # Im(rho)
RHO = (-0.5, R3)                     # rho     = e^{2 pi i/3}
RHO1 = (0.5, R3)                     # rho + 1

def arcpts(cx, cy, r, a0, a1, n=320):
    """Points on the circle of centre (cx,cy), radius r, from angle a0 to a1."""
    return [(cx + r * math.cos(a0 + (a1 - a0) * k / n),
             cy + r * math.sin(a0 + (a1 - a0) * k / n)) for k in range(n + 1)]

def shift(pts, dx):
    return [(x + dx, y) for (x, y) in pts]

D = math.radians(1.0) * 0            # (no fudge; kept explicit)

def fund_domain(ytop):
    """F = { |Re t| <= 1/2, |t| >= 1 }, clipped at Im t = ytop."""
    return ([(-0.5, ytop), RHO]
            + arcpts(0.0, 0.0, 1.0, math.radians(120), math.radians(60))
            + [(0.5, ytop)])

def S_fund_domain():
    """S.F, the curvilinear triangle with vertices rho, rho+1 and the cusp 0.
    Bounded above by |t| = 1 and below by |t-1| = 1 and |t+1| = 1."""
    return (arcpts(0.0, 0.0, 1.0, math.radians(120), math.radians(60))
            + arcpts(1.0, 0.0, 1.0, math.radians(120), math.radians(180))
            + arcpts(-1.0, 0.0, 1.0, math.radians(0), math.radians(60)))

def T_region():
    """The region T itself: top edge, then right circle, then left circle."""
    return ([(0.0, 1.0), (1.0, 1.0)]
            + arcpts(1.0, 0.5, 0.5, math.radians(90), math.radians(180))
            + arcpts(0.0, 0.5, 0.5, math.radians(0), math.radians(90)))

# --------------------------------------------------------- hyperbolic area --
def hyperbolic_area_of_T(n=2_000_000):
    """area = int_0^1 [ 1/y_lo(x) - 1 ] dx  with  y_lo = 1/2 + sqrt(1/4 - min(x,1-x)^2).
    Midpoint rule; the integrand is smooth on each half."""
    tot = 0.0
    for k in range(n):
        x = (k + 0.5) / n
        u = min(x, 1.0 - x)
        ylo = 0.5 + math.sqrt(max(0.0, 0.25 - u * u))
        tot += 1.0 / ylo - 1.0
    return tot / n

# ------------------------------------------------------------------- canvas --
XMIN, XMAX, YMAX = -0.60, 1.60, 1.30
SCALE, MARGIN = 520.0, 40.0
STRIP = 56.0                         # band above the plot, holding the key
W = (XMAX - XMIN) * SCALE + 2 * MARGIN
H = YMAX * SCALE + 2 * MARGIN + STRIP

def P(p):
    x, y = p
    return ((x - XMIN) * SCALE + MARGIN, (YMAX - y) * SCALE + MARGIN + STRIP)

def path(pts, close=True):
    d = "M " + " L ".join("%.2f,%.2f" % P(p) for p in pts)
    return d + (" Z" if close else "")

def text(p, s, size=19, anchor="middle", cls="lbl", dx=0.0, dy=0.0):
    X, Y = P(p)
    return ('<text x="%.2f" y="%.2f" text-anchor="%s" class="%s" '
            'font-size="%d">%s</text>' % (X + dx, Y + dy, anchor, cls, size, s))

def build():
    o = []
    o.append('<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" '
             'viewBox="0 0 %.0f %.0f">' % (W, H, W, H))
    o.append("""<style>
    text.lbl  { font-family: Georgia,'Times New Roman',serif; fill:#1c1c1c; }
    text.key  { font-family: Georgia,'Times New Roman',serif; fill:#1c1c1c; }
    text.cal  { font-family: Georgia,'Times New Roman',serif; font-style:italic; }
    </style>""")
    o.append('<rect width="100%" height="100%" fill="#ffffff"/>')

    # --- translates of F and of S.F, drawn first (background) ---
    for dx in (-1.0, 0.0, 1.0):
        o.append('<path d="%s" fill="#dce7f2" fill-opacity="0.85" '
                 'stroke="#7f9ab4" stroke-width="1.4"/>'
                 % path(shift(fund_domain(YMAX + 0.2), dx)))
    for dx in (-1.0, 0.0, 1.0, 2.0):
        o.append('<path d="%s" fill="#dfeedd" fill-opacity="0.9" '
                 'stroke="#7ba475" stroke-width="1.4"/>'
                 % path(shift(S_fund_domain(), dx)))

    # --- the line Im t = 1: the T-segment the whole construction hangs off ---
    o.append('<path d="%s" stroke="#444" stroke-width="1.6" '
             'stroke-dasharray="9,6" fill="none"/>'
             % path([(XMIN, 1.0), (XMAX, 1.0)], close=False))
    o.append(text((XMAX - 0.02, 1.03), "Im τ = 1", 17, "end"))

    # --- T itself ---
    o.append('<path d="%s" fill="#f0a860" fill-opacity="0.62" '
             'stroke="#a8410a" stroke-width="4.0" stroke-linejoin="round"/>'
             % path(T_region()))

    # --- real axis + ticks ---
    o.append('<path d="%s" stroke="#222" stroke-width="1.6" fill="none"/>'
             % path([(XMIN, 0.0), (XMAX, 0.0)], close=False))
    for t in (-0.5, 0.0, 0.5, 1.0, 1.5):
        o.append('<path d="%s" stroke="#222" stroke-width="1.4" fill="none"/>'
                 % path([(t, 0.0), (t, -0.018)], close=False))
        lab = {-0.5: "−½", 0.0: "0", 0.5: "½", 1.0: "1", 1.5: "3/2"}[t]
        o.append(text((t, 0.0), lab, 17, "middle", "lbl", 0, 28))

    # --- marked points ---
    for p, lab, dx, dy in [((0.0, 1.0), "i", -14, -10),
                           ((1.0, 1.0), "1+i", 26, -10),
                           ((0.5, 0.5), "(1+i)/2", 0, 26),
                           (RHO, "ρ", -16, 4),
                           (RHO1, "ρ+1", 26, 4)]:
        X, Y = P(p)
        o.append('<circle cx="%.2f" cy="%.2f" r="4.6" fill="#a8410a"/>' % (X, Y))
        o.append(text(p, lab, 18, "middle", "lbl", dx, dy))

    # --- region labels ---
    o.append(text((0.0, 1.16), "ℱ", 30, "middle", "cal"))
    o.append(text((1.0, 1.16), "Tℱ", 26, "middle", "cal"))
    o.append(text((0.0, 0.28), "Sℱ", 24, "middle", "cal"))
    o.append(text((1.0, 0.28), "TSℱ", 24, "middle", "cal"))
    o.append(text((0.5, 0.70), "\U0001d4af", 36, "middle", "cal"))

    # --- key, in its own strip above the plot ---
    cols = [("#f0a860", "#a8410a", "\U0001d4af — Definition 4.21, boundary in bold"),
            ("#dce7f2", "#7f9ab4", "ℱ = { |Re τ| ≤ ½, |τ| ≥ 1 }"),
            ("#dfeedd", "#7ba475", "Sℱ,  S: τ ↦ −1/τ")]
    x = MARGIN + 4
    for fill, stroke, lab in cols:
        o.append('<rect x="%.0f" y="%.0f" width="30" height="17" fill="%s" '
                 'stroke="%s" stroke-width="2"/>' % (x, MARGIN - 4, fill, stroke))
        o.append('<text x="%.0f" y="%.0f" class="key" font-size="18">%s</text>'
                 % (x + 40, MARGIN + 10, lab))
        x += 40 + 8.7 * len(lab) + 46

    o.append("</svg>")
    return "\n".join(o)

if __name__ == "__main__":
    area = hyperbolic_area_of_T()
    print("hyperbolic area of T = %.10f" % area)
    print("pi - 3               = %.10f" % (math.pi - 3.0))
    print("difference           = %.2e" % abs(area - (math.pi - 3.0)))
    out = __file__.rsplit("/", 1)[0] + "/T_region.svg"
    with open(out, "w") as fh:
        fh.write(build())
    print("wrote", out)
