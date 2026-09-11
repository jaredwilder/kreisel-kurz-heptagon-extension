"""kk_extension_equations.py — Phase 2: exact extension Diophantine system for a KK heptagon.

We ask for X=(x,y) in R^2 with |X-P_i| = integer for all 7 points P_i of a certified heptagon,
retaining general position. Anchors P_0=(0,0), P_1=(a,0).  For X:
    r := |X-P_0| in Z_{>0},   s := |X-P_1| in Z_{>0}
    x = (a^2 + r^2 - s^2)/(2a),     y^2 = r^2 - x^2.

Each other point P_i = (x_i, y_i*sqrt(C)) (C = characteristic = 2002; x_i, y_i in Q). Then
    t_i^2 := |X-P_i|^2 = r^2 + rho_i - 2 x_i x - 2 y_i * (sqrt(C) y),   rho_i = |P_i|^2 = D[0][i]^2.
Introduce w := sqrt(C) * y  (so w^2 = C*(r^2 - x^2)). Then t_i^2 = A_i - 2 y_i w with A_i = r^2+rho_i-2x_i x.

Eliminating w gives the exact system:
  MASTER (M), reference index m:  (A_m - t_m^2)^2 = 4 C y_m^2 (r^2 - x^2)
  LINKAGE (L_i), i != m:          y_i t_m^2 - y_m t_i^2 = (y_i-y_m) r^2 + (y_i rho_m - y_m rho_i)
                                                          - 2 x (y_i x_m - y_m x_i)
This module prints these equations with the ACTUAL coefficients (exact), for H1 and H2.
"""
from __future__ import annotations
from fractions import Fraction as F
import sympy as sp
from kk_heptagon import H1, D1, D2, CHAR1, reconstruct_full

C = 2002


def build_system(points, D, label):
    n = len(points)
    a = D[0][1]
    rho = [D[0][i] ** 2 for i in range(n)]
    r, s = sp.symbols('r s', positive=True)
    t = sp.symbols('t2 t3 t4 t5 t6', positive=True)  # t[k] = distance to P_{k+2}
    x = sp.Rational(a**2, 1)
    x = (sp.Integer(a)**2 + r**2 - s**2) / (2 * sp.Integer(a))

    def Q(fr):  # Fraction -> sympy Rational
        return sp.Rational(fr.numerator, fr.denominator)

    # choose master reference m = 2 (0-indexed), i.e. first non-anchor with y != 0
    m = 2
    xm, ym = Q(points[m][0]), Q(points[m][1])
    Am = r**2 + rho[m] - 2 * xm * x
    tm = t[0]  # t2
    M = sp.expand((Am - tm**2)**2 - 4 * C * ym**2 * (r**2 - x**2))

    L = []
    for idx, i in enumerate(range(3, n)):   # i = 3..6 (0-indexed points P_3..P_6)
        xi, yi = Q(points[i][0]), Q(points[i][1])
        ti = t[idx + 1]
        lhs = yi * tm**2 - ym * ti**2
        rhs = (yi - ym) * r**2 + (yi * rho[m] - ym * rho[i]) - 2 * x * (yi * xm - ym * xi)
        L.append(sp.expand(lhs - rhs))

    print(f"\n================= {label} =================")
    print(f"a = {a},  characteristic C = {C},  reference point index m = {m}")
    print(f"x = (a^2 + r^2 - s^2)/(2a) = {sp.nsimplify(x)}")
    print("\n--- MASTER surface M(r,s,t2) = 0  [degree-4 hypersurface] ---")
    # clear denominators for a clean integer-coefficient equation
    Mc = sp.expand(M * sp.denom(sp.together(M)))
    print(sp.simplify(Mc), "= 0")
    print("\n--- LINKAGE equations L_i = 0 (i=3..6): each LINEAR in the squares t2^2,t_i^2,r^2,s^2 ---")
    for k, Li in enumerate(L):
        Lic = sp.expand(Li * sp.denom(sp.together(Li)))
        print(f"L_{k+3}:  {sp.simplify(Lic)} = 0")
    return {"r": r, "s": s, "t": t, "x": x, "M": M, "L": L, "a": a, "rho": rho}


def analyze_master_curve(sysd, label):
    """Slice the master surface to a plane curve and report its genus (the arithmetic heart)."""
    r, s, t = sysd["r"], sysd["s"], sysd["t"]
    M = sysd["M"]
    tm = t[0]
    # Treat the master as a curve in (s, t2) for a fixed generic r -> report degree & genus attempt.
    print(f"\n--- {label}: master as a curve in (s, t2), r symbolic ---")
    Mpoly = sp.Poly(sp.expand(M * 2), [s, tm])
    print("total degree in (s,t2):", Mpoly.total_degree())
    try:
        g = sp.Curve  # placeholder; sympy has no direct genus. Report structure instead.
    except Exception:
        pass
    # Structure: M is quadratic in t2^2 (i.e. quartic in t2) and quartic in s. Report as biquadratic.
    print("degree in t2:", sp.degree(sp.expand(M), tm), " degree in s:", sp.degree(sp.expand(M), s),
          " degree in r:", sp.degree(sp.expand(M), r))


if __name__ == "__main__":
    s1 = build_system(H1, D1, "HEPTAGON H1 (diam 22270)")
    analyze_master_curve(s1, "H1")
    # H2
    from kk_heptagon import reconstruct_full
    H2 = reconstruct_full(D2, 2002)
    s2 = build_system(H2, D2, "HEPTAGON H2 (diam 66810)")
    analyze_master_curve(s2, "H2")
