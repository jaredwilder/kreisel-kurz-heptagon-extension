"""kk_local_rigidity.py — local deformation theory of the 7 KNOWN solutions on the extension surface.

Answers the operator's question directly: are the seven trivial extensions X=P_j rigid at first order,
or does the obstruction live higher up / globally?

SETUP (exact, over K=Q(sqrt2002)). Free rational coordinates (x, w) where the candidate 8th point is
X = (x, w/sqrt(2002)) in the real plane. For heptagon vertex P_i = (x_i, y_i*sqrt(2002)):

    Q_i(x,w) := 2002 (x - x_i)^2 + (w - 2002 y_i)^2       [= 2002 * |X-P_i|^2]
    extension condition:  F_i := 2002 t_i^2 - Q_i(x,w) = 0,  t_i = |X - P_i|

So V = {F_0 = ... = F_6 = 0} subset A^9 with coordinates (x, w, t_0..t_6).
A trivial solution is X = P_j: there t_j = 0 and t_i = d_ij (the published integer distances).

WHAT THIS COMPUTES
  1. Jacobian rank of V at each trivial point  (rank 7 = smooth; rank < 7 = singular).
  2. The exact tangent cone at a trivial point.
  3. Whether that tangent cone has RATIONAL points (= rational tangent directions exist
     => the 7 known solutions are NOT first-order rigid, and the obstruction is higher-order/global).

Run: python oracle/kbk/engine/kk_local_rigidity.py
"""
from __future__ import annotations
from fractions import Fraction as F
import sys
sys.path.insert(0, "oracle/kbk/engine")
import sympy as sp
from kk_heptagon import H1, D1, D2, reconstruct_full

C = 2002


def _Q(fr):
    return sp.Rational(fr.numerator, fr.denominator)


def analyse(points, D, label):
    n = len(points)
    x, w = sp.symbols('x w')
    ts = sp.symbols(f't0:{n}')
    xs = [_Q(p[0]) for p in points]
    ys = [_Q(p[1]) for p in points]

    Fs = []
    for i in range(n):
        Qi = C * (x - xs[i]) ** 2 + (w - C * ys[i]) ** 2
        Fs.append(sp.expand(C * ts[i] ** 2 - Qi))

    variables = [x, w] + list(ts)
    Jac = sp.Matrix([[sp.diff(f, v) for v in variables] for f in Fs])

    print(f"=== {label} ===")
    print(f"V = {{F_0..F_{n-1}}} in A^{len(variables)}  ->  expected dim {len(variables) - n} (a surface)")
    print()

    # --- 1. Jacobian rank at each trivial point X = P_j -------------------------------
    print("Jacobian rank at each KNOWN (trivial) solution X = P_j:")
    ranks = {}
    for j in range(n):
        sub = {x: xs[j], w: C * ys[j]}
        for i in range(n):
            sub[ts[i]] = sp.Integer(0) if i == j else sp.Integer(D[i][j])
        Jj = Jac.subs(sub)
        rk = Jj.rank()
        zero_rows = [i for i in range(n) if all(Jj[i, k] == 0 for k in range(len(variables)))]
        ranks[j] = rk
        print(f"  P{j}: rank = {rk}  (smooth needs {n})  "
              f"tangent-space dim = {len(variables) - rk}   identically-zero rows: {zero_rows}")
    print()

    # --- 2. where IS V singular? the t-block is diagonal(2*C*t_i) ---------------------
    print("Singular locus, exactly: the t-block of the Jacobian is diag(2*2002*t_i),")
    print("so rank drops iff some t_i = 0 iff Q_i(x,w) = 0. Over R, Q_i is a sum of two")
    print("squares, so Q_i = 0 <=> X = P_i.  =>  the REAL singular points of V are EXACTLY")
    print("the n known trivial solutions. The 'trivial' points are the surface's nodes.")
    print()

    # --- 3. tangent cone at a trivial point + rational points on it ------------------
    j = 2  # any non-anchor vertex
    u, v, tj = sp.symbols('u v tj')
    # substitute x = x_j + u, w = 2002 y_j + v ; F_j has NO linear part:
    Fj_local = sp.expand(Fs[j].subs({x: xs[j] + u, w: C * ys[j] + v, ts[j]: tj}))
    print(f"Tangent cone at P{j}: F_{j} in local coords (u,v,tj) = x-x_{j}, w-2002y_{j}, t_{j}:")
    print(f"  {sp.factor(Fj_local)} = 0")
    expected = sp.expand(C * tj ** 2 - C * u ** 2 - v ** 2)
    print(f"  equals 2002*tj^2 - 2002*u^2 - v^2 ?  {sp.simplify(Fj_local - expected) == 0}")
    print("  (the other F_i, i != j, have NONZERO linear parts and solve for t_i by the")
    print("   implicit function theorem, leaving exactly this quadric cone in (u,v,tj).)")
    print()

    # rational points on the cone 2002 tj^2 = 2002 u^2 + v^2
    # v^2 = 2002(tj-u)(tj+u); put m = tj-u, nn = tj+u  ->  v^2 = 2002 m nn
    # m = 2002, nn = 1  ->  v = 2002, tj = 2003/2, u = -2001/2
    tj0, u0, v0 = sp.Rational(2003, 2), sp.Rational(-2001, 2), sp.Integer(2002)
    check = sp.simplify(C * tj0 ** 2 - C * u0 ** 2 - v0 ** 2)
    print(f"RATIONAL point on the tangent cone: (u,v,tj) = ({u0}, {v0}, {tj0})   residual = {check}")
    print(f"  => the node is SPLIT over Q: rational tangent directions to V exist at P{j}.")
    print(f"  => the 7 known solutions are NOT first-order rigid. No obstruction at linear order.")
    print()
    return ranks


def scaling_check():
    """The reduction that matters: rational-distance 8-sets and integral 8-sets are the SAME problem."""
    print("=== SCALE-INVARIANCE OF THE CHARACTERISTIC (why the search space was too narrow) ===")
    # char = squarefree part of 16*Area^2 ; scaling all coords by N scales Area by N^2,
    # so 16Area^2 by N^4, so the squarefree part is UNCHANGED.
    from kk_heptagon import _squarefree
    for N in (1, 2, 3, 5, 7):
        pts = [(p[0] * N, p[1] * N) for p in H1]
        # area of triangle (P0,P1,P2) with real y = y*sqrt(2002): 2*Area = |det| * sqrt(2002)
        (x0, y0), (x1, y1), (x2, y2) = pts[0], pts[1], pts[2]
        det = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)
        area2 = (det * det * C) / 4          # Area^2  (the sqrt(2002) squares out)
        val = 16 * area2
        sf = _squarefree(val.numerator * val.denominator)
        print(f"  N={N}: char(N*H1) = squarefree(16*Area^2) = {sf}")
    print()
    print("  => scaling H1 by N keeps characteristic 2002, so scaled copies stay legal targets.")
    print("  => An 8-point set with all 28 distances RATIONAL scales (by the common denominator)")
    print("     to an integral 8-point set, and general position is scale-invariant. Hence:")
    print("       exists integral general-position 8-set  <=>  exists RATIONAL-distance one.")
    print("  => searching integer (r,s) against H1 only finds denominator-1 extensions.")
    print("     Rational extensions of H1 = integer extensions of N*H1 for some N > 1.")
    print("     MY B=100000 SEARCHES WERE THE N=1 SLICE ONLY.")


if __name__ == "__main__":
    analyse(H1, D1, "H1 (diameter 22270)")
    scaling_check()
