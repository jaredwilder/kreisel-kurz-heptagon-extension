"""kk_heptagon.py — Kreisel-Kurz heptagon reconstruction + certification (Phase 1 of the closure package).

Source: Kreisel & Kurz, "There are integral heptagons, no three points on a line, no four on a circle",
arXiv:0804.1303 (Discrete & Comput. Geom. 39:786-790, 2008). d(2,7)=22270.

Points of heptagon H1 live over the number field Q(sqrt(2002)) (characteristic 2002 = 2*7*11*13):
each point is (x_i, y_i * sqrt(2002)) with x_i, y_i in Q. Squared distance between (x_i,y_i*s) and
(x_j,y_j*s) with s^2 = 2002 is  (x_i-x_j)^2 + 2002*(y_i-y_j)^2  in Q, which must equal an integer square.

Exact arithmetic via fractions.Fraction — NO floats in the certificate path.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import isqrt
import json

CHAR1 = 2002  # squarefree; field Q(sqrt(2002))

# Exact coordinates from Figure 1 of the paper: (x, y) meaning point = (x, y*sqrt(2002)).
H1 = [
    (F(0), F(0)),
    (F(22270), F(0)),
    (F(26127018, 2227), F(932064, 2227)),
    (F(245363, 17), F(3144, 17)),
    (F(17615968, 2227), F(238464, 2227)),
    (F(56068, 17), F(3144, 17)),
    (F(19079044, 2227), F(-54168, 2227)),
]

# Published distance matrix (1), integer entries.
D1 = [
    [0, 22270, 22098, 16637, 9248, 8908, 8636],
    [22270, 0, 21488, 11397, 15138, 20698, 13746],
    [22098, 21488, 0, 10795, 14450, 13430, 20066],
    [16637, 11397, 10795, 0, 7395, 11135, 11049],
    [9248, 15138, 14450, 7395, 0, 5780, 5916],
    [8908, 20698, 13430, 11135, 5780, 0, 10744],
    [8636, 13746, 20066, 11049, 5916, 10744, 0],
]

# Published distance matrix (2) — second heptagon (coordinates NOT given in paper; we reconstruct).
D2 = [
    [0, 66810, 66555, 66294, 49928, 41238, 40290],
    [66810, 0, 32385, 64464, 32258, 25908, 52020],
    [66555, 32385, 0, 34191, 16637, 33147, 33405],
    [66294, 64464, 34191, 0, 34322, 53244, 26724],
    [49928, 32258, 16637, 34322, 0, 20066, 20698],
    [41238, 25908, 33147, 53244, 20066, 0, 32232],
    [40290, 52020, 33405, 26724, 20698, 32232, 0],
]


def sq_dist_field(pi, pj, char):
    """Exact squared distance between (x_i, y_i*sqrt(char)) and (x_j, y_j*sqrt(char)) as a Fraction."""
    dx = pi[0] - pj[0]
    dy = pi[1] - pj[1]
    return dx * dx + char * dy * dy


def is_perfect_square_fraction(q: F):
    """If q is a rational perfect square return its integer root (q must be a nonneg integer square)."""
    if q < 0:
        return None
    if q.denominator != 1:
        return None
    r = isqrt(q.numerator)
    return r if r * r == q.numerator else None


def det3(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
            - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
            + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))


def det4(m):
    tot = F(0)
    for j in range(4):
        minor = [[m[i][k] for k in range(4) if k != j] for i in range(1, 4)]
        tot += ((-1) ** j) * m[0][j] * det3(minor)
    return tot


def certify(points, char, published, label):
    """Phase 1 certificate: distances integer + match published; no 3 collinear; no 4 concyclic."""
    n = len(points)
    report = {"label": label, "char": char, "n": n}
    # 1) distance-matrix match (independent recomputation vs published)
    dist_ok = True
    computed = [[0]*n for _ in range(n)]
    mism = []
    for i in range(n):
        for j in range(i+1, n):
            q = sq_dist_field(points[i], points[j], char)
            root = is_perfect_square_fraction(q)
            if root is None:
                dist_ok = False
                mism.append((i, j, "not integer square", str(q)))
                continue
            computed[i][j] = computed[j][i] = root
            if published is not None and root != published[i][j]:
                dist_ok = False
                mism.append((i, j, root, published[i][j]))
    report["all_21_distances_integer_and_match"] = dist_ok
    report["mismatches"] = mism
    report["computed_distance_matrix"] = computed
    # 2) no three collinear: area determinant nonzero for every triple.
    #    point k has real coords (x_k, y_k*sqrt(char)); area ∝ det[[x_j-x_i, (y_j-y_i)], [x_k-x_i,(y_k-y_i)]]
    #    (the common sqrt(char) factors out of the y-column, nonzero iff the rational det is nonzero)
    coll = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                xi, yi = points[i]; xj, yj = points[j]; xk, yk = points[k]
                d = (xj-xi)*(yk-yi) - (yj-yi)*(xk-xi)
                if d == 0:
                    coll.append((i, j, k))
    report["no_three_collinear"] = (len(coll) == 0)
    report["collinear_triples"] = coll
    # 3) no four concyclic: in-circle determinant over the field. Rows [X^2+Y^2, X, Y, 1] with
    #    X=x, Y=y*sqrt(char). X^2+Y^2 = x^2 + char*y^2 (rational). The Y column carries sqrt(char);
    #    the 4x4 determinant = A + B*sqrt(char) for rationals A,B; concyclic iff det==0 iff A=B=0.
    conc = []
    s2 = char
    for quad in _combs(range(n), 4):
        # build rational A-part and B-part of the determinant
        # columns: c0 = x^2+char*y^2 (rat), c1 = x (rat), c2 = y (rat, coefficient of sqrt), c3 = 1
        rows = []
        for idx in quad:
            x, y = points[idx]
            rows.append([x*x + s2*y*y, x, y, F(1)])
        # determinant with c2 scaled by sqrt(char): expand -> rational part (terms with even # of c2)
        # Easiest exact route: det = det(matrix with c2 replaced by y)*sqrt(char)?? No — only column2 has sqrt.
        # Since exactly ONE column carries sqrt(char), the whole determinant is (that cofactor sum)*sqrt(char)
        # PLUS zero — i.e. determinant = sqrt(char) * det(rows_with_plain_y). So concyclic iff that det==0.
        d = det4(rows)
        if d == 0:
            conc.append(tuple(quad))
    report["no_four_concyclic"] = (len(conc) == 0)
    report["concyclic_quadruples"] = conc
    report["PHASE1_CERTIFIED"] = (dist_ok and not coll and not conc)
    return report


def _combs(it, r):
    from itertools import combinations
    return combinations(it, r)


def reconstruct_from_matrix(D):
    """Reconstruct exact coordinates over Q(sqrt(char)) from an integer distance matrix.
    Place P0=(0,0), P1=(D01,0). For each later point solve x from the two anchor distances; y^2 then
    is rational; the shared squarefree part is the characteristic. Returns (char, points) with points
    = (x_i, y_i) meaning (x_i, y_i*sqrt(char)). Sign of y chosen to match a reference point's config."""
    n = len(D)
    a = D[0][1]  # P0-P1 distance
    d0 = [F(D[0][i]) for i in range(n)]
    d1 = [F(D[1][i]) for i in range(n)]
    xs = [F(0), F(a)]
    ysq = [F(0), F(0)]  # y_i^2 as a rational (real y^2, before extracting sqrt)
    for i in range(2, n):
        # x from: x^2+Y^2 = d0i^2 ; (x-a)^2+Y^2 = d1i^2  => x = (a^2 + d0i^2 - d1i^2)/(2a)
        xi = (F(a)*F(a) + d0[i]*d0[i] - d1[i]*d1[i]) / (2*F(a))
        Ysq = d0[i]*d0[i] - xi*xi   # real y^2 (rational)
        xs.append(xi); ysq.append(Ysq)
    # characteristic = squarefree part of the (common) squarefree kernel of the nonzero Ysq
    def squarefree_part(q: F):
        # squarefree part of a positive rational = squarefree part of (num*den)
        m = q.numerator * q.denominator
        return _squarefree(m)
    chars = {squarefree_part(y) for y in ysq if y != 0}
    return chars, xs, ysq


def reconstruct_full(D, char):
    """Full reconstruction over Q(sqrt(char)): returns points (x_i, y_i) meaning (x_i, y_i*sqrt(char)),
    with signs chosen so that EVERY pairwise distance matches D. Raises if impossible over this field."""
    n = len(D)
    _, xs, ysq = reconstruct_from_matrix(D)
    # y_i (coefficient of sqrt(char)):  (y_i*sqrt(char))^2 = char*y_i^2 = ysq_i  => y_i = sqrt(ysq_i/char)
    ycoef_abs = []
    for i in range(n):
        val = ysq[i] / char if ysq[i] != 0 else F(0)
        root = is_perfect_square_fraction_rational(val)
        if root is None:
            raise ValueError(f"point {i}: ysq/char = {val} is not a rational square over Q(sqrt({char}))")
        ycoef_abs.append(root)
    # sign search: P0,P1 have y=0. Fix P2 sign = +. Determine each remaining sign by matching one distance.
    from itertools import product
    free = [i for i in range(2, n) if ycoef_abs[i] != 0]
    # P2 (first free) fixed +, brute force the rest (small)
    for signs in product([1, -1], repeat=len(free)):
        if free and signs[0] != 1:
            continue
        pts = [(xs[i], F(0)) for i in range(n)]
        for k, i in enumerate(free):
            pts[i] = (xs[i], signs[k] * ycoef_abs[i])
        ok = True
        for i in range(n):
            for j in range(i+1, n):
                q = sq_dist_field(pts[i], pts[j], char)
                r = is_perfect_square_fraction(q)
                if r is None or r != D[i][j]:
                    ok = False; break
            if not ok:
                break
        if ok:
            return pts
    raise ValueError("no sign assignment reproduces the distance matrix over this field")


def is_perfect_square_fraction_rational(q: F):
    """Rational sqrt of a nonneg rational q if it is a perfect rational square, else None."""
    if q < 0:
        return None
    rn = isqrt(q.numerator); rd = isqrt(q.denominator)
    if rn*rn == q.numerator and rd*rd == q.denominator:
        return F(rn, rd)
    return None


def _squarefree(m: int) -> int:
    if m == 0:
        return 0
    res = 1
    d = 2
    mm = abs(m)
    while d*d <= mm:
        cnt = 0
        while mm % d == 0:
            mm //= d; cnt += 1
        if cnt % 2 == 1:
            res *= d
        d += 1
    if mm > 1:
        res *= mm
    return res


if __name__ == "__main__":
    rep1 = certify(H1, CHAR1, D1, "H1 (diameter 22270)")
    print("=== HEPTAGON 1 CERTIFICATE ===")
    print("all 21 distances integer & match published:", rep1["all_21_distances_integer_and_match"])
    print("mismatches:", rep1["mismatches"])
    print("no three collinear:", rep1["no_three_collinear"])
    print("no four concyclic:", rep1["no_four_concyclic"])
    print("PHASE1_CERTIFIED:", rep1["PHASE1_CERTIFIED"])

    print("\n=== HEPTAGON 2 RECONSTRUCTION (from matrix 2) ===")
    chars2, xs2, ysq2 = reconstruct_from_matrix(D2)
    char2 = sorted(chars2)[0]
    print("characteristic of H2:", char2)
    H2 = reconstruct_full(D2, char2)
    print("reconstructed H2 coords (x, y*sqrt(%d)):" % char2)
    for p in H2:
        print("  (", p[0], ",", p[1], "* sqrt(%d) )" % char2)
    rep2 = certify(H2, char2, D2, "H2 (diameter 66810)")
    print("all 21 distances integer & match published:", rep2["all_21_distances_integer_and_match"])
    print("no three collinear:", rep2["no_three_collinear"])
    print("no four concyclic:", rep2["no_four_concyclic"])
    print("PHASE1_CERTIFIED:", rep2["PHASE1_CERTIFIED"])

    with open("oracle/kbk/engine/kk_heptagon_certificates.json", "w") as f:
        json.dump({
            "H1": {"coords": [[str(p[0]), str(p[1])] for p in H1], "char": CHAR1,
                   "certified": rep1["PHASE1_CERTIFIED"],
                   "distances_match": rep1["all_21_distances_integer_and_match"],
                   "computed_matrix": rep1["computed_distance_matrix"]},
            "H2": {"coords": [[str(p[0]), str(p[1])] for p in H2], "char": char2,
                   "certified": rep2["PHASE1_CERTIFIED"],
                   "distances_match": rep2["all_21_distances_integer_and_match"],
                   "computed_matrix": rep2["computed_distance_matrix"]},
        }, f, indent=2)
    print("\nwrote kk_heptagon_certificates.json")
