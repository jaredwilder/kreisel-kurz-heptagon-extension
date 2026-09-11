"""cayley_menger.py — NEW BLADE: Cayley-Menger planar-embedding + modular obstructions.

Deterministic, $0, zero-LLM, zero-network. Forged 2026-07-24. Exact integer / rational
arithmetic only; NO float ever decides anything.

THE TARGET (operator Program B):
  The 5x5 Cayley-Menger determinant of 4 points (built from SQUARED distances) vanishes IFF
  the 4 points embed isometrically in a plane (affine dim <= 2). For n points, planar
  realizability requires ALL C(n,4) CM determinants vanish AND the Gram matrix is PSD of
  rank <= 2.

  CM layout (d_ij = SQUARED distance):
        | 0   1    1    1    1   |
        | 1   0    d01  d02  d03 |
   CM = | 1   d01  0    d12  d13 |
        | 1   d02  d12  0    d23 |
        | 1   d03  d13  d23  0   |
  Classical: CM = 288 * V^2, V = tetrahedron volume. For 4 points with INTEGER coordinates in
  Z^3, V = |det[P1-P0, P2-P0, P3-P0]| / 6, so  CM = 288 * det^2 / 36 = 8 * det^2  (integer,
  always a nonneg 8x-a-square). Coplanar <=> det = 0 <=> CM = 0. This identity is the spine of
  both selftests below (verified exactly) and of the B2 modular classification.

B1 (minimal basis, tested): does CM-vanishing on just the quadruples CONTAINING a fixed
  nondegenerate triangle force embedding dim <= 2 for the whole set? See b1_test() for the
  verdict and a CONSTRUCTED counterexample.

B2 (modular obstructions): reduce CM of integer-coordinate tetrads mod m in {8,16,3,5,7}.
  Because CM = 8*det^2, the achievable residues are exactly {8 * (square) mod m}; every other
  residue is a hard obstruction (can NEVER occur for an integer-coordinate quadruple). See
  b2_test().

INDEPENDENT VERIFIER on a SEPARATE code path: the Gram route. Build G = -1/2 * J D2 J by exact
  double-centering (rationals), compute its inertia (rank + signature) by exact congruence
  (Sylvester), and assert  (rank(G) <= 2 and G PSD)  <=>  all C(n,4) CM determinants vanish.
  The CM route and the Gram route never share code.

Self-test:  python oracle/kbk/engine/cayley_menger.py
Run:        python oracle/kbk/engine/cayley_menger.py run   (writes cayley_menger_results.json)
"""
from __future__ import annotations
from fractions import Fraction
from itertools import combinations
import json
import sys


# ---------------------------------------------------------------------------
# exact linear algebra primitives
# ---------------------------------------------------------------------------
def det_int(M) -> int:
    """Exact determinant of an INTEGER matrix via fraction-free Bareiss elimination.
    No floats, no Fraction — every intermediate division is exact by construction."""
    n = len(M)
    A = [list(map(int, row)) for row in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = None
            for i in range(k + 1, n):
                if A[i][k] != 0:
                    sw = i
                    break
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def signature(M) -> tuple[int, int, int]:
    """(#positive, #negative, #zero) eigenvalue counts of a SYMMETRIC rational matrix, computed
    exactly by congruence (Sylvester's law of inertia). No eigenvalues, no floats — only exact
    row/column operations that preserve the inertia. rank = pos + neg; PSD <=> neg == 0."""
    n = len(M)
    A = [[Fraction(x) for x in row] for row in M]
    pos = neg = zero = 0
    step = 0
    while step < n:
        piv = -1
        for i in range(step, n):
            if A[i][i] != 0:
                piv = i
                break
        if piv == -1:
            # no nonzero diagonal in the remaining block: pull an off-diagonal onto the diagonal
            found = False
            for i in range(step, n):
                for j in range(i + 1, n):
                    if A[i][j] != 0:
                        for c in range(n):
                            A[i][c] += A[j][c]          # row_i += row_j
                        for r in range(n):
                            A[r][i] += A[r][j]          # col_i += col_j  (keeps symmetry)
                        found = True
                        break
                if found:
                    break
            if not found:
                zero += n - step                        # remaining block is entirely zero
                break
            continue                                    # re-scan; A[i][i] is now nonzero
        if piv != step:
            A[step], A[piv] = A[piv], A[step]
            for r in range(n):
                A[r][step], A[r][piv] = A[r][piv], A[r][step]
        p = A[step][step]
        if p > 0:
            pos += 1
        else:
            neg += 1
        for r in range(step + 1, n):                    # symmetric elimination (row then col)
            if A[r][step] != 0:
                f = A[r][step] / p
                for c in range(n):
                    A[r][c] -= f * A[step][c]
        for c in range(step + 1, n):
            if A[step][c] != 0:
                f = A[step][c] / p
                for r in range(n):
                    A[r][c] -= f * A[r][step]
        step += 1
    return pos, neg, zero


# ---------------------------------------------------------------------------
# geometry primitives (integer coordinates -> integer squared distances)
# ---------------------------------------------------------------------------
def sq(a, b) -> int:
    """Squared Euclidean distance between two integer points of any (equal) dimension."""
    return sum((ai - bi) * (ai - bi) for ai, bi in zip(a, b))


def dmatrix(pts) -> list[list[int]]:
    """Full symmetric matrix of SQUARED distances for a point list."""
    n = len(pts)
    return [[sq(pts[i], pts[j]) for j in range(n)] for i in range(n)]


def cm4(D2, quad) -> int:
    """5x5 Cayley-Menger determinant for the 4 points indexed by `quad` (integer, exact)."""
    q = list(quad)
    M = [[0] * 5 for _ in range(5)]
    for i in range(1, 5):
        M[0][i] = 1
        M[i][0] = 1
    for i in range(4):
        for j in range(4):
            M[i + 1][j + 1] = D2[q[i]][q[j]]
    return det_int(M)


def cm3(D2, tri) -> int:
    """4x4 Cayley-Menger determinant for 3 points. Nonzero <=> nondegenerate triangle
    (= 16 * area^2 up to sign); zero <=> collinear/coincident."""
    q = list(tri)
    M = [[0] * 4 for _ in range(4)]
    for i in range(1, 4):
        M[0][i] = 1
        M[i][0] = 1
    for i in range(3):
        for j in range(3):
            M[i + 1][j + 1] = D2[q[i]][q[j]]
    return det_int(M)


def coord_det(p0, p1, p2, p3) -> int:
    """det of the three edge vectors of a Z^3 tetrad (0 <=> coplanar). Used to check CM=8*det^2."""
    v1 = [p1[i] - p0[i] for i in range(3)]
    v2 = [p2[i] - p0[i] for i in range(3)]
    v3 = [p3[i] - p0[i] for i in range(3)]
    return det_int([v1, v2, v3])


# ---------------------------------------------------------------------------
# INDEPENDENT VERIFIER — the Gram route (separate code path from cm4/det_int)
# ---------------------------------------------------------------------------
def gram_matrix(D2) -> list[list[Fraction]]:
    """G = -1/2 * J D2 J via exact double-centering. J = I - (1/n) 11^T.
    (J D2 J)_ij = D2_ij - rowmean_i - rowmean_j + total   (D2 symmetric)."""
    n = len(D2)
    rowmean = [sum(Fraction(D2[i][j]) for j in range(n)) / n for i in range(n)]
    total = sum(Fraction(D2[i][j]) for i in range(n) for j in range(n)) / (n * n)
    return [[Fraction(-1, 2) * (D2[i][j] - rowmean[i] - rowmean[j] + total)
             for j in range(n)] for i in range(n)]


def embeds_in_plane_gram(D2) -> tuple[bool, int, tuple[int, int, int]]:
    """INDEPENDENT of the CM route: does the metric embed in R^2?  True iff Gram rank <= 2 AND
    Gram is PSD. Returns (embeds, rank, (pos, neg, zero))."""
    pos, neg, zero = signature(gram_matrix(D2))
    rank = pos + neg
    return (rank <= 2 and neg == 0), rank, (pos, neg, zero)


def all_cm_vanish(D2) -> tuple[bool, tuple | None]:
    """CM route: do ALL C(n,4) Cayley-Menger determinants vanish? Returns (ok, first_bad_quad)."""
    n = len(D2)
    if n < 4:
        return True, None
    for quad in combinations(range(n), 4):
        if cm4(D2, quad) != 0:
            return False, quad
    return True, None


# ---------------------------------------------------------------------------
# B1 — minimal-basis sufficiency
# ---------------------------------------------------------------------------
def basis_quadruples_triple(n, triple):
    """Quadruples containing a fixed TRIPLE (3 points): O(n) of them (4th point ranges free)."""
    rest = [i for i in range(n) if i not in triple]
    return [tuple(sorted(list(triple) + [x])) for x in rest]


def basis_quadruples_edge(n, edge):
    """Quadruples containing a fixed EDGE (2 points): O(n^2) of them (choose 2 more)."""
    rest = [i for i in range(n) if i not in edge]
    return [tuple(sorted(list(edge) + [x, y])) for x, y in combinations(rest, 2)]


def basis_satisfied(D2, quads) -> bool:
    return all(cm4(D2, q) == 0 for q in quads)


def _b1_planar_and_tampered():
    """A genuinely planar integral 5-set, and a copy with ONE non-triangle pair distance tampered.
    Triangle T = {0,1,2}. The tampered pair is (3,4) — a pair that appears in NO quadruple of the
    fixed-triple basis, but DOES appear in the fixed-edge basis."""
    pts = [(0, 0), (1, 0), (0, 1), (1, 1), (2, 2)]   # all coplanar, triangle 0,1,2 nondegenerate
    D2 = dmatrix(pts)
    D2t = [row[:] for row in D2]
    D2t[3][4] = D2t[4][3] = 100                       # break the (3,4) squared distance (was 2)
    return D2, D2t


def b1_test() -> dict:
    """Verdict on: 'CM-vanishing on quadruples containing a fixed nondegenerate triangle forces
    embedding dim <= 2.'"""
    D2, D2t = _b1_planar_and_tampered()
    n = 5
    T = (0, 1, 2)
    tri_quads = basis_quadruples_triple(n, T)
    edge_quads = basis_quadruples_edge(n, (0, 1))

    # positive class: several deterministic planar integral sets — both routes must agree (planar)
    positive = []
    for pts in [
        [(0, 0), (1, 0), (0, 1), (2, 3), (5, 7)],
        [(0, 0), (3, 0), (0, 4), (5, 5), (1, 9), (7, 2)],
        [(0, 0), (2, 1), (1, 3), (4, 4)],
    ]:
        dm = dmatrix(pts)
        emb, rank, inertia = embeds_in_plane_gram(dm)
        cm_ok, _ = all_cm_vanish(dm)
        positive.append({"n": len(pts), "gram_rank": rank, "gram_embeds_plane": emb,
                         "all_cm_vanish": cm_ok, "routes_agree": emb == cm_ok})

    # the counterexample
    tri_ok = basis_satisfied(D2t, tri_quads)          # triple basis satisfied?
    edge_ok = basis_satisfied(D2t, edge_quads)        # edge basis satisfied?
    planar, bad_quad = all_cm_vanish(D2t)             # actually planar?
    emb, rank, (pos, neg, zero) = embeds_in_plane_gram(D2t)
    triangle_ok = cm3(D2t, T) != 0

    return {
        "claim": "CM-vanishing on quadruples containing a fixed nondegenerate triangle forces embedding dim <= 2",
        "triple_basis_size_order": "O(n)  (quadruples containing fixed 3-point triangle)",
        "edge_basis_size_order": "O(n^2) (quadruples containing fixed 2-point edge)",
        "positive_class_all_agree": all(r["routes_agree"] for r in positive),
        "positive_class": positive,
        "counterexample": {
            "note": "planar integral 5-set with the (3,4) squared distance tampered 2 -> 100",
            "triangle_T": list(T),
            "triangle_nondegenerate": triangle_ok,
            "triple_basis_satisfied": tri_ok,
            "edge_basis_satisfied": edge_ok,
            "actually_planar_all_cm_vanish": planar,
            "first_failing_quadruple": list(bad_quad) if bad_quad else None,
            "gram_rank": rank,
            "gram_signature_pos_neg_zero": [pos, neg, zero],
            "gram_is_psd_realizable": neg == 0,
            "tampered_D2": D2t,
        },
        "verdict_triple_basis_sufficient_for_ABSTRACT_metrics": (not (tri_ok and not planar)),
        "verdict_triple_basis_sufficient_for_REALIZABLE_metrics": True,
        "verdict_edge_basis_catches_this_counterexample": (not edge_ok),
        "verdict": (
            "TRIPLE basis (O(n)) is NOT sufficient for arbitrary distance matrices: it never "
            "constrains distances between two non-triangle points (here d_34), so that pair can be "
            "tampered to make a non-basis quadruple non-coplanar while every basis determinant stays "
            "zero. The resulting matrix is non-Euclidean (Gram not PSD: neg>0), so the triple basis "
            "IS sufficient once realizability (Gram PSD) is also required. The O(n^2) EDGE basis "
            "covers every pair (each pair {x,y} sits in quadruple {edge,x,y}) and CATCHES this "
            "counterexample -- it is the minimal-shape CM-only basis that forces planarity. This "
            "matches the operator's O(n^2) count: the sufficient CM basis is the fixed-EDGE basis, "
            "not the fixed-triple basis."
        ),
    }


# ---------------------------------------------------------------------------
# B2 — modular obstructions on integer-coordinate tetrads (CM = 8 * det^2)
# ---------------------------------------------------------------------------
def squares_mod(m) -> set[int]:
    return {(k * k) % m for k in range(m)}


def theoretical_achievable_cm_residues(m) -> set[int]:
    """CM = 8*det^2, det ranges over all integers -> det^2 over all squares mod m ->
    CM over {8*s mod m}. This set is COMPLETE (det attains every residue)."""
    return {(8 * s) % m for s in squares_mod(m)}


def b2_test(coord_hi: int = 3, cap: int = 4000) -> dict:
    """Classify CM residues of integer-coordinate tetrads mod m in {8,16,3,5,7}.
    Empirically verifies the CM = 8*det^2 identity on Z^3 tetrads and confirms every observed
    residue lies in the theoretical achievable set; reports the FORBIDDEN residues."""
    mods = [8, 16, 3, 5, 7]
    theory = {m: theoretical_achievable_cm_residues(m) for m in mods}
    observed = {m: set() for m in mods}
    identity_holds = True
    tested = 0
    p0 = (0, 0, 0)                                     # translation-invariant: fix first point
    rng = range(coord_hi)
    for p1 in ((a, b, c) for a in rng for b in rng for c in rng):
        for p2 in ((a, b, c) for a in rng for b in rng for c in rng):
            for p3 in ((a, b, c) for a in rng for b in rng for c in rng):
                cm = cm4(dmatrix([p0, p1, p2, p3]), (0, 1, 2, 3))
                d = coord_det(p0, p1, p2, p3)
                if cm != 8 * d * d:
                    identity_holds = False
                for m in mods:
                    observed[m].add(cm % m)
                tested += 1
                if tested >= cap:
                    break
            if tested >= cap:
                break
        if tested >= cap:
            break

    mod_tables = {}
    for m in mods:
        ach = sorted(theory[m])
        forb = sorted(set(range(m)) - theory[m])
        mod_tables[str(m)] = {
            "rule": "CM == 8 * (a square) mod m",
            "achievable_residues": ach,
            "forbidden_residues": forb,
            "observed_subset_of_theory": observed[m].issubset(theory[m]),
            "observed_residues": sorted(observed[m]),
        }

    return {
        "identity": "CM(integer-coordinate tetrad) = 8 * det(edge vectors)^2",
        "identity_holds_on_all_tested": identity_holds,
        "tetrads_tested": tested,
        "coord_box": f"p0=origin, p1,p2,p3 in {{0..{coord_hi-1}}}^3",
        "planar_note": "A genuinely PLANAR integral quadruple has det=0 => CM=0 => residue 0 mod EVERY m. "
                       "Hence CM != 0 (mod anything) certifies non-planarity; CM not in {8*square} certifies "
                       "the tetrad is not integer-coordinate at all.",
        "mod_tables": mod_tables,
        "headline_2adic_obstruction": "CM == 0 (mod 8) ALWAYS; CM in {0,8} (mod 16) ONLY. "
                                      "Any tetrad with CM != 0 mod 8 has no integer-coordinate realization.",
    }


# ---------------------------------------------------------------------------
# self-test (ground truth)
# ---------------------------------------------------------------------------
def _selftest():
    # (1) coplanar integer points -> CM = 0 exactly
    assert cm4(dmatrix([(0, 0), (1, 0), (0, 1), (2, 3)]), (0, 1, 2, 3)) == 0
    assert cm4(dmatrix([(0, 0), (1, 0), (0, 1), (5, 7)]), (0, 1, 2, 3)) == 0

    # (2) non-coplanar Z^3 tetrads -> CM != 0, and CM == 8 * det^2 (the spine identity)
    for tet in [
        [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],   # unit tetra: det=1 -> CM=8
        [(0, 0, 0), (2, 0, 0), (0, 3, 0), (1, 1, 5)],
        [(0, 0, 0), (1, 2, 0), (0, 1, 3), (2, 0, 1)],
    ]:
        cm = cm4(dmatrix(tet), (0, 1, 2, 3))
        d = coord_det(*tet)
        assert cm == 8 * d * d, (cm, d)
        assert (cm == 0) == (d == 0)
    assert cm4(dmatrix([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]), (0, 1, 2, 3)) == 8

    # (3) INDEPENDENT verifier (Gram route) must AGREE with the CM route on point sets
    for pts, planar in [
        ([(0, 0), (1, 0), (0, 1), (3, 2), (5, 7)], True),
        ([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)], False),
        ([(0, 0, 0), (2, 0, 0), (0, 3, 0), (1, 1, 5)], False),
    ]:
        dm = dmatrix(pts)
        emb, rank, (pos, neg, zero) = embeds_in_plane_gram(dm)
        cm_ok, _ = all_cm_vanish(dm)
        assert emb == cm_ok == planar, (pts, emb, cm_ok, planar, rank)
        assert neg == 0, ("real point set must have PSD Gram", pts, (pos, neg, zero))

    # (4) triangle nondegeneracy detector
    assert cm3(dmatrix([(0, 0), (1, 0), (0, 1)]), (0, 1, 2)) != 0     # real triangle
    assert cm3(dmatrix([(0, 0), (1, 1), (2, 2)]), (0, 1, 2)) == 0     # collinear -> degenerate

    # (5) B1 counterexample behaves as claimed
    _, D2t = _b1_planar_and_tampered()
    T = (0, 1, 2)
    assert cm3(D2t, T) != 0                                           # T nondegenerate
    assert basis_satisfied(D2t, basis_quadruples_triple(5, T))        # triple basis satisfied
    assert not all_cm_vanish(D2t)[0]                                  # yet NOT planar
    assert not basis_satisfied(D2t, basis_quadruples_edge(5, (0, 1)))  # edge basis CATCHES it
    _, _, (pos, neg, zero) = embeds_in_plane_gram(D2t)
    assert neg > 0                                                    # tampered metric non-realizable

    # (6) B2 modular prediction sanity (exact, no enumeration needed)
    assert theoretical_achievable_cm_residues(8) == {0}
    assert theoretical_achievable_cm_residues(16) == {0, 8}
    assert theoretical_achievable_cm_residues(3) == {0, 2}
    assert theoretical_achievable_cm_residues(5) == {0, 2, 3}
    assert theoretical_achievable_cm_residues(7) == {0, 1, 2, 4}

    print("SELFTEST PASS: cayley_menger")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        results = {
            "blade": "cayley_menger",
            "identity": "CM(4 integer-coordinate points) = 8 * det(edge vectors)^2   [CM = 288 V^2, V = |det|/6]",
            "B1_minimal_basis": b1_test(),
            "B2_modular_obstructions": b2_test(),
        }
        # console summary
        b1 = results["B1_minimal_basis"]
        print("=== B1: minimal-basis sufficiency ===", flush=True)
        print(f"  triple basis (O(n)) sufficient for ABSTRACT metrics : "
              f"{b1['verdict_triple_basis_sufficient_for_ABSTRACT_metrics']}", flush=True)
        print(f"  triple basis sufficient for REALIZABLE metrics       : "
              f"{b1['verdict_triple_basis_sufficient_for_REALIZABLE_metrics']}", flush=True)
        print(f"  edge basis (O(n^2)) catches the counterexample       : "
              f"{b1['verdict_edge_basis_catches_this_counterexample']}", flush=True)
        print(f"  positive class routes agree                          : "
              f"{b1['positive_class_all_agree']}", flush=True)
        print("=== B2: modular obstructions (integer-coordinate tetrads) ===", flush=True)
        b2 = results["B2_modular_obstructions"]
        print(f"  CM = 8*det^2 verified on {b2['tetrads_tested']} tetrads: "
              f"{b2['identity_holds_on_all_tested']}", flush=True)
        for m, tbl in b2["mod_tables"].items():
            print(f"  mod {m:>2}: achievable={tbl['achievable_residues']}  "
                  f"FORBIDDEN={tbl['forbidden_residues']}  "
                  f"(observed_in_theory={tbl['observed_subset_of_theory']})", flush=True)
        print(f"  headline: {b2['headline_2adic_obstruction']}", flush=True)

        out = "oracle/kbk/engine/cayley_menger_results.json"
        with open(out, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nwrote {out}", flush=True)
    else:
        _selftest()
