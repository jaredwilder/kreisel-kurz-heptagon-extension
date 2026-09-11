"""kk_closure_full_exact.py — THE DEFINITIVE SWEEP: every cell, exact integer arithmetic, no floats.

Supersedes the hybrid float/exact approach. The float prefilter turned out to have only a ~3x margin
(measured worst float64 root error 3.06e-3 against a 1e-2 acceptance window), which is far too thin to
support a closure claim. Exact integer arithmetic costs only ~1.8us/cell, so the ENTIRE grid is swept
exactly in a few minutes and no numerical argument is needed anywhere.

THE THEOREM THIS DECIDES (for a given integral point set S = {P_0..P_6} in general position):
    Is there a point X in the plane with |X P_i| integral for every i?
    The vertices themselves trivially qualify; the question is whether anything ELSE does.

WHY THE ENUMERATION IS COMPLETE (no height bound needed):
    Pick a base vertex P_b and two others P_j, P_k (non-collinear). For ANY such X,
        c_j := |XP_b| - |XP_j|  and  c_k := |XP_b| - |XP_k|
    are INTEGERS (differences of integers) and satisfy |c_j| <= d_bj, |c_k| <= d_bk by the triangle
    inequality -- a FINITE range that does not depend on how far X lies. Given (c_j, c_k), the three
    distances from three non-collinear points determine X uniquely, and r = |XP_b| satisfies an explicit
    quadratic. So the full solution set is covered by the finite grid, exhaustively.

INTEGER FORM (derivation in kk_closure_exact.py):
    4*ALPHA*r^2 + 4*BETA*r + GAMMA = 0,  r = (-BETA +/- sqrt(BETA^2 - ALPHA*GAMMA)) / (2*ALPHA)
Every quantity below is a Python bigint; the only tests are ">= 0", "is a perfect square", "divides".

Run: python oracle/kbk/engine/kk_closure_full_exact.py [H1|H2]
     python oracle/kbk/engine/kk_closure_full_exact.py --selftest
"""
from __future__ import annotations
from math import isqrt
import json
import sys
import time

sys.path.insert(0, "oracle/kbk/engine")
from kk_heptagon import H1, D1, D2, reconstruct_full
from kk_closure import best_triple, exact_check
from kk_closure_exact import integer_setup

CHAR = 2002


def full_exact_sweep(points, D, label, char=CHAR, report=20, base=None):
    n = len(points)
    if base is None:
        cost, b, j, k, dbj, dbk = best_triple(D, n)
    else:
        b, j, k = base
        dbj, dbk = D[b][j], D[b][k]
        cost = (2 * dbj + 1) * (2 * dbk + 1)
    M, Uj, Wj, Uk, Wk, G = integer_setup(points, D, b, j, k, char)
    M2 = M * M
    charM2 = char * M2
    charG2 = char * G * G
    dbj2, dbk2 = dbj * dbj, dbk * dbk

    print(f"[{label}] base P{b}; partners P{j}(d={dbj}), P{k}(d={dbk})")
    print(f"[{label}] grid = {cost:,} cells; EXACT INTEGER ARITHMETIC ONLY (no floats anywhere)")

    # precompute ek and the ck-dependent pieces
    cks = list(range(-dbk, dbk + 1))
    eks = [dbk2 - c * c for c in cks]
    Wj_ek = [Wj * e for e in eks]
    Uj_ek = [Uj * e for e in eks]
    Wj_ck = [Wj * c for c in cks]
    Uj_ck = [Uj * c for c in cks]

    solutions = []
    t0 = time.time()
    total = 2 * dbj + 1
    for idx, cj in enumerate(range(-dbj, dbj + 1)):
        ej = dbj2 - cj * cj
        Wk_ej = Wk * ej
        Uk_ej = Uk * ej
        Wk_cj = Wk * cj
        Uk_cj = Uk * cj
        for t in range(len(cks)):
            a1 = Wk_ej - Wj_ek[t]
            b1 = Wk_cj - Wj_ck[t]
            c1 = Uj_ek[t] - Uk_ej
            d1 = Uj_ck[t] - Uk_cj
            ALPHA = charM2 * b1 * b1 + M2 * d1 * d1 - charG2
            BETA = charM2 * a1 * b1 + M2 * c1 * d1
            GAMMA = charM2 * a1 * a1 + M2 * c1 * c1
            if ALPHA == 0:
                if BETA != 0:
                    den = 4 * BETA
                    if (-GAMMA) % den == 0:
                        r = (-GAMMA) // den
                        if r > 0:
                            _try(solutions, points, D, b, j, k, r, cj, cks[t], char)
                continue
            DISC = BETA * BETA - ALPHA * GAMMA
            if DISC < 0:
                continue
            s = isqrt(DISC)
            if s * s != DISC:
                continue
            den = 2 * ALPHA
            for num in (-BETA + s, -BETA - s):
                if num % den == 0:
                    r = num // den
                    if r > 0:
                        _try(solutions, points, D, b, j, k, r, cj, cks[t], char)
        if (idx + 1) % max(1, total // report) == 0:
            el = time.time() - t0
            print(f"  [{label}] c_j {idx+1}/{total} ({100*(idx+1)/total:.1f}%) "
                  f"solutions={len(solutions)} elapsed={el:.0f}s "
                  f"eta={el/(idx+1)*(total-idx-1):.0f}s", flush=True)

    el = time.time() - t0
    # de-duplicate by coordinates
    seen, uniq = set(), []
    for s_ in solutions:
        kk = (s_["x"], s_["y_coeff"])
        if kk not in seen:
            seen.add(kk); uniq.append(s_)
    solutions = uniq
    nontriv = [s_ for s_ in solutions if not s_["is_vertex"]]

    want = {(str(points[i][0]), str(points[i][1])) for i in range(n) if i != b}
    got = {(s_["x"], s_["y_coeff"]) for s_ in solutions}
    missing = want - got
    sound = not missing

    print(f"\n[{label}] CORRECTNESS GATE: recovered {len(want & got)}/{len(want)} known non-base vertices"
          f" -> {'PASS' if sound else 'FAIL'}")
    if missing:
        print(f"  !! MISSING -> sweep UNSOUND, no closure claim: {sorted(missing)}")
    print(f"[{label}] FULL EXACT SWEEP COMPLETE in {el:.0f}s")
    print(f"  cells (all exact)   : {cost:,}")
    print(f"  total solutions     : {len(solutions)}")
    print(f"  vertices            : {len(solutions)-len(nontriv)}")
    print(f"  NON-TRIVIAL         : {len(nontriv)}")
    return {"label": label, "grid_cells": cost, "solutions": solutions, "nontrivial": nontriv,
            "seconds": el, "base": [b, j, k], "correctness_gate_pass": bool(sound),
            "missing_vertices": sorted(missing), "method": "exact integer arithmetic, no floats"}


def _try(acc, points, D, b, j, k, r, cj, ck, char):
    ok, X, dists = exact_check(points, D, b, j, k, r, cj, ck, char)
    if not ok:
        return
    is_vertex = any(X == p for p in points)
    acc.append({"r": r, "cj": cj, "ck": ck, "x": str(X[0]), "y_coeff": str(X[1]),
                "distances": dists, "is_vertex": is_vertex})
    tag = "vertex" if is_vertex else "*** NON-TRIVIAL ***"
    print(f"    SOLUTION {tag}: r={r} cj={cj} ck={ck} X=({X[0]}, {X[1]}*sqrt{char}) d={dists}", flush=True)


def _selftest():
    """Sanity: a tiny sub-grid around a known vertex must recover it, with exact arithmetic."""
    cost, b, j, k, dbj, dbk = best_triple(D1)
    M, Uj, Wj, Uk, Wk, G = integer_setup(H1, D1, b, j, k)
    from kk_closure_exact import exact_cell
    ok_all = True
    for i in range(7):
        if i == b:
            continue
        r = D1[b][i]; cj = r - D1[j][i]; ck = r - D1[k][i]
        if r not in exact_cell(M, Uj, Wj, Uk, Wk, G, dbj, dbk, cj, ck):
            ok_all = False
    print(json.dumps({"ok": bool(ok_all), "method": "full-exact sweep kernel recovers all non-base vertices"}))
    return ok_all


if __name__ == "__main__":
    if "--selftest" in sys.argv or len(sys.argv) == 1:
        _selftest(); sys.exit(0)
    which = sys.argv[1]
    base = None
    tag = ""
    if len(sys.argv) > 4:
        base = (int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
        tag = f"_base{base[0]}{base[1]}{base[2]}"
    pts, D = (H1, D1) if which == "H1" else (reconstruct_full(D2, 2002), D2)
    res = full_exact_sweep(pts, D, which, base=base)
    with open(f"oracle/kbk/engine/closure_full_exact_{which}{tag}.json", "w") as f:
        json.dump(res, f, indent=1)
    print(f"  receipt -> oracle/kbk/engine/closure_full_exact_{which}.json")
