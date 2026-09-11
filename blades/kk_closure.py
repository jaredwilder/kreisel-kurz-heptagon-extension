"""kk_closure.py — THE CLOSURE ALGORITHM: exhaustive, unconditional, no height bound.

THE IDEA THAT REMOVES THE HEIGHT BOUND. Every previous search fixed a box r,s <= B and could only ever
give a height-BOUNDED receipt, because nothing stopped the 8th point X from lying farther out. But:

    for ANY point X and any two vertices P_b, P_j, the triangle inequality gives
        | |XP_b| - |XP_j| |  <=  d_bj
    and if both distances are integers, c_j := |XP_b| - |XP_j| is an INTEGER in [-d_bj, d_bj].

That range is FINITE and does not depend on how far away X is. Choosing a base vertex P_b and two others
P_j, P_k (non-collinear), the pair (c_j, c_k) ranges over a finite grid. And three distances from three
non-collinear points determine X UNIQUELY (trilateration). So:

    { X : all distances to P_b,P_j,P_k integral }  is enumerated EXACTLY by the finite grid (c_j,c_k),
    each cell contributing at most 2 points (a quadratic in r = |XP_b|).

Hence the full solution set is finite AND effectively enumerable -- no height bound, no conjecture.

DERIVATION (exact, in Q(sqrt(char)); points are (x, v*sqrt(char)) relative to P_b at the origin).
With u_i = x_i - x_b, w_i = y_i - y_b, e_i = d_bi^2 - c_i^2, and Delta0 = u_j w_k - w_j u_k != 0:

    2 u_j x + 2*char* w_j v = e_j + 2 r c_j          (from |X-P_j|^2 = (r-c_j)^2, using x^2+char*v^2 = r^2)
    2 u_k x + 2*char* w_k v = e_k + 2 r c_k

    =>  x = A + B r,   v = C + D r      with
        A = (w_k e_j - w_j e_k)/(2*Delta0)          B = (w_k c_j - w_j c_k)/Delta0
        C = (u_j e_k - u_k e_j)/(2*char*Delta0)     D = (u_j c_k - u_k c_j)/(char*Delta0)

    and substituting into x^2 + char*v^2 = r^2 gives the QUADRATIC IN r:
        alpha r^2 + beta r + gamma = 0,
        alpha = B^2 + char*D^2 - 1,  beta = 2(A B + char*C D),  gamma = A^2 + char*C^2.

Two-stage for speed WITHOUT losing exactness: float64 vectorised prefilter (keep r within 1e-3 of an
integer -- float error here is ~1e-6, so nothing real can be missed), then EXACT Fraction verification
of every survivor. A solution is reported only if exact arithmetic confirms all 7 distances are integers.

Run:  python oracle/kbk/engine/kk_closure.py [H1|H2] [--selftest]
"""
from __future__ import annotations
from fractions import Fraction as F
from math import isqrt
import json
import sys
import time
import argparse
from pathlib import Path
import numpy as np

sys.path.insert(0, "oracle/kbk/engine")
from kk_heptagon import H1, D1, D2, reconstruct_full

CHAR = 2002


def best_triple(D, n=7):
    """Base vertex b and two others j,k minimising the (c_j,c_k) grid size."""
    best = None
    for b in range(n):
        ds = sorted((D[b][i], i) for i in range(n) if i != b)
        (dj, j), (dk, k) = ds[0], ds[1]
        cost = (2 * dj + 1) * (2 * dk + 1)
        if best is None or cost < best[0]:
            best = (cost, b, j, k, dj, dk)
    return best


def exact_check(points, D, b, j, k, r, cj, ck, char=CHAR):
    """Exact rational verification: rebuild X from (r,cj,ck) and test ALL distances are integers.
    Returns (ok, X, distances) with X = (x, v) meaning (x, v*sqrt(char)) in ABSOLUTE coordinates."""
    xb, yb = points[b]
    uj, wj = points[j][0] - xb, points[j][1] - yb
    uk, wk = points[k][0] - xb, points[k][1] - yb
    D0 = uj * wk - wj * uk
    if D0 == 0:
        return False, None, None
    ej = F(D[b][j]) ** 2 - F(cj) ** 2
    ek = F(D[b][k]) ** 2 - F(ck) ** 2
    A = (wk * ej - wj * ek) / (2 * D0)
    B = (wk * F(cj) - wj * F(ck)) / D0
    C = (uj * ek - uk * ej) / (2 * char * D0)
    Dd = (uj * F(ck) - uk * F(cj)) / (char * D0)
    x = A + B * F(r)
    v = C + Dd * F(r)
    # the defining constraint must hold EXACTLY
    if x * x + char * v * v != F(r) ** 2:
        return False, None, None
    X = (x + xb, v + yb)
    dists = []
    for i in range(len(points)):
        q = (X[0] - points[i][0]) ** 2 + char * (X[1] - points[i][1]) ** 2
        if q.denominator != 1:
            return False, None, None
        rt = isqrt(q.numerator)
        if rt * rt != q.numerator:
            return False, None, None
        dists.append(rt)
    return True, X, dists


def closure(points, D, label, char=CHAR, chunk_report=25):
    n = len(points)
    cost, b, j, k, dbj, dbk = best_triple(D, n)
    print(f"[{label}] base P{b}; partners P{j} (d={dbj}), P{k} (d={dbk})")
    print(f"[{label}] EXHAUSTIVE grid: c_j in [-{dbj},{dbj}] x c_k in [-{dbk},{dbk}] = {cost:,} cells")
    print(f"[{label}] each cell -> quadratic in r -> <=2 candidate points. NO height bound involved.")

    xb, yb = points[b]
    uj, wj = float(points[j][0] - xb), float(points[j][1] - yb)
    uk, wk = float(points[k][0] - xb), float(points[k][1] - yb)
    D0 = uj * wk - wj * uk
    dbj2, dbk2 = float(D[b][j]) ** 2, float(D[b][k]) ** 2

    ck = np.arange(-dbk, dbk + 1, dtype=np.float64)
    ek = dbk2 - ck * ck
    solutions = []
    exact_tests = 0
    t0 = time.time()
    total_cj = 2 * dbj + 1

    for idx, cj in enumerate(range(-dbj, dbj + 1)):
        ej = dbj2 - float(cj) ** 2
        A = (wk * ej - wj * ek) / (2 * D0)
        B = (wk * float(cj) - wj * ck) / D0
        C = (uj * ek - uk * ej) / (2 * char * D0)
        Dd = (uj * ck - uk * float(cj)) / (char * D0)
        alpha = B * B + char * Dd * Dd - 1.0
        beta = 2.0 * (A * B + char * C * Dd)
        gamma = A * A + char * C * C
        t1 = beta * beta
        t2 = 4.0 * alpha * gamma
        disc = t1 - t2
        # TANGENCY GUARD: disc can be EXACTLY 0 mathematically (double root, e.g. X = P_j) yet come out
        # slightly negative in float64 via catastrophic cancellation. Reject only if disc is negative
        # beyond a generous relative tolerance; otherwise clamp to 0. (Observed real error ~1e-9 on
        # magnitudes ~1e5, i.e. ~1e-14 relative; the 1e-9 relative tolerance below is a 10^5 margin.)
        scale = np.abs(t1) + np.abs(t2) + 1.0
        with np.errstate(invalid='ignore', divide='ignore'):
            ok_disc = disc >= -1e-9 * scale
            sq = np.sqrt(np.clip(disc, 0.0, None))
            # numerically stable roots: q = -(beta + sign(beta)*sq)/2 ; roots are q/alpha and gamma/q
            sgn = np.where(beta >= 0, 1.0, -1.0)
            qq = -0.5 * (beta + sgn * sq)
            roots = []
            with np.errstate(divide='ignore', invalid='ignore'):
                roots.append(np.where(np.abs(alpha) > 0, qq / np.where(alpha == 0, np.nan, alpha), np.nan))
                roots.append(np.where(np.abs(qq) > 0, gamma / np.where(qq == 0, np.nan, qq), np.nan))
                # DEGENERATE LINEAR CASE alpha ~ 0: beta*r + gamma = 0
                roots.append(np.where(np.abs(beta) > 0, -gamma / np.where(beta == 0, np.nan, beta), np.nan))
            for r in roots:
                cand = ok_disc & np.isfinite(r) & (r > 0)
                # generous near-integer window (float error ~1e-6; 1e-2 is a 10^4 margin)
                cand &= np.abs(r - np.rint(r)) < 1e-2
                # distances must be non-negative
                cand &= (r - float(cj) >= -1e-6) & (r - ck >= -1e-6)
                hits = np.flatnonzero(cand)
                for h in hits:
                    exact_tests += 1
                    rr = int(round(r[h]))
                    ok, X, dists = exact_check(points, D, b, j, k, rr, cj, int(ck[h]), char)
                    if ok:
                        is_vertex = any(X == p for p in points)
                        solutions.append({"r": rr, "cj": cj, "ck": int(ck[h]),
                                          "x": str(X[0]), "y_coeff": str(X[1]),
                                          "distances": dists, "is_vertex": is_vertex})
                        tag = "vertex" if is_vertex else "*** NON-TRIVIAL ***"
                        print(f"    SOLUTION {tag}: r={rr} cj={cj} ck={int(ck[h])} X=({X[0]}, {X[1]}*sqrt{char}) d={dists}", flush=True)
        if (idx + 1) % max(1, total_cj // chunk_report) == 0:
            el = time.time() - t0
            pct = 100 * (idx + 1) / total_cj
            print(f"  [{label}] c_j {idx+1}/{total_cj} ({pct:.1f}%)  exact_tests={exact_tests}  "
                  f"solutions={len(solutions)}  elapsed={el:.0f}s  eta={el/(idx+1)*(total_cj-idx-1):.0f}s", flush=True)

    el = time.time() - t0
    # de-duplicate (the same point can be reached by more than one root formula)
    seen, uniq = set(), []
    for s in solutions:
        kkey = (s["x"], s["y_coeff"])
        if kkey not in seen:
            seen.add(kkey); uniq.append(s)
    solutions = uniq
    nontriv = [s for s in solutions if not s["is_vertex"]]

    # ---- MANDATORY CORRECTNESS GATE -------------------------------------------------
    # The scan MUST recover every non-base vertex. If it does not, the sweep is unsound
    # and no closure claim may be made from it.
    want = {(str(points[i][0]), str(points[i][1])) for i in range(n) if i != b}
    got = {(s["x"], s["y_coeff"]) for s in solutions}
    missing = want - got
    sound = not missing
    print(f"\n[{label}] CORRECTNESS GATE: recovered {len(want & got)}/{len(want)} known non-base vertices"
          f"  -> {'PASS' if sound else 'FAIL'}")
    if missing:
        print(f"  !! MISSING (sweep is UNSOUND, closure claim INVALID): {sorted(missing)}")
    print(f"\n[{label}] CLOSURE COMPLETE in {el:.0f}s")
    print(f"  grid cells examined : {cost:,}  (EXHAUSTIVE -- no height bound)")
    print(f"  exact verifications : {exact_tests}")
    print(f"  total solutions     : {len(solutions)}  (points at integral distance from ALL {n} vertices)")
    print(f"  of which vertices   : {len(solutions)-len(nontriv)}")
    print(f"  NON-TRIVIAL         : {len(nontriv)}")
    return {"label": label, "grid_cells": cost, "exact_tests": exact_tests,
            "solutions": solutions, "nontrivial": nontriv, "seconds": el,
            "base": [b, j, k], "correctness_gate_pass": bool(sound),
            "missing_vertices": sorted(missing)}


def _selftest():
    """The algorithm MUST rediscover the known trivial solutions X = P_i (i != base)."""
    cost, b, j, k, dbj, dbk = best_triple(D1)
    found = []
    for i in range(7):
        if i == b:
            continue
        r = D1[b][i]; cj = r - D1[j][i]; ck = r - D1[k][i]
        ok, X, dists = exact_check(H1, D1, b, j, k, r, cj, ck)
        found.append((i, ok, X == H1[i] if ok else False))
    allok = all(ok and match for _, ok, match in found)
    print(json.dumps({"ok": bool(allok), "method": "kk_closure exact_check rediscovers all 6 non-base vertices"}))
    return allok


def run_spec(spec: dict) -> dict:
    """Aim the existing no-height-bound closure engine at a declared built-in instance."""
    if not isinstance(spec, dict):
        return {"ok": False, "verdict": "REFUSED", "reason": "input spec must be an object"}
    operation = spec.get("operation", "close")
    if operation == "selftest":
        ok = _selftest()
        return {"ok": bool(ok), "verdict": "SELFTEST" if ok else "FAILED"}
    which = spec.get("instance")
    if operation != "close" or which not in {"H1", "H2"}:
        return {"ok": False, "verdict": "REFUSED", "reason": "operation must be selftest or close with instance H1 or H2"}
    pts, distances = (H1, D1) if which == "H1" else (reconstruct_full(D2, 2002), D2)
    result = closure(pts, distances, which)
    sound = bool(result.get("correctness_gate_pass"))
    return {"ok": sound, "verdict": "CLOSED" if sound else "UNSOUND", "result": result}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("instance", nargs="?", choices=("H1", "H2"))
    parser.add_argument("--selftest", action="store_true")
    parser.add_argument("--input-file")
    parser.add_argument("--out")
    args = parser.parse_args()
    if args.input_file:
        try:
            payload = run_spec(json.loads(Path(args.input_file).read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError) as exc:
            payload = {"ok": False, "verdict": "REFUSED", "reason": str(exc)}
        text = json.dumps(payload)
        if args.out:
            Path(args.out).write_text(text + "\n", encoding="utf-8")
        print(text)
        sys.exit(0 if payload.get("ok") else 2)
    if args.selftest or args.instance is None:
        sys.exit(0 if _selftest() else 1)
    which = args.instance
    pts, D = (H1, D1) if which == "H1" else (reconstruct_full(D2, 2002), D2)
    res = closure(pts, D, which)
    with open(f"oracle/kbk/engine/closure_{which}.json", "w") as f:
        json.dump(res, f, indent=1)
    print(f"  receipt -> oracle/kbk/engine/closure_{which}.json")
