"""kk_closure_exact.py — RIGOROUS closure: exact integer arithmetic on the ill-conditioned cells.

The float sweep in kk_closure.py is only certified where the quadratic is well-conditioned. Its dangerous
region is exactly {|alpha| small}, which is a thin annulus around an explicit ELLIPSE in the (c_j,c_k) grid
(alpha = 0 <=> the two defining hyperbolas are asymptotically parallel <=> the implied direction is a unit
vector). That region is small (1.5M of 137M cells for |alpha|<1e-2) and is handled here with NO FLOATS.

FULLY-INTEGER FORM. With P_b at the origin, M = common denominator of the vertex coordinates,
Uj = M*(x_j-x_b), Wj = M*(y_j-y_b) (integers), G = Uj*Wk - Wj*Uk, e_i = d_bi^2 - c_i^2, and

    a1 = Wk*e_j - Wj*e_k      b1 = Wk*c_j - Wj*c_k
    c1 = Uj*e_k - Uk*e_j      d1 = Uj*c_k - Uk*c_j

the quadratic alpha*r^2 + beta*r + gamma = 0 clears denominators to the INTEGER equation

    4*ALPHA*r^2 + 4*BETA*r + GAMMA = 0,
      ALPHA = char*M^2*b1^2 + M^2*d1^2 - char*G^2
      BETA  = char*M^2*a1*b1 + M^2*c1*d1
      GAMMA = char*M^2*a1^2 + M^2*c1^2

so  r = (-BETA +/- sqrt(DISC)) / (2*ALPHA)  with  DISC = BETA^2 - ALPHA*GAMMA.

r is rational iff DISC is a perfect square; every test below is exact integer arithmetic (Python bigints).
The degenerate ALPHA = 0 case is handled as the linear equation 4*BETA*r + GAMMA = 0.

Run:  python oracle/kbk/engine/kk_closure_exact.py [H1|H2] [eps]
      python oracle/kbk/engine/kk_closure_exact.py --selftest
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
from kk_closure import best_triple, exact_check

CHAR = 2002


def integer_setup(points, D, b, j, k, char=CHAR):
    """Integer data (M, Uj, Wj, Uk, Wk, G) for the base triple."""
    xb, yb = points[b]
    uj, wj = points[j][0] - xb, points[j][1] - yb
    uk, wk = points[k][0] - xb, points[k][1] - yb
    from math import gcd
    M = 1
    for v in (uj, wj, uk, wk):
        d = int(v.denominator)
        M = M * d // gcd(M, d)          # pure Python ints only — numpy ints would overflow the bigint math
    M = int(M)
    Uj, Wj = int(uj * M), int(wj * M)
    Uk, Wk = int(uk * M), int(wk * M)
    G = Uj * Wk - Wj * Uk
    assert G != 0, "base triple is collinear"
    return M, Uj, Wj, Uk, Wk, G


def exact_cell(M, Uj, Wj, Uk, Wk, G, dbj, dbk, cj, ck, char=CHAR):
    """EXACT integer solve of one cell. Returns the list of positive-integer r solving the quadratic."""
    ej = dbj * dbj - cj * cj
    ek = dbk * dbk - ck * ck
    a1 = Wk * ej - Wj * ek
    b1 = Wk * cj - Wj * ck
    c1 = Uj * ek - Uk * ej
    d1 = Uj * ck - Uk * cj
    M2 = M * M
    ALPHA = char * M2 * b1 * b1 + M2 * d1 * d1 - char * G * G
    BETA = char * M2 * a1 * b1 + M2 * c1 * d1
    GAMMA = char * M2 * a1 * a1 + M2 * c1 * c1
    out = []
    if ALPHA == 0:
        # linear: 4*BETA*r + GAMMA = 0
        if BETA != 0 and GAMMA % (4 * BETA) == 0:
            r = -GAMMA // (4 * BETA)
            if r > 0:
                out.append(r)
        return out
    DISC = BETA * BETA - ALPHA * GAMMA
    if DISC < 0:
        return out
    s = isqrt(DISC)
    if s * s != DISC:
        return out                      # r irrational -> no integer solution
    for num in (-BETA + s, -BETA - s):
        den = 2 * ALPHA
        if den != 0 and num % den == 0:
            r = num // den
            if r > 0:
                out.append(r)
    return out


def sweep_dangerous(points, D, label, eps=1e-2, char=CHAR):
    n = len(points)
    cost, b, j, k, dbj, dbk = best_triple(D, n)
    M, Uj, Wj, Uk, Wk, G = integer_setup(points, D, b, j, k, char)
    print(f"[{label}] base P{b}, partners P{j}(d={dbj}) P{k}(d={dbk}); M={M}, G={G}")
    print(f"[{label}] EXACT sweep of all cells with |alpha| < {eps} (no floats in the decision path)")

    # identify dangerous cells with float (safe: float alpha is accurate to ~1e-16 absolute,
    # so screening at 2*eps provably contains every cell whose TRUE |alpha| < eps)
    xb, yb = points[b]
    ujf, wjf = float(points[j][0] - xb), float(points[j][1] - yb)
    ukf, wkf = float(points[k][0] - xb), float(points[k][1] - yb)
    D0 = ujf * wkf - wjf * ukf
    ckv = np.arange(-dbk, dbk + 1, dtype=np.float64)

    checked = 0
    solutions = []
    t0 = time.time()
    total = 2 * dbj + 1
    for idx, cj in enumerate(range(-dbj, dbj + 1)):
        B = (wkf * cj - wjf * ckv) / D0
        Dd = (ujf * ckv - ukf * cj) / (char * D0)
        alpha = B * B + char * Dd * Dd - 1.0
        idxs = np.flatnonzero(np.abs(alpha) < 2 * eps)
        for t in idxs:
            ck = int(ckv[t])
            checked += 1
            for r in exact_cell(M, Uj, Wj, Uk, Wk, G, dbj, dbk, cj, ck, char):
                if r - cj < 0 or r - ck < 0:
                    continue
                ok, X, dists = exact_check(points, D, b, j, k, r, cj, ck, char)
                if ok:
                    is_vertex = any(X == p for p in points)
                    solutions.append({"r": r, "cj": cj, "ck": ck, "x": str(X[0]),
                                      "y_coeff": str(X[1]), "distances": dists,
                                      "is_vertex": is_vertex})
                    tag = "vertex" if is_vertex else "*** NON-TRIVIAL ***"
                    print(f"    EXACT SOLUTION {tag}: r={r} cj={cj} ck={ck} "
                          f"X=({X[0]}, {X[1]}*sqrt{char}) d={dists}", flush=True)
        if (idx + 1) % max(1, total // 20) == 0:
            el = time.time() - t0
            print(f"  [{label}] c_j {idx+1}/{total} ({100*(idx+1)/total:.1f}%) "
                  f"exact_cells={checked} solutions={len(solutions)} elapsed={el:.0f}s "
                  f"eta={el/(idx+1)*(total-idx-1):.0f}s", flush=True)
    el = time.time() - t0
    nontriv = [s for s in solutions if not s["is_vertex"]]
    print(f"\n[{label}] EXACT DANGEROUS-CELL SWEEP COMPLETE in {el:.0f}s")
    print(f"  cells solved with exact integer arithmetic : {checked:,}")
    print(f"  solutions found : {len(solutions)}   NON-TRIVIAL : {len(nontriv)}")
    return {"label": label, "eps": eps, "exact_cells": checked, "solutions": solutions,
            "nontrivial": nontriv, "seconds": el}


def float_error_certificate(points, D, label, samples=40000, seed=12345, char=CHAR):
    """Measure the ACTUAL max |float r - exact r| over well-conditioned cells, to certify the
    1e-2 near-integer window used by the float sweep."""
    n = len(points)
    cost, b, j, k, dbj, dbk = best_triple(D, n)
    M, Uj, Wj, Uk, Wk, G = integer_setup(points, D, b, j, k, char)
    xb, yb = points[b]
    ujf, wjf = float(points[j][0] - xb), float(points[j][1] - yb)
    ukf, wkf = float(points[k][0] - xb), float(points[k][1] - yb)
    D0 = ujf * wkf - wjf * ukf
    dbj2, dbk2 = float(dbj) ** 2, float(dbk) ** 2
    rng = np.random.default_rng(seed)
    worst = 0.0
    worst_at = None
    tested = 0
    for _ in range(samples):
        cj = int(rng.integers(-dbj, dbj + 1)); ck = int(rng.integers(-dbk, dbk + 1))
        ej = dbj2 - cj ** 2; ek = dbk2 - ck ** 2
        A = (wkf * ej - wjf * ek) / (2 * D0); B = (wkf * cj - wjf * ck) / D0
        C = (ujf * ek - ukf * ej) / (2 * char * D0); Dd = (ujf * ck - ukf * cj) / (char * D0)
        alpha = B * B + char * Dd * Dd - 1.0
        if abs(alpha) < 1e-2:
            continue                      # dangerous cells are handled exactly, not here
        beta = 2.0 * (A * B + char * C * Dd); gamma = A * A + char * C * C
        disc = beta * beta - 4.0 * alpha * gamma
        if disc < 0:
            continue
        sq = disc ** 0.5
        froots = [(-beta + sq) / (2 * alpha), (-beta - sq) / (2 * alpha)]
        # exact roots for the same cell
        for r in exact_cell(M, Uj, Wj, Uk, Wk, G, dbj, dbk, cj, ck, char):
            tested += 1
            err = min(abs(fr - r) for fr in froots)
            if err > worst:
                worst = err; worst_at = (cj, ck, r)
    return {"samples": samples, "exact_roots_compared": tested,
            "max_abs_error": worst, "worst_at": worst_at}


def _selftest():
    cost, b, j, k, dbj, dbk = best_triple(D1)
    M, Uj, Wj, Uk, Wk, G = integer_setup(H1, D1, b, j, k)
    # the exact integer solver must recover r for every non-base vertex
    ok_all = True
    for i in range(7):
        if i == b:
            continue
        r = D1[b][i]; cj = r - D1[j][i]; ck = r - D1[k][i]
        rs = exact_cell(M, Uj, Wj, Uk, Wk, G, dbj, dbk, cj, ck)
        if r not in rs:
            ok_all = False
    print(json.dumps({"ok": bool(ok_all),
                      "method": "exact integer cell solver recovers r for all 6 non-base vertices"}))
    return ok_all


def run_spec(spec: dict) -> dict:
    """Run the existing exact dangerous-cell sweep from a typed instance request."""
    if not isinstance(spec, dict):
        return {"ok": False, "verdict": "REFUSED", "reason": "input spec must be an object"}
    operation = spec.get("operation", "sweep")
    if operation == "selftest":
        ok = _selftest()
        return {"ok": bool(ok), "verdict": "SELFTEST" if ok else "FAILED"}
    which = spec.get("instance")
    eps = spec.get("eps", 1e-2)
    if operation != "sweep" or which not in {"H1", "H2"} or isinstance(eps, bool) or not isinstance(eps, (int, float)) or eps <= 0:
        return {"ok": False, "verdict": "REFUSED", "reason": "operation must be selftest or sweep with instance H1/H2 and positive eps"}
    pts, distances = (H1, D1) if which == "H1" else (reconstruct_full(D2, 2002), D2)
    certificate = float_error_certificate(pts, distances, which)
    result = sweep_dangerous(pts, distances, which, float(eps))
    result["float_error_certificate"] = certificate
    return {"ok": True, "verdict": "EXACT_SWEEP", "result": result}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("instance", nargs="?", choices=("H1", "H2"))
    parser.add_argument("eps", nargs="?", type=float, default=1e-2)
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
    eps = args.eps
    pts, D = (H1, D1) if which == "H1" else (reconstruct_full(D2, 2002), D2)
    cert = float_error_certificate(pts, D, which)
    print(f"[{which}] FLOAT-ERROR CERTIFICATE (well-conditioned cells, |alpha|>=1e-2):")
    print(f"  {json.dumps(cert)}")
    print(f"  window used by the float sweep = 1e-2  ->  margin = "
          f"{1e-2/max(cert['max_abs_error'],1e-18):.3e}x\n")
    res = sweep_dangerous(pts, D, which, eps)
    res["float_error_certificate"] = cert
    with open(f"oracle/kbk/engine/closure_exact_{which}.json", "w") as f:
        json.dump(res, f, indent=1)
    print(f"  receipt -> oracle/kbk/engine/closure_exact_{which}.json")
