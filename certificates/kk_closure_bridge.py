"""kk_closure_bridge.py — route SHAPE-134 (Kreisel-Kurz / integral-distance CLOSURE) targets into
`kk_closure.py` + `kk_closure_exact.py`, fail-closed through the `gated_adapter` kernel.

THE SHAPE. Given an integral point set S = {P_0..P_{n-1}} in Q(sqrt(char)), is there ANY further point X
in the plane at integral distance from every P_i? The vertices trivially qualify; the question is whether
anything ELSE does. The enumeration is EXHAUSTIVE WITH NO HEIGHT BOUND: for a base vertex P_b and two
others P_j, P_k, c_j = |XP_b| - |XP_j| is an integer in [-d_bj, d_bj] by the triangle inequality (a finite
range independent of how far X lies), and (r, c_j, c_k) fixes X by trilateration. See kk_closure.py.

GRAMMAR
    kkclosure:run,<label>,<D>[,<max_cells>]   sweep the (c_j,c_k) grid for label H1|H2 with matrix D1|D2
    kkclosure:check,<points>,<D>              independently verify a candidate integral point set

    <label>     H1 | H2                 <D>  D1 | D2   (must be consistent with the label)
    <max_cells> a positive integer cell budget, or `full` for the unbounded exhaustive sweep.
                DEFAULT IS A BUDGET, NOT `full` — see the trap below.
    <points>    H1 | H2, or a JSON list [["x","y_coeff"],...] meaning (x, y_coeff*sqrt(char))

⛔ THE TRAP THIS ADAPTER EXISTS TO NOT FALL INTO. The entire value of this result is that it is NOT a
sample. A budgeted run that stopped early is `INCOMPLETE_BUDGET_EXHAUSTED` and is NEVER upgraded to a
closure, no matter what it found (or did not find). Every verdict therefore carries `chunks_completed`
and `chunks_total` (chunks = c_j rows) plus `cells_examined` / `cells_total`, and `claim_strength` says
in words what a partial sweep does and does not license. A truncated run is proved to come back
INCOMPLETE by selftest case `truncated_run_is_INCOMPLETE_not_CLOSED`.

If a complete banked artifact matches the run's own parameters it is surfaced under `banked_reference` —
clearly separated from `verdict`, which describes only what THIS call computed.

IN-ROUTE GROUND-TRUTH GATE (every call, before any work, via gated_adapter.run_gated).
The full exhaustive sweep is MEASURED at 153s (banked, exact integer) / 170s (this machine, float
prefilter) against a 60s gate budget, so the gate cannot BE the exhaustive sweep — that is exactly how
SHAPE-027 shipped inert. What the gate does instead, in ~1s:

  1. the banked artifact `closure_full_exact_H1.json` parses and is a COMPLETE run
     (correctness_gate_pass true, missing_vertices empty);
  2. `best_triple(D1)` re-derives the banked base triple [4,5,6];
  3. every one of the 6 banked solution records is REPRODUCED BYTE-IDENTICALLY by re-running
     `exact_check(H1, D1, b, j, k, r, cj, ck)` and re-serialising the record — not compared field-by-
     field with a tolerance, compared as bytes;
  4. an INDEPENDENT arithmetic path (`kk_closure_exact.exact_cell`, pure bigint, no Fractions)
     recovers the same r for every banked cell;
  5. the banked solution set is exactly the 6 non-base vertices and `nontrivial` is empty;
  6. INDEPENDENT RE-VERIFICATION by `oracle/kernel/checkers.check_integral_point_set(7, 2002, H1)`,
     which imports nothing from the search code;
  7. NEGATIVE CONTROL: a perturbed banked record (r+1) must FAIL exact_check. Without this the gate
     would pass against an exact_check that returned True unconditionally.

The smaller banked artifact `closure_exact_H1.json` (dangerous-cell sweep, 2,935,632 cells) is checked
for shape only — it is a sub-sweep, not a closure, and is never a source of a CLOSED verdict.

$0, deterministic, zero-LLM, NO NETWORK — every number comes from exact integer / Fraction arithmetic.
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
from fractions import Fraction

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
_REPO = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

import gated_adapter as GA

try:
    from kk_heptagon import H1, D1, D2, reconstruct_full
    from kk_closure import best_triple, exact_check, CHAR
    from kk_closure_exact import integer_setup, exact_cell
    _SOLVER_OK, _SOLVER_ERR = True, None
except Exception as exc:  # pragma: no cover - environment guard
    _SOLVER_OK, _SOLVER_ERR = False, f"{type(exc).__name__}: {exc}"
    CHAR = 2002

try:
    from oracle.kernel.checkers import check_integral_point_set
    _CHECKER_OK, _CHECKER_ERR = True, None
except Exception as exc:  # pragma: no cover - environment guard
    _CHECKER_OK, _CHECKER_ERR = False, f"{type(exc).__name__}: {exc}"

GRAMMAR = re.compile(r"^\s*kkclosure\s*:\s*(run|check)\s*,\s*(.+?)\s*$", re.IGNORECASE)

BANKED_FULL = os.path.join(_HERE, "closure_full_exact_H1.json")
BANKED_SUB = os.path.join(_HERE, "closure_exact_H1.json")

# Measured, not estimated: the full exhaustive H1 sweep is 153s banked / 170s re-run on this machine.
FULL_SWEEP_MEASURED_S = {"banked_exact_integer": 153.2, "rerun_float_prefilter": 170.1}
GATE_BUDGET_S = 60.0
DEFAULT_MAX_CELLS = 1_500_000   # a BUDGET by default; `full` must be asked for explicitly


# --------------------------------------------------------------------------------------------------
# grammar
# --------------------------------------------------------------------------------------------------
def is_kk_closure_target(query: str) -> bool:
    return bool(GRAMMAR.match(query or ""))


def _load_set(label: str):
    """(points, D) for a label. H2 is reconstructed from its distance matrix, as kk_closure does."""
    lab = label.strip().upper()
    if lab == "H1":
        return list(H1), D1
    if lab == "H2":
        return list(reconstruct_full(D2, CHAR)), D2
    return None, None


def _record(points, D, b, j, k, r, cj, ck, char=CHAR):
    """Rebuild ONE solution record in exactly the field order kk_closure/_try writes it, so a banked
    record can be compared as BYTES rather than field-by-field."""
    ok, X, dists = exact_check(points, D, b, j, k, r, cj, ck, char)
    if not ok:
        return None
    return {"r": r, "cj": cj, "ck": ck, "x": str(X[0]), "y_coeff": str(X[1]),
            "distances": dists, "is_vertex": any(X == p for p in points)}


# --------------------------------------------------------------------------------------------------
# the in-route ground-truth gate
# --------------------------------------------------------------------------------------------------
def ground_truth_gate() -> dict:
    checks: dict = {}

    def put(name, ok, **extra):
        checks[name] = {"pass": bool(ok), **extra}

    if not _SOLVER_OK:
        put("solver_importable", False, error=_SOLVER_ERR)
        return {"ok": False, "detail": checks}
    put("solver_importable", True)

    if not os.path.exists(BANKED_FULL):
        put("banked_artifact_present", False, path=BANKED_FULL)
        return {"ok": False, "detail": checks}
    with open(BANKED_FULL, encoding="utf-8") as fh:
        banked = json.load(fh)
    put("banked_artifact_present", True, path=os.path.basename(BANKED_FULL))

    # 1. the banked artifact must itself be a COMPLETE run
    put("banked_run_is_complete",
        banked.get("correctness_gate_pass") is True and banked.get("missing_vertices") == [],
        correctness_gate_pass=banked.get("correctness_gate_pass"))

    # 2. base triple re-derived
    cost, b, j, k, dbj, dbk = best_triple(D1, 7)
    put("base_triple_reproduced", [b, j, k] == banked.get("base"),
        got=[b, j, k], banked=banked.get("base"))
    put("grid_size_reproduced", cost == banked.get("grid_cells"),
        got=cost, banked=banked.get("grid_cells"))

    # 3. BYTE-IDENTICAL reproduction of every banked solution record
    sols = banked.get("solutions") or []
    byte_ok, mismatches = True, []
    for s in sols:
        rebuilt = _record(H1, D1, b, j, k, s["r"], s["cj"], s["ck"])
        if rebuilt is None or json.dumps(rebuilt, sort_keys=True) != json.dumps(s, sort_keys=True):
            byte_ok = False
            mismatches.append({"cell": [s.get("r"), s.get("cj"), s.get("ck")]})
    put("banked_solutions_reproduced_byte_identical", byte_ok and len(sols) == 6,
        n=len(sols), mismatches=mismatches)

    # 4. INDEPENDENT arithmetic path (pure bigint, no Fractions) recovers the same r
    M, Uj, Wj, Uk, Wk, G = integer_setup(H1, D1, b, j, k, CHAR)
    indep_ok = all(s["r"] in exact_cell(M, Uj, Wj, Uk, Wk, G, dbj, dbk, s["cj"], s["ck"], CHAR)
                   for s in sols)
    put("independent_integer_path_agrees", indep_ok)

    # 5. the solution set is exactly the 6 non-base vertices; nothing non-trivial
    want = {(str(H1[i][0]), str(H1[i][1])) for i in range(7) if i != b}
    got = {(s["x"], s["y_coeff"]) for s in sols}
    put("solution_set_is_exactly_the_non_base_vertices", want == got,
        missing=sorted(want - got), extra=sorted(got - want))
    put("banked_has_no_nontrivial_solution", banked.get("nontrivial") == [])

    # 6. INDEPENDENT re-verification by the kernel checker (imports nothing from the search code)
    if not _CHECKER_OK:
        put("independent_kernel_checker", False, error=_CHECKER_ERR)
    else:
        ok_pts, why = check_integral_point_set(7, CHAR, [[str(p[0]), str(p[1])] for p in H1])
        put("independent_kernel_checker", ok_pts, detail=why)

    # 7. NEGATIVE CONTROL — a perturbed cell must be REJECTED. Without this the whole gate would pass
    #    against an exact_check that answered True unconditionally.
    s0 = sols[0] if sols else None
    perturbed = _record(H1, D1, b, j, k, s0["r"] + 1, s0["cj"], s0["ck"]) if s0 else "no-solutions"
    put("perturbed_cell_is_rejected", perturbed is None or
        json.dumps(perturbed, sort_keys=True) != json.dumps(s0, sort_keys=True))

    # the sub-sweep artifact: shape only. It is NOT a closure and is never a source of one.
    if os.path.exists(BANKED_SUB):
        with open(BANKED_SUB, encoding="utf-8") as fh:
            sub = json.load(fh)
        put("sub_sweep_artifact_is_labelled_partial",
            "eps" in sub and sub.get("exact_cells", 0) < banked.get("grid_cells", 0),
            exact_cells=sub.get("exact_cells"), of=banked.get("grid_cells"))

    return {"ok": all(v["pass"] for v in checks.values()), "detail": checks}


# --------------------------------------------------------------------------------------------------
# the work
# --------------------------------------------------------------------------------------------------
def _sweep(points, D, label, max_cells, char=CHAR, base=None, max_chunks=None) -> dict:
    """Chunked EXACT-INTEGER sweep, composed from the solver's own primitives
    (`integer_setup` + `exact_cell` + `exact_check`). Chunks are c_j rows.

    Stops when the cell budget or the chunk cap is hit, and says so. A stopped sweep is INCOMPLETE.
    """
    n = len(points)
    if base is None:
        cost, b, j, k, dbj, dbk = best_triple(D, n)
    else:
        b, j, k = base
        dbj, dbk = D[b][j], D[b][k]
        cost = (2 * dbj + 1) * (2 * dbk + 1)
    M, Uj, Wj, Uk, Wk, G = integer_setup(points, D, b, j, k, char)
    cks = list(range(-dbk, dbk + 1))
    chunks_total = 2 * dbj + 1
    if max_chunks is not None:
        chunks_total = min(chunks_total, max_chunks) if max_chunks >= 0 else chunks_total

    solutions, cells, chunks_done = [], 0, 0
    t0 = time.perf_counter()
    stopped = None
    for cj in range(-dbj, dbj + 1):
        if max_chunks is not None and chunks_done >= max_chunks:
            stopped = "chunk_cap_reached"
            break
        if max_cells is not None and cells + len(cks) > max_cells:
            stopped = "cell_budget_exhausted"
            break
        for ck in cks:
            cells += 1
            for r in exact_cell(M, Uj, Wj, Uk, Wk, G, dbj, dbk, cj, ck, char):
                if r - cj < 0 or r - ck < 0:
                    continue
                rec = _record(points, D, b, j, k, r, cj, ck, char)
                if rec is not None:
                    solutions.append(rec)
        chunks_done += 1
    el = time.perf_counter() - t0

    seen, uniq = set(), []
    for s in solutions:
        key = (s["x"], s["y_coeff"])
        if key not in seen:
            seen.add(key)
            uniq.append(s)
    solutions = uniq
    nontrivial = [s for s in solutions if not s["is_vertex"]]

    full_chunks = 2 * dbj + 1
    complete = stopped is None and chunks_done == full_chunks and cells == cost

    want = {(str(points[i][0]), str(points[i][1])) for i in range(n) if i != b}
    got = {(s["x"], s["y_coeff"]) for s in solutions}
    missing = sorted(want - got)
    correctness_gate = complete and not missing

    return {"label": label, "base": [b, j, k], "grid_cells": cost,
            "cells_examined": cells, "cells_total": cost,
            "chunks_completed": chunks_done, "chunks_total": full_chunks,
            "coverage_fraction": round(cells / cost, 9) if cost else 0.0,
            "complete": complete, "stopped_because": stopped,
            "solutions": solutions, "nontrivial": nontrivial,
            "missing_vertices": missing,
            "correctness_gate_pass": bool(correctness_gate),
            "seconds": round(el, 3),
            "method": "exact integer arithmetic (exact_cell) + exact Fraction verification (exact_check)"}


def _verdict_for(sw: dict) -> dict:
    """The ONE place a verdict is assigned. INCOMPLETE is never upgraded."""
    if not sw["complete"]:
        return {"verdict": "INCOMPLETE_BUDGET_EXHAUSTED",
                "claim_strength": (
                    "NONE — this is a PARTIAL sweep "
                    f"({sw['chunks_completed']}/{sw['chunks_total']} chunks, "
                    f"{sw['cells_examined']:,}/{sw['cells_total']:,} cells). It licenses NO closure "
                    "claim whatsoever: the point that breaks the closure could lie in any unexamined "
                    "cell. Re-run with max_cells=full for the unconditional result.")}
    if not sw["correctness_gate_pass"]:
        return {"verdict": "UNSOUND_SWEEP",
                "claim_strength": ("NONE — the completed sweep failed to recover known non-base "
                                   f"vertices {sw['missing_vertices']}; the enumeration is unsound "
                                   "and no closure claim may be made from it.")}
    if sw["nontrivial"]:
        return {"verdict": "EXTENSION_FOUND",
                "claim_strength": ("PROOF (existence) — the exhaustive no-height-bound enumeration "
                                   "produced a point at integral distance from every vertex, verified "
                                   "in exact Fraction arithmetic.")}
    return {"verdict": "CLOSED_EXHAUSTIVE",
            "claim_strength": ("PROOF (non-existence, UNCONDITIONAL) — every one of the "
                               f"{sw['cells_total']:,} grid cells was solved in exact integer "
                               "arithmetic, with no height bound anywhere in the argument, and the "
                               "only points at integral distance from all vertices are the vertices "
                               "themselves.")}


def _parse_run(args: str) -> dict:
    toks = [t.strip() for t in args.split(",") if t.strip()]
    if len(toks) not in (2, 3):
        return {"error": "run_needs_<label>,<D>[,<max_cells>]"}
    label, dlab = toks[0].upper(), toks[1].upper()
    if label not in ("H1", "H2"):
        return {"error": f"unknown_label:{toks[0]!r} (expected H1|H2)"}
    if dlab not in ("D1", "D2"):
        return {"error": f"unknown_matrix:{toks[1]!r} (expected D1|D2)"}
    if {"H1": "D1", "H2": "D2"}[label] != dlab:
        return {"error": f"label_matrix_mismatch:{label}_requires_{ {'H1':'D1','H2':'D2'}[label] }"}
    max_cells: int | None = DEFAULT_MAX_CELLS
    if len(toks) == 3:
        t = toks[2].lower()
        if t == "full":
            max_cells = None
        elif re.fullmatch(r"\d+", t) and int(t) > 0:
            max_cells = int(t)
        else:
            return {"error": f"max_cells_must_be_a_positive_int_or_full:{toks[2]!r}"}
    return {"label": label, "dlab": dlab, "max_cells": max_cells}


def _do_check(args: str) -> dict:
    """`kkclosure:check,<points>,<D>` — INDEPENDENT verification of a candidate integral point set."""
    m = re.match(r"^(.*),\s*(D1|D2)\s*$", args.strip(), re.IGNORECASE)
    if not m:
        return {"ok": False, "trusted": False, "reason": "check_needs_<points>,<D1|D2>"}
    pts_tok, dlab = m.group(1).strip(), m.group(2).upper()
    if pts_tok.upper() in ("H1", "H2"):
        points, _ = _load_set(pts_tok)
        raw = [[str(p[0]), str(p[1])] for p in points]
    else:
        try:
            raw = json.loads(pts_tok)
            points = [(Fraction(str(p[0])), Fraction(str(p[1]))) for p in raw]
            raw = [[str(a), str(b)] for a, b in points]
        except Exception as exc:
            return {"ok": False, "trusted": False, "reason": f"unparsable_points: {exc}"}
    if not _CHECKER_OK:
        return GA.tool_missing("oracle.kernel.checkers_unavailable", _CHECKER_ERR or "")
    ok_pts, why = check_integral_point_set(len(raw), CHAR, raw)
    D = D1 if dlab == "D1" else D2
    matches_D = None
    if len(points) == len(D):
        matches_D = all(
            (points[a][0] - points[b][0]) ** 2 + CHAR * (points[a][1] - points[b][1]) ** 2
            == Fraction(D[a][b]) ** 2
            for a in range(len(D)) for b in range(len(D)))
    return {"ok": True, "trusted": True, "mode": "check", "n": len(raw), "char": CHAR,
            "verdict": "INTEGRAL_POINT_SET_VERIFIED" if ok_pts else "NOT_AN_INTEGRAL_POINT_SET",
            "claim_strength": ("PROOF — verified from the problem statement in exact Fraction/isqrt "
                               "arithmetic by oracle/kernel/checkers.py, which imports nothing from "
                               "the search code" if ok_pts else
                               "REJECTED — the candidate fails the integral-general-position definition"),
            "detail": why, "matches_distance_matrix": matches_D,
            "chunks_completed": None, "chunks_total": None,
            "note": "check mode verifies a point SET; it makes no closure claim (see run mode)"}


def attack(query: str, time_limit_s: float = 45.0) -> dict:
    m = GRAMMAR.match(query or "")
    if not m:
        return GA.off_grammar("kkclosure")
    if not _SOLVER_OK:
        return GA.tool_missing("kk_closure_solver_unavailable", _SOLVER_ERR or "")

    mode, args = m.group(1).lower(), m.group(2)

    if mode == "check":
        return GA.run_gated(
            name="kk_closure_bridge",
            gate_fn=ground_truth_gate,
            work_fn=lambda: _do_check(args),
            gate_description=("banked closure_full_exact_H1.json reproduced byte-identically "
                              "(+ independent bigint path, kernel checker, perturbed-cell control)"),
            gate_budget_s=GATE_BUDGET_S,
            extra={"mode": "check", "shape": "SHAPE-134"},
        )

    parsed = _parse_run(args)
    if "error" in parsed:
        return {"ok": False, "trusted": False, "reason": parsed["error"]}

    def work():
        points, D = _load_set(parsed["label"])
        sw = _sweep(points, D, parsed["label"], parsed["max_cells"])
        sw.update(_verdict_for(sw))
        sw["budget"] = {"max_cells": parsed["max_cells"] if parsed["max_cells"] is not None else "full",
                        "default_is_a_budget_not_full": parsed["max_cells"] == DEFAULT_MAX_CELLS}
        sw["full_sweep_measured_seconds"] = FULL_SWEEP_MEASURED_S
        # a banked COMPLETE result is reference material, never this run's verdict
        if parsed["label"] == "H1" and os.path.exists(BANKED_FULL):
            with open(BANKED_FULL, encoding="utf-8") as fh:
                bk = json.load(fh)
            sw["banked_reference"] = {
                "artifact": os.path.basename(BANKED_FULL),
                "grid_cells": bk.get("grid_cells"), "seconds": bk.get("seconds"),
                "nontrivial": len(bk.get("nontrivial") or []),
                "note": ("a separately banked COMPLETE sweep, reproduced byte-identically by this "
                         "call's gate; it is NOT the verdict of this call"),
            }
        return sw

    out = GA.run_gated(
        name="kk_closure_bridge",
        gate_fn=ground_truth_gate,
        work_fn=work,
        gate_description=("banked closure_full_exact_H1.json reproduced byte-identically "
                          "(+ independent bigint path, kernel checker, perturbed-cell control)"),
        gate_budget_s=GATE_BUDGET_S,
        extra={"mode": "run", "shape": "SHAPE-134"},
    )
    if out.get("ok"):
        res = out["result"]
        out["verdict"] = res["verdict"]
        out["claim_strength"] = res["claim_strength"]
        out["chunks_completed"] = res["chunks_completed"]
        out["chunks_total"] = res["chunks_total"]
    return out


# --------------------------------------------------------------------------------------------------
# selftest fixture: the 3-4-5 triangle over Q(sqrt(1)) — a genuine integral point set whose ENTIRE
# grid is small enough to sweep exhaustively inside a selftest. H1's smallest grid is 136,801,313
# cells (measured below), so the complete-run branch cannot be exercised on H1 itself.
# --------------------------------------------------------------------------------------------------
_TOY_PTS = [(Fraction(0), Fraction(0)), (Fraction(5), Fraction(0)),
            (Fraction(16, 5), Fraction(12, 5))]
_TOY_D = [[0, 5, 4], [5, 0, 3], [4, 3, 0]]


# --------------------------------------------------------------------------------------------------
# selftest
# --------------------------------------------------------------------------------------------------
def selftest() -> dict:
    checks = []
    t = {}

    checks.append(("recognises_run_grammar", is_kk_closure_target("kkclosure:run,H1,D1"), True))
    checks.append(("recognises_check_grammar",
                   is_kk_closure_target('kkclosure:check,H1,D1'), True))
    checks.append(("rejects_natural_language",
                   is_kk_closure_target("is the KK heptagon closed?"), False))
    checks.append(("rejects_other_grammar", is_kk_closure_target("intrel:constant,pi"), False))
    checks.append(("refuses_off_grammar", attack("nope").get("reason"),
                   "query_not_in_kkclosure_grammar"))

    # ---- the gate, both directions -----------------------------------------------------
    t0 = time.perf_counter()
    g = ground_truth_gate()
    t["gate_s"] = round(time.perf_counter() - t0, 3)
    checks.append(("gate_is_green", g["ok"], True))
    checks.append(("gate_within_budget", t["gate_s"] <= GATE_BUDGET_S, True))
    checks.append(("gate_reproduces_banked_bytes",
                   g["detail"]["banked_solutions_reproduced_byte_identical"]["pass"], True))
    checks.append(("gate_has_negative_control",
                   g["detail"]["perturbed_cell_is_rejected"]["pass"], True))
    checks.append(("gate_uses_independent_kernel_checker",
                   g["detail"]["independent_kernel_checker"]["pass"], True))

    # a red gate must WITHHOLD and must not run the work
    ran = []
    red = GA.run_gated(name="kk_closure_bridge", gate_fn=lambda: {"ok": False, "detail": {"x": {"pass": False}}},
                       work_fn=lambda: ran.append(1) or {}, gate_description="forced-red")
    checks.append(("red_gate_withholds_verdict", red["trusted"], False))
    checks.append(("red_gate_does_not_run_sweep", ran, []))

    # ---- THE TRAP: a truncated run is INCOMPLETE, never CLOSED --------------------------
    t1 = time.perf_counter()
    trunc = _sweep(H1, D1, "H1", max_cells=None, max_chunks=3)
    trunc.update(_verdict_for(trunc))
    t["truncated_sweep_s"] = round(time.perf_counter() - t1, 3)
    checks.append(("truncated_run_is_INCOMPLETE_not_CLOSED",
                   trunc["verdict"], "INCOMPLETE_BUDGET_EXHAUSTED"))
    checks.append(("truncated_run_is_not_complete", trunc["complete"], False))
    checks.append(("truncated_run_reports_chunk_counts",
                   trunc["chunks_completed"] == 3 and trunc["chunks_total"] > 3, True))
    checks.append(("truncated_run_claim_strength_is_NONE",
                   trunc["claim_strength"].startswith("NONE"), True))
    checks.append(("truncated_run_correctness_gate_is_false",
                   trunc["correctness_gate_pass"], False))

    # budgeted route call: same property end-to-end through attack()
    t2 = time.perf_counter()
    a = attack("kkclosure:run,H1,D1,20000")
    t["budgeted_route_s"] = round(time.perf_counter() - t2, 3)
    checks.append(("budgeted_route_is_trusted", a.get("trusted"), True))
    checks.append(("budgeted_route_is_INCOMPLETE", a.get("verdict"), "INCOMPLETE_BUDGET_EXHAUSTED"))
    checks.append(("budgeted_route_reports_chunks",
                   a.get("chunks_completed") is not None and a.get("chunks_total") is not None, True))
    checks.append(("budgeted_route_banked_is_not_the_verdict",
                   a["result"]["banked_reference"]["note"].endswith("NOT the verdict of this call"), True))
    checks.append(("default_is_a_budget_not_full",
                   attack("kkclosure:run,H1,D1")["result"]["complete"], False))

    # ---- the COMPLETE direction, on a sweep that genuinely finishes ---------------------
    # MEASURED FACT, not a convenience: the SMALLEST base triple on H1 is 136,801,313 cells (the
    # banked one, [4,5,6]) — there is no cheap base, and a complete H1 sweep is 153-170s. So the
    # complete-run branch is exercised on a TOY integral point set (the 3-4-5 triangle over
    # Q(sqrt(1))), whose whole grid is swept here, while H1's completeness rests on the banked
    # artifact reproduced byte-identically by the gate above. Both facts are reported, neither
    # is dressed up as the other.
    t3 = time.perf_counter()
    full = _sweep(_TOY_PTS, _TOY_D, "toy345", max_cells=None, char=1)
    full.update(_verdict_for(full))
    t["complete_toy_sweep_s"] = round(time.perf_counter() - t3, 3)
    checks.append(("complete_sweep_reports_complete", full["complete"], True))
    checks.append(("complete_sweep_examines_every_cell",
                   full["cells_examined"] == full["cells_total"], True))
    checks.append(("complete_sweep_all_chunks",
                   full["chunks_completed"] == full["chunks_total"], True))
    checks.append(("complete_sweep_recovers_all_non_base_vertices", full["missing_vertices"], []))
    checks.append(("complete_sweep_correctness_gate_passes", full["correctness_gate_pass"], True))
    checks.append(("complete_sweep_gets_a_terminal_verdict",
                   full["verdict"] in ("CLOSED_EXHAUSTIVE", "EXTENSION_FOUND"), True))
    checks.append(("complete_sweep_is_never_INCOMPLETE",
                   full["verdict"] == "INCOMPLETE_BUDGET_EXHAUSTED", False))
    # same grid, one chunk short -> must flip to INCOMPLETE. This is the trap, isolated.
    near = _sweep(_TOY_PTS, _TOY_D, "toy345", max_cells=None, char=1,
                  max_chunks=full["chunks_total"] - 1)
    near.update(_verdict_for(near))
    checks.append(("one_chunk_short_is_INCOMPLETE", near["verdict"], "INCOMPLETE_BUDGET_EXHAUSTED"))
    checks.append(("h1_has_no_cheap_base_triple",
                   min((2 * D1[b][j] + 1) * (2 * D1[b][k] + 1)
                       for b in range(7) for j in range(7) for k in range(7)
                       if len({b, j, k}) == 3) == 136_801_313, True))

    # ---- check mode --------------------------------------------------------------------
    t4 = time.perf_counter()
    c = attack("kkclosure:check,H1,D1")
    t["check_route_s"] = round(time.perf_counter() - t4, 3)
    checks.append(("check_route_is_trusted", c.get("trusted"), True))
    checks.append(("check_verifies_H1", c["result"]["verdict"], "INTEGRAL_POINT_SET_VERIFIED"))
    checks.append(("check_confirms_distance_matrix", c["result"]["matches_distance_matrix"], True))
    bad = attack('kkclosure:check,[["0","0"],["1","0"],["2","0"]],D1')
    checks.append(("check_rejects_collinear_junk",
                   bad["result"]["verdict"], "NOT_AN_INTEGRAL_POINT_SET"))

    # ---- grammar guards ----------------------------------------------------------------
    checks.append(("rejects_label_matrix_mismatch",
                   attack("kkclosure:run,H1,D2").get("reason") is not None, True))
    checks.append(("rejects_bad_budget",
                   attack("kkclosure:run,H1,D1,-5").get("reason") is not None, True))

    detail = {name: {"pass": got == want, "got": got} for name, got, want in checks}
    return {"ok": all(v["pass"] for v in detail.values()),
            "method": "KK integral-distance closure bridge — SHAPE-134",
            "timing_s": t,
            "full_sweep_measured_seconds": FULL_SWEEP_MEASURED_S,
            "gate_budget_s": GATE_BUDGET_S,
            "detail": detail}


if __name__ == "__main__":
    out = selftest()
    print(json.dumps(out, indent=2, ensure_ascii=False, default=str))
    sys.exit(0 if out["ok"] else 1)
