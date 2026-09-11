"""kk_hexagon_octagon.py — can a HEXAGON of H1/H2 be extended by TWO points to an octagon?

WHY THIS IS A GENUINELY NEW TEST. The closure sweep proved "H1 + X" has no solution: no eighth point can be
added to the full heptagon. It says NOTHING about dropping one vertex of H1 and adding TWO different points.
That 8-set contains a hexagon of H1 but not H1 itself, so it is completely untouched by the earlier result.

And crucially this test has NO DIAMETER CEILING: the triangle-inequality parametrisation enumerates ALL
points at integral distance from the hexagon, however far away they lie. So a negative answer here is a
statement about octagons of ANY diameter -- unlike the d(2,8) > 30000 bound, which inherits the 30000 ceiling
from Kreisel-Kurz's exhaustive range.

METHOD. For each 6-subset Hex of a certified heptagon:
  1. exhaustively enumerate E(Hex) = { X : |X P| integral for all P in Hex }   (exact integer arithmetic)
  2. for every pair X != Y in E(Hex) \\ Hex with |XY| integral, test whether Hex + {X,Y} is a genuine
     integral octagon in general position (no 3 collinear, no 4 concyclic).
A single surviving pair IS AN OCTAGON and settles Erdos's n=8 question affirmatively.

Run: python oracle/kbk/engine/kk_hexagon_octagon.py [H1|H2] [omit_index|all]
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import combinations
from math import isqrt
import json
import sys
import time

sys.path.insert(0, "oracle/kbk/engine")
from kk_heptagon import H1, D1, D2, reconstruct_full
from kk_closure_full_exact import full_exact_sweep

CHAR = 2002


def dist_int(a, b, char=CHAR):
    """Exact integer distance, or None."""
    q = (a[0] - b[0]) ** 2 + char * (a[1] - b[1]) ** 2
    if q < 0 or q.denominator != 1:
        return None
    r = isqrt(q.numerator)
    return r if r * r == q.numerator else None


def collinear(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) == 0


def concyclic(a, b, c, d, char=CHAR):
    def row(p):
        return [p[0] * p[0] + char * p[1] * p[1], p[0], p[1], F(1)]
    m = [row(a), row(b), row(c), row(d)]
    def d3(mm):
        return (mm[0][0] * (mm[1][1] * mm[2][2] - mm[1][2] * mm[2][1])
                - mm[0][1] * (mm[1][0] * mm[2][2] - mm[1][2] * mm[2][0])
                + mm[0][2] * (mm[1][0] * mm[2][1] - mm[1][1] * mm[2][0]))
    return sum(((-1) ** j) * m[0][j] * d3([[m[i][k] for k in range(4) if k != j] for i in range(1, 4)])
               for j in range(4)) == 0


def general_position(pts, char=CHAR):
    for a, b, c in combinations(pts, 3):
        if collinear(a, b, c):
            return False, "collinear"
    for a, b, c, d in combinations(pts, 4):
        if concyclic(a, b, c, d, char):
            return False, "concyclic"
    return True, "ok"


def analyse(points, D, label, omit, char=CHAR):
    hex_pts = [points[i] for i in range(len(points)) if i != omit]
    hex_D = [[D[i][j] for j in range(len(points)) if j != omit]
             for i in range(len(points)) if i != omit]
    print(f"\n===== {label}: hexagon = heptagon minus P{omit} =====", flush=True)
    res = full_exact_sweep(hex_pts, hex_D, f"{label}-drop{omit}", char=char, report=10)
    if not res["correctness_gate_pass"]:
        print("  !! GATE FAILED -- result discarded"); return None

    ext = [(F(s["x"]), F(s["y_coeff"])) for s in res["solutions"] if not s["is_vertex"]]
    print(f"  E(Hex) non-hexagon extension points: {len(ext)}")
    for e in ext:
        known = "= dropped vertex P%d" % omit if e == points[omit] else "NEW"
        print(f"    ({e[0]}, {e[1]}*sqrt{char})   {known}")

    octagons = []
    for X, Y in combinations(ext, 2):
        dxy = dist_int(X, Y, char)
        if dxy is None:
            continue
        cand = hex_pts + [X, Y]
        ok, why = general_position(cand, char)
        print(f"    PAIR at integral distance {dxy}: general position -> {ok} ({why})", flush=True)
        if ok:
            octagons.append({"X": [str(X[0]), str(X[1])], "Y": [str(Y[0]), str(Y[1])], "dXY": dxy})
            print("    *** OCTAGON FOUND *** ", flush=True)
    print(f"  integral-distance extension PAIRS: "
          f"{sum(1 for X,Y in combinations(ext,2) if dist_int(X,Y,char))}; OCTAGONS: {len(octagons)}")
    return {"label": label, "omit": omit, "n_extensions": len(ext),
            "extensions": [[str(e[0]), str(e[1])] for e in ext],
            "octagons": octagons, "seconds": res["seconds"]}


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "H1"
    sel = sys.argv[2] if len(sys.argv) > 2 else "all"
    pts, D = (H1, D1) if which == "H1" else (reconstruct_full(D2, 2002), D2)
    omits = range(len(pts)) if sel == "all" else [int(sel)]
    out = []
    t0 = time.time()
    for o in omits:
        r = analyse(pts, D, which, o)
        if r:
            out.append(r)
        with open(f"oracle/kbk/engine/hexagon_octagon_{which}.json", "w") as f:
            json.dump(out, f, indent=1)
    tot = sum(len(r["octagons"]) for r in out)
    print(f"\n===== {which} SUMMARY ({time.time()-t0:.0f}s) =====")
    for r in out:
        print(f"  drop P{r['omit']}: {r['n_extensions']} extension points, {len(r['octagons'])} octagons")
    print(f"  TOTAL OCTAGONS FOUND: {tot}")
