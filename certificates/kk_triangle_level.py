"""kk_triangle_level.py — THE BOTTOM RUNG: exclude every triangle of H1 (and H2) from any octagon.

WHY THIS SUBSUMES EVERYTHING. Every sub-configuration of size >= 3 contains a triangle. So

    "no integral octagon in general position contains any triangle of H"

implies at once: no octagon contains any quadrilateral / pentagon / hexagon of H, and H itself is maximal.
It is the strongest statement of the exclusion ladder and it carries NO DIAMETER CEILING, because the
triangle-inequality parametrisation enumerates E(T) exhaustively however far the extension points lie.

METHOD per triangle T:
  1. exhaustively enumerate E(T) exactly (integer arithmetic, correctness-gated);
  2. PRUNE by the general-position necessary condition (measured to remove ~89%);
  3. build the integral-distance compatibility graph on the survivors;
  4. enumerate 5-cliques by walking the graph;
  5. an octagon exists iff some 5-clique is in general position together with T.

Runs cost-ordered (cheapest grids first) with incremental receipts, so partial results are usable.

Run: python oracle/kbk/engine/kk_triangle_level.py [H1|H2] [max_cells]
"""
from __future__ import annotations
from itertools import combinations
import json
import sys
import time

sys.path.insert(0, "oracle/kbk/engine")
from kk_heptagon import H1, D1, D2, reconstruct_full
from kk_ladder import exclude_subset
from kk_closure_full_exact import best_triple

CHAR = 2002


def grid_cost(tri, D):
    """Cells the exact sweep will examine for this triangle."""
    sub = [[D[i][j] for j in tri] for i in tri]
    cost, b, j, k, dbj, dbk = best_triple(sub, 3)
    return cost


def run(points, D, label, max_cells=None):
    tris = list(combinations(range(len(points)), 3))
    ordered = sorted(tris, key=lambda t: grid_cost(t, D))
    print(f"[{label}] triangle level: {len(ordered)} triangles, cost-ordered")
    print(f"[{label}] grid sizes: min {grid_cost(ordered[0],D):,}  max {grid_cost(ordered[-1],D):,}")
    if max_cells:
        ordered = [t for t in ordered if grid_cost(t, D) <= max_cells]
        print(f"[{label}] restricted to grids <= {max_cells:,} cells -> {len(ordered)} triangles")
    out, t0 = [], time.time()
    for n, tri in enumerate(ordered, 1):
        c = grid_cost(tri, D)
        print(f"\n[{label}] ({n}/{len(ordered)}) triangle {tri}  grid {c:,} cells ...", flush=True)
        r = exclude_subset(tri, 5, points, D, CHAR)
        r["grid_cost"] = c
        out.append(r)
        if not r.get("gate"):
            print(f"    !! GATE FAIL -> discarded", flush=True)
        else:
            print(f"    E={r['n_ext_raw']} -> pruned {r['n_pruned_by_general_position']} "
                  f"({100*r['n_pruned_by_general_position']/max(r['n_ext_raw'],1):.1f}%) "
                  f"-> {r['n_ext']} survivors, {r['n_edges']} edges, "
                  f"{r['cliques_found']} 5-cliques, OCTAGONS={len(r['octagons'])}", flush=True)
            if r["octagons"]:
                print(f"    *** OCTAGON FOUND *** {json.dumps(r['octagons'][0])}", flush=True)
        with open(f"oracle/kbk/engine/triangle_level_{label}.json", "w") as f:
            json.dump(out, f, indent=1)
        el = time.time() - t0
        print(f"    [{el:.0f}s elapsed, {el/n*(len(ordered)-n):.0f}s eta]", flush=True)
    ok = [r for r in out if r.get("gate")]
    tot = sum(len(r["octagons"]) for r in ok)
    print(f"\n===== [{label}] TRIANGLE LEVEL SUMMARY =====")
    print(f"  triangles completed : {len(ok)}/{len(ordered)}")
    print(f"  total 5-cliques examined : {sum(r['cliques_found'] for r in ok)}")
    print(f"  TOTAL OCTAGONS FOUND : {tot}")
    if tot == 0 and len(ok) == len(tris):
        print(f"  => THEOREM: no integral octagon in general position contains ANY triangle of {label}.")
        print(f"     (subsumes: no octagon contains any sub-configuration of {label}; {label} is maximal)")
    return out


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "H1"
    max_cells = int(sys.argv[2]) if len(sys.argv) > 2 else None
    pts, D = (H1, D1) if which == "H1" else (reconstruct_full(D2, 2002), D2)
    run(pts, D, which, max_cells)
