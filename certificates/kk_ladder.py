"""kk_ladder.py — the SUB-CONFIGURATION EXCLUSION LADDER (closure-wedge attack on n=8).

WEDGE SHAPE. The diameter bound d(2,8) > 30000 inherits Kreisel-Kurz's 30000 exhaustiveness ceiling, so it
is a residual bounded only from BELOW -- not closure. Sub-configuration exclusion has no such ceiling:

    "no integral octagon in general position contains <configuration S>"

is proved by exhaustively enumerating E(S) = { X : |X P| integral for all P in S } -- which the
triangle-inequality parametrisation makes FINITE regardless of how far X lies -- and showing E(S) contains no
(8-|S|)-clique that is mutually at integral distance and in general position with S.

THE LADDER. An octagon contains 28 hexagons, 56 pentagons, 70 quadrilaterals, 56 triangles. Excluding a
whole level is strictly stronger than excluding the level above it:

    |S|=6 -> need a 2-clique in E(S)      (DONE for H1: 0 octagons)
    |S|=5 -> need a 3-clique
    |S|=4 -> need a 4-clique
    |S|=3 -> need a 5-clique

BOUNDARY PROBE FIRST (doctrine 11.2): the campaign is only viable if |E(S)| stays small as |S| shrinks.
This module measures that growth before committing compute to a full level.

Run: python oracle/kbk/engine/kk_ladder.py probe        # measure |E(S)| growth
     python oracle/kbk/engine/kk_ladder.py level <m>    # exclude all |S|=m subsets of H1
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import combinations
import json
import sys
import time

sys.path.insert(0, "oracle/kbk/engine")
from kk_heptagon import H1, D1
from kk_closure_full_exact import full_exact_sweep
from kk_hexagon_octagon import dist_int, general_position

CHAR = 2002


def extensions_of(subset_idx, points=H1, D=D1, char=CHAR, quiet=True):
    """E(S) beyond S itself, for S = the given index subset. Exhaustive + exact."""
    S = [points[i] for i in subset_idx]
    DS = [[D[i][j] for j in subset_idx] for i in subset_idx]
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf) if quiet else contextlib.nullcontext():
        res = full_exact_sweep(S, DS, f"S{''.join(map(str,subset_idx))}", char=char, report=1)
    if not res["correctness_gate_pass"]:
        return None, res
    ext = [(F(s["x"]), F(s["y_coeff"])) for s in res["solutions"] if not s["is_vertex"]]
    return ext, res


def exclude_subset(subset_idx, need, points=H1, D=D1, char=CHAR):
    """Is there a `need`-clique in E(S), mutually integral + general position with S?"""
    S = [points[i] for i in subset_idx]
    ext_all, res = extensions_of(subset_idx, points, D, char)
    if ext_all is None:
        return {"subset": list(subset_idx), "gate": False}

    # ---- PRUNING THEOREM (the chord-forcing invariant, used as a NECESSARY condition) ----
    # Every point of an octagon must INDIVIDUALLY be in general position with S. So any X in E(S)
    # for which S u {X} already has 3 collinear or 4 concyclic points cannot occur in ANY octagon
    # containing S, and is discarded before the clique search. This is sound, not heuristic -- and
    # the 66x chord enrichment predicts it removes most of E(S).
    ext = [X for X in ext_all if general_position(S + [X], char)[0]]
    n_pruned = len(ext_all) - len(ext)

    # compatibility graph on the SURVIVORS
    n = len(ext)
    adj = {i: set() for i in range(n)}
    for a, b in combinations(range(n), 2):
        if dist_int(ext[a], ext[b], char) is not None:
            adj[a].add(b); adj[b].add(a)
    n_edges = sum(len(v) for v in adj.values()) // 2
    # graph-walking clique enumeration -- never touches C(n, need)
    found = []
    cliques = 0

    def extend(cur, cand_set):
        nonlocal cliques
        if len(cur) == need:
            cliques += 1
            pts = S + [ext[i] for i in cur]
            ok, _ = general_position(pts, char)
            if ok:
                found.append([[str(ext[i][0]), str(ext[i][1])] for i in cur])
            return
        for v in sorted(cand_set):
            if v <= (cur[-1] if cur else -1):
                continue
            extend(cur + [v], cand_set & adj[v])

    extend([], set(range(n)))
    return {"subset": list(subset_idx), "gate": True,
            "n_ext_raw": len(ext_all), "n_pruned_by_general_position": n_pruned,
            "n_ext": n, "n_edges": n_edges,
            "cliques_needed": need, "cliques_found": cliques, "octagons": found,
            "seconds": res["seconds"], "cells": res["grid_cells"]}


def probe():
    """Boundary probe: how does |E(S)| grow as |S| shrinks? Keep P4,P5,P6 so the grid stays small."""
    print("BOUNDARY PROBE — |E(S)| vs |S| (subsets of H1 retaining P4,P5,P6 for a small grid)\n")
    tests = [
        (0, 1, 2, 3, 4, 5, 6),        # the whole heptagon
        (1, 2, 3, 4, 5, 6),           # hexagon
        (2, 3, 4, 5, 6),              # pentagon
        (3, 4, 5, 6),                 # quadrilateral
        (4, 5, 6),                    # triangle
    ]
    rows = []
    for sub in tests:
        t0 = time.time()
        ext, res = extensions_of(sub)
        if ext is None:
            print(f"  |S|={len(sub)} {sub}: GATE FAIL"); continue
        need = 8 - len(sub)
        print(f"  |S|={len(sub)} {sub}: |E(S)| = {len(ext):5d}   need a {need}-clique   "
              f"({res['grid_cells']:,} cells, {time.time()-t0:.0f}s)", flush=True)
        rows.append({"size": len(sub), "subset": list(sub), "n_ext": len(ext),
                     "need": need, "cells": res["grid_cells"], "seconds": time.time() - t0})
        with open("oracle/kbk/engine/ladder_probe.json", "w") as f:
            json.dump(rows, f, indent=1)
    print("\nVIABILITY: the ladder is climbable while |E(S)| stays small enough to enumerate cliques.")
    return rows


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "probe"
    if mode == "probe":
        probe()
    else:
        m = int(sys.argv[2])
        need = 8 - m
        out = []
        for sub in combinations(range(7), m):
            r = exclude_subset(sub, need)
            out.append(r)
            tag = f"{len(r.get('octagons',[]))} OCTAGONS" if r.get("gate") else "GATE FAIL"
            print(f"  S={sub}: |E|={r.get('n_ext','?')} need {need}-clique -> {tag}", flush=True)
            with open(f"oracle/kbk/engine/ladder_level{m}.json", "w") as f:
                json.dump(out, f, indent=1)
        tot = sum(len(r.get("octagons", [])) for r in out)
        print(f"\nLEVEL |S|={m}: TOTAL OCTAGONS FOUND = {tot}")
