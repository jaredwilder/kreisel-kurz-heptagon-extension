"""kk_octagon_hunt.py — FIND AN OCTAGON. Constructive search, not exclusion.

Every computation in this lane so far has been an EXCLUSION run on H1/H2 -- configurations that are already
the most-searched objects in the literature. If an octagon exists it is somewhere nobody has looked, and an
exclusion-only programme is structurally guaranteed to miss it. This module hunts for the object.

THE HUNT. Both known heptagons have characteristic 2002 = 2*7*11*13 -- that is the arithmetic that supports
seven points in general position. So: generate integral triangles of characteristic 2002 that are NOT
sub-triangles of H1/H2, and drive each one toward eight points.

THE KEY EFFICIENCY. For any S containing T, E(S) is a SUBSET of E(T). So ONE exhaustive sweep per triangle
suffices: an octagon containing T is exactly a 5-clique inside E(T) that is mutually at integral distance
and in general position together with T. That is one sweep + one graph search per triangle.

A HIT HERE SETTLES ERDOS n=8 AFFIRMATIVELY.

Run: python oracle/kbk/engine/kk_octagon_hunt.py <max_side> [max_cells]
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import combinations
from math import isqrt
import json
import sys
import time

sys.path.insert(0, "oracle/kbk/engine")
from kk_heptagon import H1, _squarefree
from kk_closure_full_exact import full_exact_sweep
from kk_hexagon_octagon import dist_int, general_position

CHAR = 2002


def heron16(a, b, c):
    """16 * Area^2 for a triangle with integer sides."""
    return 2*a*a*b*b + 2*b*b*c*c + 2*c*c*a*a - a**4 - b**4 - c**4


def rational_sqrt(q: F):
    if q < 0:
        return None
    rn, rd = isqrt(q.numerator), isqrt(q.denominator)
    return F(rn, rd) if rn*rn == q.numerator and rd*rd == q.denominator else None


def triangle_points(a, b, c, char=CHAR):
    """Points for a triangle with |P0P1|=a, |P0P2|=b, |P1P2|=c, in (x, y_coeff) form
    where the real point is (x, y_coeff*sqrt(char)). None if it does not live over Q(sqrt(char))."""
    x = F(a*a + b*b - c*c, 2*a)
    y2 = F(b)*b - x*x                      # real y^2
    if y2 <= 0:
        return None                        # degenerate / collinear
    q2 = y2 / char
    q = rational_sqrt(q2)
    if q is None:
        return None
    return [(F(0), F(0)), (F(a), F(0)), (x, q)]


def gen_triangles(max_side, char=CHAR):
    """Integral triangles with characteristic `char` and sides <= max_side (a <= b <= c).

    FAST TEST: squarefree(h) == char  <=>  h % char == 0 AND h/char is a perfect square.
    No factorisation needed. Vectorised over c with numpy (int64 holds h ~ max_side^4 up to ~55000)."""
    import numpy as np
    out = []
    for a in range(1, max_side + 1):
        a2 = a * a
        for b in range(a, max_side + 1):
            lo, hi = b, min(max_side, a + b - 1)
            if lo > hi:
                continue
            c = np.arange(lo, hi + 1, dtype=np.int64)
            u = c * c
            # h = -u^2 + 2u(a^2+b^2) - (a^2-b^2)^2
            h = -u * u + 2 * u * (a2 + b * b) - (a2 - b * b) ** 2
            m = (h > 0) & (h % char == 0)
            if not m.any():
                continue
            k = h[m] // char
            r = np.sqrt(k.astype(np.float64)).astype(np.int64)
            for cand in (-1, 0, 1):                      # guard float sqrt rounding
                ok = (r + cand) ** 2 == k
                for cc in c[m][ok]:
                    out.append((a, b, int(cc)))
    return sorted(set(out))


def hunt_triangle(a, b, c, char=CHAR):
    """One exhaustive sweep + 5-clique search. Returns the octagons found (normally none)."""
    pts = triangle_points(a, b, c, char)
    if pts is None:
        return None
    D = [[0]*3 for _ in range(3)]
    D[0][1] = D[1][0] = a
    D[0][2] = D[2][0] = b
    D[1][2] = D[2][1] = c
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        res = full_exact_sweep(pts, D, f"T{a}_{b}_{c}", char=char, report=1)
    if not res["correctness_gate_pass"]:
        return {"tri": (a, b, c), "gate": False}
    ext_all = [(F(s["x"]), F(s["y_coeff"])) for s in res["solutions"] if not s["is_vertex"]]
    # sound pruning: every octagon point is individually in general position with T
    ext = [X for X in ext_all if general_position(pts + [X], char)[0]]
    n = len(ext)
    adj = {i: set() for i in range(n)}
    for i, j in combinations(range(n), 2):
        if dist_int(ext[i], ext[j], char) is not None:
            adj[i].add(j); adj[j].add(i)
    found, cliques = [], 0

    def extend(cur, cand):
        nonlocal cliques
        if len(cur) == 5:
            cliques += 1
            cand_pts = pts + [ext[i] for i in cur]
            if general_position(cand_pts, char)[0]:
                found.append([[str(ext[i][0]), str(ext[i][1])] for i in cur])
            return
        for v in sorted(cand):
            if v <= (cur[-1] if cur else -1):
                continue
            extend(cur + [v], cand & adj[v])

    extend([], set(range(n)))
    return {"tri": (a, b, c), "gate": True, "cells": res["grid_cells"],
            "n_ext_raw": len(ext_all), "n_ext": n, "cliques": cliques,
            "octagons": found, "seconds": res["seconds"]}


def gen_triangles_thin(a_max, b_max, char=CHAR):
    """Char-`char` triangles with SMALLEST side <= a_max and the others <= b_max.

    THE LEVER: a sweep costs (2a+1)(2b+1), set by the two smallest sides. A thin triangle (one small
    side, two large) is therefore CHEAP while reaching a LARGE-diameter configuration -- the regime a
    square a,b <= B search never touches. a=100 reaches b ~ 50,000, far beyond H1's scale."""
    import numpy as np
    out = []
    for a in range(1, a_max + 1):
        a2 = a * a
        for b in range(a, b_max + 1):
            lo, hi = b, min(b_max, a + b - 1)
            if lo > hi:
                continue
            c = np.arange(lo, hi + 1, dtype=np.int64)
            u = c * c
            h = -u * u + 2 * u * (a2 + b * b) - (a2 - b * b) ** 2
            m = (h > 0) & (h % char == 0)
            if not m.any():
                continue
            k = h[m] // char
            r = np.sqrt(k.astype(np.float64)).astype(np.int64)
            for cand in (-1, 0, 1):
                for cc in c[m][(r + cand) ** 2 == k]:
                    out.append((a, b, int(cc)))
    return sorted(set(out))


if __name__ == "__main__":
    if sys.argv[1] == "thin":
        a_max = int(sys.argv[2]); b_max = int(sys.argv[3])
        max_cells = int(sys.argv[4]) if len(sys.argv) > 4 else 20_000_000
        print(f"generating THIN characteristic-{CHAR} triangles: min side <= {a_max}, "
              f"others <= {b_max} ...", flush=True)
        tris = gen_triangles_thin(a_max, b_max)
        print(f"  {len(tris)} found", flush=True)
        tris = [t for t in tris if (2*t[0]+1)*(2*t[1]+1) <= max_cells]
        tris.sort(key=lambda t: (2*t[0]+1)*(2*t[1]+1))
        print(f"  {len(tris)} within the {max_cells:,}-cell budget", flush=True)
        if tris:
            print(f"  largest diameter reachable: {max(t[2] for t in tris)}", flush=True)
        out, t0, best = [], time.time(), 0
        for i, (a, b, c) in enumerate(tris, 1):
            r = hunt_triangle(a, b, c)
            if r is None or not r.get("gate"):
                continue
            out.append(r)
            best = max(best, r["n_ext"])
            if r["octagons"]:
                print(f"\n*** OCTAGON FOUND *** triangle {(a,b,c)}", flush=True)
                print(json.dumps(r["octagons"][0], indent=1), flush=True)
            if i % 200 == 0 or r["n_ext"] >= 8:
                print(f"  [{i}/{len(tris)}] {(a,b,c)}: E={r['n_ext_raw']}->{r['n_ext']} gp, "
                      f"{r['cliques']} 5-cliques  [{time.time()-t0:.0f}s, max_gp={best}]", flush=True)
            with open("oracle/kbk/engine/octagon_hunt_thin.json", "w") as f:
                json.dump(out, f, indent=1)
        print(f"\n===== THIN HUNT SUMMARY ({time.time()-t0:.0f}s) =====")
        print(f"  triangles hunted : {len(out)}   max |E_gp| : {best}")
        print(f"  OCTAGONS FOUND : {sum(len(r['octagons']) for r in out)}")
        sys.exit(0)

    max_side = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    max_cells = int(sys.argv[2]) if len(sys.argv) > 2 else 60_000_000
    print(f"generating characteristic-{CHAR} triangles with sides <= {max_side} ...", flush=True)
    tris = gen_triangles(max_side)
    # cost of a sweep is ~ (2a+1)(2b+1); hunt cheapest first
    tris = [t for t in tris if (2*t[0]+1)*(2*t[1]+1) <= max_cells]
    tris.sort(key=lambda t: (2*t[0]+1)*(2*t[1]+1))
    print(f"  {len(tris)} triangles of characteristic {CHAR} within the cell budget", flush=True)
    if tris:
        print(f"  cheapest {tris[0]} ({(2*tris[0][0]+1)*(2*tris[0][1]+1):,} cells)"
              f"  dearest {tris[-1]} ({(2*tris[-1][0]+1)*(2*tris[-1][1]+1):,} cells)", flush=True)
    out, t0, best = [], time.time(), 0
    for i, (a, b, c) in enumerate(tris, 1):
        r = hunt_triangle(a, b, c)
        if r is None:
            continue
        out.append(r)
        if r.get("gate"):
            if r["n_ext"] > best:
                best = r["n_ext"]
            if r["octagons"]:
                print(f"\n*** OCTAGON FOUND *** triangle {(a,b,c)}", flush=True)
                print(json.dumps(r["octagons"][0], indent=1), flush=True)
            if i % 25 == 0 or r["n_ext"] >= 5:
                el = time.time() - t0
                print(f"  [{i}/{len(tris)}] {(a,b,c)}: E={r['n_ext_raw']}->{r['n_ext']} gp, "
                      f"{r['cliques']} 5-cliques, oct={len(r['octagons'])}  "
                      f"[{el:.0f}s, max_gp_seen={best}]", flush=True)
        with open("oracle/kbk/engine/octagon_hunt.json", "w") as f:
            json.dump(out, f, indent=1)
    tot = sum(len(r.get("octagons", [])) for r in out)
    print(f"\n===== HUNT SUMMARY ({time.time()-t0:.0f}s) =====")
    print(f"  triangles hunted : {len(out)}")
    print(f"  max |E_gp(T)| seen : {best}   (an octagon needs >= 5)")
    print(f"  OCTAGONS FOUND : {tot}")
