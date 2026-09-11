"""kk_maxset.py — GROW every seed to its MAXIMUM. Not the binary 5-clique question.

What was avoided: the hunt only ever asked "is there a 5-clique?" (=> octagon) and reported 0. It never
asked HOW BIG the set actually gets, so there was no distance-to-goal metric at all.

  max total size 8  => OCTAGON. Erdos n=8 settled affirmatively.
  max total size 7  => a NEW INTEGRAL HEPTAGON. Only two are known (Kreisel-Kurz 2008); they explicitly
                       pose "further examples or an infinite family" as open. A third closes that.
  max total size <7 => a real measured distance to the frontier, per seed.

Given a seed T, E(T) is enumerated exactly once (finite, no height bound), pruned by the sound
general-position necessary condition, and then we find the MAXIMUM clique that is mutually at integral
distance AND in general position with T -- branch and bound with incremental checks.

Run: python oracle/kbk/engine/kk_maxset.py <mode> ...
     mode=char  <max_side> <characteristic>      grow every seed of that characteristic
     mode=scale <max_side> <char> <q>            same, on q-scaled seeds (uses the rational<=>integral thm)
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import combinations
import json
import sys
import time

sys.path.insert(0, "oracle/kbk/engine")
from kk_closure_full_exact import full_exact_sweep
from kk_hexagon_octagon import dist_int, general_position
from kk_octagon_hunt import gen_triangles, triangle_points


def grow(pts, D, char, label=""):
    """Largest integral general-position set containing the seed. Returns (max_total, witness)."""
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        res = full_exact_sweep(pts, D, label or "seed", char=char, report=1)
    if not res["correctness_gate_pass"]:
        return None, None, 0
    ext_all = [(F(s["x"]), F(s["y_coeff"])) for s in res["solutions"] if not s["is_vertex"]]
    ext = [X for X in ext_all if general_position(pts + [X], char)[0]]
    n = len(ext)
    adj = {i: set() for i in range(n)}
    for i, j in combinations(range(n), 2):
        if dist_int(ext[i], ext[j], char) is not None:
            adj[i].add(j); adj[j].add(i)

    best = {"k": 0, "set": []}

    def extend(cur, cand):
        if len(cur) > best["k"]:
            best["k"] = len(cur); best["set"] = list(cur)
        if len(cur) + len(cand) <= best["k"]:
            return                      # bound: cannot beat the incumbent
        for idx, v in enumerate(sorted(cand)):
            if len(cur) + (len(cand) - idx) <= best["k"]:
                return
            trial = pts + [ext[i] for i in cur] + [ext[v]]
            if not general_position(trial, char)[0]:
                continue
            extend(cur + [v], (cand & adj[v]) - {v})

    extend([], set(range(n)))
    total = len(pts) + best["k"]
    wit = [[str(ext[i][0]), str(ext[i][1])] for i in best["set"]]
    return total, wit, len(ext)


def run(max_side, char, q=1):
    tris = gen_triangles(max_side, char)
    tris = [t for t in tris if (2 * t[0] * q + 1) * (2 * t[1] * q + 1) <= 20_000_000]
    tris.sort(key=lambda t: (2 * t[0] + 1) * (2 * t[1] + 1))
    tag = f"char{char}" + (f"_q{q}" if q > 1 else "")
    print(f"[{tag}] {len(tris)} seeds (sides<={max_side}, scale q={q})", flush=True)
    out, t0, best_seen = [], time.time(), 0
    for i, (a, b, c) in enumerate(tris, 1):
        pts = triangle_points(a * q, b * q, c * q, char)
        if pts is None:
            continue
        D = [[0] * 3 for _ in range(3)]
        D[0][1] = D[1][0] = a * q; D[0][2] = D[2][0] = b * q; D[1][2] = D[2][1] = c * q
        total, wit, nx = grow(pts, D, char, f"{tag}-{a}_{b}_{c}")
        if total is None:
            continue
        out.append({"tri": [a * q, b * q, c * q], "char": char, "max_total": total,
                    "n_ext_gp": nx, "witness": wit})
        if total > best_seen:
            best_seen = total
            print(f"  NEW BEST max_total={total} at {(a*q,b*q,c*q)} (|E_gp|={nx})  [{time.time()-t0:.0f}s]", flush=True)
            if total >= 7:
                print(f"  *** {'OCTAGON' if total>=8 else 'NEW HEPTAGON'} *** {json.dumps(wit)}", flush=True)
        if i % 250 == 0:
            print(f"  [{i}/{len(tris)}] best_so_far={best_seen}  [{time.time()-t0:.0f}s]", flush=True)
        with open(f"oracle/kbk/engine/maxset_{tag}.json", "w") as f:
            json.dump(out, f, indent=1)
    import collections
    dist = collections.Counter(r["max_total"] for r in out)
    print(f"\n===== [{tag}] MAX-SET DISTRIBUTION ({time.time()-t0:.0f}s) =====")
    for k in sorted(dist):
        print(f"  max total {k}: {dist[k]} seeds")
    print(f"  BEST: {best_seen}   (7 = new heptagon, 8 = OCTAGON)")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "char":
        run(int(sys.argv[2]), int(sys.argv[3]))
    elif mode == "scale":
        run(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
