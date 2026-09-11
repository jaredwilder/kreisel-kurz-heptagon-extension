"""kk_multichar.py — sweep CHARACTERISTICS. The axis never tried.

Every hunt so far fixed characteristic 2002 because both known heptagons have it. That is two data
points. If an octagon exists there is no reason its characteristic is 2002. This sweeps c over many
squarefree values and, for each, measures the largest integral general-position set reachable from a
seed of that characteristic.

max total 8 => OCTAGON.  max total 7 => a NEW HEPTAGON (open per Kreisel-Kurz).

Run: python oracle/kbk/engine/kk_multichar.py <max_side> <n_chars> [cell_budget]
"""
from __future__ import annotations
import json
import sys
import time

sys.path.insert(0, "oracle/kbk/engine")
from kk_maxset import grow
from kk_octagon_hunt import gen_triangles, triangle_points


def squarefree_values(limit):
    """Squarefree integers >= 2, ascending."""
    out = []
    for n in range(2, limit + 1):
        m, ok, d = n, True, 2
        while d * d <= m:
            if m % (d * d) == 0:
                ok = False; break
            d += 1
        if ok:
            out.append(n)
    return out


if __name__ == "__main__":
    max_side = int(sys.argv[1]) if len(sys.argv) > 1 else 220
    n_chars = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    budget = int(sys.argv[3]) if len(sys.argv) > 3 else 6_000_000

    # characteristics: small squarefree values + the two known-good ones
    chars = squarefree_values(120)[:n_chars]
    for extra in (2002, 1001, 4004):
        if extra not in chars:
            chars.append(extra)
    print(f"sweeping {len(chars)} characteristics, seeds with sides <= {max_side}", flush=True)

    results, t0, global_best = [], time.time(), (0, None)
    for ci, char in enumerate(chars, 1):
        tris = gen_triangles(max_side, char)
        tris = [t for t in tris if (2 * t[0] + 1) * (2 * t[1] + 1) <= budget]
        if not tris:
            print(f"  [{ci}/{len(chars)}] c={char}: no seeds", flush=True)
            continue
        tris.sort(key=lambda t: (2 * t[0] + 1) * (2 * t[1] + 1))
        best, best_tri = 0, None
        for (a, b, c) in tris:
            pts = triangle_points(a, b, c, char)
            if pts is None:
                continue
            D = [[0] * 3 for _ in range(3)]
            D[0][1] = D[1][0] = a; D[0][2] = D[2][0] = b; D[1][2] = D[2][1] = c
            total, wit, nx = grow(pts, D, char, f"c{char}")
            if total and total > best:
                best, best_tri = total, (a, b, c)
                if total >= 7:
                    print(f"  *** c={char} {(a,b,c)} reaches {total} — "
                          f"{'OCTAGON' if total >= 8 else 'NEW HEPTAGON'} ***", flush=True)
                    print(f"      {json.dumps(wit)}", flush=True)
        results.append({"char": char, "n_seeds": len(tris), "best_total": best, "best_tri": best_tri})
        if best > global_best[0]:
            global_best = (best, char)
        print(f"  [{ci}/{len(chars)}] c={char}: {len(tris)} seeds, best max-set = {best} "
              f"(tri {best_tri})  [{time.time()-t0:.0f}s]", flush=True)
        with open("oracle/kbk/engine/multichar.json", "w") as f:
            json.dump(results, f, indent=1)

    print(f"\n===== MULTI-CHARACTERISTIC SUMMARY ({time.time()-t0:.0f}s) =====")
    for r in sorted(results, key=lambda r: -r["best_total"])[:15]:
        print(f"  c={r['char']:6d}: {r['n_seeds']:5d} seeds, best {r['best_total']}")
    print(f"  GLOBAL BEST: max-set {global_best[0]} at characteristic {global_best[1]}")
    print(f"  (7 = NEW HEPTAGON, open per Kreisel-Kurz; 8 = OCTAGON, closes Erdos n=8)")
