"""kk_scale_corrected_hunt.py — the hunt at the scale the evidence says it must run.

WHY THIS EXISTS (min-distance scale separation, 2026-07-25). Measured minimum distance of the extremal
integral general-position n-set:

    n         3    4    5    6      7
    min dist  1    5    26   68     5780        (H1; H2's is 16637)
    growth         5.0x 5.2x 2.6x   85x

Every step is 2.6-5.2x except n=6 -> 7, which is 85x. Every edge of both known heptagons is >= 5780.
All earlier hunts used seeds with sides <= 260 (sweep) or <= 900 (contract) -- 6 to 22x BELOW the
smallest edge any heptagon has. Those seeds could not have been sub-triangles of a 7-set, so the
observed cap at 5-6 was a scale error, not an obstruction.

COROLLARY (a compute floor, not an optimisation target): a sweep costs (2a+1)(2b+1) in the two
smallest sides, so `gen_triangles_thin` makes seeds cheap by making one side small -- but a small side
is exactly what the wedge forbids. The cheap seeds are precisely the useless ones. Any viable seed
costs >= (2*5780)^2 ~ 1.3e8 cells. That floor is inherent.

This hunts char-2002 triangles with ALL sides >= min_side, which is where a 7-set (and therefore any
8-set containing one) can actually live.

Run: python oracle/kbk/engine/kk_scale_corrected_hunt.py <min_side> <max_side> <n_seeds>
"""
from __future__ import annotations
import json
import sys
import time
from fractions import Fraction as F
from math import isqrt

sys.path.insert(0, "oracle/kbk/engine")
from kk_maxset import grow
from kk_octagon_hunt import triangle_points

CHAR = 2002


def char_of(a: int, b: int, c: int) -> int | None:
    """squarefree(16*Area^2) via Heron, or None if degenerate. 16A^2 = 2a2b2+2b2c2+2c2a2-a4-b4-c4."""
    a2, b2, c2 = a * a, b * b, c * c
    h = 2 * a2 * b2 + 2 * b2 * c2 + 2 * c2 * a2 - a2 * a2 - b2 * b2 - c2 * c2
    if h <= 0:
        return None
    # squarefree part
    m, res, d = h, 1, 2
    while d * d <= m:
        e = 0
        while m % d == 0:
            m //= d; e += 1
        if e % 2:
            res *= d
        d += 1
    return res * m if m > 1 else res


def gen_seeds(min_side: int, max_side: int, char: int, want: int):
    """Char-`char` triangles with min_side <= a <= b <= c <= max_side.

    FAST TEST (from gen_triangles): squarefree(h) == char  <=>  h % char == 0 AND h/char is a perfect
    square, where h = 16*Area^2. No factorisation. Vectorised over c with numpy + Python ints for the
    perfect-square check (h ~ max_side^4 overflows int64 above ~55000, so we use object dtype)."""
    import numpy as np
    found = []
    for a in range(min_side, max_side + 1):
        a2 = a * a
        for b in range(a, max_side + 1):
            lo, hi = max(b, min_side), min(max_side, a + b - 1)
            if lo > hi:
                continue
            b2 = b * b
            c = np.arange(lo, hi + 1, dtype=object)
            c2 = c * c
            h = 2 * a2 * b2 + 2 * b2 * c2 + 2 * c2 * a2 - a2 * a2 - b2 * b2 - c2 * c2
            for idx in range(len(c)):
                hv = int(h[idx])
                if hv <= 0 or hv % char:
                    continue
                q = hv // char
                r = isqrt(q)
                if r * r == q:
                    found.append((a, b, int(c[idx])))
                    if len(found) >= want:
                        return found
    return found


if __name__ == "__main__":
    min_side = int(sys.argv[1]) if len(sys.argv) > 1 else 5780
    max_side = int(sys.argv[2]) if len(sys.argv) > 2 else 12000
    want = int(sys.argv[3]) if len(sys.argv) > 3 else 40

    print(f"[scale-corrected] char={CHAR}, sides in [{min_side},{max_side}] -- at or above the "
          f"5780 heptagon-edge floor", flush=True)
    t0 = time.time()
    seeds = gen_seeds(min_side, max_side, CHAR, want)
    print(f"[scale-corrected] {len(seeds)} viable seeds found in {time.time()-t0:.0f}s", flush=True)
    if not seeds:
        print("no char-2002 triangles in this band"); sys.exit(0)

    seeds.sort(key=lambda t: (2 * t[0] + 1) * (2 * t[1] + 1))
    out, best = [], 0
    for i, (a, b, c) in enumerate(seeds, 1):
        pts = triangle_points(a, b, c, CHAR)
        if pts is None:
            continue
        D = [[0] * 3 for _ in range(3)]
        D[0][1] = D[1][0] = a; D[0][2] = D[2][0] = b; D[1][2] = D[2][1] = c
        cells = (2 * a + 1) * (2 * b + 1)
        t1 = time.time()
        total, wit, nx = grow(pts, D, CHAR, f"sc-{a}_{b}_{c}")
        if total is None:
            continue
        out.append({"tri": [a, b, c], "char": CHAR, "max_total": total,
                    "n_ext_gp": nx, "cells": cells, "witness": wit})
        tag = ""
        if total > best:
            best = total
            tag = "  <== NEW BEST"
            if total >= 7:
                print(f"  *** {'OCTAGON' if total >= 8 else 'SEVEN-POINT SET'} at {(a,b,c)} ***",
                      flush=True)
                print(f"      {json.dumps(wit)}", flush=True)
        print(f"  [{i}/{len(seeds)}] {(a,b,c)} cells={cells:,} |E_gp|={nx} max={total} "
              f"({time.time()-t1:.0f}s){tag}", flush=True)
        with open("oracle/kbk/engine/scale_corrected_hunt.json", "w") as f:
            json.dump(out, f, indent=1)
    print(f"\n[scale-corrected] BEST max-set = {best} over {len(out)} seeds "
          f"({time.time()-t0:.0f}s)   (7 = new heptagon, 8 = OCTAGON)")
