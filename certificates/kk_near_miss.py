"""kk_near_miss.py — isolate the NEAR-MISS candidates: X satisfying >= 3 of the 5 linkage conditions.

The linkage-statistics scan found, in the band r <= 40000 for H1, a candidate hitting 4 of the 5
conditions (alongside the 5 genuine trivial embeddings which hit 5 of 5). A 4-of-5 point is an
"almost octagon": seven of the eight required distances are simultaneously integral and only one fails.
Those are the concrete objects worth studying -- exactly one arithmetic condition away from settling n=8.

This isolates them exactly, reports WHICH condition fails and by how much, and re-verifies each.

Run: python oracle/kbk/engine/kk_near_miss.py <B> [H1|H2] [minhits]
"""
from __future__ import annotations
from fractions import Fraction as F
from math import isqrt
import json
import sys
sys.path.insert(0, "oracle/kbk/engine")
from kk_heptagon import H1, D1, D2, reconstruct_full

C = 2002


def perfect_sqrt_frac(q: F):
    if q < 0:
        return None
    rn = isqrt(q.numerator); rd = isqrt(q.denominator)
    return F(rn, rd) if rn * rn == q.numerator and rd * rd == q.denominator else None


def _qr_mods(mods):
    return [(m, {(v * v) % m for v in range(m)}) for m in mods]

_MODS = _qr_mods([64, 81, 25, 49, 121, 169, 289, 361, 529])


def find(points, D, B, label, minhits=3):
    a = D[0][1]
    xi = [p[0] for p in points]; yi = [p[1] for p in points]
    a2 = F(a) * F(a)
    n = len(points)
    # the trivial solutions, as (x, w) pairs, to label them
    trivial_xw = {(xi[j], C * yi[j]) for j in range(n)}
    found = []
    for r in range(1, B + 1):
        s0 = max(1, abs(r - a) + 1)
        s1 = min(B, r + a - 1)
        if s0 > s1:
            continue
        for s in range(s0, s1 + 1):
            N16 = (r + s + a) * (-r + s + a) * (r - s + a) * (r + s - a)
            g = C * N16
            if g <= 0:
                continue
            skip = False
            for m, qs in _MODS:
                if (g % m) not in qs:
                    skip = True; break
            if skip:
                continue
            gr = isqrt(g)
            if gr * gr != g:
                continue
            x = (a2 + F(r) * F(r) - F(s) * F(s)) / (2 * F(a))
            y2 = F(N16, (2 * a) ** 2)
            for w in ({F(gr, 2 * a), F(-gr, 2 * a)} if gr != 0 else {F(0)}):
                sq, fail = [], []
                for i in range(2, n):
                    ti2 = (x - xi[i]) ** 2 + y2 - 2 * yi[i] * w + C * yi[i] ** 2
                    root = perfect_sqrt_frac(ti2)
                    (sq if root is not None else fail).append((i, ti2, root))
                if len(sq) >= minhits:
                    is_trivial = (x, w) in trivial_xw
                    found.append({
                        "r": r, "s": s, "x": str(x), "w": str(w),
                        "hits": len(sq), "trivial_vertex": is_trivial,
                        "satisfied": [(i, str(root)) for i, _, root in sq],
                        "failed": [(i, str(t2)) for i, t2, _ in fail],
                    })
                    tag = "TRIVIAL (X = P_j)" if is_trivial else "*** NON-TRIVIAL NEAR-MISS ***"
                    print(f"  {len(sq)}/5 at r={r} s={s}  {tag}", flush=True)
                    if not is_trivial:
                        print(f"      x={x}  w={w}")
                        print(f"      satisfied: {[(i, str(rt)) for i, _, rt in sq]}")
                        for i, t2, _ in fail:
                            # how far from a square is the failing one?
                            nearest = None
                            if t2 >= 0 and t2.denominator == 1:
                                rt = isqrt(t2.numerator)
                                nearest = f"between {rt}^2 and {rt+1}^2 (gap to lower {t2.numerator - rt*rt})"
                            print(f"      FAILED i={i}: t_{i}^2 = {t2}   {nearest or '(non-integer square)'}")
    return found


if __name__ == "__main__":
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 40000
    which = sys.argv[2] if len(sys.argv) > 2 else "H1"
    minhits = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    pts, D = (H1, D1) if which == "H1" else (reconstruct_full(D2, 2002), D2)
    print(f"[{which}] scanning r<={B} for candidates with >= {minhits} of 5 linkage conditions satisfied")
    res = find(pts, D, B, which, minhits)
    nontrivial = [f for f in res if not f["trivial_vertex"]]
    print(f"\ntotal >= {minhits}/5 : {len(res)}   of which NON-TRIVIAL: {len(nontrivial)}")
    with open(f"oracle/kbk/engine/near_miss_{which}_B{B}.json", "w") as f:
        json.dump(res, f, indent=1)
    print(f"receipt -> oracle/kbk/engine/near_miss_{which}_B{B}.json")
