"""kk_scaled_search.py — search the DENOMINATOR axis: rational-distance extensions of a KK heptagon.

WHY THIS EXISTS (the gap in the earlier program). An 8-point set with all 28 distances RATIONAL scales
by the common denominator to an INTEGRAL 8-point set, and general position is scale-invariant. So

    exists integral general-position 8-set  <=>  exists rational-distance general-position 8-set.

An integer-(r,s) search against H1 therefore only finds DENOMINATOR-1 extensions. A rational extension
of H1 with denominator q is exactly an INTEGER extension of the scaled heptagon q*H1 (whose characteristic
is still 2002 -- verified, since 16*Area^2 scales by q^4). The earlier B=100000 runs were the q=1 slice only.

This driver runs the same verified exact machinery against q*H1 for q = 2,3,... -- a genuinely different
arithmetic direction, not a larger box.

Run: python oracle/kbk/engine/kk_scaled_search.py <B> <q> [H1|H2]
"""
from __future__ import annotations
import json
import sys
sys.path.insert(0, "oracle/kbk/engine")
from kk_heptagon import H1, D1, D2, reconstruct_full
from kk_deep_search import deep_search


def scaled(points, D, q):
    pts = [(p[0] * q, p[1] * q) for p in points]
    Dq = [[d * q for d in row] for row in D]
    return pts, Dq


if __name__ == "__main__":
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    q = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    which = sys.argv[3] if len(sys.argv) > 3 else "H1"
    base_pts, base_D = (H1, D1) if which == "H1" else (reconstruct_full(D2, 2002), D2)
    pts, Dq = scaled(base_pts, base_D, q)
    a = Dq[0][1]
    print(f"[{which} x{q}] scaled diameter = {a}; searching integer (r,s) with r <= {B}")
    print(f"  (an integer extension of {q}*{which} == a denominator-{q} RATIONAL extension of {which})")
    if B <= a // 2:
        print(f"  WARNING: B={B} <= diameter/2={a//2}; triangle inequality r+s>a makes the band nearly empty.")
    res = deep_search(pts, Dq, B, f"{which}x{q}")
    res["scale_q"] = q
    res["scaled_diameter"] = a
    with open(f"oracle/kbk/engine/scaled_search_{which}_q{q}_B{B}.json", "w") as f:
        json.dump(res, f, indent=1)
    print(f"  receipt -> oracle/kbk/engine/scaled_search_{which}_q{q}_B{B}.json")
