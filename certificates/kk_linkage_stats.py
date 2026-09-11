"""kk_linkage_stats.py — how independent ARE the 5 linkage square-conditions? (invariant/density hunt)

For every candidate X that passes the field gate (X in Q(sqrt2002), i.e. r,s integers AND 2002*Heron16 a
perfect square), the two anchor distances t_0=r, t_1=s are integers BY CONSTRUCTION. The remaining five
conditions -- t_i^2 = Q_i(x,w)/2002 a perfect square for i=2..6 -- are what must hold simultaneously.

This measures the JOINT distribution: over all gate-passing X in a band, how many of the 5 are squares?
  * If the max ever seen is 1-2, the conditions behave independently -> a density heuristic applies and
    predicts the height at which an extension could first appear. That quantifies whether ANY bounded
    search can succeed.
  * If some X hits 4/5, that is a genuine near-miss and a target for closer study.

Deliberately reports the EMPIRICAL per-condition square rate rather than assuming the textbook 1/t density.

Run: python oracle/kbk/engine/kk_linkage_stats.py <B> [H1|H2]
"""
from __future__ import annotations
from fractions import Fraction as F
from math import isqrt
import json
import sys
import time
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


def collect(points, D, B, label):
    a = D[0][1]
    xi = [p[0] for p in points]; yi = [p[1] for p in points]
    a2 = F(a) * F(a)
    n = len(points)
    hist = {k: 0 for k in range(0, n - 1)}   # how many of the 5 linkage conds are squares
    best = {"count": -1, "at": None}
    gate_passes = 0
    per_condition_hits = [0] * n
    t0 = time.time()
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
            gate_passes += 1
            x = (a2 + F(r) * F(r) - F(s) * F(s)) / (2 * F(a))
            y2 = F(N16, (2 * a) ** 2)
            for w in ({F(gr, 2 * a), F(-gr, 2 * a)} if gr != 0 else {F(0)}):
                cnt = 0
                for i in range(2, n):
                    ti2 = (x - xi[i]) ** 2 + y2 - 2 * yi[i] * w + C * yi[i] ** 2
                    if perfect_sqrt_frac(ti2) is not None:
                        cnt += 1
                        per_condition_hits[i] += 1
                hist[cnt] = hist.get(cnt, 0) + 1
                if cnt > best["count"]:
                    best = {"count": cnt, "at": (r, s, str(w))}
        if r % max(1, B // 10) == 0:
            print(f"  [{label}] r<={r}/{B} gate={gate_passes} best_simultaneous={best['count']} "
                  f"elapsed={time.time()-t0:.0f}s", flush=True)
    total_evals = sum(hist.values())
    print(f"\n[{label}] B={B}: gate_passes={gate_passes}, X-evaluations(with both w signs)={total_evals}")
    print(f"  histogram of (# of the 5 linkage conditions that are perfect squares):")
    for k in sorted(hist):
        if hist[k]:
            print(f"    {k} of 5 : {hist[k]}")
    print(f"  MAX simultaneous = {best['count']} at (r,s,w)={best['at']}")
    print(f"  per-condition square counts i=2..6: {per_condition_hits[2:]}")
    if total_evals:
        rate = sum(per_condition_hits[2:]) / (5 * total_evals)
        print(f"  EMPIRICAL per-condition square rate p = {rate:.6g}")
        if rate > 0:
            import math
            print(f"  independence heuristic: P(all 5) ~ p^5 = {rate**5:.3g};")
            print(f"    expected X-evals needed for one hit ~ {1/rate**5:.3g}")
            print(f"    (observed {total_evals} evals for B={B} => evals grow ~linearly in B, so the")
            print(f"     heuristic height needed is astronomically beyond any feasible B.)")
    return {"B": B, "gate_passes": gate_passes, "hist": {str(k): v for k, v in hist.items()},
            "max_simultaneous": best["count"], "best_at": best["at"],
            "per_condition_hits": per_condition_hits[2:], "total_evals": total_evals}


if __name__ == "__main__":
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 40000
    which = sys.argv[2] if len(sys.argv) > 2 else "H1"
    pts, D = (H1, D1) if which == "H1" else (reconstruct_full(D2, 2002), D2)
    res = collect(pts, D, B, which)
    with open(f"oracle/kbk/engine/linkage_stats_{which}_B{B}.json", "w") as f:
        json.dump(res, f, indent=1)
    print(f"  receipt -> oracle/kbk/engine/linkage_stats_{which}_B{B}.json")
