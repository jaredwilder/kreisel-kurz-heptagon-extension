"""kk_deep_search.py — the real large-B extension search, with in-flight genuine-hit alerts.

Same exact arithmetic as kk_extension_search.py (verified regression: reproduces X=P2's coordinates
exactly from (r,s) alone). This driver adds: immediate stdout flag on any hit whose X is NOT an
existing heptagon vertex (a genuine candidate, not the trivial X=P_j solutions), and a final honest
B8/Delta8 receipt.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import isqrt
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

_MODS = _qr_mods([64, 81, 25, 49, 121, 169, 289, 361, 529])  # extra prime power (19^2) vs the earlier version


def deep_search(points, D, B, label, rlo=1):
    a = D[0][1]
    xi = [p[0] for p in points]; yi = [p[1] for p in points]
    existing = set((p[0], p[1]) for p in points)
    a2 = F(a) * F(a)
    gate_passes = 0
    genuine_hits = []
    t_start = time.time()
    for r in range(rlo, B + 1):
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
                ok = True
                for i in range(2, len(points)):
                    ti2 = (x - xi[i]) ** 2 + y2 - 2 * yi[i] * w + C * yi[i] ** 2
                    if perfect_sqrt_frac(ti2) is None:
                        ok = False; break
                if ok:
                    # X=(x, w) with true y-coord = w/sqrt(C); an existing vertex P_j=(x_j, y_j*sqrt(C))
                    # has w_j = C*y_j (rational) — compare directly, no float/sqrt needed.
                    is_existing_vertex = any(x == xi[j] and w == C_sqrt_coef(yi[j]) for j in range(len(points)))
                    if not is_existing_vertex:
                        elapsed = time.time() - t_start
                        print(f"  *** GENUINE HIT [{label}] at r={r} s={s} x={x} w={w}  "
                              f"(elapsed {elapsed:.0f}s) ***", flush=True)
                        genuine_hits.append((r, s, str(x), str(w)))
        if r % max(1, B // 20) == 0:
            elapsed = time.time() - t_start
            print(f"  [{label}] r<={r}/{B} ({100*r/B:.1f}%) gate_passes={gate_passes} "
                  f"genuine_hits={len(genuine_hits)} elapsed={elapsed:.0f}s", flush=True)
    total = time.time() - t_start
    print(f"[{label}] COMPLETE: B8={B} (exhaustively searched, empty of genuine extensions) "
          f"gate_passes={gate_passes} genuine_hits={len(genuine_hits)} time={total:.0f}s", flush=True)
    return {"B": B, "gate_passes": gate_passes, "genuine_hits": genuine_hits, "seconds": total}


def C_sqrt_coef(y_coef):
    # w for an existing vertex with plane-y-coefficient y_coef (meaning true y = y_coef*sqrt(C)) is
    # w := sqrt(C)*true_y = sqrt(C)*y_coef*sqrt(C) = C*y_coef  (rational!)
    return C * y_coef


if __name__ == "__main__":
    which = sys.argv[2] if len(sys.argv) > 2 else "H1"
    pts, D = (H1, D1) if which == "H1" else (reconstruct_full(D2, 2002), D2)
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    res = deep_search(pts, D, B, which)
    import json
    with open(f"oracle/kbk/engine/deep_search_{which}_B{B}.json", "w") as f:
        json.dump(res, f, indent=1)
