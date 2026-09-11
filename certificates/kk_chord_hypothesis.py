"""kk_chord_hypothesis.py — TEST the chord-concurrence hypothesis on every partial near-miss.

HYPOTHESIS (from the single 4-of-5 near-miss found at r=17816, s=22794, which turned out to lie on
THREE chords simultaneously: P0P5, P1P3, P4P6):

    Non-trivial candidates that satisfy MANY of the linkage square-conditions are chord points --
    they lie on lines through pairs of heptagon vertices. Integrality is inherited from the
    collinearity relations, and collinearity is exactly what general position forbids. I.e. the
    arithmetic conditions and the general-position conditions are ANTI-CORRELATED.

If true, this is a real structural obstruction mechanism (and a pruning rule: chord points can be
skipped, and they are precisely the arithmetically-rich ones).

This scans for ALL non-trivial candidates with >= `minhits` of 5 conditions and, for each, counts how
many heptagon chords it lies on. Reports the correlation between linkage-hits and chord-membership,
including the base rate over ordinary gate-passing candidates (so the claim is measured, not assumed).

Run: python oracle/kbk/engine/kk_chord_hypothesis.py <B> [H1|H2] [minhits]
"""
from __future__ import annotations
from fractions import Fraction as F
from math import isqrt
from itertools import combinations
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


def chords_through(pt, points):
    """Which vertex-pairs (i,j) have pt collinear with P_i, P_j? pt and points in (x, y-coeff) form."""
    out = []
    for i, j in combinations(range(len(points)), 2):
        a, b = points[i], points[j]
        if (b[0] - a[0]) * (pt[1] - a[1]) - (b[1] - a[1]) * (pt[0] - a[0]) == 0:
            out.append((i, j))
    return out


def scan(points, D, B, label, minhits=2):
    a = D[0][1]
    xi = [p[0] for p in points]; yi = [p[1] for p in points]
    a2 = F(a) * F(a)
    n = len(points)
    trivial_xw = {(xi[j], C * yi[j]) for j in range(n)}
    rows = []
    base_chord_count = 0     # chord-membership among ALL gate-passing non-trivial candidates
    base_total = 0
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
                if (x, w) in trivial_xw:
                    continue
                hits = 0
                for i in range(2, n):
                    ti2 = (x - xi[i]) ** 2 + y2 - 2 * yi[i] * w + C * yi[i] ** 2
                    if perfect_sqrt_frac(ti2) is not None:
                        hits += 1
                pt = (x, w / C)      # convert w back to a y-coefficient
                ch = chords_through(pt, points)
                base_total += 1
                if ch:
                    base_chord_count += 1
                if hits >= minhits:
                    rows.append({"r": r, "s": s, "x": str(x), "w": str(w),
                                 "linkage_hits": hits, "n_chords": len(ch), "chords": ch})
                    print(f"  hits={hits}/5  r={r} s={s}  chords={len(ch)} {ch}", flush=True)
    print(f"\n[{label}] B={B}")
    print(f"  BASE RATE: non-trivial gate-passing candidates on >=1 chord: "
          f"{base_chord_count}/{base_total}"
          f"{' = %.4g%%' % (100*base_chord_count/base_total) if base_total else ''}")
    if rows:
        onch = sum(1 for r_ in rows if r_["n_chords"] > 0)
        print(f"  HIGH-LINKAGE (>= {minhits}/5) candidates on >=1 chord: {onch}/{len(rows)}"
              f" = {100*onch/len(rows):.4g}%")
        for k in sorted({r_["linkage_hits"] for r_ in rows}, reverse=True):
            grp = [r_ for r_ in rows if r_["linkage_hits"] == k]
            g_on = sum(1 for r_ in grp if r_["n_chords"] > 0)
            print(f"    {k}/5 : {g_on}/{len(grp)} on a chord "
                  f"(mean #chords {sum(r_['n_chords'] for r_ in grp)/len(grp):.2f})")
    return {"B": B, "base_chord": base_chord_count, "base_total": base_total, "rows": rows}


if __name__ == "__main__":
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 40000
    which = sys.argv[2] if len(sys.argv) > 2 else "H1"
    minhits = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    pts, D = (H1, D1) if which == "H1" else (reconstruct_full(D2, 2002), D2)
    print(f"[{which}] chord-hypothesis test, r<={B}, reporting non-trivial candidates with >={minhits}/5")
    res = scan(pts, D, B, which, minhits)
    with open(f"oracle/kbk/engine/chord_hypothesis_{which}_B{B}.json", "w") as f:
        json.dump(res, f, indent=1)
    print(f"receipt -> oracle/kbk/engine/chord_hypothesis_{which}_B{B}.json")
