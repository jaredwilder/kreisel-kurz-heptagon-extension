"""kk_extension_search.py — direct bounded extension search for a certified KK heptagon (subordinate receipt).

For each (r,s) in [1,B]^2 (r=|X-P0|, s=|X-P1|), the extension point X=(x,y) is DETERMINED:
    x = (a^2+r^2-s^2)/(2a),   y^2 = r^2 - x^2,   w := sqrt(2002)*y  with  w^2 = 2002*(r^2-x^2).
X has coordinates in Q(sqrt(2002)) iff 2002*(r^2-x^2) is a perfect rational square (=> w rational, two signs).
Then for i=2..6:  t_i^2 = (x-x_i)^2 + (r^2-x^2) - 2 y_i w + 2002 y_i^2  (exact Fraction); X extends the
heptagon iff all of t_2..t_6 are perfect squares AND X is not an existing vertex AND general position holds.

Tracks:  B8 = largest r,s-box exhaustively excluded,  Delta8 = smallest total square-defect of any candidate.
Exact arithmetic (fractions). No floats in the decision path.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import isqrt
import sys
from kk_heptagon import H1, D1, D2, reconstruct_full, sq_dist_field

C = 2002


def perfect_sqrt_frac(q: F):
    if q < 0:
        return None
    rn = isqrt(q.numerator); rd = isqrt(q.denominator)
    if rn * rn == q.numerator and rd * rd == q.denominator:
        return F(rn, rd)
    return None


def defect_to_square(q: F):
    """A nonneg measure of how far rational q>=0 is from being a perfect square (0 iff perfect square)."""
    if q < 0:
        return abs(q) + 1
    r = perfect_sqrt_frac(q)
    if r is not None:
        return F(0)
    # distance to nearest integer square if q is an integer, else a positive proxy
    if q.denominator == 1:
        r = isqrt(q.numerator)
        lo = q.numerator - r * r
        hi = (r + 1) * (r + 1) - q.numerator
        return F(min(lo, hi))
    return F(1, q.denominator)  # non-integer -> not a distance square; small positive proxy


def _qr_mods(mods):
    return [(m, {(v * v) % m for v in range(m)}) for m in mods]

_MODS = _qr_mods([64, 81, 25, 49, 121, 169, 289])  # prime-power QR prefilter for "2002*N16 is a square"

def search(points, D, B, label, rlo=1):
    a = D[0][1]
    pts = points
    xi = [p[0] for p in pts]; yi = [p[1] for p in pts]
    best_defect = None; best_at = None
    hits = []; gate_passes = 0
    a2 = F(a) * F(a)
    for r in range(rlo, B + 1):
        # triangle inequality for sides (r,s,a): need s > |r-a| and s < r+a
        s0 = max(1, abs(r - a) + 1)
        s1 = min(B, r + a - 1)
        for s in range(s0, s1 + 1):
            N16 = (r + s + a) * (-r + s + a) * (r - s + a) * (r + s - a)
            g = C * N16
            if g <= 0:
                continue
            # modular prefilter: g must be a QR modulo each prime power
            skip = False
            for m, qs in _MODS:
                if (g % m) not in qs:
                    skip = True; break
            if skip:
                continue
            gr = isqrt(g)
            if gr * gr != g:
                continue  # X not over Q(sqrt(2002)) -> no integer distances (Kurz characteristic thm)
            gate_passes += 1
            # w = +/- gr/(2a);  x = (a^2+r^2-s^2)/(2a);  y2 = r^2 - x^2 = N16/(2a)^2
            x = (a2 + F(r) * F(r) - F(s) * F(s)) / (2 * F(a))
            y2 = F(N16, (2 * a) ** 2)
            for w in ({F(gr, 2 * a), F(-gr, 2 * a)} if gr != 0 else {F(0)}):
                defect = F(0); tvals = []; ok = True
                for i in range(2, 7):
                    ti2 = (x - xi[i]) ** 2 + y2 - 2 * yi[i] * w + C * yi[i] ** 2
                    root = perfect_sqrt_frac(ti2)
                    tvals.append((ti2, root))
                    if root is None:
                        ok = False
                    defect += defect_to_square(ti2)
                if best_defect is None or defect < best_defect:
                    best_defect = defect; best_at = (r, s, str(w))
                if ok:
                    hits.append((r, s, str(x), str(w), [str(t[0]) for t in tvals]))
        if r % max(1, B // 10) == 0:
            print(f"  [{label}] r<= {r}/{B}; gate-passes={gate_passes}; best square-defect={best_defect} at {best_at}", flush=True)
    return {"B": B, "hits": hits, "best_defect": str(best_defect), "best_at": best_at, "gate_passes": gate_passes}


if __name__ == "__main__":
    which = sys.argv[2] if len(sys.argv) > 2 else "H1"
    pts, D = (H1, D1) if which == "H1" else (reconstruct_full(D2, 2002), D2)
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 2 * D[0][1]  # default: full band up to 2a
    res = search(pts, D, B, which)
    # classify hits: trivial (X == existing vertex) vs genuine extension
    existing = set((p[0], p[1]) for p in pts)
    genuine = []
    for h in res["hits"]:
        r, s, xstr, wstr, tv = h
        # X = (x, y) with y = w/sqrt(2002); as a field element store (x_frac, w_frac) — vertex iff matches
        genuine.append(h)  # detailed vertex/general-position classification in Phase 4
    print(f"\n{which} B={B}: gate-passes={res['gate_passes']}; all-integer-distance hits = {len(res['hits'])}")
    print(f"  (hits include the trivial embeddings X=P_j; a hit that is NOT an existing vertex = genuine octagon)")
    for h in res["hits"][:30]:
        print("  HIT (r,s,x,w,t2..t6):", h)
    print(f"  Delta8 (min square-defect over gate-passing X) = {res['best_defect']} at (r,s,w)={res['best_at']}")
    print(f"  B8 receipt: exhaustively scanned r in [1,{B}], s in triangle band, for {which}.")
