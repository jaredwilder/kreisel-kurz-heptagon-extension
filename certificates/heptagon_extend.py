"""heptagon_extend.py — NEW BLADE: integral point-set extension / bounded-maximality search in R^2.

Deterministic, $0, zero-LLM, NO NETWORK. Exact rational arithmetic (fractions.Fraction) — no float
ever decides anything. Forged 2026-07-24 to attack the Erdos integral-point-set program: does a given
integral point set in GENERAL POSITION admit a further point at integer distance from all of it?
(Kreisel-Kurz found the smallest general-position integral heptagon; whether it extends — or is
"maximal" — is exactly a bounded extension search.)

THE TARGET (operator vector: heptagon extension / maximality):
  Input  P: a set of points in R^2 whose pairwise distances are all integers, in GENERAL POSITION
            (no 3 collinear, no 4 concyclic), NORMALIZED so P[0]=(0,0) and P[1]=(a,0), a a positive int.
  Ask    : is there a point X with integer distance to EVERY point of P, all distances <= a bound B,
            keeping general position?

TWO-ANCHOR REDUCTION (the exact-arithmetic engine):
  Fix anchors P1=(0,0), P2=(a,0). A candidate X with |X-P1|=r, |X-P2|=s (r,s positive integers) has
      x = (a^2 + r^2 - s^2) / (2a)          -- rational, denominator divides 2a
      y^2 = r^2 - x^2                        -- rational
  X exists as a real point iff y^2 >= 0, and has RATIONAL coordinates iff y^2 is a rational square.
  We enumerate integer (r,s) in [1,B]^2 (finite), keep those giving a rational X off the P1-P2 line,
  then INDEPENDENTLY re-check that the remaining distances |X-P_i| are all integers and that general
  position is preserved. Every accepted X is at all-integer distance to P by construction+recheck.

  SCOPE / HONESTY: this engine works inside Q^2 (rational coordinates). Integer-distance sets do NOT
  all admit rational coordinates (e.g. a bare 2-3-4 triangle does not), so a "none up to B" verdict is
  a claim about RATIONAL-coordinate extensions with all new distances <= B. It is a BOUNDED,
  COMPUTATIONAL statement — never an unconditional "this set is maximal". The bound B is always stated.

  A maximality verdict on an INVALID input is meaningless, so verify_input() runs FIRST: it re-checks
  integrality + general position from scratch (separate code path from the search) before any claim.

GROUND TRUTH (self-test):
  * Genuine integral 4-set in general position, rational coords:
        A=(0,0) B=(8,0) C=(4,3) D=(4,-3)   distances AB=8 AC=AD=BC=BD=5 CD=6, non-concyclic.
    Dropping D leaves triangle T={A,B,C}; the search over T MUST re-discover D=(4,-3). (found-extension)
  * Over the same T with a tight bound B=4, the triangle inequality (|X-A|+|X-B| >= |A-B|=8) forces
    r=s=4 with X on line AB (collinear) -> NO general-position extension. (none-up-to-B)

Self-test:  python oracle/kbk/engine/heptagon_extend.py
Config run: python oracle/kbk/engine/heptagon_extend.py run   (writes heptagon_extend_results.json)
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from itertools import combinations
import json
import math
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# exact rational primitives (no float ever decides)
# ---------------------------------------------------------------------------
def _F(v) -> Fraction:
    return v if isinstance(v, Fraction) else Fraction(v)


def as_pt(pt) -> tuple[Fraction, Fraction]:
    return (_F(pt[0]), _F(pt[1]))


def rational_sqrt(q: Fraction):
    """Exact square root of a nonnegative rational, or None if not a rational square."""
    if q < 0:
        return None
    if q == 0:
        return Fraction(0)
    n, d = q.numerator, q.denominator  # d > 0, lowest terms
    rn, rd = math.isqrt(n), math.isqrt(d)
    if rn * rn == n and rd * rd == d:
        return Fraction(rn, rd)
    return None


def integer_distance(a, b):
    """Integer distance between two rational points, or None if the distance is not an integer."""
    dx = _F(a[0]) - _F(b[0])
    dy = _F(a[1]) - _F(b[1])
    d2 = dx * dx + dy * dy  # Fraction
    if d2.denominator != 1:
        return None
    n = d2.numerator
    r = math.isqrt(n)
    return r if r * r == n else None


def collinear(a, b, c) -> bool:
    """Three rational points collinear iff the affine cross product is exactly 0."""
    cross = (_F(b[0]) - _F(a[0])) * (_F(c[1]) - _F(a[1])) \
        - (_F(b[1]) - _F(a[1])) * (_F(c[0]) - _F(a[0]))
    return cross == 0


def _det4(m) -> Fraction:
    """Exact 4x4 determinant over the rationals via cofactor expansion (no float)."""
    def det3(mm):
        return (
            mm[0][0] * (mm[1][1] * mm[2][2] - mm[1][2] * mm[2][1])
            - mm[0][1] * (mm[1][0] * mm[2][2] - mm[1][2] * mm[2][0])
            + mm[0][2] * (mm[1][0] * mm[2][1] - mm[1][1] * mm[2][0])
        )
    total = Fraction(0)
    for j in range(4):
        minor = [[m[i][k] for k in range(4) if k != j] for i in range(1, 4)]
        total += ((-1) ** j) * m[0][j] * det3(minor)
    return total


def cocircular(a, b, c, d) -> bool:
    """Four rational points concyclic (or all collinear) iff the lift determinant is 0.
    Rows [x^2+y^2, x, y, 1]. Collinearity is forbidden separately, so det==0 => a genuine circle."""
    def row(pt):
        x, y = _F(pt[0]), _F(pt[1])
        return [x * x + y * y, x, y, Fraction(1)]
    return _det4([row(a), row(b), row(c), row(d)]) == 0


# ---------------------------------------------------------------------------
# independent verifier #1: is the INPUT a genuine general-position integral set?
# (separate code path from the search — no maximality claim is valid without this)
# ---------------------------------------------------------------------------
def verify_input(P) -> tuple[bool, str]:
    pts = [as_pt(p) for p in P]
    if len(set(pts)) != len(pts):
        return False, "duplicate point"
    if len(pts) < 2:
        return False, "need >= 2 points to anchor"
    # normalization required by the two-anchor reduction
    if pts[0] != (Fraction(0), Fraction(0)):
        return False, "not normalized: P[0] must be (0,0)"
    a = pts[1]
    if a[1] != 0 or a[0] <= 0 or a[0].denominator != 1:
        return False, "not normalized: P[1] must be (a,0) with a a positive integer"
    for x, y in combinations(pts, 2):
        if integer_distance(x, y) is None:
            return False, f"non-integer distance between {x},{y}"
    for x, y, z in combinations(pts, 3):
        if collinear(x, y, z):
            return False, f"collinear triple {x},{y},{z}"
    for x, y, z, w in combinations(pts, 4):
        if cocircular(x, y, z, w):
            return False, f"cocircular quadruple {x},{y},{z},{w}"
    return True, "ok"


# ---------------------------------------------------------------------------
# independent verifier #2: does candidate X genuinely extend P within bound B?
# (separate code path — re-derives every distance and every general-position constraint)
# ---------------------------------------------------------------------------
def verify_extension(P, X, B: int) -> tuple[bool, str]:
    pts = [as_pt(p) for p in P]
    X = as_pt(X)
    if X in pts:
        return False, "X duplicates an existing point"
    for p in pts:
        d = integer_distance(X, p)
        if d is None:
            return False, f"non-integer distance X->{p}"
        if d > B:
            return False, f"distance X->{p} = {d} exceeds bound B={B}"
        if d == 0:
            return False, "zero distance"
    for a, b in combinations(pts, 2):
        if collinear(a, b, X):
            return False, f"X collinear with {a},{b}"
    for a, b, c in combinations(pts, 3):
        if cocircular(a, b, c, X):
            return False, f"X cocircular with {a},{b},{c}"
    return True, "ok"


# ---------------------------------------------------------------------------
# the extension search (two-anchor reduction). Exhaustive over rational X whose
# distances to the two anchors are integers in [1,B]. Every hit is re-checked by
# the independent verifier above before it is accepted.
# ---------------------------------------------------------------------------
def extension_search(P, B: int):
    """Return the sorted list of all rational points X that extend P with every new distance an
    integer <= B, keeping general position. [] means: maximal among such X, up to bound B."""
    pts = [as_pt(p) for p in P]
    a = pts[1][0]                       # P[1] = (a, 0), a a positive integer (checked by verify_input)
    ai = a.numerator
    found = []
    seen = set()
    for r in range(1, B + 1):
        for s in range(1, B + 1):
            x = Fraction(ai * ai + r * r - s * s, 2 * ai)
            y2 = Fraction(r * r) - x * x
            if y2 < 0:
                continue
            ys = rational_sqrt(y2)
            if ys is None:                # X not rational -> outside this engine's Q^2 scope
                continue
            for y in ({ys, -ys}):         # both branches; y==0 is on the anchor line (rejected below)
                X = (x, y)
                if X in seen:
                    continue
                ok, _ = verify_extension(pts, X, B)
                if ok:
                    seen.add(X)
                    found.append(X)
    return sorted(found)


# ---------------------------------------------------------------------------
# serialization helpers (exact — Fractions become strings like "4" or "7/2")
# ---------------------------------------------------------------------------
def _pt_str(pt):
    return [str(_F(pt[0])), str(_F(pt[1]))]


def analyze(name: str, P, B: int) -> dict:
    ok_in, reason = verify_input(P)
    row = {
        "config": name,
        "n_points": len(P),
        "points": [_pt_str(p) for p in P],
        "verified_input": ok_in,
        "input_reason": reason,
        "bound_B": B,
    }
    if not ok_in:
        row["extension_search"] = "SKIPPED (invalid input — no maximality claim possible)"
        return row
    exts = extension_search(P, B)
    row["n_extensions_found"] = len(exts)
    row["extensions"] = [
        {"X": _pt_str(X), "distances": [integer_distance(X, p) for p in P]} for X in exts
    ]
    if exts:
        row["verdict"] = f"EXTENDABLE within bound B={B}"
    else:
        row["verdict"] = (
            f"NO rational-coordinate extension with all new distances <= B={B} "
            f"(BOUNDED/COMPUTATIONAL maximality; NOT an unconditional maximality claim)"
        )
    return row


# ---------------------------------------------------------------------------
# self-test  (ground truth: a set that extends, a set that does not)
# ---------------------------------------------------------------------------
def _selftest():
    # --- primitive sanity ---
    assert collinear((0, 0), (1, 1), (2, 2)) and not collinear((0, 0), (1, 0), (0, 1))
    assert integer_distance((0, 0), (3, 4)) == 5
    assert integer_distance((0, 0), (1, 1)) is None       # sqrt(2) not integer
    assert rational_sqrt(Fraction(9, 4)) == Fraction(3, 2)
    assert rational_sqrt(Fraction(2)) is None
    # rectangle vertices are concyclic
    assert cocircular((0, 0), (4, 0), (4, 3), (0, 3))
    assert not cocircular((0, 0), (8, 0), (4, 3), (4, -3))  # the general-position quad

    # --- verify_input accepts the genuine general-position integral 4-set ---
    quad = [(0, 0), (8, 0), (4, 3), (4, -3)]
    ok, reason = verify_input(quad)
    assert ok, f"genuine quad rejected: {reason}"
    # ...and REJECTS a concyclic set (rectangle) and a collinear set
    assert not verify_input([(0, 0), (4, 0), (4, 3), (0, 3)])[0]     # concyclic
    assert not verify_input([(0, 0), (2, 0), (5, 0)])[0]             # collinear
    assert not verify_input([(1, 0), (0, 0), (8, 0), (4, 3)])[0]     # unnormalized P[0]

    # --- GROUND TRUTH 1: triangle T extends; search MUST re-find D=(4,-3) ---
    T = [(0, 0), (8, 0), (4, 3)]
    assert verify_input(T)[0]
    exts = extension_search(T, 8)
    assert (Fraction(4), Fraction(-3)) in exts, "search failed to re-discover the known extension"
    for X in exts:                                   # every reported extension independently re-verified
        assert verify_extension(T, X, 8)[0], f"search returned invalid extension {X}"
    # independent confirmation of the planted point on its own path
    assert verify_extension(T, (Fraction(4), Fraction(-3)), 8)[0]
    assert not verify_extension(T, (Fraction(1), Fraction(1)), 8)[0]   # non-integer distances

    # --- GROUND TRUTH 2: same triangle, tight bound B=4 -> no extension (triangle inequality) ---
    assert extension_search(T, 4) == [], "expected NO extension for T within bound B=4"

    print("SELFTEST PASS: heptagon_extend")


def run_spec(spec: dict) -> dict:
    if not isinstance(spec, dict) or not isinstance(spec.get("points"), list) or not isinstance(spec.get("bound"), int):
        return {"ok": False, "verdict": "REFUSED", "reason": "points list and integer bound are required"}
    if spec["bound"] < 1:
        return {"ok": False, "verdict": "REFUSED", "reason": "bound must be a positive integer"}
    try:
        result = analyze(str(spec.get("label", "unnamed")), spec["points"], spec["bound"])
    except (IndexError, TypeError, ValueError, ZeroDivisionError) as exc:
        return {"ok": False, "verdict": "REFUSED", "reason": f"invalid point-set spec: {exc}"}
    return {"ok": bool(result["verified_input"]),
            "verdict": "ANALYZED" if result["verified_input"] else "INVALID_INPUT",
            "result": result}


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Exact bounded rational integral point-set extension search")
    parser.add_argument("run", nargs="?")
    parser.add_argument("--input-file")
    parser.add_argument("--out")
    args = parser.parse_args()
    if args.input_file:
        try:
            result = run_spec(json.loads(Path(args.input_file).read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError) as exc:
            result = {"ok": False, "verdict": "REFUSED", "reason": f"cannot read input spec: {exc}"}
        text = json.dumps(result, indent=2)
        if args.out:
            Path(args.out).write_text(text + "\n", encoding="utf-8")
        print(text)
        raise SystemExit(0 if result.get("ok") else 2)
    if args.run == "run":
        results = {
            "blade": "heptagon_extend",
            "arithmetic": "exact rational (fractions.Fraction) — zero float decisions",
            "scope_note": (
                "Extension search runs inside Q^2 (rational coordinates). 'None up to B' means: no "
                "rational-coordinate point at integer distance <= B from all of P. Bounded, computational."
            ),
            "configs": [],
        }

        # genuine, hand-verified general-position integral configs
        triangle = [(0, 0), (8, 0), (4, 3)]                 # extends to the quad below
        quad = [(0, 0), (8, 0), (4, 3), (4, -3)]            # non-concyclic integral 4-set

        for name, P, B in [
            ("triangle_(0,0)(8,0)(4,3)", triangle, 20),
            ("quad_(0,0)(8,0)(4,3)(4,-3)", quad, 100),
        ]:
            row = analyze(name, P, B)
            results["configs"].append(row)
            print(f"{name}: verified_input={row['verified_input']}  "
                  f"B={B}  extensions_found={row.get('n_extensions_found')}  "
                  f"-> {row.get('verdict', row.get('extension_search'))}", flush=True)

        # honest growth attempt: try to grow the quad by one genuine, independently-verified point
        grow_B = 200
        current = list(quad)
        grown = None
        exts = extension_search(current, grow_B)
        if exts:
            cand = current + [exts[0]]
            if verify_input(cand)[0]:                        # independent re-check of the grown set
                grown = cand
        results["growth_from_quad"] = {
            "bound_B": grow_B,
            "grew": grown is not None,
            "grown_set": [_pt_str(p) for p in grown] if grown else None,
            "note": ("added an independently re-verified 5th point" if grown
                     else f"no rational integral 5th point found up to B={grow_B} (bounded result)"),
        }
        print(f"growth_from_quad (B={grow_B}): grew={grown is not None}", flush=True)

        # the real Kreisel-Kurz general-position integral HEPTAGON — DATA DEPENDENCY, NOT FABRICATED
        results["kk_heptagon"] = {
            "status": "PENDING REAL COORDINATES",
            "reason": (
                "Plugging the genuine Kreisel-Kurz minimal general-position integral heptagon into this "
                "engine requires its VERIFIED coordinates (normalized to P[0]=(0,0), P[1]=(a,0)). "
                "These are NOT fabricated here. Once supplied, run: verify_input(H) must pass, then "
                "extension_search(H, B) gives the bounded extension/maximality verdict for that heptagon."
            ),
        }
        print("kk_heptagon: PENDING REAL COORDINATES (not fabricated)", flush=True)

        out = "oracle/kbk/engine/heptagon_extend_results.json"
        with open(out, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nwrote {out}", flush=True)
    else:
        _selftest()
