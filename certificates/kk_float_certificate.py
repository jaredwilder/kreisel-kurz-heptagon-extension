"""kk_float_certificate.py — measure the ACTUAL float64 error of the closure sweep's root computation.

The earlier certificate was VACUOUS: it only compared float roots against roots that happened to be
integers, and in 40000 random cells there were none, so it measured nothing (0 comparisons) and reported
a meaningless margin. This version compares the float64 root against the SAME root computed in 60-digit
decimal arithmetic from EXACT rational coefficients, for every sampled cell that has a real root --
regardless of whether the root is an integer. That is a real measurement.

What it certifies: the closure sweep keeps a cell if its float root lies within 1e-2 of an integer. A true
solution has an exactly-integer root, so it can only be missed if the float error exceeds 1e-2. The number
reported here is the worst observed float error over well-conditioned cells (|alpha| >= 1e-2, the regime the
float sweep is responsible for; smaller |alpha| is swept exactly by kk_closure_exact.py).

Run: python oracle/kbk/engine/kk_float_certificate.py [H1|H2] [samples]
"""
from __future__ import annotations
from fractions import Fraction as F
from decimal import Decimal, getcontext
import json
import sys
import numpy as np

sys.path.insert(0, "oracle/kbk/engine")
from kk_heptagon import H1, D1, D2, reconstruct_full
from kk_closure import best_triple

CHAR = 2002
getcontext().prec = 60


def _dec(fr: F) -> Decimal:
    return Decimal(fr.numerator) / Decimal(fr.denominator)


def certificate(points, D, label, samples=200000, seed=7, char=CHAR, alpha_floor=1e-2):
    n = len(points)
    cost, b, j, k, dbj, dbk = best_triple(D, n)
    xb, yb = points[b]
    uj, wj = points[j][0] - xb, points[j][1] - yb      # exact Fractions
    uk, wk = points[k][0] - xb, points[k][1] - yb
    D0 = uj * wk - wj * uk
    ujf, wjf, ukf, wkf, D0f = float(uj), float(wj), float(uk), float(wk), float(D0)
    dbj2, dbk2 = dbj * dbj, dbk * dbk

    rng = np.random.default_rng(seed)
    worst = 0.0
    worst_at = None
    compared = 0
    for _ in range(samples):
        cj = int(rng.integers(-dbj, dbj + 1)); ck = int(rng.integers(-dbk, dbk + 1))
        # ---- float64 path (exactly what the sweep does) ----
        ejf = float(dbj2 - cj * cj); ekf = float(dbk2 - ck * ck)
        Af = (wkf * ejf - wjf * ekf) / (2 * D0f); Bf = (wkf * cj - wjf * ck) / D0f
        Cf = (ujf * ekf - ukf * ejf) / (2 * char * D0f); Ddf = (ujf * ck - ukf * cj) / (char * D0f)
        af = Bf * Bf + char * Ddf * Ddf - 1.0
        if abs(af) < alpha_floor:
            continue                                  # handled by the exact sweep, not the float one
        bf = 2.0 * (Af * Bf + char * Cf * Ddf); gf = Af * Af + char * Cf * Cf
        discf = bf * bf - 4.0 * af * gf
        if discf < 0:
            continue
        sqf = discf ** 0.5
        froots = [(-bf + sqf) / (2 * af), (-bf - sqf) / (2 * af)]
        # ---- exact-coefficient, 60-digit path ----
        ej = F(dbj2 - cj * cj); ek = F(dbk2 - ck * ck)
        A = (wk * ej - wj * ek) / (2 * D0); B = (wk * F(cj) - wj * F(ck)) / D0
        C = (uj * ek - uk * ej) / (2 * char * D0); Dd = (uj * F(ck) - uk * F(cj)) / (char * D0)
        a = B * B + char * Dd * Dd - 1; bq = 2 * (A * B + char * C * Dd); g = A * A + char * C * C
        disc = bq * bq - 4 * a * g
        if disc < 0:
            continue
        sq = _dec(disc).sqrt()
        ad, bd = _dec(a), _dec(bq)
        droots = [(-bd + sq) / (2 * ad), (-bd - sq) / (2 * ad)]
        for dr in droots:
            err = min(abs(Decimal(repr(fr)) - dr) for fr in froots)
            compared += 1
            if float(err) > worst:
                worst = float(err); worst_at = (cj, ck, str(dr)[:40])
    return {"label": label, "samples": samples, "roots_compared": compared,
            "alpha_floor": alpha_floor, "max_abs_float_error": worst, "worst_at": worst_at,
            "sweep_window": 1e-2,
            "margin": (1e-2 / worst) if worst > 0 else None}


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "H1"
    samples = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
    pts, D = (H1, D1) if which == "H1" else (reconstruct_full(D2, 2002), D2)
    cert = certificate(pts, D, which, samples)
    print(json.dumps(cert, indent=1))
    if cert["roots_compared"] == 0:
        print("!! VACUOUS CERTIFICATE -- measured nothing. Do not use.")
    else:
        print(f"\nworst float64 error over {cert['roots_compared']:,} real roots = {cert['max_abs_float_error']:.3e}")
        print(f"sweep keeps anything within 1e-2 of an integer -> margin {cert['margin']:.3e}x")
    with open(f"oracle/kbk/engine/float_certificate_{which}.json", "w") as f:
        json.dump(cert, f, indent=1)
