# kreisel-kurz-heptagon-extension

**An exact Diophantine reduction of the integral-octagon existence question**, with both Kreisel–Kurz heptagons reconstructed over `Q(sqrt 2002)` and the extension problem reduced to four simultaneous square conditions over a genus-zero base.

Author: Jared Wilder. First public timestamp: 2026-09-10.

Companion to [integral-point-sets](https://github.com/jaredwilder/integral-point-sets), which carries the theorem `d(2,8) > 30000`.

## Exact reconstruction of the two Kreisel–Kurz heptagons

Kreisel and Kurz (arXiv:0804.1303) determined `d(2,7) = 22270`. Both their heptagons, H1 of diameter 22270 and H2 of diameter 66810, live over the number field **Q(sqrt 2002)**, where `2002 = 2·7·11·13`.

Both are reconstructed here with **all 21 distances integer and matching the published matrix exactly**, no three collinear, no four concyclic, certified. H2's coordinates were recovered by sign-solving against all 21 distances, for example `P2 = (7690545/131, (91800/131) sqrt 2002)`.

Machine-readable certificates are in `certificates/kk_heptagon_certificates.json`. All arithmetic is exact.

## Exact extension system

With anchors `P0 = (0,0)` and `P1 = (a,0)`, an extension point `X = (x,y)` with `r = |XP0|` and `s = |XP1|` integral satisfies

`x = (a^2 + r^2 - s^2) / 2a`,

`y^2 = r^2 - x^2`.

Writing `w = sqrt(2002) * y`, each remaining distance obeys the exact identity `t_i^2 = A_i - 2 y_i w` with `A_i = r^2 + rho_i - 2 x_i x`.

Eliminating `w` gives a **master quartic surface**, quadratic in the squares `(r^2, s^2, t2^2)`. For H1 its leading coefficient is `123988225 = (a/2)^2` exactly. The full coefficient list is included in the package.

Four **linkage equations** `L3` through `L6` are each **linear in the squares**, expressing `t_i^2` as an affine form in `(r^2, s^2, t2^2)`.

> **Exact reduction.** An extension exists if and only if there is an integer triple `(r,s,t2)` on the master quadric for which four explicit affine forms are simultaneously perfect squares.

The three-point subcondition—integer distance to `P0`, `P1` and `P2` alone—is the master discriminant `W^2 = D(r^2,s^2)`, a **smooth conic, hence genus zero and rational**. The base problem is therefore rational; the arithmetic difficulty is concentrated in the four stacked square conditions.

A further field constraint is proved directly: in the anchor embedding, any extension point lies in `Q(sqrt 2002)`.

## Arithmetic consequences and remaining surface

The full system has no plain congruence obstruction: the seven original vertices `X = P_j` are exact integral solutions of the full distance equations and therefore survive reduction modulo every prime. This was checked for all primes through 103 on both heptagons.

That observation localizes the next arithmetic target correctly: the **non-trivial quotient locus**, with the seven tautological points removed. Local methods applied to that quotient remain available; the computation above does not rule them out.

So the current mathematical status is an **exact reduction** plus exact reconstruction and field constraints. The companion repository separately proves the diameter bound `d(2,8) > 30000`.

## Also here

- `certificates/cayley_menger_results.json` — vanishing of the Cayley–Menger determinant on quadruples containing a fixed nondegenerate triangle forces embedding dimension at most 2, verified at `n = 4,5,6` with Gram rank 2;
- `certificates/certificate_n*_d*.json` — six enumeration certificates for planar integral point sets under diameter constraints, including `n = 5, d = 73` and `n = 6, d = 174`, the exact minima;
- `certificates/closure_exact_H1.json` — **2,935,632 exact cells** enumerated, with a float-error certificate over 40,000 samples reporting maximum error `0.0`;
- `blades/` — the twenty-two exact-arithmetic programs that produced the package.

## Scope

This repository reduces the existence question to a precise Diophantine system; it does not assert existence or nonexistence of an integral octagon. Its theorem-level contribution is the exact reduction, reconstruction, field localization, and accompanying certified finite arithmetic.

## License

Apache-2.0.