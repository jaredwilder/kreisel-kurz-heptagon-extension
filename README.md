# Kreisel–Kurz heptagons and the integral-octagon extension problem

**An exact Diophantine reduction of the integral-octagon existence question**, with both Kreisel–Kurz heptagons reconstructed over `Q(sqrt 2002)` and the extension problem reduced to four simultaneous square conditions over a genus-zero base.

Author: Jared Wilder. First public timestamp: 2026-09-10.

Companion to [integral-point-sets](https://github.com/jaredwilder/integral-point-sets), which contains the theorem `d(2,8) > 30000`.

## Exact reconstruction of the two Kreisel–Kurz heptagons

Kreisel and Kurz (arXiv:0804.1303) determined `d(2,7)=22270`. Both of their heptagons—H1 of diameter 22270 and H2 of diameter 66810—lie in the number field **`Q(sqrt 2002)`**, where `2002=2·7·11·13`.

Both are reconstructed here with **all 21 pairwise distances integral and exactly matching the published matrices**, with no three points collinear and no four concyclic. H2's coordinates were recovered by solving the sign choices against all 21 distances; for example,

```text
P2 = (7690545/131, (91800/131) sqrt 2002).
```

Machine-readable certificates are in `certificates/kk_heptagon_certificates.json`. All arithmetic is exact.

## Exact extension system

With anchors `P0=(0,0)` and `P1=(a,0)`, an extension point `X=(x,y)` with integral distances `r=|XP0|` and `s=|XP1|` satisfies

```text
x = (a^2 + r^2 - s^2) / 2a,
y^2 = r^2 - x^2.
```

Writing `w=sqrt(2002)·y`, each remaining distance satisfies

```text
t_i^2 = A_i - 2 y_i w,
A_i = r^2 + rho_i - 2 x_i x.
```

Eliminating `w` gives a **master quartic surface**, quadratic in the squares `(r^2,s^2,t2^2)`. For H1 its leading coefficient is exactly

```text
123988225 = (a/2)^2.
```

Four linkage equations `L3` through `L6` are each **linear in the squares**, expressing `t_i^2` as an affine form in `(r^2,s^2,t2^2)`.

> **Exact reduction.** An extension exists if and only if there is an integer triple `(r,s,t2)` on the master quadric for which four explicit affine forms are simultaneously perfect squares.

The three-point subproblem—integral distance to `P0`, `P1`, and `P2`—is the discriminant equation

```text
W^2 = D(r^2,s^2),
```

a **smooth conic**, hence genus zero and rational. The base geometry is therefore rational; the arithmetic difficulty is concentrated in the four simultaneous square conditions.

A further field constraint is proved directly: in this embedding, every extension point lies in `Q(sqrt 2002)`.

## Arithmetic consequences and remaining question

The full distance system has no ordinary congruence obstruction, because the seven original vertices `X=P_j` are themselves exact integral solutions and therefore survive reduction modulo every prime. This was checked for all primes through 103 for both heptagons.

The natural next arithmetic object is therefore the **non-trivial quotient locus**, after removing those seven tautological points. Local methods applied to that quotient remain available.

Thus the contribution here is an exact reconstruction, field localization, and Diophantine reduction of the extension problem. The companion repository separately proves the diameter bound `d(2,8)>30000`.

## Additional certificates and programs

- `certificates/cayley_menger_results.json` — vanishing of the Cayley–Menger determinant on quadruples containing a fixed nondegenerate triangle forces embedding dimension at most 2, verified at `n=4,5,6` with Gram rank 2;
- `certificates/certificate_n*_d*.json` — six enumeration certificates for planar integral point sets under diameter constraints, including `n=5,d=73` and `n=6,d=174`, the exact minima;
- `certificates/closure_exact_H1.json` — **2,935,632 exact cells** enumerated, with a floating-point cross-check over 40,000 samples reporting maximum error `0.0`;
- `blades/` — historical directory containing the twenty-two exact-arithmetic programs used to produce the package.

## Scope

This repository reduces the existence question to a precise Diophantine system; it does not decide whether an integral octagon exists. Its theorem-level contribution is the exact reduction, reconstruction, field localization, and accompanying certified finite arithmetic.

## License

Apache-2.0.