# kreisel-kurz-heptagon-extension

An exact Diophantine reduction of the integral-octagon existence question, with both Kreisel-Kurz
heptagons reconstructed over Q(sqrt 2002) and the extension problem reduced to four stacked
square-conditions over a genus-zero base.

Author: Jared Wilder. First public timestamp: 2026-09-10.

Companion to [integral-point-sets](https://github.com/jaredwilder/integral-point-sets), which
carries the theorem d(2,8) > 30000.

## Phase 1: the heptagons, reconstructed exactly

Kreisel and Kurz (arXiv:0804.1303) determined d(2,7) = 22270. Both their heptagons, H1 of diameter
22270 and H2 of diameter 66810, live over the number field **Q(sqrt 2002)**, where
2002 = 2 * 7 * 11 * 13.

Both are reconstructed here with **all 21 distances integer and matching the published matrix
exactly**, no three collinear, no four concyclic, certified. H2's coordinates were recovered by
sign-solving against all 21 distances, for example P2 = (7690545/131, (91800/131) sqrt 2002).

Machine-readable certificates in `certificates/kk_heptagon_certificates.json`. All arithmetic exact.

## Phase 2: the extension system

With anchors P0 = (0,0) and P1 = (a,0), an extension point X = (x,y) with r = |XP0| and s = |XP1|
integral satisfies x = (a^2 + r^2 - s^2) / 2a and y^2 = r^2 - x^2. Writing w = sqrt(2002) * y, each
distance obeys the exact identity t_i^2 = A_i - 2 y_i w with A_i = r^2 + rho_i - 2 x_i x.

Eliminating w gives a **master quartic surface**, quadratic in the squares (r^2, s^2, t2^2). For H1
its leading coefficient is 123988225 = (a/2)^2 exactly. The full coefficient list is in the package.

Four **linkage equations** L3 through L6 are each **linear in the squares**, expressing t_i^2 as an
affine form in (r^2, s^2, t2^2).

**The reduction:** an extension exists if and only if there is an integer triple (r, s, t2) on the
master quadric with four affine forms simultaneously perfect squares.

The three-point sub-condition, integer distance to P0, P1 and P2 alone, is the master discriminant
W^2 = D(r^2, s^2), which is a **smooth conic, hence genus zero and rational**. So the base of the
problem is easy and **the entire difficulty sits in the four stacked square-conditions.**

## Phase 3: what the arithmetic attack established, and exactly what it did not

Verdict: **REDUCED.** Not closed.

**No congruence obstruction exists on the full system.** The seven points X = P_j are themselves
exact integral solutions of the full system, so they survive reduction modulo every prime. Verified
for all primes up to 103, on both heptagons.

The package then draws its own boundary, and this is quoted rather than paraphrased:

> This is the exact statement, and no more: it rules out a congruence/reduction obstruction *on the
> full system*. It does NOT claim "every local method fails" - a local method applied to the
> QUOTIENT variety (trivial points removed) is untouched by this and remains open. The real object
> is the *non-trivial* rational locus.

A field constraint is also **proven rather than invoked**: in the anchor embedding, any extension
point must lie in Q(sqrt 2002).

## Also here

- `certificates/cayley_menger_results.json` — vanishing of the Cayley-Menger determinant on
  quadruples containing a fixed nondegenerate triangle forces embedding dimension at most 2,
  verified at n = 4, 5, 6 with Gram rank 2.
- `certificates/certificate_n*_d*.json` — six enumeration certificates for n-point planar sets
  under diameter constraints, including n = 5 at d = 73 and n = 6 at d = 174, the exact minima.
- `certificates/closure_exact_H1.json` — 2,935,632 exact cells enumerated, with a float-error
  certificate over 40,000 samples reporting maximum error 0.0.
- `blades/` — the twenty-two exact-arithmetic programs that produced all of the above.

## Scope

**Whether an integral octagon in general position exists is open.** This reduces the question; it
does not answer it. The companion repository proves only that any such octagon has diameter greater
than 30000.

## License

Apache-2.0.
