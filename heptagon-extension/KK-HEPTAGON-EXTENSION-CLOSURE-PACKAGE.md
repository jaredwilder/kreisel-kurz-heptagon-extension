# Kreisel–Kurz Heptagon Extension Closure Package — receipt

Source objects: Kreisel & Kurz, arXiv:0804.1303, `d(2,7)=22270`. Two heptagons H1 (diam 22270) and
H2 (diam 66810), both over the number field **Q(√2002)** (characteristic 2002 = 2·7·11·13).
All computation exact ($0, deterministic, zero-LLM). Blades: `kk_heptagon.py`, `kk_extension_equations.py`,
`kk_extension_sieve.py`, `kk_extension_search.py`.

## Phase 1 — reconstruction certificates (COMPLETE)

| heptagon | field | 21 distances integer & match published matrix | no 3 collinear | no 4 concyclic | CERTIFIED |
|---|---|---|---|---|---|
| H1 | Q(√2002) | ✓ (exact, Figure-1 coords) | ✓ | ✓ | **TRUE** |
| H2 | Q(√2002) | ✓ (coords reconstructed from matrix (2)) | ✓ | ✓ | **TRUE** |

Machine-readable: `kk_heptagon_certificates.json`. H2 exact coordinates recovered (sign-solved to match all
21 distances), e.g. P2 = (7690545/131, 91800/131·√2002).

## Phase 2 — exact extension Diophantine system (COMPLETE)

Anchors P0=(0,0), P1=(a,0). Extension X=(x,y), r=|XP0|, s=|XP1| ∈ Z, x=(a²+r²−s²)/(2a), y²=r²−x².
With w := √2002·y (so w²=2002(r²−x²)) and points P_i=(x_i, y_i√2002), rho_i=|P_i|²=D[0][i]²:
each distance obeys the **identity** t_i² = A_i − 2 y_i w, A_i = r²+rho_i−2x_i x (validated exactly).
Eliminating w:

- **MASTER surface** M(r,s,t2)=0 — a quartic; quadratic in the squares (r²,s²,t2²). For H1:
  `115433536 r⁴ − 113525712 r²s² − 117341360 r²t2² − 60318681457512960 r² + 122080401 s⁴`
  `− 130635090 s²t2² − 57300321130741440 s² + 123988225 t2⁴ − 56303406090964800 t2² + 27956215004362898680857600 = 0`
  (leading coeff 123988225 = (a/2)² exactly.) H2 analogue printed by the blade.
- **LINKAGES** L_3..L_6 — each LINEAR in the squares, giving t_i² as an affine form ℓ_i(r²,s²,t2²). For H1:
  `L3: 133104 r² + 387096 s² + 411864 t2² − 932064 t3² − 135117693635040 = 0`  (and L4,L5,L6 similarly).

**Reduction:** an extension ⇔ a triple (r,s,t2)∈Z³ on the master quadric with ℓ_3,ℓ_4,ℓ_5,ℓ_6 simultaneously
perfect squares. The **3-point sub-condition** (integer distance to P0,P1,P2) is the master discriminant
`W² = D(r²,s²)`, a **smooth conic ⇒ genus 0 (rational)** — so the difficulty is entirely the four stacked
square-conditions, not the base.

## Phase 3 — arithmetic attack (status: local methods provably insufficient; verdict REDUCED)

1. **No congruence obstruction on the full system (precise claim).** The 7 points X=P_j (j=0..6) are exact
   integral solutions of the full system, so they reduce to survivors mod every prime — hence **no congruence
   obstruction on the full distance system can be empty** (verified survivors p ≤ 103, both heptagons). This is
   the exact statement, and no more: it rules out a congruence/reduction obstruction *on the full system*. It
   does NOT claim "every local method fails" — a local method applied to the QUOTIENT variety (trivial points
   removed) is untouched by this and remains open. The real object is the *non-trivial* rational locus.
2. **Field constraint — PROVEN (not merely invoked).** In the anchor embedding P0=(0,0), P1=(a,0): any point
   P at integral distance from both anchors has x=(a²+r²−s²)/(2a) ∈ Q and y²=r²−x² ∈ Q. The triangle
   △P0P1P has 16·Area² = 4a²y², so its characteristic (squarefree part) = squarefree(y²). By the
   characteristic theorem this equals the set's characteristic 2002. Hence y² = 2002·(rational square), i.e.
   **y = q√2002 with q ∈ Q, so X ∈ Q(√2002).** (Verified: all 7 heptagon points have squarefree(y²)=2002.)
   The subtlety flagged in review is real — the characteristic theorem alone does NOT give a single field;
   the anchor placement is what closes it, and the two together give the 4-line proof above. Equivalent gate:
   `2002 · Heron16(r,s,a) = perfect square`.
3. **Triangle-inequality band.** Since a is the diameter, any X satisfies r+s>a: the search band is
   r,s ∈ (a−r, a+r), i.e. distances comparable to the diameter — NOT small heights.
4. **Bounded receipt (subordinate, COMPLETE for BOTH heptagons, updated 2026-07-25).** Exact search
   (`kk_deep_search.py`, regression-verified: correctly reproduces + suppresses the known trivial X=P_j
   solutions), band r,s ∈ [1, B]:

   | heptagon | B | B / diameter | gate-passes | genuine (non-vertex) extensions | wall time |
   |---|---|---|---|---|---|
   | H1 (d=22270) | 100000 | 4.49× | 2088 | **0** | 1117s |
   | H2 (d=66810) | 100000 | 1.50× | 1945 | **0** | 1972s |

   Zero genuine octagons for either heptagon within these bands (supersedes the earlier partial B=44540/1236-
   gate-pass H1-only run). This does NOT prove MAXIMAL — X could lie farther out (no proven height bound; see
   §3.1) — only curve-descent or an effective height bound can close that. It IS a real, honest, exhaustive
   floor in the correct universe (quadratic-field, not lattice), now real for H2 as well as H1.

   **Denominator axis (q ≥ 2), per the §4c correction** — `kk_scaled_search.py`, integer extensions of q·H1
   (⟺ denominator-q RATIONAL extensions of H1):

   | target | scaled diameter | B | gate-passes | genuine extensions | wall time |
   |---|---|---|---|---|---|
   | 2·H1 | 44540 | 100000 | 2086 | **0** | 2025s |
   | 3·H1 | 66810 | 120000 | 2379 | **0** | 2953s |

   So H1 admits **no rational extension of denominator 1, 2, or 3** in these bands (in original units: distances
   up to 4.49×, 2.25×, and 1.80× the diameter respectively). Three denominator slices now exhaustively empty,
   plus H2 at q=1. Still height-bounded, still not a maximality proof — but the denominator axis, which was
   entirely unexplored before this correction, is now covered for q ≤ 3.

## Phase 4b — LOCAL DEFORMATION THEORY of the seven known solutions (2026-07-25)

Prompted by the critique "the seven known points are not trivial, they are the entire mystery." Computed
exactly (`kk_local_rigidity.py`, sympy over ℚ), on V = {F_i := 2002 t_i² − Q_i(x,w) = 0}_{i=0..6} ⊂ 𝔸⁹ with
Q_i(x,w) = 2002(x−x_i)² + (w−2002y_i)²:

1. **The seven known solutions are exactly the nodes of the extension surface.** At X=P_j the Jacobian has
   an identically-zero row (row j) → **rank 6, not 7** (verified all j=0..6); tangent-space dim 3 > 2 = dim V.
   The t-block of the Jacobian is diag(2·2002·t_i), so rank drops **iff some t_i = 0 iff Q_i(x,w)=0**, and Q_i
   is a sum of two squares, so over ℝ that happens **iff X = P_i**. So Sing(V)(ℝ) = the 7 trivial points,
   exactly. They are not noise in the search — they are the singular locus.
2. **Exact tangent cone at each node:** in local coords (u,v,t_j) = (x−x_j, w−2002y_j, t_j), F_j has no linear
   part and equals **2002 t_j² − 2002u² − v² = 0** — an A₁ quadric cone. (The other six F_i have nonzero linear
   parts and solve for t_i by the implicit function theorem, leaving exactly this cone.)
3. **NOT first-order rigid.** That cone has rational points, e.g. (u,v,t_j) = (−2001/2, 2002, 2003/2) (residual
   exactly 0). So each node is **split over ℚ**: rational tangent directions to V exist at every known
   solution. **The obstruction to an 8th point is therefore neither first-order nor local at the nodes** — it
   is higher-order/global. This matches the independent earlier finding that no congruence obstruction on the
   full system can be empty.

## Phase 4c — THE SEARCH SPACE WAS TOO NARROW (correction, 2026-07-25)

General position and "all distances rational" are both **scale-invariant**, and distances scale linearly, so
multiplying an 8-point rational-distance set by the common denominator of its 28 distances gives an integral
one. Hence:

> **∃ integral general-position 8-set ⟺ ∃ rational-distance general-position 8-set.**

A rational extension of H1 with denominator q is exactly an **integer** extension of the scaled heptagon q·H1,
whose characteristic is still 2002 (verified for q=1,2,3,5,7: 16·Area² scales by q⁴, so the squarefree part is
unchanged). **Therefore all the earlier integer-(r,s) searches — including the B=100000 runs — were the q=1
slice only**, and could not have seen a denominator-q rational extension. This also widens the target from
"contains H1 exactly" to "contains any scaled copy of H1", which partly answers the earlier adversarial note
that an 8-set's 7-subsets need not be H1 or H2 themselves. New driver: `kk_scaled_search.py` (regression-
verified: correctly locates and suppresses the trivial points of 2·H1).

## Phase 4d — ⭐ THE CHORD-FORCING INVARIANT (2026-07-25) — the session's real result

Prompted by the critique "hunt for a hidden invariant that distinguishes genuine vertices from every
rejected extension." Found one, and it is **measured, quantified, and yields an unconditional proof strategy**.

**Step 1 — a genuine near-miss exists.** `kk_linkage_stats.py` histogrammed, over every gate-passing X in
r ≤ 40000 for H1, how many of the 5 non-anchor linkage square-conditions hold (t₀=r, t₁=s are integral by
construction). Result: 0/5 → 2174, 1/5 → 56, 2/5 → 16, **4/5 → 1**, 5/5 → 5. The five 5/5's are exactly the
five known non-anchor trivial embeddings (instrument validated). The **4/5 is a genuine non-trivial near-miss**:

> **X = (x, w) = (112136/17, 12588576/17)**, i.e. y-coefficient 6288/17 = **2 × y(P₃) = 2 × y(P₅)**.
> Integral distances to **six of seven** vertices: |XP₀|=17816, |XP₁|=22794, |XP₃|=11397, |XP₄|=11832,
> |XP₅|=8908, |XP₆|=17748. Only |XP₂| fails: t₂² = 31113124, which is **960 short of 5578²**.

So {P₀,P₁,P₃,P₄,P₅,P₆,X} has **all 21 distances integral** with diameter 22794 — *below* Kreisel–Kurz's proven
uniqueness bound of 30000. It must therefore fail general position, and it does, informatively:
**three collinear triples (P₀P₅X, P₁P₃X, P₄P₆X) and three concyclic quadruples.** X is the concurrence point
of three chords of the heptagon. (No contradiction with K–K. Verified exactly.)

**Step 2 — the mechanism is real, not anecdotal.** `kk_chord_hypothesis.py` measured chord-membership for
every non-trivial gate-passing candidate in the band:

| population | on ≥1 chord |
|---|---|
| all non-trivial gate-passing candidates (base rate) | 28 / 2247 = **1.246 %** |
| candidates with ≥ 2 of 5 linkage conditions | 14 / 17 = **82.35 %** |
| the single 4/5 candidate | on **3** chords (expected 0.0125) |

**Enrichment 66.1×, binomial p = 1.4 × 10⁻²⁴.** Arithmetic richness and chord-membership are strongly, in fact
decisively, correlated — and chord-membership is *exactly what general position forbids*. The integrality
conditions like to be simultaneously satisfied precisely where the geometry excludes the point.

**Step 3 — the proof strategy this buys (no Bombieri–Lang).**

> **Chord-Forcing Lemma (CONJECTURE, unproven).** Any X ∈ ℚ(√2002) at integral distance from all seven
> vertices of H1 lies on a chord P_iP_j.
>
> **Corollary (unconditional if the Lemma holds).** H1 admits no general-position integral extension.
> *Proof.* Let X extend H1 in general position; then X is at integral distance from all 7 vertices, so by the
> Lemma X, P_i, P_j are collinear for some i≠j — contradicting no-3-collinear. ∎

This is consistent with the trivial solutions (X = P_j lies on every chord through P_j), and it would give
MAXIMAL **unconditionally** — replacing the Bombieri–Lang-conditional route entirely.

**Step 4 — THE MECHANISM, PROVEN (why chord points are arithmetically rich).**

> **Freebie Lemma (PROVEN — trivial, but load-bearing).** If X is collinear with vertices P_i, P_j and
> d_ij = |P_iP_j| ∈ ℤ, then |XP_j| = | |XP_i| ± d_ij |. Hence **|XP_i| ∈ ℤ ⟺ |XP_j| ∈ ℤ**: each chord
> through X converts **two** integrality conditions into **one**. Verified on all three chords through the
> near-miss (P₀P₅: 17816 = 8908 + 8908; P₁P₃: 22794 = 11397 + 11397; P₄P₆: 17748 = 11832 + 5916).

This *quantitatively* explains the enrichment. With the measured per-condition square rate p ≈ 0.0104, a generic
X needs 5 independent conditions (p⁵ ≈ 1.2×10⁻¹⁰), while an X on 3 chords covering 6 vertices needs only ~2–3
(p² … p³ ≈ 10⁻⁴ … 10⁻⁶) — a **10³–10⁵× boost**, more than enough to produce the observed 66×.

**And the near-miss is generated by a hidden affine relation inside H1 itself:**

> **P₀ − P₁ = 2(P₅ − P₃)** — verified exactly. Segment P₃P₅ is parallel to P₁P₀ with exactly **half** its
> length (|P₃P₅| = 11135 = 22270/2 = |P₀P₁|/2; both horizontal).

That single relation makes **X := 2P₅ − P₀ = 2P₃ − P₁ = 3P₄ − 2P₆** one well-defined point (all three
identities verified exactly) lying on three chords, with **six of seven distances forced integral** by
reflection (P₅ is the midpoint of P₀X, P₃ the midpoint of P₁X). Only |XP₂| is unconstrained by any relation —
and that is exactly the one that fails.

**This answers "why exactly these seven?"** The heptagon's own internal affine coincidences — the same
arithmetic that lets it be integral at all — generate arithmetically-rich extension candidates, and those
candidates are *automatically collinear with vertices*, hence excluded by general position. Arithmetic richness
and admissibility are structurally, not accidentally, opposed.

**Honest status of the Lemma.** The *mechanism* (Freebie Lemma + affine generation) is proven and explains the
enrichment. The **converse — that a FULL 7-of-7 extension must be a chord point — remains unproven**, and that
converse is the whole load-bearing gap: the Freebie Lemma shows chord points are hugely favoured, but does not
preclude an off-chord solution. Supporting evidence is the 66× enrichment plus the 4/5 point sitting
on three chords. Against over-reading: 3 of the 17 high-linkage candidates are *not* on a chord (all at the
2/5 level), so the Lemma cannot hold at a low threshold; and there is no non-trivial 5/5 candidate in existence
to test directly (that is the open problem). One band, one heptagon. It is a well-posed, falsifiable target —
which is strictly better than the conjectural surface-theoretic ceiling it would replace.

## ⛳ PHASE 5 — **CLOSED. H1 IS MAXIMAL.** (2026-07-25, machine-verified, unconditional)

> **THEOREM (exhaustively verified).** Let H be either Kreisel–Kurz heptagon (H1, diameter 22270; H2,
> diameter 66810). There is **no point X in the plane, other than the seven vertices themselves, at integral
> distance from all seven vertices of H.** Hence neither heptagon admits any integral extension — in
> particular no octagon. **BOTH KNOWN KREISEL–KURZ HEPTAGONS ARE MAXIMAL.**

**No height bound. No Bombieri–Lang. No Chabauty. No floating point.** The whole surface-theoretic apparatus
turned out to be the wrong altitude; an elementary finite parametrisation closes it.

**Why the enumeration is complete (the key move).** Pick a base vertex P_b and two others P_j, P_k
(non-collinear). For *any* candidate X with all distances integral,
  c_j := |XP_b| − |XP_j| and c_k := |XP_b| − |XP_k|
are **integers**, and by the triangle inequality |c_j| ≤ d_bj, |c_k| ≤ d_bk. **That range is finite and does
not depend on how far away X lies** — which is exactly what every earlier bounded search lacked. Given
(c_j, c_k), three distances from three non-collinear points determine X uniquely, and r = |XP_b| satisfies an
explicit quadratic. Clearing denominators gives the all-integer form

    4·ALPHA·r² + 4·BETA·r + GAMMA = 0,   r = (−BETA ± √(BETA²−ALPHA·GAMMA)) / (2·ALPHA)

with ALPHA, BETA, GAMMA explicit bigints. Every test is exact: `≥ 0`, `is a perfect square`, `divides`.

**Results — two independent exhaustive sweeps, different base triples, different grids:**

| heptagon | base triple | grid cells (all exact) | solutions | vertices | **non-trivial** | gate | time |
|---|---|---|---|---|---|---|---|
| H1 (d=22270) | (P₄; P₅, P₆) | 136,801,313 | 6 | 6 | **0** | 6/6 PASS | 153s |
| H1 (d=22270) | (P₅; P₄, P₀) | 205,982,337 | 6 | 6 | **0** | 6/6 PASS | 223s |
| **H2 (d=66810)** | (P₄; P₂, P₅) | **1,335,425,575** | 6 | 6 | **0** | 6/6 PASS | 1530s |

**Total: 1.678 billion cells swept in exact integer arithmetic. Zero non-trivial solutions.**

The **correctness gate** is mandatory and non-negotiable: the sweep must recover all six known non-base
vertices, or the run is declared unsound and no claim is made from it. Both runs pass 6/6.

**Two bugs were caught and fixed on the way — both would have produced a false "closure":**
1. *Tangency/float*: for X = P_j the discriminant is exactly 0, but float64 returned −2.4×10⁻⁹ via
   catastrophic cancellation, silently dropping two known vertices. Caught by the correctness gate.
2. *Precision margin*: a proper 60-digit certificate showed the float sweep's worst root error was
   3.06×10⁻³ against a 10⁻² acceptance window — a margin of only **3.3×**, far too thin. (An earlier
   "certificate" reporting a 10¹⁶× margin was **vacuous** — it compared 0 roots.) Resolved by discarding
   floats entirely: exact integer arithmetic costs only ~1.8 µs/cell, so the whole grid is swept exactly.

**Corollary — the Chord-Forcing Lemma is now PROVEN for H1** (vacuously): the only X at integral distance
from all seven vertices are the vertices, and each vertex lies on chords through it. The 66× enrichment and
the Freebie Lemma explained *why* near-misses cluster on chords; the exhaustive sweep shows no candidate ever
completes.

### Scope — stated precisely, no overreach
This closes the **extension question for the known heptagons**. It does **not** settle Erdős's n = 8 problem:
an 8-point integral general-position set need not contain H1 or H2 — its seven-point subsets could be entirely
different heptagons (see the adversarial-review note). What is now settled is that **these** heptagons are
dead ends, so any octagon must come from a genuinely different heptagon.

## Verdict (Phase 5)

**H1, H2: REDUCED.** Exact object: the master quartic surface M(r,s,t2)=0 with the four linear linkages, over
Q(√2002). **The one exact unresolved arithmetic statement whose proof settles maximality:**

> Does { M(r,s,t2)=0 ; ℓ_i(r²,s²,t2²) a perfect square, i=3..6 } admit a solution (r,s,t2)∈Z³ with the
> induced X ∉ {P_0,…,P_6} and general position — i.e. a **non-trivial rational point** beyond the 7 known?

**The path to a maximality PROOF is a chain, and only its first link is currently solid — do not collapse it:**

1. ✅ construct the exact extension equations *(done, validated as identities)*
2. ⬜ identify the correct variety / its irreducible components *(in progress; it is a SURFACE, dim 2 — see §4)*
3. ⬜ prove a birational model
4. ⬜ compute genus (of the fibration's curve fibers) — **with a CAS, not by eye**
5. ⬜ study the rational points (Mordell-Weil / Chabauty–Coleman)
6. ⬜ deduce maximality

Links 2–6 are goals, not results. "MAXIMAL ⇔ curve has only trivial rational points" is the *destination* of
this chain, not an established equivalence today.

## Phase 4 — the extension variety (algebraic geometry; the real object)

Corrected framing after review. The extension locus is **not a curve — it is a surface**:

- **Dimension = 2, proven.** In coordinates X=(x,w) (w := √2002·y, both rational — see the PROVEN field
  containment above), each point P_i gives `E_i: 2002(x−x_i)² + (w−2002 y_i)² = 2002 t_i²`. Each rational-
  distance condition adds a variable (t_i) AND an equation → net 0. So V ⊂ A⁹ is 2-dimensional. Confirmed by
  Lang–Weil point count: |V(𝔽₁₀₀₉)| = 1,033,600 ≈ p² (ratio 1.015; a curve would give ≈ p).
- **Structure:** V is the fiber product of **7 conic-branched double covers** of the (x,w)-plane (branch conics
  D_i = 0). The 7 trivial embeddings X=P_j are rational points on it.
- **Fibration & fiber genus (Singular-COMPUTED, not eyeballed).** Project V to the x-line; the fiber over a
  generic x₀ is a curve whose geometric genus (over 𝔽₃₂₀₀₃, char-independent) is:

  | # points | 2 | 3 | 4 | 5 |
  |---|---|---|---|---|
  | fiber genus | 1 | 5 | 17 | 49 |

  (A naive fiber-product formula predicted 1,3,9,25 and was WRONG — Singular gives 1,5,17,49. Compute, never
  guess.) The 2-point fiber is elliptic (genus 1); by 3 points the fiber is genus 5 → **Faltings ⇒ each fiber
  has finitely many rational points.** The four COMPUTED values fit the recurrence gₖ = 2gₖ₋₁ + 2ᵏ⁻¹ − 1,
  which *predicts* g₆=129, g₇=321 (full-heptagon fiber) — but that is a fit, NOT a computation: the genus-~130
  space-curve normalization exceeded Singular's reach here, so 129/321 remain unconfirmed and are not banked.

### Kodaira dimension — COMPUTED (Sage 9.5 stood up in WSL): κ(V) = 2, GENERAL TYPE

Key structural correction (caught by computing, not assuming): the branch conics D_i = 2002(x−x_i)²+(w−2002y_i)²
are **degenerate** — a sum of two squares factors over ℂ into two conjugate lines, and every "circle" passes
through the SAME two circular points at infinity [1:±i√2002:0]. So the naive "7 smooth conics ⇒ 4H ⇒ general
type" was right on the class but wrong on the geometry. In isotropic coordinates p=√2002·x+iw, q=√2002·x−iw the
branch locus becomes a **grid**: 7 vertical lines {p=a_i} + 7 horizontal {q=b_i}, with a_i=√2002 x_i+2002 i y_i
(7 distinct, verified) and t_i² = (p−a_i)(q−b_i).

Hence **V is birational to (C_p × C_q)/(ℤ/2)⁷** (diagonal action), where C_p, C_q are the (ℤ/2)⁷ covers of the
p- and q-lines. Genera:
- Sage: the hyperelliptic quotient y²=∏(p−a_i) (7 distinct roots) has **genus 3** ⇒ g(C_p) ≥ 3.
- Riemann–Hurwitz (verified vs Sage at n=2→0, n=3→1): **g(C_p)=g(C_q)=129**.

C_p × C_q = product of two genus-129 curves ⇒ κ = 1+1 = 2. V is a finite (ℤ/2)⁷ quotient ⇒ **κ(V)=2** (Kodaira
dimension is invariant under finite quotients; the fixed-point singularities are canonical). Blades:
`kk_kodaira.py`-equivalent scripts + Sage `/tmp/kod3.sage`.

### Consequence for maximality — stated honestly (conditional vs unconditional)

- **Unconditional (Faltings):** V is fibered over the x-line with fibers of genus ≥ 5 (Singular: 1,5,17,49,…),
  so each fiber has finitely many rational points ⇒ for any bounded x-window only finitely many extensions
  (the B8 sweep is one instance).
- **Conditional (Bombieri–Lang, OPEN):** a surface of general type has rational points NOT Zariski-dense — they
  lie on a finite union of curves + finitely many points. With the 7 trivial points known, this says only
  finitely many sporadic extensions can exist, confined to special curves. This is strong structural support
  for MAXIMALITY, but Bombieri–Lang is a major open conjecture — it is NOT a proof that H1/H2 are maximal.

**Honest tool state now:** Singular (genus/decomposition), gp (elliptic/hyperelliptic rank), **Sage 9.5 (WSL)**.

### The organizing principle — general theorem (heptagon-independent), 2026-07-24

The κ=2 computation never used any special property of H1 — only "7 distinct points over ℚ(√c)". Verified on
H1, H2, and a synthetic heptagon over ℚ(√6): identical grid structure, 7 distinct a_i. This upgrades the
result from "H1's extension is general type" to a **classification theorem**:

> **Theorem (extension-variety type).** For n distinct points over ℚ(√c), the rational-distance extension
> variety is birational to (C_p × C_q)/(ℤ/2)ⁿ, where C_p, C_q are the (ℤ/2)ⁿ covers of the p-, q-lines with
> genus **g_n = 2ⁿ⁻²(n−3) + 1** (Riemann–Hurwitz; small cases n=2→0, n=3→1 independently checked; n=7→129).
> Hence the Kodaira dimension is:
> - **n = 2:** rational (κ = −∞)
> - **n = 3:** κ = 0 — **K3 / abelian type** (this is exactly the K3 case anticipated in math724)
> - **n ≥ 4:** **general type (κ = 2)** — since g_n ≥ 5 ≥ 2.

For heptagons (n=7) the general-type conclusion is Sage-backed (genus-3 hyperelliptic quotient, valid for any
7 distinct a_i). For n=4 it rests on the R-H genus g₄=5 (no hyperelliptic quotient reaches genus 2 there).

**Consequence for "must every 8-set contain H1 or H2?" — the question becomes moot.** Answer: **no** (an
8-set's seven 7-subsets can be other heptagons, and none is known to exist while n=8 is open). But we do NOT
need it: *every* integral heptagon — whatever appears as a 7-subset — has a **general-type** extension
variety by the theorem above. So the structural analysis already covers all possible 7-subsets uniformly,
without identifying them. This is the search-space collapse: not "search millions of 8th points," but "every
extension locus is general type ⇒ (Bombieri–Lang) finitely many rational points, uniformly."

**Honest scope (unchanged):** general type ⇒ maximality only *conditional on Bombieri–Lang*; unconditionally it
gives Faltings-finiteness on each genus-≥2 fiber but not global non-existence. The theorem strengthens the
organizing principle; it does not by itself close Erdős n=8.

### Chabauty–Coleman assessment — TWO walls, one structural, one tooling (2026-07-24)

Attempted per the plan; here is the honest result.

1. **There is no genus-2 fiber.** The fibration over the x-line has fibers of genus 1,5,17,49,… — 2 never
   occurs. The genus-2 object is a *quotient*: over a fixed rational x₀, the 3-point sub-condition {P₀,P₂,P₃}
   gives the quadratic-twist curve `C_{x₀}: 2002·Y² = g₀(w)g₂(w)g₃(w)` (g_i = 2002(x₀−x_i)²+(w−2002y_i)²),
   degree 6, **genus 2 over ℚ (Sage-confirmed)**, whose ℚ-points ↔ rational-distance points to P₀,P₂,P₃ on
   that fiber. Reduced model computed (Sage `C.reduce()`).
2. **Structural wall (fatal regardless of tooling): Chabauty is per-curve; this is a surface.** C_{x₀} depends
   on the free parameter x₀, and there is **no height bound on x₀** (that claim was already retracted). So even
   a complete Chabauty run on one C_{x₀} only clears extensions with that single x-coordinate — an infinite
   family of fibers remains. Chabauty–Coleman on the fibers therefore **cannot** yield unconditional global
   maximality.
3. **Tooling wall:** Sage 9.5 stock has **no** genus-2 Jacobian rank / 2-descent / analytic-rank / Chabauty
   (probed: all absent), and **Magma is not installed** (only the interface stub; `magma -n` → not found).
   gp does elliptic/hyperelliptic genus but not Chabauty–Coleman. So the method is not executable here even for
   a single fiber.

**Consequence / honest redirect.** The real bottleneck for *unconditional* maximality is NOT the Chabauty step
— it is an **effective height bound on the 8th point** (a Vojta-type statement; open). WITH such a bound, the
surface reduces to finitely many fibers and the existing exhaustive fiber-search (`kk_extension_search.py`)
closes maximality directly — no Chabauty needed. WITHOUT it, no per-curve method (Chabauty included) suffices.
This is precisely the wall Kreisel–Kurz hit (their Ω(d³) search reached only diameter 70000). So the frontier
is the height bound, and the strongest current statement remains: **κ=2 (general type) ⇒ maximality holds
conditional on Bombieri–Lang.**

## Next compute (specified, not hand-waved)
- Extract the explicit genus of the 4-, 5-, 6-, 7-point curves (resultant tower on the master surface).
- Run the B8 band sweep to completion (modular-prefiltered) → largest diameter exhaustively excluded.
- Where a curve is elliptic/hyperelliptic: Mordell–Weil rank + Chabauty–Coleman on the non-trivial locus.
