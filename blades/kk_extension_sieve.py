"""kk_extension_sieve.py — Phase 3 (congruence sieve + local solubility) for the KK heptagon extension.

Master M(r,s,t2)=0 and linkages L_i (t_i^2 linear in r^2,s^2,t2^2). We sieve: for prime p (not dividing the
relevant leading coeffs), enumerate (r,s,t2) in (Z/p)^3; keep those with M==0; for each, test whether every
linkage admits a t_i (i.e. the required value is a QR mod p). A modulus with ZERO survivors would be a rigorous
MAXIMAL certificate. We ALSO inject the 7 trivial embeddings X=P_j to test whether they survive every modulus
(if so, no congruence obstruction can exist — the real object is the non-trivial locus).
"""
from __future__ import annotations
from fractions import Fraction as F
import sympy as sp
from kk_heptagon import H1, D1, D2, reconstruct_full

C = 2002


def get_polys(points, D):
    n = len(points); a = D[0][1]
    rho = [D[0][i] ** 2 for i in range(n)]
    r, s = sp.symbols('r s'); tt = sp.symbols('t2 t3 t4 t5 t6')
    x = (sp.Integer(a) ** 2 + r ** 2 - s ** 2) / (2 * sp.Integer(a))
    def Q(fr): return sp.Rational(fr.numerator, fr.denominator)
    m = 2; xm, ym = Q(points[m][0]), Q(points[m][1])
    Am = r ** 2 + rho[m] - 2 * xm * x
    M = sp.expand((Am - tt[0] ** 2) ** 2 - 4 * C * ym ** 2 * (r ** 2 - x ** 2))
    M = sp.expand(M * sp.denom(sp.together(M)))
    Ls = []
    for idx, i in enumerate(range(3, n)):
        xi, yi = Q(points[i][0]), Q(points[i][1])
        ti = tt[idx + 1]
        Li = yi * tt[0] ** 2 - ym * ti ** 2 - ((yi - ym) * r ** 2 + (yi * rho[m] - ym * rho[i]) - 2 * x * (yi * xm - ym * xi))
        Ls.append(sp.expand(Li * sp.denom(sp.together(Li))))
    return (r, s, tt, M, Ls, a, rho)


def qr_table(p):
    return {(v * v) % p for v in range(p)}  # includes 0


def sieve_prime(polys, p):
    r, s, tt, M, Ls, a, rho = polys
    # integer-coefficient lambdas mod p
    Mf = sp.lambdify((r, s, tt[0]), M, 'math')
    # For each L_i solve for t_i^2 coefficient: L_i = c_i * t_i^2 + rest(r,s,t2). Extract.
    lin = []
    for k, Li in enumerate(Ls):
        ti = tt[k + 1]
        ci = sp.Poly(Li, ti).coeff_monomial(ti ** 2)
        rest = sp.expand(Li - ci * ti ** 2)  # rest(r,s,t2); L_i=0 => c_i t_i^2 = -rest
        restf = sp.lambdify((r, s, tt[0]), rest, 'math')
        lin.append((int(ci) % p, restf))
    QR = qr_table(p)
    survivors = 0
    for rr in range(p):
        for ss in range(p):
            for zz in range(p):
                if int(Mf(rr, ss, zz)) % p != 0:
                    continue
                ok = True
                for ci, restf in lin:
                    need = (-int(restf(rr, ss, zz))) % p  # = c_i * t_i^2 mod p
                    if ci % p != 0:
                        inv = pow(ci % p, p - 2, p)
                        val = (need * inv) % p
                        if val not in QR:
                            ok = False; break
                    else:
                        # c_i ≡ 0: need must be ≡ 0 for a solution (t_i free), else none
                        if need % p != 0:
                            ok = False; break
                if ok:
                    survivors += 1
    return survivors


def trivial_points(points, D):
    """The 7 embeddings X=P_j: (r,s,t2,t3,t4,t5,t6) = distances from P_j to (P0,P1,P2,...,P6)."""
    n = len(points)
    out = []
    for j in range(n):
        r = D[0][j]; s = D[1][j]
        ts = [D[i][j] for i in range(2, n)]
        out.append((r, s, ts))
    return out


def check_trivial_survive(polys, points, D, primes):
    r, s, tt, M, Ls, a, rho = polys
    Mf = sp.lambdify((r, s, tt[0]), M, 'math')
    Lfs = [sp.lambdify((r, s, tt[0], tt[k + 1]), Ls[k], 'math') for k in range(len(Ls))]
    triv = trivial_points(points, D)
    allgood = True
    for j, (rr, ss, ts) in enumerate(triv):
        t2 = ts[0]
        mval = int(Mf(rr, ss, t2))
        lvals = [int(Lfs[k](rr, ss, t2, ts[k + 1])) for k in range(len(Ls))]
        exact_zero = (mval == 0 and all(v == 0 for v in lvals))
        allgood = allgood and exact_zero
        print(f"  trivial X=P_{j}: M=0? {mval==0}  all L=0? {all(v==0 for v in lvals)}  -> exact solution: {exact_zero}")
    return allgood


if __name__ == "__main__":
    for label, pts, D in [("H1", H1, D1), ("H2", reconstruct_full(D2, 2002), D2)]:
        print(f"\n================= {label} =================")
        polys = get_polys(pts, D)
        print("Do the 7 trivial embeddings X=P_j satisfy the exact system?")
        triv_ok = check_trivial_survive(polys, pts, D, None)
        print(f"ALL 7 trivial points are exact integral solutions: {triv_ok}")
        print("\nCongruence sieve (survivor counts; 0 would prove MAXIMAL):")
        obstruction = None
        for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103]:
            try:
                surv = sieve_prime(polys, p)
            except Exception as e:
                print(f"  p={p}: error {e}"); continue
            tag = "  <-- OBSTRUCTION (no survivors)" if surv == 0 else ""
            if surv == 0 and obstruction is None:
                obstruction = p
            print(f"  p={p:3d}: survivors mod p = {surv}{tag}")
        print(f"\n{label} verdict: " + ("MAXIMAL via congruence obstruction at p=%d" % obstruction if obstruction
              else "NO congruence obstruction found — local methods cannot prove MAXIMAL (trivial embeddings survive)."))
