from __future__ import annotations

import sympy as sp

from RGE.SMEFTWeinbergFlavorRGE import (
    beta_weinberg_matrix,
    g2,
    lambdaH,
    one_generation_reduction,
    symbolic_complex_matrix,
    symbolic_symmetric_matrix,
)


def main() -> int:
    ye, yu, yd = sp.symbols("ye yu yd")

    expected_1g = (
        -3 * g2**2
        + 2 * lambdaH
        + 6 * yu * sp.conjugate(yu)
        + 6 * yd * sp.conjugate(yd)
        - ye * sp.conjugate(ye)
    )

    reduced = one_generation_reduction()
    diff_1g = sp.simplify(reduced - expected_1g)

    print("=" * 72)
    print("FULL-FLAVOR SMEFT WEINBERG RGE REGRESSION")
    print("=" * 72)
    print()
    print("One-generation reduction:")
    print(reduced)
    print()
    print("Expected:")
    print(expected_1g)
    print()
    print("Difference:")
    print(diff_1g)
    print()

    if diff_1g != 0:
        print("FAIL: one-generation reduction")
        return 1

    K = symbolic_symmetric_matrix("k", 3)
    Ye = symbolic_complex_matrix("ye", 3, 3)
    Yu = symbolic_complex_matrix("yu", 3, 3)
    Yd = symbolic_complex_matrix("yd", 3, 3)

    beta = beta_weinberg_matrix(K, Ye, Yu, Yd)

    symmetry_difference = (
        beta - beta.T
    ).applyfunc(sp.simplify)

    print("Symmetry check beta_K - beta_K^T:")
    print(symmetry_difference)
    print()

    if any(entry != 0 for entry in symmetry_difference):
        print("FAIL: beta_K is not symmetric")
        return 1

    # Diagonal-Yukawa check.  For diagonal Ye, each K_pq component must receive
    # the expected flavor-dependent factor
    #
    #   -3/2 (|ye_p|^2 + |ye_q|^2) K_pq
    #
    # in addition to the universal trace/gauge/quartic part.
    e1, e2, e3 = sp.symbols("e1 e2 e3")
    u1, u2, u3 = sp.symbols("u1 u2 u3")
    d1, d2, d3 = sp.symbols("d1 d2 d3")

    Ye_diag = sp.diag(e1, e2, e3)
    Yu_diag = sp.diag(u1, u2, u3)
    Yd_diag = sp.diag(d1, d2, d3)

    beta_diag = beta_weinberg_matrix(
        K,
        Ye_diag,
        Yu_diag,
        Yd_diag,
    )

    T = (
        sum(x * sp.conjugate(x) for x in (e1, e2, e3))
        + 3 * sum(x * sp.conjugate(x) for x in (u1, u2, u3))
        + 3 * sum(x * sp.conjugate(x) for x in (d1, d2, d3))
    )

    universal = 2 * lambdaH - 3 * g2**2 + 2 * T
    charged = (e1, e2, e3)

    for p in range(3):
        for q in range(3):
            expected_component = (
                universal
                - sp.Rational(3, 2) * (
                    charged[p] * sp.conjugate(charged[p])
                    + charged[q] * sp.conjugate(charged[q])
                )
            ) * K[p, q]

            difference = sp.simplify(
                beta_diag[p, q] - expected_component
            )

            if difference != 0:
                print(
                    f"FAIL: diagonal flavor component ({p + 1},{q + 1})"
                )
                print("Difference:", difference)
                return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
