from __future__ import annotations

"""
Validate the mixed EFT1 Wilson-coefficient flavor convention.

For a Majorana heavy fermion F, integrating F out at tree level gives

    L_eff ~ (1 / (2 M_F)) J J + h.c.

with

    J_r = sum_p [ y1^*_{pr} (L_p S1) + y2^*_{pr} (L_p S2) ].

The mixed S1-S2 piece is therefore

    C12_{pq}
      = (1 / (2 M_F)) sum_r [
            y1^*_{pr} y2^*_{qr}
          + y2^*_{pr} y1^*_{qr}
        ].

This is symmetric in p <-> q and reduces in one generation to

    C12 = y1^* y2^* / M_F,

which is exactly the normalization already used by the one-generation
Matchete EFT1 seed.

This file performs symbolic and non-diagonal numerical checks of that
convention.  It does not modify matching output.
"""

import argparse
import json
from pathlib import Path
from typing import Any

import sympy as sp


def build_mixed_c12(
    y1: sp.MatrixBase,
    y2: sp.MatrixBase,
    masses: list[sp.Expr] | tuple[sp.Expr, ...],
) -> sp.Matrix:
    """Return the symmetric mixed Wilson matrix C12^{pq}."""
    if y1.shape != y2.shape:
        raise ValueError("y1 and y2 must have the same shape.")

    n_lepton, n_heavy = y1.shape
    if len(masses) != n_heavy:
        raise ValueError("Need one heavy-fermion mass per heavy generation.")

    result = sp.MutableDenseMatrix.zeros(n_lepton, n_lepton)

    for p in range(n_lepton):
        for q in range(n_lepton):
            value = sp.S.Zero
            for r in range(n_heavy):
                value += sp.Rational(1, 2) / masses[r] * (
                    sp.conjugate(y1[p, r]) * sp.conjugate(y2[q, r])
                    + sp.conjugate(y2[p, r]) * sp.conjugate(y1[q, r])
                )
            result[p, q] = sp.simplify(value)

    return sp.Matrix(result)


def build_ordered_kernel(
    y1: sp.MatrixBase,
    y2: sp.MatrixBase,
    masses: list[sp.Expr] | tuple[sp.Expr, ...],
) -> sp.Matrix:
    """Return the unsymmetrized ordered y1*_pr y2*_qr / M_r kernel."""
    n_lepton, n_heavy = y1.shape
    result = sp.MutableDenseMatrix.zeros(n_lepton, n_lepton)

    for p in range(n_lepton):
        for q in range(n_lepton):
            result[p, q] = sp.simplify(
                sum(
                    sp.conjugate(y1[p, r])
                    * sp.conjugate(y2[q, r])
                    / masses[r]
                    for r in range(n_heavy)
                )
            )

    return sp.Matrix(result)


def symbolic_check() -> dict[str, Any]:
    y1p, y1q, y2p, y2q, mf = sp.symbols(
        "y1p y1q y2p y2q MF",
        nonzero=True,
    )

    mixed_pq = sp.Rational(1, 2) / mf * (
        sp.conjugate(y1p) * sp.conjugate(y2q)
        + sp.conjugate(y2p) * sp.conjugate(y1q)
    )
    mixed_qp = sp.Rational(1, 2) / mf * (
        sp.conjugate(y1q) * sp.conjugate(y2p)
        + sp.conjugate(y2q) * sp.conjugate(y1p)
    )

    y1, y2 = sp.symbols("y1 y2")
    one_generation = sp.simplify(
        sp.Rational(1, 2) / mf * (
            sp.conjugate(y1) * sp.conjugate(y2)
            + sp.conjugate(y2) * sp.conjugate(y1)
        )
    )
    expected_1g = sp.conjugate(y1) * sp.conjugate(y2) / mf

    return {
        "pq_qp_symmetry": sp.simplify(mixed_pq - mixed_qp) == 0,
        "one_generation_reduction": sp.simplify(
            one_generation - expected_1g
        ) == 0,
        "mixed_pq": str(mixed_pq),
        "one_generation": str(one_generation),
    }


def nondiagonal_check() -> dict[str, Any]:
    """Use deliberately non-diagonal complex Yukawas."""
    y1 = sp.Matrix(
        [
            [1 + sp.I, 2],
            [3, 1 - 2 * sp.I],
            [2 - sp.I, -1],
        ]
    )
    y2 = sp.Matrix(
        [
            [2, -sp.I],
            [1 + sp.I, 4],
            [-2, 3 + sp.I],
        ]
    )
    masses = [sp.Integer(5), sp.Integer(7)]

    c12 = build_mixed_c12(y1, y2, masses)
    ordered = build_ordered_kernel(y1, y2, masses)

    symmetry_residual = sp.simplify(c12 - c12.T)
    ordered_antisymmetric = sp.simplify(ordered - ordered.T)

    # The symmetrized coefficient must equal 1/2(K + K^T).
    reconstruction = sp.simplify(c12 - sp.Rational(1, 2) * (ordered + ordered.T))

    return {
        "c12_is_symmetric": symmetry_residual == sp.zeros(3),
        "ordered_kernel_is_generically_nonsymmetric": (
            ordered_antisymmetric != sp.zeros(3)
        ),
        "equals_half_ordered_plus_transpose": reconstruction == sp.zeros(3),
        "c12": [[str(sp.simplify(x)) for x in c12.row(i)] for i in range(3)],
        "ordered_antisymmetric_residual": [
            [str(sp.simplify(x)) for x in ordered_antisymmetric.row(i)]
            for i in range(3)
        ],
    }


def run_validation(output: Path | None = None) -> dict[str, Any]:
    symbolic = symbolic_check()
    numeric = nondiagonal_check()

    passed = all(
        [
            symbolic["pq_qp_symmetry"],
            symbolic["one_generation_reduction"],
            numeric["c12_is_symmetric"],
            numeric["ordered_kernel_is_generically_nonsymmetric"],
            numeric["equals_half_ordered_plus_transpose"],
        ]
    )

    result = {
        "status": "Success" if passed else "Failed",
        "convention": (
            "C12[p,q] = 1/2 sum_r "
            "(conjugate(y1[p,r])*conjugate(y2[q,r]) + "
            "conjugate(y2[p,r])*conjugate(y1[q,r]))/MF_r"
        ),
        "derivation": (
            "The Majorana tree-level EFT is (1/2) J M_F^{-1} J. "
            "Expanding J = y1^*(L S1) + y2^*(L S2) gives two mixed "
            "cross terms. Combining them in the symmetric psi^2 phi^2 "
            "Wilson convention yields the 1/2(y1*y2 + y2*y1) flavor kernel."
        ),
        "symbolic_checks": symbolic,
        "nondiagonal_checks": numeric,
        "implication": (
            "The previous ordered C12[p,q] = sum_r y1^*_{pr} y2^*_{qr}/MF_r "
            "is sufficient only in one generation. Full flavor must use the "
            "symmetrized kernel above."
        ),
    }

    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, indent=2), encoding="utf-8")

    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("c12_flavor_symmetry_validation.json"),
    )
    args = parser.parse_args()

    result = run_validation(args.output)
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "Success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
