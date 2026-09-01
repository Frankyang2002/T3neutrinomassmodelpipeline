from __future__ import annotations

"""
Full three-generation one-loop SMEFT RGE for the Weinberg coefficient.

Convention:
    V(H) = (lambdaH / 2) (H^\dagger H)^2

With this convention,

    16 pi^2 beta_K =
        (2 lambdaH - 3 g2^2 + 2 T) K
        - 3/2 [Ye Ye^\dagger K + K (Ye Ye^\dagger)^T],

where

    T = Tr(
        Ye Ye^\dagger
        + 3 Yu Yu^\dagger
        + 3 Yd Yd^\dagger
    ).

K is a complex symmetric 3 x 3 matrix.
"""

import sympy as sp


g2 = sp.Symbol("g2")
lambdaH = sp.Symbol("lambdaH")


def symbolic_complex_matrix(prefix: str, rows: int, cols: int) -> sp.Matrix:
    """Return a matrix with independent complex symbolic entries."""

    return sp.Matrix(
        rows,
        cols,
        lambda i, j: sp.Symbol(f"{prefix}{i + 1}{j + 1}"),
    )


def symbolic_symmetric_matrix(prefix: str, size: int) -> sp.Matrix:
    """Return a symmetric matrix with independent upper-triangle entries."""

    entries: dict[tuple[int, int], sp.Symbol] = {}

    for i in range(size):
        for j in range(i, size):
            entries[(i, j)] = sp.Symbol(f"{prefix}{i + 1}{j + 1}")

    return sp.Matrix(
        size,
        size,
        lambda i, j: entries[(min(i, j), max(i, j))],
    )


def dagger(matrix: sp.MatrixBase) -> sp.Matrix:
    """Hermitian conjugate."""

    return sp.conjugate(matrix.T)


def yukawa_trace(
    Ye: sp.MatrixBase,
    Yu: sp.MatrixBase,
    Yd: sp.MatrixBase,
) -> sp.Expr:
    """Return T = Tr(Ye Ye^dag + 3 Yu Yu^dag + 3 Yd Yd^dag)."""

    return sp.expand(
        sp.trace(Ye * dagger(Ye))
        + 3 * sp.trace(Yu * dagger(Yu))
        + 3 * sp.trace(Yd * dagger(Yd))
    )


def beta_weinberg_matrix(
    K: sp.MatrixBase,
    Ye: sp.MatrixBase,
    Yu: sp.MatrixBase,
    Yd: sp.MatrixBase,
    *,
    higgs_quartic: sp.Expr = lambdaH,
    weak_coupling: sp.Expr = g2,
    expand_result: bool = True,
) -> sp.Matrix:
    """
    Return 16 pi^2 beta_K for the full three-generation SMEFT.

    K must be square and symmetric.  The Yukawa matrices must have compatible
    three-generation flavor dimensions.
    """

    if K.rows != K.cols:
        raise ValueError("K must be square.")

    n = K.rows

    for name, matrix in (("Ye", Ye), ("Yu", Yu), ("Yd", Yd)):
        if matrix.rows != n or matrix.cols != n:
            raise ValueError(
                f"{name} must be {n}x{n} to match K."
            )

    T = yukawa_trace(Ye, Yu, Yd)
    YeYeDag = Ye * dagger(Ye)

    universal = (
        2 * higgs_quartic
        - 3 * weak_coupling**2
        + 2 * T
    )

    beta = (
        universal * K
        - sp.Rational(3, 2) * (
            YeYeDag * K
            + K * YeYeDag.T
        )
    )

    if expand_result:
        return beta.applyfunc(sp.expand)

    return beta


def one_generation_reduction() -> sp.Expr:
    """Reduce the matrix equation to the validated one-generation benchmark."""

    kappa = sp.Symbol("kappa")
    ye, yu, yd = sp.symbols("ye yu yd")

    K = sp.Matrix([[kappa]])
    Ye = sp.Matrix([[ye]])
    Yu = sp.Matrix([[yu]])
    Yd = sp.Matrix([[yd]])

    beta = beta_weinberg_matrix(K, Ye, Yu, Yd)[0, 0]

    return sp.factor(sp.simplify(beta / kappa))


def build_symbolic_three_generation_system() -> dict:
    """Construct symbolic K, Ye, Yu, Yd and their full beta matrix."""

    K = symbolic_symmetric_matrix("k", 3)
    Ye = symbolic_complex_matrix("ye", 3, 3)
    Yu = symbolic_complex_matrix("yu", 3, 3)
    Yd = symbolic_complex_matrix("yd", 3, 3)

    beta = beta_weinberg_matrix(K, Ye, Yu, Yd)

    return {
        "K": K,
        "Ye": Ye,
        "Yu": Yu,
        "Yd": Yd,
        "T": yukawa_trace(Ye, Yu, Yd),
        "beta": beta,
    }


if __name__ == "__main__":
    expected = (
        -3 * g2**2
        + 2 * lambdaH
        + 6 * sp.Symbol("yu") * sp.conjugate(sp.Symbol("yu"))
        + 6 * sp.Symbol("yd") * sp.conjugate(sp.Symbol("yd"))
        - sp.Symbol("ye") * sp.conjugate(sp.Symbol("ye"))
    )

    reduced = one_generation_reduction()

    print("One-generation reduction:")
    print(reduced)
    print()
    print("Expected:")
    print(expected)
    print()
    print("Difference:")
    print(sp.simplify(reduced - expected))
