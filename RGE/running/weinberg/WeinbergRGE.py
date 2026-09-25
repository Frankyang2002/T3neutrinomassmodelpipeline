"""SMEFT Weinberg-operator beta functions.

This module defines the one-loop RGE model for the final SM + Weinberg EFT.
It contains beta-function algebra only. Numerical state packing, ODE integration,
trajectory storage, and benchmark substitution live under ``Numerical/``.

Project convention:
    L_EFT contains (1/2) C5 O5 + h.c.
    m_nu = -(v^2/2) C5.
"""

from __future__ import annotations

import numpy as np
import sympy as sp


g2 = sp.Symbol("g2")
lambdaH = sp.Symbol("lambdaH")
LOOP_FACTOR = 16.0 * np.pi**2


def symbolic_complex_matrix(prefix: str, rows: int, cols: int) -> sp.Matrix:
    """Return a symbolic complex matrix with independent named entries."""
    return sp.Matrix(
        rows,
        cols,
        lambda i, j: sp.Symbol(f"{prefix}{i + 1}{j + 1}"),
    )


def symbolic_symmetric_matrix(prefix: str, size: int) -> sp.Matrix:
    """Return a symbolic symmetric matrix."""
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
    return sp.conjugate(matrix.T)


def yukawa_trace(
    Ye: sp.MatrixBase,
    Yu: sp.MatrixBase,
    Yd: sp.MatrixBase,
) -> sp.Expr:
    """Return Tr(Ye Ye† + 3 Yu Yu† + 3 Yd Yd†)."""
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
    """Return 16 pi^2 beta_K for the full three-generation SMEFT."""

    if K.rows != K.cols:
        raise ValueError("K must be square.")

    n = K.rows
    for name, matrix in (("Ye", Ye), ("Yu", Yu), ("Yd", Yd)):
        if matrix.rows != n or matrix.cols != n:
            raise ValueError(f"{name} must be {n}x{n} to match K.")

    T = yukawa_trace(Ye, Yu, Yd)
    YeYeDag = Ye * dagger(Ye)

    beta = (
        (2 * higgs_quartic - 3 * weak_coupling**2 + 2 * T) * K
        - sp.Rational(3, 2)
        * (YeYeDag * K + K * YeYeDag.T)
    )

    if expand_result:
        return beta.applyfunc(sp.expand)
    return beta


def one_generation_reduction() -> sp.Expr:
    """Return the one-generation beta_K/K reduction."""
    kappa = sp.Symbol("kappa")
    ye, yu, yd = sp.symbols("ye yu yd")
    K = sp.Matrix([[kappa]])
    Ye = sp.Matrix([[ye]])
    Yu = sp.Matrix([[yu]])
    Yd = sp.Matrix([[yd]])
    beta = beta_weinberg_matrix(K, Ye, Yu, Yd)[0, 0]
    return sp.factor(sp.simplify(beta / kappa))


def _trace_yy_dagger(matrix: np.ndarray) -> float:
    value = np.trace(matrix @ matrix.conj().T)
    return float(np.real_if_close(value, tol=1000).real)


def _trace_yyyy(matrix: np.ndarray) -> float:
    yy = matrix @ matrix.conj().T
    value = np.trace(yy @ yy)
    return float(np.real_if_close(value, tol=1000).real)


def sm_weinberg_beta(
    *,
    gY: float,
    g2_value: float,
    g3: float,
    lambdaH_value: float,
    ye: np.ndarray,
    yu: np.ndarray,
    yd: np.ndarray,
    K: np.ndarray,
) -> tuple[
    float,
    float,
    float,
    float,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Return the one-loop beta numerators for the final SM + Weinberg EFT.

    The returned quantities satisfy

        16*pi^2 dX/dln(mu) = beta_X.

    This function contains the RGE model only; it performs no numerical
    integration and stores no trajectory.
    """

    ye = np.asarray(ye, dtype=complex)
    yu = np.asarray(yu, dtype=complex)
    yd = np.asarray(yd, dtype=complex)
    K = np.asarray(K, dtype=complex)

    ye_ye_dag = ye @ ye.conj().T
    yu_yu_dag = yu @ yu.conj().T
    yd_yd_dag = yd @ yd.conj().T

    T = (
        _trace_yy_dagger(ye)
        + 3.0 * _trace_yy_dagger(yu)
        + 3.0 * _trace_yy_dagger(yd)
    )
    H4 = (
        _trace_yyyy(ye)
        + 3.0 * _trace_yyyy(yu)
        + 3.0 * _trace_yyyy(yd)
    )

    beta_gY = (41.0 / 6.0) * gY**3
    beta_g2 = (-19.0 / 6.0) * g2_value**3
    beta_g3 = -7.0 * g3**3

    beta_yu = (
        (
            T
            - (17.0 / 12.0) * gY**2
            - (9.0 / 4.0) * g2_value**2
            - 8.0 * g3**2
        )
        * yu
        + 1.5 * (yu_yu_dag @ yu - yd_yd_dag @ yu)
    )

    beta_yd = (
        (
            T
            - (5.0 / 12.0) * gY**2
            - (9.0 / 4.0) * g2_value**2
            - 8.0 * g3**2
        )
        * yd
        + 1.5 * (yd_yd_dag @ yd - yu_yu_dag @ yd)
    )

    beta_ye = (
        (
            T
            - (15.0 / 4.0) * gY**2
            - (9.0 / 4.0) * g2_value**2
        )
        * ye
        + 1.5 * (ye_ye_dag @ ye)
    )

    beta_lambdaH = (
        12.0 * lambdaH_value**2
        + lambdaH_value
        * (
            -3.0 * gY**2
            - 9.0 * g2_value**2
            + 4.0 * T
        )
        + (3.0 / 4.0) * gY**4
        + (3.0 / 2.0) * gY**2 * g2_value**2
        + (9.0 / 4.0) * g2_value**4
        - 4.0 * H4
    )

    beta_K = (
        (
            2.0 * lambdaH_value
            - 3.0 * g2_value**2
            + 2.0 * T
        )
        * K
        - 1.5
        * (
            ye_ye_dag @ K
            + K @ ye_ye_dag.T
        )
    )

    return (
        float(beta_gY),
        float(beta_g2),
        float(beta_g3),
        float(beta_lambdaH),
        np.asarray(beta_ye, dtype=complex),
        np.asarray(beta_yu, dtype=complex),
        np.asarray(beta_yd, dtype=complex),
        np.asarray(beta_K, dtype=complex),
    )
