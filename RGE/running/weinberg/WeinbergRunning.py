"""
Full-flavour one-loop SM + Weinberg-operator running.

Conventions
-----------
- t = ln(mu)
- gY is the ordinary SM hypercharge coupling.
- V(H) = (lambdaH / 2) (H^\dagger H)^2.
- ye, yu, yd are full complex 3x3 Yukawa matrices.
- K is the full complex symmetric 3x3 Weinberg coefficient.
- Project convention: L_EFT contains (1/2) K O5 + h.c., hence
      m_nu = -(v^2/2) K.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp


g2 = sp.Symbol("g2")
lambdaH = sp.Symbol("lambdaH")
LOOP = 16.0 * np.pi**2


def symbolic_complex_matrix(prefix: str, rows: int, cols: int) -> sp.Matrix:
    return sp.Matrix(
        rows,
        cols,
        lambda i, j: sp.Symbol(f"{prefix}{i + 1}{j + 1}"),
    )


def symbolic_symmetric_matrix(prefix: str, size: int) -> sp.Matrix:
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
    kappa = sp.Symbol("kappa")
    ye, yu, yd = sp.symbols("ye yu yd")
    K = sp.Matrix([[kappa]])
    Ye = sp.Matrix([[ye]])
    Yu = sp.Matrix([[yu]])
    Yd = sp.Matrix([[yd]])
    beta = beta_weinberg_matrix(K, Ye, Yu, Yd)[0, 0]
    return sp.factor(sp.simplify(beta / kappa))


@dataclass(frozen=True)
class SMInitialConditions:
    gY: float
    g2: float
    g3: float
    lambdaH: float
    ye: np.ndarray
    yu: np.ndarray
    yd: np.ndarray
    K: np.ndarray

    def validated(self) -> "SMInitialConditions":
        ye = np.asarray(self.ye, dtype=complex)
        yu = np.asarray(self.yu, dtype=complex)
        yd = np.asarray(self.yd, dtype=complex)
        K = np.asarray(self.K, dtype=complex)

        for name, arr in (("ye", ye), ("yu", yu), ("yd", yd)):
            if arr.shape != (3, 3):
                raise ValueError(f"{name} must be a full complex 3x3 matrix.")
            if not np.all(np.isfinite(arr.real)) or not np.all(
                np.isfinite(arr.imag)
            ):
                raise ValueError(f"{name} must contain only finite entries.")

        if K.shape != (3, 3):
            raise ValueError("K must be a 3x3 matrix.")
        if not np.all(np.isfinite(K.real)) or not np.all(np.isfinite(K.imag)):
            raise ValueError("K must contain only finite entries.")
        if not np.allclose(K, K.T, rtol=1e-10, atol=1e-14):
            raise ValueError("K must be symmetric.")

        scalars = {
            "gY": self.gY,
            "g2": self.g2,
            "g3": self.g3,
            "lambdaH": self.lambdaH,
        }
        for name, value in scalars.items():
            if not np.isfinite(float(value)):
                raise ValueError(f"{name} must be finite.")

        return SMInitialConditions(
            gY=float(self.gY),
            g2=float(self.g2),
            g3=float(self.g3),
            lambdaH=float(self.lambdaH),
            ye=ye.copy(),
            yu=yu.copy(),
            yd=yd.copy(),
            K=K.copy(),
        )


@dataclass(frozen=True)
class NumericalRGEResult:
    mu_initial: float
    mu_final: float
    gY: float
    g2: float
    g3: float
    lambdaH: float
    ye: np.ndarray
    yu: np.ndarray
    yd: np.ndarray
    K: np.ndarray
    solver_success: bool
    solver_message: str
    nfev: int


def _pack_complex_matrix(matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=complex)
    return np.concatenate(
        [matrix.real.reshape(-1), matrix.imag.reshape(-1)]
    )


def _unpack_complex_matrix(
    y: np.ndarray,
    start: int,
) -> tuple[np.ndarray, int]:
    real = y[start : start + 9].reshape(3, 3)
    imag = y[start + 9 : start + 18].reshape(3, 3)
    return real + 1j * imag, start + 18


def _pack(state: SMInitialConditions) -> np.ndarray:
    state = state.validated()
    return np.concatenate(
        [
            np.array(
                [state.gY, state.g2, state.g3, state.lambdaH],
                dtype=float,
            ),
            _pack_complex_matrix(state.ye),
            _pack_complex_matrix(state.yu),
            _pack_complex_matrix(state.yd),
            _pack_complex_matrix(state.K),
        ]
    )


def _unpack(y: np.ndarray) -> tuple[
    float,
    float,
    float,
    float,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    y = np.asarray(y, dtype=float)
    if y.shape != (76,):
        raise ValueError(
            f"Full-flavour Weinberg ODE vector must have length 76; got {y.shape}."
        )

    gY, g2_value, g3, lambdaH_value = y[:4]
    index = 4
    ye, index = _unpack_complex_matrix(y, index)
    yu, index = _unpack_complex_matrix(y, index)
    yd, index = _unpack_complex_matrix(y, index)
    K, index = _unpack_complex_matrix(y, index)

    return (
        float(gY),
        float(g2_value),
        float(g3),
        float(lambdaH_value),
        ye,
        yu,
        yd,
        K,
    )


def _trace_yy_dagger(matrix: np.ndarray) -> float:
    value = np.trace(matrix @ matrix.conj().T)
    return float(np.real_if_close(value, tol=1000).real)


def _trace_yyyy(matrix: np.ndarray) -> float:
    yy = matrix @ matrix.conj().T
    value = np.trace(yy @ yy)
    return float(np.real_if_close(value, tol=1000).real)


def _beta(
    _t: float,
    y: np.ndarray,
) -> np.ndarray:
    """One-loop full-flavour SM beta functions plus beta_K."""

    gY, g2_value, g3, lambdaH_value, ye, yu, yd, K = _unpack(y)

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

    derivative = np.concatenate(
        [
            np.array(
                [beta_gY, beta_g2, beta_g3, beta_lambdaH],
                dtype=float,
            ),
            _pack_complex_matrix(beta_ye),
            _pack_complex_matrix(beta_yu),
            _pack_complex_matrix(beta_yd),
            _pack_complex_matrix(beta_K),
        ]
    )

    return derivative / LOOP


def evolve_weinberg(
    initial: SMInitialConditions,
    mu_initial: float,
    mu_final: float,
    *,
    rtol: float = 1e-8,
    atol: float = 1e-11,
) -> NumericalRGEResult:
    if mu_initial <= 0 or mu_final <= 0:
        raise ValueError("RGE scales must be positive.")

    initial = initial.validated()
    y0 = _pack(initial)

    solution = solve_ivp(
        _beta,
        t_span=(np.log(mu_initial), np.log(mu_final)),
        y0=y0,
        method="DOP853",
        rtol=rtol,
        atol=atol,
    )

    if not solution.success:
        raise RuntimeError(
            f"Numerical RGE integration failed: {solution.message}"
        )

    (
        gY_value,
        g2_out,
        g3_out,
        lambdaH_out,
        ye_out,
        yu_out,
        yd_out,
        K_out,
    ) = _unpack(solution.y[:, -1])

    K_out = 0.5 * (K_out + K_out.T)

    return NumericalRGEResult(
        mu_initial=float(mu_initial),
        mu_final=float(mu_final),
        gY=float(gY_value),
        g2=float(g2_out),
        g3=float(g3_out),
        lambdaH=float(lambdaH_out),
        ye=ye_out.copy(),
        yu=yu_out.copy(),
        yd=yd_out.copy(),
        K=K_out,
        solver_success=bool(solution.success),
        solver_message=str(solution.message),
        nfev=int(solution.nfev),
    )


def neutrino_mass_matrix(
    K: np.ndarray,
    *,
    vev_gev: float = 246.22,
) -> np.ndarray:
    """Return m_nu = -(v^2/2) K in the project convention."""

    K = np.asarray(K, dtype=complex)
    if K.shape != (3, 3):
        raise ValueError("K must be a 3x3 matrix.")
    return -(float(vev_gev) ** 2 / 2.0) * K
