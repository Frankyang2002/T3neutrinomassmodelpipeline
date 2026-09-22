"""
Deals with after all BSM fields are integrated out
We run the weinberg from matching threshold to low energy
Note this is the flavour matrix of coefficient matrix, the weak structure is extracted previously
"""

from __future__ import annotations
from scipy.integrate import solve_ivp
from dataclasses import dataclass
import sympy as sp
import numpy as np



g2 = sp.Symbol("g2")
lambdaH = sp.Symbol("lambdaH")


def symbolic_complex_matrix(prefix: str, rows: int, cols: int) -> sp.Matrix:
    """Return a matrix with independent complex symbolic entries.
    it is very general, with no assumptions of symmetries"""

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

    We get the weinberg equation with the 2lambda-3g^2+2T-3/2(YeYedag K+K(YeYedag)^T)
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
    """Reduce the matrix equation to one generation so we can check if its valid"""

    kappa = sp.Symbol("kappa")
    ye, yu, yd = sp.symbols("ye yu yd")

    K = sp.Matrix([[kappa]])
    Ye = sp.Matrix([[ye]])
    Yu = sp.Matrix([[yu]])
    Yd = sp.Matrix([[yd]])

    beta = beta_weinberg_matrix(K, Ye, Yu, Yd)[0, 0]

    return sp.factor(sp.simplify(beta / kappa))


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


# ---------------------------------------------------------------------------
# Numerical Weinberg running
# Consolidated from NumericalWeinbergRGE.py
# ---------------------------------------------------------------------------

r"""
Numerical one-loop SM evolution of the full 3x3 Weinberg coefficient.
- t = ln(mu)
- gY is the ordinary SM hypercharge coupling
- V(H) = (lambdaH / 2) (H^\dagger H)^2.
- We use diagonal yukawa for our evolution, so we have D=diag(y_e^2,y_mu^2,y_tau^2)

- K is the full complex symmetric 3x3 Weinberg coefficient.

The gauge, diagonal Yukawa and Higgs-quartic couplings are evolved
simultaneously at one loop.
"""



LOOP = 16.0 * np.pi**2


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
        '''Package everything to start numerical running'''
        ye = np.asarray(self.ye, dtype=float)
        yu = np.asarray(self.yu, dtype=float)
        yd = np.asarray(self.yd, dtype=float)
        K = np.asarray(self.K, dtype=complex)

        for name, arr in (("ye", ye), ("yu", yu), ("yd", yd)):
            if arr.shape != (3,):
                raise ValueError(f"{name} must contain exactly three entries.")

        if K.shape != (3, 3):
            raise ValueError("K must be a 3x3 matrix.")

        if not np.allclose(K, K.T, rtol=1e-10, atol=1e-14):
            raise ValueError("K must be symmetric.")

        return SMInitialConditions(
            gY=float(self.gY),
            g2=float(self.g2),
            g3=float(self.g3),
            lambdaH=float(self.lambdaH),
            ye=ye,
            yu=yu,
            yd=yd,
            K=K,
        )


@dataclass(frozen=True)
class NumericalRGEResult:
    '''stored state for after numerical evolution finishes'''
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


def _pack(state: SMInitialConditions) -> np.ndarray:
    """Pack real SM couplings and complex K into one real ODE vector.
    Prepared for scipy.solve.ivp
    We also split K into real components for this"""

    return np.concatenate(
        [
            np.array(
                [state.gY, state.g2, state.g3, state.lambdaH],
                dtype=float,
            ),
            state.ye,
            state.yu,
            state.yd,
            state.K.real.reshape(-1),
            state.K.imag.reshape(-1),
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
    """Unpack the ODE vector.
    It reconstructs our weinberg"""

    gY, g2, g3, lambdaH = y[:4]

    ye = y[4:7]
    yu = y[7:10]
    yd = y[10:13]

    K_real = y[13:22].reshape(3, 3)
    K_imag = y[22:31].reshape(3, 3)
    K = K_real + 1j * K_imag

    return gY, g2, g3, lambdaH, ye, yu, yd, K


def _beta(
    _t: float,
    y: np.ndarray,
) -> np.ndarray:
    """
    One-loop SM beta functions plus the Weinberg-coefficient beta function.
    We evaluate one loop beta functinos for all our parameters
    g_y23, y_eud, lambda_H, K etc
    
    We evaluate the beta function we know, 
    beta = (2lambda-3g^2_2+2T)K-3/2(DK+KD^T) 
    D = diag(ye^2,ymu^2ytau^2)   """

    gY, g2, g3, lambdaH, ye, yu, yd, K = _unpack(y)

    ye2 = ye**2
    yu2 = yu**2
    yd2 = yd**2

    T = (
        np.sum(ye2)
        + 3.0 * np.sum(yu2)
        + 3.0 * np.sum(yd2)
    )

    H4 = (
        np.sum(ye2**2)
        + 3.0 * np.sum(yu2**2)
        + 3.0 * np.sum(yd2**2)
    )

    beta_gY = (41.0 / 6.0) * gY**3
    beta_g2 = (-19.0 / 6.0) * g2**3
    beta_g3 = -7.0 * g3**3

    beta_yu = yu * (
        1.5 * (yu2 - yd2)
        + T
        - (17.0 / 12.0) * gY**2
        - (9.0 / 4.0) * g2**2
        - 8.0 * g3**2
    )

    beta_yd = yd * (
        1.5 * (yd2 - yu2)
        + T
        - (5.0 / 12.0) * gY**2
        - (9.0 / 4.0) * g2**2
        - 8.0 * g3**2
    )

    beta_ye = ye * (
        1.5 * ye2
        + T
        - (15.0 / 4.0) * gY**2
        - (9.0 / 4.0) * g2**2
    )

    beta_lambdaH = (
        12.0 * lambdaH**2
        + lambdaH * (
            -3.0 * gY**2
            - 9.0 * g2**2
            + 4.0 * T
        )
        + (3.0 / 4.0) * gY**4
        + (3.0 / 2.0) * gY**2 * g2**2
        + (9.0 / 4.0) * g2**4
        - 4.0 * H4
    )

    De = np.diag(ye2)

    beta_K = (
        (2.0 * lambdaH - 3.0 * g2**2 + 2.0 * T) * K
        - 1.5 * (
            De @ K
            + K @ De.T
        )
    )

    return np.concatenate(
        [
            np.array(
                [
                    beta_gY,
                    beta_g2,
                    beta_g3,
                    beta_lambdaH,
                ],
                dtype=float,
            ),
            beta_ye,
            beta_yu,
            beta_yd,
            beta_K.real.reshape(-1),
            beta_K.imag.reshape(-1),
        ]
    ) / LOOP


def evolve_weinberg(
    initial: SMInitialConditions,
    mu_initial: float,
    mu_final: float,
    *,
    rtol: float = 1e-8,
    atol: float = 1e-11,
) -> NumericalRGEResult:
    """Evolve all one-loop SM couplings and K between two positive scales from mu initial to mu final"""

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
        gY,
        g2,
        g3,
        lambdaH,
        ye,
        yu,
        yd,
        K,
    ) = _unpack(solution.y[:, -1])

    # Numerical integration can introduce tiny antisymmetric roundoff pieces.
    K = 0.5 * (K + K.T)

    return NumericalRGEResult(
        mu_initial=float(mu_initial),
        mu_final=float(mu_final),
        gY=float(gY),
        g2=float(g2),
        g3=float(g3),
        lambdaH=float(lambdaH),
        ye=ye.copy(),
        yu=yu.copy(),
        yd=yd.copy(),
        K=K,
        solver_success=bool(solution.success),
        solver_message=str(solution.message),
        nfev=int(solution.nfev),
    )


def neutrino_mass_matrix(
    K: np.ndarray,
    *,
    vev_gev: float = 246.22,
) -> np.ndarray:
    """Return m_nu=-v^2 K in the matched Matchete convention."""

    K = np.asarray(K, dtype=complex)

    if K.shape != (3, 3):
        raise ValueError("K must be a 3x3 matrix.")

    return -vev_gev**2 * K
