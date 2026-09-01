from __future__ import annotations

"""
Numerical one-loop SM evolution of the full 3x3 Weinberg coefficient.

Conventions
-----------
- t = ln(mu)
- gY is the ordinary SM hypercharge coupling, not GUT-normalized g1.
- V(H) = (lambdaH / 2) (H^\dagger H)^2.
- The charged-lepton, up-quark and down-quark Yukawa matrices are taken
  diagonal during the numerical evolution.
- K is the full complex symmetric 3x3 Weinberg coefficient.

The evolved equations are

    16 pi^2 dK/dt =
        (2 lambdaH - 3 g2^2 + 2 T) K
        - 3/2 [De K + K De^T],

where De = diag(|ye_i|^2) and

    T = sum_i |ye_i|^2
        + 3 sum_i |yu_i|^2
        + 3 sum_i |yd_i|^2.

The gauge, diagonal Yukawa and Higgs-quartic couplings are evolved
simultaneously at one loop.
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


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
    """Pack real SM couplings and complex K into one real ODE vector."""

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
    """Unpack the ODE vector."""

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

    The independent variable is t=ln(mu), so the right-hand side is d/dt.
    """

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
        6.0 * lambdaH**2
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
    """Evolve all one-loop SM couplings and K between two positive scales."""

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
    """Return m_nu=-(v^2/2)K in GeV."""

    K = np.asarray(K, dtype=complex)

    if K.shape != (3, 3):
        raise ValueError("K must be a 3x3 matrix.")

    return -vev_gev**2 * K
