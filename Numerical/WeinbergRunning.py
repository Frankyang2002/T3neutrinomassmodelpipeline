"""Numerical integration of the final SM + Weinberg EFT.

The RGE equations themselves live in
``RGE.running.weinberg.WeinbergRGE``. This module owns numerical state
validation, packing/unpacking, ``solve_ivp`` integration, and endpoint results.

Conventions
-----------
- t = ln(mu)
- gY is the ordinary SM hypercharge coupling.
- V(H) = (lambdaH / 2) (H^dagger H)^2.
- ye, yu, yd are full complex 3x3 Yukawa matrices.
- K is the full complex symmetric 3x3 Weinberg coefficient.
- Project convention: L_EFT contains (1/2) K O5 + h.c., hence
      m_nu = -(v^2/2) K.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from RGE.running.weinberg.WeinbergRGE import (
    LOOP_FACTOR,
    sm_weinberg_beta,
)


# Historical public constant retained for tests and callers.
LOOP = LOOP_FACTOR


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


def _beta(
    _t: float,
    y: np.ndarray,
) -> np.ndarray:
    """Return d(state-vector)/dln(mu) from the final-EFT RGE model."""

    gY, g2_value, g3, lambdaH_value, ye, yu, yd, K = _unpack(y)

    (
        beta_gY,
        beta_g2,
        beta_g3,
        beta_lambdaH,
        beta_ye,
        beta_yu,
        beta_yd,
        beta_K,
    ) = sm_weinberg_beta(
        gY=gY,
        g2_value=g2_value,
        g3=g3,
        lambdaH_value=lambdaH_value,
        ye=ye,
        yu=yu,
        yd=yd,
        K=K,
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
    """Numerically evolve the final SM + Weinberg EFT between two scales."""

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
