"""Canonical numerical parameter states for the T3 UV theory.

This module defines data only. It deliberately does not integrate RGEs.

The state mirrors the renormalisable couplings exported by
``RGE/running/rgbeta/T3RGBetaModel.wl`` on the current pipeline.

The Higgs quadratic mass is intentionally absent: the current RGBeta T3
front end does not export a beta function for it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from common.T3Model import (
    valid_shared_scalar_dimensions,
    valid_t3_dimensions,
)

N_SM_FLAVOR = 3
N_HEAVY_FLAVOR = 3


def _finite_real(name: str, value: Any) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite.")
    return result


def _finite_complex(name: str, value: Any) -> complex:
    result = complex(value)
    if not (np.isfinite(result.real) and np.isfinite(result.imag)):
        raise ValueError(f"{name} must be finite.")
    return result


def _positive_scale(value: Any) -> float:
    result = _finite_real("mu_gev", value)
    if result <= 0.0:
        raise ValueError("mu_gev must be positive.")
    return result


def _complex_matrix(name: str, value: Any, shape: tuple[int, int]) -> np.ndarray:
    matrix = np.asarray(value, dtype=complex)
    if matrix.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {matrix.shape}.")
    if not np.all(np.isfinite(matrix.real)) or not np.all(np.isfinite(matrix.imag)):
        raise ValueError(f"{name} must contain only finite entries.")
    return matrix.copy()


def _require_optional_real(name: str, value: float | None, *, required: bool) -> float | None:
    if required:
        if value is None:
            raise ValueError(f"{name} is required for this SU(2) representation.")
        return _finite_real(name, value)
    if value is not None:
        raise ValueError(f"{name} is not present for this SU(2) representation.")
    return None


def _require_optional_complex(
    name: str,
    value: complex | None,
    *,
    required: bool,
) -> complex | None:
    if required:
        if value is None:
            raise ValueError(
                f"{name} is required for this SU(2)/hypercharge assignment."
            )
        return _finite_complex(name, value)
    if value is not None:
        raise ValueError(
            f"{name} is not present for this SU(2)/hypercharge assignment."
        )
    return None


@dataclass(frozen=True)
class T3Representation:
    """Discrete representation data selecting one numerical theory."""

    d_s1: int
    d_s2: int
    d_f: int
    alpha: int
    shared_scalar: bool = False

    def validated(self) -> "T3Representation":
        d_s1 = int(self.d_s1)
        d_s2 = int(self.d_s2)
        d_f = int(self.d_f)
        alpha = int(self.alpha)

        if self.shared_scalar:
            if d_s1 != d_s2:
                raise ValueError("Shared-scalar mode requires d_s1 == d_s2.")
            if alpha != -1:
                raise ValueError("Shared-scalar mode requires alpha = -1.")
            if not valid_shared_scalar_dimensions(d_s1, d_f):
                raise ValueError(
                    "Unsupported shared-scalar dimensions. "
                    "Current support is dS=2 with dF=1 or 3."
                )
        elif not valid_t3_dimensions(d_s1, d_s2, d_f):
            raise ValueError(
                "Unsupported ordinary T3 dimensions "
                f"(dS1,dS2,dF)=({d_s1},{d_s2},{d_f})."
            )

        return T3Representation(
            d_s1=d_s1,
            d_s2=d_s2,
            d_f=d_f,
            alpha=alpha,
            shared_scalar=bool(self.shared_scalar),
        )

    @property
    def y_s1(self) -> float:
        return self.alpha / 2.0

    @property
    def y_s2(self) -> float:
        return (self.alpha + 2) / 2.0

    @property
    def y_f(self) -> float:
        return (self.alpha + 1) / 2.0

    @property
    def majorana_fermion(self) -> bool:
        return self.y_f == 0.0


@dataclass(frozen=True)
class SMNumericalState:
    """SM couplings common to all stages."""

    gY: float
    g2: float
    g3: float
    lambdaH: float
    yu: np.ndarray
    yd: np.ndarray
    ye: np.ndarray

    def validated(self) -> "SMNumericalState":
        return SMNumericalState(
            gY=_finite_real("gY", self.gY),
            g2=_finite_real("g2", self.g2),
            g3=_finite_real("g3", self.g3),
            lambdaH=_finite_real("lambdaH", self.lambdaH),
            yu=_complex_matrix("yu", self.yu, (N_SM_FLAVOR, N_SM_FLAVOR)),
            yd=_complex_matrix("yd", self.yd, (N_SM_FLAVOR, N_SM_FLAVOR)),
            ye=_complex_matrix("ye", self.ye, (N_SM_FLAVOR, N_SM_FLAVOR)),
        )


@dataclass(frozen=True)
class T3UVState:
    """Canonical ordinary split-scalar T3 numerical state at one scale."""

    mu_gev: float
    representation: T3Representation
    sm: SMNumericalState

    y1: np.ndarray
    y2: np.ndarray
    MF: np.ndarray

    mS1Sq: float
    mS2Sq: float

    lambdaS1: float
    lambdaS2: float
    lambdaH1: float
    lambdaH2: float
    lambda12: float
    lambdaT3: complex

    lambdaH1Adj: float | None = None
    lambdaH2Adj: float | None = None
    lambdaS1Adj: float | None = None
    lambdaS2Adj: float | None = None
    lambda12Adj: float | None = None
    lambda12Cross: float | None = None

    lambdaHHdagS2S2: complex | None = None
    lambdaHHdagS1barS1bar: complex | None = None
    lambdaS1bar2S2bar2: complex | None = None
    lambdaS1barS2S2bar2: complex | None = None
    lambdaS1S1bar2S2bar: complex | None = None
    lambdaHHdagS1barS2barCross: complex | None = None

    def validated(self) -> "T3UVState":
        rep = self.representation.validated()
        if rep.shared_scalar:
            raise ValueError(
                "T3UVState is for split scalars; use SharedT3UVState instead."
            )

        special_doublet_case = (
            rep.d_s1 == 2 and rep.d_s2 == 2 and rep.alpha == -1
        )

        mf = _complex_matrix(
            "MF",
            self.MF,
            (N_HEAVY_FLAVOR, N_HEAVY_FLAVOR),
        )
        if rep.majorana_fermion and not np.allclose(
            mf, mf.T, rtol=1.0e-10, atol=1.0e-14
        ):
            raise ValueError(
                "MF must be complex symmetric when the heavy fermion is Majorana."
            )

        return T3UVState(
            mu_gev=_positive_scale(self.mu_gev),
            representation=rep,
            sm=self.sm.validated(),
            y1=_complex_matrix("y1", self.y1, (N_SM_FLAVOR, N_HEAVY_FLAVOR)),
            y2=_complex_matrix("y2", self.y2, (N_SM_FLAVOR, N_HEAVY_FLAVOR)),
            MF=mf,
            mS1Sq=_finite_real("mS1Sq", self.mS1Sq),
            mS2Sq=_finite_real("mS2Sq", self.mS2Sq),
            lambdaS1=_finite_real("lambdaS1", self.lambdaS1),
            lambdaS2=_finite_real("lambdaS2", self.lambdaS2),
            lambdaH1=_finite_real("lambdaH1", self.lambdaH1),
            lambdaH2=_finite_real("lambdaH2", self.lambdaH2),
            lambda12=_finite_real("lambda12", self.lambda12),
            lambdaT3=_finite_complex("lambdaT3", self.lambdaT3),
            lambdaH1Adj=_require_optional_real(
                "lambdaH1Adj", self.lambdaH1Adj, required=rep.d_s1 > 1
            ),
            lambdaH2Adj=_require_optional_real(
                "lambdaH2Adj", self.lambdaH2Adj, required=rep.d_s2 > 1
            ),
            lambdaS1Adj=_require_optional_real(
                "lambdaS1Adj", self.lambdaS1Adj, required=rep.d_s1 == 3
            ),
            lambdaS2Adj=_require_optional_real(
                "lambdaS2Adj", self.lambdaS2Adj, required=rep.d_s2 == 3
            ),
            lambda12Adj=_require_optional_real(
                "lambda12Adj",
                self.lambda12Adj,
                required=rep.d_s1 > 1 and rep.d_s2 > 1,
            ),
            lambda12Cross=_require_optional_real(
                "lambda12Cross",
                self.lambda12Cross,
                required=rep.d_s1 == 3 and rep.d_s2 == 3,
            ),
            lambdaHHdagS2S2=_require_optional_complex(
                "lambdaHHdagS2S2",
                self.lambdaHHdagS2S2,
                required=special_doublet_case,
            ),
            lambdaHHdagS1barS1bar=_require_optional_complex(
                "lambdaHHdagS1barS1bar",
                self.lambdaHHdagS1barS1bar,
                required=special_doublet_case,
            ),
            lambdaS1bar2S2bar2=_require_optional_complex(
                "lambdaS1bar2S2bar2",
                self.lambdaS1bar2S2bar2,
                required=special_doublet_case,
            ),
            lambdaS1barS2S2bar2=_require_optional_complex(
                "lambdaS1barS2S2bar2",
                self.lambdaS1barS2S2bar2,
                required=special_doublet_case,
            ),
            lambdaS1S1bar2S2bar=_require_optional_complex(
                "lambdaS1S1bar2S2bar",
                self.lambdaS1S1bar2S2bar,
                required=special_doublet_case,
            ),
            lambdaHHdagS1barS2barCross=_require_optional_complex(
                "lambdaHHdagS1barS2barCross",
                self.lambdaHHdagS1barS2barCross,
                required=special_doublet_case,
            ),
        )


@dataclass(frozen=True)
class SharedT3UVState:
    """Canonical shared-scalar/scotogenic numerical UV state."""

    mu_gev: float
    representation: T3Representation
    sm: SMNumericalState

    h: np.ndarray
    MF: np.ndarray
    mSSq: float

    lambdaS: float
    lambda3: float
    lambda4: float
    lambda5: float

    def validated(self) -> "SharedT3UVState":
        rep = self.representation.validated()
        if not rep.shared_scalar:
            raise ValueError(
                "SharedT3UVState requires representation.shared_scalar=True."
            )

        mf = _complex_matrix(
            "MF",
            self.MF,
            (N_HEAVY_FLAVOR, N_HEAVY_FLAVOR),
        )
        if not np.allclose(mf, mf.T, rtol=1.0e-10, atol=1.0e-14):
            raise ValueError(
                "MF must be complex symmetric in the supported shared-scalar models."
            )

        return SharedT3UVState(
            mu_gev=_positive_scale(self.mu_gev),
            representation=rep,
            sm=self.sm.validated(),
            h=_complex_matrix("h", self.h, (N_SM_FLAVOR, N_HEAVY_FLAVOR)),
            MF=mf,
            mSSq=_finite_real("mSSq", self.mSSq),
            lambdaS=_finite_real("lambdaS", self.lambdaS),
            lambda3=_finite_real("lambda3", self.lambda3),
            lambda4=_finite_real("lambda4", self.lambda4),
            lambda5=_finite_real("lambda5", self.lambda5),
        )
