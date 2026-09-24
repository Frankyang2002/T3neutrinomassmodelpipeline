"""Canonical numerical states for the intermediate EFT after F is removed.

This file mirrors ``T3RGBetaBuildEFT1`` and
``T3RGBetaSharedEFT1OneLoopBetas`` in the current symbolic pipeline.

The intermediate renormalisable theories are

ordinary T3:
    SM + S1 + S2

shared-scalar T3:
    SM + S

The heavy fermion F, its Yukawa couplings, and MF are absent.  Higher-
dimensional Wilson coefficients generated at the fermion threshold are not
stored here; they remain part of the existing EFT1 Wilson/matching layer.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from Numerical.State import (
    N_SM_FLAVOR,
    SMNumericalState,
    T3Representation,
)


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


def _require_optional_real(
    name: str,
    value: float | None,
    *,
    required: bool,
) -> float | None:
    if required:
        if value is None:
            raise ValueError(
                f"{name} is required for this SU(2) representation."
            )
        return _finite_real(name, value)

    if value is not None:
        raise ValueError(
            f"{name} is not present for this SU(2) representation."
        )

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
class T3EFT1State:
    """Renormalisable ordinary-T3 state after integrating out F."""

    mu_gev: float
    representation: T3Representation
    sm: SMNumericalState

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

    def validated(self) -> "T3EFT1State":
        rep = self.representation.validated()

        if rep.shared_scalar:
            raise ValueError(
                "T3EFT1State is for split scalars; use SharedT3EFT1State."
            )

        special_doublet_case = (
            rep.d_s1 == 2
            and rep.d_s2 == 2
            and rep.alpha == -1
        )

        return T3EFT1State(
            mu_gev=_positive_scale(self.mu_gev),
            representation=rep,
            sm=self.sm.validated(),
            mS1Sq=_finite_real("mS1Sq", self.mS1Sq),
            mS2Sq=_finite_real("mS2Sq", self.mS2Sq),
            lambdaS1=_finite_real("lambdaS1", self.lambdaS1),
            lambdaS2=_finite_real("lambdaS2", self.lambdaS2),
            lambdaH1=_finite_real("lambdaH1", self.lambdaH1),
            lambdaH2=_finite_real("lambdaH2", self.lambdaH2),
            lambda12=_finite_real("lambda12", self.lambda12),
            lambdaT3=_finite_complex("lambdaT3", self.lambdaT3),
            lambdaH1Adj=_require_optional_real(
                "lambdaH1Adj",
                self.lambdaH1Adj,
                required=rep.d_s1 > 1,
            ),
            lambdaH2Adj=_require_optional_real(
                "lambdaH2Adj",
                self.lambdaH2Adj,
                required=rep.d_s2 > 1,
            ),
            lambdaS1Adj=_require_optional_real(
                "lambdaS1Adj",
                self.lambdaS1Adj,
                required=rep.d_s1 == 3,
            ),
            lambdaS2Adj=_require_optional_real(
                "lambdaS2Adj",
                self.lambdaS2Adj,
                required=rep.d_s2 == 3,
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
class SharedT3EFT1State:
    """Renormalisable shared-scalar/scotogenic state after F is removed."""

    mu_gev: float
    representation: T3Representation
    sm: SMNumericalState

    mSSq: float
    lambdaS: float
    lambda3: float
    lambda4: float
    lambda5: float

    def validated(self) -> "SharedT3EFT1State":
        rep = self.representation.validated()

        if not rep.shared_scalar:
            raise ValueError(
                "SharedT3EFT1State requires shared_scalar=True."
            )

        return SharedT3EFT1State(
            mu_gev=_positive_scale(self.mu_gev),
            representation=rep,
            sm=self.sm.validated(),
            mSSq=_finite_real("mSSq", self.mSSq),
            lambdaS=_finite_real("lambdaS", self.lambdaS),
            lambda3=_finite_real("lambda3", self.lambda3),
            lambda4=_finite_real("lambda4", self.lambda4),
            lambda5=_finite_real("lambda5", self.lambda5),
        )
