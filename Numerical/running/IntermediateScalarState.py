"""Scalar-only EFT numerical state, including native general RGBeta scalar couplings."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Mapping
import numpy as np
from Numerical.core.State import (SMNumericalState, T3Representation, _finite_real,
    _finite_complex, _positive_scale, _require_optional_real, _require_optional_complex)
from Numerical.core.GeneralCouplingLayout import spec_for

@dataclass(frozen=True)
class T3IntermediateScalarState:
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
    general_couplings: Mapping[str, Any] = field(default_factory=dict)

    def validated(self) -> "T3IntermediateScalarState":
        rep = self.representation.validated()
        if rep.shared_scalar:
            raise ValueError("Use SharedT3IntermediateScalarState for shared scalars")
        is_doublet = rep.d_s1 == rep.d_s2 == 2 and rep.alpha == -1
        real_required = {'lambdaH1Adj':rep.d_s1>1,'lambdaH2Adj':rep.d_s2>1,
            'lambdaS1Adj':rep.d_s1==3,'lambdaS2Adj':rep.d_s2==3,
            'lambda12Adj':rep.d_s1>1 and rep.d_s2>1,
            'lambda12Cross':rep.d_s1==3 and rep.d_s2==3}
        complex_required = ('lambdaHHdagS2S2','lambdaHHdagS1barS1bar',
            'lambdaS1bar2S2bar2','lambdaS1barS2S2bar2','lambdaS1S1bar2S2bar',
            'lambdaHHdagS1barS2barCross')
        real_optionals = {k:_require_optional_real(k,getattr(self,k),required=req)
            for k,req in real_required.items()}
        complex_optionals = {k:_require_optional_complex(k,getattr(self,k),required=is_doublet)
            for k in complex_required}
        extras = {}
        for key,value in self.general_couplings.items():
            if not isinstance(key,str) or not (key.startswith('lambdaT3GeneralQ') or
                key in ('muT3GeneralS1Sq','muT3GeneralS2Sq','muT3GeneralS1S2')):
                raise ValueError(f"Not an EFT1 general scalar coupling: {key!r}")
            extras[key] = spec_for(key).validate(value)
        return T3IntermediateScalarState(
            mu_gev=_positive_scale(self.mu_gev), representation=rep, sm=self.sm.validated(),
            mS1Sq=_finite_real('mS1Sq',self.mS1Sq),mS2Sq=_finite_real('mS2Sq',self.mS2Sq),
            lambdaS1=_finite_real('lambdaS1',self.lambdaS1),
            lambdaS2=_finite_real('lambdaS2',self.lambdaS2),
            lambdaH1=_finite_real('lambdaH1',self.lambdaH1),
            lambdaH2=_finite_real('lambdaH2',self.lambdaH2),
            lambda12=_finite_real('lambda12',self.lambda12),
            lambdaT3=_finite_complex('lambdaT3',self.lambdaT3),
            **real_optionals,**complex_optionals,general_couplings=extras)

@dataclass(frozen=True)
class SharedT3IntermediateScalarState:
    mu_gev: float
    representation: T3Representation
    sm: SMNumericalState
    mSSq: float
    lambdaS: float
    lambda3: float
    lambda4: float
    lambda5: float
    def validated(self) -> "SharedT3IntermediateScalarState":
        rep=self.representation.validated()
        if not rep.shared_scalar: raise ValueError('Shared scalar representation required')
        return SharedT3IntermediateScalarState(
            mu_gev=_positive_scale(self.mu_gev),representation=rep,sm=self.sm.validated(),
            mSSq=_finite_real('mSSq',self.mSSq),lambdaS=_finite_real('lambdaS',self.lambdaS),
            lambda3=_finite_real('lambda3',self.lambda3),lambda4=_finite_real('lambda4',self.lambda4),
            lambda5=_finite_real('lambda5',self.lambda5))
T3EFT1State=T3IntermediateScalarState
SharedT3EFT1State=SharedT3IntermediateScalarState
