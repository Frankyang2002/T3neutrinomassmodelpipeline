"""Canonical T3 numerical states; ordinary models also accept RGBeta general couplings."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Mapping
import numpy as np
from common.T3Model import valid_shared_scalar_dimensions, valid_t3_dimensions
from Numerical.core.GeneralCouplingLayout import spec_for

N_SM_FLAVOR = N_HEAVY_FLAVOR = 3

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
    result = _finite_real('mu_gev', value)
    if result <= 0:
        raise ValueError('mu_gev must be positive.')
    return result

def _complex_matrix(name: str, value: Any, shape: tuple[int, int]) -> np.ndarray:
    matrix = np.asarray(value, dtype=complex)
    if matrix.shape != shape or not (np.isfinite(matrix.real).all() and np.isfinite(matrix.imag).all()):
        raise ValueError(f'{name} must contain finite entries with shape {shape}.')
    return matrix.copy()

def _require_optional_real(name: str, value: float | None, *, required: bool) -> float | None:
    if required:
        if value is None:
            raise ValueError(f'{name} is required for this SU(2) representation.')
        return _finite_real(name, value)
    if value is not None:
        raise ValueError(f'{name} is not present for this SU(2) representation.')
    return None

def _require_optional_complex(name: str, value: complex | None, *, required: bool) -> complex | None:
    if required:
        if value is None:
            raise ValueError(f'{name} is required for this SU(2)/hypercharge assignment.')
        return _finite_complex(name, value)
    if value is not None:
        raise ValueError(f'{name} is not present for this SU(2)/hypercharge assignment.')
    return None

@dataclass(frozen=True)
class T3Representation:
    d_s1: int
    d_s2: int
    d_f: int
    alpha: int
    shared_scalar: bool = False

    def validated(self) -> 'T3Representation':
        d1, d2, df, a = map(int, (self.d_s1, self.d_s2, self.d_f, self.alpha))
        if self.shared_scalar:
            if d1 != d2 or a != -1 or not valid_shared_scalar_dimensions(d1, df):
                raise ValueError('Unsupported shared-scalar representation.')
        elif not valid_t3_dimensions(d1, d2, df):
            raise ValueError(f'Unsupported ordinary T3 dimensions (dS1,dS2,dF)=({d1},{d2},{df}).')
        return T3Representation(d1, d2, df, a, bool(self.shared_scalar))

    @property
    def y_s1(self) -> float: return self.alpha / 2.0
    @property
    def y_s2(self) -> float: return (self.alpha + 2) / 2.0
    @property
    def y_f(self) -> float: return (self.alpha + 1) / 2.0
    @property
    def majorana_fermion(self) -> bool: return self.y_f == 0.0

@dataclass(frozen=True)
class SMNumericalState:
    gY: float
    g2: float
    g3: float
    lambdaH: float
    yu: np.ndarray
    yd: np.ndarray
    ye: np.ndarray

    def validated(self) -> 'SMNumericalState':
        return SMNumericalState(
            gY=_finite_real('gY', self.gY), g2=_finite_real('g2', self.g2),
            g3=_finite_real('g3', self.g3), lambdaH=_finite_real('lambdaH', self.lambdaH),
            yu=_complex_matrix('yu', self.yu, (3, 3)),
            yd=_complex_matrix('yd', self.yd, (3, 3)),
            ye=_complex_matrix('ye', self.ye, (3, 3)))

_OPTIONAL_REAL = ('lambdaH1Adj','lambdaH2Adj','lambdaS1Adj','lambdaS2Adj','lambda12Adj','lambda12Cross')
_OPTIONAL_COMPLEX = ('lambdaHHdagS2S2','lambdaHHdagS1barS1bar','lambdaS1bar2S2bar2',
                     'lambdaS1barS2S2bar2','lambdaS1S1bar2S2bar','lambdaHHdagS1barS2barCross')

@dataclass(frozen=True)
class T3UVState:
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
    # Additional couplings are present only when general RGBeta mode is selected.
    # Mapping entries are copied and validated, never silently filtered.
    general_couplings: Mapping[str, Any] = field(default_factory=dict)

    def validated(self) -> 'T3UVState':
        rep = self.representation.validated()
        if rep.shared_scalar:
            raise ValueError('T3UVState is for split scalars; use SharedT3UVState instead.')
        mf = _complex_matrix('MF', self.MF, (3, 3))
        if rep.majorana_fermion and not np.allclose(mf, mf.T, rtol=1e-10, atol=1e-14):
            raise ValueError('MF must be complex symmetric when the heavy fermion is Majorana.')
        special = rep.d_s1 == 2 and rep.d_s2 == 2 and rep.alpha == -1
        rflags = (rep.d_s1 > 1, rep.d_s2 > 1, rep.d_s1 == 3,
                  rep.d_s2 == 3, rep.d_s1 > 1 and rep.d_s2 > 1,
                  rep.d_s1 == 3 and rep.d_s2 == 3)
        optional_real = {name: _require_optional_real(name, getattr(self,name), required=req)
                         for name,req in zip(_OPTIONAL_REAL, rflags)}
        optional_complex = {name: _require_optional_complex(name, getattr(self,name), required=special)
                            for name in _OPTIONAL_COMPLEX}
        base = set(self.__dataclass_fields__) | {'gY','g2','g3','lambdaH','yu','yd','ye'}
        extras = {}
        for name,value in self.general_couplings.items():
            if not isinstance(name,str) or name in base:
                raise ValueError(f'General coupling collides with baseline or has invalid name: {name!r}')
            extras[name] = spec_for(name).validate(value)
        return T3UVState(
            mu_gev=_positive_scale(self.mu_gev), representation=rep, sm=self.sm.validated(),
            y1=_complex_matrix('y1',self.y1,(3,3)), y2=_complex_matrix('y2',self.y2,(3,3)), MF=mf,
            mS1Sq=_finite_real('mS1Sq',self.mS1Sq), mS2Sq=_finite_real('mS2Sq',self.mS2Sq),
            lambdaS1=_finite_real('lambdaS1',self.lambdaS1),lambdaS2=_finite_real('lambdaS2',self.lambdaS2),
            lambdaH1=_finite_real('lambdaH1',self.lambdaH1),lambdaH2=_finite_real('lambdaH2',self.lambdaH2),
            lambda12=_finite_real('lambda12',self.lambda12),lambdaT3=_finite_complex('lambdaT3',self.lambdaT3),
            **optional_real, **optional_complex, general_couplings=extras)

@dataclass(frozen=True)
class SharedT3UVState:
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

    def validated(self) -> 'SharedT3UVState':
        rep = self.representation.validated()
        if not rep.shared_scalar:
            raise ValueError('SharedT3UVState requires representation.shared_scalar=True.')
        mf = _complex_matrix('MF',self.MF,(3,3))
        if not np.allclose(mf,mf.T,rtol=1e-10,atol=1e-14):
            raise ValueError('MF must be complex symmetric in the supported shared-scalar models.')
        return SharedT3UVState(
            mu_gev=_positive_scale(self.mu_gev),representation=rep,sm=self.sm.validated(),
            h=_complex_matrix('h',self.h,(3,3)),MF=mf,mSSq=_finite_real('mSSq',self.mSSq),
            lambdaS=_finite_real('lambdaS',self.lambdaS),lambda3=_finite_real('lambda3',self.lambda3),
            lambda4=_finite_real('lambda4',self.lambda4),lambda5=_finite_real('lambda5',self.lambda5))
