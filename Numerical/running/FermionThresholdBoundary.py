"""Project existing ordinary/shared T3 UV states into scalar-only EFT1 at F threshold."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from Numerical.core.State import T3UVState,SharedT3UVState
from Numerical.running.IntermediateScalarState import T3IntermediateScalarState,SharedT3IntermediateScalarState
UVState=T3UVState|SharedT3UVState
IntermediateScalarState=T3IntermediateScalarState|SharedT3IntermediateScalarState
@dataclass(frozen=True)
class FermionThresholdDiagnostic:
    mu_gev:float
    masses_gev:np.ndarray
    @property
    def minimum_mass_gev(self)->float:return float(np.min(self.masses_gev))
    @property
    def maximum_mass_gev(self)->float:return float(np.max(self.masses_gev))
def fermion_threshold_masses(state:UVState)->FermionThresholdDiagnostic:
    state=state.validated();m=np.sort(np.asarray(np.linalg.svd(np.asarray(state.MF,dtype=complex),compute_uv=False),dtype=float))
    if not np.all(np.isfinite(m)) or np.any(m<0):raise ValueError('Invalid fermion singular masses')
    return FermionThresholdDiagnostic(float(state.mu_gev),m)
def build_intermediate_scalar_boundary(uv_state:UVState,*,matching_scale_gev:float|None=None)->IntermediateScalarState:
    uv_state=uv_state.validated()
    scale=float(uv_state.mu_gev if matching_scale_gev is None else matching_scale_gev)
    if not np.isfinite(scale) or scale<=0:raise ValueError('matching_scale_gev must be finite and positive')
    if not np.isclose(scale,uv_state.mu_gev,rtol=1e-12,atol=0):
        raise ValueError('UV state must be evolved to matching scale before projecting EFT1')
    if isinstance(uv_state,SharedT3UVState):
        return SharedT3IntermediateScalarState(mu_gev=scale,representation=uv_state.representation,
            sm=uv_state.sm,**{n:getattr(uv_state,n) for n in ('mSSq','lambdaS','lambda3','lambda4','lambda5')}).validated()
    optional=('lambdaH1Adj','lambdaH2Adj','lambdaS1Adj','lambdaS2Adj','lambda12Adj','lambda12Cross',
        'lambdaHHdagS2S2','lambdaHHdagS1barS1bar','lambdaS1bar2S2bar2',
        'lambdaS1barS2S2bar2','lambdaS1S1bar2S2bar','lambdaHHdagS1barS2barCross')
    general={k:v for k,v in uv_state.general_couplings.items()
        if k.startswith('lambdaT3GeneralQ') or k in
        ('muT3GeneralS1Sq','muT3GeneralS2Sq','muT3GeneralS1S2')}
    return T3IntermediateScalarState(mu_gev=scale,representation=uv_state.representation,sm=uv_state.sm,
        **{n:getattr(uv_state,n) for n in ('mS1Sq','mS2Sq','lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12','lambdaT3')},
        **{n:getattr(uv_state,n) for n in optional},general_couplings=general).validated()
fermion_singular_masses=fermion_threshold_masses
project_after_fermion_threshold=build_intermediate_scalar_boundary
EFT1State=IntermediateScalarState
