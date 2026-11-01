"""Existing EFT1 ODE vector layout extended to include general scalar couplings."""
from __future__ import annotations
import numpy as np
from Numerical.core.State import SMNumericalState
from Numerical.core.StateVector import (StateVectorLayout, _VectorBuilder,
    _pack_sm, _read_real_scalar, _read_complex_scalar, _read_complex_matrix, _check_vector)
from Numerical.running.IntermediateScalarState import (
    T3IntermediateScalarState,SharedT3IntermediateScalarState)
IntermediateScalarState=T3IntermediateScalarState|SharedT3IntermediateScalarState
_REAL_OPTIONALS=('lambdaH1Adj','lambdaH2Adj','lambdaS1Adj','lambdaS2Adj','lambda12Adj','lambda12Cross')
_COMPLEX_OPTIONALS=('lambdaHHdagS2S2','lambdaHHdagS1barS1bar','lambdaS1bar2S2bar2',
    'lambdaS1barS2S2bar2','lambdaS1S1bar2S2bar','lambdaHHdagS1barS2barCross')

def pack_intermediate_scalar_state(state:IntermediateScalarState)->tuple[np.ndarray,StateVectorLayout]:
    state=state.validated();b=_VectorBuilder();_pack_sm(b,state.sm)
    if isinstance(state,SharedT3IntermediateScalarState):
        for n in ('mSSq','lambdaS','lambda3','lambda4','lambda5'):
            b.add_real_scalar(n,getattr(state,n))
    else:
        for n in ('mS1Sq','mS2Sq','lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12'):
            b.add_real_scalar(n,getattr(state,n))
        b.add_complex_scalar('lambdaT3',state.lambdaT3)
        for n in _REAL_OPTIONALS:
            v=getattr(state,n)
            if v is not None:b.add_real_scalar(n,v)
        for n in _COMPLEX_OPTIONALS:
            v=getattr(state,n)
            if v is not None:b.add_complex_scalar(n,v)
        for n,v in sorted(state.general_couplings.items()):
            b.add_complex_scalar(n,v)
    return b.finish(representation=state.representation,reference_mu_gev=state.mu_gev,
        shared_scalar=state.representation.shared_scalar)

def unpack_intermediate_scalar_state(vector:np.ndarray,layout:StateVectorLayout,*,mu_gev:float|None=None)->IntermediateScalarState:
    vector=_check_vector(vector,layout); blocks={b.name:b for b in layout.blocks}
    def real(n):return _read_real_scalar(vector,blocks[n])
    def comp(n):return _read_complex_scalar(vector,blocks[n])
    def mat(n):return _read_complex_matrix(vector,blocks[n])
    sm=SMNumericalState(gY=real('gY'),g2=real('g2'),g3=real('g3'),lambdaH=real('lambdaH'),
        yu=mat('yu'),yd=mat('yd'),ye=mat('ye'))
    scale=layout.reference_mu_gev if mu_gev is None else float(mu_gev)
    if layout.shared_scalar:
        return SharedT3IntermediateScalarState(mu_gev=scale,representation=layout.representation,
            sm=sm,**{n:real(n) for n in ('mSSq','lambdaS','lambda3','lambda4','lambda5')}).validated()
    extra={}
    for n in blocks:
        if n.startswith('lambdaT3GeneralQ') or n in ('muT3GeneralS1Sq','muT3GeneralS2Sq','muT3GeneralS1S2'):
            extra[n]=comp(n)
    return T3IntermediateScalarState(mu_gev=scale,representation=layout.representation,sm=sm,
        **{n:real(n) for n in ('mS1Sq','mS2Sq','lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12')},
        lambdaT3=comp('lambdaT3'),
        **{n:real(n) if n in blocks else None for n in _REAL_OPTIONALS},
        **{n:comp(n) if n in blocks else None for n in _COMPLEX_OPTIONALS},
        general_couplings=extra).validated()
EFT1State=IntermediateScalarState
pack_eft1_state=pack_intermediate_scalar_state
unpack_eft1_state=unpack_intermediate_scalar_state
