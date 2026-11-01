"""Evaluate standard and general scalar-only EFT1 RGBeta payloads."""
from __future__ import annotations
from typing import Any,Mapping
import numpy as np
from Numerical.core.BetaVector import beta_values_to_derivative,validate_beta_keys
from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
from Numerical.running.IntermediateScalarState import (
    SharedT3IntermediateScalarState,T3IntermediateScalarState)
from Numerical.running.IntermediateScalarStateVector import (
    IntermediateScalarState,pack_intermediate_scalar_state)
EXPECTED_CONVENTION='16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'

def intermediate_scalar_state_environment(state:IntermediateScalarState)->dict[str,Any]:
    state=state.validated()
    env={'gY':state.sm.gY,'g2':state.sm.g2,'g3':state.sm.g3,
         'yu':state.sm.yu,'yd':state.sm.yd,'ye':state.sm.ye,'lambdaH':state.sm.lambdaH}
    if isinstance(state,SharedT3IntermediateScalarState):
        env.update({n:getattr(state,n) for n in ('mSSq','lambdaS','lambda3','lambda4','lambda5')})
        return env
    for n in ('mS1Sq','mS2Sq','lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12','lambdaT3',
        'lambdaH1Adj','lambdaH2Adj','lambdaS1Adj','lambdaS2Adj','lambda12Adj','lambda12Cross',
        'lambdaHHdagS2S2','lambdaHHdagS1barS1bar','lambdaS1bar2S2bar2',
        'lambdaS1barS2S2bar2','lambdaS1S1bar2S2bar','lambdaHHdagS1barS2barCross'):
        v=getattr(state,n)
        if v is not None:env[n]=v
    env.update(state.general_couplings)
    return env

def _validate_intermediate_scalar_metadata(payload:Mapping[str,Any],state:IntermediateScalarState)->None:
    meta=payload.get('metadata',{})
    if not isinstance(meta,Mapping):raise ValueError('EFT1 RGBeta metadata must be a mapping')
    rep=state.representation.validated()
    expected={'SharedScalar':rep.shared_scalar,'dS1':rep.d_s1,'dS2':rep.d_s2,
        'dF':rep.d_f,'alpha':rep.alpha,'IntegratedField':'F'}
    actual={**meta,'SharedScalar':bool(meta.get('SharedScalar',False))}
    for k,v in expected.items():
        if actual.get(k)!=v:raise ValueError(f'EFT1 RGBeta {k}: payload={actual.get(k)!r}, state={v!r}')
    fields=['S'] if rep.shared_scalar else ['S1','S2']
    if list(meta.get('ActiveBSMFields',[]))!=fields:
        raise ValueError('Unexpected EFT1 ActiveBSMFields')
    if meta.get('ReportBetaConvention')!=EXPECTED_CONVENTION:
        raise ValueError('Unexpected EFT1 RGBeta convention')

def evaluate_intermediate_scalar_rgbeta_payload(payload:Mapping[str,Any],state:IntermediateScalarState)->dict[str,Any]:
    fmt=payload.get('format'); status=payload.get('status')
    if fmt=='general-t3-quadratic-rgbeta-experimental-v1':
        if status!='SuccessCandidateMassBeta' or payload.get('metadata',{}).get('GeneralScalarStage')!='EFT1':
            raise ValueError('General EFT1 RGBeta payload stage/status mismatch')
        if not isinstance(state,T3IntermediateScalarState):
            raise ValueError('General EFT1 requires ordinary split scalar state')
    elif status!='Success':raise ValueError('EFT1 RGBeta status is not Success')
    state=state.validated();_validate_intermediate_scalar_metadata(payload,state)
    betas=payload.get('report_betas')
    if not isinstance(betas,Mapping):raise ValueError('EFT1 report_betas missing')
    _,layout=pack_intermediate_scalar_state(state)
    validate_beta_keys({n:None for n in betas},layout)
    env=intermediate_scalar_state_environment(state)
    result={str(n):evaluate_inputform_expression(str(expr),env) for n,expr in betas.items()}
    validate_beta_keys(result,layout)
    return result

def intermediate_scalar_derivative_from_payload(payload:Mapping[str,Any],state:IntermediateScalarState):
    state=state.validated();_,layout=pack_intermediate_scalar_state(state)
    result=evaluate_intermediate_scalar_rgbeta_payload(payload,state)
    return beta_values_to_derivative(result,layout,divide_by_loop_factor=True),layout
EFT1State=IntermediateScalarState
eft1_state_environment=intermediate_scalar_state_environment
evaluate_eft1_rgbeta_payload=evaluate_intermediate_scalar_rgbeta_payload
eft1_derivative_from_payload=intermediate_scalar_derivative_from_payload
