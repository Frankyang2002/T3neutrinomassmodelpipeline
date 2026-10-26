"""Generic model-dependent RGBeta state extension for general scalar quartics.

This is experimental plumbing. It requires a complete RG beta export and does
not identify RGBeta couplings with Matchete's potential coefficients. The eight
per-model quadratic mixing structures have not been inserted in mass matrices.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Mapping
import json

import numpy as np
from scipy.integrate import solve_ivp

from Numerical.core.BetaVector import beta_values_to_derivative
from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression, state_environment
from Numerical.core.State import T3UVState
from Numerical.core.StateVector import StateVectorLayout, VectorBlock, pack_uv_state, unpack_uv_state
from Numerical.running.IntermediateScalarRGBetaEvaluator import intermediate_scalar_state_environment
from Numerical.running.IntermediateScalarState import T3IntermediateScalarState
from Numerical.running.IntermediateScalarStateVector import (
    pack_intermediate_scalar_state, unpack_intermediate_scalar_state,
)

Status = 'CandidateGeneralScalarRGBetaOnly'
_EXPECTED_CONVENTION = '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'


def _model_name(rep: Any) -> str:
    models={(1,3,2):'A',(2,2,1):'B',(2,2,3):'C',(3,1,2):'D',(3,3,2):'E'}
    key=(int(rep.d_s1),int(rep.d_s2),int(rep.d_f))
    if key not in models:
        raise ValueError(f'Not one of the certified 16 ordinary models: {key}')
    a=int(rep.alpha)
    return f'{models[key]}_{a:+d}'


def registry_model(registry:Mapping[str,Any], rep:Any) -> dict[str,Any]:
    name=_model_name(rep)
    model=registry.get('models',{}).get(name)
    if not isinstance(model,dict):
        raise ValueError(f'Model {name!r} missing from general RGBeta registry.')
    if list(model.get('dimensions',[]))!=[rep.d_s1,rep.d_s2,rep.d_f] or model.get('alpha')!=rep.alpha:
        raise ValueError('General RGBeta registry representation mismatch.')
    if model.get('quadratics_not_implemented',0) > 0:
        # Note this in reports; RG state still safely represents only quartics.
        pass
    if len(model['new_rgbeta_symbols'])+len(model['preexisting_rgbeta_symbols'])!=model['total_additional_quartics']:
        raise ValueError('General scalar quartic inventory is incomplete.')
    return model


@dataclass(frozen=True)
class GeneralT3ScalarRGState:
    baseline: T3UVState | T3IntermediateScalarState
    extra: Mapping[str,complex]

    def validated(self, registry:Mapping[str,Any]) -> 'GeneralT3ScalarRGState':
        base=self.baseline.validated()
        model=registry_model(registry,base.representation)
        expected={item['name'] for item in model['new_rgbeta_symbols']}
        received=set(self.extra)
        if expected!=received:
            raise ValueError(f'Additional RGBeta state couplings differ: missing={sorted(expected-received)} extra={sorted(received-expected)}')
        for item in model['preexisting_rgbeta_symbols']:
            attr=item['existing_name']
            if getattr(base,attr,None) is None:
                raise ValueError(f'Existing RGBeta coupling {attr} missing in baseline state.')
        values={}
        for name in sorted(expected):
            val=complex(self.extra[name])
            if not (np.isfinite(val.real) and np.isfinite(val.imag)):
                raise ValueError(f'{name}: nonfinite complex coupling')
            values[name]=val
        return GeneralT3ScalarRGState(base,values)


def _pack_baseline(baseline):
    if isinstance(baseline,T3UVState):
        return pack_uv_state(baseline)
    if isinstance(baseline,T3IntermediateScalarState):
        return pack_intermediate_scalar_state(baseline)
    raise TypeError('Only ordinary T3 UV/EFT1 states are supported.')


def pack_general_scalar_state(state:GeneralT3ScalarRGState, registry:Mapping[str,Any]):
    state=state.validated(registry)
    vector,base_layout=_pack_baseline(state.baseline)
    blocks=list(base_layout.blocks)
    values=[vector]
    end=base_layout.size
    for name in sorted(state.extra):
        z=state.extra[name]
        blocks.append(VectorBlock(name=name,kind='complex_scalar',start=end,stop=end+2,shape=()))
        values.append(np.array([z.real,z.imag],dtype=float))
        end+=2
    layout=StateVectorLayout(representation=base_layout.representation,
                             reference_mu_gev=base_layout.reference_mu_gev,
                             shared_scalar=False,blocks=tuple(blocks))
    return np.concatenate(values),layout


def unpack_general_scalar_state(vector,layout,registry,*,uv:bool,mu_gev:float):
    v=np.asarray(vector,dtype=float)
    if v.shape!=(layout.size,) or not np.all(np.isfinite(v)):
        raise ValueError('Invalid general scalar ODE vector.')
    model=registry_model(registry,layout.representation)
    extra_names=sorted(x['name'] for x in model['new_rgbeta_symbols'])
    end=layout.size-2*len(extra_names)
    base_blocks=tuple(b for b in layout.blocks if b.stop<=end)
    base_layout=replace(layout,blocks=base_blocks)
    baseline=(unpack_uv_state(v[:end],base_layout,mu_gev=mu_gev) if uv else
              unpack_intermediate_scalar_state(v[:end],base_layout,mu_gev=mu_gev))
    extra={name:complex(float(v[end+2*i]),float(v[end+2*i+1]))
           for i,name in enumerate(extra_names)}
    return GeneralT3ScalarRGState(baseline,extra).validated(registry)


def evaluate_general_scalar_derivative(payload, state, registry):
    """Evaluate all declared betas, fail closed for missing/extra keys."""
    state=state.validated(registry)
    meta=payload.get('metadata',{})
    if payload.get('status')!='Success' or not isinstance(meta,dict):
        raise ValueError('Invalid general RGBeta payload status/metadata.')
    rep=state.baseline.representation
    model=_model_name(rep)
    stage='UV' if isinstance(state.baseline,T3UVState) else 'EFT1'
    expected_meta={'GeneralScalarModel':model,'GeneralScalarStage':stage,
                   'GeneralScalarBasisStatus':'CandidateNativeRGBetaUnmappedToMatchete',
                   'dS1':rep.d_s1,'dS2':rep.d_s2,'dF':rep.d_f,'alpha':rep.alpha,
                   'ReportBetaConvention':_EXPECTED_CONVENTION}
    for k,v in expected_meta.items():
        if meta.get(k)!=v:
            raise ValueError(f'General RGBeta metadata mismatch {k}: expected {v!r}, got {meta.get(k)!r}')
    expected_extras=sorted(state.extra)
    if sorted(meta.get('GeneralScalarNewBetaKeys',[]))!=expected_extras:
        raise ValueError('General RGBeta additional beta keys mismatch registry.')
    betas=payload.get('report_betas')
    if not isinstance(betas,dict):
        raise ValueError('Missing report_betas dictionary.')
    env=(state_environment(state.baseline) if stage=='UV' else
         intermediate_scalar_state_environment(state.baseline))
    env.update(state.extra)
    values={name:evaluate_inputform_expression(expr,env) for name,expr in betas.items()}
    _,layout=pack_general_scalar_state(state,registry)
    dy=beta_values_to_derivative(values,layout,divide_by_loop_factor=True)
    return dy,layout


def run_general_scalar_segment(initial:GeneralT3ScalarRGState,payload,registry,mu_final_gev:float,
                               *,rtol=1e-8,atol=1e-11):
    """Generic experimental UV or scalar-EFT solve_ivp evolution."""
    initial=initial.validated(registry)
    y0,layout=pack_general_scalar_state(initial,registry)
    uv=isinstance(initial.baseline,T3UVState)
    mu0=float(initial.baseline.mu_gev)
    mu1=float(mu_final_gev)
    if not np.isfinite(mu1) or mu1<=0 or mu1==mu0:
        raise ValueError('Final RG scale must be positive and different.')
    test,other_layout=evaluate_general_scalar_derivative(payload,initial,registry)
    if layout!=other_layout or test.shape!=y0.shape:
        raise RuntimeError('General scalar derivative/layout mismatch.')
    def rhs(t,y):
        yy=np.array(y,dtype=float,copy=True)
        if uv and layout.representation.majorana_fermion:
            block=layout.block('MF')
            n=int(np.prod(block.shape))
            mreal=yy[block.start:block.start+n].reshape(block.shape)
            mimag=yy[block.start+n:block.stop].reshape(block.shape)
            mreal[:]=0.5*(mreal+mreal.T)
            mimag[:]=0.5*(mimag+mimag.T)
        state=unpack_general_scalar_state(yy,layout,registry,uv=uv,mu_gev=float(np.exp(t)))
        return evaluate_general_scalar_derivative(payload,state,registry)[0]
    sol=solve_ivp(rhs,(np.log(mu0),np.log(mu1)),y0,method='DOP853',rtol=rtol,atol=atol)
    if not sol.success:raise RuntimeError('General scalar RGE integration failed: '+str(sol.message))
    final_vector=np.array(sol.y[:,-1],copy=True)
    if uv and layout.representation.majorana_fermion:
        block=layout.block('MF')
        n=int(np.prod(block.shape))
        mr=final_vector[block.start:block.start+n].reshape(block.shape)
        mi=final_vector[block.start+n:block.stop].reshape(block.shape)
        mr[:]=0.5*(mr+mr.T)
        mi[:]=0.5*(mi+mi.T)
    final=unpack_general_scalar_state(final_vector,layout,registry,uv=uv,mu_gev=mu1)
    return {'status':Status,'solver_success':True,'mu_initial_gev':mu0,
            'mu_final_gev':mu1,'nfev':int(sol.nfev), 'state':final,
            'physical_neutrino_mass_prediction':False,
            'quartic_matchete_mapping_verified':False}


def load_registry(path:Path|str):
    report=json.loads(Path(path).read_text(encoding='utf-8'))
    if report.get('status')!='GeneratedCandidateBasisNotWolframVerified' or report['totals']['matchete_extra_quartics']!=50:
        raise ValueError('Unexpected general RGBeta registry status/completeness.')
    return report


def project_general_after_fermion(state:GeneralT3ScalarRGState,registry):
    """Tree-level projection F -> scalar EFT, retaining *all* extra quartics."""
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    st=state.validated(registry)
    if not isinstance(st.baseline,T3UVState):
        raise TypeError('Fermion projection requires a UV state.')
    eft=build_intermediate_scalar_boundary(st.baseline,matching_scale_gev=st.baseline.mu_gev)
    return GeneralT3ScalarRGState(eft,dict(st.extra)).validated(registry)


def run_general_scalar_two_segments(initial, uv_payload, eft1_payload, registry,
                                    *,mu_fermion_gev:float,mu_scalar_gev:float):
    """Renormalisable UV + EFT1 flow only; no full general-scalar matching."""
    from Numerical.running.FermionThresholdBoundary import fermion_threshold_masses
    from Numerical.running.ScalarThresholdBoundary import scalar_threshold_masses,build_final_sm_boundary
    initial=initial.validated(registry)
    if not isinstance(initial.baseline,T3UVState):
        raise TypeError('Two-stage flow requires a UV state.')
    mu0=float(initial.baseline.mu_gev)
    mu_f=float(mu_fermion_gev)
    mu_s=float(mu_scalar_gev)
    if not (np.isfinite(mu_f) and np.isfinite(mu_s) and mu0>mu_f>mu_s>0):
        raise ValueError('Require mu_UV > mu_F > mu_S > 0')
    uv=run_general_scalar_segment(initial,uv_payload,registry,mu_f)
    at_f=uv['state']
    masses_f=fermion_threshold_masses(at_f.baseline)
    if min(masses_f.masses_gev)<=mu_f:
        raise ValueError('Grouped F-first benchmark invalid: active fermion lighter than mu_F.')
    after_f=project_general_after_fermion(at_f,registry)
    if dict(at_f.extra)!=dict(after_f.extra):
        raise RuntimeError('Extra quartic discontinuity at tree-level fermion threshold')
    eft=run_general_scalar_segment(after_f,eft1_payload,registry,mu_s)
    at_s=eft['state']
    masses_s=scalar_threshold_masses(at_s.baseline)
    if max(masses_s.masses_gev)>=mu_f or min(masses_s.masses_gev)<=mu_s:
        raise ValueError('Grouped scalar threshold ordering invalid for current running masses.')
    sm=build_final_sm_boundary(at_s.baseline,matching_scale_gev=mu_s)
    model=registry_model(registry,initial.baseline.representation)
    return {'status':'GeneralScalarRenormalisableTwoSegmentCandidate',
            'model':model['model'],'uv':uv,'eft1':eft,
            'final_sm_boundary':sm,'uv_nfev':uv['nfev'],'eft1_nfev':eft['nfev'],
            'mass_mixing_quadratics_not_included':model['quadratics_not_implemented'],
            'all_extra_quartics_transported':True,
            'full_general_scalar_neutrino_mass_prediction':False,
            'general_scalar_wilson_matching_complete':False}
