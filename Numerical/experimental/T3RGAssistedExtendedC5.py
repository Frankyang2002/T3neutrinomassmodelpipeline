"""Six-model UV-to-extended-C5 numerical bridge with explicit, UNCERTIFIED coupling maps.

The exported RGBeta UV theory is integrated only above all massive thresholds.
Its one-step Matchete C5 is evaluated at the UV endpoint under declared
complex-linear coupling and heavy-mass-basis hypotheses. No EFT1 evolution,
sequential matching, or certified RGBeta/Matchete vertex map is implied.

This is independent of the successful 16-model ODE implementation: it reuses
that solver and the exact source-checked six-model C5 evaluator unchanged.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any, Mapping
import copy
import json
import math
import re

import numpy as np

from Numerical.experimental import T3GeneralYukawaCoupledUV as UV
from Numerical.experimental import T3ExtendedC5ExpressionReduction as C5

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'Reports/output/workstream_a/rg_assisted_extended_c5'
EXTRA_QUARTIC = {
    'A_p0':('lambdaT3ExtraQ021010I1','lambdaT3GeneralQ021010I1'),
    'A_m2':('lambdaT3ExtraQ020101I1','lambdaT3GeneralQ020101I1'),
    'D_p0':('lambdaT3ExtraQ021010I1','lambdaT3GeneralQ021010I1'),
    'D_m2':('lambdaT3ExtraQ020101I1','lambdaT3GeneralQ020101I1'),
    'E_p0':('lambdaT3ExtraQ021010I1','lambdaT3GeneralQ021010I1'),
    'E_m2':('lambdaT3ExtraQ020101I1','lambdaT3GeneralQ020101I1'),
}
SCHEMA='t3-conditional-rgbeta-matchete-bridge-scenario-v1'
STATUS='ConditionalUVEndpointOneStepExtendedC5_NotPhysicalPrediction'


def check(condition: bool, msg: str) -> None:
    if not condition:
        raise ValueError(msg)


def zvalue(v:Any, label:str)->complex:
    check(isinstance(v,dict) and set(v)=={'real','imag'},label+' must have real,imag')
    z=complex(float(v['real']),float(v['imag']))
    check(bool(np.isfinite([z.real,z.imag]).all()),label+' is not finite')
    return z


def encode(z:complex)->dict[str,float]:
    z=complex(z)
    return {'real':float(z.real),'imag':float(z.imag)}


def mat_json(a:np.ndarray)->dict[str,Any]:
    x=np.asarray(a,dtype=complex)
    check(x.shape==(3,3) and bool(np.isfinite(x).all()),'Invalid 3x3 matrix')
    return {'real':x.real.tolist(),'imag':x.imag.tolist()}


def from_json(a:Mapping[str,Any])->np.ndarray:
    check(isinstance(a,Mapping) and set(a)=={'real','imag'},'Missing complex matrix real/imag')
    r=np.asarray(a['real'],float); i=np.asarray(a['imag'],float)
    check(r.shape==(3,3) and i.shape==r.shape and bool(np.isfinite(r).all())
          and bool(np.isfinite(i).all()),'Invalid complex matrix components')
    return r+1j*i


def read_scenario(config:Mapping[str,Any])->dict[str,Any]:
    check(isinstance(config,Mapping) and config.get('format')==SCHEMA,'Wrong bridge scenario format')
    check(config.get('physical_coupling_dictionary_certified') is False
          and config.get('hierarchical_threshold_matching_certified') is False,
          'Scenario MUST explicitly declare both physical gates false')
    check(config.get('basis_hypothesis')=='complex-linear-unverified-source-components',
          'Only declared complex-linear unverified maps are implemented')
    m=config.get('matchete_from_rgbeta')
    required={'y1','y2','extra_yukawa','lambdaT3','extra_quartic'}
    check(isinstance(m,dict) and set(m)==required,'Five conversion multipliers required')
    scales={k:zvalue(v,k) for k,v in m.items()}
    check(all(abs(z)>1e-14 for z in scales.values()),'Conversion cannot have zero magnitude')
    phase=config.get('extra_quartic_relative_phase_radians')
    check(isinstance(phase,list) and 1<=len(phase)<=24,'Require 1-24 phases')
    phases=[float(t) for t in phase]
    check(all(math.isfinite(t) and abs(t)<=2*math.pi for t in phases),'Invalid scan phase')
    buffer=float(config.get('endpoint_buffer',1.5))
    check(math.isfinite(buffer) and 1.05<=buffer<=5.,'Endpoint buffer must be between 1.05 and 5')
    n_save=int(config.get('n_save',7))
    check(2<=n_save<=100,'Save count out of range')
    seed_factor=float(config.get('seed_factor',0.3))
    check(math.isfinite(seed_factor) and 0<=seed_factor<=5,'Invalid seed_factor')
    ordering=config.get('ordering','AUTO')
    check(ordering in ('AUTO','NO','IO'),'Unknown neutrino mass ordering')
    vev=float(config.get('vev_gev',246.22))
    check(math.isfinite(vev) and vev>0,'Invalid Higgs vev')
    return {'scales':scales,'phases':phases,'endpoint_buffer':buffer,
            'n_save':n_save,'seed_factor':seed_factor,'ordering':ordering,'vev_gev':vev,'raw':config}


def check_rgbeta_quartic(model:str, root:Path, payload:Mapping[str,Any])->dict[str,Any]:
    check(model in EXTRA_QUARTIC,'Unsupported model')
    mat_name,rg_name=EXTRA_QUARTIC[model]
    extra=UV.EXTRAS[model]
    betas=payload.get('report_betas',{})
    for key in ('y1','y2',extra,'MF','mS1Sq','mS2Sq','lambdaT3',rg_name):
        check(key in betas,f'{model}: missing required RGBeta beta {key}')
    alpha=UV.SPECS[model][3]
    source=root/'RGE/general/generated'/f'{model[0]}_{alpha:+d}_extra_quartics.wl'
    check(source.is_file(),f'Missing generated RGBeta quartic declaration: {source}')
    content=source.read_text(encoding='utf-8-sig')
    # This proves field-content identification, NOT cross-package tensor orientation.
    body=re.escape('{Bar @ H, Bar @ H, S1, S2}' if alpha==0 else
                   '{Bar @ H, Bar @ H, Bar @ S1, Bar @ S2}')
    pattern=r'AddQuartic\[\s*'+re.escape(rg_name)+r'\s*,\s*'+body+r'\s*,\s*GroupInvariant\s*->'
    check(re.search(pattern,content,re.S) is not None,
          f'{model}: RGBeta extra quartic field content/source changed')
    return {'matchete_name':mat_name,'rgbeta_name':rg_name,
            'rgbeta_extension_source':str(source.relative_to(root)),
            'source_sha256':sha256(source.read_bytes()).hexdigest(),
            'field_monomial':('Hdag Hdag S1 S2' if alpha==0 else 'Hdag Hdag S1dag S2dag'),
            'field_content_same_coupling_not_certified':True}


@dataclass(frozen=True)
class DiracBasis:
    masses:np.ndarray
    left:np.ndarray
    right:np.ndarray
    residual:float


def vectorlike_basis(mf:np.ndarray)->DiracBasis:
    """Biunitary mass decomposition as in existing FinalC5TrajectoryAdapter.

    Assumes y1/y1Conj contract with the same left heavy slot and
    y2/y2Conj with the right heavy slot. This action-level slot map is not
    independently certified here, so all numerical predictions remain conditional.
    """
    mf=UV.matrix(mf)
    u,m,vh=np.linalg.svd(mf)
    order=np.argsort(m)
    left=u[:,order];right=vh.conj().T[:,order];m=m[order]
    check(bool(np.all(m>0)), 'Non-positive heavy Dirac singular value')
    delta=left.conj().T@mf@right-np.diag(m)
    residual=float(np.linalg.norm(delta)/max(float(np.linalg.norm(mf)),1e-300))
    check(residual<=5e-11,'Dirac SVD mass-basis residual too large')
    return DiracBasis(m,left,right,residual)


def to_matchete_inputs(model:str,uv:Mapping[str,Any], scenario:Mapping[str,Any],
                       *,extra_phase:float=0.,include_extra:bool=True)->dict[str,Any]:
    check(model in EXTRA_QUARTIC,'Model not supported by extra-C5 calculation')
    mat_name,rgb_name=EXTRA_QUARTIC[model]
    extra=UV.EXTRAS[model]
    s=scenario['scales']
    basis=vectorlike_basis(uv['MF'])
    # Both C-conjugate Yukawas share the F-side (left) heavy rotation,
    # while the un-conjugated-side Yukawas use the right rotation.
    def heavy_rotate(name:str)->np.ndarray:
        m=UV.matrix(uv[name])
        return m@ (basis.left.conj() if name in ('y1','y1Conj') else basis.right)
    lambdas={
        'lambdaT3':s['lambdaT3']*complex(uv['lambdaT3']),
        mat_name:s['extra_quartic']*np.exp(1j*extra_phase)*complex(uv[rgb_name]),
    }
    p={'MF_gev':basis.masses.tolist(),
       'MS1_gev':math.sqrt(float(uv['mS1Sq'])),
       'MS2_gev':math.sqrt(float(uv['mS2Sq'])),
       'y1':mat_json(s['y1']*heavy_rotate('y1')),
       'y2':mat_json(s['y2']*heavy_rotate('y2')),
       extra:mat_json(s['extra_yukawa']*heavy_rotate(extra) if include_extra
                      else np.zeros((3,3),complex)),
       'lambdaT3':encode(lambdas['lambdaT3']),
       mat_name:encode(lambdas[mat_name]),
       'vev_gev':scenario['vev_gev'],
       'ye':mat_json(UV.matrix(uv['ye'])),
       'ordering':scenario['ordering'],
    }
    check(p['MS1_gev']>0 and p['MS2_gev']>0,'Nonpositive scalar mass')
    result={'model':model,
            'basis_declaration':'unverified-rgbeta-to-matchete-transport-hypothesis',
            'parameters':p}
    return result


def evaluate_hypothetically_transported_c5(model:str, input_data:Mapping[str,Any])->dict[str,Any]:
    """Evaluate verified one-step kernel WITHOUT bypassing the existing gate.

    The old Matchete-only evaluator intentionally refuses RGBeta transport;
    this separate evaluator takes the unverified provenance explicitly.
    It uses the exact kernel/formula established by C5.verify_model.
    """
    check(input_data.get('model')==model and
          input_data.get('basis_declaration')=='unverified-rgbeta-to-matchete-transport-hypothesis',
          'Hypothetical transport provenance missing')
    p=input_data['parameters'];extra=C5.MODELS[model][0];quartic=C5.MODELS[model][1]
    masses=np.asarray(p['MF_gev'],float)
    check(masses.shape==(3,) and bool(np.isfinite(masses).all()) and
          bool(np.all(masses>0)),'Invalid mass eigenvalues')
    m1=float(p['MS1_gev']);m2=float(p['MS2_gev'])
    check(min(m1,m2)>0 and bool(np.isfinite([m1,m2]).all()),'Invalid scalar masses')
    loops=np.array([m*C5.three_mass_loop(m*m,m1*m1,m2*m2) for m in masses])
    coeff=complex(C5.MODEL_COEFFICIENT[model].evalf())/(16*math.pi**2)
    y1=from_json(p['y1']);y2=from_json(p['y2']);ye=from_json(p['ye']);yx=from_json(p[extra])
    lam=zvalue(p['lambdaT3'],'lambdaT3');lex=zvalue(p[quartic],quartic)
    def pair(a:np.ndarray,b:np.ndarray,scalar:complex)->np.ndarray:
        z=coeff*scalar*(np.conjugate(a)*loops[np.newaxis,:])@np.conjugate(b).T
        return z+z.T
    base=pair(y1,y2,lam)
    direct=pair(yx if extra=='y1Conj' else y1,
                y2 if extra=='y1Conj' else yx,np.conjugate(lex))
    total=base+direct
    check(bool(np.isfinite(total).all()) and
          bool(np.allclose(total,total.T,rtol=1e-12,atol=1e-30)),
          'Extended one-step C5 failed symmetry/finite guard')
    from physics.NeutrinoTrajectory import charged_lepton_mass_basis_matrix
    from physics.NeutrinoObservables import calculate_neutrino_observables
    v=float(p['vev_gev'])
    mnu=-.5*v*v*total
    obs=calculate_neutrino_observables(charged_lepton_mass_basis_matrix(mnu,ye),
                                      ordering=p['ordering'])
    check(float(obs.takagi_residual)<1e-6,'Takagi factorization failed')
    return {
        'status':'HypotheticallyTransportedOneStepC5_NotPhysicalPrediction',
        'basis_declaration':input_data['basis_declaration'],
        'C5_baseline_gev_inverse':mat_json(base),
        'C5_extra_gev_inverse':mat_json(direct),
        'C5_total_gev_inverse':mat_json(total),
        'conditional_observables':{
            'masses_ev':obs.masses_ev.tolist(),
            'delta_m21_sq_ev2':float(obs.delta_m21_sq_ev2),
            'delta_m31_sq_ev2':float(obs.delta_m31_sq_ev2),
            'pmns_abs':obs.pmns_abs.tolist(),
            'takagi_residual':float(obs.takagi_residual),
            'ordering':str(obs.ordering),
        },
        'absolute_coupling_dictionary_certified':False,
        'physical_neutrino_mass_prediction_certified':False,
    }


def _matrix_from_result(c5:Mapping[str,Any],name:str)->np.ndarray:
    return from_json(c5[name])


def _obs(c5:Mapping[str,Any])->dict[str,Any]:
    # The conditional evaluator computes observables only after Ye is provided.
    return c5['conditional_observables']


def validate_endpoint(model:str,env:Mapping[str,Any],mu:float)->dict[str,float]:
    heavy=float(max(np.linalg.svd(UV.matrix(env['MF']),compute_uv=False)))
    a=float(env['mS1Sq']);b=float(env['mS2Sq'])
    check(a>0 and b>0,'Tachyonic/zero bare scalar square mass')
    ms1,ms2=math.sqrt(a),math.sqrt(b)
    largest=max(heavy,ms1,ms2)
    check(mu>1.02*largest,f'{model}: endpoint {mu:g} GeV already crosses heavy mass threshold {largest:g} GeV')
    return {'max_fermion_singular_gev':heavy,'scalar1_bare_gev':ms1,
            'scalar2_bare_gev':ms2,'max_heavy_mass_gev':largest,
            'first_threshold_is_not_crossed':True}


def compare_solution(model:str,endpoint:Mapping[str,Any],control:Mapping[str,Any],
                     scenario:Mapping[str,Any],mu:float)->dict[str,Any]:
    validate_endpoint(model,endpoint,mu)
    validate_endpoint(model,control,mu)
    results=[]
    for phase in scenario['phases']:
        a=to_matchete_inputs(model,endpoint,scenario,extra_phase=phase)
        b=to_matchete_inputs(model,control,scenario,extra_phase=phase)
        full=evaluate_hypothetically_transported_c5(model,a)
        baseline=evaluate_hypothetically_transported_c5(model,b)
        ct=_matrix_from_result(full,'C5_total_gev_inverse')
        c0=_matrix_from_result(baseline,'C5_total_gev_inverse')
        direct=_matrix_from_result(full,'C5_extra_gev_inverse')
        total_change=ct-c0
        check(bool(np.isfinite(ct).all()) and bool(np.isfinite(total_change).all()),'Nonfinite C5')
        # Baseline/extra separation must reproduce the sum numerically.
        check(np.allclose(ct,_matrix_from_result(full,'C5_baseline_gev_inverse')+direct,
                          rtol=2e-12,atol=1e-32),'C5 decomposition broken')
        results.append({
            'extra_relative_phase_radians':phase,
            'C5_total_gev_inverse':mat_json(ct),
            'C5_direct_extra_gev_inverse':mat_json(direct),
            'C5_seeded_minus_zero_seed_gev_inverse':mat_json(total_change),
            'frobenius_C5_total_gev_inverse':float(np.linalg.norm(ct)),
            'frobenius_C5_extra_gev_inverse':float(np.linalg.norm(direct)),
            'frobenius_seeded_minus_zero_seed_gev_inverse':float(np.linalg.norm(total_change)),
            'fractional_extra_direct_to_total':float(np.linalg.norm(direct)/max(float(np.linalg.norm(ct)),1e-300)),
            'conditional_observables':_obs(full),
            'zero_seed_conditional_observables':_obs(baseline),
        })
    mass_points=np.array([r['conditional_observables']['masses_ev'] for r in results],float)
    return {'scenarios':results,'masses_ev_min_by_order':mass_points.min(axis=0).tolist(),
            'masses_ev_max_by_order':mass_points.max(axis=0).tolist(),
            'n_phases':len(results),
            'all_predictions_are_hypothesis_dependent':True}


def run_model(model:str, scenario:Mapping[str,Any], *, root:Path=ROOT,
              native_source:Path|None=None)->dict[str,Any]:
    check(model in EXTRA_QUARTIC,'Unsupported model')
    source=Path(native_source) if native_source else root/'Reports/output/workstream_a/extended_c5_native/native'
    native=C5.verify_model(source,model)
    check(native['exact_reduction_of_all_three_sectors_after_identification'] is True,
          'Native one-step coefficient not reduced')
    payload,config=UV.load_case(model,root/'Reports/output/workstream_a/integrated_general_yukawa_uv/native',
                                root/'Numerical/benchmarks/fitted_reference')
    mapping=check_rgbeta_quartic(model,root,payload)
    settings=scenario
    initial=UV.prepare_state(model,payload,config,settings['seed_factor'])
    uv_mu,endpoint_mu=UV.endpoint_for(config,initial,factor=settings['endpoint_buffer'])
    # Avoid a retrospective lower-scale prediction: the endpoint must be UV-only.
    validate_endpoint(model,initial,uv_mu)
    layout,solution,_=UV.integrate(payload,initial,uv_mu,endpoint_mu,n_save=settings['n_save'])
    final=layout.unpack(solution.y[:,-1])
    gate=validate_endpoint(model,final,endpoint_mu)
    extra=UV.EXTRAS[model]
    zero={key:(value.copy() if isinstance(value,np.ndarray) else value)
          for key,value in initial.items()}
    zero[extra]=np.zeros((3,3),complex)
    _,control,_=UV.integrate(payload,zero,uv_mu,endpoint_mu,n_save=settings['n_save'])
    ctrl=layout.unpack(control.y[:,-1])
    compared=compare_solution(model,final,ctrl,settings,endpoint_mu)
    report={
        'status':STATUS,'model':model,'extra_yukawa':extra,
        'uv_start_gev':uv_mu,'uv_endpoint_gev':endpoint_mu,
        'uv_rhs_evaluations':int(solution.nfev),'uv_control_rhs_evaluations':int(control.nfev),
        'uv_full_n_coupling_blocks':len(layout.blocks),
        'one_step_native_evidence':{'status':native['status'],
                                    'coefficient_source_sha256':native['source_sha256']},
        'source_quartic_field_content':mapping,
        'threshold_ordering':gate,
        'coupling_map_hypotheses':{k:encode(v) for k,v in settings['scales'].items()},
        'heavy_slot_rotation_hypothesis':'Dirac biunitary SVD, y1/y1Conj via conjugate left rotation, y2/y2Conj via right rotation',
        'complex_linear_source_map_verified':False,
        'absolute_physical_matchete_rgbeta_dictionary_certified':False,
        'scalar_bare_mass_is_physical_eigenmass_certified':False,
        'sequential_threshold_matching_performed':False,
        'EFT1_running_performed':False,
        'one_step_matchete_expression_evaluated_with_hypothetical_transport':True,
        'physical_neutrino_mass_prediction_certified':False,
        'status_explanation':'UV RGE endpoint is above all heavy masses; one-step Matchete C5 formula evaluated using declared unverified coupling conversions and bare scalar masses; neither sequential matching nor low energy running is included.',
        'phase_scan':compared,
    }
    return report


def self_test()->None:
    check(len(EXTRA_QUARTIC)==6 and set(EXTRA_QUARTIC)==set(C5.MODELS)==set(UV.EXTRAS),
          'Six-model inventory mismatched')
    a=np.array([[4,1,0],[0,7,0],[0,0,12]],complex)
    bb=vectorlike_basis(a)
    check(abs(bb.residual)<1e-12 and len(bb.masses)==3,'SVD mass basis failed')
    for name in EXTRA_QUARTIC:
        check(EXTRA_QUARTIC[name][0]==C5.MODELS[name][1],
              'Matchete quartic/source model mismatch')
        check(EXTRA_QUARTIC[name][1].replace('General','Extra')==EXTRA_QUARTIC[name][0],
              'RG quartic suffix mismatch')
    payload={'report_betas':{key:'x' for key in ('y1','y2','y1Conj','MF','mS1Sq','mS2Sq','lambdaT3',
                                               'lambdaT3GeneralQ021010I1')}}
    # Self-test doesn't require native Wolfram package: sources gate tested separately.
    try:validate_endpoint('A_p0',{'MF':np.eye(3)*10,'mS1Sq':100.,'mS2Sq':16.},10.)
    except ValueError:pass
    else:raise AssertionError('Threshold crossing was accepted')
    z=1.2+0.5j; check(abs(zvalue(encode(z),'test')-z)<1e-14,'Complex scalar JSON broken')
    try:read_scenario({'format':SCHEMA,'physical_coupling_dictionary_certified':True,
                       'hierarchical_threshold_matching_certified':False})
    except ValueError:pass
    else:raise AssertionError('Physical overclaim accepted')
    print('PASS: six native-to-RGBeta names, Dirac SVD, threshold gate, complex JSON and no-certification controls')
