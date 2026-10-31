"""Conditional six-model one-step matching + SM/Weinberg low-energy study.

This is a deliberately explicit **sudden-decoupling approximation**:
  UV RGBeta -> matching at mu_match ABOVE all heavy pole-mass proxies
  -> integrate only SM+Weinberg below mu_match.

It is NOT the physical hierarchical EFT1 scalar-matching calculation.
All cross-package coupling conversions remain scenario hypotheses.

Reuses the already verified native C5 reduction, UV ODE and SM+Weinberg
solver. Does not edit their production implementations.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any, Callable, Mapping
import json
import math

import numpy as np

from Numerical.experimental import T3RGAssistedExtendedC5 as BR
from Numerical.experimental import T3GeneralYukawaCoupledUV as UV
from Numerical.experimental import T3ExtendedC5ExpressionReduction as C5

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'Reports/output/workstream_a/extended_c5_low_energy'
STATUS = 'ConditionalSuddenDecouplingSMWeinbergCompleted_NotHierarchicalPrediction'
MODELS = tuple(BR.EXTRA_QUARTIC)


def need(cond: bool, msg: str) -> None:
    if not cond:
        raise ValueError(msg)


def finite(x: Any, label: str, positive: bool = False) -> float:
    v = float(x)
    need(math.isfinite(v) and (v > 0 if positive else True), f'{label} must be finite'+(' and positive' if positive else ''))
    return v


def load_study(obj: Mapping[str, Any]) -> dict[str, Any]:
    need(isinstance(obj, Mapping) and obj.get('format')=='t3-extended-c5-conditional-low-energy-v1', 'Incorrect study format')
    need(obj.get('physical_dictionary_certified') is False and obj.get('hierarchical_threshold_matching_certified') is False,
         'Both physical-certification flags must be explicitly false')
    vals = obj.get('endpoint_buffers')
    need(isinstance(vals, list) and 1<=len(vals)<=6, 'Require 1-6 matching scale buffers')
    buffers = [finite(x, 'endpoint_buffer', positive=True) for x in vals]
    need(all(1.05<=x<=5.0 for x in buffers) and len(set(buffers))==len(buffers), 'Invalid or duplicate endpoint buffers')
    low = finite(obj.get('mu_low_gev'), 'mu_low_gev', positive=True)
    need(1.0<=low<=1000.0, 'Low scale out of study range')
    rtol=finite(obj.get('rtol',3e-8), 'rtol', positive=True)
    atol=finite(obj.get('atol',1e-18), 'atol', positive=True)
    need(rtol<=1e-4 and atol<=1e-10, 'Low-energy tolerances too loose')
    return {'buffers':buffers,'mu_low_gev':low,'rtol':rtol,'atol':atol,
            'allow_unverified_sudden_decoupling':obj.get('allow_unverified_sudden_decoupling') is True,
            'raw':dict(obj)}


def complex_matrix_json(m: np.ndarray) -> dict[str, Any]:
    return BR.mat_json(np.asarray(m,complex))


def observe(c5: np.ndarray, ye: np.ndarray, *, vev_gev:float, ordering:str) -> dict[str, Any]:
    from physics.NeutrinoTrajectory import charged_lepton_mass_basis_matrix
    from physics.NeutrinoObservables import calculate_neutrino_observables
    mass = -.5*vev_gev**2*np.asarray(c5,complex)
    obs = calculate_neutrino_observables(charged_lepton_mass_basis_matrix(mass,ye),ordering=ordering)
    need(obs.takagi_residual<1e-6,'Bad low-energy Takagi residual')
    return {'masses_ev':obs.masses_ev.tolist(), 'delta_m21_sq_ev2':float(obs.delta_m21_sq_ev2),
            'delta_m31_sq_ev2':float(obs.delta_m31_sq_ev2),
            'delta_m32_sq_ev2':float(getattr(obs,'delta_m32_sq_ev2',obs.delta_m31_sq_ev2-obs.delta_m21_sq_ev2)),
            'ordering':str(obs.ordering),'pmns_abs':obs.pmns_abs.tolist(),
            'takagi_residual':float(obs.takagi_residual),
            'mnu_matrix_norm_ev':float(np.linalg.norm(mass)*1e9)}


def threshold_diagnostic(env: Mapping[str, Any], mu: float) -> dict[str, Any]:
    masses=np.linalg.svd(UV.matrix(env['MF']),compute_uv=False)
    s1=math.sqrt(finite(env['mS1Sq'],'mS1Sq',positive=True));s2=math.sqrt(finite(env['mS2Sq'],'mS2Sq',positive=True))
    heavy=sorted([float(x) for x in masses]+[s1,s2],reverse=True)
    need(mu>1.02*heavy[0],'Matching endpoint crosses or touches a heavy threshold')
    # F-first hierarchical path only applies if all heavy fermions > both scalars.
    f_first=bool(float(min(masses))>max(s1,s2))
    return {'mu_matching_gev':float(mu),'heavy_mass_proxies_gev_descending':heavy,
            'fermion_mass_singular_values_gev':sorted(map(float,masses)),
            'scalar_bare_mass_proxies_gev':[s1,s2],
            'mu_over_largest_mass':float(mu/heavy[0]),
            'max_to_min_heavy_mass_ratio':float(heavy[0]/heavy[-1]),
            'fermion_first_hierarchy_in_mass_proxies':f_first,
            'mass_proxies_not_pole_masses':True,
            'separate_heavy_thresholds_integrated_out':False,
            'intermediate_scalar_eft_run':False}


def _sm_initial(env: Mapping[str, Any], k: np.ndarray):
    from Numerical.running.SMWeinbergEvolution import SMWeinbergInitialConditions
    return SMWeinbergInitialConditions(
        gY=float(env['gY']),g2=float(env['g2']),g3=float(env['g3']),
        lambdaH=float(env['lambdaH']),
        ye=UV.matrix(env['ye']),yu=UV.matrix(env['yu']),yd=UV.matrix(env['yd']),
        K=np.asarray(k,dtype=complex)).validated()


def transport(c5: np.ndarray, env: Mapping[str, Any], mu_match:float, low:float,
              rtol:float,atol:float, *, evolve: Callable|None=None):
    if evolve is None:
        from Numerical.running.SMWeinbergEvolution import evolve_sm_weinberg
        evolve=evolve_sm_weinberg
    need(mu_match>low and low>0,'SM/Weinberg segment requires mu_match > mu_low > 0')
    k=np.asarray(c5,complex)
    need(k.shape==(3,3) and np.isfinite(k).all(),'Bad C5 shape/values')
    need(np.allclose(k,k.T,rtol=1e-10,atol=1e-29),'Non-symmetric matched C5')
    result=evolve(_sm_initial(env,k),mu_match,low,rtol=rtol,atol=atol)
    need(bool(result.solver_success),'Final SM+Weinberg ODE solver did not succeed')
    final=np.asarray(result.K,complex)
    need(np.allclose(final,final.T,rtol=1e-10,atol=1e-29),'SM+Weinberg C5 not symmetric')
    need(np.isfinite(final).all(),'Nonfinite final Weinberg coefficient')
    return result


def evaluate_endpoint(model:str,initial:Mapping[str,Any],endpoint:Mapping[str,Any],control:Mapping[str,Any],
                      mu:float, scenario:Mapping[str,Any], study:Mapping[str,Any],*,evolve:Callable|None=None)->dict[str,Any]:
    need(model in MODELS,'Unknown extended model')
    diag=threshold_diagnostic(endpoint,mu)
    threshold_diagnostic(control,mu)
    low=study['mu_low_gev'];need(mu>low,'Low scale is not below matching scale')
    # Both seeded and zero-seeded UV backgrounds are evolved independently.
    # Even initially-zero extra Yukawas can be radiatively generated.
    phases=[]
    for phase in scenario['phases']:
        full=BR.evaluate_hypothetically_transported_c5(model,BR.to_matchete_inputs(model,endpoint,scenario,extra_phase=phase))
        zero=BR.evaluate_hypothetically_transported_c5(model,BR.to_matchete_inputs(model,control,scenario,extra_phase=phase))
        kfull=BR.from_json(full['C5_total_gev_inverse'])
        kbase=BR.from_json(full['C5_baseline_gev_inverse'])
        kextra=BR.from_json(full['C5_extra_gev_inverse'])
        kzero=BR.from_json(zero['C5_total_gev_inverse'])
        need(np.allclose(kfull,kbase+kextra,rtol=1e-11,atol=1e-30),'Invalid C5 decomposition')
        # Use the same running SM state for baseline and extra decomposition:
        # the Weinberg RGE is linear in C5 for a fixed SM background.
        r_full=transport(kfull,endpoint,mu,low,study['rtol'],study['atol'],evolve=evolve)
        r_base=transport(kbase,endpoint,mu,low,study['rtol'],study['atol'],evolve=evolve)
        r_extra=transport(kextra,endpoint,mu,low,study['rtol'],study['atol'],evolve=evolve)
        r_zero=transport(kzero,control,mu,low,study['rtol'],study['atol'],evolve=evolve)
        resid=float(np.linalg.norm(r_full.K-r_base.K-r_extra.K)/max(float(np.linalg.norm(r_full.K)),1e-300))
        need(resid<5e-6,f'SM+Weinberg C5 linearity check failed: {resid}')
        low_obs=observe(r_full.K,r_full.ye,vev_gev=scenario['vev_gev'],ordering=scenario['ordering'])
        zero_obs=observe(r_zero.K,r_zero.ye,vev_gev=scenario['vev_gev'],ordering=scenario['ordering'])
        phases.append({
            'extra_relative_phase_radians':float(phase),
            'C5_matching_gev_inverse':complex_matrix_json(kfull),
            'C5_low_gev_inverse':complex_matrix_json(r_full.K),
            'C5_low_baseline_gev_inverse':complex_matrix_json(r_base.K),
            'C5_low_direct_extra_gev_inverse':complex_matrix_json(r_extra.K),
            'C5_low_zero_seed_gev_inverse':complex_matrix_json(r_zero.K),
            'C5_low_seeded_minus_zero_seed_gev_inverse':complex_matrix_json(r_full.K-r_zero.K),
            'C5_low_norm_gev_inverse':float(np.linalg.norm(r_full.K)),
            'direct_extra_fraction_at_low':float(np.linalg.norm(r_extra.K)/max(float(np.linalg.norm(r_full.K)),1e-300)),
            'seeded_minus_zero_seed_fraction_at_low':float(np.linalg.norm(r_full.K-r_zero.K)/max(float(np.linalg.norm(r_full.K)),1e-300)),
            'c5_linearity_residual':resid,
            'low_energy_observables':low_obs,
            'zero_seed_low_energy_observables':zero_obs,
            'SM_Weinberg_solver_nfev_full':int(r_full.nfev),
        })
    return {'status':STATUS,'model':model,'mu_match_gev':float(mu),'mu_low_gev':float(low),
            'threshold_diagnostics':diag,'phase_scenarios':phases,
            'one_step_sudden_decoupling_approximation':True,
            'coupling_transport_physical_dictionary_certified':False,
            'hierarchical_EFT1_matching_performed':False,
            'heavy_flavor_full_threshold_matching_certified':False,
            'low_energy_predictions_physical_certified':False}


def run_model(model:str,scenario:Mapping[str,Any],study:Mapping[str,Any], *, root:Path=ROOT,
              buffers:list[float]|None=None,native_source:Path|None=None)->dict[str,Any]:
    need(study['allow_unverified_sudden_decoupling'],'Study must opt into the explicitly UNVERIFIED sudden-decoupling approximation')
    root=Path(root)
    native=Path(native_source) if native_source is not None else root/'Reports/output/workstream_a/extended_c5_native/native'
    proof=C5.verify_model(native,model)
    need(proof.get('exact_reduction_of_all_three_sectors_after_identification') is True,'Native expression proof unavailable')
    payload,benchmark=UV.load_case(model,root/'Reports/output/workstream_a/integrated_general_yukawa_uv/native',root/'Numerical/benchmarks/fitted_reference')
    rg_quartic=BR.check_rgbeta_quartic(model,root,payload)
    initial=UV.prepare_state(model,payload,benchmark,scenario['seed_factor'])
    zero={key:(val.copy() if isinstance(val,np.ndarray) else val) for key,val in initial.items()}
    zero[UV.EXTRAS[model]]=np.zeros((3,3),complex)
    out=[]
    for factor in (buffers if buffers is not None else study['buffers']):
        factor=finite(factor,'endpoint_buffer',positive=True)
        match_start,match_end=UV.endpoint_for(benchmark,initial,factor=factor)
        layout,solver,_=UV.integrate(payload,initial,match_start,match_end,n_save=scenario['n_save'])
        control_layout,control,_=UV.integrate(payload,zero,match_start,match_end,n_save=scenario['n_save'])
        need(layout.names==control_layout.names,'Seeded/control dynamic beta layout changed')
        endpoint=layout.unpack(solver.y[:,-1])
        zero_endpoint=control_layout.unpack(control.y[:,-1])
        compared=evaluate_endpoint(model,initial,endpoint,zero_endpoint,match_end,scenario,study)
        compared['endpoint_buffer']=factor
        compared['mu_uv_start_gev']=float(match_start)
        compared['full_UV_RHS_evaluations']=int(solver.nfev)
        compared['zero_seed_UV_RHS_evaluations']=int(control.nfev)
        compared['full_RGBeta_beta_blocks']=len(layout.names)
        out.append(compared)
    report={
        'model':model,
        'status':STATUS,
        'extra_yukawa':UV.EXTRAS[model],
        'one_step_source_status':proof['status'],
        'one_step_source_sha256':proof['source_sha256'],
        'RGBeta_extra_quartic_source':rg_quartic,
        'hypothesized_matchete_from_rgbeta':{k:BR.encode(v) for k,v in scenario['scales'].items()},
        'native_c5_extra_term_verified':True,
        'matching_scale_scan':out,
        'matching_scale_dependence_not_a_certified_uncertainty_band':True,
        'physical_rgbeta_matchete_coupling_dictionary_certified':False,
        'UV_to_EFT1_sequential_matching_performed':False,
        'intermediate_scalar_eft_running_performed':False,
        'final_SM_Weinberg_RGE_performed_under_sudden_decoupling':True,
        'physical_low_energy_prediction_certified':False,
        'note':'Full coupled UV -> one-step C5 at a common endpoint above all heavy masses -> SM+Weinberg running to low scale. This does NOT implement sequential fermion/scalar thresholds, heavy matching logarithm resummation, full flavour matching certification, pole scalar masses or the absolute RGBeta/Matchete convention.'
    }
    return report


def self_test()->None:
    need(len(MODELS)==6,'Six-model registry mismatch')
    assert set(MODELS)==set(C5.MODELS)==set(UV.EXTRAS)
    scenario={'format':'t3-extended-c5-conditional-low-energy-v1',
              'physical_dictionary_certified':False,'hierarchical_threshold_matching_certified':False,
              'allow_unverified_sudden_decoupling':True,'endpoint_buffers':[1.15,1.5],
              'mu_low_gev':91.1876}
    cfg=load_study(scenario)
    need(cfg['buffers']==[1.15,1.5],'Study parsing failure')
    for change in ({'physical_dictionary_certified':True},
                   {'hierarchical_threshold_matching_certified':True},
                   {'endpoint_buffers':[1.0]}, {'mu_low_gev':0}):
        try:load_study({**scenario,**change})
        except ValueError:pass
        else:raise AssertionError('Invalid high-impact study configuration accepted '+str(change))
    state={'MF':np.diag([1000.,1100.,1250.]).astype(complex),
           'mS1Sq':600.**2,'mS2Sq':700.**2}
    assert threshold_diagnostic(state,1800.)['fermion_first_hierarchy_in_mass_proxies']
    try:threshold_diagnostic(state,1250.)
    except ValueError:pass
    else:raise AssertionError('UV threshold-crossing negative control accepted')
    print('PASS: six-model study configuration, physical gates, mass ordering, and threshold crossing controls')
