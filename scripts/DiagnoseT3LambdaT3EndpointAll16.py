"""Fitted EFT1 lambdaT3(muS) and coupling-normalised LLSS -> Weinberg response.

Project path: scripts/DiagnoseT3LambdaT3EndpointAll16.py

Reconstructs the U=1,E=1 renormalisable trajectories using the SAME
implementations/settings as RunT3GeneralSequentialFitAblation. No refitting,
LLSS matching, Weinberg running, or SMEFT evolution is performed. Reuses the
saved source-resolved leading logs and previously extracted C_LLSS norms.

B_-1/C_-1: preserves the existing muS scalar basis approximation; numerical
normalisation of their direct quartic channels is not independently certified.
The ratio kappa/|lambdaT3| is not an anomalous-dimension eigenvalue.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import traceback
import numpy as np

from scripts import RunT3GeneralSequentialFitAblation as seq
from scripts import RunT3GeneralAll16FromReference as general

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / 'Reports/output/thesis_qualitative/eft1_llss_contributions'
CAMPAIGN = ROOT / 'output/studies/general_scalar_extended/corrected_llss_eft1_freeze'


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def as_complex(v):
    if isinstance(v,dict) and {'real','imag'} <= set(v):
        return complex(float(v['real']),float(v['imag']))
    return complex(v)


def encode(v):
    z = as_complex(v)
    if not np.isfinite([z.real,z.imag]).all():
        raise ValueError('Nonfinite lambdaT3')
    return {'real':z.real,'imag':z.imag}


def reconstructed_nonmixed(model, saved, canonical, wolframscript):
    from Numerical.orchestration.GeneralRunningBenchmark import initial_state, _state_values
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment

    uv_beta, eft_beta, _, _, _, _ = general.nonmixed_inputs(
        model, canonical, generate_general_c5=False, wolframscript=wolframscript)
    overrides = read(saved/'general_uv_overrides.json')
    _, initial = initial_state(uv_beta, mu_uv=seq.SCALES[0],
                              mu_f=seq.SCALES[1],mu_s=seq.SCALES[2],
                              overrides=overrides)
    uv_solution = run_uv_segment(initial,uv_beta,seq.SCALES[1],
                                 rtol=1e-7,atol=1e-10)
    uv_end = uv_solution.final_state
    eft_start = build_intermediate_scalar_boundary(
        uv_end,matching_scale_gev=seq.SCALES[1])
    from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state
    _,layout = pack_intermediate_scalar_state(eft_start)
    if set(layout.names)!=set(eft_beta.get('report_betas',{})):
        raise ValueError('EFT1 coupling inventory mismatch')
    eft_solution = run_intermediate_scalar_segment(
        eft_start,eft_beta,seq.SCALES[2],rtol=1e-7,atol=1e-10)
    start_value=as_complex(_state_values(eft_start)['lambdaT3'])
    end_value=as_complex(_state_values(eft_solution.final_state)['lambdaT3'])
    return start_value,end_value,{
        'UV_nfev':int(uv_solution.nfev),
        'EFT1_nfev':int(eft_solution.nfev),
        'fit_overrides_sha256':sha(saved/'general_uv_overrides.json'),
        'state_origin':'GeneralRunningBenchmark._state_values',
        'mixed_muS_basis_assumption':False,
    }


def reconstructed_mixed(model,saved):
    from scripts import RunT3BCGeneralFromReference as runner
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc
    from scripts import RunT3GeneralScalarBCPilot as pilot
    from scripts import RunT3ZeroExtraBoundaryBenchmark as core
    from Numerical.experimental.GeneralT3ScalarRGState import (
        unpack_general_scalar_state,project_general_after_fermion,
        pack_general_scalar_state)
    from Numerical.orchestration.GeneralRunningBenchmark import _state_values

    fitted=read(saved/'result.json')
    source=runner.build_reference_source(model)
    source['cfg']=read(saved/'fitted_config.json')
    registry,new,reused=source['registry'],source['new'],source['reused']
    q,m=seq.mixed_fit_fractions(fitted)
    beta_bundle=bc.beta_bundle(model,source,registry)
    initial,uv_layout,_,_=pilot.prepare_initial(
        source,registry,'general_quartics_and_mixing' if m>0 else 'general_quartics',q,m)
    uv_sol=core.integrate_segment(
        initial,lambda t,x:core.rhs_at(t,x,layout=uv_layout,
             registry=registry,betas=beta_bundle['UV'],stage='UV'),
        seq.SCALES[0],seq.SCALES[1],n_points=3)
    uv_vec=uv_sol.y[:,-1].copy()
    uv=unpack_general_scalar_state(uv_vec[:uv_layout.size],uv_layout,
                                   registry,uv=True,mu_gev=seq.SCALES[1])
    eft=project_general_after_fermion(uv,registry)
    eft_vec,eft_layout=pack_general_scalar_state(eft,registry)
    eft_vec=np.concatenate((eft_vec,uv_vec[uv_layout.size:]))
    eft_beta=core.check_payload(
        beta_bundle['payloads']['EFT1'],model=model,stage='EFT1',
        rep=eft_layout.representation,
        base_names=tuple(eft_layout.names[:-len(new)]),
        new_names=new,reused_names=reused)
    eft_sol=core.integrate_segment(
        eft_vec,lambda t,x:core.rhs_at(t,x,layout=eft_layout,
             registry=registry,betas=eft_beta,stage='EFT1'),
        seq.SCALES[1],seq.SCALES[2],n_points=3)
    # This follows the exact endpoint construction of the corrected campaign.
    first=unpack_general_scalar_state(
        eft_vec[:eft_layout.size],eft_layout,registry,
        uv=False,mu_gev=seq.SCALES[1]).baseline
    last=unpack_general_scalar_state(
        eft_sol.y[:,-1][:eft_layout.size],eft_layout,registry,
        uv=False,mu_gev=seq.SCALES[2]).baseline
    start=as_complex(_state_values(first)['lambdaT3'])
    end=as_complex(_state_values(last)['lambdaT3'])
    return start,end,{
        'UV_nfev':int(uv_sol.nfev),
        'EFT1_nfev':int(eft_sol.nfev),
        'fitted_config_sha256':sha(saved/'fitted_config.json'),
        'fitted_result_sha256':sha(saved/'result.json'),
        'state_origin':'GeneralT3ScalarRGState baseline -> _state_values',
        'mixed_muS_basis_assumption':True,
    }


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-run',default='thesis16_merged_v1')
    parser.add_argument('--run-id',default='corrected_llss16_v4')
    parser.add_argument('--models',default=None)
    parser.add_argument('--wolframscript',default='wolframscript')
    parser.add_argument('--output',type=Path,default=None)
    args=parser.parse_args(argv)
    models=tuple(s.strip() for s in args.models.split(',') if s.strip()) if args.models else seq.MODELS
    if not models or len(set(models))!=len(models) or set(models)-set(seq.MODELS):
        parser.error('Invalid model IDs')
    source_root=general.OUTPUT/args.source_run
    settings=read(source_root/'settings.json')
    canonical=settings.get('canonical_run_id','full16_general_v1')
    case_root=CAMPAIGN/args.run_id
    threshold_root=REPORT/args.run_id/'fermion_threshold_all16'
    out=args.output or REPORT/args.run_id/'lambdaT3_endpoint_normalized_all16'
    out=out.resolve()
    if out==case_root.resolve() or out==threshold_root.resolve() or case_root.resolve() in out.parents or threshold_root.resolve() in out.parents:
        raise ValueError('Refusing to write into a source directory')
    rows=[]; details={};errors=[]
    for model in models:
        try:
            fit_dir=source_root/seq.tag(model)
            case_path=case_root/seq.tag(model)/'cases.json'
            threshold_path=threshold_root/seq.tag(model)/'summary.json'
            case=read(case_path)
            threshold=read(threshold_path)
            if case.get('model')!=model or case.get('status')!='Computed' or threshold.get('model')!=model:
                raise ValueError('Saved model/status mismatch')
            if threshold.get('reference_cases_sha256') and threshold['reference_cases_sha256']!=sha(case_path):
                raise ValueError('Previously extracted LLSS norm uses different saved cases')
            branch=case['cases']['1111']
            if not branch['uv_on'] or not branch['eft1_on'] or not branch['wilson_log_on']:
                raise ValueError('1111 switch mismatch')
            if model in ('B_-1','C_-1'):
                initial,end,provenance=reconstructed_mixed(model,fit_dir)
            else:
                initial,end,provenance=reconstructed_nonmixed(model,fit_dir,canonical,args.wolframscript)
            if abs(end)==0 or not math.isfinite(abs(end)):
                raise ValueError('Cannot divide by zero/nonfinite lambdaT3(muS)')
            raw=float(threshold['C_LLSS_formal_frobenius_gev_inverse'])
            if raw<=0 or not math.isfinite(raw):
                raise ValueError('Invalid saved LLSS norm')
            quartics=[x for x in branch['per_general_interaction'] if x['category']=='quartic_decomposition']
            t3=[x for x in quartics if x['name']=='native:lambdaT3']
            if len(t3)!=1:
                raise ValueError('Expected one lambdaT3 channel')
            def mat(raw):
                z=np.asarray(raw['real'],float)+1j*np.asarray(raw['imag'],float)
                if z.shape!=(3,3) or not np.isfinite(z).all():
                    raise ValueError('Invalid Weinberg matrix')
                return z
            contribution=mat(t3[0]['c5_leading_log'])
            total=mat(branch['c5_induced_wilson_log'])
            if np.linalg.norm(sum((mat(x['c5_leading_log']) for x in quartics),np.zeros((3,3),complex))-total)>1e-7*max(np.linalg.norm(total),1e-40):
                raise ValueError('Quartic decomposition fails additivity')
            kappa=float(np.linalg.norm(contribution)/raw)
            prefactor=math.log(seq.SCALES[1]/seq.SCALES[2])/(16*math.pi**2)
            row={
                'model':model,'lambdaT3_muF_real':initial.real,
                'lambdaT3_muF_imag':initial.imag,'lambdaT3_muS_real':end.real,
                'lambdaT3_muS_imag':end.imag,'abs_lambdaT3_muS':abs(end),
                'relative_lambdaT3_running':abs(end-initial)/max(abs(initial),1e-40),
                'LLSS_norm_GeV_inv':raw,
                'C5_lambdaT3_log_norm_GeV_inv':float(np.linalg.norm(contribution)),
                'C5_total_log_norm_GeV_inv':float(np.linalg.norm(total)),
                'kappa_T3':kappa,
                'kappa_T3_over_abs_lambdaT3_muS':kappa/abs(end),
                'kappa_T3_over_looplog_and_abs_lambdaT3':kappa/(prefactor*abs(end)),
                'signed_T3_projection_percent':100*float(np.vdot(total,contribution).real/np.linalg.norm(total)**2),
                'UV_nfev':provenance['UV_nfev'],'EFT1_nfev':provenance['EFT1_nfev'],
                'mixed_muS_basis_assumption':provenance['mixed_muS_basis_assumption'],
            }
            rows.append(row)
            details[model]={
                'lambdaT3_muF':encode(initial),'lambdaT3_muS':encode(end),
                'endpoint_provenance':provenance,
                'cases_sha256':sha(case_path),
                'threshold_summary_sha256':sha(threshold_path),
            }
            print(f'{model}: lambdaT3(muS)={end.real:.9e}{end.imag:+.9e}j; '
                  f'kappa_T3={kappa:.6e}; kappa/|lambdaT3|={row["kappa_T3_over_abs_lambdaT3_muS"]:.6e}; '
                  f'RG change={row["relative_lambdaT3_running"]:.3e}',flush=True)
        except Exception as exc:
            errors.append({'model':model,'error':f'{type(exc).__name__}: {exc}',
                           'traceback':traceback.format_exc(limit=6)})
            print(f'{model}: FAILED {type(exc).__name__}: {exc}',flush=True)
    out.mkdir(parents=True,exist_ok=True)
    columns=['model','lambdaT3_muF_real','lambdaT3_muF_imag','lambdaT3_muS_real',
             'lambdaT3_muS_imag','abs_lambdaT3_muS','relative_lambdaT3_running',
             'LLSS_norm_GeV_inv','C5_lambdaT3_log_norm_GeV_inv',
             'C5_total_log_norm_GeV_inv','kappa_T3',
             'kappa_T3_over_abs_lambdaT3_muS',
             'kappa_T3_over_looplog_and_abs_lambdaT3',
             'signed_T3_projection_percent','UV_nfev','EFT1_nfev',
             'mixed_muS_basis_assumption']
    with (out/'all16_lambdaT3_endpoint_normalized.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)
    (out/'endpoint_details.json').write_text(json.dumps(details,indent=2)+'\n',encoding='utf-8')
    (out/'manifest.json').write_text(json.dumps({
        'source_run':args.source_run,'reference_run':args.run_id,
        'requested_models':list(models),'completed_models':[x['model'] for x in rows],
        'errors':errors,'muF_GeV':seq.SCALES[1],'muS_GeV':seq.SCALES[2],
        'running':'same renormalisable UV/EFT1 numerical trajectory as corrected campaign; no refit',
        'C5_and_CLLSS':'read saved previous campaigns and diagnostic; not recomputed',
        'lambdaT3_value':'reconstructed EFT1 renormalisable endpoint at muS, not saved in cases',
        'absolute_spinor_normalization_certified':False,
        'mixed_source_muS_basis_assumption':'B_-1/C_-1 source includes scalar mixing U(muS)',
        'BC_channel_convention_certified':False,
        'coupling_normalized_kappa_not_anomalous_dimension':True,
        'full_hard_fermion_C5_matching_calculated':False,
    },indent=2)+'\n',encoding='utf-8')
    print('Output:',out)
    return 0 if not errors else 2


if __name__=='__main__':
    raise SystemExit(main())
