"""T3 six self-bilinears: candidate anomalous-scalar HH-quartic C5 continuation.

Uses saved continuous_v1, neutrino_v1 and self_bilinear_v1 results. No UV/EFT1
re-evolution, Wolfram run, new fit or source modification. The sole one-HH
quartic is tensor-matched to the existing T3MixCG and its loop is formed by
anomalous spectral subtraction of the existing model-specific C5 kernel.

This is conditional on lambda_Matchete = sign * lambda_RGBeta and the native
Matchete-to-RGBeta coupling/bilinear-basis identification. No sequential
F-first C12/C5 Wilson matching is performed.

python -m scripts.RunT3SelfBilinearAnomalousC5 --all --run-id continuous_v1
python -m scripts.RunT3SelfBilinearAnomalousC5 --all --run-id continuous_v1 --resume
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
from pathlib import Path

import numpy as np

from Numerical.experimental.T3SelfBilinearAnomalousKernel import (
    MODELS, EXTRA, model_alpha, tensor_ratio, anomalous_phase,
    anomalous_from_unit_kernels, unit_topology_kernels,
)
from scripts import RunT3FullScalarContinuous as continuous
from scripts import RunT3FullScalarNeutrino as neutrino
from scripts import RunT3SelfBilinearNeutrino as previous

ROOT=continuous.ROOT
FORMAT='t3-self-bilinear-anomalous-HH-C5-candidate-v1'
STATUS='CandidateSelfBilinearNormalPlusAnomalousHHSourceTensorC5'
STAGE22=ROOT/'output/general_t3_stage22/matchete_tensors.json'
DEFAULT_OUTPUT_ID='self_anomalous_v1'


def require(ok:bool,message:str)->None:
    if not ok:
        raise ValueError('T3 self anomalous C5: '+message)


def fingerprint(path:Path)->str:
    p=Path(path)
    require(p.is_file(),'missing required input '+str(p))
    return sha256(p.read_bytes()).hexdigest()


def group_path(model:str)->Path:
    require(model in MODELS,'unexpected model')
    return ROOT/f'output/group_factors/T3_{model[0]}_uv_cg_registry.json'


def native_matching(model:str,report:dict)->dict:
    rep=report['eft1_scalar_endpoint_state']['representation']
    ds1,ds2=rep['d_s1'],rep['d_s2']
    expected={'A':(1,3),'D':(3,1),'E':(3,3)}[model[0]]
    require((ds1,ds2)==expected,'model dimensions changed')
    group=neutrino.load(group_path(model))
    meta=group.get('metadata',{})
    require(meta.get('dS1')==ds1 and meta.get('dS2')==ds2 and
            meta.get('dF')==rep['d_f'] and meta.get('alpha')==0,
            'class-level native topology CG registry incompatible')
    return tensor_ratio(neutrino.load(STAGE22),group,model,ds1,ds2)


def endpoint_extra(model:str,source:dict)->complex:
    """Read and cross-check the saved full coupled EFT1 scalar quartic."""
    alpha=model_alpha(model)
    symbol=EXTRA[alpha]
    rows=source.get('eft1_trajectory')
    require(isinstance(rows,list) and len(rows)>1,'missing coupled EFT1 trajectory')
    raw=rows[-1].get('extra_quartics',{}).get(symbol)
    require(raw is not None,'missing evolved extra HH coupling '+symbol)
    val=neutrino.complex_pair(raw,symbol)
    return val


def relative_delta(new:np.ndarray,base:np.ndarray)->float|None:
    d=float(np.linalg.norm(base))
    return float(np.linalg.norm(new-base)/d) if d>0 else None


def calculate(model:str,*,source_root:Path,output_root:Path,
              previous_output_id:str='self_bilinear_v1',
              control_output_id:str='neutrino_v1',
              extra_sign:int=-1,resume:bool=False)->dict:
    require(model in MODELS,'only six A/D/E self-bilinears supported')
    require(extra_sign in (-1,0,1),'--extra-sign must be -1, 0 or +1')
    tag=continuous.model_tag(model)
    record_path=source_root/tag/'continuous_trajectory.json'
    source=neutrino.verify_continuous(model,record_path)
    old_path=source_root/previous_output_id/(tag+'_self_bilinear_c5.json')
    old=neutrino.load(old_path)
    require(old.get('format')==previous.FORMAT and old.get('status')==previous.STATUS and
            old.get('source_continuous_sha256')==fingerprint(record_path) and
            old.get('normal_channel_scalar_bilinear_included') is True and
            old.get('anomalous_prop_matching_included') is False and
            old.get('full_general_scalar_c5_certified') is False,
            'previous normal-channel certificate missing/incompatible')
    basis=previous.candidate_basis(model,source)
    uv=neutrino.decode_uv(source['uv_fermion_threshold_state'])
    eft=neutrino.decode_eft(source['eft1_scalar_endpoint_state'])
    from scripts.RunT3GeneralScalar48Campaign import matching
    kernel,kernel_path,kind=matching(model)
    require(fingerprint(kernel_path)==old['source_matching_kernel_sha256'],
            'matching source no longer agrees with earlier normal result')
    signatures={
        'continuous':fingerprint(record_path),
        'normal_c5':fingerprint(old_path),
        'matching':fingerprint(kernel_path),
        'stage22_extra_tensor':fingerprint(STAGE22),
        'native_class_topology_tensor':fingerprint(group_path(model)),
        'anomalous_runner':fingerprint(Path(__file__).resolve()),
        'anomalous_kernel':fingerprint(ROOT/'Numerical/experimental/T3SelfBilinearAnomalousKernel.py'),
    }
    target=output_root/(tag+'_anomalous_hh_c5.json')
    if target.exists():
        data=neutrino.load(target)
        require(resume and data.get('format')==FORMAT and data.get('status')==STATUS and
                data.get('input_sha256')==signatures and data.get('matchete_extra_sign')==extra_sign,
                'existing result cannot be reused; use --resume with unchanged sources or a new --output-id')
        return data
    metric=native_matching(model,source)
    normal=neutrino.complex_matrix(old['normal_propagator_resummed_c5_at_1tev'],'normal topology C5')
    saved_low=neutrino.complex_matrix(old['c5_at_1tev_gev_inverse'],'normal reported C5')
    require(np.allclose(normal,saved_low,rtol=1e-12,atol=1e-20),
            'previous normal report has inconsistent C5 entries')
    plus,minus,parity=unit_topology_kernels(kernel,uv,eft,basis,normal)
    rgbeta=endpoint_extra(model,source)
    matchete=complex(extra_sign)*rgbeta  # explicit candidate map, not certified
    anom=anomalous_from_unit_kernels(plus,minus,mu_sq=basis.mu_squared,
            alpha=model_alpha(model),matchete_extra=matchete,
            tensor_factor=metric['ratio'])
    combined=neutrino.validated_c5(normal+anom)
    normal_low=neutrino.low_energy(model,normal,eft)
    require(np.allclose(
          neutrino.complex_matrix(normal_low['c5_at_100gev_gev_inverse'],'normal low C5'),
          neutrino.complex_matrix(old['c5_at_100gev_gev_inverse'],'saved normal low C5'),
          rtol=4e-9,atol=1e-20),
          'previous normal low-energy running does not reproduce')
    low=neutrino.low_energy(model,combined,eft)
    norm_mnu=neutrino.complex_matrix(normal_low['mnu_at_100gev_ev'],'normal low-energy mnu')
    new_mnu=neutrino.complex_matrix(low['mnu_at_100gev_ev'],'new low-energy mnu')
    payload={
        'format':FORMAT,'status':STATUS,'model':model,
        'source_continuous_sha256':signatures['continuous'],
        'source_normal_c5_sha256':signatures['normal_c5'],
        'input_sha256':signatures,
        'matching_kernel_kind':kind,
        'quadratic_kind':source['quadratic_kind'],
        'self_species':basis.species,'self_component_dimension':basis.dim,
        'mu2_at_scalar_endpoint_gev2':continuous.cjson(basis.mu_squared),
        'split_scalar_mass_squared_gev2':list(basis.eigen_masses_squared),
        'anomalous_internal_scalar_phase':continuous.cjson(anomalous_phase(basis.mu_squared,model_alpha(model))),
        'native_matchete_tensor_ratio':continuous.cjson(metric['ratio']),
        'native_tensor_proportionality_residual':metric['relative_tensor_residual'],
        'native_mass_metric_delta_residual':metric['native_self_metric_residual'],
        'source_matchete_hh_symbol':metric['matchete_extra_symbol'],
        'source_rgbeta_hh_symbol':metric['rgbeta_extra_symbol'],
        'class_level_SU2_topology_CG_from_alpha_zero':True,
        'rgbeta_hh_coupling_at_scalar_endpoint':continuous.cjson(rgbeta),
        'matchete_extra_sign':extra_sign,
        'assumed_matchete_hh_coupling':continuous.cjson(matchete),
        'matchete_hermitian_conjugate_HH_coupling':continuous.cjson(matchete.conjugate()),
        'original_lambdaT3_at_scalar_endpoint':continuous.cjson(complex(eft.lambdaT3)),
        'unit_topology_kernel_linearity_relative_error':parity,
        'original_topology_normal_c5_at_1tev':neutrino.matrix_payload(normal),
        'extra_quartic_anomalous_c5_at_1tev':neutrino.matrix_payload(anom),
        'normal_plus_anomalous_c5_at_1tev':neutrino.matrix_payload(combined),
        'anomalous_vs_normal_c5_amplitude_norm_ratio':float(np.linalg.norm(anom)/np.linalg.norm(normal)) if np.linalg.norm(normal)>0 else None,
        'relative_combined_c5_change':relative_delta(combined,normal),
        'normal_baseline_neutrino_masses_ev_descending':normal_low['takagi_masses_ev_descending'],
        'relative_low_energy_mnu_change':relative_delta(new_mnu,norm_mnu),
        'anomalous_self_propagator_and_one_HH_extra_quartic_included':True,
        'normal_propagator_quadratic_mass_resummation_included':True,
        'all_additional_quartic_direct_one_loop_contributions_certified':False,
        'native_Matchete_RGBeta_bilinear_mapping_certified':False,
        'native_Matchete_RGBeta_quartic_mapping_certified':False,
        'full_sequential_wilson_matching_performed':False,
        'full_general_scalar_matching_physically_certified':False,
        'matching_scope':('One-loop original topology normal propagator plus single anomalous '
            'self-scalar propagator with the one independent Hdag Hdag extra quartic. '
            'Source native Matchete tensor proportionality checked; RGBeta sign and '
            'basis identification conditional. One-step spectral reweighting of existing '
            'hierarchical/legacy source kernel is NOT a new explicit sequential EFT match.'),
        **low,
    }
    neutrino.atomic_json(target,payload)
    return payload


def merge_sixteen(*,source_root:Path,output_root:Path,
                  control_output_id:str='neutrino_v1')->dict:
    """Keep ten existing results and replace six restricted controls in a NEW table."""
    require(set(MODELS).issubset(set(continuous.MODELS)) and len(continuous.MODELS)==16,
            'unexpected 16-model list')
    rows=[]
    for model in continuous.MODELS:
        tag=continuous.model_tag(model)
        expected_hash=fingerprint(source_root/tag/'continuous_trajectory.json')
        if model in MODELS:
            path=output_root/(tag+'_anomalous_hh_c5.json')
            report=neutrino.load(path)
            require(report.get('format')==FORMAT and report.get('status')==STATUS and
                    report.get('source_continuous_sha256')==expected_hash and
                    report.get('anomalous_self_propagator_and_one_HH_extra_quartic_included') is True and
                    report.get('full_general_scalar_matching_physically_certified') is False,
                    model+': anomalous candidate is incompatible')
            category='self-bilinear-normal-plus-anomalous-HH-candidate'
        else:
            path=source_root/control_output_id/(tag+'_neutrino.json')
            report=neutrino.load(path)
            require(report.get('format')==neutrino.FORMAT and
                    report.get('status') in ('CandidateNoAdditionalQuadraticMatchedC5',
                        'CandidateMassEigenBCEvolvedYukawasSMWeinberg') and
                    report.get('source_continuous_sha256')==expected_hash and
                    report.get('full_general_scalar_physical_c5_certified') is False and
                    report.get('full_sequential_threshold_matching_performed') is False and
                    report.get('additional_quadratic_c5_omitted') is False,
                    model+': other candidate is incompatible')
            category='existing-candidate-one-loop-matching'
        masses=report.get('takagi_masses_ev_descending')
        require(isinstance(masses,list) and len(masses)==3 and
                all(np.isfinite(z) and z>=0 for z in masses),
                model+': missing finite neutrino eigenvalues')
        rows.append({'model':model,'status':report['status'],
            'matching_category':category,'source_report':path.relative_to(ROOT).as_posix(),
            'source_sha256':fingerprint(path),
            'largest_neutrino_mass_ev':masses[0],
            'middle_neutrino_mass_ev':masses[1],
            'smallest_neutrino_mass_ev':masses[2],
            **report['oscillation_observables_normal_ordering']})
    payload={'format':'t3-16-model-candidate-anomalous-HH-neutrino-summary-v1',
             'models_processed':len(rows),
             'models_anomalous_self_bilinear_included':len(MODELS),
             'other_models_existing_candidate_results':len(rows)-len(MODELS),
             'all_16_numerical_results_available':True,
             'physical_full_general_scalar_matching_certified':False,
             'sequential_wilson_matching_performed':False,
             'rows':rows}
    output_root.mkdir(parents=True,exist_ok=True)
    neutrino.atomic_json(output_root/'combined_16_summary.json',payload)
    columns=['model','status','matching_category','largest_neutrino_mass_ev',
             'middle_neutrino_mass_ev','smallest_neutrino_mass_ev',
             'delta_m21_sq_ev2','delta_m3l_sq_ev2',
             'sin2_theta12','sin2_theta13','sin2_theta23']
    tmp=output_root/'combined_16_summary.csv.tmp'
    with tmp.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)
    tmp.replace(output_root/'combined_16_summary.csv')
    return payload


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    selector=ap.add_mutually_exclusive_group(required=True)
    selector.add_argument('--all',action='store_true')
    selector.add_argument('--model',choices=MODELS)
    ap.add_argument('--run-id',required=True,type=continuous.run_id)
    ap.add_argument('--source-root',type=Path,default=continuous.OUTPUT)
    ap.add_argument('--previous-output-id',default='self_bilinear_v1')
    ap.add_argument('--control-output-id',default='neutrino_v1')
    ap.add_argument('--output-id',default=DEFAULT_OUTPUT_ID)
    ap.add_argument('--extra-sign',type=int,choices=(-1,0,1),default=-1,
                    help='Explicit conditional lambda_M=sign*lambda_RGBeta; default -1')
    ap.add_argument('--resume',action='store_true')
    ap.add_argument('--combine-only',action='store_true')
    args=ap.parse_args()
    for x in (args.previous_output_id,args.control_output_id,args.output_id):
        continuous.run_id(x)
    require(len({args.previous_output_id,args.control_output_id,args.output_id})==3,
            'source and output IDs must be distinct')
    folder=args.source_root/args.run_id
    require(folder.is_dir(),'missing completed continuous trajectory '+str(folder))
    output=folder/args.output_id
    if args.combine_only:
        require(args.all,'--combine-only requires --all')
        summary=merge_sixteen(source_root=folder,output_root=output,
                               control_output_id=args.control_output_id)
        print('Combined:',output/'combined_16_summary.csv',
              'models=',summary['models_processed'])
        return 0
    requested=MODELS if args.all else (args.model,)
    rows=[]
    for model in requested:
        try:
            data=calculate(model,source_root=folder,output_root=output,
                previous_output_id=args.previous_output_id,
                control_output_id=args.control_output_id,
                extra_sign=args.extra_sign,resume=args.resume)
            rows.append({'model':model,'status':data['status'],
                'masses_ev':data['takagi_masses_ev_descending'],
                'relative_anomalous_amplitude':data['anomalous_vs_normal_c5_amplitude_norm_ratio']})
            amplitude_ratio=data['anomalous_vs_normal_c5_amplitude_norm_ratio']
            ratio_text=(f'{amplitude_ratio:.6g}' if amplitude_ratio is not None
                        else 'undefined-zero-baseline')
            print(f'{model}: normal+anomalous HH C5; masses={data["takagi_masses_ev_descending"]}; '
                  f'|C_anom|/|C_normal|={ratio_text}',flush=True)
        except (ValueError,KeyError,TypeError,RuntimeError,OSError,OverflowError,
                FloatingPointError,np.linalg.LinAlgError) as exc:
            rows.append({'model':model,'status':'Blocked','error':str(exc)[:700]})
            print(f'{model}: BLOCKED {exc}',flush=True)
        neutrino.atomic_json(output/'summary.json',{
            'format':FORMAT,'requested_models':list(requested),
            'completed':sum(z['status']==STATUS for z in rows),
            'blocked':sum(z['status']=='Blocked' for z in rows),
            'matchete_extra_sign':args.extra_sign,
            'physical_matching_certified':False,'rows':rows})
    if args.all and all(z['status']!= 'Blocked' for z in rows):
        summary=merge_sixteen(source_root=folder,output_root=output,
                               control_output_id=args.control_output_id)
        print('Combined candidate table:',output/'combined_16_summary.csv',
              'models=',summary['models_processed'])
    print('Report:',output/'summary.json')
    return 2 if any(z['status']=='Blocked' for z in rows) else 0

if __name__=='__main__':
    raise SystemExit(main())
