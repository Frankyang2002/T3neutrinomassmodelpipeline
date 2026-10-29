"""Corrected six-model A/D/E self-bilinear C5 using mass-dependent one-scale source.

The previous hierarchical kernel has NO MS1/MS2 dependence and therefore
forces an identically zero anomalous spectral difference. This runner uses
saved mass-dependent Matchete one-scale c5_coefficient.txt instead; it does
NOT certify physical sequential matching or the RGBeta/Matchete sign map.

    python -m scripts.RunT3SelfBilinearMassDependentC5 --all --run-id continuous_v1
    python -m scripts.RunT3SelfBilinearMassDependentC5 --all --run-id continuous_v1 --resume
"""
from __future__ import annotations

import argparse
import csv
import math
from hashlib import sha256
from pathlib import Path

import numpy as np

from Numerical.experimental.T3SelfBilinearAnomalousKernel import (
    MODELS, tensor_ratio, model_alpha, anomalous_phase,
)
from Numerical.experimental.T3SelfBilinearMassDependentKernel import (
    source_path, verify_legacy_source, split_unit_kernels, spectral_c5,
)
from scripts import RunT3FullScalarContinuous as continuous
from scripts import RunT3FullScalarNeutrino as neutrino
from scripts import RunT3SelfBilinearNeutrino as prior

ROOT = continuous.ROOT
FORMAT = 't3-self-bilinear-massdependent-one-scale-HH-C5-candidate-v1'
STATUS = 'CandidateMassDependentOneScaleNormalPlusAnomalousHH'
STAGE22 = ROOT/'output/general_t3_stage22/matchete_tensors.json'
DEFAULT_OUTPUT_ID = 'self_anomalous_massdep_v1'


def require(value: bool, message: str) -> None:
    if not value:
        raise ValueError('Self-bilinear mass-dependent one-scale C5: '+message)


def sha(path: Path) -> str:
    require(Path(path).is_file(), 'missing '+str(path))
    return sha256(Path(path).read_bytes()).hexdigest()


def cg_path(model: str) -> Path:
    require(model in MODELS, 'unsupported model')
    return ROOT/f'output/group_factors/T3_{model[0]}_uv_cg_registry.json'


def relative(a: np.ndarray, b: np.ndarray) -> float | None:
    denom=float(np.linalg.norm(b))
    return float(np.linalg.norm(a-b)/denom) if denom>0 else None


def calculate(model: str, *, source_root: Path, output_root: Path,
              sign: int=-1, resume: bool=False) -> dict:
    require(model in MODELS, 'unsupported model')
    require(sign in (-1,0,1), 'extra sign must be -1, 0 or +1')
    tag=continuous.model_tag(model)
    trajectory=source_root/tag/'continuous_trajectory.json'
    source=neutrino.verify_continuous(model,trajectory)
    basis=prior.candidate_basis(model,source)
    uv=neutrino.decode_uv(source['uv_fermion_threshold_state'])
    eft=neutrino.decode_eft(source['eft1_scalar_endpoint_state'])
    matching_source=source_path(ROOT,model)
    witness=verify_legacy_source(matching_source)
    from scripts.RunT3GeneralScalar48Campaign import LegacyKernel
    kernel=LegacyKernel(matching_source)

    stage=neutrino.load(STAGE22)
    group=neutrino.load(cg_path(model))
    rep=source['eft1_scalar_endpoint_state']['representation']
    expect={'A':(1,3),'D':(3,1),'E':(3,3)}[model[0]]
    require((rep['d_s1'],rep['d_s2'])==expect, 'unexpected representation')
    metric=tensor_ratio(stage,group,model,*expect)
    alpha=model_alpha(model)
    extra_name=metric['rgbeta_extra_symbol']
    rows=source.get('eft1_trajectory')
    require(isinstance(rows,list) and len(rows)>=2, 'missing scalar EFT1 trajectory')
    extra_data=rows[-1].get('extra_quartics',{})
    require(extra_name in extra_data, 'missing evolved HH extra quartic '+extra_name)
    extra_rgbeta=neutrino.complex_pair(extra_data[extra_name],extra_name)
    extra_matchete=sign*extra_rgbeta

    inputs={
        'continuous':sha(trajectory),
        'one_scale_matchete':sha(matching_source),
        'native_extra_tensors':sha(STAGE22),
        'native_topology_cg':sha(cg_path(model)),
        'runner':sha(Path(__file__).resolve()),
        'massdependent_kernel':sha(ROOT/'Numerical/experimental/T3SelfBilinearMassDependentKernel.py'),
        'tensor_kernel':sha(ROOT/'Numerical/experimental/T3SelfBilinearAnomalousKernel.py'),
    }
    target=output_root/f'{tag}_massdependent_c5.json'
    if target.is_file():
        old=neutrino.load(target)
        require(resume and old.get('format')==FORMAT and old.get('status')==STATUS
                and old.get('input_sha256')==inputs and old.get('extra_sign')==sign,
                'output exists or inputs changed; choose a new --output-id or --resume')
        return old

    plus,minus,sensitivity=split_unit_kernels(kernel,uv,eft,basis)
    channels=spectral_c5(plus,minus,lam_t3=eft.lambdaT3,
            mu_sq=basis.mu_squared,alpha=alpha,
            extra_matchete=extra_matchete,tensor_factor=metric['ratio'])
    # The old zero resulted from a mass-insensitive kernel. Do not silently
    # accept another zero when BOTH the source quartic and mass split are nonzero.
    if sign!=0 and basis.mu_squared!=0 and extra_rgbeta!=0:
        require(np.linalg.norm(channels['anomalous'])>0,
                'unexpected zero with nonzero extra quartic and mass splitting')
    normal_low=neutrino.low_energy(model,channels['normal'],eft)
    low=neutrino.low_energy(model,channels['combined'],eft)
    m0=neutrino.complex_matrix(normal_low['mnu_at_100gev_ev'],'normal low mnu')
    m1=neutrino.complex_matrix(low['mnu_at_100gev_ev'],'new low mnu')
    out={
        'format':FORMAT,'status':STATUS,'model':model,'input_sha256':inputs,
        'matching_kernel_type':'mass-dependent legacy one-scale Matchete C5',
        'source_mass_dependence_check':witness,
        'previous_hierarchical_mass_independent_kernel_rejected':True,
        'quadratic_kind':source['quadratic_kind'],
        'self_species':basis.species,
        'mu2_at_scalar_endpoint_gev2':continuous.cjson(basis.mu_squared),
        'split_mass_squared_gev2':list(basis.eigen_masses_squared),
        'scalar_mass_response_of_unit_kernel':sensitivity,
        'extra_rgbeta_hh_symbol':extra_name,
        'extra_rgbeta_at_scalar_endpoint':continuous.cjson(extra_rgbeta),
        'extra_sign':sign,'extra_matchete_assumed':continuous.cjson(extra_matchete),
        'assumed_source_sign_not_certified':True,
        'native_matchete_tensor_ratio':continuous.cjson(metric['ratio']),
        'native_tensor_relative_residual':metric['relative_tensor_residual'],
        'anomalous_spectral_phase':continuous.cjson(anomalous_phase(basis.mu_squared,alpha)),
        'normal_one_scale_c5_at_1tev':neutrino.matrix_payload(channels['normal']),
        'anomalous_one_scale_c5_at_1tev':neutrino.matrix_payload(channels['anomalous']),
        'combined_one_scale_c5_at_1tev':neutrino.matrix_payload(channels['combined']),
        'anomalous_over_normal_c5':channels['anomalous_over_normal'],
        'relative_low_energy_mnu_change_vs_normal':relative(m1,m0),
        'normal_one_scale_takagi_masses_ev_descending':normal_low['takagi_masses_ev_descending'],
        'physical_full_general_scalar_matching_certified':False,
        'full_sequential_wilson_matching_performed':False,
        'renormalisation_scheme_compatible_with_other_ten_models_certified':False,
        'scope':('Mass-sensitive legacy Matchete one-scale topology with normal and '
                 'one anomalous self-scalar line and the independent HH quartic. '
                 'RGBeta-Matchete coupling map is conditional. This one-scale result '
                 'MUST NOT be interpreted as physically certified sequential EFT matching '
                 'or directly scheme-equivalent to the other ten candidate results.'),
        **low,
    }
    neutrino.atomic_json(target,out)
    return out


def merge(*,source_root:Path,output_root:Path,
          other_output_id:str='neutrino_v1') -> dict:
    require(set(MODELS).issubset(continuous.MODELS) and len(continuous.MODELS)==16,
            'expected 16 ordinary models')
    rows=[]
    for model in continuous.MODELS:
        tag=continuous.model_tag(model)
        continuous_sha=sha(source_root/tag/'continuous_trajectory.json')
        if model in MODELS:
            path=output_root/f'{tag}_massdependent_c5.json'
            report=neutrino.load(path)
            require(report.get('format')==FORMAT and report.get('status')==STATUS and
                    report.get('input_sha256',{}).get('continuous')==continuous_sha and
                    report.get('physical_full_general_scalar_matching_certified') is False,
                    model+': missing/incompatible mass-dependent one-scale result')
            category='self-bilinear-one-scale-massdependent-conditional'
        else:
            path=source_root/other_output_id/(tag+'_neutrino.json')
            report=neutrino.load(path)
            require(report.get('format')==neutrino.FORMAT and
                    report.get('status') in ('CandidateNoAdditionalQuadraticMatchedC5',
                       'CandidateMassEigenBCEvolvedYukawasSMWeinberg') and
                    report.get('source_continuous_sha256')==continuous_sha and
                    report.get('full_general_scalar_physical_c5_certified') is False,
                    model+': missing/incompatible existing candidate')
            category='existing-candidate-matching-different-scheme'
        masses=report['takagi_masses_ev_descending']
        require(isinstance(masses,list) and len(masses)==3 and
                all(isinstance(x,(int,float)) and math.isfinite(x) and x>=0 for x in masses),
                'invalid neutrino masses '+model)
        rows.append({
           'model':model,'status':report['status'],'matching_category':category,
           'mass1_ev_descending':masses[0],'mass2_ev_descending':masses[1],
           'mass3_ev_descending':masses[2],
           'anomalous_over_normal_c5':report.get('anomalous_over_normal_c5'),
           'source_report':path.relative_to(ROOT).as_posix(),
           'source_sha256':sha(path),
           **report['oscillation_observables_normal_ordering'],
        })
    data={
        'format':'t3-self-bilinear-massdependent-candidate-16-model-summary-v1',
        'status':'16CandidateNeutrinoResultsMixedMatchingApproximations',
        'models_processed':16,'models_new_massdependent':6,
        'results_compatible_as_single_physical_scheme':False,
        'physical_full_general_scalar_matching_certified':False,
        'full_sequential_wilson_matching_performed':False,
        'rows':rows,
    }
    neutrino.atomic_json(output_root/'combined_16_summary.json',data)
    fields=['model','status','matching_category','mass1_ev_descending','mass2_ev_descending',
            'mass3_ev_descending','anomalous_over_normal_c5','delta_m21_sq_ev2',
            'delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23']
    tmp=output_root/'combined_16_summary.csv.tmp'
    with tmp.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore')
        w.writeheader();w.writerows(rows)
    tmp.replace(output_root/'combined_16_summary.csv')
    return data


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    selector=ap.add_mutually_exclusive_group(required=True)
    selector.add_argument('--all',action='store_true')
    selector.add_argument('--model',choices=MODELS)
    ap.add_argument('--run-id',required=True,type=continuous.run_id)
    ap.add_argument('--source-root',type=Path,default=continuous.OUTPUT)
    ap.add_argument('--output-id',default=DEFAULT_OUTPUT_ID)
    ap.add_argument('--extra-sign',type=int,choices=(-1,0,1),default=-1)
    ap.add_argument('--resume',action='store_true')
    ap.add_argument('--combine-only',action='store_true')
    args=ap.parse_args()
    continuous.run_id(args.output_id)
    require(args.output_id not in ('neutrino_v1','self_anomalous_v1','self_bilinear_v1'),
            '--output-id must not overwrite previous calculations')
    folder=args.source_root/args.run_id
    require(folder.is_dir(),'missing completed continuous run '+str(folder))
    out=folder/args.output_id
    selected=MODELS if args.all else (args.model,)
    if args.combine_only:
        require(args.all,'--combine-only requires --all')
    failed=[]
    if not args.combine_only:
        for model in selected:
            try:
                row=calculate(model,source_root=folder,output_root=out,
                              sign=args.extra_sign,resume=args.resume)
                print(f'{model}: mass-dependent anomalous |C_A|/|C_N|='
                      f'{row["anomalous_over_normal_c5"]:.8g}; '
                      f'mass sensitivity={row["scalar_mass_response_of_unit_kernel"]:.8g}; '
                      f'masses={row["takagi_masses_ev_descending"]}',flush=True)
            except (ValueError,TypeError,KeyError,OSError,RuntimeError,
                    FloatingPointError,OverflowError,np.linalg.LinAlgError) as exc:
                failed.append(model)
                print(f'{model}: BLOCKED {exc}',flush=True)
    if failed:
        print('Blocked models:',', '.join(failed))
        return 2
    if args.all:
        combined=merge(source_root=folder,output_root=out)
        print('Combined candidate results:',out/'combined_16_summary.csv',
              'models=',combined['models_processed'],flush=True)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
