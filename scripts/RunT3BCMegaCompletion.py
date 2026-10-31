"""Large local B/C integration: complex heavy RGE covariance + source matching
+ full three-generation SM-only low-energy reference + provenance gates.

This is NOT a physical UV->EFT1->EFT2->SM Wilson matching calculation.
New source reports are never overwritten except via explicit --resume reuse.

From project root:
  python -u -X utf8 -m scripts.RunT3BCMegaCompletion --all
  python -u -X utf8 -m scripts.RunT3BCMegaCompletion --all --resume
"""
from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path
import sys
import traceback

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.experimental import T3BCMegaPhysicalBridge as core
from Numerical.experimental import T3BCComplexUVStress as uv

DEFAULT_OUTPUT=ROOT/'output/studies/t3_bc_mega_completion_v1'
ONE_SCALE=ROOT/'output/studies/t3_bc_running_one_scale_triangle_v1'
FIVE_MASS=ROOT/'output/studies/t3_full_stack_overnight_v1/bc_full_mass'
MIXED_V2=ROOT/'output/studies/t3_bc_mixed_quadratic_kernel_v2'
UV_TOPOLOGY=ROOT/'Lagrangian/interactions/T3Topology.wl'
NATIVE=ROOT/'output/general_t3_stage22/matchete_tensors.json'


def put_new(path:Path,value:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    core.require(not path.exists(),'refusing output overwrite '+str(path))
    tmp=path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(path)


def one_pair_paths(model:str,onescale:Path=ONE_SCALE,five:Path=FIVE_MASS)->tuple[Path,Path]:
    tag=model.replace('-','m')
    return (onescale/f'{tag}_running_one_scale.json',
            five/f'{tag}_uv_to_eft1_coupled_five_mass_v1.json')


def get_kernel(model:str,report_path:Path,kernel_root:Path=MIXED_V2)->dict:
    tag=model.replace('-','m')
    path=kernel_root/f'{tag}_mixed_kernel.json'
    report=core.load(path)
    core.require(report.get('status')=='FullRestrictedSourceKernelChecksPassed' and
                 report.get('model')==model and
                 report.get('physical_C5_certified') is False and
                 report.get('full_sequential_EFT_matching_performed') is False,
                 'mixed kernel source no longer matches restricted scope')
    core.require(report.get('input_running_one_scale_sha256')==core.digest(report_path),
                 'mixed-kernel v2 built from a different source C5')
    nc=report.get('native_quartic_source',{})
    core.require(nc.get('status')=='NativeSevenQuarticsRotatedAtRunningEigenbasis' and
                 nc.get('quartic_count')==7,
                 'native seven-tensor source is missing')
    core.require(nc.get('native_matchete_tensors_sha256')==core.digest(NATIVE),
                 'native source quartics changed since v2')
    source=report.get('full_mass_vs_strict_source_coefficient',{})
    core.require(float(source.get('source_matchete_full_mass_parity_relative',1))<2e-8,
                 'mixed-v2 full-mass source parity failure')
    require_stress=report.get('complex_majorana_spinor_convention_stress',{})
    core.require(require_stress.get('status')=='SyntheticComplexHeavyMassPhaseSensitivityCalculated',
                 'complex fermion phase stress absent in v2')
    return {'status':'VerifiedExistingNativeMixedKernelReportReused',
            'input_sha256':core.digest(path),
            'native_quartics_rotated':nc['quartic_count'],
            'two_contraction_110101':nc['two_independent_110101_contractions'],
            'max_native_quartic_polynomial_relative_error':max(
                x['max_relative_polynomial_and_hc_error']
                for x in nc['source_vertices'].values()),
            'mass_resummation_eigenbasis_relative':report['max_original_basis_eigenbasis_relative_error'],
            'max_mass_insertion_zero_error':report['largest_zero_insertion_kernel_error'],
            'full_vs_strict_C5_relative':source['leading_1_over_MF_restricted_C5_relative_error'],
            'synthetic_spinor_prescription_C5_relative':require_stress['alternative_y1_V_vs_Vconj_relative_C5_difference'],
            'complete_physical_matching':False}


def optional_native_record(model:str,run_file:Path,*,existing_dir:Path=MIXED_V2)->dict:
    try:
        return get_kernel(model,run_file,existing_dir)
    except (ValueError,OSError,KeyError,TypeError) as exc:
        # No physics claim: this does NOT silently substitute a missing
        # source.  Caller records failure to establish native tensor gate.
        return {'status':'NativeMixedSourceProvenanceGateBlocked',
                'reason':str(exc),'complete_physical_matching':False}


def source_hash_gate(running:dict,model:str)->dict:
    """Hash exact archived direct expression, without importing experimental filename."""
    from pathlib import Path
    tag=model.replace('-','m')
    source=ROOT/'output/general_t3_mass_eigen_matchete/one_loop_probe'/(
         tag+'_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    record=ROOT/'output/general_t3_mass_eigen_matchete/one_loop_probe'/(
         tag+'_final_one_loop_probe.json')
    old=running['frozen_source_parity']
    core.require(core.digest(source)==old['native_direct_source_sha256'],
                 'native Matchete one-loop expression changed')
    core.require(core.digest(record)==old['native_direct_matching_report_sha256'],
                 'native Matchete one-loop report changed')
    # The archived evaluation used the exact v22 source-hash inventory,
    # and the complex Wilson coefficient is already verified for that source.
    return {'status':'NativeMatcheteSourceHashUnchanged',
            'direct_expression_sha256':core.digest(source),
            'native_matching_report_sha256':core.digest(record),
            'nonzero_mixed_quadratic_inserted_in_native_matchete_source':False,
            'full_finite_EFT_match':False}


def one_model(model:str,*,output_dir:Path,source_one:Path=ONE_SCALE,
              source_five:Path=FIVE_MASS,source_kernel:Path=MIXED_V2,
              complex_uv:bool=True,proxy:bool=True,rtol:float=1e-8)->dict:
    candidate_path,uv_path=one_pair_paths(model,source_one,source_five)
    candidate=core.load(candidate_path); source=core.load(uv_path)
    evidence,sm=core.check_running_pair(model,source,candidate,
                                       source_path=uv_path,running_path=candidate_path)
    topo=core.check_topology_orientation(UV_TOPOLOGY)
    source_gate=source_hash_gate(candidate,model)
    mixed=optional_native_record(model,candidate_path,existing_dir=source_kernel)
    alg=core.complex_spinor_action_test(seed=17+(0 if model=='B_-1' else 1))
    response={'format':core.FORMAT,'model':model,
              'status':'PartialSourceSupportedStudy',
              'same_scale_matchete_parity':evidence,
              'native_lagrangian_spinor_source':topo,
              'direct_source_provenance':source_gate,
              'mixed_scalar_full_mass_and_native_quartic_gate':mixed,
              'complex_field_redefinition_action':alg,
              'physical_C5_certified':False,
              'full_sequential_EFT_matching_performed':False}
    if complex_uv:
        try:
            response['native_RGBeta_heavy_basis_stress']=uv.run_complex_rge_covariance(
                model,seed=135+(0 if model=='B_-1' else 1),integrate=True,points=5)
        except Exception as e:
            response['native_RGBeta_heavy_basis_stress']={
                'status':'AttemptedNativeComplexRGECovarianceBlocked',
                'error_type':type(e).__name__,'reason':str(e),
                'physical_basis_comparison_not_proven':True}
    else:
        response['native_RGBeta_heavy_basis_stress']={'status':'SkippedByExplicitOption'}
    if proxy:
        a=core.run_sm_only_proxy(model,source,candidate,low_scale=100.,rtol=rtol)
        b=core.run_sm_only_proxy(model,source,candidate,low_scale=100.,rtol=rtol/10)
        A=core.unpack_matrix(a['reference_C5_at_low_gev_inverse'])
        B=core.unpack_matrix(b['reference_C5_at_low_gev_inverse'])
        response['SM_only_low_energy_proxy']=a
        response['SM_only_solver_tolerance_check']={
            'rtol_a':rtol,'rtol_b':rtol/10,
            'C5_relative_difference':core.relative(A,B),
            'nfev_b':b['solver_nfev'],
            'passed_numeric_convergence':core.relative(A,B)<1e-5,
            'full_sequential_matching':False}
        core.require(response['SM_only_solver_tolerance_check']['passed_numeric_convergence'],
                     'SM-only proxy solver not converged')
    else:
        response['SM_only_low_energy_proxy']={'status':'SkippedByExplicitOption'}
    rge_status=response['native_RGBeta_heavy_basis_stress']['status']
    if rge_status=='AttemptedNativeComplexRGECovarianceBlocked':
        response['status']='CompletedMatchingProxyWithNativeComplexUVStressBlocked'
    elif rge_status=='NativeRGBetaComplexHeavyBasisCovarianceComputed' and not response['native_RGBeta_heavy_basis_stress']['source_redefinition_is_covariant_numerically']:
        response['status']='CompletedRunningWithComplexRGBetaBasisCovarianceMismatch'
    elif mixed['status']!='VerifiedExistingNativeMixedKernelReportReused':
        response['status']='CompletedSourceSubsetNativeQuarticProvenanceBlocked'
    elif proxy:
        response['status']='CompletedCandidateIntegrationPhysicalMatchingNotCertified'
    else:
        response['status']='CompletedSourceChecksWithSMProxyExplicitlySkipped'
    return response


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--all',action='store_true')
    ap.add_argument('--model',choices=core.MODELS)
    ap.add_argument('--output-id',default='t3_bc_mega_completion_v1')
    ap.add_argument('--resume',action='store_true')
    ap.add_argument('--skip-complex-uv',action='store_true')
    ap.add_argument('--skip-sm-proxy',action='store_true')
    ap.add_argument('--onescale-dir',type=Path,default=ONE_SCALE)
    ap.add_argument('--five-mass-dir',type=Path,default=FIVE_MASS)
    ap.add_argument('--kernel-v2-dir',type=Path,default=MIXED_V2)
    ap.add_argument('--rtol',type=float,default=1e-8)
    ns=ap.parse_args()
    ap.error('use --all or --model') if not (ns.all or ns.model) or (ns.all and ns.model) else None
    core.require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,70}',ns.output_id)),
                 'unsafe output ID')
    output=ROOT/'output/studies'/ns.output_id
    core.require(not output.exists() or ns.resume,
                 'output already exists; use --resume or a fresh --output-id')
    core.require(1e-11<ns.rtol<1e-4,'invalid tolerance')
    models=core.MODELS if ns.all else (ns.model,)
    results={};errors={}
    for model in models:
        tag=model.replace('-','m');filename=output/f'{tag}_mega.json'
        if ns.resume and filename.is_file():
            row=core.load(filename)
            core.require(row.get('format')==core.FORMAT and row.get('model')==model,
                         'resume report schema mismatch')
            results[model]={'status':row['status'],'report':filename.name,
                            'resumed':True}
            print(model,'RESUME',row['status'],flush=True)
            continue
        try:
            row=one_model(model,output_dir=output,source_one=ns.onescale_dir,
                source_five=ns.five_mass_dir,source_kernel=ns.kernel_v2_dir,
                complex_uv=not ns.skip_complex_uv,proxy=not ns.skip_sm_proxy,
                rtol=ns.rtol)
            put_new(filename,row)
            results[model]={'status':row['status'],'report':filename.name,
                            'resumed':False}
            m=row['SM_only_low_energy_proxy']
            desc=(' low-energy SM-only proxy eV='+str(m['SM_only_100gev_reference']['masses_descending_ev'])
                  if 'SM_only_100gev_reference' in m else '')
            c=row['native_RGBeta_heavy_basis_stress']
            print(model,'PASS',row['status'],desc,
                  'complex-RGE='+c['status'],flush=True)
        except Exception as e:
            errors[model]={'error_type':type(e).__name__,'reason':str(e),
                           'traceback':traceback.format_exc(limit=8)}
            print(model,'BLOCKED',type(e).__name__,str(e),flush=True)
    observed=[r['status'] for r in results.values()]
    if errors:
        campaign_status='PartialWithErrors'
    elif any('Mismatch' in status for status in observed):
        campaign_status='CompletedWithDetectedNativeBetaBasisMismatch'
    elif any('Blocked' in status for status in observed):
        campaign_status='CompletedWithUnresolvedSourceOrNativeBasisBlockers'
    else:
        campaign_status='CompletedConditionalBCTargets'
    summary={'format':'t3-bc-mega-source-and-numerical-batch-v1',
             'status':campaign_status,
             'models':results,'errors':errors,
             'scope':core.blocked_physical_gates(),
             'physical_C5_certified':False,
             'same_trajectory_sequential_EFT_matching_performed':False}
    output.mkdir(parents=True,exist_ok=True)
    # Keep a progress on failure, and a stable final report on completed campaign.
    summary_path=output/'summary.json'
    core.require(ns.resume or not summary_path.exists(), 'refusing to overwrite prior summary')
    tmp=output/'summary.json.tmp'
    tmp.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(summary_path)
    print('SUMMARY',summary['status'],summary_path,flush=True)
    return 0 if not errors else 3


if __name__=='__main__':raise SystemExit(main())
