"""A(-4) diagnostic: MF LLSS Ward projection -> existing EFT1 -> C -> SMEFT.

No UV/Stream1/Stream2 original input file is modified. The projector is an
uncertified source-boundary approximation; all output classifications remain
WardProjectedBoundaryDiagnosticOnly even if C5 at MS is gauge invariant.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

import numpy as np

from Numerical.parallel_t3.wilson.rough_stream2 import load_canonical,write_c_endpoint
from Numerical.parallel_t3.wilson.workflow import run_wilson
from Numerical.parallel_t3.wilson.ward_singlet_v7 import (
    STATUS,project_a4_boundary,sha256,_gauge_norm)
from Numerical.parallel_t3.wilson.adm import beta,total_beta
from Numerical.parallel_t3.threshold.rough_adapter import prepare_pair
from Numerical.parallel_t3.threshold.preflight import preflight_pair
from Numerical.parallel_t3.threshold.rough_numerics import run_rough_case


def _write_json(path:Path,object:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(object,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def execute(project:Path, *, output:Path, steps:int=6, rtol:float=1e-6,
            atol:float=1e-21, overwrite:bool=False, ward_tol:float=1e-7,
            skip_c:bool=False, project_only:bool=False) -> dict:
    project=Path(project).resolve()
    output=Path(output)
    if not output.is_absolute():output=project/output
    # A private result directory, never the genuine all16 benchmark tree.
    if output.exists() and any(output.iterdir()):
        if not overwrite:raise FileExistsError(f'{output}: --overwrite required')
        record=output/'ward_projection_summary.json'
        if not record.is_file():
            raise ValueError(f'Refusing to erase unrecognised output directory: {output}')
        prior=json.loads(record.read_text(encoding='utf-8-sig'))
        if prior.get('format')!='t3-ward-A4-diagnostic-run-v7':
            raise ValueError('Refusing to erase directory with unknown provenance')
    tmp=output.with_name(output.name+'.partial')
    if tmp.exists():
        record=tmp/'ward_projection_summary.json'
        if not record.is_file():
            raise FileExistsError(f'Refusing to erase an unknown partial folder {tmp}')
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    report={'format':'t3-ward-A4-diagnostic-run-v7','model':'A_-4',
            'status':'BLOCKED','classification':STATUS,
            'not_physical_matching':True,
            'original_source_arrays_unchanged_and_projected_copy_isolated':True,
            'unverified_native_scalar_map':True,
            'threshold_order':'UV -> ALL fermions at MF -> BOTH scalars at MS -> SMEFT'}
    _write_json(tmp/'ward_projection_summary.json',report)
    try:
        folder=project/'output/parallel_t3/rough_inputs/rough_fermion_first_v1/a/A_m4'
        case=load_canonical(project,folder,'A_-4',c5_boundary='zero')
        report['source']={name:sha256(folder/name) for name in
            ('background.npz','background.json','boundary.npz','boundary.json')}
        report['original_MF_Ward']=_gauge_norm(case.initial,case.basis)
        original=np.asarray(case.initial,complex)
        restored,pdiag=project_a4_boundary(original)
        report['boundary_projection']=pdiag
        for label,x in [('Yukawa',case.background.y[-1]),
                        ('quartic',case.background.quartic[-1])]:
            # Existing v6 preflight formula, no inferred field transformations.
            from Numerical.parallel_t3.threshold.rough_ward_gate_v6 import gauge_ward_report
            result=gauge_ward_report(x,case.basis,sector='yukawa' if label=='Yukawa' else 'quartic')
            report['MF_'+label+'_Ward']=result['ward_relative']
            if result['passed'] is not True:
                raise ValueError(f'MF {label} background violates gauge Ward: {result["ward_relative"]}')
        b0,backend=beta(restored,*case.background.at(np.log(case.mf)),
                        variant='fonseca_eq22_six_term')
        ward_beta=_gauge_norm(total_beta(b0),case.basis)
        report['MF_projected_boundary_beta_Ward']=ward_beta
        report['beta_backend']=backend
        _write_json(tmp/'ward_projection_summary.json',report)
        if max(ward_beta.values())>ward_tol:
            raise ValueError('Existing EFT1 ADM does not preserve the projected invariant MF boundary: '
                             +repr(ward_beta))
        if project_only:
            report['status']=STATUS
            report['stage']='projected_boundary_and_beta_only'
            _write_json(tmp/'ward_projection_summary.json',report)
            if output.exists():shutil.rmtree(output)
            tmp.rename(output)
            print('A_-4:',STATUS,'projection-only; retained_fraction=',
                  round(pdiag['kept_LLSS_frobenius_fraction'],6))
            return report
        result=run_wilson(case.basis,restored,case.background,case.mf,case.ms,
                          adm_variant='fonseca_eq22_six_term',solver='DOP853',
                          steps=steps,rtol=rtol,atol=atol)
        report['MS_full_Ward']=_gauge_norm(result['C_final'],case.basis)
        report['C5_at_MS_norm_GeV_inverse']=float(np.linalg.norm(result['C5_final']))
        if max(report['MS_full_Ward'].values())>ward_tol:
            raise ValueError('The evolved full EFT1 Wilson has lost gauge invariance: '
                             +repr(report['MS_full_Ward']))
        # Explicitly mark the preprojected source in the original input
        # metadata; these modifications are to an in-memory source record.
        case.source['ward_boundary_projected_diagnostic']=True
        case.source['ward_boundary_projection_record']='ward_projection_summary.json'
        write_c_endpoint(case,result,tmp/'wilson'/'A_m4',
                         adm_variant='fonseca_eq22_six_term',rtol=rtol,
                         atol=atol,steps=steps,overwrite=False)
        wpath=tmp/'wilson'/'A_m4'/'wilson.json'
        wm=json.loads(wpath.read_text(encoding='utf-8'))
        wm['ward_projected_boundary_diagnostic']=True
        wm['output_status']=STATUS
        wm['not_physical_matching']=True
        wm['original_source_boundary_sha256']=report['source']
        wm['original_MF_SU2_Ward_residual']=report['original_MF_Ward']['SU2']
        wm['projected_MF_SU2_Ward_residual']=pdiag['projected_Ward']['SU2']
        wm['boundary_projection_frobenius_fraction_kept']=pdiag['kept_LLSS_frobenius_fraction']
        _write_json(wpath,wm)
        if skip_c:
            report['status']=STATUS
            report['stage']='EFT1_only_C_not_requested'
        else:
            bpath=project/'output/parallel_t3/rough_inputs/rough_fermion_first_v1/c/A_m4/background.json'
            pair=prepare_pair('A_-4',wpath,bpath,tmp/'adapted'/'A_m4')
            pf=preflight_pair(pair.wilson,pair.background,expected_model='A_-4')
            report['C_preflight']=pf
            if pf['status']!='ready_conditional':
                raise ValueError('Existing strict C preflight blocked restored A4: '
                                 +repr(pf['blockers']))
            # Stage C under a plainly marked diagnostic directory. The base
            # helper labels its result based only on LLHH purity, so rewrite
            # the stage manifest BEFORE publishing final diagnostic outputs.
            stage_out=tmp/'C_result_preclassification'/'A_m4'
            info=run_rough_case('A_-4',pair,stage_out,
                                allow_projection_residual=False)
            mpath=stage_out/'manifest.json'; opath=stage_out/'observables.json'
            cm=json.loads(mpath.read_text(encoding='utf-8'))
            co=json.loads(opath.read_text(encoding='utf-8'))
            cm['status']=STATUS;cm['classification']=STATUS
            cm['uncertified_upstream_boundary_Ward_projection']=pdiag
            cm['not_a_genuine_RoughNumericalBenchmark']=True
            co['status']=STATUS;co['not_a_genuine_RoughNumericalBenchmark']=True
            co['boundary_LLSS_orthogonal_projection_kept_fraction']=pdiag['kept_LLSS_frobenius_fraction']
            _write_json(mpath,cm);_write_json(opath,co)
            dest=tmp/'results'/'A_m4';dest.parent.mkdir(parents=True,exist_ok=True)
            stage_out.rename(dest)
            report['masses_eV_DIAGNOSTIC_ONLY']=info['masses_eV']
            report['dm21_sq_eV2_DIAGNOSTIC_ONLY']=info['dm21_sq_eV2']
            report['status']=STATUS
            report['stage']='EFT1_single_MS_matching_SMEFT_low'
        _write_json(tmp/'ward_projection_summary.json',report)
        if output.exists():shutil.rmtree(output)
        tmp.rename(output)
        report['summary']=str(output/'ward_projection_summary.json')
        print('A_-4:',STATUS, 'kept_LLSS_fraction=',round(pdiag['kept_LLSS_frobenius_fraction'],6))
        if 'masses_eV_DIAGNOSTIC_ONLY' in report:
            print('DIAGNOSTIC_ONLY masses_eV=',report['masses_eV_DIAGNOSTIC_ONLY'])
        return report
    except Exception as exc:
        report['blocker']=f'{type(exc).__name__}: {exc}'
        _write_json(tmp/'ward_projection_summary.json',report)
        # Preserve failed diagnostics in a private directory for inspection.
        failed=output.with_name(output.name+'.blocked')
        if failed.exists():
            existing=failed/'ward_projection_summary.json'
            if not existing.is_file():raise
            shutil.rmtree(failed)
        tmp.rename(failed)
        raise


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--output',type=Path,
        default=Path('output/parallel_t3/ward_projected_diagnostic_v7'))
    p.add_argument('--steps',type=int,default=6)
    p.add_argument('--rtol',type=float,default=1e-6)
    p.add_argument('--atol',type=float,default=1e-21)
    p.add_argument('--ward-tolerance',type=float,default=1e-7)
    p.add_argument('--skip-c',action='store_true')
    p.add_argument('--project-only',action='store_true',help='Validate source Ward projection and MF beta only; no ODE')
    p.add_argument('--overwrite',action='store_true')
    args=p.parse_args(argv)
    try:
        execute(args.project_root,output=args.output,steps=args.steps,
                rtol=args.rtol,atol=args.atol,ward_tol=args.ward_tolerance,
                skip_c=args.skip_c,project_only=args.project_only,overwrite=args.overwrite)
    except Exception as exc:
        print('A_-4 Ward-projected diagnostic BLOCKED:',type(exc).__name__,str(exc),file=sys.stderr)
        return 1
    return 0

if __name__=='__main__':sys.exit(main())
