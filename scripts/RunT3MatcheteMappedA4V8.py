"""A(-4): source-exact Matchete T3MixCG triplet frame -> existing EFT1 -> C.

All fermions at MF -> BOTH scalar multiplets in EFT1 -> BOTH at MS -> SMEFT.
This uses an orthogonal map fixed by the archived Matchete CG, NOT a Ward
singlet projection. Other source/spinor/convention assumptions remain explicitly
uncertified. Input files and regular all-16 output directories are unchanged.
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
from Numerical.parallel_t3.wilson.adm import beta,total_beta
from Numerical.parallel_t3.wilson.ward_singlet_v7 import _gauge_norm
from Numerical.parallel_t3.wilson.matchete_a4_map_v8 import (
    CLASSIFICATION, SOURCE_REL, sha256,source_heavy_from_canonical,
    convert_a4_boundary)
from Numerical.parallel_t3.threshold.rough_ward_gate_v6 import gauge_ward_report
from Numerical.parallel_t3.threshold.rough_adapter import prepare_pair
from Numerical.parallel_t3.threshold.preflight import preflight_pair
from Numerical.parallel_t3.threshold.rough_numerics import run_rough_case

REPORT_FORMAT='t3-source-CG-mapped-A4-run-v8'
DEFAULT_OUTPUT=Path('output/parallel_t3/source_cg_mapped_a4_v8')


def _json(path:Path,record:dict):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def execute(project:Path,*,output:Path=DEFAULT_OUTPUT,boundary_only:bool=False,
            skip_c:bool=False,overwrite:bool=False,steps:int=6,rtol:float=1e-6,
            atol:float=1e-21,ward_tolerance:float=1e-7) -> dict:
    project=Path(project).resolve();output=Path(output)
    if not output.is_absolute():output=project/output
    # Never allow writing into standard genuine or old projected result folders.
    forbidden=[project/'output/parallel_t3/rough_results',
               project/'output/parallel_t3/rough_wilson',
               project/'output/parallel_t3/ward_projected_diagnostic_v7',
               project/'output/parallel_t3/rough_inputs']
    if any(output==p or p in output.parents for p in forbidden):
        raise ValueError('v8 must use its own source-CG-mapped diagnostic output')
    record_path=output/'source_cg_map_summary.json'
    if output.exists() and any(output.iterdir()):
        if not overwrite:raise FileExistsError(f'{output}: --overwrite required')
        if not record_path.is_file() or json.loads(record_path.read_text()).get('format')!=REPORT_FORMAT:
            raise ValueError('Refusing to overwrite an unrecognised output folder')
    tmp=output.with_name(output.name+'.partial')
    if tmp.exists():
        fp=tmp/'source_cg_map_summary.json'
        if not fp.is_file() or json.loads(fp.read_text()).get('format')!=REPORT_FORMAT:
            raise FileExistsError('Refusing to delete unknown partial output')
        shutil.rmtree(tmp)
    tmp.parent.mkdir(parents=True,exist_ok=True);tmp.mkdir()
    report={'format':REPORT_FORMAT,'model':'A_-4','status':'BLOCKED',
            'classification':CLASSIFICATION,'not_physical_matching':True,
            'numerical_source_cg_map_verified':False,
            'EFT_order':'UV -> all fermions out at MF -> BOTH scalars out at MS -> SMEFT',
            'original_source_files_not_modified':True,
            'no_orthogonal_singlet_projection':True}
    _json(tmp/'source_cg_map_summary.json',report)
    try:
        folder=project/'output/parallel_t3/rough_inputs/rough_fermion_first_v1/a/A_m4'
        case=load_canonical(project,folder,'A_-4',c5_boundary='zero')
        source_file=project/SOURCE_REL
        meta=json.loads((folder/'boundary.json').read_text(encoding='utf-8-sig'))
        native_declared=str(meta.get('native_source_path','')).replace('\\','/')
        if native_declared!=SOURCE_REL:
            raise ValueError('Original Stream1 manifest does not name the corrected A(-4) Matchete source')
        if meta.get('source_sha256',{}).get(SOURCE_REL)!=sha256(source_file):
            raise ValueError('Source seed does not match Stream1 boundary manifest SHA256')
        report['source']={k:sha256(folder/k) for k in
                          ('background.npz','background.json','boundary.npz','boundary.json')}
        R,source=source_heavy_from_canonical(source_file)
        report['source_CG']=source
        with np.load(folder/'boundary.npz',allow_pickle=False) as d:
            if 'C_LLSS_formal_muF' not in d:
                raise ValueError('Require original source-resolved Matchete formal LLSS array')
            raw=np.asarray(d['C_LLSS_formal_muF'],complex).copy()
            oldR=np.asarray(d['heavy_source_from_canonical']) if 'heavy_source_from_canonical' in d else None
            if oldR is None or oldR.shape!=(8,8) or not np.allclose(oldR,np.eye(8),atol=1e-12):
                raise ValueError('Original Stream1 baseline was not the explicitly assumed identity map')
        fac=meta.get('gamma_cc_to_weyl_signed_factor')
        if fac is None:raise ValueError('Missing explicit original GammaCC signed factor')
        if isinstance(fac,dict) and set(fac)=={'real','imag'}:
            fac=complex(fac['real'],fac['imag'])
        mapped,diagnosis=convert_a4_boundary(case.initial,raw,complex(fac),R,
                                           ward_tolerance=ward_tolerance)
        report['MF_source_boundary_map']=diagnosis
        for label,arr,sector in [
                ('Yukawa',case.background.y[-1],'yukawa'),
                ('quartic',case.background.quartic[-1],'quartic')]:
            ward=gauge_ward_report(arr,case.basis,sector=sector)
            report['MF_'+label+'_Ward']=ward['ward_relative']
            if not ward['passed']:
                raise ValueError(f'MF {label} source background violates gauge Ward: {ward["ward_relative"]}')
        b0,backend=beta(mapped,*case.background.at(np.log(case.mf)),
                        variant='fonseca_eq22_six_term')
        report['MF_mapped_boundary_beta_Ward']=_gauge_norm(total_beta(b0),case.basis)
        report['MF_beta_backend']=backend
        _json(tmp/'source_cg_map_summary.json',report)
        if max(report['MF_mapped_boundary_beta_Ward'].values())>ward_tolerance:
            raise ValueError('Candidate ADM violates gauge Ward on source-CG-mapped boundary')
        report['numerical_source_cg_map_verified']=True
        if boundary_only:
            report['status']=CLASSIFICATION
            report['stage']='source_CG_boundary_and_MF_beta_only'
        else:
            result=run_wilson(case.basis,mapped,case.background,case.mf,case.ms,
                              adm_variant='fonseca_eq22_six_term',solver='DOP853',
                              steps=steps,rtol=rtol,atol=atol)
            report['MS_full_Ward']=_gauge_norm(result['C_final'],case.basis)
            if max(report['MS_full_Ward'].values())>ward_tolerance:
                raise ValueError('EFT1 evolved Wilson violates gauge Ward')
            case.source['Matchete_A4_source_derived_triplet_map']=source
            case.source['not_physical_matching']=True
            write_c_endpoint(case,result,tmp/'wilson'/'A_m4',
                             adm_variant='fonseca_eq22_six_term',rtol=rtol,
                             atol=atol,steps=steps,overwrite=False)
            wpath=tmp/'wilson'/'A_m4'/'wilson.json'
            wm=json.loads(wpath.read_text(encoding='utf-8'))
            wm.update({'output_status':CLASSIFICATION,
                       'source_CG_derived_triplet_map_verified':True,
                       'source_CG_evidence':source,
                       'not_physical_matching':True})
            _json(wpath,wm)
            if skip_c:
                report['status']=CLASSIFICATION
                report['stage']='EFT1_evolved_only'
            else:
                bpath=project/'output/parallel_t3/rough_inputs/rough_fermion_first_v1/c/A_m4/background.json'
                pair=prepare_pair('A_-4',wpath,bpath,tmp/'adapted'/'A_m4')
                pf=preflight_pair(pair.wilson,pair.background,expected_model='A_-4')
                report['C_preflight']=pf
                if pf['status']!='ready_conditional':
                    raise ValueError('Strict scalar-threshold preflight failed: '+repr(pf['blockers']))
                stage=tmp/'C_unclassified'/'A_m4'
                info=run_rough_case('A_-4',pair,stage,allow_projection_residual=False)
                manifest=stage/'manifest.json';observables=stage/'observables.json'
                m=json.loads(manifest.read_text());o=json.loads(observables.read_text())
                m.update({'status':CLASSIFICATION,'classification':CLASSIFICATION,
                          'source_CG_map':source,'not_a_certified_prediction':True})
                o.update({'status':CLASSIFICATION,'not_a_certified_prediction':True})
                _json(manifest,m);_json(observables,o)
                dest=tmp/'results'/'A_m4';dest.parent.mkdir(parents=True,exist_ok=True)
                stage.rename(dest)
                report['masses_eV_CONDITIONAL']=info['masses_eV']
                report['dm21_sq_eV2_CONDITIONAL']=info['dm21_sq_eV2']
                report['status']=CLASSIFICATION
                report['stage']='EFT1_one_MS_scalar_matching_SMEFT_low'
        _json(tmp/'source_cg_map_summary.json',report)
        if output.exists():shutil.rmtree(output)
        tmp.rename(output)
        report['summary_path']=str(output/'source_cg_map_summary.json')
        print('A_-4:',CLASSIFICATION,'stage=',report['stage'],
              'source_LLSS_fraction_retained=',round(diagnosis['LLSS_norm_fraction_retained'],10),
              'mapped_SU2_Ward=',format(diagnosis['source_cg_mapped_ward']['SU2'],'.5g'))
        if 'masses_eV_CONDITIONAL' in report:
            print('SOURCE_CG_MAPPED_CONDITIONAL masses_eV=',report['masses_eV_CONDITIONAL'])
        return report
    except Exception as exc:
        report['blocker']=f'{type(exc).__name__}: {exc}'
        _json(tmp/'source_cg_map_summary.json',report)
        failed=output.with_name(output.name+'.blocked')
        if failed.exists():
            check=failed/'source_cg_map_summary.json'
            if not check.is_file() or json.loads(check.read_text()).get('format')!=REPORT_FORMAT:
                raise
            shutil.rmtree(failed)
        tmp.rename(failed)
        raise


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    p.add_argument('--boundary-only',action='store_true',
                   help='Source CG basis map and MF beta, no EFT1 ODE or scalar matching')
    p.add_argument('--skip-c',action='store_true',help='Evolve EFT1 only')
    p.add_argument('--overwrite',action='store_true')
    p.add_argument('--steps',type=int,default=6)
    p.add_argument('--rtol',type=float,default=1e-6)
    p.add_argument('--atol',type=float,default=1e-21)
    p.add_argument('--ward-tolerance',type=float,default=1e-7)
    args=p.parse_args(argv)
    try:
        execute(args.project_root,output=args.output,boundary_only=args.boundary_only,
                skip_c=args.skip_c,overwrite=args.overwrite,steps=args.steps,
                rtol=args.rtol,atol=args.atol,ward_tolerance=args.ward_tolerance)
    except Exception as e:
        print('A(-4) source CG map BLOCKED:',type(e).__name__,str(e),file=sys.stderr)
        return 1
    return 0

if __name__=='__main__':sys.exit(main())
