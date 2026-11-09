"""One guarded all-16 T3 source-CG rough numerical benchmark (conditional).

UV -> ALL fermions out at common MF -> EFT1 with BOTH scalars and full
coupled Wilson -> BOTH scalars out at common MS -> SMEFT -> observables.

Creates a separate, source-audited Stream1 input set from saved Workstream B
and one batch output; never rewrites the original sources/rough results.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

import numpy as np

from Numerical.parallel_t3.rough_inputs import builder
from Numerical.parallel_t3.background.contracts import tag
from Numerical.parallel_t3.wilson.basis import ALL_MODELS
from Numerical.parallel_t3.wilson.coordinates import embed_formal_llss
from Numerical.parallel_t3.wilson.rough_stream2 import load_canonical, write_c_endpoint
from Numerical.parallel_t3.wilson.workflow import run_wilson
from Numerical.parallel_t3.wilson.adm import beta, total_beta
from Numerical.parallel_t3.wilson.ward_singlet_v7 import _gauge_norm
from Numerical.parallel_t3.wilson.matchete_a4_map_v8 import (
    convert_a4_boundary, SOURCE_REL)
from Numerical.parallel_t3.wilson.source_audited_all16 import (
    _BASE_QUARTIC_UNITS, assert_source_rgbeta_definition, hash_file,
    source_invariant_units, source_scalar_map, mass_ward_report,
    validated_wilson_map,covariance_corrected_mass)
from Numerical.parallel_t3.threshold.rough_ward_gate_v6 import gauge_ward_report
from Numerical.parallel_t3.threshold.rough_adapter import prepare_pair
from Numerical.parallel_t3.threshold.preflight import preflight_pair
from Numerical.parallel_t3.threshold.rough_numerics import run_rough_case

FORMAT = 't3-source-audited-conditional-all16-v1'
DELIVERY = 'source-corrected-all16-batch-v2'
RUN_ID = 'source_audited_all16_v1'
OUTPUT = Path('output/parallel_t3/source_audited_all16_v1')
INPUT = Path('output/parallel_t3/rough_inputs') / RUN_ID
CLASSIFICATION = 'SourceCGAndInvariantMappedNumericalCandidateNotCertified'


def _json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False,sort_keys=True)+'\n',encoding='utf-8')


def _npz(path: Path, **arrays) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.partial.npz')
    np.savez_compressed(tmp,**arrays)
    tmp.replace(path)


def _load_npz(path: Path) -> dict:
    with np.load(path,allow_pickle=False) as d:
        return {k:np.array(d[k]) for k in d.files}


def _update_boundary(project: Path, folder: Path, model: str,
                     tol: float) -> tuple[dict,dict]:
    """Build source-CG map, norm-check source LLSS, atomically save new input."""
    old = load_canonical(project, folder, model, c5_boundary='zero')
    bp=folder/'boundary.json'; npz=folder/'boundary.npz'
    meta=json.loads(bp.read_text(encoding='utf-8-sig'))
    raw_data=_load_npz(npz)
    raw=raw_data['C_LLSS_formal_muF']
    native=meta.get('native_source_path','').replace('\\','/')
    if not native or Path(native).is_absolute() or '..' in Path(native).parts:
        raise ValueError('Missing/unsafe original Matchete source path')
    source=project/native
    if meta.get('source_sha256',{}).get(native) != hash_file(source):
        raise ValueError('Original Matchete source changed relative to Stream1 manifest')
    r, evidence=source_scalar_map(project, model, source)
    fac=meta.get('gamma_cc_to_weyl_signed_factor')
    if fac is None:
        raise ValueError('Missing explicit signed GammaCC-to-Weyl benchmark factor')
    if isinstance(fac,dict):
        fac=complex(fac['real'],fac['imag'])
    c, report=validated_wilson_map(old.basis,old.initial,raw,complex(fac),r,tol)
    if model == 'A_-4':
        # Actual regression to v8 source-CG conversion, not merely its expected mass.
        if native != SOURCE_REL:
            raise ValueError('A(-4) source differs from proven v8 native seed')
        v8, v8_report=convert_a4_boundary(old.initial,raw,complex(fac),r,
                                         ward_tolerance=tol)
        if not np.allclose(c,v8,atol=1e-24,rtol=1e-12):
            raise ValueError('A(-4) source bridge deviated from verified v8 tensor')
        report['A4_v8_regression']={
            'source_LLSS_fraction_retained':v8_report['LLSS_norm_fraction_retained'],
            'source_cg_mapped_ward':v8_report['source_cg_mapped_ward']}
    # The source-referenced Matchete formal LLSS stays UNCHANGED. Only its
    # explicit real orthogonal map and the clearly labelled spare are replaced.
    if not np.allclose(raw_data['heavy_source_from_canonical'],
                       np.eye(len(r)),atol=1e-12):
        raise ValueError('New isolated Stream1 source unexpectedly non-identity')
    raw_data['heavy_source_from_canonical']=r
    raw_data['C_full_provisional']=c
    _npz(npz,**raw_data)
    meta['npz_sha256']=hash_file(npz)
    meta['source_audited_heavy_map']=evidence
    meta['source_audited_orthogonal_not_projected']=True
    meta['Matchete_RGBeta_scalar_map_certified']=False
    _json(bp,meta)
    # This input has been rebuilt from native source and *conditionally* correct
    # scalar components. Scalar mass and quartic are already in the canonical
    # highest-weight frame -- never apply r a second time to these tensors.
    assumptions=folder/'assumptions.json'
    assumptions_doc=json.loads(assumptions.read_text(encoding='utf-8-sig'))
    assumptions_doc.update({
        'source_CG_Wilson_map':evidence,
        'corrected_lambdaT3_source_registered_invariant':True,
        'corrected_E_lambda12Cross_same_variance_triplet_metric':True,
        'native_RGBeta_coefficient_normalisations_certified':False,
        'scalar_mass_background_already_canonical_and_NOT_double_rotated':True,
    })
    _json(assumptions,assumptions_doc)
    return evidence, report


def _update_mass(project:Path,model:str)->dict:
    """Correct only the canonical scalar mass input for the common MS gate.

    Source B trajectories and existing old C inputs stay unchanged; this
    updates the isolated v1-run copy in place with updated SHA256 manifest.
    """
    folder=project/INPUT/'c'/tag(model)
    source=folder/'background.npz'
    meta_path=folder/'background.json'
    arrays=_load_npz(source)
    original=np.asarray(arrays['mass_squared_gev2_ms'],float)
    converted,report=covariance_corrected_mass(project,model,original)
    arrays['mass_squared_gev2_ms']=converted
    _npz(source,**arrays)
    meta=json.loads(meta_path.read_text(encoding='utf-8-sig'))
    meta['tensor_npz_sha256']=hash_file(source)
    meta['source_audited_mass_transform']=report
    meta['scalar_eigenvalues_squared_gev2']=report['corrected_mass_eigenvalues_sq_gev2']
    meta['gauge_ward_certified']=False
    _json(meta_path,meta)
    return report


def _preflight_source_case(project:Path, model:str, folder:Path,
                           tol:float, *, unit_source:dict) -> tuple[object,dict]:
    case=load_canonical(project,folder,model,c5_boundary='zero')
    b=case.basis
    report={'source_invariant_definition':unit_source}
    report['MF_Wilson_Ward']=_gauge_norm(case.initial,b)
    if max(report['MF_Wilson_Ward'].values())>tol:
        raise ValueError('MF source Wilson not gauge invariant: '+repr(report['MF_Wilson_Ward']))
    # The underlying native trajectories remain those saved by Workstream B.
    for label,tensor,sector in (
        ('MF_Yukawa',case.background.y[-1],'yukawa'),
        ('MS_Yukawa',case.background.y[0],'yukawa'),
        ('MF_quartic',case.background.quartic[-1],'quartic'),
        ('MS_quartic',case.background.quartic[0],'quartic')):
        w=gauge_ward_report(tensor,b,sector=sector)
        report[label+'_Ward']=w['ward_relative']
        if not w['passed'] or w['max_relative']>tol:
            raise ValueError(f'{label} source component gauge Ward failed: {w["ward_relative"]}')
    # Check full coupled scalar quadratic mixing, including B/C special X/Y.
    sc=project/INPUT/'c'/tag(model)/'background.npz'
    if not sc.is_file():
        raise ValueError('Common-scalar matching background missing')
    mass=_load_npz(sc)['mass_squared_gev2_ms']
    mw=mass_ward_report(mass,b)
    report['MS_full_heavy_mass_Ward']=mw
    if max(mw.values())>tol:
        raise ValueError('MS scalar mass mixing not canonical gauge invariant: '+repr(mw))
    if np.min(np.linalg.eigvalsh(mass)) <= 0:
        raise ValueError('MS scalar mass matrix has nonpositive mass-squared eigenvalues')
    b0,backend=beta(case.initial,*case.background.at(np.log(case.mf)),
                    variant='fonseca_eq22_six_term')
    report['MF_beta_backend']=backend
    report['MF_beta_Ward']=_gauge_norm(total_beta(b0),b)
    if max(report['MF_beta_Ward'].values())>tol:
        raise ValueError('Candidate full ADM fails source-corrected MF beta Ward')
    return case,report


def run_one(project:Path, out:Path, model:str, *, tol:float,steps:int,
            rtol:float,atol:float,gamma_factor:float,skip_ode:bool=False) -> dict:
    item={'model':model,'slug':tag(model),'status':'BLOCKED',
          'classification':'BlockedUnresolvedNativeSourceConvention',
          'source_map':'BLOCKED','boundary_Ward':'BLOCKED','quartic_Ward':'BLOCKED',
          'EFT1':'NOT_RUN','scalar_matching':'NOT_RUN','neutrino_masses':'NOT_RUN'}
    folder=project/INPUT/'a'/tag(model)
    try:
        opt=builder.Options(root=project,B_run_id=builder.B_RUN,out_run_id=RUN_ID,
                            gamma_factor=gamma_factor, assume_formal_identity=True,
                            assume_sm_weyl_yukawa=True,reduced_potential=True,
                            hard_c5_zero_benchmark=True,overwrite=True)
        # Local in-process injection into the EXISTING Stream1 builder. Not a
        # new RGE, not a replacement of the saved Workstream B trajectory.
        original=builder.quartic_unit_tensors
        if original is not _BASE_QUARTIC_UNITS:
            raise ValueError('Unexpected runtime quartic builder already modified')
        try:
            builder.quartic_unit_tensors=source_invariant_units
            built=builder.build_model(model,opt)
        finally:
            builder.quartic_unit_tensors=original
        item['source_builder']=built['status']
        src,correction=_update_boundary(project,folder,model,tol)
        item['source_map']='VERIFIED_ORTHOGONAL_OR_WARD_GATED_IDENTITY'
        item['boundary_Ward']=correction['ward']
        item['source_CG']=src
        item['source_boundary_norm']=correction
        item['scalar_mass_transform']=_update_mass(project,model)
        case,pre=_preflight_source_case(project,model,folder,tol,
                   unit_source=assert_source_rgbeta_definition(project))
        item['quartic_Ward']={
            'MF':pre['MF_quartic_Ward'],'MS':pre['MS_quartic_Ward']}
        item['preflight_Ward']=pre
        item['stage']='SOURCE_AND_MF_VALIDATED'
        if skip_ode:
            item.update(status='BOUNDARY_ONLY',classification=CLASSIFICATION)
            return item
        result=run_wilson(case.basis,case.initial,case.background,
                          case.mf,case.ms,adm_variant='fonseca_eq22_six_term',
                          solver='DOP853',steps=steps,rtol=rtol,atol=atol)
        item['MS_Wilson_Ward']=_gauge_norm(result['C_final'],case.basis)
        if max(item['MS_Wilson_Ward'].values())>tol:
            raise ValueError('Full coupled EFT1 Wilson fails MS Ward')
        item['EFT1']='PASS_FULL_COUPLED_DOP853'
        case.source['source_CG_evidence']=src
        case.source['RGBeta_invariant_status']='conditional_native_normalisation'
        dest=out/'wilson'/tag(model)
        write_c_endpoint(case,result,dest,adm_variant='fonseca_eq22_six_term',
                         rtol=rtol,atol=atol,steps=steps,overwrite=False)
        wpath=dest/'wilson.json'
        wm=json.loads(wpath.read_text(encoding='utf-8'))
        wm.update({'output_status':CLASSIFICATION,'source_CG_map':src,
                   'not_a_certified_prediction':True})
        _json(wpath,wm)
        bpath=project/INPUT/'c'/tag(model)/'background.json'
        pair=prepare_pair(model,wpath,bpath,out/'adapted'/tag(model))
        pf=preflight_pair(pair.wilson,pair.background,expected_model=model)
        item['C_preflight']=pf
        if pf['status']!='ready_conditional':
            raise ValueError('Common MS scalar matching preflight: '+repr(pf['blockers']))
        stage=out/'_unclassified'/tag(model)
        info=run_rough_case(model,pair,stage,allow_projection_residual=False)
        manifest=stage/'manifest.json'
        observables=stage/'observables.json'
        m=json.loads(manifest.read_text(encoding='utf-8'))
        o=json.loads(observables.read_text(encoding='utf-8'))
        m.update({'status':CLASSIFICATION,'classification':CLASSIFICATION,
                  'source_CG_map':src,'not_a_certified_prediction':True})
        o.update({'status':CLASSIFICATION,'not_a_certified_prediction':True})
        _json(manifest,m);_json(observables,o)
        target=out/'results'/tag(model)
        target.parent.mkdir(parents=True,exist_ok=True)
        stage.rename(target)
        item['scalar_matching']='PASS_ONE_COMMON_MS_COMPLETE_MASS_MATRIX'
        item['neutrino_masses']='CONDITIONAL_NUMERICAL_ONLY'
        item['masses_eV_CONDITIONAL']=info['masses_eV']
        item['dm21_sq_eV2_CONDITIONAL']=info['dm21_sq_eV2']
        item['classification']=CLASSIFICATION
        item['status']='NUMERICAL_CANDIDATE_NOT_CERTIFIED'
        item['stage']='FULL_EFT_CHAIN'
    except Exception as exc:
        item['blocker']=f'{type(exc).__name__}: {exc}'
        item['blocked_at_stage']=item.get('stage','SOURCE_INPUT_OR_MAPPING')
        item['status']='BLOCKED'
        print(model,'BLOCKED',item['blocker'],file=sys.stderr)
    return item


def execute(project:Path,*,output:Path=OUTPUT,models=ALL_MODELS,
            overwrite:bool=False,boundary_only:bool=False,steps:int=6,
            rtol:float=1e-6,atol:float=1e-21,ward_tol:float=1e-7,
            gamma_factor:float=1.,resume:bool=False)->dict:
    project=Path(project).resolve()
    out=output if output.is_absolute() else project/output
    if out.resolve()==project or not (project/'Numerical/parallel_t3/rough_inputs/builder.py').is_file():
        raise ValueError('Project root / output path is not a valid local checkout')
    if overwrite and resume:
        raise ValueError('--overwrite and --resume cannot be combined')
    present=out.exists() and any(out.iterdir())
    if present and not (overwrite or resume):
        raise FileExistsError('Output exists; pass --resume to retry blocked models, or --overwrite to rerun all')
    old_summary=None
    if present:
        check=out/'summary.json'
        if not check.is_file():
            raise FileExistsError('Existing output lacks summary; cannot safely resume/overwrite')
        old_summary=json.loads(check.read_text(encoding='utf-8'))
        if old_summary.get('format')!=FORMAT:
            raise FileExistsError('Refusing an unrecognised summary output directory')
    if overwrite and present:
        shutil.rmtree(out)
    out.mkdir(parents=True,exist_ok=True)
    summary=(old_summary if resume and old_summary else {
             'format':FORMAT,'classification':CLASSIFICATION,
             'run_id':RUN_ID,'mf_gev':100000.,'ms_gev':1000.,
             'models_requested':list(models),'results':{},
             'eft_order':'UV -> ALL fermions -> EFT1 with BOTH scalars -> ALL scalars -> SMEFT',
             'candidate_ADM':'fonseca_eq22_six_term',
             'full_NPZ_outputs_from_genuine_saved_Workstream_B':True,
             'source_mapping_convention_status':'CG-derived, RGBeta quartic coefficients conditional',
             'gamma_cc_to_weyl_factor_certified':False,
             'scalar_coupling_normalisations_certified':False,
             'candidate_ADM_certified':False,
             'not_certified_physical_predictions':True})
    summary['delivery_revision']=DELIVERY
    summary['resume_last_requested']=bool(resume)
    summary['models_requested']=list(ALL_MODELS) if resume else list(models)
    _json(out/'summary.json',summary)
    if resume:
        selected=[model for model in models if summary['results'].get(model,{}).get('status')
                  not in ('NUMERICAL_CANDIDATE_NOT_CERTIFIED',)]
        print('Resume selected:',', '.join(selected) if selected else '(none)')
    else:
        selected=list(models)
    for model in selected:
        if resume:
            # Clear only this model's previous partial outputs. Completed models
            # are never touched; rerunning a blocked model is deterministic.
            for kind in ('wilson','adapted','results','_unclassified'):
                owned=out/kind/tag(model)
                if owned.is_dir():shutil.rmtree(owned)
        summary['results'][model]=run_one(project,out,model,tol=ward_tol,
            steps=steps,rtol=rtol,atol=atol,gamma_factor=gamma_factor,
            skip_ode=boundary_only)
        _json(out/'summary.json',summary)
        r=summary['results'][model]
        if 'masses_eV_CONDITIONAL' in r:
            print(model,r['status'],'masses_eV:',r['masses_eV_CONDITIONAL'])
        elif r['status']=='BOUNDARY_ONLY':
            print(model,'MF SOURCE AND WARD PASS')
    summary['executed_this_call']=list(selected)
    summary['complete_count']=sum(v['status']=='NUMERICAL_CANDIDATE_NOT_CERTIFIED'
                                  for v in summary['results'].values())
    summary['blocked_count']=sum(v['status']=='BLOCKED' for v in summary['results'].values())
    summary['boundary_only_count']=sum(v['status']=='BOUNDARY_ONLY'
                                    for v in summary['results'].values())
    _json(out/'summary.json',summary)
    return summary


def show(output:Path) -> None:
    data=json.loads((output/'summary.json').read_text(encoding='utf-8'))
    if data.get('format')!=FORMAT:
        raise ValueError('Unrecognized all-16 summary')
    print('Model   Status                              Conditional masses (eV) or exact blocker')
    print('-'*120)
    for m,record in data['results'].items():
        if 'masses_eV_CONDITIONAL' in record:
            message=', '.join(f'{x:.8g}' for x in record['masses_eV_CONDITIONAL'])
        else:
            message=record.get('blocker',record.get('stage','not run'))
        print(f'{m:7} {record["status"]:35} {message}')
    print('Conditional numerical, NOT certified physical masses.')


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--output',type=Path,default=OUTPUT)
    p.add_argument('--overwrite',action='store_true')
    p.add_argument('--resume',action='store_true',help='Only retry blocked/missing models from existing batch')
    p.add_argument('--boundary-only',action='store_true')
    p.add_argument('--show',action='store_true')
    p.add_argument('--model',choices=ALL_MODELS,help='Optional isolated regression case')
    p.add_argument('--steps',type=int,default=6)
    p.add_argument('--rtol',type=float,default=1e-6)
    p.add_argument('--atol',type=float,default=1e-21)
    p.add_argument('--ward-tolerance',type=float,default=1e-7)
    p.add_argument('--gamma-factor',type=float,default=1.)
    args=p.parse_args(argv)
    try:
        output=args.output if args.output.is_absolute() else args.project_root/args.output
        if args.show:
            show(output)
            return 0
        summary=execute(args.project_root,output=args.output,
                        models=[args.model] if args.model else ALL_MODELS,
                        overwrite=args.overwrite,boundary_only=args.boundary_only,
                        steps=args.steps,rtol=args.rtol,atol=args.atol,
                        ward_tol=args.ward_tolerance,gamma_factor=args.gamma_factor,
                        resume=args.resume)
        print('Batch summary:',output/'summary.json',
              'complete=',summary['complete_count'],
              'blocked=',summary['blocked_count'])
        return 0 if summary['blocked_count']==0 else 2
    except Exception as e:
        print('ALL16 SOURCE-AUDITED RUNNER ERROR:',type(e).__name__,str(e),file=sys.stderr)
        return 1

if __name__=='__main__':
    raise SystemExit(main())
