"""Read-only muS source-HHSS versus Batch3B candidate comparison (A_-4/D_-2).

Uses existing GeneralT3SequentialLeadingLog.scalar_tensors and source release.
The D_-2 source TensorBundle is expressed in the signed mass-basis release;
undo the source's documented R map before comparing with the fixed native
Batch3B tensor. This rotation is a *candidate basis identification*, not
an independent proof of native CG or physical normalization.

Does not rerun UV/fermion matching, modify RGEs, or overwrite historical output.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import numpy as np

MODELS=('A_-4','D_-2')
SCHEMA='t3-worker-b-source-hhss-comparison-v1'


def fixed_native_hhss(source_hhss, rotation=None):
    """Boundary coordinates b map to native n=R b; Q_native=R Q_boundary R^T."""
    h=np.asarray(source_hhss,float)
    if h.shape!=(4,4,8,8) or not np.isfinite(h).all():
        raise ValueError('source HHSS must be finite (4,4,8,8)')
    if rotation is None:
        return h.copy()
    r=np.asarray(rotation,float)
    if r.shape!=(8,8) or not np.isfinite(r).all() or not np.allclose(r.T@r,np.eye(8),rtol=1e-10,atol=1e-10):
        raise ValueError('source original-from-boundary matrix must be orthogonal (8,8)')
    return np.einsum('ia,ABab,jb->ABij',r,h,r,optimize=True)


def compute_source_hhss(project, model, state, *, release_dir, v4_dir, seed_dir, registry_path):
    """Re-evaluate Matchete native + completed RGBeta portal + extra source HHSS."""
    if model not in MODELS: raise ValueError('unsupported model')
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_release_manifest, load_source
    from Numerical.running.GeneralT3SequentialLeadingLog import scalar_tensors
    from Numerical.lo_certification.couplings.extra_hhss_geometry import model_extra_geometry
    from scripts import RunT3GeneralSequentialFitAblation as original
    from Numerical.joint_scalar_matching.special_mass_v6 import build_special_basis
    registry=json.loads(Path(registry_path).read_text(encoding='utf-8-sig'))
    release=load_release_manifest(release_dir)
    boundary=load_source(release[model],expected_model=model)
    channels=original._load_channels(model,boundary,v4=Path(v4_dir),seeds=Path(seed_dir))
    geometry=model_extra_geometry(model,registry['models'][model],include_components=True)
    # snapshot's source_path is provenance only; running inputs come exclusively
    # from the independently selected muS endpoint state supplied by caller.
    snapshot=original.snapshot(model,dict(state),Path(registry_path))
    bundle=scalar_tensors(model,boundary,channels,snapshot,registry,geometry=geometry)
    h_boundary=np.asarray(bundle.quartic,float)
    if model=='D_-2':
        # This is the exact same basis map used *inside* scalar_tensors.
        basis=build_special_basis(model,snapshot,boundary)
        r=basis['R_original_from_boundary']
        h=fixed_native_hhss(h_boundary,r)
        # Round trip is required; never compare tensors in different bases.
        restored=np.einsum('ia,ABij,jb->ABab',r,h,r,optimize=True)
        roundtrip=float(np.linalg.norm(restored-h_boundary)/max(np.linalg.norm(h_boundary),1e-100))
        if roundtrip>1e-10: raise ValueError('D_-2 source basis roundtrip failed')
        basis_note='source signed mass basis -> original real component basis using R; native/CG certification open'
    else:
        h=fixed_native_hhss(h_boundary)
        roundtrip=0.
        basis_note='source real coordinates equal original nonmixed real components'
    from Numerical.running.T3CoupledQuarticContext import digest
    files={
        'release_boundary':Path(boundary.path),
        'quartic_seed':Path(channels.seed_path),
        'rgbeta_registry':Path(registry_path),
    }
    provenance={label:{'path':str(path),'sha256':digest(path)} for label,path in files.items()}
    return h, {'basis_conversion':basis_note,
               'basis_roundtrip_relative':roundtrip,
               'source_channel_count':len(bundle.names),
               'source_channel_names':list(bundle.names),
               'source_provenance':provenance,
               'source_metadata':bundle.provenance,
               'absolute_physical_normalization_certified':False}


def run_one(project, model, response_root, endpoint_path, release_dir, v4_dir, seed_dir, registry_path):
    from Numerical.running.T3CoupledQuarticContext import CandidateQuarticContext
    from scripts.RunT3WorkerBEndpointAudit import _extract_endpoint
    from scripts.RunT3WorkerBQuarticAudit import audit_tensors
    mu_s,state,state_hash=_extract_endpoint(endpoint_path,model)
    context=CandidateQuarticContext.load(project,response_root,model)
    h, source_meta=compute_source_hhss(project,model,state,release_dir=release_dir,
                                      v4_dir=v4_dir,seed_dir=seed_dir,registry_path=registry_path)
    result=audit_tensors(context,state,hhss=h)
    # HHSS *disagreement* is measured and reported, not equated with invalid
    # tensors: native-to-RGBeta quartic normalization is still uncertified.
    result.update(source_meta)
    result.update({'mu_s_gev':mu_s,'endpoint_file':str(endpoint_path),
                   'endpoint_sha256':state_hash,'source_hhss_evaluated':True,
                   'native_rgbeta_absolute_normalization_certified':False})
    return h,result


def _pairs(items, models):
    out={}
    for item in items:
        model,sep,path=item.partition('=')
        if not sep or model not in models or not path or model in out:
            raise ValueError('use --state-file MODEL=PATH once for each model')
        out[model]=Path(path)
    return out


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--models',default='A_-4,D_-2')
    parser.add_argument('--response-root',type=Path,required=True)
    parser.add_argument('--state-file',action='append',default=[],metavar='MODEL=PATH')
    parser.add_argument('--release-dir',type=Path,default=Path('output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'))
    parser.add_argument('--v4-dir',type=Path,default=Path('output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'))
    parser.add_argument('--seed-dir',type=Path,default=Path('output/full16_symbolic'))
    parser.add_argument('--registry',type=Path,default=Path('RGE/general/generated/general_rgbeta_quartic_registry.json'))
    parser.add_argument('--out-dir',type=Path,default=Path('output/audits/workerb_source_hhss_v1'))
    args=parser.parse_args(argv)
    project=args.project.resolve()
    def absolute(p):return (p if p.is_absolute() else project/p).resolve()
    models=tuple(x.strip() for x in args.models.split(',') if x.strip())
    if not models or len(models)!=len(set(models)) or not set(models)<=set(MODELS):
        parser.error('models must be unique subset of A_-4,D_-2')
    state_paths=_pairs(args.state_file,models)
    if set(state_paths)!=set(models):
        parser.error('supply one --state-file MODEL=PATH for each selected model')
    out=absolute(args.out_dir)
    if out.exists():parser.error(f'refusing to overwrite existing output directory: {out}')
    root=absolute(args.response_root)
    sources=(absolute(args.release_dir),absolute(args.v4_dir),absolute(args.seed_dir),absolute(args.registry))
    for path in [root,*sources,*map(absolute,state_paths.values())]:
        if not path.exists():parser.error(f'missing input: {path}')
    report={'schema':SCHEMA,'status':'Computed_SourceHHSSComparison_PhysicalNormalizationOpen',
            'response_root':str(root),'results':{},'source_hhss_evaluated':True,
            'native_rgbeta_normalization_certified':False,'joint_scalar_matching_performed':False}
    arrays={}
    for model in models:
        path=absolute(state_paths[model]);h,result=run_one(project,model,root,path,*sources)
        arrays[model]=h
        report['results'][model]=result
        if result['status']!='PassedAvailableTensorChecks':
            report['status']='FailedTensorInternalChecks'
    out.mkdir(parents=True,exist_ok=False)
    for model,h in arrays.items():
        np.savez_compressed(out/f'{model}_source_hhss_native.npz',hhss=h)
        report['results'][model]['source_hhss_npz']=str(out/f'{model}_source_hhss_native.npz')
    report_path=out/'report.json'
    report_path.write_text(json.dumps(report,indent=2,allow_nan=False,default=str)+'\n',encoding='utf-8')
    for model,res in report['results'].items():
        print(f'{model}: candidate/source HHSS relative={res["source_hhss_vs_candidate_relative"]:.8e}; '
              f'Bose residual={res["hybrid_bose_symmetry_relative"]:.3e}')
    print('Report:',report_path)
    return 0 if report['status'].startswith('Computed_') else 1

if __name__=='__main__':raise SystemExit(main())
