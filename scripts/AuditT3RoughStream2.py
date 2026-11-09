"""Read-only numerical audit of Stream 2's saved full EFT1 Wilson outputs.

Run from the project root:
    python -m scripts.AuditT3RoughStream2

Does NOT rerun RGEs, alter model files, add a LLSS->C5 log or perform scalar
matching. Checks exact saved arrays, projected C5, 16pi2 beta-channel sums,
trajectory endpoints, same-source IDs and local input/source file SHA256.
Numerical consistency is not an independent certification of the candidate ADM.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path('output/parallel_t3/rough_wilson/rough_fermion_first_v1')
MODELS = tuple(f'{family}_{alpha:+d}' for family, alphas in (
    ('A',(-4,-2,0)), ('B',(-3,-1,1)), ('C',(-3,-1,1)),
    ('D',(-2,0,2)), ('E',(-4,-2,0,2))) for alpha in alphas)
NAMES = ('yukawa_wavefunction','scalar_pair','mixed_yukawa_gauge',
         'crossed_yukawa','conjugate_coefficient_yukawa',
         'scalar_anomalous_dimension','fermion_anomalous_dimension')
CONVENTION = 'C_ijab_full_sum_no_factorial_C5_1over2_H_R_interleaved_v1'
W_FORMAT = 't3-parallel-c-wilson-input-v1'


def slug(model: str) -> str:
    family, a = model.split('_')
    a = int(a)
    return f'{family}_{"m" if a < 0 else "p"}{abs(a)}'


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def near(a, b, *, rtol=1e-8, atol=1e-25) -> bool:
    return bool(np.allclose(a,b,rtol=rtol,atol=atol))


def _max(a):
    return float(np.max(np.abs(a))) if np.size(a) else 0.


def _symmetric(c, *, label):
    scale=max(_max(c),1e-100)
    residual=max(_max(c-c.swapaxes(0,1)),_max(c-c.swapaxes(2,3)))/scale
    if not np.isfinite(residual) or residual>1e-8:
        raise ValueError(f'{label} pair-symmetry relative residual {residual:.3g}')
    return float(residual)


def _indices(meta, nferm, nscalar):
    lp=meta.get('lepton_indices')
    h=meta.get('higgs_source_indices')
    heavy=meta.get('heavy_source_indices')
    if not isinstance(lp,list) or len(lp)!=3 or any(
            not isinstance(pair,list) or len(pair)!=2 for pair in lp):
        raise ValueError('Lepton index pairs absent')
    lf=[int(q) for pair in lp for q in pair]
    if len(set(lf))!=6 or any(x<0 or x>=nferm for x in lf):
        raise ValueError('Invalid lepton indices')
    if h!=[0,1,2,3] or heavy!=list(range(4,nscalar)):
        raise ValueError('Invalid master interleaved scalar indices')
    return lf,[int(pair[0]) for pair in lp],int(h[2]),[int(x) for x in heavy]


def _projection(c, nus, hreal):
    q=4*np.asarray(c[np.ix_(nus,nus,[hreal],[hreal])][:,:,0,0])
    return (q+q.T)/2


def _paths(meta, project_root):
    """Input files recorded by Stream 2: source paths may be absolute on Windows."""
    refs=meta.get('input',{}).get('file_sha256',{})
    if not isinstance(refs,dict) or not refs:
        raise ValueError('Input file SHA256 provenance missing')
    source_manifests=[]
    checked=[]
    for key, value in refs.items():
        if not isinstance(value,dict) or not value.get('sha256') or not value.get('path'):
            raise ValueError(f'Malformed saved source record: {key}')
        p=Path(value['path'])
        if not p.is_absolute():
            p=project_root/p
        if not p.is_file() or sha256(p)!=value['sha256']:
            raise ValueError(f'Stream 1 input missing/changed since computation: {key}: {p}')
        checked.append(key)
        if p.suffix.lower()=='.json':
            source_manifests.append(json.loads(p.read_text(encoding='utf-8-sig')))
    for name in ('matching_context_id','matched_seed_id','lagrangian_source_id',
                 'scalar_basis_id'):
        value=meta.get(name)
        if not isinstance(value,str) or not value.strip():
            raise ValueError(f'Missing inherited Stream 1 identifier: {name}')
        found=[m[name] for m in source_manifests if name in m]
        if not found or any(x!=value for x in found):
            raise ValueError(f'{name} not inherited unchanged from source manifest(s)')
    return checked


def audit_model(path: Path, project_root: Path, *, require_existing_engine=False) -> dict:
    """Audit *one* saved, candidate numerical C interface, with no source mutation."""
    meta_path=path/'wilson.json'
    if not meta_path.is_file():
        raise FileNotFoundError(str(meta_path))
    meta=json.loads(meta_path.read_text(encoding='utf-8-sig'))
    model=meta.get('model')
    if model not in MODELS or path.name!=slug(model):
        raise ValueError('Model/directory mismatch')
    if meta.get('format')!=W_FORMAT or meta.get('tensor_npz')!='wilson_at_MS.npz':
        raise ValueError('Not the Workstream C v1 Wilson contract')
    if meta.get('canonical_wilson_convention')!=CONVENTION:
        raise ValueError('Unknown Wilson/Weinberg normalization')
    if meta.get('scheme')!='MSbar' or meta.get('units')!='GeV':
        raise ValueError('Scheme or units mismatch')
    if meta.get('synthetic_fixture') is not False:
        raise ValueError('Production data marked synthetic or not explicitly marked genuine')
    if meta.get('wilson_evolved_to_ms') is not True or meta.get('includes_eft1_weinberg') is not True:
        raise ValueError('EFT1 evolution missing')
    if meta.get('legacy_direct_llss_to_c5_log_added') is not False:
        raise ValueError('Duplicate legacy LLSS->C5 logarithm indicated')
    method=meta.get('solver',{}).get('method')
    if method!='DOP853':
        raise ValueError(f'Unexpected solver: {method}')
    backend=meta.get('solver',{}).get('backends',[])
    if not backend:
        raise ValueError('Actual numerical ADM backend not identified')
    if require_existing_engine and backend!=['existing_repository_matrix_free_engine']:
        raise ValueError('Expected repository-native matrix-free beta backend')
    if meta.get('adm_variant') not in ('fonseca_eq22_six_term','ward_covariant_candidate'):
        raise ValueError('Unknown/noncovariant ADM selected')
    if meta.get('strict_1loop_endpoint_at_ms') is True:
        raise ValueError('This audit does not validate optional strict 1loop endpoints')
    npz_path=path/'wilson_at_MS.npz'
    if sha256(npz_path)!=meta.get('tensor_npz_sha256'):
        raise ValueError('Evolved NPZ SHA256 mismatch')
    if not near(float(meta.get('mu_f_gev',float('nan'))),1e5,atol=1e-6) or not near(
            float(meta.get('mu_s_gev',float('nan'))),1e3,atol=1e-6):
        raise ValueError('MF/MS do not match 100000/1000 GeV')
    with np.load(npz_path,allow_pickle=False) as raw:
        needed=('C_psi2phi2_ms','C_psi2phi2_tree_mf','C5_EFT1_ms','C5_EFT1_mf',
                'C_LLSS_ms','C_LLSS_mf','beta_C_at_MF_16pi2','beta_C_at_MS_16pi2',
                'trajectory_mu_gev','trajectory_C5_norm_gev_inverse',
                'trajectory_full_norm_gev_inverse')
        if any(k not in raw for k in needed):
            raise ValueError('NPZ missing core full-tensor diagnostic arrays')
        c1=np.asarray(raw['C_psi2phi2_ms'])
        c0=np.asarray(raw['C_psi2phi2_tree_mf'])
        expected_ns=16 if model.startswith('E_') else 12
        if c1.ndim!=4 or c1.shape!=c0.shape or c1.shape[0:2]!=(45,45) or \
                c1.shape[2]!=c1.shape[3] or c1.shape[2]!=expected_ns:
            raise ValueError(f'Wrong full 45-Weyl tensor geometry: {c1.shape}')
        if not all(np.isfinite(v).all() for v in (c0,c1)):
            raise ValueError('Nonfinite Wilson endpoint')
        symmetry=max(_symmetric(c1,label='MS'),_symmetric(c0,label='MF'))
        lf,nus,hr,heavy=_indices(meta,c1.shape[0],c1.shape[2])
        c50,c51=_projection(c0,nus,hr),_projection(c1,nus,hr)
        if not near(c50,raw['C5_EFT1_mf'],atol=1e-25) or not near(
                c51,raw['C5_EFT1_ms'],atol=1e-25):
            raise ValueError('C5 extracted from master tensor disagrees with saved C5')
        ll0=c0[np.ix_(lf,lf,heavy,heavy)]
        ll1=c1[np.ix_(lf,lf,heavy,heavy)]
        if not near(ll0,raw['C_LLSS_mf'],atol=1e-25) or not near(
                ll1,raw['C_LLSS_ms'],atol=1e-25):
            raise ValueError('Saved LLSS differs from full tensor')
        other=c0.copy()
        other[np.ix_(lf,lf,heavy,heavy)]=0
        n_other=float(np.linalg.norm(other))
        beta_max=0.
        for endpoint,full in [('MF',raw['beta_C_at_MF_16pi2']),
                              ('MS',raw['beta_C_at_MS_16pi2'])]:
            if full.shape!=c0.shape or not np.isfinite(full).all():
                raise ValueError(f'Invalid {endpoint} full beta')
            channels=[]
            for name in NAMES:
                key=f'beta_C5_{name}_at_{endpoint}_16pi2'
                if key not in raw or raw[key].shape!=(3,3):
                    raise ValueError(f'Missing projected beta channel: {key}')
                channels.append(np.asarray(raw[key]))
            projected=_projection(full,nus,hr)
            err=_max(np.sum(channels,axis=0)-projected)
            if err>max(1e-26,1e-8*_max(projected)):
                raise ValueError(f'Named beta channels do not sum to full {endpoint} beta: {err}')
            beta_max=max(beta_max,err)
        mu=np.asarray(raw['trajectory_mu_gev'],float)
        trajectory_c5=np.asarray(raw['trajectory_C5_norm_gev_inverse'],float)
        trajectory_full=np.asarray(raw['trajectory_full_norm_gev_inverse'],float)
        if mu.ndim!=1 or len(mu)<2 or any(len(x)!=len(mu) for x in (
                trajectory_c5,trajectory_full)) or not all(np.isfinite(x).all()
                for x in (mu,trajectory_c5,trajectory_full)):
            raise ValueError('Missing/invalid sampled scale trajectory')
        if not (np.all(np.diff(mu)<0) and np.isclose(mu[0],1e5) and
                np.isclose(mu[-1],1e3)):
            raise ValueError('Non-monotone/incorrect MF->MS trajectory')
        if not near(trajectory_c5[[0,-1]],[np.linalg.norm(c50),np.linalg.norm(c51)],
                    atol=1e-25) or not near(
                    trajectory_full[[0,-1]],[np.linalg.norm(c0),np.linalg.norm(c1)],
                    atol=1e-25):
            raise ValueError('Trajectory endpoint does not agree with master tensors')
        n0=float(np.linalg.norm(c50)); n1=float(np.linalg.norm(c51))
        nl0=float(np.linalg.norm(ll0)); nl1=float(np.linalg.norm(ll1))
        atol=float(meta['solver']['atol'])
        significant=n1>max(10*atol,1e-14*float(np.linalg.norm(c0)))
        pure_llss=n_other<=1e-12*max(nl0,1e-100)
        generated=bool(n0<=1e-100 and significant and nl0>0 and pure_llss)
        if bool(meta.get('C5_generated_from_LLSS_only_boundary'))!=generated:
            raise ValueError('Inconsistent generated-from-LLSS assertion')
        if not near(n1,meta['C5_EFT1_ms_norm_gev_inverse'],atol=1e-25):
            raise ValueError('Manifest MS C5 norm differs from tensor')
        if not near(nl0,meta['LLSS_mf_norm_gev_inverse'],atol=1e-25):
            raise ValueError('Manifest MF LLSS norm differs from tensor')
        if not near(n_other,meta['non_LLSS_tree_mf_norm_gev_inverse'],atol=1e-25):
            raise ValueError('Manifest MF non-LLSS norm differs from tensor')
    checked=_paths(meta,project_root)
    code_checked=[]
    hashes=meta.get('local_runtime',{}).get('local_source_sha256',{})
    if hashes is not None:
        if not isinstance(hashes,dict):
            raise ValueError('Malformed recorded runtime source hashes')
        for rel,digest in hashes.items():
            if digest is None:continue  # Module was absent during the original run.
            relpath=Path(rel)
            if relpath.is_absolute() or '..' in relpath.parts:
                raise ValueError('Unsafe recorded runtime source path')
            path=project_root/relpath
            if not path.is_file() or sha256(path)!=digest:
                raise ValueError(f'Runtime ADM/source changed since numerical run: {rel}')
            code_checked.append(rel)
    return {
        'status':'passed_numerical_endpoint_audit',
        'model':model,'adm_variant':meta['adm_variant'],'backend':backend,
        'n_weyl':45,'n_real_scalars':c0.shape[2],
        'C5_MF_norm_GeV_inv':n0,'C5_MS_norm_GeV_inv':n1,
        'LLSS_MF_norm_GeV_inv':nl0,'LLSS_MS_norm_GeV_inv':nl1,
        'non_LLSS_MF_norm_GeV_inv':n_other,
        'generated_from_pure_LLSS':generated,
        'significant_C5_MS':bool(significant),
        'candidate_ADM_not_certified':not bool(meta.get('wilson_adm_certified',False)),
        'pair_symmetry_relative_residual_max':symmetry,
        'beta_channel_sum_max_abs_residual':beta_max,
        'input_files_sha256_rechecked':checked,
        'runtime_code_sha256_rechecked':code_checked,
        'boundary_basis_map':meta.get('input',{}).get('boundary',{}).get('basis_map'),
        'context_ids':{name:meta[name] for name in (
            'matching_context_id','matched_seed_id','lagrangian_source_id','scalar_basis_id')},
        'numerical_trajectory_samples':int(len(mu)),
        'physical_scalar_matching_performed':False,
    }


def run(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=Path('.'))
    p.add_argument('--output-root',type=Path,default=ROOT)
    p.add_argument('--models',nargs='+',choices=MODELS,default=list(MODELS))
    p.add_argument('--require-existing-engine',action='store_true',
                   help='Fail unless the repository-native ADM backend ran')
    p.add_argument('--report',type=Path,
                   help='JSON destination (default: output-root/stream2_audit.json)')
    args=p.parse_args(argv)
    project=args.project_root.resolve()
    root=args.output_root
    if not root.is_absolute():root=project/root
    report=args.report if args.report else root/'stream2_audit.json'
    if not report.is_absolute():report=project/report
    results={}
    for model in args.models:
        try:
            record=audit_model(root/slug(model),project,
                        require_existing_engine=args.require_existing_engine)
            results[model]=record
            print(f'{model}: PASS  ||C5(MS)||={record["C5_MS_norm_GeV_inv"]:.6e} '
                  f'from_pure_LLSS={record["generated_from_pure_LLSS"]} '
                  f'backend={record["backend"]}',flush=True)
        except (ValueError,KeyError,TypeError,OSError,RuntimeError) as e:
            results[model]={'status':'failed_audit','reason':f'{type(e).__name__}: {e}'}
            print(f'{model}: FAIL  {type(e).__name__}: {e}',file=sys.stderr,flush=True)
    success=[v for v in results.values() if v['status']=='passed_numerical_endpoint_audit']
    summary={
        'format':'t3-rough-stream2-endpoint-audit-v1',
        'candidate_not_physical_certification':True,
        'output_root':str(root),
        'selected_models':list(args.models),
        'passed_count':len(success),
        'failed_count':len(results)-len(success),
        'generated_from_pure_LLSS_count':sum(bool(x['generated_from_pure_LLSS']) for x in success),
        'nonzero_C5_MS_models':[x['model'] for x in success if x['significant_C5_MS']],
        'results':results,
    }
    report.parent.mkdir(parents=True,exist_ok=True)
    report.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(f'Audit: {len(success)}/{len(results)} passed; '
          f'{summary["generated_from_pure_LLSS_count"]} generated from pure LLSS. '
          f'Report: {report}',flush=True)
    return 0 if not summary['failed_count'] else 2


if __name__=='__main__':
    raise SystemExit(run())
