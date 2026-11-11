"""Read-only Worker B cached-quartic endpoint audit for A_-4 and D_-2.

This audits *existing* Batch3B affine response data; no invariant is inferred.
The independent source HHSS and direct quartic are explicit inputs. Without
those inputs, their tests are reported as NOT_EVALUATED, never as passed.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping
import numpy as np

def symmetry_residual(q: np.ndarray) -> float:
    """Relative Bose exchange residual of a real rank-four tensor."""
    q=np.asarray(q,float)
    norm=max(float(np.linalg.norm(q)),1e-100)
    return max(float(np.linalg.norm(q-q.transpose(perm)))/norm for perm in
               ((1,0,2,3),(0,2,1,3),(0,1,3,2)))

MODELS = ('A_-4', 'D_-2')
SCHEMA = 't3-worker-b-quartic-endpoint-audit-v1'
SECTORS = {
    'HHHH': (slice(0,4),)*4,
    'HHSS': (slice(0,4),slice(0,4),slice(4,12),slice(4,12)),
    'SSSS': (slice(4,12),)*4,
}


def sha(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def rel(actual: np.ndarray, reference: np.ndarray) -> float:
    return float(np.linalg.norm(actual-reference)/max(float(np.linalg.norm(reference)),1e-100))


def load_npz_array(path: Path, key: str, shape: tuple[int,...]) -> np.ndarray:
    with np.load(path, allow_pickle=False) as z:
        if key not in z:
            raise ValueError(f'{path}: missing array {key!r}; found {z.files}')
        raw=np.asarray(z[key])
    if raw.shape != shape or not np.isfinite(raw).all():
        raise ValueError(f'{path}: {key} must be finite, shape {shape}')
    if np.iscomplexobj(raw) and np.max(abs(raw.imag)) > 1e-11*max(1.,float(np.max(abs(raw.real)))):
        raise ValueError(f'{path}: {key} is not real')
    return np.asarray(raw.real, float)


def audit_tensors(context: Any, state: Mapping[str,Any], *, hhss:np.ndarray|None=None,
                  direct:np.ndarray|None=None, tolerance:float=2e-9) -> dict:
    """Audit affine endpoint and six-orbit source replacement without mutation."""
    from Numerical.running.T3GeneralScalarBackground import replace_source_hhss
    if context.model not in MODELS: raise ValueError('model unsupported')
    q=np.asarray(context.full_reference(state),float)
    if q.shape != (12,)*4 or not np.isfinite(q).all():
        raise ValueError('cached quartic has invalid shape or nonfinite values')
    dirs=tuple(map(tuple,context.directions))
    if not dirs or len(dirs)!=len(set(dirs)):
        raise ValueError('empty or duplicated cached coupling directions')
    base=np.asarray(context.reference_tensor,float)
    cols=np.asarray(context.columns,float)
    if base.shape!=(12,)*4 or cols.shape!=(len(dirs),)+(12,)*4:
        raise ValueError('invalid cached affine columns')
    residual=symmetry_residual(q)
    column_residual=max(symmetry_residual(x) for x in cols)
    # Partition by exact field content. Every 4-index orbit is accounted for.
    mask=np.indices((12,)*4)<4
    higgs_count=np.sum(mask,axis=0)
    sector_norms={f'H{h}S{4-h}':float(np.linalg.norm(q[higgs_count==h])) for h in range(5)}
    names=sorted({name for name,_ in dirs})
    missing=[name for name in names if name not in state]
    if missing: raise ValueError(f'missing endpoint couplings {missing}')
    result={
        'model':context.model,'response_directions':[list(x) for x in dirs],
        'distinct_couplings':names,
        'candidate_reference_norm':float(np.linalg.norm(base)),
        'candidate_endpoint_norm':float(np.linalg.norm(q)),
        'candidate_sector_norms':sector_norms,
        'candidate_bose_symmetry_relative':residual,
        'maximum_response_column_bose_symmetry_relative':column_residual,
        'source_hhss_evaluated':hhss is not None,
        'direct_quartic_evaluated':direct is not None,
        'physical_quartic_normalization_certified':False,
        'full_operator_inventory_certified':False,
    }
    if hhss is not None:
        h=np.asarray(hhss,float)
        if h.shape!=(4,4,8,8) or not np.isfinite(h).all():
            raise ValueError('independent source HHSS must be finite (4,4,8,8)')
        candidate=q[:4,:4,4:,4:]
        hybrid=replace_source_hhss(q,h)
        changed=(higgs_count==2)
        result.update({
            'source_hhss_vs_candidate_relative':rel(candidate,h),
            'source_hhss_norm':float(np.linalg.norm(h)),
            'hybrid_bose_symmetry_relative':symmetry_residual(hybrid),
            'hybrid_nonhhss_preservation_max_abs':float(np.max(abs(hybrid[~changed]-q[~changed]))),
            'hybrid_hhss_exact_max_abs':float(np.max(abs(hybrid[:4,:4,4:,4:]-h))),
        })
    if direct is not None:
        d=np.asarray(direct,float)
        if d.shape!=(12,)*4 or not np.isfinite(d).all():
            raise ValueError('direct quartic must be finite (12,)*4')
        result['independent_direct_vs_cached_relative']=rel(q,d)
        result['independent_direct_bose_symmetry_relative']=symmetry_residual(d)
        result['independent_direct_sector_relatives']={
            f'H{h}S{4-h}':rel(q[higgs_count==h],d[higgs_count==h]) for h in range(5)
        }
    result['status']='PassedAvailableTensorChecks' if (
        residual<=tolerance and column_residual<=tolerance and
        (hhss is None or (result['hybrid_bose_symmetry_relative']<=tolerance and
                              result['hybrid_nonhhss_preservation_max_abs']==0. and
                              result['hybrid_hhss_exact_max_abs']==0.)) and
        (direct is None or result['independent_direct_vs_cached_relative']<=tolerance)
    ) else 'FailedTensorChecks'
    return result


def parse_pairs(records:list[str], label:str, allowed:set[str]) -> dict[str,Path]:
    out={}
    for item in records:
        name,sep,raw=item.partition('=')
        if not sep or name not in allowed or not raw or name in out:
            raise ValueError(f'{label}: use MODEL=PATH once per model')
        out[name]=Path(raw)
    return out


def find_response_roots(project:Path, explicit:Path|None) -> list[Path]:
    if explicit is not None: return [explicit.resolve()]
    roots=[]
    for base in [project/'output',project/'OUTPUT_B3B']:
        if not base.is_dir():continue
        for report in base.glob('**/quartic_response.npz'):
            if report.parent.name in MODELS and (report.parent/'report.json').is_file():
                roots.append(report.parent.parent.resolve())
    return sorted(set(roots))


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--response-root',type=Path,help='parent directory with A_-4/ and D_-2/ Batch3B reports')
    p.add_argument('--models',default='A_-4,D_-2')
    p.add_argument('--state-file',action='append',default=[],metavar='MODEL=PATH')
    p.add_argument('--source-hhss',action='append',default=[],metavar='MODEL=NPZ')
    p.add_argument('--direct-quartic',action='append',default=[],metavar='MODEL=NPZ')
    p.add_argument('--hhss-key',default='hhss')
    p.add_argument('--direct-key',default='quartic')
    p.add_argument('--find-responses',action='store_true')
    p.add_argument('--out',type=Path,default=Path('output/audits/workerb_quartic_endpoint_v1/report.json'))
    args=p.parse_args(argv)
    project=args.project.resolve()
    models=tuple(x.strip() for x in args.models.split(',') if x.strip())
    if not models or len(set(models))!=len(models) or not set(models)<=set(MODELS):
        p.error('models must be unique subset of A_-4,D_-2')
    roots=find_response_roots(project,args.response_root)
    if args.find_responses:
        if not roots: print('No Batch3B quartic_response.npz/report.json pairs found under output/.')
        for root in roots:
            print('Response root:',root)
            for model in models:
                print(f'  {model}:', 'FOUND' if (root/model/'report.json').is_file() else 'MISSING')
        return 0
    if len(roots)!=1:
        p.error(f'{len(roots)} Batch3B response roots found; use --find-responses and --response-root PATH')
    state_paths=parse_pairs(args.state_file,'state-file',set(models))
    hhss_paths=parse_pairs(args.source_hhss,'source-hhss',set(models))
    direct_paths=parse_pairs(args.direct_quartic,'direct-quartic',set(models))
    from Numerical.running.T3CoupledQuarticContext import CandidateQuarticContext
    # Reuse the audited endpoint discovery and schema validation.
    from scripts.RunT3WorkerBEndpointAudit import _discover, _extract_endpoint
    report={'schema':SCHEMA,'status':'PassedAvailableTensorChecks',
            'batch3b_response_root':str(roots[0]),'results':{},
            'source_hhss_certified':False,'physical_quartic_certified':False}
    for model in models:
        paths=_discover(project,model)
        if model in state_paths:
            raw=state_paths[model]
            path=(raw if raw.is_absolute() else project/raw).resolve()
        elif len(paths)==1: path=paths[0]
        else: p.error(f'{model}: {len(paths)} endpoint inputs; supply --state-file {model}=PATH')
        mu_s, state, state_hash=_extract_endpoint(path,model)
        context=CandidateQuarticContext.load(project,roots[0],model)
        hhss=load_npz_array(hhss_paths[model],args.hhss_key,(4,4,8,8)) if model in hhss_paths else None
        direct=load_npz_array(direct_paths[model],args.direct_key,(12,)*4) if model in direct_paths else None
        data=audit_tensors(context,state,hhss=hhss,direct=direct)
        data.update({'mu_s_gev':mu_s,'endpoint_file':str(path),'endpoint_sha256':state_hash,
                     'source_hhss_file':str(hhss_paths[model]) if model in hhss_paths else None,
                     'independent_direct_file':str(direct_paths[model]) if model in direct_paths else None})
        report['results'][model]=data
        if data['status']!='PassedAvailableTensorChecks':report['status']='FailedTensorChecks'
    report['all_source_hhss_supplied']=all(m in hhss_paths for m in models)
    report['all_independent_direct_supplied']=all(m in direct_paths for m in models)
    report['status']=report['status'] if report['all_source_hhss_supplied'] else report['status']+'_SourceHHSSNotEvaluated'
    target=args.out if args.out.is_absolute() else project/args.out
    if target.exists():p.error(f'output exists, refusing overwrite: {target}')
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2,allow_nan=False))
    return 0 if all(r['status']=='PassedAvailableTensorChecks' for r in report['results'].values()) else 1

if __name__=='__main__':
    raise SystemExit(main())
