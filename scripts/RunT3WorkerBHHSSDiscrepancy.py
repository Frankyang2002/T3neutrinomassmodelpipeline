"""Read-only diagnostic of A_-4/D_-2 Batch3B vs independently saved source HHSS.

Does NOT rescale, sign-flip, alter, or certify the physical quartic.
Inputs are the same saved running endpoint and source_hhss_native.npz used by
RunT3WorkerBSourceHHSSAudit. Diagnostics are descriptive, not a correction.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any
import numpy as np

MODELS = ('A_-4', 'D_-2')
SCHEMA = 't3-worker-b-hhss-discrepancy-diagnostics-v1'
DIMS = {'A_-4': (1, 3), 'D_-2': (3, 1)}


def _finite_array(value: Any, shape: tuple[int, ...], label: str) -> np.ndarray:
    a = np.asarray(value)
    if a.shape != shape or not np.isfinite(a).all() or np.iscomplexobj(a):
        raise ValueError(f'{label}: expected finite real array with shape {shape}')
    return np.array(a, dtype=float, copy=True)


def _relative(a: np.ndarray, b: np.ndarray) -> float | None:
    norm = float(np.linalg.norm(b))
    return float(np.linalg.norm(a - b) / norm) if norm > 0 else None


def _metrics(candidate: np.ndarray, source: np.ndarray) -> dict:
    """Projection c ~= alpha*s and shape residual; no parameter is modified."""
    if candidate.shape != source.shape:
        raise ValueError('candidate/source HHSS shapes disagree')
    c, s = candidate.reshape(-1), source.reshape(-1)
    cc = float(np.dot(c, c)); ss = float(np.dot(s, s)); cs = float(np.dot(c, s))
    nc, ns = np.sqrt(cc), np.sqrt(ss)
    if ns == 0:
        return {'candidate_norm':nc, 'source_norm':ns, 'relative_to_source':None,
                'cosine':None, 'least_squares_candidate_over_source':None,
                'residual_after_best_real_scale_relative_to_source':None,
                'residual_if_candidate_sign_flipped_relative_to_source':None,
                'status':'ZeroSourceNorm'}
    alpha = cs / ss
    residual_sq = max(cc - 2 * alpha * cs + alpha * alpha * ss, 0.)
    residual = np.sqrt(residual_sq) / ns
    return {'candidate_norm':nc, 'source_norm':ns,
            'relative_to_source':_relative(candidate, source),
            'cosine':cs/(nc*ns) if nc else None,
            'least_squares_candidate_over_source':alpha,
            'residual_after_best_real_scale_relative_to_source':float(residual),
            'residual_if_candidate_sign_flipped_relative_to_source':_relative(-candidate, source),
            'status':'Compared'}


def diagnose(model: str, candidate: np.ndarray, source: np.ndarray, *, top_n: int = 12) -> dict:
    """Break down mismatch by Higgs components and S1/S2 field content.

    S1 and S2 are canonical real interleaved multiplet components in the fixed
    native 8-state basis. Blocks cover the full HHSS index set exactly once.
    """
    if model not in MODELS: raise ValueError('unsupported model')
    c = _finite_array(candidate, (4, 4, 8, 8), 'candidate')
    s = _finite_array(source, (4, 4, 8, 8), 'source')
    if not isinstance(top_n, int) or not 0 <= top_n <= 100:
        raise ValueError('top_n must be between 0 and 100')
    d1, _ = DIMS[model]
    heavy = {'S1':slice(0, 2*d1), 'S2':slice(2*d1, 8)}
    blocks = {}
    for p in ('S1','S2'):
        for q in ('S1','S2'):
            blocks[f'{p}{q}'] = _metrics(c[:,:,heavy[p],heavy[q]],s[:,:,heavy[p],heavy[q]])
    higgs = {f'H{a}H{b}':_metrics(c[a,b],s[a,b])
             for a in range(4) for b in range(a,4)}
    difference = np.abs(c-s)
    flat = difference.ravel()
    take = sorted(range(flat.size), key=lambda k:(-float(flat[k]),k))[:top_n]
    largest = []
    for k in take:
        a,b,i,j = map(int,np.unravel_index(k,difference.shape))
        largest.append({'indices_zero_based':[a,b,i,j],
                        'candidate':float(c[a,b,i,j]),'source':float(s[a,b,i,j]),
                        'absolute_difference':float(difference[a,b,i,j]),
                        'heavy_species':[('S1' if i<2*d1 else 'S2'),
                                         ('S1' if j<2*d1 else 'S2')]})
    overall = _metrics(c,s)
    return {'model':model,'overall':overall,
            'heavy_species_blocks':blocks,'higgs_component_blocks':higgs,
            'largest_absolute_entry_differences':largest,
            'source_symmetry':{
                'higgs_swap_absolute_max':float(np.max(abs(s-s.swapaxes(0,1)))),
                'heavy_swap_absolute_max':float(np.max(abs(s-s.swapaxes(2,3))))},
            'candidate_symmetry':{
                'higgs_swap_absolute_max':float(np.max(abs(c-c.swapaxes(0,1)))),
                'heavy_swap_absolute_max':float(np.max(abs(c-c.swapaxes(2,3))))},
            'interpretation':'Projection and sign diagnostics only; no automatic tensor correction',
            'physical_normalization_certified':False}


def _pair_args(items: list[str], models: tuple[str, ...], opt: str) -> dict[str, Path]:
    out = {}
    for item in items:
        key, eq, value = item.partition('=')
        if not eq or key not in models or not value or key in out:
            raise ValueError(f'{opt}: supply MODEL=PATH once per selected model')
        out[key] = Path(value)
    return out


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--models',default='A_-4,D_-2')
    p.add_argument('--response-root',type=Path,required=True)
    p.add_argument('--state-file',action='append',default=[],metavar='MODEL=PATH')
    p.add_argument('--source-dir',type=Path,default=Path('output/audits/workerb_source_hhss_v1'))
    p.add_argument('--source-hhss',action='append',default=[],metavar='MODEL=NPZ')
    p.add_argument('--out',type=Path,default=Path('output/audits/workerb_hhss_diagnostic_v1/report.json'))
    p.add_argument('--top',type=int,default=12)
    args=p.parse_args(argv)
    models=tuple(x.strip() for x in args.models.split(',') if x.strip())
    if not models or len(set(models))!=len(models) or not set(models)<=set(MODELS):
        p.error('only unique A_-4,D_-2 models are supported')
    project=args.project.resolve()
    def absolute(x):
        x=Path(x)
        return (x if x.is_absolute() else project/x).resolve()
    state_paths=_pair_args(args.state_file,models,'--state-file')
    source_paths=_pair_args(args.source_hhss,models,'--source-hhss')
    if set(state_paths)!=set(models):p.error('explicit --state-file required for each model')
    response_root=absolute(args.response_root)
    source_dir=absolute(args.source_dir)
    out=absolute(args.out)
    if out.exists():p.error(f'refusing to overwrite existing result: {out}')
    from Numerical.running.T3CoupledQuarticContext import CandidateQuarticContext
    from scripts.RunT3WorkerBEndpointAudit import _extract_endpoint
    from scripts.RunT3WorkerBQuarticAudit import load_npz_array
    report={'schema':SCHEMA,'status':'Computed_NotCertified',
            'response_root':str(response_root),'source_dir':str(source_dir),
            'results':{},'physical_normalization_certified':False,
            'source_hhss_tensor_modified':False}
    for model in models:
        state_path=absolute(state_paths[model])
        mu,state,state_hash=_extract_endpoint(state_path,model)
        ctx=CandidateQuarticContext.load(project,response_root,model)
        candidate=ctx.full_reference(state)[:4,:4,4:,4:]
        source_path=absolute(source_paths[model]) if model in source_paths else source_dir/f'{model}_source_hhss_native.npz'
        source=load_npz_array(source_path,'hhss',(4,4,8,8))
        diag=diagnose(model,candidate,source,top_n=args.top)
        diag.update({'mu_s_gev':mu,'endpoint_file':str(state_path),
                     'endpoint_sha256':state_hash,'source_hhss_file':str(source_path),
                     'source_hhss_sha256':_hash(source_path)})
        report['results'][model]=diag
        overall=diag['overall']
        print(f'{model}: relative={overall["relative_to_source"]:.9g}, '
              f'best_scale={overall["least_squares_candidate_over_source"]:.9g}, '
              f'best_scale_residual={overall["residual_after_best_real_scale_relative_to_source"]:.9g}, '
              f'sign_flip_residual={overall["residual_if_candidate_sign_flipped_relative_to_source"]:.9g}')
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Report:',out)
    return 0


if __name__=='__main__':raise SystemExit(main())
