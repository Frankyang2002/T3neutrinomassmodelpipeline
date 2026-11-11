"""Batch Worker B: HHSS direction geometry, fixed-basis hybrid export, and audits.

A_-4/D_-2 only.  Evaluate the existing source-aware HHSS independently,
compare it with the hashed Batch3B affine candidate, and replace *only* the
six HHSS placements in the fixed-basis candidate quartic.  This is a
conditional EFT1 scalar background, not a physical normalization certificate.

No renormalisation-group solve, no fermion matching, no physical C5, no
modification of historical files.  All output goes to a fresh directory.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import numpy as np

from scripts.RunT3WorkerBHHSSDirectionAudit import (
    SourceEvaluator, audit_responses, number, put, check_symmetry,
)

MODELS = ('A_-4', 'D_-2')
SCHEMA = 't3-worker-b-batch-hhss-reconciliation-v1'


def rel(a: np.ndarray, b: np.ndarray) -> float | None:
    a = np.asarray(a, float); b = np.asarray(b, float)
    d = float(np.linalg.norm(b))
    return float(np.linalg.norm(a-b) / d) if d else None


def _sha(path: Path) -> str:
    h = sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def _sector(a: int, b: int, split: int) -> str:
    return ('S1S1' if a < split and b < split else
            'S2S2' if a >= split and b >= split else 'S1S2')


def component_diagnostic(candidate: np.ndarray, source: np.ndarray, *, model: str, top=16) -> dict:
    """Report highest absolute component residuals in fixed real coordinates."""
    if model not in MODELS:raise ValueError('invalid model for component labels')
    split=2 if model=='A_-4' else 6
    c = np.asarray(candidate, float); s = np.asarray(source, float)
    check_symmetry(c, 'component candidate')
    check_symmetry(s, 'component source')
    diff = c-s
    n = min(max(int(top), 1), diff.size)
    indices = np.argpartition(np.abs(diff).ravel(), -n)[-n:]
    indices = sorted(indices, key=lambda j: (-abs(diff.ravel()[j]), int(j)))
    rows=[]
    for j in indices:
        A, B, a, b = map(int, np.unravel_index(int(j), diff.shape))
        rows.append({'indices_zero_based_HHSS':[A,B,a,b],
                     'candidate':float(c[A,B,a,b]),
                     'source':float(s[A,B,a,b]),
                     'difference':float(diff[A,B,a,b]),
                     'heavy_species_pair':_sector(a,b,split)})
    return {'maximum_absolute_difference':float(np.max(np.abs(diff))),
            'largest_components':rows}


def species_blocks(candidate: np.ndarray, source: np.ndarray, *, model: str) -> dict:
    """A/D: one singlet (2 real), one triplet (6 real), ordered by S1/S2."""
    if model not in MODELS:raise ValueError('invalid model for species blocks')
    split=2 if model=='A_-4' else 6
    c=np.asarray(candidate,float);s=np.asarray(source,float)
    check_symmetry(c,'block candidate');check_symmetry(s,'block source')
    out={}
    for key, aa, bb in [('S1S1',slice(0,split),slice(0,split)),
                        ('S1S2',slice(0,split),slice(split,8)),
                        ('S2S1',slice(split,8),slice(0,split)),
                        ('S2S2',slice(split,8),slice(split,8))]:
        x=c[:,:,aa,bb];y=s[:,:,aa,bb]
        ns=float(np.linalg.norm(y));nc=float(np.linalg.norm(x))
        out[key]={'source_norm':ns,'candidate_norm':nc,
                  'relative_difference':rel(x,y),
                  'best_scale_candidate_over_source':float(np.vdot(y,x).real/ns**2) if ns else None}
    return out


def direction_geometry(directions: list[tuple[str,str]], candidate_columns: np.ndarray,
                       source_columns: np.ndarray, *, rtol=1e-9) -> dict:
    """Cross-direction candidate/source subspace comparison by SVD.

    A direction-wise scale can conceal mixing between distinct invariant
    basis vectors. Compare whole directional spans; do NOT fit or apply a
    correction to the production tensor.
    """
    d=tuple(tuple(v) for v in directions)
    c=np.asarray(candidate_columns,float);s=np.asarray(source_columns,float)
    if c.shape!=(len(d),4,4,8,8) or s.shape!=c.shape or not np.isfinite(c).all() or not np.isfinite(s).all():
        raise ValueError('direction column geometry has invalid shape')
    if len(d)!=len(set(d)):raise ValueError('duplicate quartic directions')
    for x in [*c,*s]:check_symmetry(x,'direction geometry')
    C=c.reshape(len(d),-1).T
    S=s.reshape(len(d),-1).T
    sv=np.linalg.svd(S,compute_uv=False)
    tol=(sv[0]*rtol if sv.size else 0.)
    rank=int(np.count_nonzero(sv>tol))
    X=np.linalg.lstsq(S,C,rcond=rtol)[0] if rank else np.zeros((len(d),len(d)))
    residual=np.linalg.norm(C-S@X)/max(float(np.linalg.norm(C)),1e-100)
    norms=np.linalg.norm(S,axis=0)
    gram=S.T@C
    cross=[]
    for j, (name,axis) in enumerate(d):
        strengths=[]
        for k, label in enumerate(d):
            denom=float(np.linalg.norm(C[:,j])*norms[k])
            cosine=float(gram[k,j]/denom) if denom else None
            if cosine is not None:strengths.append((abs(cosine),k,cosine))
        strengths.sort(reverse=True)
        cross.append({'candidate_direction':[name,axis],
                      'top_source_overlaps':[
                          {'source_direction':list(d[k]),'cosine':float(v)}
                          for _,k,v in strengths[:5]],
                      'candidate_column_outside_source_span_relative':
                          float(np.linalg.norm(C[:,j]-S@X[:,j])/
                                max(np.linalg.norm(C[:,j]),1e-100))})
    return {'source_column_rank':rank,
            'source_column_singular_values':sv.tolist(),
            'candidate_direction_subspace_relative_residual':float(residual),
            'best_linear_direction_map_source_to_candidate':X.tolist(),
            'direction_order':[list(x) for x in d],
            'candidate_top_cross_overlaps':cross,
            'note':'Diagnostic least-squares map only. No fitted rescaling or mixing is applied.'}


def hybrid_tensor(candidate: np.ndarray, source_hhss: np.ndarray) -> tuple[np.ndarray,dict]:
    """Replace source HHSS once, retaining all original non-HHSS entries."""
    from Numerical.running.T3CoupledQuarticContext import _put_hhss
    from Numerical.running.T3ScalarQuarticLinearResponse import symmetry_residual
    q=np.asarray(candidate,float)
    h=np.asarray(source_hhss,float)
    if q.shape!=(12,)*4 or not np.isfinite(q).all():
        raise ValueError('full quartic must be finite rank-4 12-real tensor')
    check_symmetry(h,'source to hybrid')
    result=_put_hhss(q.copy(),h)
    sym=symmetry_residual(result)
    if sym>1e-10:raise ValueError(f'hybrid quartic Bose symmetry failure {sym}')
    mask=np.zeros(q.shape,bool)
    for ax1 in range(4):
        for ax2 in range(ax1+1,4):
            bounds=[slice(4,None)]*4
            bounds[ax1]=slice(None,4);bounds[ax2]=slice(None,4)
            mask[tuple(bounds)]=True
    if not np.array_equal(result[~mask],q[~mask]):
        raise ValueError('hybrid modification outside HHSS positions')
    if not np.allclose(result[:4,:4,4:,4:],h,atol=1e-12,rtol=1e-12):
        raise ValueError('source HHSS not reproduced in hybrid')
    return result,{'full_bose_symmetry_relative':float(sym),
                   'non_hhss_unchanged_exactly':True,
                   'source_hhss_replaced_exactly':True,
                   'candidate_hhss_relative_to_source':rel(q[:4,:4,4:,4:],h),
                   'hybrid_rank4_shape':list(result.shape)}


def source_direction_columns(state, direct_source, direction_order):
    """Independent finite difference of source evaluator at every cached axis."""
    base=np.asarray(direct_source(state),float)
    check_symmetry(base,'source base')
    columns=[]
    for name,axis in direction_order:
        if name not in state or axis not in ('real','imag'):
            raise ValueError('invalid source response direction')
        z=number(state[name]); is_complex=(name,'imag') in direction_order
        probe=put(state,name,z+(1j if axis=='imag' else 1.),is_complex)
        col=np.asarray(direct_source(probe),float)-base
        check_symmetry(col,'source '+name+'.'+axis)
        columns.append(col)
    return base,np.stack(columns)


def reconcile(model, ctx, state, source_evaluator, *, top=16):
    """All computations share one endpoint, candidate fixed scalar basis."""
    from Numerical.running.T3ScalarQuarticLinearResponse import symmetry_residual
    if model not in MODELS or model!=ctx.model:
        raise ValueError('invalid model/context pairing')
    full=ctx.full_reference(state)
    candidate_h=full[:4,:4,4:,4:].copy()
    cc=np.asarray(ctx.columns[:,:4,:4,4:,4:],float)
    source, sc=source_direction_columns(state,source_evaluator,ctx.directions)
    diag=audit_responses(state,candidate_h,ctx.directions,cc,source_evaluator)
    hybrid,verification=hybrid_tensor(full,source)
    geometry=direction_geometry(list(ctx.directions),cc,sc)
    diag.update({'component_diagnostic':component_diagnostic(candidate_h,source,model=model,top=top),
                 'heavy_species_blocks':species_blocks(candidate_h,source,model=model),
                 'cross_direction_geometry':geometry,
                 'hybrid_verification':verification,
                 'candidate_bose_symmetry_relative':float(symmetry_residual(full)),
                 'physical_normalization_certified':False,
                 'full_source_invariant_inventory_certified':False,
                 'non_hhss_candidate_basis_certified':False})
    return {'source_hhss':source,'source_direction_columns':sc,
            'candidate_hhss':candidate_h,'hybrid_full_tensor':hybrid},diag


def _pairs(items,models):
    out={}
    for item in items:
        model,sep,path=item.partition('=')
        if not sep or model not in models or not path or model in out:
            raise ValueError('provide one --state-file MODEL=PATH per selected model')
        out[model]=Path(path)
    return out


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--models',default='A_-4,D_-2')
    p.add_argument('--response-root',type=Path,required=True)
    p.add_argument('--state-file',action='append',default=[],metavar='MODEL=PATH')
    p.add_argument('--release-dir',type=Path,default=Path('output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'))
    p.add_argument('--v4-dir',type=Path,default=Path('output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'))
    p.add_argument('--seed-dir',type=Path,default=Path('output/full16_symbolic'))
    p.add_argument('--registry',type=Path,default=Path('RGE/general/generated/general_rgbeta_quartic_registry.json'))
    p.add_argument('--out-dir',type=Path,default=Path('output/audits/workerb_batch_reconcile_v1'))
    p.add_argument('--top',type=int,default=16)
    args=p.parse_args(argv)
    models=tuple(x.strip() for x in args.models.split(',') if x.strip())
    if not models or len(set(models))!=len(models) or not set(models)<=set(MODELS):
        p.error('unique subset A_-4,D_-2 required')
    if args.top<1 or args.top>100:p.error('top must be in 1..100')
    project=args.project.resolve()
    def absolute(x):
        q=Path(x)
        return (q if q.is_absolute() else project/q).resolve()
    paths=_pairs(args.state_file,models)
    if set(paths)!=set(models):p.error('explicit --state-file needed for each model')
    out=absolute(args.out_dir)
    if out.exists():p.error(f'refusing to overwrite existing results at {out}')
    response_root=absolute(args.response_root)
    srcpaths=[absolute(args.release_dir),absolute(args.v4_dir),absolute(args.seed_dir),absolute(args.registry)]
    inputs=[response_root,*srcpaths,*[absolute(paths[m]) for m in models]]
    for x in inputs:
        if not x.exists():p.error(f'missing input: {x}')
    from Numerical.running.T3CoupledQuarticContext import CandidateQuarticContext
    from scripts.RunT3WorkerBEndpointAudit import _extract_endpoint
    records={};arrays={}
    for model in models:
        path=absolute(paths[model]);mu,state,state_sha=_extract_endpoint(path,model)
        ctx=CandidateQuarticContext.load(project,response_root,model)
        evaluator=SourceEvaluator(model,*srcpaths,state)
        data,result=reconcile(model,ctx,state,evaluator,top=args.top)
        result.update({'model':model,'mu_s_gev':mu,'endpoint_path':str(path),
                       'endpoint_sha256':state_sha,
                       'source_evaluator_calls':evaluator.evaluations,
                       'source_release_path':str(evaluator.boundary.path),
                       'source_release_sha256':_sha(evaluator.boundary.path),
                       'native_quartic_seed_path':str(evaluator.channels.seed_path),
                       'native_quartic_seed_sha256':_sha(evaluator.channels.seed_path)})
        arrays[model]=data;records[model]=result
        print(f'{model}: source affine={result["source_affine_interpolation_relative"]:.3e}; '
              f'HHSS candidate/source={result["endpoint_comparison"]["candidate_vs_source_relative"]:.6g}; '
              f'candidate outside source direction span={result["cross_direction_geometry"]["candidate_direction_subspace_relative_residual"]:.6g}')
        for record in result['cross_direction_geometry']['candidate_top_cross_overlaps']:
            if record['candidate_direction'][0]=='lambdaT3':
                print('  lambdaT3.'+record['candidate_direction'][1]+': '+
                      str(record['top_source_overlaps'][:3]))
    out.mkdir(parents=True,exist_ok=False)
    for model,data in arrays.items():
        npz=out/f'{model}_reconciliation.npz'
        np.savez_compressed(npz,**data)
        records[model]['npz_path']=str(npz)
        records[model]['npz_sha256']=_sha(npz)
    report={'schema':SCHEMA,'status':'ComputedSourceHHSSHybrid_CandidateNormalizationOpen',
            'response_root':str(response_root),
            'results':records,
            'muF_integrate_out_all_fermions_gev':1e5,
            'muS_joint_integrate_out_both_scalars_gev':1e3,
            'hybrid_used_for_physical_C5':False,
            'production_solver_changed':False,
            'absolute_physical_normalization_certified':False}
    report_path=out/'report.json'
    report_path.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Report:',report_path)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
