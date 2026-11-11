"""Worker B: constrained, NON-CERTIFYING scalar-coordinate basis audit.

Reads the earlier Batch Reconcile NPZ and hashed Batch3B response. Tests
species-preserving complex-component permutations, uniform phases 1/i/-1/-i,
and species-wide conjugation. The direction-space least-squares fit is DIAGNOSTIC,
not a physical coupling redefinition. Neither sources nor production code change.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path
import numpy as np

MODELS=('A_-4','D_-2')
SCHEMA='t3-worker-b-scalar-coordinate-basis-search-v1'


def sha(path):
    h=sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def real_transform(n_complex, perm, phase, conjugate):
    """Canonical [Re z0, Im z0,...] coordinates: old = R @ new.

    old z[k] = phase * (conj(new z[perm[k]]) if conjugate else new z[perm[k]]).
    """
    if sorted(perm)!=list(range(n_complex)) or phase not in (0,1,2,3):
        raise ValueError('Invalid complex permutation or quadrant phase')
    if type(conjugate) is not bool:raise ValueError('conjugation must be bool')
    u=(1j)**phase
    r=np.zeros((2*n_complex,2*n_complex))
    for k,j in enumerate(perm):
        # multiplication by u composed with conjugation
        block=np.array([[u.real,(-u.imag if not conjugate else u.imag)],
                        [u.imag,(u.real if not conjugate else -u.real)]])
        r[2*k:2*k+2,2*j:2*j+2]=block
    if not np.array_equal(r.T@r,np.eye(2*n_complex)):
        raise ValueError('Nonorthogonal scalar transformation')
    return r


def candidate_rotations(model):
    """384 transformations: 6 triplet permutations x 4^2 phases x 2^2 conjugations."""
    if model not in MODELS:raise ValueError('Unsupported scalar model')
    dims=(1,3) if model=='A_-4' else (3,1)
    total=0
    for perm1 in permutations(range(dims[0])):
        for perm2 in permutations(range(dims[1])):
            for phase1,phase2,conj1,conj2 in product(range(4),range(4),(False,True),(False,True)):
                r1=real_transform(dims[0],perm1,phase1,conj1)
                r2=real_transform(dims[1],perm2,phase2,conj2)
                R=np.zeros((8,8));a=2*dims[0]
                R[:a,:a]=r1;R[a:,a:]=r2
                total+=1
                yield {'S1_permutation':list(perm1),'S2_permutation':list(perm2),
                       'S1_phase_i_power':phase1,'S2_phase_i_power':phase2,
                       'S1_conjugate':conj1,'S2_conjugate':conj2},R
    assert total==384


def rotate_hhss(h,R):
    h=np.asarray(h,float);R=np.asarray(R,float)
    if h.shape!=(4,4,8,8) or R.shape!=(8,8) or not np.isfinite(h).all():
        raise ValueError('Invalid HHSS tensor / scalar rotation')
    if not np.allclose(R.T@R,np.eye(8),atol=1e-12):
        raise ValueError('Rotation must be orthogonal')
    return np.einsum('ABij,ia,jb->ABab',h,R,R,optimize=True)


def fit_metrics(candidate,source):
    """Coordinate-level equality, best single scale, and 2-axis source span."""
    c=np.asarray(candidate,float);s=np.asarray(source,float)
    if c.shape!=s.shape or c.ndim!=2 or not np.isfinite(c).all() or not np.isfinite(s).all():
        raise ValueError('Expected finite flattened direction matrices')
    nc=float(np.linalg.norm(c));ns=float(np.linalg.norm(s))
    sc=float(np.sum(c*s)/max(ns**2,1e-100))
    direct=float(np.linalg.norm(c-s)/max(nc,1e-100))
    sign=float(np.linalg.norm(c+s)/max(nc,1e-100))
    scaled=float(np.linalg.norm(c-sc*s)/max(nc,1e-100))
    # Column span: allow a real 2x2 redefinition of (Re lambdaT3, Im lambdaT3)
    X=np.linalg.lstsq(s,c,rcond=1e-11)[0]
    projected=float(np.linalg.norm(c-s@X)/max(nc,1e-100))
    return {'identity_relative_to_candidate':direct,'signflip_relative_to_candidate':sign,
            'best_global_scale':sc,'after_global_scale_relative_to_candidate':scaled,
            'after_2x2_coupling_mix_relative_to_candidate':projected,
            'coupling_mix_source_to_candidate':X.tolist(),
            'candidate_norm':nc,'source_norm':ns}


def search(model,candidate_columns,source_columns,directions,*,top=8):
    """Rank candidates using the lambdaT3 2D span; do not fit physics values."""
    if model not in MODELS:raise ValueError('Invalid model')
    d=[tuple(x) for x in directions]
    if len(set(d))!=len(d):raise ValueError('Duplicate response directions')
    c=np.asarray(candidate_columns,float);s=np.asarray(source_columns,float)
    if c.shape!=s.shape or c.shape!=(len(d),4,4,8,8):
        raise ValueError('Expected matching HHSS response arrays')
    idx=[]
    for axis in ('real','imag'):
        if ('lambdaT3',axis) not in d:raise ValueError('Missing lambdaT3.'+axis)
        idx.append(d.index(('lambdaT3',axis)))
    C=c[idx].reshape(2,-1).T
    results=[]
    for description,R in candidate_rotations(model):
        S=np.stack([rotate_hhss(s[i],R) for i in idx]).reshape(2,-1).T
        m=fit_metrics(C,S)
        # Report basis residual independent of simple lambdaT3 phase redefinitions.
        results.append((m['after_2x2_coupling_mix_relative_to_candidate'],
                        m['after_global_scale_relative_to_candidate'],
                        m['identity_relative_to_candidate'],description,m,R))
    results.sort(key=lambda x:(x[0],x[1],x[2]))
    base=fit_metrics(C,s[idx].reshape(2,-1).T)
    records=[]
    for _,_,_,desc,m,_ in results[:top]:records.append({'rotation':desc,**m})
    best_R=results[0][-1]
    # A best-fit 2x2 coupling map can hide a large mismatch in other directions.
    S_all=np.stack([rotate_hhss(q,best_R) for q in s]).reshape(len(d),-1).T
    C_all=c.reshape(len(d),-1).T
    X=np.linalg.lstsq(S_all,C_all,rcond=1e-11)[0]
    all_resid=float(np.linalg.norm(C_all-S_all@X)/max(float(np.linalg.norm(C_all)),1e-100))
    # Additional exact fixed-label check for lambdaH2, where a pure minus sign was observed.
    portal=[]
    if ('lambdaH2','real') in d:
        j=d.index(('lambdaH2','real'))
        p=rotate_hhss(s[j],best_R)
        pc=c[j]; denom=max(float(np.linalg.norm(pc)),1e-100)
        portal=[{'direction':['lambdaH2','real'],
                 'best_rotation_identity_residual':float(np.linalg.norm(pc-p)/denom),
                 'best_rotation_signflip_residual':float(np.linalg.norm(pc+p)/denom)}]
    return {'model':model,'searched_rotations':len(results),
            'baseline_lambdaT3':base,'best_lambdaT3':records[0],
            'top_rotations':records,'all_direction_span_residual_at_best_rotation':all_resid,
            'lambdaH2_check':portal,
            'interpretation':'Finite constrained basis hypothesis tests only. No scalar-basis/coupling identification certified.',
            'physical_normalization_certified':False}


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--models',default='A_-4,D_-2')
    p.add_argument('--response-root',type=Path,required=True)
    p.add_argument('--reconcile-dir',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--top',type=int,default=8)
    args=p.parse_args(argv)
    models=tuple(x.strip() for x in args.models.split(',') if x.strip())
    if not models or len(set(models))!=len(models) or not set(models)<=set(MODELS):
        p.error('Select A_-4,D_-2')
    if args.top<1 or args.top>40:p.error('top must be between 1 and 40')
    project=args.project.resolve()
    def full(q):return (q if q.is_absolute() else project/q).resolve()
    out=full(args.out)
    if out.exists():p.error('Refusing to overwrite existing report: '+str(out))
    response_root=full(args.response_root);reconcile_dir=full(args.reconcile_dir)
    from Numerical.running.T3CoupledQuarticContext import CandidateQuarticContext
    report_path=reconcile_dir/'report.json'
    report=json.loads(report_path.read_text(encoding='utf-8-sig'))
    if report.get('schema')!='t3-worker-b-batch-hhss-reconciliation-v1':
        raise ValueError('Wrong source reconciliation report schema')
    records={}
    for model in models:
        ctx=CandidateQuarticContext.load(project,response_root,model)
        rec=report['results'][model]
        npz=reconcile_dir/f'{model}_reconciliation.npz'
        if rec.get('npz_sha256')!=sha(npz):raise ValueError(f'{model}: altered reconciliation tensor')
        with np.load(npz,allow_pickle=False) as f:
            src=f['source_direction_columns']
        records[model]=search(model,ctx.columns[:,:4,:4,4:,4:],src,ctx.directions,top=args.top)
        b=records[model]['best_lambdaT3']
        print(f'{model}: baseline lambdaT3 2x2 span={records[model]["baseline_lambdaT3"]["after_2x2_coupling_mix_relative_to_candidate"]:.6g}; '
              f'best constrained basis span={b["after_2x2_coupling_mix_relative_to_candidate"]:.6g}; '
              f'all-directions span={records[model]["all_direction_span_residual_at_best_rotation"]:.6g}')
    out.parent.mkdir(parents=True,exist_ok=True)
    output={'schema':SCHEMA,'status':'ConstrainedScalarBasisHypothesesEvaluated',
            'response_root':str(response_root),'reconciliation_root':str(reconcile_dir),
            'results':records,'applied_to_production':False,
            'physical_basis_identification_certified':False}
    out.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Report:',out)
    return 0


if __name__=='__main__':raise SystemExit(main())
