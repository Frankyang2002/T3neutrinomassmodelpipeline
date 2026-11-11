"""Read-only orthogonal-basis-invariant audit for T3 lambdaT3 HHSS tensors.

Input is the hash-bound Worker B reconciliation NPZ and a hash-validated Batch3B
response context.  For HHSS q_ABab, singular values of the 16x64 flattening
are invariant under ANY orthogonal O(4) Higgs and O(8) heavy scalar basis
change: M -> (O_H kron O_H)^T M (O_S kron O_S).  An incompatible
singular-value spectrum rules out orthogonal field rotations plus a scalar
normalization for the *fixed* coupling direction.  An angular scan over the
source lambdaT3 real/imag plane is DIAGNOSTIC, not a proof about arbitrary
GL(2,R) coupling redefinitions.

No coupling correction, new matching run, production patch, or physical C5.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import numpy as np

MODELS=('A_-4','D_-2')
SCHEMA='t3-worker-b-hhss-spectral-obstruction-v1'


def sha(path):
    h=sha256()
    with Path(path).open('rb') as f:
        for buf in iter(lambda:f.read(1048576),b''):
            h.update(buf)
    return h.hexdigest()


def check_hhss(h):
    q=np.asarray(h,float)
    if q.shape!=(4,4,8,8) or not np.isfinite(q).all():
        raise ValueError('HHSS must be finite 4x4x8x8')
    norm=max(float(np.linalg.norm(q)),1e-100)
    e=max(float(np.linalg.norm(q-q.swapaxes(0,1)))/norm,
          float(np.linalg.norm(q-q.swapaxes(2,3)))/norm)
    if e>1e-9:raise ValueError(f'HHSS Bose symmetry violated: {e}')
    return q


def singular_fingerprint(h):
    q=check_hhss(h)
    singular=np.linalg.svd(q.reshape(16,64),compute_uv=False)
    n=float(np.linalg.norm(singular))
    if n<=1e-14:
        return {'norm':n,'normalized_singular_values':None,'rank':0,
                'stable_rank':None}
    sn=singular/n
    return {'norm':n,'normalized_singular_values':sn.tolist(),
            'rank':int(np.count_nonzero(singular>1e-10*singular[0])),
            'stable_rank':float(np.sum(singular**2)/(singular[0]**2))}


def compare_fingerprints(candidate,source):
    """Spectral obstruction for the given fixed coupling direction.

    The L2 distance of UNIT-normalized singular values is a LOWER BOUND
    on best relative error after any orthogonal left/right multiplication
    and an arbitrary nonzero real scale, when both tensors have unit norm.
    It is a lower bound under the larger O(16)xO(64) class and therefore
    also under physically admissible O(4)xO(8). Never assert sufficiency.
    """
    a=singular_fingerprint(candidate); b=singular_fingerprint(source)
    if a['normalized_singular_values'] is None or b['normalized_singular_values'] is None:
        return {'status':'ZeroTensor','candidate':a,'source':b,
                'normalized_spectral_lower_bound':None}
    u=np.asarray(a['normalized_singular_values']);v=np.asarray(b['normalized_singular_values'])
    return {'status':'Compared','candidate':a,'source':b,
            'normalized_spectral_lower_bound':float(np.linalg.norm(u-v)),
            'rank_matches':a['rank']==b['rank'],
            'exact_orthogonal_rotation_plus_scale_excluded':bool(np.linalg.norm(u-v)>1e-8)}


def rotate_hhss(h,heavy,higgs=None):
    q=check_hhss(h);S=np.asarray(heavy,float)
    H=np.eye(4) if higgs is None else np.asarray(higgs,float)
    for r,n in ((S,8),(H,4)):
        if r.shape!=(n,n) or not np.allclose(r.T@r,np.eye(n),atol=1e-10):
            raise ValueError('Nonorthogonal field transformation')
    return np.einsum('ABij,AC,BD,ik,jl->CDkl',q,H,H,S,S,optimize=True)


def species_sector_fingerprints(h,model):
    """Pure S1/S1 and S2/S2 block singular values, invariant under
    species-preserving orthogonal rotations. Cross block likewise for
    species-preserving transformations, but not under S1-S2 mixing.
    """
    if model not in MODELS:raise ValueError('Unknown model')
    q=check_hhss(h); n=2 if model=='A_-4' else 6
    blocks={'S1S1':q[:,:,:n,:n], 'S1S2':q[:,:,:n,n:],
            'S2S1':q[:,:,n:,:n], 'S2S2':q[:,:,n:,n:]}
    out={}
    for k,t in blocks.items():
        singular=np.linalg.svd(t.reshape(16,-1),compute_uv=False)
        norm=float(np.linalg.norm(singular))
        out[k]={'norm':norm,'normalized_singular_values':(singular/norm).tolist() if norm>1e-14 else None,
                'rank':int(np.count_nonzero(singular>singular[0]*1e-10)) if norm>1e-14 else 0}
    return out


def scan_source_phase(candidate,source_real,source_imag,*,samples=256):
    """A directional scan over normalized linear combinations of the
    source lambdaT3 real and imaginary tensor responses.  Every sampled
    candidate-source comparison allows independent scalar normalization.
    This is NOT an exhaustive GL2 coupling-basis search.
    """
    if samples<8 or samples>4096:raise ValueError('samples must be 8..4096')
    c=check_hhss(candidate);r=check_hhss(source_real);i=check_hhss(source_imag)
    fingerprint=singular_fingerprint(c)
    if fingerprint['normalized_singular_values'] is None:
        raise ValueError('Candidate lambdaT3 response is zero')
    cs=np.asarray(fingerprint['normalized_singular_values'])
    best=None
    for j in range(samples):
        theta=np.pi*j/samples # opposite signs have same singular spectrum
        mixed=np.cos(theta)*r+np.sin(theta)*i
        src=singular_fingerprint(mixed)
        if src['normalized_singular_values'] is None:continue
        d=float(np.linalg.norm(cs-np.asarray(src['normalized_singular_values'])))
        if best is None or d<best['spectral_distance']:
            best={'theta_rad':float(theta),'theta_deg':float(np.degrees(theta)),
                  'spectral_distance':d,'source_mix_norm':src['norm'],
                  'source_mix_rank':src['rank']}
    if best is None:raise ValueError('Both source lambdaT3 columns vanish')
    return {'samples':samples,'best_sampled_phase':best,
            'interpretation':'Finite angular scan over source lambdaT3 plane, allowing overall amplitude; not a proven global lower bound for arbitrary GL2 coupling maps.'}


def diagnose(model, directions, candidate, source,*,samples=256):
    if model not in MODELS:raise ValueError('Unknown model')
    names=tuple(tuple(x) for x in directions)
    c=np.asarray(candidate,float);s=np.asarray(source,float)
    if len(set(names))!=len(names) or c.shape!=s.shape or c.shape!=(len(names),4,4,8,8):
        raise ValueError('Direction inventory or response shape mismatch')
    indexes=[]
    for name in ('lambdaT3','lambdaH2'):
        axes=['real','imag'] if name=='lambdaT3' else ['real']
        for axis in axes:
            if (name,axis) not in names:raise ValueError(f'Missing direction {name}.{axis}')
            indexes.append((name,axis,names.index((name,axis))))
    rows={}
    for name,axis,j in indexes:
        v=compare_fingerprints(c[j],s[j]);
        v['source_species_sectors']=species_sector_fingerprints(s[j],model)
        v['candidate_species_sectors']=species_sector_fingerprints(c[j],model)
        # Identity, opposite sign and optimal scalar-fit residuals are basis-dependent diagnostics.
        nc=max(float(np.linalg.norm(c[j])),1e-100);ns2=float(np.sum(s[j]**2))
        alpha=float(np.sum(s[j]*c[j])/ns2) if ns2>1e-24 else None
        v.update({'identity_rel_to_candidate':float(np.linalg.norm(c[j]-s[j])/nc),
                  'signflip_rel_to_candidate':float(np.linalg.norm(c[j]+s[j])/nc),
                  'best_fixed_basis_scale':alpha,
                  'best_fixed_basis_scale_residual':float(np.linalg.norm(c[j]-alpha*s[j])/nc) if alpha is not None else None})
        rows[f'{name}.{axis}']=v
    ir=names.index(('lambdaT3','real'));ii=names.index(('lambdaT3','imag'))
    phase={axis:scan_source_phase(c[j],s[ir],s[ii],samples=samples)
           for axis,j in [('real',ir),('imag',ii)]}
    return {'model':model,'direction_count':len(names), 'directions':rows,
            'lambdaT3_source_phase_scan':phase,
            'orthogonal_invariance_note':'Fixed-direction singular spectra invariant under independent real orthogonal Higgs and heavy-field rotations. Nonzero singular-spectrum distance forbids exact rotation+scale for that fixed direction only.',
            'physical_normalization_certified':False,'coupling_basis_identification_certified':False}


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--models',default='A_-4,D_-2')
    p.add_argument('--response-root',type=Path,required=True)
    p.add_argument('--reconcile-dir',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--phase-samples',type=int,default=256)
    args=p.parse_args(argv)
    models=tuple(x.strip() for x in args.models.split(',') if x.strip())
    if not models or len(set(models))!=len(models) or not set(models)<=set(MODELS):p.error('Select A_-4,D_-2')
    if not 8<=args.phase_samples<=4096:p.error('phase-samples must be between 8 and 4096')
    project=args.project.resolve()
    def full(q):return (q if q.is_absolute() else project/q).resolve()
    output=full(args.out)
    if output.exists():p.error('Refusing to overwrite existing report: '+str(output))
    resp=full(args.response_root);reconcile=full(args.reconcile_dir)
    from Numerical.running.T3CoupledQuarticContext import CandidateQuarticContext
    report_file=reconcile/'report.json'
    report=json.loads(report_file.read_text(encoding='utf-8-sig'))
    if report.get('schema')!='t3-worker-b-batch-hhss-reconciliation-v1':
        raise ValueError('Wrong reconciliation report schema')
    results={}
    for model in models:
        ctx=CandidateQuarticContext.load(project,resp,model)
        row=report['results'][model]
        file=reconcile/f'{model}_reconciliation.npz'
        if row.get('npz_sha256')!=sha(file):raise ValueError(f'{model}: reconciliation cache hash differs')
        with np.load(file,allow_pickle=False) as f:
            source=np.asarray(f['source_direction_columns'],float)
        results[model]=diagnose(model,ctx.directions,ctx.columns[:,:4,:4,4:,4:],source,samples=args.phase_samples)
        t=results[model]['directions']['lambdaT3.real']
        print(f'{model}: lambdaT3.real fixed-axis spectral obstruction={t["normalized_spectral_lower_bound"]:.6g}; '
              f'best scanned source phase spectral distance={results[model]["lambdaT3_source_phase_scan"]["real"]["best_sampled_phase"]["spectral_distance"]:.6g}')
    output.parent.mkdir(parents=True,exist_ok=True)
    doc={'schema':SCHEMA,'status':'SpectralInvariantsEvaluated',
         'response_root':str(resp),'reconciliation_root':str(reconcile),
         'results':results,'no_production_changes':True,'physical_quartic_certified':False}
    output.write_text(json.dumps(doc,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Report:',output)
    return 0


if __name__=='__main__':raise SystemExit(main())
