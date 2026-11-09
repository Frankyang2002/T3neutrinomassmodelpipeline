"""Read-only first-cause gauge Ward gate for the fermion-first rough T3 pipeline.

Read saved full Wilson at MF/MS, saved canonical Stream-1 background at MF/MS,
and optionally recompute the EXISTING candidate beta at MF (not a new RGE).
Use the project's actual Weyl/real-scalar generators. No field transformations,
SU(2) averaging, scalar matching, or artificial Wilson boundary is performed.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import json
import math
import numpy as np

from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.wilson.adm import beta, total_beta
from Numerical.parallel_t3.threshold.models import normalize_model
from Numerical.parallel_t3.threshold.pipeline import _load_pair, W_FORMAT
from Numerical.parallel_t3.threshold.rough_origin import source_charge_analysis
from Numerical.parallel_t3.threshold.matching import extract_lepton_sectors

FORMAT='t3-rough-stream3-first-cause-ward-gate-v6'
GAUGES=('SU2','U1','SU3')
GATE_TOLERANCE=1e-8


def file_sha(path:Path)->str:
    return sha256(Path(path).read_bytes()).hexdigest()


def _nonzero_norm(a):
    return float(np.linalg.norm(np.asarray(a)))


def ward_numerator(tensor:np.ndarray, tf:np.ndarray | None,
                   ts:np.ndarray)->np.ndarray:
    """Covariant variation of Y_ija, C_ijab, or lambda_abcd.

    Source convention: all tensor indices co-transform with the same signed
    Hermitian generators used by Wilson ADM. Do not conjugate C or Y here.
    """
    x=np.asarray(tensor)
    s=np.asarray(ts)
    if x.ndim==3:
        if tf is None or x.shape[:2]!=(len(tf),len(tf)) or x.shape[2]!=len(s):
            raise ValueError('Yukawa axes must be [nW,nW,nReal]')
        return (np.einsum('ki,kja->ija',tf,x,optimize=True)
                +np.einsum('kj,ika->ija',tf,x,optimize=True)
                +np.einsum('ca,ijc->ija',s,x,optimize=True))
    if x.ndim==4 and tf is not None:
        if x.shape[:2]!=(len(tf),len(tf)) or x.shape[2:]!=(len(s),len(s)):
            raise ValueError('Wilson axes must be [nW,nW,nReal,nReal]')
        return (np.einsum('ki,kjab->ijab',tf,x,optimize=True)
                +np.einsum('kj,ikab->ijab',tf,x,optimize=True)
                +np.einsum('ca,ijcb->ijab',s,x,optimize=True)
                +np.einsum('cb,ijac->ijab',s,x,optimize=True))
    if x.ndim==4 and tf is None:
        if x.shape!=(len(s),)*4:raise ValueError('Quartic axes must be [nReal]^4')
        return (np.einsum('ka,kbcd->abcd',s,x,optimize=True)
                +np.einsum('kb,akcd->abcd',s,x,optimize=True)
                +np.einsum('kc,abkd->abcd',s,x,optimize=True)
                +np.einsum('kd,abck->abcd',s,x,optimize=True))
    raise ValueError('Supported tensors: 3-index Y, 4-index Wilson, 4-index quartic')


def gauge_ward_report(x:np.ndarray, basis, *, sector:str)->dict:
    """Full tensor Frobenius Ward residuals; zero tensors remain uninformative."""
    norm=_nonzero_norm(x)
    out={'tensor_norm':norm,'status':'zero_tensor' if norm<1e-100 else 'checked',
         'ward_relative':{}, 'ward_absolute':{}}
    if sector not in ('wilson','yukawa','quartic'):raise ValueError('Unknown sector')
    for name in GAUGES:
        measurements=[]
        for tf,ts in basis.generators()[name]:
            value=ward_numerator(x,tf if sector!='quartic' else None,ts)
            measurements.append(_nonzero_norm(value))
        # Report max generator (not RMS), consistent with prior Ward table.
        absolute=max(measurements,default=0.)
        out['ward_relative'][name]=None if norm<1e-100 else absolute/norm
        out['ward_absolute'][name]=absolute
    out['max_relative']=None if norm<1e-100 else max(out['ward_relative'].values())
    out['passed']=None if norm<1e-100 else out['max_relative']<GATE_TOLERANCE
    return out


def read_background(npz:Path, meta:Path, basis, *, source_sha:str | None=None)->dict:
    npz,meta=Path(npz),Path(meta)
    document=json.loads(meta.read_text(encoding='utf-8-sig'))
    digest=file_sha(npz)
    if document.get('model')!=basis.model:
        raise ValueError('Stream1 canonical background model mismatch')
    if document.get('npz_sha256')!=digest and document.get('tensor_npz_sha256')!=digest:
        raise ValueError('Stream1 canonical background NPZ checksum differs from manifest')
    if source_sha is not None and digest!=source_sha:
        raise ValueError('Stream1 background no longer matches the background used by Stream2')
    if tuple(document.get('fermion_labels',())) != basis.fermion_labels or \
       tuple(document.get('scalar_labels',())) != basis.scalar_labels:
        raise ValueError('Stream1 background basis labels differ from current Wilson basis')
    with np.load(npz,allow_pickle=False) as d:
        arrays={key:np.asarray(d[key]).copy() for key in d.files}
    mu=np.asarray(arrays['mu_gev'],float)
    n=len(mu)
    expected=(n,basis.nf,basis.nf,basis.ns)
    if n<2 or arrays['Y'].shape!=expected or arrays['lambda4'].shape!=(n,)+(basis.ns,)*4:
        raise ValueError('Canonical Stream1 background shape mismatch')
    if not np.isfinite(mu).all() or len(np.unique(mu))!=n or np.any(mu<=0):
        raise ValueError('Stream1 background scales are invalid')
    for name in ('Y','lambda4','g2','gY','g3'):
        if not np.isfinite(arrays[name]).all():raise ValueError('Nonfinite '+name)
    iMF=int(np.argmax(mu));iMS=int(np.argmin(mu))
    if not np.isclose(mu[iMF],1e5,rtol=1e-9) or not np.isclose(mu[iMS],1e3,rtol=1e-9):
        raise ValueError('This rough benchmark expects canonical endpoints MF=1e5,MS=1e3 GeV')
    return {'meta':document,'sha256':digest,'arrays':arrays,'iMF':iMF,'iMS':iMS}


def _llhh(c,model:str,basis):
    indices=np.asarray([list(pair) for pair in basis.leptons],int)
    h,s=extract_lepton_sectors(c,indices,np.arange(4,dtype=int),
        np.arange(4,basis.ns,dtype=int),number_of_real_scalars=basis.ns)
    return source_charge_analysis(h)


def _origin_flags(c0,y0,q0,b0, *, tol=GATE_TOLERANCE)->list[str]:
    bad=[]
    for label,item in (('MF_wilson_boundary',c0),('MF_SM_Yukawa',y0),('MF_full_quartic',q0),
                       ('MF_full_beta',b0)):
        if item is None:continue
        ward=item['ward_relative']
        for group in GAUGES:
            residual=ward[group]
            if residual is not None and residual>tol:
                bad.append(f'{label}_{group}_WardViolation')
    return bad


def check_portal_units(model:str,basis=None)->dict:
    """Independently check source-coded, conditional Stream1 quartic units.

    No attempt to replace their Clebsch tensors or fix the native dictionary.
    For a complex portal test two real combinations (phase=1 and i).
    """
    from Numerical.parallel_t3.rough_inputs.builder import quartic_unit_tensors
    basis=basis or build_basis(model)
    source=quartic_unit_tensors(model)
    results={}
    for name,val in source.items():
        v=np.asarray(val)
        phases=(('real',1.),('imaginary',1j)) if np.max(abs(v.imag))>1e-14 else (('real',1.),)
        part={}
        for phase_label,phase in phases:
            realtensor=np.asarray(np.real(phase*v),float)
            if _nonzero_norm(realtensor)<1e-100: continue
            part[phase_label]=gauge_ward_report(realtensor,basis,sector='quartic')
        results[name]=part
    bad=[name for name,parts in results.items() if any(
        part['passed'] is False for part in parts.values())]
    return {'status':'SourceNativeCGWardTestNotNormalizationCertification',
            'quartic_component_units':results,'noninvariant_unit_names':bad,
            'lambdaT3_present':'lambdaT3' in results}


def audit_model(root:Path,model:str,*,with_beta:bool=False,with_units:bool=True,
                input_run_id:str='rough_fermion_first_v1')->dict:
    root=Path(root)
    model=normalize_model(model)
    slug=model.replace('_-','_m').replace('_+','_p')
    basis=build_basis(model)
    wpath=root/'output/parallel_t3/rough_wilson'/input_run_id/slug/'wilson.json'
    bpath=root/'output/parallel_t3/rough_inputs'/input_run_id/'a'/slug/'background.json'
    w,wa,wjsonsha,wnpzsha=_load_pair(wpath,W_FORMAT)
    if w.get('model')!=model:raise ValueError('Wilson manifest model mismatch')
    c0=np.asarray(wa['C_psi2phi2_tree_mf'],complex)
    cMS=np.asarray(wa['C_psi2phi2_ms'],complex)
    expected=(basis.nf,basis.nf,basis.ns,basis.ns)
    if c0.shape!=expected or cMS.shape!=expected:
        raise ValueError('Stream2 full Wilson shape/basis mismatch')
    src_ref=(w.get('input',{}).get('file_sha256',{}).get('background_npz',{})
             .get('sha256'))
    bg=read_background(bpath.with_suffix('.npz'),bpath,basis,source_sha=src_ref)
    arr=bg['arrays']
    a=arr['Y'];lam=arr['lambda4']
    i0,i1=bg['iMF'],bg['iMS']
    if with_units:
        # Unit inventory exposes source-code CG failure independently of
        # the values/zeros used in the saved benchmark.
        units=check_portal_units(model,basis)
    else:units={'status':'not_requested'}
    mf={'Wilson':gauge_ward_report(c0,basis,sector='wilson'),
        'Yukawa':gauge_ward_report(a[i0],basis,sector='yukawa'),
        'quartic':gauge_ward_report(lam[i0],basis,sector='quartic')}
    ms={'Wilson':gauge_ward_report(cMS,basis,sector='wilson'),
        'Yukawa':gauge_ward_report(a[i1],basis,sector='yukawa'),
        'quartic':gauge_ward_report(lam[i1],basis,sector='quartic')}
    beta_saved = (gauge_ward_report(np.asarray(wa['beta_C_at_MF_16pi2']),basis,
                                    sector='wilson') if 'beta_C_at_MF_16pi2' in wa else None)
    beta_saved_ms = (gauge_ward_report(np.asarray(wa['beta_C_at_MS_16pi2']),basis,
                                       sector='wilson') if 'beta_C_at_MS_16pi2' in wa else None)
    beta_info=None
    if with_beta:
        couplings=tuple((float(np.real(arr[name][i0])),basis.generators()[group])
                        for name,group in [('g2','SU2'),('gY','U1'),('g3','SU3')])
        channels,backend=beta(c0,a[i0],lam[i0],couplings,
            variant=w.get('adm_variant','fonseca_eq22_six_term'))
        summed=total_beta(channels)
        prior=np.asarray(wa['beta_C_at_MF_16pi2'])
        mismatch=float(np.linalg.norm(summed-prior)/max(np.linalg.norm(prior),
                    np.linalg.norm(summed),1e-100))
        beta_info={'backend':backend,'variant':w.get('adm_variant'),
                   'stored_beta_relative_disagreement':mismatch,
                   'total':gauge_ward_report(summed,basis,sector='wilson'),
                   'total_LLHH':_llhh(summed,model,basis),
                   'by_channel':{key:{'full_Ward':gauge_ward_report(value,basis,sector='wilson'),
                          'LLHH':_llhh(value,model,basis)} for key,value in channels.items()}}
    flags=_origin_flags(mf['Wilson'],mf['Yukawa'],mf['quartic'],
                        beta_info['total'] if beta_info else beta_saved)
    # Note that even perfect MF Ward identities do not certify the native
    # Matchete <-> RGBeta dictionary. They only identify this failure stage.
    return {'format':FORMAT,'model':model,'status':'SourceWardFirstCauseDiagnosticOnly',
       'source_sha256':{'Stream2_manifest':wjsonsha,'Stream2_full_wilson':wnpzsha,
                        'Stream1_background':bg['sha256']},
       'Stream1_matching_provenance':{k:w.get(k) for k in
                       ('matching_context_id','matched_seed_id','lagrangian_source_id','scalar_basis_id')},
       'MF':mf,'MS':ms,'MF_LLHH':_llhh(c0,model,basis),
       'MS_LLHH':_llhh(cMS,model,basis),
       'saved_candidate_beta_MF':beta_saved,'saved_candidate_beta_MS':beta_saved_ms,
       'candidate_beta_MF':beta_info,'source_quartic_units':units,
       'first_cause_candidates':flags,
       'not_a_proof_of_native_source_normalization':True,
       'no_source_arrays_or_code_modified':True,
       'all_fermions_out_before_both_scalars':True}
