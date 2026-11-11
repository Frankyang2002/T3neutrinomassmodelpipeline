"""Read-only A(-4)/D(-2) fermion-hard source audit; NOT a physical EFT1 Wilson match.

Run at repository root: python scripts/WorkerAFiniteChannelAudit.py --repo . --output ...
Requires existing Batch1 and Batch2 outputs, does not modify their source files.
"""
from __future__ import annotations
import argparse
from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
import numpy as np

BASE = Path('output/parallel_t3/full_general_batch2')
B1 = Path('output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1')


def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))


def decode(x):
    if isinstance(x, dict) and set(x) == {'real', 'imag'}:
        return complex(x['real'], x['imag'])
    if isinstance(x,list):
        return [decode(a) for a in x]
    return x


def complex_json(x):
    a=np.asarray(x,complex)
    return {'real':a.real.tolist(), 'imag':a.imag.tolist()}


def fingerprint(path):
    return sha256(path.read_bytes()).hexdigest()


def channel(t):
    s=t['raw_inputform']
    assert t['sector']=='CandidateLeptonNumberViolatingLLHH'
    assert s.count('GammaCC')==1
    yuk='y2Conj' if 'Coupling[y2Conj,' in s else 'y2'
    quartic='lambdaT3ExtraQ020101I1' if 'Coupling[lambdaT3ExtraQ020101I1,' in s else 'lambdaT3'
    chirality='holomorphic_HH' if 'Proj[-1]' in s and 'Proj[1]' not in s else 'anti_holomorphic_HdagHdag'
    assert chirality in ('holomorphic_HH','anti_holomorphic_HdagHdag')
    kind='log' if 'Log[' in s else ('pole' if 'ϵ' in s or 'ε' in s else 'constant')
    assert s.count('Coupling[y1,')==1 and s.count('Coupling['+yuk+',')==1
    return (yuk,quartic,chirality,kind)


def audit(root, model):
    raw=root/BASE/'thesis_full_general_batch2_raw_v1'/model
    laurent=root/BASE/'thesis_full_general_batch2_laurent_v2'/model
    fpath=root/B1/model/'fermion_matching_input.json'
    spath=root/B1/model/'scalar_matching_input.json'
    ledger_path=raw/'one_loop_operator_ledger.json'
    lpath=laurent/'laurent_audit.json'
    raw_source=raw/'native_matching/data/eft1_after_F_one_loop_inputform.txt'
    for p in (fpath,spath,ledger_path,lpath,raw_source):
        if not p.is_file(): raise FileNotFoundError(p)
    led=read(ledger_path); lau=read(lpath); f=read(fpath); s=read(spath)
    assert led['status']=='UnrenormalisedOneLoopTermsEnumerated'
    assert led['source_sha256']==fingerprint(raw_source)
    assert lau['source_sha256']==fingerprint(raw_source)
    assert lau['explicit_laurent_reconstruction_exact'] is True
    assert lau['uv_hard_minus_eft1_soft_subtraction_performed'] is False
    assert lau['renormalised_MSbar_Wilson_boundary_computed'] is False
    assert f['model']==s['model']==model
    assert f['all_heavy_fermions_to_integrate_out_together'] is True
    assert s['both_heavy_scalar_multiplets_to_integrate_out_together'] is True
    assert f['scale_GeV']==1e5 and s['scale_GeV']==1e3
    terms=[t for t in led['terms'] if t['sector']=='CandidateLeptonNumberViolatingLLHH']
    inv=Counter(channel(t) for t in terms)
    baseline={(y,q,c,k) for y,q in [('y2','lambdaT3')] for c in ('holomorphic_HH','anti_holomorphic_HdagHdag') for k in ('constant','pole','log')}
    expected=set(baseline)
    if model=='D_-2':
        expected.update((y,q,c,k) for y,q in [('y2Conj','lambdaT3ExtraQ020101I1')] for c in ('holomorphic_HH','anti_holomorphic_HdagHdag') for k in ('constant','pole','log'))
    if set(inv)!=expected or any(v!=1 for v in inv.values()):
        raise RuntimeError(f'{model}: LLHH channel inventory changed: {inv}')
    v=f['native_uv_rgbeta_couplings_at_muF']
    M=np.asarray(decode(v['MF']),complex);Y1=np.asarray(decode(v['y1']),complex);Y2=np.asarray(decode(v['y2']),complex)
    assert M.shape==Y1.shape==Y2.shape==(3,3)
    offdiag=float(np.linalg.norm(M-np.diag(np.diag(M)))/np.linalg.norm(M))
    if offdiag>1e-8:raise RuntimeError(f'{model}: MF mass-basis unverified: {offdiag}')
    masses=np.real(np.diag(M))
    if np.any(masses<=0) or np.max(abs(np.imag(np.diag(M))))>1e-8*np.max(masses):
        raise RuntimeError('non-positive or complex native fermion masses')
    lam=complex(decode(v['lambdaT3']))
    g=2/(math.sqrt(3)*16*math.pi**2)
    w=(Y1.conj()/masses[None,:])@Y2.conj().T
    pole=g*lam*w
    formal_hard_regular=g*lam*((Y1.conj()*(1+np.log((1e5/masses)**2))[None,:]/masses[None,:])@Y2.conj().T)
    record={
        'model':model,'status':'source_verified_direct_channel_diagnostic_NOT_finite_matching',
        'source_sha256':{'raw_matchete_one_loop':fingerprint(raw_source),'ledger':fingerprint(ledger_path),'laurent_audit':fingerprint(lpath),'fermion_input':fingerprint(fpath),'scalar_input':fingerprint(spath)},
        'all_fermions_integrated_at_GeV':1e5,'both_scalars_retained_through_EFT1':True,
        'total_one_loop_operators':led['term_count'],'sectors':led['sector_counts'],
        'direct_llhh_terms':len(terms),'restricted_Y1_Y2_llhh_terms':len(baseline),
        'excluded_Y2Conj_llhh_terms':len(terms)-len(baseline),
        'hard_mass_eigenvalues_diagonal_approx_GeV':masses.tolist(),
        'MF_offdiagonal_relative_norm':offdiag,
        'restricted_Y2Conj_norm_at_F':float(np.linalg.norm(np.asarray(decode(v.get('y2Conj',np.zeros((3,3)))),complex))),
        'raw_direct_hard_pole_ordered_GeVinv':complex_json(pole),
        'raw_direct_hard_regular_ordered_GeVinv':complex_json(formal_hard_regular),
        'raw_direct_hard_pole_Frobenius_GeVinv':float(np.linalg.norm(pole)),
        'raw_direct_hard_regular_Frobenius_GeVinv':float(np.linalg.norm(formal_hard_regular)),
        'hard_regular_is_NOT_MSbar_boundary':True,
        'source_eft1_soft_subtraction_performed':False,
        'gammaCC_to_Weyl_sign_and_normalization_certified':False,
        'complete_operator_basis_renormalization_certified':False,
        'finite_fermion_boundary_available':False,
    }
    if model=='A_-4':
        sm=s['scalar_quadratic_terms_GeV2'];a=float(sm['mS1Sq']);b=float(sm['mS2Sq'])
        # Existing a4_direct_hard_soft.py formula: this is a formal direct HARD+SOFT
        # comparison at fixed fermion-scale Yukawas, NOT a boundary Wilson coefficient.
        if a<=0 or b<=0 or a==b:raise RuntimeError('invalid diagonal scalar masses')
        logs=(a*np.log(a/masses**2)-b*np.log(b/masses**2))/(a-b)
        candidate=g*lam*((Y1.conj()*logs[None,:]/masses[None,:])@Y2.conj().T)
        record['A4_formal_direct_hard_soft_remainder_GeVinv']=complex_json(candidate)
        record['A4_formal_direct_hard_soft_remainder_norm_GeVinv']=float(np.linalg.norm(candidate))
        record['A4_formal_remainder_provenance']='Formula in Numerical/parallel_t3/complete/a4_direct_hard_soft.py; uses F Yukawas and S scalar masses; not an EFT1 Wilson boundary'
    else:
        record['D2_extra_channel']='Removed by setting independent y2Conj=0; no physical finite subtraction inferred'
        record['D2_mass_basis']='additional scalar self-bilinear requires mass-basis-aware full soft matching; not done'
    return record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path('.'))
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    root=args.repo.resolve()
    report={'schema':'worker-a-direct-channel-audit-v1',
            'one_loop_MSbar_boundary_completed':False,
            'models':{m:audit(root,m) for m in ('A_-4','D_-2')},
            'unresolved':['Independent full-basis EFT1 soft subtraction and counterterms',
                          'GammaCC-to-Weyl sign/normalization and CG dictionary',
                          'D(-2) self-bilinear mass-basis soft terms'],
            'note':'Source-backed numerical diagnostics only; never feed any raw or formal remainder into coupled solver.'}
    target=args.output
    if not target.is_absolute():target=root/target
    target.parent.mkdir(parents=True,exist_ok=True)
    with target.open('x',encoding='utf-8') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    for m,v in report['models'].items():
        print(f"{m}: LLHH {v['direct_llhh_terms']}, |raw pole|={v['raw_direct_hard_pole_Frobenius_GeVinv']:.8e} GeV^-1, |raw regular|={v['raw_direct_hard_regular_Frobenius_GeVinv']:.8e} GeV^-1")
        if m=='A_-4':print(f"  |formal A4 hard+soft remainder|={v['A4_formal_direct_hard_soft_remainder_norm_GeVinv']:.8e} GeV^-1 (NOT physical C5)")
    print('Wrote',target)

if __name__=='__main__':main()
