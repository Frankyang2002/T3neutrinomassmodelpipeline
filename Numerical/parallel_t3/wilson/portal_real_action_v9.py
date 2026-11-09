"""Workstream A v9: 58 native-portal real-action consistency audit (CONDITIONAL).

Compares 1/2 PlusHc[phase*MatcheteCG] against a constructed real-scalar
RGBeta 6*Matrix / 4! action, using native v5 registrations and the *abstract*
fundamental Fierz / S2 Casimir adjoint orientation of v6/v8.

NO native RGBeta Lam/projector evaluated; NO numerical native RGBeta S2 frame;
NO full scalar action, mass eigenbasis, spinor, or Wilson certification.
"""
from __future__ import annotations

import itertools
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any

import numpy as np

FORMAT = 't3-parallel-a-higgs-portal-real-action-conditional-v9'
STATUS = ('NativePortal58RealActionAlgebraVerified_'
          'NativeRGBetaLamProjectorsAndFullScalarNormalizationUnverifiedConditional')
TOL = 2.e-10
EXPECTED_MODELS = frozenset({
    'A_-4','A_-2','A_+0','B_-3','B_-1','B_+1','C_-3','C_-1','C_+1',
    'D_-2','D_+0','D_+2','E_-4','E_-2','E_+0','E_+2',
})
DIMENSIONS = {'A':(1,3,2),'B':(2,2,1),'C':(2,2,3),
              'D':(3,1,2),'E':(3,3,2)}
PERMS = tuple(itertools.permutations(range(4)))


class V9EvidenceError(ValueError):
    """Missing/stale/contradictory native evidence or failed action test."""


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise V9EvidenceError('Workstream A v9: '+reason)


def relative(a: Any, b: Any) -> float:
    a, b = np.asarray(a), np.asarray(b)
    return float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(a),np.linalg.norm(b)))


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def read(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj,dict),str(path)+' must contain a JSON object')
    return obj


def identity(d: int) -> np.ndarray:
    return np.einsum('ij,ab->ijab',np.eye(2),np.eye(d)).astype(complex)


def fundamental_fierz() -> np.ndarray:
    """Independently source-grounded RGBeta fundamental completeness identity."""
    return np.array([
        (int(i==v and u==j)-.5*int(i==j and u==v))/2
        for i in range(2) for j in range(2)
        for u in range(2) for v in range(2)
    ],dtype=complex).reshape(2,2,2,2)


def _spectral_values(k: np.ndarray) -> np.ndarray:
    mat = k.transpose(0,2,1,3).reshape(6,6)
    require(relative(mat,mat.conj().T)<TOL,'triplet adjoint is non-Hermitian')
    return np.linalg.eigvalsh(mat)


def spectral_triplet_adjoint(columns: list[np.ndarray]) -> tuple[np.ndarray,int,float]:
    """Select the sign without reading v6/v7's Matchete-fitted map."""
    require(len(columns)==2 and all(x.shape==(2,2,3,3) for x in columns),
            'triplet portal requires two native tensors')
    ident=identity(3)
    orth=[]
    for x in columns:
        a=np.vdot(ident,x)/np.vdot(ident,ident)
        require(abs(a.imag)<TOL,'triplet singlet coefficient is complex')
        orth.append(x-a.real*ident)
    pivot=next((x for x in orth if np.linalg.norm(x)>1.e-9),None)
    require(pivot is not None,'triplet traceless portal missing')
    candidate=np.sqrt(3.)*pivot/np.linalg.norm(pivot)
    expected=np.array([-1.,-1.,.5,.5,.5,.5])
    errors={sign:relative(_spectral_values(sign*candidate),expected)
            for sign in (1,-1)}
    accepted=[s for s in (1,-1) if errors[s]<TOL]
    require(len(accepted)==1,'S2 total-Casimir spectrum did not fix a unique sign')
    sign=accepted[0]
    require(errors[-sign]>.1,'S2 opposite-sign negative control failed')
    return sign*candidate,sign,errors[sign]


def _real_projection(d: int) -> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    """x=[H_1R,H_1I,...,S_1R,S_1I,...], z=(R+iI)/sqrt(2).

    The 1/sqrt(2) factors are *not* in these matrices: RGBeta's complex
    scalar coupling prefactor 6 / 4! cancels four such factors exactly.
    """
    n=4+2*d
    h=np.zeros((2,n),complex);s=np.zeros((d,n),complex)
    for i in range(2):
        h[i,2*i]=1.;h[i,2*i+1]=1j
    for i in range(d):
        s[i,4+2*i]=1.;s[i,5+2*i]=1j
    return h.conj(),h,s.conj(),s


def canonical_lam(t: np.ndarray) -> np.ndarray:
    """Conditional RGBeta real Lam *constructed* from source 6*Matrix rule.

    Quartic tensor is fully symmetric in real axes, action is Lam*x^4/24.
    RGBeta's native Lam/SDelS projector is NOT called by this Python code.
    """
    require(t.ndim==4 and t.shape[:2]==(2,2) and t.shape[2]==t.shape[3]
            and t.shape[2] in (1,2,3),'invalid portal scalar dimensions')
    d=t.shape[2]
    require(relative(t,t.conj().transpose(1,0,3,2))<TOL,
            'non-Hermitian portal tensor')
    hb,h,sb,s=_real_projection(d)
    unsym=np.einsum('ijab,ip,jq,ar,bs->pqrs',t,hb,h,sb,s,optimize=True)
    sym=sum(np.transpose(unsym,p) for p in PERMS)/24.
    lam=6*sym
    require(float(np.linalg.norm(lam.imag))<TOL,
            'real scalar Lam has residual imaginary component')
    return lam.real


def complex_fields(x: np.ndarray,d: int) -> tuple[np.ndarray,np.ndarray]:
    require(x.shape==(4+2*d,),'real scalar vector wrong dimension')
    return ((x[0:4:2]+1j*x[1:4:2])/np.sqrt(2),
            (x[4::2]+1j*x[5::2])/np.sqrt(2))


def complex_action(t: np.ndarray,x: np.ndarray) -> float:
    h,s=complex_fields(x,t.shape[2])
    y=np.einsum('ijab,i,j,a,b->',t,h.conj(),h,s.conj(),s,optimize=True)
    require(abs(y.imag)<TOL*max(1.,abs(y)),'complex portal action not real')
    return float(y.real)


def lam_action(lam: np.ndarray,x: np.ndarray) -> float:
    require(lam.shape==(len(x),)*4,'Lam real axis mismatch')
    return float(np.einsum('ijkl,i,j,k,l->',lam,x,x,x,x,optimize=True)/24.)


def source_half_plus_hc(t_raw: np.ndarray,phase: complex,x: np.ndarray) -> float:
    h,s=complex_fields(x,t_raw.shape[2])
    raw=phase*np.einsum('ijab,i,j,a,b->',t_raw,h.conj(),h,s.conj(),s,optimize=True)
    return float(raw.real) # 1/2 (raw + conjugate)


def _source_record_ok(r: dict,which: int,d:int,k:int) -> None:
    require(r.get('which')==which and r.get('invariant')==k,
            'native Matchete portal index/order changed')
    require(r.get('coupling')==f'lambdaH{which}Inv{k}',
            'Matchete portal symbol mismatch')
    require(r.get('dimensions')==[2,2,d,d] and
            r.get('phase_inside_plus_hc') in ('1','I'),
            'Matchete portal dimensions/PlusHc phase unsupported')


def _rgb_record_ok(r: dict,which:int,d:int,kind:str) -> None:
    symbol=f'lambdaH{which}'+('Adj' if kind=='adjoint' else '')
    require(r.get('symbol')==symbol and r.get('which')==which and
            r.get('kind')==kind and r.get('dimensions')==[2,2,d,d],
            'native RGBeta portal registration mismatch')
    require(r.get('native_self_conjugate') is True and
            r.get('native_permutation_count')==24,
            'native RGBeta unique arrangements or self-conjugacy unsupported')
    norm=re.sub(r'\s+','',r.get('native_coupling_inputform',''))
    require(norm==f'6*Matrix[{symbol}][##1]&',
            'RGBeta portal 6*Matrix coupling factor missing')
    inv=r.get('native_invariant_inputform','')
    require(isinstance(inv,str) and
            (('tGen[' in inv) if kind=='adjoint' else ('del[' in inv)),
            'RGBeta native invariant symbolic form changed')
    require(r.get('native_numerical_generator_components_exported') is False,
            'unexpected numerical RGBeta generator claim')


def _expected_v6_matrix(v6_sector:dict,d:int) -> np.ndarray:
    kind={1:'singlet',2:'doublet',3:'triplet'}[d]
    require(v6_sector.get('sector')==kind,'v6 sector type mismatch')
    key=('candidate_rgbeta_from_matchete_real_matrix' if d==3
         else 'rgbeta_from_matchete_real_matrix')
    require(key in v6_sector,'v6 coupling dictionary schema changed')
    result=np.asarray(v6_sector[key],dtype=float)
    n=1 if d==1 else 2
    require(result.shape==(n,n) and np.all(np.isfinite(result)),
            'v6 coupling transformation malformed')
    return result


def audit_sector(raw: list[np.ndarray], phased: list[np.ndarray],
                 mat_records: list[dict],rgb_records: list[dict],
                 v6_sector: dict,*, tag:str,which:int,d:int,
                 v8_sector:dict|None=None,
                 seed:int=13241,samples:int=5) -> dict:
    """Native CG component action vs constructed RGBeta action, 1/2 PlusHc.

    Inputs raw/phased must originate from the v5 native Matchete accessor;
    full provenance and re-execution checks are performed by audit_from_paths.
    """
    n=1 if d==1 else 2
    require(len(raw)==len(phased)==len(mat_records)==len(rgb_records)==n,
            f'{tag}/S{which}: portal sector cardinality mismatch')
    for k,r in enumerate(mat_records,1):_source_record_ok(r,which,d,k)
    for k,r in enumerate(rgb_records):
        _rgb_record_ok(r,which,d,'singlet' if k==0 else 'adjoint')
    for r,p,q in zip(mat_records,raw,phased):
        require(p.shape==q.shape==(2,2,d,d),f'{tag}/S{which}: CG tensor shape')
        phase=1j if r['phase_inside_plus_hc']=='I' else 1.
        require(relative(q,(phase*p+(phase*p).conj().transpose(1,0,3,2))/2)<TOL,
                f'{tag}/S{which}: Matchete 1/2 PlusHc projection mismatch')
    ident=identity(d)
    if d==1:
        generator=None; sign=None; spectral=None; target=[ident]
    elif d==2:
        generator=fundamental_fierz();sign=None;spectral=None;target=[ident,generator]
    else:
        generator,sign,spectral=spectral_triplet_adjoint(phased)
        target=[ident,generator]
        require(v8_sector is not None and
                v8_sector.get('source_CG_sign_selected_without_using_v6_map') is True and
                v8_sector.get('sign_relative_to_first_native_traceless_CG')==sign and
                float(v8_sector.get('selected_spectral_checks',{}).get('max_residual',1.))<TOL and
                v8_sector.get('RGBeta_numeric_triplet_frame_certified') is False,
                f'{tag}/S{which}: saved v8 independent spectral sign mismatch')
    gram=np.asarray([[np.vdot(x,y) for y in target] for x in target])
    require(relative(gram,np.diag(np.diag(gram)))<TOL,
            'singlet and adjoint not orthogonal')
    convert=np.array([[np.vdot(t,x)/np.vdot(t,t) for x in phased]
                      for t in target],dtype=complex)
    require(float(np.linalg.norm(convert.imag))<TOL,
            f'{tag}/S{which}: complex couplings in real portal dictionary')
    convert=convert.real
    expected=_expected_v6_matrix(v6_sector,d)
    require(relative(convert,expected)<TOL,
            f'{tag}/S{which}: source action contradicts verified v6 transformation')
    require(np.linalg.matrix_rank(convert,tol=1.e-10)==n,
            f'{tag}/S{which}: rank-deficient portal map')
    tensor_error=max(relative(x,sum(convert[k,j]*target[k] for k in range(n)))
                     for j,x in enumerate(phased))
    require(tensor_error<TOL,'Matchete portal tensor not in RGBeta invariant span')
    rng=np.random.default_rng(seed)
    max_errors={'source_vs_complex':0.,'native_CG_vs_constructed_Lam':0.,
                'RGBeta_conditional_vs_constructed_Lam':0.,
                'native_CG_vs_RGBeta_conditional':0.}
    for _ in range(samples):
        x=rng.normal(size=4+2*d)
        coefficients=rng.uniform(-1.7,1.7,size=n)  # real source coupling
        rgb_coeff=convert@coefficients
        mat_t=sum(c*t for c,t in zip(coefficients,phased))
        rgb_t=sum(c*t for c,t in zip(rgb_coeff,target))
        native_source=sum(c*source_half_plus_hc(r,1j if rec['phase_inside_plus_hc']=='I' else 1.,x)
                          for c,rec,r in zip(coefficients,mat_records,raw))
        complex_mat=complex_action(mat_t,x)
        rgb_complex=complex_action(rgb_t,x)
        lam_mat=lam_action(canonical_lam(mat_t),x)
        lam_rgb=lam_action(canonical_lam(rgb_t),x)
        scores={'source_vs_complex':relative(native_source,complex_mat),
                'native_CG_vs_constructed_Lam':relative(native_source,lam_mat),
                'RGBeta_conditional_vs_constructed_Lam':relative(rgb_complex,lam_rgb),
                'native_CG_vs_RGBeta_conditional':relative(native_source,rgb_complex)}
        for key,value in scores.items():max_errors[key]=max(max_errors[key],value)
    worst=max(tensor_error,*max_errors.values(),relative(convert,expected),
              spectral or 0.)
    require(worst<TOL,f'{tag}/S{which}: action normalisation residual {worst:.4g}')
    return {
        'model':tag,'which':which,'dimension':d,'native_matchete_couplings':[x['coupling'] for x in mat_records],
        'rgbeta_registered_couplings':[x['symbol'] for x in rgb_records],
        'native_Matchete_phase_inside_half_PlusHc':[x['phase_inside_plus_hc'] for x in mat_records],
        'conditional_rgbeta_from_matchete_real_matrix':convert.tolist(),
        'source_CG_tensor_to_invariant_span_residual':tensor_error,
        'constructed_real_Lam_action_residuals':max_errors,
        'max_relative_residual':worst,
        'triplet_spectral_sign_relative_to_native_CG':sign,
        'triplet_spectral_residual':spectral,
        'rgbeta_quartic_prefactor_from_native_registration':'6*Matrix[coupling]',
        'rgbeta_real_action_convention':'constructed Lam_real[x^4]/4!; z=(R+iI)/sqrt(2)',
        'native_RGBeta_Lam_projector_executed':False,
        'native_RGBeta_S2_component_frame_certified':False,
        'physical_coupling_dictionary_certified':False,
    }


def audit(mat:dict,rgb:dict,v6:dict,v8:dict,*,samples:int=5)->dict:
    """Uses verified source exports; caller revalidates their provenance."""
    require(set(mat.get('models',{}))==set(rgb.get('models',{}))==
            set(v6.get('models',{}))==EXPECTED_MODELS,'16 native portal models missing')
    require(v6.get('triplet_sectors_Pauli_algebra_conditional')==14 and
            v8.get('triplet_sectors_basis_invariant_spectrum_verified')==14,
            'v6/v8 triplet audit not complete')
    require(v6.get('scope',{}).get('native_RGBeta_triplet_intertwiner_certified') is False,
            'unexpected v6 native S2 certification')
    require(v8.get('certification',{}).get('native_RGBeta_numerical_S2_frame_certified') is False,
            'unexpected v8 native numerical S2 certification')
    count=0;sector_count=0;anti=0;triplets=0;max_error=0.;models={}
    # Import source-validated CG decoder only here (Windows repository).
    from Numerical.parallel_t3.wilson.portal_tensor_audit import _component_tensor
    from Numerical.parallel_t3.wilson.portal_orientation_v6 import _projected_portals
    for tag in sorted(EXPECTED_MODELS):
        d1,d2,df=DIMENSIONS[tag[0]]
        require(mat['models'][tag].get('dimensions')==[d1,d2,df] and
                rgb['models'][tag].get('dimensions')==[d1,d2,df],
                tag+': representation dimensions mismatch')
        v6s=v6['models'][tag].get('sectors')
        require(isinstance(v6s,list) and len(v6s)==2,tag+': v6 sectors not two')
        v8s=v8.get('model_sectors',{}).get(tag,[])
        rows=[]
        for which,d in ((1,d1),(2,d2)):
            m=sorted((p for p in mat['models'][tag]['portals'] if p.get('which')==which),
                     key=lambda p:p.get('invariant',-1))
            r=sorted((p for p in rgb['models'][tag]['portals'] if p.get('which')==which),
                     key=lambda p:0 if p.get('kind')=='singlet' else 1)
            raw=[_component_tensor(rec,which,d,k) for k,rec in enumerate(m,1)]
            phased=_projected_portals(m,which,d)
            v8row=next((v for v in v8s if v.get('which')==which),None)
            row=audit_sector(raw,phased,m,r,v6s[which-1],tag=tag,which=which,d=d,
                             v8_sector=v8row,seed=26100+10*sector_count,samples=samples)
            rows.append(row)
            count+=len(m);sector_count+=1
            anti+=sum(t.get('phase_inside_plus_hc')=='I' for t in m)
            triplets+=int(d==3)
            max_error=max(max_error,row['max_relative_residual'])
        models[tag]=rows
    require(count==58 and sector_count==32 and triplets==14 and anti==14,
            'incorrect 58/32/14/14 native portal census')
    return {
        'format':FORMAT,'status':STATUS,
        'models_checked':16,'portal_sectors_checked':sector_count,
        'native_Matchete_portal_coupling_slots':count,
        'triplet_sectors_spectral_normalization_verified':triplets,
        'anti_Hermitian_Matchete_phase_slots_checked':anti,
        'max_constructed_real_action_relative_residual':max_error,
        'models':models,
        'certification':{
            'native_Matchete_portal_CG_source_guarded':True,
            'native_RGBeta_registered_portal_6Matrix_and_permutations_source_guarded':True,
            'conditional_real_Lam_action_algebra_verified':True,
            'native_RGBeta_Lam_projector_numerically_executed':False,
            'independent_RGBeta_to_Matchete_S2_frame_certified':False,
            'complete_scalar_potential_dictionary_certified':False,
            'BCminus1_mixed_scalar_mass_diagonalization_certified':False,
            'signed_GammaCC_to_Weyl_certified':False,
            'physical_Wilson_evolution_authorized':False,
            'EFT_order':'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT',
        },
    }


def audit_from_paths(mat_path:Path,rgb_path:Path,v5_path:Path,
                     fierz_path:Path,v6_path:Path,v7_path:Path,
                     native_v8_path:Path,v8_path:Path,root:Path,
                     *,samples:int=5)->dict:
    """Re-run prior native source contracts and compare saved v8 row evidence."""
    from Numerical.parallel_t3.wilson.s2_spectral_v8 import (
        audit_from_paths as v8_recheck,STATUS_REDUCED,
    )
    require(1<=samples<=15,'invalid sample count')
    actual_v8=v8_recheck(mat_path,rgb_path,v5_path,fierz_path,v6_path,
                         v7_path,native_v8_path,root)
    saved_v8=read(v8_path)
    require(saved_v8.get('format')==actual_v8['format'] and
            saved_v8.get('status')==actual_v8['status']==STATUS_REDUCED,
            'saved v8 not the native-reduction-observed result')
    require(saved_v8.get('input_sha256')==actual_v8.get('input_sha256'),
            'saved v8 provenance hashes modified')
    require(saved_v8.get('triplet_sectors_basis_invariant_spectrum_verified')==14 and
            set(saved_v8.get('model_sectors',{}))==EXPECTED_MODELS,
            'saved v8 model/triplet inventory changed')
    for tag in EXPECTED_MODELS:
        old=saved_v8['model_sectors'][tag];new=actual_v8['model_sectors'][tag]
        require(len(old)==len(new),tag+': saved v8 sector count changed')
        for a,b in zip(old,new):
            require(a.get('which')==b['which'] and
                    a.get('sign_relative_to_first_native_traceless_CG')==
                        b['sign_relative_to_first_native_traceless_CG'] and
                    relative(a.get('selected_spectral_checks',{}).get('eigenvalues',[]),
                             b['selected_spectral_checks']['eigenvalues'])<TOL,
                    tag+': saved v8 spectral evidence changed')
    result=audit(read(mat_path),read(rgb_path),read(v6_path),actual_v8,samples=samples)
    result['input_sha256']={name:sha(path) for name,path in (
        ('matchete_native_v5',mat_path),('rgbeta_native_v5',rgb_path),
        ('portal_dictionary_v5',v5_path),('native_fierz_v6',fierz_path),
        ('portal_orientation_v6',v6_path),('s2_intertwiner_v7',v7_path),
        ('native_s2_reduction_v8',native_v8_path),('s2_spectrum_v8',v8_path),
    )}
    return result
