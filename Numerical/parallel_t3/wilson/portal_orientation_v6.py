"""Workstream A v6: native SU(2) doublet Fierz and conditional triplet Lie-sign audit.

Input: CURRENT v5 Matchete/RGBeta native portal exports, unchanged v5 audit,
and a live native RGBeta SU(2) Fierz witness. Does not modify any coupling.

Doublet: compare with RGBeta's exact fundamental completeness identity.
Triplet: impose Pauli-frame [G1,G2] = i G3, C2=2 and covariance. This is
NOT a native RGBeta S2 -> Matchete Cartesian intertwiner certificate.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import numpy as np

from .matchete_rgbeta_dictionary import DictionaryEvidenceError
from .portal_tensor_audit import (
    MAT_FORMAT, RGB_FORMAT, OUT_FORMAT as V5_FORMAT, STATUS as V5_STATUS,
    REP, _component_tensor, audit_from_paths as v5_audit_from_paths,
)

FORMAT = 't3-parallel-a-portal-orientation-fierz-lie-v6'
STATUS = ('DoubletNativeRGBetaFierzFixed_TripletPauliLieSignFixedConditional_'
          'NativeS2IntertwinerUnverified')
FIERZ_FORMAT = 't3-parallel-a-native-rgbeta-fierz-witness-v6'
FIERZ_STATUS = 'NativeSU2FundamentalFierzExact_S2GeneratorsStillAbstract'
SOURCES_FIERZ = ('RGE/running/rgbeta/T3RGBetaModel.wl',
                 'scripts/RunT3ParallelAExportNativeFierz.wl')
TOL = 2.e-10


def _require(ok: bool, why: str) -> None:
    if not ok:
        raise DictionaryEvidenceError('Portal v6: ' + why)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _read(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    _require(isinstance(obj, dict), f'{path}: JSON is not an object')
    return obj


def _rel(left: np.ndarray | float | complex,
         right: np.ndarray | float | complex) -> float:
    a = np.asarray(left)
    b = np.asarray(right)
    return float(np.linalg.norm(a-b)/max(1., np.linalg.norm(a), np.linalg.norm(b)))


def _pauli() -> np.ndarray:
    return np.asarray(([[0., .5], [.5, 0.]],
                       [[0., -.5j], [.5j, 0.]],
                       [[.5, 0.], [0., -.5]]), dtype=complex)


def _identity(d: int) -> np.ndarray:
    return np.einsum('ij,ab->ijab', np.eye(2), np.eye(d)).astype(complex)


def _fundamental_fierz() -> np.ndarray:
    """RGBeta SU(2) Fierz: 1/2(delta_il delta_kj - delta_ij delta_kl/2).

    Axes: {Hbar_i,H_j,Sbar_k,S_l}; both H and S transform as SU2L[fund].
    """
    return np.asarray([
        .5*(int(i == b and a == j) - .5*int(i == j and a == b))
        for i in range(2) for j in range(2)
        for a in range(2) for b in range(2)
    ], dtype=complex).reshape(2,2,2,2)


def _projected_portals(records: list[dict], which: int, d: int) -> list[np.ndarray]:
    out = []
    for index, r in enumerate(records, 1):
        raw = _component_tensor(r, which, d, index)
        phase = 1j if r['phase_inside_plus_hc'] == 'I' else 1.
        value = phase*raw
        out.append((value + value.conj().transpose(1,0,3,2))/2)
    return out


def _doublet_sector(projected: list[np.ndarray], v5: dict, model: str,
                    which: int, rng: np.random.Generator) -> dict:
    _require(len(projected) == 2, f'{model}/S{which}: wrong doublet portal count')
    ident = _identity(2)
    k = _fundamental_fierz()
    pauli = _pauli()
    direct = np.einsum('aij,akl->ijkl',pauli,pauli)
    _require(_rel(k,direct)<TOL, 'Fundamental Fierz differs from Pauli completeness')
    _require(_rel(np.vdot(ident,k),0) < TOL, 'Fierz tensor not traceless')
    _require(_rel(np.vdot(k,k).real,.75)<TOL, 'Fierz generator normalization changed')
    transform = np.zeros((2,2),dtype=float)
    worst = 0.
    for col,t in enumerate(projected):
        a = np.vdot(ident,t)/np.vdot(ident,ident)
        b = np.vdot(k,t)/np.vdot(k,k)
        _require(max(abs(a.imag),abs(b.imag))<TOL,
                 f'{model}/S{which}: non-real doublet coupling coefficients')
        transform[:,col] = [a.real,b.real]
        residual=_rel(t,a.real*ident+b.real*k)
        _require(residual<TOL,f'{model}/S{which}: Fierz does not span native portals')
        worst=max(worst,residual)
    _require(np.linalg.matrix_rank(transform,tol=1.e-10)==2,
             f'{model}/S{which}: singular Fierz mapping')
    old = np.asarray(v5['candidate_rgbeta_from_native_matchete_real_matrix'],float)
    _require(old.shape==(2,2) and _rel(transform[0],old[0])<TOL,
             f'{model}/S{which}: previous singlet map changed')
    plus=_rel(transform[1],old[1]);minus=_rel(transform[1],-old[1])
    _require(min(plus,minus)<TOL,
             f'{model}/S{which}: v5 oriented adjoint differs from native Fierz')
    sign=1 if plus<minus else -1
    for _ in range(9):
        h=rng.normal(size=2)+1j*rng.normal(size=2)
        s=rng.normal(size=2)+1j*rng.normal(size=2)
        l=rng.normal(size=2)
        raw=sum(c*np.einsum('ijab,i,j,a,b->',t,h.conj(),h,s.conj(),s)
                for c,t in zip(l,projected))
        rhs=transform@l
        expectation=rhs[0]*np.vdot(h,h)*np.vdot(s,s)+rhs[1]*sum(
            (h.conj() @ g @ h)*(s.conj() @ g @ s) for g in pauli)
        worst=max(worst,_rel(raw,expectation))
        _require(_rel(raw,expectation)<TOL,
                 f'{model}/S{which}: doublet potential holdout failed')
    return {
        'sector':'doublet', 'model':model, 'which':which,
        'native_generator_orientation':'SU2L[fund]_Fierz_same_component_axes',
        'matchete_source_couplings':v5['native_matchete_coupling_order'],
        'rgbeta_coupling_order':[f'lambdaH{which}',f'lambdaH{which}Adj'],
        'rgbeta_from_matchete_real_matrix':transform.tolist(),
        'orientation_relative_to_v5_maybe_oriented_adjoint':sign,
        'max_relative_component_and_potential_residual':worst,
        'native_RGBeta_fundamental_fierz_verified':True,
        'relative_doublet_component_frame_assumption':
            'Matchete H and S doublets use common native SU2L[fund] indices',
        'native_RGBeta_complex_scalar_coordinate_bridge_certified':False,
        'physical_complete_potential_map_certified':False,
    }


def _extract_triplet_generators(k: np.ndarray) -> np.ndarray:
    # Since Tr(T_H^a T_H^b)=delta_ab/2, multiplying by 2 gives G_R^a.
    # Matrix index orientation is crucial: trace contracts T_H^a[j,i] K[i,j,a,b].
    return np.asarray([2*np.einsum('ji,ijab->ab',g,k) for g in _pauli()])


def _triplet_sign_checks(k: np.ndarray) -> dict:
    """Orientation check in an EXPLICIT Pauli frame, not native RGBeta S2.

    A shared adjoint index A is used for both H and R. Reversing sign on only
    the R generator violates [G1,G2]=iG3, while leaving Casimir unchanged.
    """
    g=_extract_triplet_generators(k)
    ident=np.eye(3,dtype=complex)
    herm=max(_rel(x,x.conj().T) for x in g)
    traces=max(abs(np.trace(x)) for x in g)
    gram=np.asarray([[np.trace(x@y) for y in g] for x in g])
    index=_rel(gram,2*np.eye(3))
    cas=_rel(sum(x@x for x in g),2*ident)
    comm=max(_rel(g[a]@g[b]-g[b]@g[a],1j*sum(
        (1 if (a,b,c) in ((0,1,2),(1,2,0),(2,0,1)) else
         -1 if (a,b,c) in ((1,0,2),(2,1,0),(0,2,1)) else 0)*g[c]
        for c in range(3))) for a,b in ((0,1),(1,2),(2,0)))
    k_matrix=k.transpose(0,2,1,3).reshape(6,6)
    covariance=0.0
    for a in range(3):
        total=np.kron(_pauli()[a],ident)+np.kron(np.eye(2),g[a])
        covariance=max(covariance,_rel(total@k_matrix,k_matrix@total))
    return {'generators':g, 'hermiticity_residual':herm,
            'trace_residual':float(traces), 'trace_normalization_residual':index,
            'quadratic_casimir_residual':cas,'lie_bracket_residual':comm,
            'simultaneous_SU2_covariance_residual':covariance,
            'max_residual':max(herm,traces,index,cas,comm,covariance)}


def _triplet_sector(projected: list[np.ndarray],v5: dict, model: str,
                    which: int,rng:np.random.Generator) -> dict:
    _require(len(projected)==2, f'{model}/S{which}: wrong triplet count')
    old=np.asarray(v5['candidate_rgbeta_from_native_matchete_real_matrix'],float)
    _require(old.shape==(2,2),f'{model}/S{which}: v5 map incorrect dimensions')
    ident=_identity(3)
    ortho=[]
    singlet=[]
    for t in projected:
        x=np.vdot(ident,t)/np.vdot(ident,ident)
        _require(abs(x.imag)<TOL,f'{model}/S{which}: complex triplet singlet')
        singlet.append(float(x.real))
        ortho.append(t-x.real*ident)
    pivot=next((x for x in ortho if np.linalg.norm(x)>1.e-10),None)
    _require(pivot is not None, f'{model}/S{which}: absent adjoint')
    # SU2L[S2] has trace index T(S2)=2: ||K_HR||^2=3*(1/2)*2=3.
    candidate=(3.**.5)*pivot/np.linalg.norm(pivot)
    yes=_triplet_sign_checks(candidate)
    no=_triplet_sign_checks(-candidate)
    _require(min(yes['max_residual'],no['max_residual'])<TOL,
             f'{model}/S{which}: neither triplet orientation closes SU2')
    _require(max(yes['lie_bracket_residual'],no['lie_bracket_residual'])>1.e-4,
             f'{model}/S{which}: ambiguous triplet SU2 bracket orientation')
    sign=1 if yes['max_residual']<no['max_residual'] else -1
    selected=candidate*sign
    check=yes if sign==1 else no
    _require(check['max_residual']<TOL,
             f'{model}/S{which}: selected triplet algebra/covariance failed')
    target=np.zeros((2,2),dtype=float)
    worst=0.0
    for col,t in enumerate(projected):
        a=np.vdot(ident,t)/np.vdot(ident,ident)
        b=np.vdot(selected,t)/np.vdot(selected,selected)
        _require(abs(a.imag)<TOL and abs(b.imag)<TOL,
                 f'{model}/S{which}: complex triplet projection')
        target[:,col]=[a.real,b.real]
        worst=max(worst,_rel(t,target[0,col]*ident+target[1,col]*selected))
    _require(_rel(target[0],old[0])<TOL and
             min(_rel(target[1],old[1]),_rel(target[1],-old[1]))<TOL,
             f'{model}/S{which}: v5 triplet map changed')
    gen=_extract_triplet_generators(selected)
    pauli=_pauli()
    for _ in range(9):
        h=rng.normal(size=2)+1j*rng.normal(size=2)
        s=rng.normal(size=3)+1j*rng.normal(size=3)
        l=rng.normal(size=2)
        action=sum(c*np.einsum('ijab,i,j,a,b->',t,h.conj(),h,s.conj(),s)
                   for c,t in zip(l,projected))
        mapped=target@l
        expectation=mapped[0]*np.vdot(h,h)*np.vdot(s,s)+mapped[1]*sum(
            (h.conj() @ x @ h)*(s.conj() @ y @ s)
            for x,y in zip(pauli,gen))
        worst=max(worst,_rel(action,expectation))
        _require(_rel(action,expectation)<TOL,
                 f'{model}/S{which}: triplet potential holdout failed')
    return {
        'sector':'triplet','model':model,'which':which,
        'matchete_source_couplings':v5['native_matchete_coupling_order'],
        'candidate_rgbeta_coupling_order':[f'lambdaH{which}',
                                          f'lambdaH{which}Adj_PauliFrameConditional'],
        'candidate_rgbeta_from_matchete_real_matrix':target.tolist(),
        'selected_Pauli_frame_sign_relative_to_v5_adjoint':sign,
        'selected_generator_checks':{k:v for k,v in check.items() if k!='generators'},
        'wrong_sign_lie_bracket_residual':float((no if sign==1 else yes)['lie_bracket_residual']),
        'max_relative_component_and_potential_residual':worst,
        'component_frame':'H Pauli/2; R reconstructed in registered Matchete Cartesian CG axes',
        'native_RGBeta_S2_intertwiner_certified':False,
        'native_RGBeta_symbolic_tGen_component_matrix_extracted':False,
        'physical_complete_potential_map_certified':False,
    }


def _check_fierz_witness(root:Path, witness:dict) -> dict:
    _require(witness.get('format')==FIERZ_FORMAT and
             witness.get('status')==FIERZ_STATUS,
             'Invalid native RGBeta Fierz witness format/status')
    flags=('native_fundamental_Fierz_exact','native_symmetric_S2_rep_dimension_three',
           'native_S2_trace_normalization_two','native_triplet_generators_symbolic',
           'native_barred_generator_minus_transpose')
    _require(all(witness.get(f) is True for f in flags),
             'Incomplete native Fierz/Casimir/conjugation check')
    _require(witness.get('native_S2_component_intertwiner_certified') is False,
             'Unsupported native S2 certification claim')
    hashes=witness.get('source_sha256')
    _require(isinstance(hashes,dict) and set(hashes)==set(SOURCES_FIERZ),
             'Missing native Fierz source hashes')
    for name in SOURCES_FIERZ:
        _require((root/name).is_file() and hashes[name]==_sha(root/name),
                 f'Native RGBeta Fierz witness stale against {name}')
    return {'format':witness['format'],'status':witness['status'],
            'source_sha256':hashes}


def audit(mat:dict,rgb:dict,v5:dict,witness:dict,root:Path,
          *,input_sha256:dict[str,str]) -> dict:
    _check_fierz_witness(root,witness)
    _require(mat.get('format')==MAT_FORMAT and rgb.get('format')==RGB_FORMAT,
             'Wrong native v5 exporter schema')
    _require(v5.get('format')==V5_FORMAT and v5.get('status')==V5_STATUS,
             'Wrong v5 result schema/status')
    _require(v5.get('input_sha256')=={
        'matchete_portals':input_sha256['matchete'],
        'rgbeta_portals':input_sha256['rgbeta']},
        'v5 result not derived from THESE exact native portal exports')
    _require(v5.get('portal_coupling_slots')==58 and
             v5.get('native_anti_Hermitian_CGs_given_I_inside_PlusHc')==14,
             'Missing v5 58 slots or 14 anti-Hermitian corrections')
    rng=np.random.default_rng(66117)
    models={};doublet=triplet=singlet=0;worst=0.0
    for tag in sorted(v5['models']):
        cls=tag[0];d1,d2,df,alphas=REP[cls]
        _require(tag in mat['models'] and tag in rgb['models'] and
                 int(tag.split('_')[1]) in alphas,
                 f'{tag}: mismatched native source model')
        sectors=[]
        for which,d in ((1,d1),(2,d2)):
            v5sector=v5['models'][tag]['sectors'][which-1]
            source=sorted([r for r in mat['models'][tag]['portals']
                          if r.get('which')==which],key=lambda r:r['invariant'])
            projected=_projected_portals(source,which,d)
            if d==2:
                row=_doublet_sector(projected,v5sector,tag,which,rng)
                doublet+=1
            elif d==3:
                row=_triplet_sector(projected,v5sector,tag,which,rng)
                triplet+=1
            else:
                _require(len(projected)==1 and
                         _rel(projected[0],_identity(1))<TOL,
                         f'{tag}/S{which}: singlet identity changed')
                row={'sector':'singlet','model':tag,'which':which,
                     'rgbeta_from_matchete_real_matrix':[[1.]],
                     'orientation_sign':'not_applicable'}
                singlet+=1
            worst=max(worst,row.get('max_relative_component_and_potential_residual',0.))
            sectors.append(row)
        models[tag]={'dimensions':[d1,d2,df],'sectors':sectors,
                     'matching_mass_basis_ready':v5['models'][tag][
                         'matchete_native_build_validation']['matching_mass_basis_ready']}
    _require(set(models)==set(v5['models']) and (doublet,triplet,singlet)==(12,14,6),
             'Unexpected all-16 SU2 portal representation census')
    _require(worst<TOL,'Large portal orientation reconstruction error')
    return {'format':FORMAT,'status':STATUS,
            'input_sha256':input_sha256,
            'native_fierz_witness':witness['status'],
            'doublet_sectors_Fierz_verified':doublet,
            'triplet_sectors_Pauli_algebra_conditional':triplet,
            'singlet_sectors':singlet,
            'original_matchete_portal_slots':58,
            'max_relative_component_and_potential_residual':worst,
            'models':models,
            'convention':{
                'doublet_Fierz':'sum_A tFund[A,i,j] tFund[A,k,l] = '
                                  '(delta[i,l]delta[k,j] - delta[i,j]delta[k,l]/2)/2',
                'Pauli_frame':'tH=(sigma1,sigma2,sigma3)/2, '
                              '[tR1,tR2]=i tR3, Tr[tR^a tR^b]=T(R)delta_ab',
                'triplet':'Matchete registered Cartesian CG axes; native RGBeta '
                          'S2 basis is abstract and has NOT been rotated into these axes'},
            'scope':{
                'native_RGBeta_doublet_Fierz_identity_verified':True,
                'doublet_relative_adjoint_sign_fixed_in_common_fundamental_axes':True,
                'triplet_sign_fixed_only_in_explicit_Pauli_Matchete_frame':True,
                'native_RGBeta_triplet_intertwiner_certified':False,
                'physical_complete_scalar_coupling_dictionary_certified':False,
                'BCminus1_mixed_scalar_mass_diagonalization_done':False,
                'fermion_GammaCC_to_Weyl_normalization_certified':False,
                'physical_Wilson_background_certified':False,
            }}


def audit_from_paths(mat_path:Path,rgb_path:Path,v5_path:Path,
                     fierz_path:Path,root:Path)->dict:
    mat=_read(mat_path);rgb=_read(rgb_path);v5=_read(v5_path)
    witness=_read(fierz_path)
    evidence={
        'matchete':_sha(mat_path),'rgbeta':_sha(rgb_path),
        'v5_audit':_sha(v5_path),'native_fierz':_sha(fierz_path)}
    # Revalidate all source hashes, built-model statuses and the v5 tensor
    # identities independently; no trust in a stale/edited v5 JSON.
    fresh=audit_from_paths_v5(mat_path,rgb_path,root)
    _require(_rel(fresh['max_relative_tensor_and_potential_residual'],
                  v5.get('max_relative_tensor_and_potential_residual',-1))<TOL,
             'Previous v5 audit residual is inconsistent')
    _require(fresh['models_requiring_mixed_scalar_mass_diagonalization']==
             v5.get('models_requiring_mixed_scalar_mass_diagonalization'),
             'Previous v5 mixed-mass inventory is inconsistent')
    for model,entry in fresh['models'].items():
        _require(model in v5.get('models',{}),'Previous v5 model missing')
        for a,b in zip(entry['sectors'],v5['models'][model]['sectors']):
            _require(_rel(np.asarray(a['candidate_rgbeta_from_native_matchete_real_matrix']),
                          np.asarray(b['candidate_rgbeta_from_native_matchete_real_matrix']))<TOL,
                     f'Previous v5 coefficients edited for {model}')
    return audit(mat,rgb,v5,witness,root,input_sha256=evidence)


def audit_from_paths_v5(mat_path:Path,rgb_path:Path,root:Path)->dict:
    return v5_audit_from_paths(mat_path,rgb_path,root)
