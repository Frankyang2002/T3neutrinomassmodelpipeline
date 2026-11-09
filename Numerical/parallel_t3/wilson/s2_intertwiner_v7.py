"""Workstream A v7: S2 representation consistency, NOT native RGBeta frame certification.

The SU2 symmetric-square generators are induced from Pauli/2 in an independently
DECLARED weight basis.  The comparison intertwiner is solved using generators
reconstructed from Matchete portal CGs.  It is NOT extracted from RGBeta.

Inputs are the unchanged v5/v6 exports and optional archived provenance witnesses.
No physical couplings, Wilson coefficients, source models, or RG evolution are changed.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import numpy as np

FORMAT = 't3-parallel-a-native-s2-frame-intertwiner-audit-v7'
STATUS = 'CanonicalS2RepresentationAndMatcheteComparisonVerified_NativeRGBetaS2FrameUnspecifiedConditional'
TOL = 2.0e-10
EXPECTED_TRIPLETS = 14


class S2EvidenceError(ValueError):
    pass


def require(ok: bool, message: str) -> None:
    if not ok:
        raise S2EvidenceError('Workstream A v7: ' + message)


def read_json(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding='utf-8-sig'))
    except (OSError, ValueError) as exc:
        raise S2EvidenceError(f'Cannot read {path}: {exc}') from exc
    require(isinstance(data, dict), f'{path}: object expected')
    return data


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def rel(a: Any, b: Any) -> float:
    a, b = np.asarray(a), np.asarray(b)
    return float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(a),np.linalg.norm(b)))


def pauli_half() -> np.ndarray:
    return np.asarray(([[0., .5], [.5, 0.]],
                       [[0., -.5j], [.5j, 0.]],
                       [[.5, 0.], [0., -.5]]),dtype=np.complex128)


def symmetric_square_cg() -> np.ndarray:
    """A CHOSEN basis, not a measurement of RGBeta's delS2 coordinates."""
    c=np.zeros((3,2,2),dtype=np.complex128)
    c[0,0,0]=1.
    c[1,0,1]=c[1,1,0]=1/np.sqrt(2)
    c[2,1,1]=1.
    return c


def induced_s2_generators(c: np.ndarray) -> np.ndarray:
    require(c.shape==(3,2,2), 'S2 CG shape must be (3,2,2)')
    return np.asarray([
        np.einsum('uij,ik,vkj->uv',c.conj(),t,c) +
        np.einsum('uij,jk,vik->uv',c.conj(),t,c)
        for t in pauli_half()], dtype=np.complex128)


def generator_checks(g: np.ndarray) -> dict[str,float]:
    require(g.shape==(3,3,3), 'triplet generator shape must be (3,3,3)')
    gram=np.asarray([[np.trace(a@b) for b in g] for a in g])
    bracket=max(rel(g[a]@g[b]-g[b]@g[a],1j*g[c])
                for a,b,c in ((0,1,2),(1,2,0),(2,0,1)))
    out={
        'hermiticity':max(rel(a,a.conj().T) for a in g),
        'tracelessness':max(float(abs(np.trace(a))) for a in g),
        'trace_normalization':rel(gram,2*np.eye(3)),
        'casimir':rel(sum(a@a for a in g),2*np.eye(3)),
        'lie_bracket':bracket,
    }
    out['maximum']=max(out.values())
    return out


def canonical_checks() -> tuple[np.ndarray,np.ndarray,dict]:
    c=symmetric_square_cg()
    g=induced_s2_generators(c)
    sym=np.zeros((2,2,2,2),dtype=float)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    sym[i,j,k,l]=(int(i==k and j==l)+int(i==l and j==k))/2
    gram=rel(np.einsum('aij,bij->ab',c.conj(),c),np.eye(3))
    projector=rel(np.einsum('aij,akl->ijkl',c.conj(),c),sym)
    symmetry=rel(c,c.transpose(0,2,1))
    equiv=max(rel(np.einsum('vu,vij->uij',g[a],c),
                  np.einsum('ik,ukj->uij',pauli_half()[a],c)+
                  np.einsum('jk,uik->uij',pauli_half()[a],c))
              for a in range(3))
    # The equality uses the transpose on the triplet CG row index.  This is
    # not a claim that RGBeta uses this particular numerical CG frame.
    info={'S2_isometry':gram,'symmetric_projector':projector,
          'CG_exchange_symmetry':symmetry,'fundamental_induced_action':equiv,
          'generators':generator_checks(g)}
    info['maximum']=max(gram,projector,symmetry,equiv,info['generators']['maximum'])
    require(info['maximum']<TOL,'Declared S2 frame fails its exact identities')
    return c,g,info


def adjoint_tensor(g: np.ndarray) -> np.ndarray:
    return np.einsum('aij,auv->ijuv',pauli_half(),g)


def intertwiner_from_generators(g_target:np.ndarray,g_source:np.ndarray)->tuple[np.ndarray,dict]:
    """Solve G_target^a U = U G_source^a, without using portal coefficients.

    Important: if g_target itself was derived from a Matchete portal CG,
    the result is a Matchete comparison, not an independent RGBeta witness.
    """
    require(g_target.shape==(3,3,3) and g_source.shape==(3,3,3),
            'Generator dimensions incorrect')
    system=np.concatenate([
        np.kron(np.eye(3),g_target[a])-
        np.kron(g_source[a].T,np.eye(3)) for a in range(3)],axis=0)
    _,s,vh=np.linalg.svd(system,full_matrices=False)
    require(s[-1]<TOL and s[-2]>1.e-5,
            'SU(2) intertwiner nullspace is not one dimensional')
    u=vh[-1].conj().reshape(3,3,order='F')
    u*=np.sqrt(3.)/np.linalg.norm(u)
    # Global phase is arbitrary; select a deterministic *reporting* convention.
    pivot=np.unravel_index(np.abs(u).argmax(),u.shape)
    u*=np.exp(-1j*np.angle(u[pivot]))
    checks={'intertwining':max(rel(g_target[a]@u,u@g_source[a]) for a in range(3)),
            'unitarity':rel(u.conj().T@u,np.eye(3)),
            'nullspace_gap':float(s[-2]),
            'smallest_singular_value':float(s[-1])}
    require(max(checks['intertwining'],checks['unitarity'])<TOL,
            'Computed comparison intertwiner not unitary/covariant')
    return u,checks


def frame_ambiguity(g_c:np.ndarray,c:np.ndarray,u:np.ndarray,g_m:np.ndarray)->dict:
    """Two distinguishable frames satisfy all native abstract S2 normalizations."""
    phases=np.exp(1j*np.asarray([0.37,-0.61,1.07]))
    w=np.diag(phases)
    g_rot=np.einsum('ux,axy,vy->auv',w,g_c,w.conj())
    c_rot=np.einsum('ux,xij->uij',w,c)
    u_rot=u@w.conj().T
    checks={
        'rotated_cg_isometry':rel(np.einsum('aij,bij->ab',c_rot.conj(),c_rot),np.eye(3)),
        'rotated_generator_algebra':generator_checks(g_rot)['maximum'],
        'rotated_intertwining':max(rel(g_m[a]@u_rot,u_rot@g_rot[a]) for a in range(3)),
        'distinct_generator_frame_distance':rel(g_c,g_rot),
        'transported_adjoint_tensor_invariance':rel(
            np.einsum('ux,ijxy,vy->ijuv',u_rot,adjoint_tensor(g_rot),u_rot.conj()),
            np.einsum('ux,ijxy,vy->ijuv',u,adjoint_tensor(g_c),u.conj()))}
    require(max(v for k,v in checks.items() if k!='distinct_generator_frame_distance')<TOL,
            'Unitary change of S2 frame did not preserve abstract identities')
    require(checks['distinct_generator_frame_distance']>0.01,
            'Nontrivial S2 alternative frame unexpectedly identical')
    return checks


def _complex_matrix(m:np.ndarray)->list[list[list[float]]]:
    return [[[float(v.real),float(v.imag)] for v in row] for row in m]


def conditional_sector(c:np.ndarray,g_c:np.ndarray,mat_portals:list[dict],
                       v6_sector:dict,which:int,model:str)->dict:
    """Compare with Matchete portal-derived generators; never call them RGBeta-native."""
    from Numerical.parallel_t3.wilson.portal_orientation_v6 import (
        _identity,_projected_portals,_extract_triplet_generators)
    mat_portals=sorted((r for r in mat_portals if r.get('which')==which),
                       key=lambda r:r['invariant'])
    projected=_projected_portals(mat_portals,which,3)
    # v6 triplet sectors use this key. The similar '...native_matchete...'
    # name belongs to the v5 sector schema and must not be used here.
    require(v6_sector.get('sector')=='triplet',
            f'{model}/S{which}: expected a v6 triplet sector')
    require('candidate_rgbeta_from_matchete_real_matrix' in v6_sector,
            f'{model}/S{which}: v6 triplet coupling matrix missing')
    a=np.asarray(v6_sector['candidate_rgbeta_from_matchete_real_matrix'],float)
    require(len(projected)==2 and a.shape==(2,2),f'{model}/S{which}: triplet columns missing')
    ident=_identity(3)
    adj=[]
    for j in range(2):
        if abs(a[1,j])>1.e-8:
            adj.append((projected[j]-a[0,j]*ident)/a[1,j])
    require(len(adj)>0,f'{model}/S{which}: no adjoint direction')
    k=adj[0]
    span=max([rel(k,z) for z in adj]+[
        rel(projected[j],a[0,j]*ident+a[1,j]*k) for j in range(2)])
    require(span<TOL,f'{model}/S{which}: v6 Matchete portal reconstruction changed')
    g_m=_extract_triplet_generators(k)
    g_check=generator_checks(g_m)
    require(g_check['maximum']<TOL,f'{model}/S{which}: portal generators do not close SU2')
    require(rel(k,adjoint_tensor(g_m))<TOL,
            f'{model}/S{which}: generator reconstruction not inverse to K')
    u,uc=intertwiner_from_generators(g_m,g_c)
    transported=np.einsum('ux,ijxy,vy->ijuv',u,adjoint_tensor(g_c),u.conj())
    component=rel(k,transported)
    require(component<TOL,f'{model}/S{which}: SU2 transport fails')
    alternate=frame_ambiguity(g_c,c,u,g_m)
    return {'model':model,'which':which,
            'comparison_basis':'Matchete-native registered portal CG axes',
            'source_generator_origin':'RECONSTRUCTED_FROM_MATCHEte_PORTAL_CG_NOT_INDEPENDENT',
            'canonical_S2_basis_origin':'DECLARED_HIGHEST_WEIGHT_NOT_NATIVE_RGBeta',
            'intertwiner_can_to_matchete_real_imag':_complex_matrix(u),
            'intertwiner_checks':uc,
            'matchete_portal_generator_checks':g_check,
            'matchete_generator_from_portal_T3_real_imag':_complex_matrix(g_m[2]),
            'adjoint_tensor_transport_residual':component,
            'v6_component_span_residual':span,
            'alternative_frame_checks':alternate,
            'independent_native_RGBeta_to_Matchete_intertwiner_certified':False}


def audit(mat:dict,v6:dict, native_rgbeta:dict|None=None,
          generator_candidate:dict|None=None)->dict:
    from Numerical.parallel_t3.wilson.portal_orientation_v6 import STATUS as V6_STATUS
    require(v6.get('status')==V6_STATUS and
            v6.get('triplet_sectors_Pauli_algebra_conditional')==14,
            'v6 all-16 conditional evidence missing')
    require(v6.get('scope',{}).get('native_RGBeta_triplet_intertwiner_certified') is False,
            'v6 native S2 claim contradicts original evidence')
    c,g,canonical=canonical_checks()
    if native_rgbeta is not None:
        require(native_rgbeta.get('status')==
                'NATIVE_SYMBOLIC_GENERATOR_RULES_CAPTURED_NUMERIC_AXES_PENDING',
                'Native RGBeta witness status changed: requires NEW verification')
        for k in ('RGBeta_numeric_eps_components_certified',
                  'RGBeta_numeric_delS2_components_certified',
                  'RGBeta_to_Matchete_basis_identified',
                  'physical_coupling_dictionary_certified'):
            require(native_rgbeta.get(k) is False,
                    f'Native RGBeta witness {k} unexpectedly changed')
    if generator_candidate is not None:
        require(generator_candidate.get('generator_basis_is_live_export') is False and
                generator_candidate.get('status')==
                'candidate_from_repository_test_fixtures_not_live_export',
                'Unknown generator candidate provenance: cannot use as native')
    models={};num=0;worst=0.
    for tag in sorted(v6['models']):
        require(tag in mat.get('models',{}),f'{tag}: missing native Matchete CG inventory')
        rows=[]
        for sector in v6['models'][tag]['sectors']:
            if sector['sector']!='triplet':continue
            which=sector['which']
            row=conditional_sector(c,g,mat['models'][tag]['portals'],sector,which,tag)
            rows.append(row);num+=1
            worst=max(worst,row['adjoint_tensor_transport_residual'],
                      row['v6_component_span_residual'],
                      row['matchete_portal_generator_checks']['maximum'],
                      row['intertwiner_checks']['intertwining'],
                      row['intertwiner_checks']['unitarity'])
        models[tag]=rows
    require(num==EXPECTED_TRIPLETS,'Triplet census is not 14')
    return {
        'format':FORMAT,'status':STATUS,
        'triplet_sectors_conditional_comparison_verified':num,
        'canonical_symmetric_square':canonical,
        'maximum_comparison_residual':worst,
        'native_RGBeta_current_S2_numeric_component_frame_exposed':False,
        'independent_native_RGBeta_S2_to_Matchete_intertwiner_certified':False,
        'matchete_generator_fixture_is_live_native':False,
        'native_witness_status':native_rgbeta.get('status') if native_rgbeta else 'not_provided',
        'native_witness_is_historical_diagnostic_only':native_rgbeta is not None,
        'matchete_candidate_status':generator_candidate.get('status') if generator_candidate else 'not_provided',
        'models':models,
        'limitations':{
            'comparison_intertwiner_used_Matchete_portal_CG':True,
            'comparison_does_not_measure_RGBeta_triplet_frame':True,
            'native_RGBeta_S2_frame_is_abstract':True,
            'unitary_S2_basis_ambiguity_demonstrated':True,
            'complete_scalar_potential_normalization_certified':False,
            'BCminus1_mixed_scalar_mass_diagonalization_certified':False,
            'signed_GammaCC_to_Weyl_certified':False,
            'physical_Wilson_evolution_authorized':False,
            'physical_EFT_order':'UV -> all fermions -> all scalars -> SMEFT',
        },
    }


def audit_from_paths(mat_path:Path,rgb_path:Path,v5_path:Path,fierz_path:Path,
                     v6_path:Path,root:Path,
                     native_rgbeta_path:Path|None=None,
                     generator_candidate_path:Path|None=None)->dict:
    from Numerical.parallel_t3.wilson.portal_orientation_v6 import (
        audit_from_paths as recheck_v6)
    original=read_json(v6_path)
    fresh=recheck_v6(mat_path,rgb_path,v5_path,fierz_path,root)
    require(original.get('format')==fresh['format'] and
            original.get('status')==fresh['status'] and
            original.get('input_sha256')==fresh['input_sha256'],
            'Saved v6 does not match freshly revalidated v5/v6 source evidence')
    require(set(original.get('models',{}))==set(fresh['models']),
            'Saved v6 model inventory changed')
    for model in fresh['models']:
        saved_sectors=original['models'][model]['sectors']
        verified_sectors=fresh['models'][model]['sectors']
        require(len(saved_sectors)==len(verified_sectors),
                f'{model}: saved v6 sector count changed')
        for old, new in zip(saved_sectors,verified_sectors):
            require(old.get('sector')==new['sector'] and
                    old.get('which')==new['which'],
                    f'{model}: saved v6 sector identity changed')
            # The triplet matrix's canonical field name in the *v6* result
            # drops "native"; v5's similar name is NOT a v6 field.
            key=('candidate_rgbeta_from_matchete_real_matrix'
                 if new['sector']=='triplet' else
                 'rgbeta_from_matchete_real_matrix')
            require(key in old and key in new,
                    f'{model}: saved v6 coupling map field missing ({key})')
            require(rel(np.asarray(old[key]),np.asarray(new[key]))<TOL,
                    f'{model}: saved v6 coupling map differs from native recheck')
    other=read_json(native_rgbeta_path) if native_rgbeta_path else None
    candidate=read_json(generator_candidate_path) if generator_candidate_path else None
    report=audit(read_json(mat_path),fresh,other,candidate)
    report['input_sha256']={
        'matchete_v5':digest(mat_path),'rgbeta_v5':digest(rgb_path),
        'portal_dictionary_v5':digest(v5_path),
        'native_fierz_v6':digest(fierz_path),'portal_orientation_v6':digest(v6_path)}
    if native_rgbeta_path:report['input_sha256']['historical_native_symbolic_witness']=digest(native_rgbeta_path)
    if generator_candidate_path:report['input_sha256']['fixture_generator_candidate']=digest(generator_candidate_path)
    return report
