"""Workstream A v5: corrected native Matchete Higgs-portal action and RGBeta metadata.

The real-coupling map is expressed in a *Matchete-oriented* adjoint frame.
Its overall adjoint sign relative to RGBeta's symbolic tGen is NOT certified.
Only data registered by native Wolfram exporters are accepted. This module
never changes a source Lagrangian, coupling, or Wilson evolution result.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import re

import numpy as np

from .matchete_rgbeta_dictionary import (
    ALL16, DictionaryEvidenceError, matchete_number,
)

MAT_FORMAT = 't3-parallel-a-native-matchete-portal-cg-phases-v5'
RGB_FORMAT = 't3-parallel-a-native-rgbeta-portal-registrations-v5'
OUT_FORMAT = 't3-parallel-a-native-portal-singlet-adjoint-conditional-v5'
STATUS = 'NativeMatchetePortalPhasesAndFullRank_Verified_RGBetaAdjointOrientationUncertified'
REP = {
    'A': (1, 3, 2, (-4, -2, 0)),
    'B': (2, 2, 1, (-3, -1, 1)),
    'C': (2, 2, 3, (-3, -1, 1)),
    'D': (3, 1, 2, (-2, 0, 2)),
    'E': (3, 3, 2, (-4, -2, 0, 2)),
}
SOURCE_MAT = ('Lagrangian/model/T3ModelCatalog.wl',
              'Lagrangian/model/LagrangianBuilder.wl',
              'Lagrangian/model/SU2Invariants.wl',
              'Lagrangian/model/T3Fields.wl',
              'Lagrangian/interactions/T3Topology.wl',
              'Lagrangian/interactions/ScalarPotential.wl',
              'Lagrangian/interactions/GeneralScalarOperators.wl')
SOURCE_RGB = ('RGE/running/rgbeta/T3RGBetaModel.wl',
              'scripts/RGBetaNativeAccess.wl')


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise DictionaryEvidenceError('Portal v5: ' + message)


def _clean(s: str) -> str:
    _require(isinstance(s, str), 'Missing RGBeta InputForm')
    return re.sub(r'\s+', '', s)


def _sha_files(root: Path, evidence: dict, files: tuple[str, ...]) -> dict:
    hashes = evidence.get('source_sha256')
    _require(isinstance(hashes, dict) and set(hashes) == set(files),
             'Incomplete native source-hash inventory')
    out = {}
    for path in files:
        src = root / path
        _require(src.is_file(), f'Missing live source {path}')
        h = sha256(src.read_bytes()).hexdigest()
        _require(hashes[path] == h, f'Stale source evidence: {path}')
        out[path] = h
    return out


def _component_tensor(rec: dict, which: int, dim: int, number: int) -> np.ndarray:
    _require(rec.get('which') == which and rec.get('invariant') == number,
             'Matchete portal order or numbering changed')
    _require(rec.get('coupling') == f'lambdaH{which}Inv{number}',
             'Matchete portal coupling does not match native operator')
    _require(rec.get('cg') == f'T3HiggsPortal{which}CGInv{number}',
             'Unexpected native Matchete CG symbol')
    _require(rec.get('native_tensor_accessor') == 'Matchete`GetCGTensor' and
             rec.get('assembled_lagrangian_contains_cg_and_coupling') is True,
             'Missing native registered CG/operator evidence')
    _require(rec.get('dimensions') == [2, 2, dim, dim],
             'Matchete four-index portal dimensions changed')
    _require(rec.get('native_dimensions') == ([2, 2, dim, dim] if dim > 1 else [2, 2]),
             'Matchete native singleton axis contract changed')
    entries = rec.get('components_inputform')
    _require(isinstance(entries, list) and len(entries) == 4*dim*dim,
             'Incomplete registered native CG components')
    tensor = np.array([matchete_number(v) for v in entries],
                      dtype=np.complex128).reshape(2, 2, dim, dim)
    _require(float(np.linalg.norm(tensor)) > 1e-12, 'Zero native portal tensor')
    hc = tensor.conj().transpose(1, 0, 3, 2)
    scale = max(1., float(np.linalg.norm(tensor)))
    is_herm = np.linalg.norm(tensor-hc) < 1e-10*scale
    is_anti = np.linalg.norm(tensor+hc) < 1e-10*scale
    _require(is_herm != is_anti, 'Native CG is neither pure Hermitian nor anti-Hermitian')
    expected_phase = '1' if is_herm else 'I'
    _require(rec.get('phase_inside_plus_hc') == expected_phase,
             f'Incorrect PlusHc phase for H{which}Inv{number}: need {expected_phase}')
    _require(rec.get('native_hermiticity') ==
             ('Hermitian' if is_herm else 'AntiHermitian'),
             'Native Hermiticity metadata disagrees with components')
    return tensor


def _native_rgbeta(rec: dict, which: int, dim: int, kind: str) -> dict:
    symbol = f'lambdaH{which}' + ('Adj' if kind == 'adjoint' else '')
    _require(rec.get('symbol') == symbol and rec.get('which') == which and
             rec.get('kind') == kind, 'Unexpected RGBeta portal coupling')
    _require(rec.get('dimensions') == [2, 2, dim, dim],
             'RGBeta portal field dimensions changed')
    expected = f'{{Bar[H],H,Bar[S{which}],S{which}}}'
    _require(_clean(rec.get('fields_inputform')) == expected,
             f'RGBeta portal field order changed: {symbol}')
    coup = _clean(rec.get('native_coupling_inputform'))
    _require(re.fullmatch(r'6\*Matrix\['+symbol+r'\]\[##1\]&', coup) is not None,
             f'RGBeta native quartic coupling prefactor changed: {symbol}')
    _require(rec.get('native_self_conjugate') is True,
             f'RGBeta {symbol} unexpectedly complex')
    inv = _clean(rec.get('native_invariant_inputform'))
    _require('tGen[' in inv if kind == 'adjoint' else 'del[' in inv,
             f'RGBeta {symbol} registered invariant changed')
    _require(rec.get('native_permutation_count', 0) in range(1, 25),
             'Invalid RGBeta permutation count')
    return {'symbol': symbol,
            'registered_invariant_inputform': rec['native_invariant_inputform'],
            'native_coupling_inputform': rec['native_coupling_inputform'],
            'native_permutation_count': rec['native_permutation_count'],
            'native_generator_components_numerically_certified': False}


def _relative_error(left: np.ndarray | complex, right: np.ndarray | complex) -> float:
    diff = float(np.linalg.norm(np.asarray(left)-np.asarray(right)))
    return diff/max(1.,float(np.linalg.norm(left)),float(np.linalg.norm(right)))


def _portal_model_sector(model: str, which: int, d: int,
                         mrecords: list[dict], rrecords: list[dict],
                         rng: np.random.Generator, samples: int) -> dict:
    n = 1 if d == 1 else 2
    _require(len(mrecords) == n and len(rrecords) == n,
             f'{model}/S{which}: incomplete portal invariant sector')
    raw = [_component_tensor(rec, which, d, k)
           for k, rec in enumerate(mrecords, 1)]
    native = [_native_rgbeta(rec, which, d, kind)
              for rec, kind in zip(rrecords,
                                   (['singlet'] if d == 1 else ['singlet', 'adjoint']))]
    # A *real* lambda multiplies 1/2 PlusHc[phase * CG * (Hbar H)(Sbar S)].
    # This source factor is represented by the Hermitian projection.
    phases = [1j if rec['phase_inside_plus_hc'] == 'I' else 1.
              for rec in mrecords]
    projected = []
    raw_vectors = []
    for t, phase in zip(raw, phases):
        v = phase*t
        t_hc = v.conj().transpose(1,0,3,2)
        projected.append((v+t_hc)/2)
        raw_vectors.append(t.reshape(-1))
    raw_rank = int(np.linalg.matrix_rank(np.stack(raw_vectors,axis=1),tol=1e-10))
    _require(raw_rank == n, f'{model}/S{which}: registered raw CG rank {raw_rank}/{n}')
    mat = np.stack([np.r_[h.real.reshape(-1),h.imag.reshape(-1)]
                    for h in projected], axis=1)
    projected_rank = int(np.linalg.matrix_rank(mat,tol=1e-10))
    _require(projected_rank == n,
             f'{model}/S{which}: corrected real-Hermitian portal rank {projected_rank}/{n}; '
             'source phase or operator basis may be stale')
    ident = np.einsum('ij,ab->ijab',np.eye(2),np.eye(d))
    den = float(np.vdot(ident,ident).real)
    singlet = np.asarray([np.vdot(ident,t)/den for t in projected])
    _require(max(abs(singlet.imag)) < 1e-10, f'{model}/S{which}: complex singlet coefficients')
    singlet = singlet.real
    orth = [t-x*ident for t,x in zip(projected,singlet)]
    # A traceless SU(2) adjoint bilinear cannot contribute to either partial
    # trace. This checks the invariant decomposition independently of the fit.
    worst_trace = 0.
    for t in orth:
        tr_h = np.einsum('iiab->ab',t)
        tr_s = np.einsum('ijaa->ij',t)
        rel = max(_relative_error(tr_h,np.zeros_like(tr_h)),
                  _relative_error(tr_s,np.zeros_like(tr_s)))
        worst_trace = max(worst_trace, rel)
        _require(rel < 1e-10, f'{model}/S{which}: nontraceless adjoint remainder')
    if d == 1:
        oriented_adj = None
        transform = singlet.reshape(1,-1)
        target_symbol = [f'lambdaH{which}']
        adj_norm = None
        pivot = None
    else:
        pivot = next((i for i,t in enumerate(orth) if np.linalg.norm(t)>1e-10),None)
        _require(pivot is not None, f'{model}/S{which}: missing independent adjoint')
        # SU(2): tr(T_H^a T_H^b)=1/2 delta_ab; tr(T_S^a T_S^b)=
        # T(R) delta_ab, with T(2)=1/2 and T(3)=2. Thus the exact physical
        # source-defined RGBeta adjoint tensor norm is sqrt(3*T(H)*T(S)).
        adj_norm = float(np.sqrt(3*.5*({2:.5,3:2.}[d])))
        oriented_adj = adj_norm*orth[pivot]/float(np.linalg.norm(orth[pivot]))
        _require(_relative_error(np.vdot(ident,oriented_adj),0) < 1e-10,
                 f'{model}/S{which}: singlet-adjoint overlap')
        adj_row = []
        for t in orth:
            v = np.vdot(oriented_adj,t)/adj_norm**2
            _require(abs(v.imag)<1e-10,
                     f'{model}/S{which}: complex real-coupling adjoint coefficient')
            adj_row.append(float(v.real))
        transform=np.stack([singlet,np.asarray(adj_row)],axis=0)
        target_symbol=[f'lambdaH{which}',f'lambdaH{which}Adj_MatcheteOriented']
    worst=0.
    for k,t in enumerate(projected):
        reconstruct = transform[0,k]*ident
        if oriented_adj is not None:
            reconstruct=reconstruct+transform[1,k]*oriented_adj
        residual=_relative_error(t,reconstruct)
        worst=max(worst,residual)
        _require(residual < 1e-10, f'{model}/S{which}: portal tensor not in singlet+adjoint span')
    for _ in range(samples):
        h = rng.normal(size=2)+1j*rng.normal(size=2)
        s = rng.normal(size=d)+1j*rng.normal(size=d)
        coefficients = rng.normal(size=n)  # real Matchete native couplings
        actual=sum(z*np.einsum('ijab,i,j,a,b->',t,h.conj(),h,s.conj(),s)
                   for z,t in zip(coefficients,projected))
        as_rgb=transform@coefficients
        expect=as_rgb[0]*np.vdot(h,h)*np.vdot(s,s)
        if oriented_adj is not None:
            expect+=as_rgb[1]*np.einsum('ijab,i,j,a,b->',
                                       oriented_adj,h.conj(),h,s.conj(),s)
        worst=max(worst,_relative_error(actual,expect))
        _require(_relative_error(actual,expect)<1e-10,
                 f'{model}/S{which}: complex-field real-action holdout failed')
    return {
        'which': which, 'scalar_dimension':d,'raw_CG_rank':raw_rank,
        'corrected_real_Hermitian_rank':projected_rank,
        'native_matchete_coupling_order':[r['coupling'] for r in mrecords],
        'candidate_rgbeta_coupling_order':target_symbol,
        'candidate_rgbeta_from_native_matchete_real_matrix':transform.tolist(),
        'native_matchete_phases_inside_half_PlusHc':[r['phase_inside_plus_hc']
                                                      for r in mrecords],
        'traceless_adjoint_partial_trace_relative_residual':worst_trace,
        'expected_RGBeta_adjoint_generator_tensor_frobenius_norm':adj_norm,
        'matchete_oriented_adjoint_reference_pivot_1based':None if pivot is None else pivot+1,
        'native_rgbeta_registered_quartics':native,
        'max_tensor_and_real_complex_potential_relative_residual':worst,
        'candidate_orientation_sign_relative_to_native_RGBeta_tGen':
            None if d == 1 else 'UNDETERMINED_PLUS_OR_MINUS',
        'native_RGBeta_component_generators_extracted':False,
        'physical_coupling_dictionary_certified':False,
    }


def _matchete_build_status(model: str, row: dict) -> dict:
    """Allow B/C(-1) *symbolic* CG evidence, not unrotated mass-basis matching.

    The native builder explicitly labels these two general-quadratic models
    RequiresMixedScalarMassDiagonalization. Never relabel them as Success.
    All other unsuccessful native builds remain hard errors.
    """
    status = row.get('native_build_status')
    full = row.get('full_lagrangian_validation')
    ready = row.get('matching_mass_basis_ready')
    mixed = row.get('mixed_scalar_mass_basis_blocked')
    blocked = row.get('additional_quadratics_heavy_mass_basis_blocked')
    inventory = row.get('portal_cg_operator_inventory_verified')
    _require(inventory is True, f'{model}: native portal assembly not verified')
    if status == 'Success':
        _require(full is True and ready is True and mixed is False and
                 blocked == 0,
                 f'{model}: Success but mass-basis/full-validation metadata inconsistent')
    elif status == 'RequiresMixedScalarMassDiagonalization':
        _require(model in {'B_-1','C_-1'} and full is False and
                 ready is False and mixed is True and blocked == 1,
                 f'{model}: unsupported mixed-mass tensor-only exception')
    else:
        _require(False, f'{model}: native Lagrangian status {status!r} is blocked')
    return {
        'native_build_status': status,
        'full_lagrangian_validation': full,
        'matching_mass_basis_ready': ready,
        'additional_quadratics_heavy_mass_basis_blocked': blocked,
        'portal_CG_extraction_only': status != 'Success',
    }


def audit(mdoc: dict, rdoc: dict, *, sha_m: str, sha_r: str,
          seed: int=5261, samples: int=5) -> dict:
    _require(mdoc.get('format') == MAT_FORMAT and rdoc.get('format') == RGB_FORMAT,
             'Unexpected native-portal evidence format')
    _require(set(mdoc.get('models',{}))==ALL16 and set(rdoc.get('models',{}))==ALL16,
             'Missing 16-model evidence')
    _require(isinstance(samples,int) and 1<=samples<=25, 'Invalid sample count')
    rng = np.random.default_rng(seed)
    models={};total=0;worst=0.;anti=0
    mixed_models=[]
    for model in sorted(ALL16):
        cls = model[0]; d1,d2,df,alphas=REP[cls]
        alpha = int(model.split('_')[1])
        _require(alpha in alphas, 'Invalid model alpha')
        mr = mdoc['models'][model]; rr=rdoc['models'][model]
        build = _matchete_build_status(model, mr)
        if build['portal_CG_extraction_only']:
            mixed_models.append(model)
        _require(mr.get('dimensions')==[d1,d2,df] and
                 rr.get('dimensions')==[d1,d2,df] and
                 mr.get('alpha')==alpha and rr.get('alpha')==alpha,
                 f'{model}: native exporter representation mismatch')
        mrecords=mr.get('portals');rrecords=rr.get('portals')
        _require(isinstance(mrecords,list) and isinstance(rrecords,list),
                 f'{model}: absent portal lists')
        rows=[]
        for which,d in ((1,d1),(2,d2)):
            m=[t for t in mrecords if t.get('which')==which]
            r=[t for t in rrecords if t.get('which')==which]
            n=1 if d==1 else 2
            _require(len(m)==n and len(r)==n, f'{model}/S{which}: portal cardinality')
            m.sort(key=lambda x:x.get('invariant',-1))
            r.sort(key=lambda x:0 if x.get('kind')=='singlet' else 1)
            if model=='A_-4' and which==2:
                _require(m[0].get('phase_inside_plus_hc')=='I',
                         'A_-4 anti-Hermitian portal correction absent')
            row=_portal_model_sector(model,which,d,m,r,rng,samples)
            rows.append(row); total+=n
            anti+=sum(rec.get('phase_inside_plus_hc')=='I' for rec in m)
            worst=max(worst,row['max_tensor_and_real_complex_potential_relative_residual'])
        _require(len(mrecords)==sum(1 if d==1 else 2 for d in (d1,d2)) and
                 len(rrecords)==len(mrecords),f'{model}: extra unexpected portal')
        models[model]={'dimensions':[d1,d2,df],'alpha':alpha,
                       'matchete_native_build_validation':build,'sectors':rows}
    _require(total==58,'Expected 58 all-model Matchete portal coupling slots')
    _require(sorted(mdoc.get('models_requiring_mixed_scalar_mass_diagonalization',[]))
             ==sorted(mixed_models),
             'Native mixed-mass status inventory disagrees with model evidence')
    _require(mdoc.get('mixed_scalar_mass_basis_matching_authorized') is False,
             'Native mixed-mass matching unexpectedly authorized')
    return {
        'format':OUT_FORMAT, 'status':STATUS,
        'model_count':16,'portal_coupling_slots':total,
        'models_requiring_mixed_scalar_mass_diagonalization':mixed_models,
        'mixed_scalar_mass_basis_matching_authorized':False,
        'native_anti_Hermitian_CGs_given_I_inside_PlusHc':anti,
        'max_relative_tensor_and_potential_residual':worst,
        'models':models,
        'input_sha256':{'matchete_portals':sha_m,'rgbeta_portals':sha_r},
        'source_normalization':{
            'Matchete':'real lambdaH{which}InvN /2 PlusHc[phase*CG*Hbar H Sbar S]',
            'RGBeta':'AddQuartic[lambdaH{which} or Adj, {Bar H,H,Bar S,S}], '
                     'self-conjugate; Lam/24 and native 6*Matrix coefficient',
            'RGBeta_source_SU2_generator_trace_indices':{'fund':.5,'S2':2.},
        },
        'scope':{
            'native_registered_Matchete_portal_CGs_and_source_Hermitian_phase':True,
            'RGBeta_native_registered_symbolic_portal_functions':True,
            'real_coupling_subspace_and_conditional_generator_norm':True,
            'physical_RGBeta_adjoint_orientation_sign_certified':False,
            'native_RGBeta_numeric_generator_components_extracted':False,
            'native_RGBeta_real_scalar_coordinate_intertwiner_certified':False,
            'all_scalar_quartic_families_certified':False,
            'GammaCC_to_Weyl_normalization_certified':False,
            'Wilson_physical_inputs_certified':False,
            'B_C_minus1_mixed_scalar_mass_diagonalization_completed':False,
            'mixed_mass_model_portal_tensors_are_not_matching_certificates':True,
        },
    }


def audit_from_paths(mat: Path, rgb: Path, root: Path, **options) -> dict:
    mp=mat.read_bytes();rp=rgb.read_bytes()
    mdoc=json.loads(mp.decode('utf-8-sig'))
    rdoc=json.loads(rp.decode('utf-8-sig'))
    _sha_files(root, mdoc, SOURCE_MAT)
    _sha_files(root, rdoc, SOURCE_RGB)
    return audit(mdoc,rdoc,sha_m=sha256(mp).hexdigest(),
                 sha_r=sha256(rp).hexdigest(),**options)
