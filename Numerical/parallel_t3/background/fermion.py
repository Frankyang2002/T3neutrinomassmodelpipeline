"""Same-run source-resolved Matchete *formal* LLSS boundary at one common MF.

No signed GammaCC/Weyl conversion, no hard C5 invention, no scalar decoupling.
Only the nonmixed source/CG backend is used. Mixed B_-1/C_-1 fail closed.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .contracts import (MIXED, ContractError, check_finite, hash_file,
                        read_json, tag)


def source_directory(root: Path, batch: Path, model: str) -> Path:
    if model == 'A_-4':
        return root / 'output/general_observable/fermion_first_hermitian_portal_v1/A_m4'
    return batch / 'matchete' / tag(model)


def evaluate_formal_llss(model: str, uv_at_mf: Any, *, root: Path,
                         batch: Path, record: Mapping[str, Any], mf: float
                         ) -> tuple[dict[str, np.ndarray], dict]:
    """Use Matchete's ACTUAL F-tree term and CG factors with the evolved UV state.

    Returns (arrays, manifest), with arrays['C_LLSS_formal'] ordered
    (lepton flavour 3, lepton flavour 3, L SU2 2, L SU2 2,
     real heavy S component 2(dS1+dS2), real heavy S component).
    Tensor values are *formal*, NOT canonical signed Weyl coefficients.
    """
    if model in MIXED:
        raise ContractError(f'{model}: mixed scalar requires dedicated native X/Y '
                            'source and Workstream C mass eigenbasis at MF; '
                            'restricted/nonmixed source is forbidden')
    source = source_directory(root, batch, model)
    native_path = source / 'native_status.json'
    seed_path = source / 'data/eft1_after_F_wilson_seed.json'
    if not native_path.is_file() or not seed_path.is_file():
        raise ContractError(f'{model}: native Matchete F-tree seed absent at {source}')
    native = read_json(native_path)
    if native.get('status') != 'FThresholdTreeWilsonExported':
        raise ContractError(f'{model}: native Matchete matching not successful')
    seed = read_json(seed_path)
    if seed.get('IntegratedFields') != ['F'] or seed.get('ActiveHeavyFields') != ['S1', 'S2']:
        raise ContractError(f'{model}: Matchete seed is not ALL F out, both scalars active')
    evidence = record.get('f_threshold', {})
    saved_seed = evidence.get('data/eft1_after_F_wilson_seed.json', {})
    if not evidence.get('fermions_out_scalars_kept') or hash_file(seed_path) != saved_seed.get('sha256'):
        raise ContractError(f'{model}: Matchete source hash/inventory mismatches Batch 1')
    if model == 'A_-4':
        qseed = source / 'data/eft1_after_F_scalar_quartic_seed.json'
        qref = evidence.get('data/eft1_after_F_scalar_quartic_seed.json', {})
        if not qseed.is_file() or hash_file(qseed) != qref.get('sha256'):
            raise ContractError('A_-4: corrected 9-quartic portal seed hash mismatch')
        portal = evidence.get('a4_portal', {})
        if not portal.get('both_validated') or portal.get('matched_quartics_missing') != []:
            raise ContractError('A_-4: corrected two-portal matched inventory not established')

    # Import the already implemented source-exact factoriser, heavy SVD/Takagi,
    # Matchete SU2/real scalar CG expander and polynomial orbit projector.
    # Do NOT replace with a guessed y1 M^-1 y2^T tensor or scalar proxy.
    from scripts.RunT3GeneralWilsonBoundary import (
        _validate_seed, flavour_matrix, encode_matrix)
    from Numerical.matching.GeneralT3FlavourC5 import prepare_heavy_flavour
    from Numerical.running.GeneralT3OrderedTensor import project_ordered_term, scalar_basis
    from Numerical.running.GeneralT3DynamicFermionWilson import _dense_operator

    terms = _validate_seed(seed, model)
    needed = {f['name'] for term in terms for f in term['factors']}
    masses, yukawas, heavy_diagnostics = prepare_heavy_flavour(
        uv_at_mf, needed | {'y1','y2'})
    if abs(float(uv_at_mf.mu_gev)-mf) > 1e-7 * mf:
        raise ContractError(f'{model}: Matchete boundary must use UV state at actual MF')
    operators = []
    for term in terms:
        term['flavour_matrix_axes_external_indices'] = [
            f['external_index'] for f in term['factors']]
        term['ordered_flavour_matrix_gev_inv'] = encode_matrix(
            flavour_matrix(term, masses, yukawas))
        operators.append(project_ordered_term(term, seed['CGRegistry'], model))
    by_source, indices, audit = _dense_operator(operators, terms, model)
    tensor = np.asarray(by_source.sum(axis=0), dtype=np.complex128)
    check_finite(tensor, by_source)
    n = 2*sum({'A': (1,3), 'B': (2,2), 'C': (2,2),
               'D': (3,1), 'E': (3,3)}[model[0]])
    if tensor.shape != (3, 3, 2, 2, n, n):
        raise ContractError(f'{model}: wrong formal LLSS tensor shape {tensor.shape}')
    if not np.allclose(tensor, tensor.transpose(1,0,3,2,4,5), atol=1e-28, rtol=1e-9):
        raise ContractError(f'{model}: failed lepton/spinor exchange symmetry')
    if not np.allclose(tensor, tensor.swapaxes(4,5), atol=1e-28, rtol=1e-9):
        raise ContractError(f'{model}: failed bosonic scalar exchange symmetry')
    scalar = scalar_basis(model)
    heavy_labels = list(getattr(scalar, 'labels', ()))[4:]
    # If the scalar basis API exposes no labels, keep an explicitly indexed
    # source convention; do not guess the relative complex-field order.
    axes = {
        'flavour_1': ['e','mu','tau'], 'flavour_2': ['e','mu','tau'],
        'lepton_su2_1': ['1','2'], 'lepton_su2_2': ['1','2'],
        'scalar_1': heavy_labels if len(heavy_labels) == n else
                    [f'Matchete_real_scalar_slot_{i+5}' for i in range(n)],
        'scalar_2': heavy_labels if len(heavy_labels) == n else
                    [f'Matchete_real_scalar_slot_{i+5}' for i in range(n)]}
    signature = hashlib.sha256(seed_path.read_bytes() + tensor.tobytes()).hexdigest()
    return {'C_LLSS_formal': tensor, 'C_LLSS_by_source_formal': by_source,
            'heavy_mass_eigenvalues_gev': np.asarray(masses, dtype=np.float64),
            'source_indices': np.asarray(indices, dtype=np.int64)}, {
        'status': 'FormalSourceExactFermionTreeLLSS', 'model':model, 'mu_F_gev': mf,
        'source_matchete_seed': str(seed_path.relative_to(root)),
        'source_matchete_seed_sha256': hash_file(seed_path),
        'source_matchete_native_status_sha256': hash_file(native_path),
        'source_signature_sha256': signature,
        'flavour_basis': 'three lepton flavours; heavy basis rotated by verified SVD/Takagi',
        'axes_C_LLSS_formal': axes,
        'source_operator_indices': [int(i) for i in indices],
        'source_operator_count': len(indices),
        'heavy_mass_basis_diagnostics': str(heavy_diagnostics),
        'required_yukawa_channels': sorted(needed),
        'heavy_fermions_integrated_out_together': 3,
        'active_BSM_fields_below_MF': ['S1', 'S2'],
        'gammaCC_to_weyl_signed_conversion_certified': False,
        'physical_normalization_certified': False,
        'other_dimension_five_hard_operators_computed': False,
        'hard_C5_status': 'hard_C5_unavailable',
        'complete_one_loop_matching': False,
        'audit': audit,
    }


def evaluate_mixed_native_candidate(model: str, uv_at_mf: Any, eft1_at_mf: Any,
                                    *, root: Path, mf: float
                                    ) -> tuple[dict[str, np.ndarray], dict]:
    """Dedicated B/C X/Y native CG interface, never restricted matching.

    Reuses the *existing* GeneralT3MixedExactFermionWilson source driver and
    its exact Majorana heavy-flavour/SU(2) factors. Provides the evolved
    same-run scalar mass matrix at MF to its native X/Y diagonaliser. Any
    X/Y-to-RGBeta coupling/basis dictionary remains a conditional assumption.
    The archived mixed native source itself is fixed at MF=1e5 GeV.
    """
    if model not in MIXED:
        raise ContractError('Dedicated mixed adapter accepts B_-1/C_-1 only')
    if abs(mf-1e5) > 1e-9:
        raise ContractError(f'{model}: archived mixed Matchete source hardcodes '
                            'MF=1e5 GeV; regenerate exact source for other MF')
    from Numerical.running.GeneralT3MixedExactFermionWilson import mixed_exact_boundary
    from Numerical.running.IntermediateScalarRGBetaEvaluator import (
        intermediate_scalar_state_environment)
    scalars = intermediate_scalar_state_environment(eft1_at_mf)
    required = ('mS1Sq','mS2Sq','muT3GeneralS1S2')
    missing = set(required)-set(scalars)
    if missing:
        raise ContractError(f'{model}: missing native mixed scalar mass coordinates {sorted(missing)}')
    source, audit = mixed_exact_boundary(model, uv_at_mf,
                                        {k:scalars[k] for k in required},
                                        repo_root=root, mu_f=mf)
    tensor = np.asarray(source.raw, dtype=np.complex128)
    by_source = np.asarray(source.by_source, dtype=np.complex128)
    check_finite(tensor, by_source)
    if tensor.shape != (3,3,2,2,8,8):
        raise ContractError(f'{model}: wrong native X/Y tensor shape')
    if not np.allclose(tensor, tensor.transpose(1,0,3,2,4,5), rtol=1e-9, atol=1e-28):
        raise ContractError(f'{model}: mixed source exchange symmetry fails')
    if not np.allclose(tensor, tensor.swapaxes(4,5), rtol=1e-9, atol=1e-28):
        raise ContractError(f'{model}: mixed bosonic exchange symmetry fails')
    n = 8
    labels = list(source.labels)
    if len(labels) != n:
        labels = [f'Matchete_mixed_real_scalar_slot_{i+5}' for i in range(n)]
    return {'C_LLSS_formal': tensor, 'C_LLSS_by_source_formal': by_source,
            'source_indices': np.asarray(source.source_indices, dtype=np.int64),
            'heavy_mass_eigenvalues_gev': np.asarray(audit['heavy_mass_eigenvalues_gev'])}, {
        'status': 'FormalNativeMixedXYFermionTreeCandidate',
        'model':model, 'mu_F_gev':mf,
        'source_matchete_seed':str(source.path.relative_to(root)),
        'source_matchete_seed_sha256':hash_file(source.path),
        'source_signature_sha256':source.hash,
        'source_operator_indices':list(source.source_indices),
        'source_operator_count':len(source.source_indices),
        'axes_C_LLSS_formal': {
            'flavour_1':['e','mu','tau'],'flavour_2':['e','mu','tau'],
            'lepton_su2_1':['1','2'],'lepton_su2_2':['1','2'],
            'scalar_1':labels, 'scalar_2':labels},
        'heavy_fermions_integrated_out_together':3,
        'active_BSM_fields_below_MF':['S1','S2'],
        'mixed_scalar_mass_matrix_scale':'MF',
        'mixed_scalar_mass_matrix_input_coordinates':required,
        'mixed_scalar_eigenbasis_status':'native_XY_diagonaliser_at_MF; '
              'source_to_RGBeta_basis_identification_candidate_not_certified',
        'physical_normalization_certified':False,
        'gammaCC_to_weyl_signed_conversion_certified':False,
        'other_dimension_five_hard_operators_computed':False,
        'hard_C5_status':'hard_C5_unavailable',
        'complete_one_loop_matching':False,
        'mixed_native_diagnostics':audit}


def check_external_hard_matching(path: Path, manifest_path: Path, *, model: str,
                                 mf: float, uv_state_sha: str) -> tuple[np.ndarray, dict]:
    """Accept ONLY a documented independent calculation in the matching basis.

    An absent hard matching is NOT a zero coefficient; it is absent from output.
    External artifacts must carry current same-run UV state hash and exact MF.
    """
    if not path.is_file() or not manifest_path.is_file():
        raise ContractError('Optional hard matching requested but files are missing')
    meta = read_json(manifest_path)
    if (meta.get('model') != model or meta.get('mu_F_gev') != mf or
        meta.get('same_run_uv_state_sha256') != uv_state_sha or
        meta.get('method') != 'independent_one_loop_fermion_threshold_calculation' or
        meta.get('hard_c5_file_sha256') != hash_file(path)):
        raise ContractError('Hard matching is not independently documented for this UV state')
    if not isinstance(meta.get('basis'), str) or not meta['basis']:
        raise ContractError('Hard matching must state its operator/basis normalization')
    with np.load(path, allow_pickle=False) as f:
        if 'C5_hard' not in f:
            raise ContractError('Hard matching file lacks C5_hard')
        c5 = np.asarray(f['C5_hard'], dtype=np.complex128)
    if c5.shape != (3,3) or not np.allclose(c5, c5.T, rtol=1e-10, atol=1e-25):
        raise ContractError('Hard C5 must be 3x3 symmetric in flavour')
    check_finite(c5)
    return c5, meta
