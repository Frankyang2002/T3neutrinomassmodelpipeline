"""Fail-closed, read-only validation of saved Workstream B -> A/C inputs.

This module NEVER evolves couplings, converts Matchete to Weyl, constructs
physical scalar tensors, or assumes that a missing hard C5 is zero.  The
only optional write is a small, explicitly provisional handoff index in B's
own output directory.  Designed to consume existing completed Batch B runs.
"""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import math
import traceback
from typing import Any

import numpy as np

from .api import Background, load_formal_boundary
from .contracts import (
    FAMILIES, HEAVY, MIXED, MODELS, SELF, ContractError,
    hash_file, read_json, safe_run_id, save_json, tag,
)

HANDOFF_SCHEMA = 't3-parallel-B-handoff-audit-v1'
BACKGROUND_STATUS = 'CandidateRenormalisableBackgroundSuccess'
FORMAL_STATUS = 'FormalSourceExactFermionTreeLLSS'
MIXED_STATUS = 'FormalNativeMixedXYFermionTreeCandidate'
BOUNDARY_SCHEMA = 't3-parallel-B-fermion-boundary-v1'


def _require(condition: bool, message: str) -> None:
    if not bool(condition):
        raise ContractError(message)


def _nearly(a: float, b: float, *, rtol: float = 2e-12) -> bool:
    return bool(np.isclose(float(a), float(b), rtol=rtol, atol=1e-10))


def _relative(root: Path, path: Path) -> str:
    return path.resolve().relative_to(root.resolve()).as_posix()


def _scalar_names(names: list[str]) -> list[str]:
    return [n for n in names if n.startswith(('lambda', 'mS', 'muT3'))]


def _source_drift(root: Path, meta: dict) -> list[str]:
    """Only warn on post-run edits; archived SHA provenance remains immutable.

    A/C can change their own files after B ran without invalidating the saved
    B trajectories. The current on-disk differences must still be visible.
    """
    warnings = []
    source = meta.get('input_provenance', {}).get('local_sources', {})
    for rel, row in source.items():
        if not isinstance(row, dict) or not row.get('required_for_B_reuse'):
            continue
        if not row.get('same_as_batch1'):
            raise ContractError('Saved B run did not pass Batch1 provenance for '+rel)
        candidate = root/rel
        if not candidate.is_file():
            warnings.append(f'Current source absent after B run: {rel}')
        elif hash_file(candidate) != row.get('sha256'):
            warnings.append(f'Current source changed after B run: {rel}')
    return warnings


def verify_one(root: Path, run_id: str, model: str) -> dict[str, Any]:
    """Validate saved B data for one of the 16 ordinary T3 model points.

    Raises ContractError if input cannot be consumed safely.  All 16
    requested models are treated independently by :func:`audit_all`.
    """
    root = Path(root).resolve()
    safe_run_id(run_id)
    _require(model in MODELS, f'Unknown T3 model: {model}')
    folder = root/'output'/'parallel_t3'/'background'/run_id/tag(model)
    _require(folder.is_dir(), f'{model}: saved B output directory missing: {folder}')
    bg = Background.load(folder)  # verifies NPZ SHA, complete keys, no F in EFT1, continuity
    bm = bg.manifest
    fm = read_json(folder/'fermion_boundary_manifest.json')
    _require(bm.get('model') == fm.get('model') == model,
             f'{model}: B and F manifests have inconsistent model IDs')
    _require(bm.get('status') == BACKGROUND_STATUS,
             f'{model}: renormalisable background is not successful')
    _require(fm.get('schema') == BOUNDARY_SCHEMA,
             f'{model}: unknown formal boundary schema')
    expected = MIXED_STATUS if model in MIXED else FORMAL_STATUS
    _require(fm.get('status') == expected and fm.get('formal_LLsS_tree_available') is True,
             f'{model}: missing expected evaluated formal F-tree status')
    _require(bm.get('full_general_symmetry_certified') is False and
             bm.get('matchete_rgbeta_coupling_dictionary_certified') is False,
             f'{model}: unexpected unsupported physical-general/basis certification')
    _require(bm.get('dimension_five_running_performed') is False and
             bm.get('scalar_matching_performed') is False,
             f'{model}: saved background incorrectly claims A/C work')
    _require(fm.get('gammaCC_to_weyl_signed_conversion_certified') is False and
             fm.get('complete_one_loop_matching') is False,
             f'{model}: false claim of physical spinor conversion/full hard matching')
    _require(fm.get('all_three_heavy_fermions_decoupled_at_single_MF') is True and
             fm.get('active_BSM_fields_below_MF') == ['S1','S2'],
             f'{model}: incorrect all-fermions-first EFT order')
    _require(fm.get('heavy_fermions_integrated_out_together') == 3,
             f'{model}: not all three heavy fermions matched together')
    _require(fm.get('same_run_uv_state_sha256') == bm.get('UV_endpoint_sha256'),
             f'{model}: formal LLSS boundary has a different UV endpoint')
    scales = bm.get('input_scales_gev', {})
    for name in ('UV','F','S'):
        _require(name in scales and math.isfinite(float(scales[name])) and float(scales[name])>0,
                 f'{model}: malformed matching scale {name}')
    _require(float(scales['UV'])>float(scales['F'])>float(scales['S']),
             f'{model}: wrong UV/F/S scale ordering')
    _require(_nearly(bg.arrays['mu_uv_gev'][0],scales['UV']) and
             _nearly(bg.arrays['mu_uv_gev'][-1],scales['F']) and
             _nearly(bg.arrays['mu_eft1_gev'][0],scales['F']) and
             _nearly(bg.arrays['mu_eft1_gev'][-1],scales['S']) and
             _nearly(fm.get('mu_F_gev'),scales['F']),
             f'{model}: saved trajectories and F boundary do not meet at declared scales')
    for stage in ('uv','eft1'):
        mus = bg.arrays[f'mu_{stage}_gev']
        t = bg.arrays[f't_{stage}']
        _require(t.shape == mus.shape and np.allclose(t,np.log(mus),atol=5e-11,rtol=2e-12),
                 f'{model}: t samples disagree with ln(mu) in {stage}')
    names_uv = bm['beta_symbol_order_uv']
    names_eft = bm['beta_symbol_order_eft1']
    _require(len(names_uv)==len(set(names_uv)) and len(names_eft)==len(set(names_eft)),
             f'{model}: repeated beta coordinate labels')
    _require(set(names_uv)-HEAVY==set(names_eft),
             f'{model}: EFT1 not the exact nonfermion UV beta inventory')
    _require(set(('MF','y1','y2'))<=set(names_uv) and not (HEAVY & set(names_eft)),
             f'{model}: invalid heavy fields active at threshold')
    scalars = _scalar_names(names_eft)
    _require({'mS1Sq','mS2Sq','lambdaS1','lambdaS2'}<=set(scalars),
             f'{model}: missing scalar quadratic/quartic RGEs in EFT1')
    _require(set(bm.get('scalar_coordinates_at_MS',{})) == set(scalars),
             f'{model}: saved MS coordinate inventory incomplete')
    _require(bm.get('common_state_continuity_at_MF') is True and
             bm.get('last_eft1_sample_is_MS') is True,
             f'{model}: saved continuity/end-scale flags not set')

    arrays, f = load_formal_boundary(folder)
    c = arrays['C_LLSS_formal']
    d1, d2, _ = FAMILIES[model[0]][0]
    nreal = 2*(d1+d2)
    _require(c.shape == (3,3,2,2,nreal,nreal),
             f'{model}: unexpected real scalar axes in LLSS {c.shape}')
    _require(c.dtype == np.complex128 and
             np.allclose(c, c.transpose(1,0,3,2,4,5), rtol=2e-9,atol=1e-28) and
             np.allclose(c,c.swapaxes(4,5),rtol=2e-9,atol=1e-28),
             f'{model}: formal LLSS field exchange symmetry broken')
    by = arrays.get('C_LLSS_by_source_formal')
    idx = arrays.get('source_indices')
    masses = arrays.get('heavy_mass_eigenvalues_gev')
    _require(by is not None and by.ndim==7 and by.shape[1:]==c.shape and
             by.dtype==np.complex128,
             f'{model}: incomplete by-source formal LLSS components')
    _require(idx is not None and idx.shape == (by.shape[0],) and
             np.issubdtype(idx.dtype,np.integer) and len(set(map(int,idx)))==len(idx),
             f'{model}: source operator labels absent, repeated or invalid')
    _require(int(fm.get('source_operator_count',-1))==len(idx) and
             list(map(int,idx)) == fm.get('source_operator_indices'),
             f'{model}: source-by-source Matchete inventory disagrees with manifest')
    _require(np.allclose(by.sum(axis=0),c,rtol=2e-10,atol=1e-28),
             f'{model}: source LLSS components do not reconstruct total')
    _require(masses is not None and masses.shape==(3,) and np.isfinite(masses).all()
             and np.min(masses)>0,
             f'{model}: invalid heavy fermion singular/Takagi mass eigenvalues')
    axes = fm.get('axes_C_LLSS_formal')
    _require(isinstance(axes,dict) and
             all(len(axes.get(k,[])) == s for k,s in (
                ('flavour_1',3),('flavour_2',3),('lepton_su2_1',2),
                ('lepton_su2_2',2),('scalar_1',nreal),('scalar_2',nreal))),
             f'{model}: missing formal flavour/SU2/scalar axis labels')
    _require(axes['scalar_1']==axes['scalar_2'],
             f'{model}: mismatched heavy scalar pair axes')
    source_path = root/fm['source_matchete_seed']
    _require(source_path.is_file() and hash_file(source_path)==fm.get('source_matchete_seed_sha256'),
             f'{model}: source Matchete F-tree file no longer matches B boundary')
    if model=='A_-4':
        _require('fermion_first_hermitian_portal_v1' in source_path.as_posix(),
                 'A_-4: stale pre-portal Matchete source selected')
    if model in MIXED:
        _require(fm.get('mixed_scalar_mass_matrix_scale')=='MF' and
                 'candidate' in fm.get('mixed_scalar_eigenbasis_status',''),
                 f'{model}: missing conditional native X/Y scalar mixing evidence')
    if model in SELF:
        _require(bm.get('neutral_self_bilinear_complex_phase') ==
                 'not_certified_by_native_real_splitting',
                 f'{model}: unsupported neutral-phase certification')
    hard = fm.get('hard_C5_status')
    _require(hard in ('hard_C5_unavailable','hard_C5_requested_but_blocked',
                      'external_independent_one_loop_candidate'),
             f'{model}: unrecognized hard C5 status')
    has_hard = 'C5_hard_external' in arrays
    _require(has_hard == (hard=='external_independent_one_loop_candidate'),
             f'{model}: inconsistent optional hard C5 artifact/status')
    if has_hard:
        h=arrays['C5_hard_external']
        _require(h.shape==(3,3) and h.dtype==np.complex128 and
                 np.isfinite(h).all() and np.allclose(h,h.T,rtol=1e-10,atol=1e-25),
                 f'{model}: malformed independently calculated hard C5 matrix')
    _require(fm.get('hard_C5_is_physically_zero') is None,
             f'{model}: missing hard Weinberg was reinterpreted as zero')

    bg_path = folder/'background.npz'
    f_path = folder/'fermion_boundary_formal.npz'
    receipt: dict[str,Any] = {
        'schema':HANDOFF_SCHEMA, 'model':model, 'tag':tag(model),
        'status':'Verified_ProvisionalBInputsForAandC',
        'field_content_eft1':['SM','S1','S2'],
        'scales_gev':{'UV':float(scales['UV']),'F':float(scales['F']),
                      'S':float(scales['S'])},
        'files':{n:{'path':_relative(root, folder/n),'sha256':hash_file(folder/n)}
                 for n in ('background.npz','background_manifest.json',
                           'fermion_boundary_formal.npz','fermion_boundary_manifest.json')},
        'same_run_uv_state_sha256':bm['UV_endpoint_sha256'],
        'beta_coordinate_inventory':{
            'UV':list(names_uv),'EFT1':list(names_eft),'EFT1_scalar':scalars,
            'fermion_keys_removed_at_MF':sorted(set(names_uv)-set(names_eft)),
        },
        'sample_counts':{'UV':int(bg.arrays['mu_uv_gev'].size),
                         'EFT1':int(bg.arrays['mu_eft1_gev'].size)},
        'formal_fermion_boundary':{
            'status':fm['status'], 'axes_shape':list(c.shape),
            'axis_labels':axes, 'source_operator_count':int(by.shape[0]),
            'norm_gev_inverse':float(np.linalg.norm(c)),
            'source_matchete_seed_sha256':fm['source_matchete_seed_sha256'],
            'heavy_fermion_eigenmasses_at_MF_gev':np.asarray(masses,float).tolist(),
            'formal_matchete_coordinates_only':True,
            'signed_gammaCC_Weyl_factor_known':False,
        },
        'hard_C5_status':hard,
        'hard_C5_is_physically_zero':None,
        'conditional_special_mass':(
            'mixed_XY_source_to_RGBeta_eigenbasis_not_certified' if model in MIXED
            else 'neutral_self_bilinear_complex_phase_not_certified' if model in SELF
            else None),
        'workstream_A':{
            'native_eft1_background_ready':True,
            'formal_LLSS_source_ready':True,
            'canonical_LLSS_and_C5_renormalization_completed':False,
            'requires_signed_operator_map_and_closed_ADM':True,
        },
        'workstream_C':{
            'native_scalar_coordinate_trajectory_ready':True,
            'scalar_mass_and_quartic_component_tensor_certified':False,
            'simultaneous_scalar_matching_completed':False,
            'requires_matchete_RGBeta_mass_and_quartic_basis_map':True,
        },
        'source_drift_warnings':_source_drift(root,bm),
        'physical_neutrino_mass_claim':False,
    }
    return receipt


def audit_all(root: Path, run_id: str, models: tuple[str,...]=MODELS,
              *, output: Path | None = None) -> dict[str,Any]:
    """Independently audit models. Output only a B-owned receipt index, no physics edits."""
    root = Path(root).resolve()
    safe_run_id(run_id)
    requested=tuple(models)
    _require(len(set(requested))==len(requested) and set(requested)<=set(MODELS),
             'Requested handoff model inventory invalid')
    results={}
    for model in requested:
        try:
            results[model]=verify_one(root,run_id,model)
        except Exception as exc:
            results[model]={
                'schema':HANDOFF_SCHEMA,'model':model,'tag':tag(model),
                'status':'Blocked_'+type(exc).__name__,
                'error':str(exc),'traceback_tail':traceback.format_exc(limit=3),
                'workstream_A_ready':False,'workstream_C_background_ready':False,
            }
    passed=sum(row['status']=='Verified_ProvisionalBInputsForAandC' for row in results.values())
    out={
        'schema':'t3-parallel-B-handoff-index-v1',
        'created_utc':datetime.now(timezone.utc).isoformat(),
        'run_id':run_id,'models_requested':list(requested),
        'passed_count':passed,'requested_count':len(requested),
        'physical_C5_certified_count':0,
        'coupled_dimension_five_running_performed':False,
        'simultaneous_scalar_matching_performed':False,
        'model_results':results,
    }
    if output is not None:
        destination=Path(output)
        if not destination.is_absolute():destination=root/destination
        allowed=root/'output'/'parallel_t3'/'background'/run_id
        _require(destination.parent.resolve()==allowed.resolve(),
                 'Handoff audit may only write directly into its own B run directory')
        _require(destination.name=='handoff_index.json',
                 'Handoff audit writes only handoff_index.json (not saved physics)')
        save_json(destination,out)
    return out
