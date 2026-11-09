"""Standalone B handoff tests: no Wolfram, GitHub or local Batch 1 files needed."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from Numerical.parallel_t3.background.contracts import hash_file, save_json, tag, HEAVY
from Numerical.parallel_t3.background.handoff import audit_all, verify_one


def fixture_run(root: Path, model: str, *, hard=False) -> Path:
    run='B_all16_qualitative_v1'
    folder=root/'output/parallel_t3/background'/run/tag(model)
    folder.mkdir(parents=True,exist_ok=True)
    uv=np.geomspace(1e6,1e5,4)
    eft=np.geomspace(1e5,1e3,6)
    common=['gY','g2','g3','yu','yd','ye','lambdaH',
            'mS1Sq','mS2Sq','lambdaS1','lambdaS2']
    heavy=['MF','y1','y2']
    arrays={'mu_uv_gev':uv,'mu_eft1_gev':eft,
            't_uv':np.log(uv),'t_eft1':np.log(eft)}
    for name in common+heavy:
        base=np.eye(3, dtype=np.complex128) * (
            1e5 if name=='MF' else .02) if name in ('yu','yd','ye','MF','y1','y2') else np.asarray(.08+0j)
        if name in ('mS1Sq','mS2Sq'):
            base=np.asarray(1.2e6+0j)
        arrays[f'uv__{name}']=np.broadcast_to(base,(len(uv),)+base.shape).copy()
        if name in common:
            arrays[f'eft1__{name}']=np.broadcast_to(base,(len(eft),)+base.shape).copy()
    np.savez_compressed(folder/'background.npz',**arrays)
    bm={
        'schema':'t3-parallel-B-background-v1', 'model':model,
        'status':'CandidateRenormalisableBackgroundSuccess',
        'input_scales_gev':{'UV':1e6,'F':1e5,'S':1e3},
        'beta_symbol_order_uv':common+heavy,
        'beta_symbol_order_eft1':common,
        'UV_endpoint_sha256':'fedcba111122223333444455556666778899aaaabbbbccccddddeeeeffff0000',
        'background_npz_sha256':hash_file(folder/'background.npz'),
        'full_general_symmetry_certified':False,
        'matchete_rgbeta_coupling_dictionary_certified':False,
        'dimension_five_running_performed':False,
        'scalar_matching_performed':False,
        'scalar_coordinates_at_MS':{k:None for k in common if k.startswith(('lambda','mS','muT3'))},
        'common_state_continuity_at_MF':True,
        'last_eft1_sample_is_MS':True,
        'neutral_self_bilinear_complex_phase':
            'not_certified_by_native_real_splitting',
        'input_provenance':{'local_sources':{}},
    }
    save_json(folder/'background_manifest.json',bm)
    n = 12 if model.startswith('E') else 8
    c = np.zeros((3,3,2,2,n,n), dtype=np.complex128)
    z=1e-10*(.25+.17j)
    for (p,q,i,j,a,b) in ((0,1,0,1,0,1),(1,0,1,0,1,0),
                          (0,1,0,1,1,0),(1,0,1,0,0,1)):
        c[p,q,i,j,a,b]=z
    by=np.stack((c*.4,c*.6))
    b_arrays={
        'C_LLSS_formal':c,
        'C_LLSS_by_source_formal':by,
        'source_indices':np.asarray([101,103],dtype=np.int64),
        'heavy_mass_eigenvalues_gev':np.asarray([9.1e4,1e5,1.09e5],float),
    }
    if hard:
        b_arrays['C5_hard_external']=np.diag([1e-13+2e-14j,0j,2e-13]).astype(np.complex128)
    np.savez_compressed(folder/'fermion_boundary_formal.npz',**b_arrays)
    src=('output/general_observable/fermion_first_hermitian_portal_v1/A_m4/'
         'data/eft1_after_F_wilson_seed.json' if model=='A_-4' else
         f'output/general_observable/ordered_wilson_tensor_v4/{tag(model)}/source_flavour.json')
    source=root/src
    source.parent.mkdir(parents=True,exist_ok=True)
    source.write_text('test fixture; not physics',encoding='utf-8')
    fm={
        'schema':'t3-parallel-B-fermion-boundary-v1', 'model':model,
        'status':'FormalNativeMixedXYFermionTreeCandidate' if model in ('B_-1','C_-1')
            else 'FormalSourceExactFermionTreeLLSS',
        'formal_LLsS_tree_available':True,
        'mu_F_gev':1e5,'same_run_uv_state_sha256':bm['UV_endpoint_sha256'],
        'fermion_boundary_npz_sha256':hash_file(folder/'fermion_boundary_formal.npz'),
        'hard_C5_status':'external_independent_one_loop_candidate' if hard else 'hard_C5_unavailable',
        'hard_C5_is_physically_zero':None,
        'complete_one_loop_matching':False,
        'gammaCC_to_weyl_signed_conversion_certified':False,
        'all_three_heavy_fermions_decoupled_at_single_MF':True,
        'heavy_fermions_integrated_out_together':3,
        'active_BSM_fields_below_MF':['S1','S2'],
        'source_operator_count':2,
        'source_operator_indices':[101,103],
        'source_matchete_seed':src,
        'source_matchete_seed_sha256':hash_file(source),
        'axes_C_LLSS_formal':{
            'flavour_1':['e','mu','tau'], 'flavour_2':['e','mu','tau'],
            'lepton_su2_1':['1','2'], 'lepton_su2_2':['1','2'],
            'scalar_1':[str(i) for i in range(n)],
            'scalar_2':[str(i) for i in range(n)],
        },
        'mixed_scalar_mass_matrix_scale':'MF',
        'mixed_scalar_eigenbasis_status':'source_to_RGBeta_candidate',
    }
    save_json(folder/'fermion_boundary_manifest.json',fm)
    return folder


def test_nonmixed_readonly_handoff(tmp_path):
    path=fixture_run(tmp_path,'A_-4')
    before={f.name:hash_file(f) for f in path.iterdir() if f.is_file()}
    row=verify_one(tmp_path,'B_all16_qualitative_v1','A_-4')
    after={f.name:hash_file(f) for f in path.iterdir() if f.is_file()}
    assert before==after
    assert row['status']=='Verified_ProvisionalBInputsForAandC'
    assert row['workstream_A']['native_eft1_background_ready']
    assert row['workstream_C']['native_scalar_coordinate_trajectory_ready']
    assert row['workstream_C']['scalar_mass_and_quartic_component_tensor_certified'] is False
    assert row['hard_C5_is_physically_zero'] is None
    assert row['formal_fermion_boundary']['axes_shape']==[3,3,2,2,8,8]
    assert len(row['beta_coordinate_inventory']['fermion_keys_removed_at_MF'])==3
    assert set(HEAVY).isdisjoint(row['beta_coordinate_inventory']['EFT1'])


def test_mixed_and_independent_hard_candidate(tmp_path):
    fixture_run(tmp_path,'B_-1',hard=True)
    x=verify_one(tmp_path,'B_all16_qualitative_v1','B_-1')
    assert x['hard_C5_status']=='external_independent_one_loop_candidate'
    assert x['conditional_special_mass']=='mixed_XY_source_to_RGBeta_eigenbasis_not_certified'
    assert x['workstream_A']['canonical_LLSS_and_C5_renormalization_completed'] is False


def test_source_and_same_run_hash_fail_closed(tmp_path):
    folder=fixture_run(tmp_path,'A_-4')
    fm=folder/'fermion_boundary_manifest.json'
    from Numerical.parallel_t3.background.contracts import read_json
    f=read_json(fm)
    f['same_run_uv_state_sha256']='someone_else'
    save_json(fm,f)
    row=audit_all(tmp_path,'B_all16_qualitative_v1',('A_-4',))['model_results']['A_-4']
    assert row['status'].startswith('Blocked_')
    assert 'different UV endpoint' in row['error']


def test_by_source_construction_fail_closed(tmp_path):
    folder=fixture_run(tmp_path,'A_-4')
    path=folder/'fermion_boundary_formal.npz'
    with np.load(path,allow_pickle=False) as f:
        data={k:np.asarray(f[k]) for k in f.files}
    data['C_LLSS_by_source_formal'][0,0,1,0,1,0,1] *= 9
    np.savez_compressed(path,**data)
    from Numerical.parallel_t3.background.contracts import read_json
    mf=folder/'fermion_boundary_manifest.json'
    meta=read_json(mf)
    meta['fermion_boundary_npz_sha256']=hash_file(path)
    save_json(mf,meta)
    row=audit_all(tmp_path,'B_all16_qualitative_v1',('A_-4',))['model_results']['A_-4']
    assert row['status'].startswith('Blocked_')
    assert 'reconstruct' in row['error']


def test_missing_c5_cannot_be_zero(tmp_path):
    folder=fixture_run(tmp_path,'A_-4')
    from Numerical.parallel_t3.background.contracts import read_json
    mf=folder/'fermion_boundary_manifest.json'
    meta=read_json(mf)
    meta['hard_C5_is_physically_zero']=True
    save_json(mf,meta)
    row=audit_all(tmp_path,'B_all16_qualitative_v1',('A_-4',))['model_results']['A_-4']
    assert row['status'].startswith('Blocked_')
    assert 'reinterpreted as zero' in row['error']


def test_partial_batch_and_B_owned_writes(tmp_path):
    fixture_run(tmp_path,'A_-4')
    out=tmp_path/'output/parallel_t3/background/B_all16_qualitative_v1/handoff_index.json'
    x=audit_all(tmp_path,'B_all16_qualitative_v1',('A_-4','B_-3'),output=out)
    assert x['passed_count']==1 and x['requested_count']==2
    assert x['model_results']['B_-3']['status'].startswith('Blocked_')
    assert out.is_file()
    with pytest.raises(ValueError,match='only write'):
        audit_all(tmp_path,'B_all16_qualitative_v1',('A_-4',),
                  output=tmp_path/'somewhere_else.json')


def test_provenance_source_drift_warning_not_physics_overwrite(tmp_path):
    folder=fixture_run(tmp_path,'A_-4')
    from Numerical.parallel_t3.background.contracts import read_json
    mf=folder/'background_manifest.json'
    meta=read_json(mf)
    code=tmp_path/'Lagrangian/interactions/ScalarPotential.wl'
    code.parent.mkdir(parents=True,exist_ok=True)
    code.write_text('original source',encoding='utf-8')
    meta['input_provenance']['local_sources']={
        'Lagrangian/interactions/ScalarPotential.wl': {
            'required_for_B_reuse':True,'same_as_batch1':True,'sha256':hash_file(code)}}
    save_json(mf,meta)
    code.write_text('modified later',encoding='utf-8')
    row=verify_one(tmp_path,'B_all16_qualitative_v1','A_-4')
    assert len(row['source_drift_warnings'])==1
    assert row['status']=='Verified_ProvisionalBInputsForAandC'


def test_all16_synthetic_saved_backgrounds_handoff(tmp_path):
    from Numerical.parallel_t3.background.contracts import MODELS
    for model in MODELS:
        fixture_run(tmp_path,model)
    result=audit_all(tmp_path,'B_all16_qualitative_v1')
    assert result['passed_count']==16
    assert len(result['model_results'])==16
    assert all(result['model_results'][m]['status']=='Verified_ProvisionalBInputsForAandC'
               for m in MODELS)
    for model in ('B_-1','C_-1'):
        assert result['model_results'][model]['conditional_special_mass'].startswith('mixed_XY')
    for model in ('A_-2','A_+0','D_-2','D_+0','E_-2','E_+0'):
        assert result['model_results'][model]['conditional_special_mass'].startswith('neutral_self')
