"""New Stream3 tests: producer schemas, one-scale matching and observables.

All generated input tensors are explicitly SYNTHETIC fixtures. This suite does
NOT certify full physical A/B models or silently create genuine model results.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from Numerical.parallel_t3.threshold.fixtures import make_fixture,make_all_fixtures
from Numerical.parallel_t3.threshold.models import ALL_MODELS,n_real
from Numerical.parallel_t3.threshold.preflight import preflight_pair
from Numerical.parallel_t3.threshold.rough_adapter import AdapterBlocked,prepare_pair
from Numerical.parallel_t3.threshold.rough_numerics import neutrino_observables,run_rough_case

ROOT=Path(__file__).resolve().parents[1]


def _write_json(p:Path,meta:dict):
    p.write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')


def _make_legacy_A(original:Path):
    """Repackage real tensors from a synthetic C fixture into A's OLD schema."""
    meta=json.loads(original.read_text())
    with np.load(original.parent/'wilson_at_MS.npz',allow_pickle=False) as z:
        data={k:z[k].copy() for k in z.files}
    source={**{k:meta[k] for k in (
        'model','matching_context_id','matched_seed_id','lagrangian_source_id',
        'mu_f_gev','mu_s_gev','scheme','units','synthetic_fixture')},
        'format':'t3-parallel-a-wilson-output-v1',
        'all_F_integrated_out_at_MF':True,
        'both_S1_S2_remain_until_MS':True,
        'no_posthoc_fixed_order_leading_log':True,
        'C5_mutual_mixing_present':True,
        'C5_projection_definition':'C5[p,q]=4*C[L[p,1],L[q,1],H[2,R],H[2,R]]',
        'basis':{'basis_id':meta['scalar_basis_id'],
            'fermion_labels':[f'L[{i},{j}]' for i in (1,2,3) for j in (1,2)]+['ec[1]','ec[2]'],
            'scalar_labels':['H[1,R]','H[1,I]','H[2,R]','H[2,I]']+meta['scalar_source_real_labels']}}
    np.savez_compressed(original.parent/'wilson_at_MS.npz',
                        C_full=data['C_psi2phi2_ms'],
                        C5_EFT1_MS=data['C5_EFT1_ms'])
    _write_json(original,source)


def test_takagi_pmns_and_mass_splittings():
    # Choose a known, nontrivial complex symmetric Majorana matrix.
    angles=np.array([.58,.15,.79])
    c12,s12=np.cos(angles[0]),np.sin(angles[0])
    c13,s13=np.cos(angles[1]),np.sin(angles[1])
    c23,s23=np.cos(angles[2]),np.sin(angles[2])
    R12=np.array([[c12,s12,0],[-s12,c12,0],[0,0,1.]])
    R13=np.array([[c13,0,s13],[0,1,0],[-s13,0,c13]])
    R23=np.array([[1,0,0],[0,c23,s23],[0,-s23,c23]])
    U=(R23@R13@R12).astype(complex)@np.diag(np.exp(1j*np.array([.2,-.3,.7])))
    masses_eV=np.array([.008,.012,.052])
    m=U.conj()@np.diag(masses_eV/1e9)@U.conj().T
    Ye=np.diag([3e-6,6e-4,.01]).astype(complex)
    ans,mat=neutrino_observables(m,Ye)
    np.testing.assert_allclose(ans['masses_eV'],masses_eV,rtol=1e-9)
    np.testing.assert_allclose([ans['sin2_theta12'],ans['sin2_theta13'],ans['sin2_theta23']],
                                np.sin(angles)**2,atol=1e-8)
    assert ans['mixing_defined']
    np.testing.assert_allclose(mat['Takagi_U'].T @ m @ mat['Takagi_U'],
                                np.diag(masses_eV/1e9),atol=1e-19)


def test_near_degenerate_mass_does_not_invent_mixing():
    m=np.diag(np.array([.01,.01,.03])/1e9)
    vals,mats=neutrino_observables(m,np.diag([3e-6,.0006,.01]))
    assert not vals['mixing_defined']
    assert vals['theta12_deg'] is None and mats=={}


def test_canonical_A_B_adapter_roundtrip_and_synthetic_gate(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'fixture')
    with pytest.raises(AdapterBlocked,match='Synthetic input forbidden'):
        prepare_pair('A_-4',a,b,tmp_path/'adapted')
    pair=prepare_pair('A_-4',a,b,tmp_path/'adapted',allow_synthetic=True)
    pf=preflight_pair(pair.wilson,pair.background,permit_synthetic=True)
    assert pf['status']=='ready_conditional',pf['blockers']
    wm=json.loads(pair.wilson.read_text());bm=json.loads(pair.background.read_text())
    assert wm['model']==bm['model']=='A_-4'
    assert wm['scheme']==bm['scheme']=='MSbar'
    assert wm['units']==bm['units']=='GeV'
    with np.load(a.parent/'wilson_at_MS.npz',allow_pickle=False) as original, \
         np.load(pair.wilson.parent/'wilson_at_MS.npz',allow_pickle=False) as adapted:
        np.testing.assert_array_equal(original['C_psi2phi2_ms'],adapted['C_psi2phi2_ms'])


def test_legacy_A_full_C_is_adapted_by_field_labels_not_array_rename(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'fixture')
    _make_legacy_A(a)
    pair=prepare_pair('A_-4',a,b,tmp_path/'bridge',allow_synthetic=True)
    pf=preflight_pair(pair.wilson,pair.background,permit_synthetic=True)
    assert pf['status']=='ready_conditional',pf['blockers']
    new=json.loads(pair.wilson.read_text())
    assert new['model']=='A_-4' and new['units']=='GeV' and new['scheme']=='MSbar'
    assert new['lepton_indices']==[[0,1],[2,3],[4,5]]
    assert new['canonical_wilson_convention'].startswith('C_ijab_full_sum')


def test_B_native_trajectories_extract_genuine_MS_endpoint(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'fixture')
    bm=json.loads(b.read_text())
    with np.load(b.parent/'background.npz',allow_pickle=False) as f:
        arrays={k:f[k].copy() for k in f.files}
    source={k:bm[k] for k in (
        'model','matching_context_id','matched_seed_id','lagrangian_source_id',
        'mu_f_gev','mu_s_gev','scheme','units','synthetic_fixture',
        'quartic_potential_convention','quartic_inventory_status','wilson_scalar_basis_id',
        'background_scalar_basis_id','scalar_background_real_labels',
        'single_common_scalar_threshold','scalar_mass_origin','quartic_channel_names')}
    source.update(schema='t3-parallel-B-background-v1',
                  quartic_scalar_order='H_RI_interleaved_then_heavy_background_v1')
    for key in ('gY','g2','g3','lambdaH'):
        v=bm['sm_couplings_at_ms'][key]
        arrays['eft1__'+key]=np.array([v*.99,v],complex)
    for key in ('ye','yu','yd'):
        v=arrays.pop(key)
        arrays['eft1__'+key]=np.stack([v*.98,v])
    arrays['mu_eft1_gev']=np.array([100000.,1000.])
    np.savez_compressed(b.parent/'background.npz',**arrays)
    _write_json(b,source)
    pair=prepare_pair('A_-4',a,b,tmp_path/'bridge',allow_synthetic=True)
    with np.load(pair.background.parent/'background.npz',allow_pickle=False) as z:
        np.testing.assert_allclose(z['ye'],arrays['eft1__ye'][-1])
    pf=preflight_pair(pair.wilson,pair.background,permit_synthetic=True)
    assert pf['status']=='ready_conditional',pf['blockers']


def test_native_beta_trajectories_without_quartic_component_tensor_block(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'fixture')
    with np.load(b.parent/'background.npz',allow_pickle=False) as f:
        arrays={k:f[k] for k in f.files if k!='quartic_real_full_ms'}
    np.savez_compressed(b.parent/'background.npz',**arrays)
    # Recompute checksum on this artificially damaged fixture to expose the
    # real missing component rather than just a hash failure.
    import hashlib
    info=json.loads(b.read_text())
    info['tensor_npz_sha256']=hashlib.sha256((b.parent/'background.npz').read_bytes()).hexdigest()
    _write_json(b,info)
    with pytest.raises(AdapterBlocked,match='full scalar quartic'):
        prepare_pair('A_-4',a,b,tmp_path/'bridge',allow_synthetic=True)


def test_different_source_state_hash_blocks_no_spurious_results(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'fixture')
    meta=json.loads(b.read_text());meta['matched_seed_id']='different'
    _write_json(b,meta)
    with pytest.raises(AdapterBlocked,match='matched_seed_id'):
        prepare_pair('A_-4',a,b,tmp_path/'bridge',allow_synthetic=True)
    assert not (tmp_path/'bridge'/'wilson_at_MS.npz').exists()


def test_missing_16_model_data_are_blockers_not_filled(tmp_path):
    script=ROOT/'scripts'/'RunT3RoughAll16.py'
    out=tmp_path/'results'
    p=subprocess.run([sys.executable,str(script),
                      '--background-root',str(tmp_path/'noB'),
                      '--wilson-root',str(tmp_path/'noA'),
                      '--output',str(out)],cwd=ROOT,text=True,capture_output=True)
    assert p.returncode==1,(p.stdout,p.stderr)
    summary=json.loads((out/'summary.json').read_text())
    assert summary['counts']['requested']==16
    assert summary['counts']['blocked']==16
    assert summary['counts']['successful_real']==0
    assert all(x['masses_eV'] is None for x in summary['models'] if 'masses_eV' in x)
    assert not list(out.rglob('C5_smeft_low.npz'))


def test_synthetic_smoke_A4_one_common_scalar_RGE(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'inputs')
    pair=prepare_pair('A_-4',a,b,tmp_path/'bridge',allow_synthetic=True)
    res=run_rough_case('A_-4',pair,tmp_path/'results',allow_synthetic=True)
    assert res['status']=='SyntheticFixtureSmokeTest'
    assert len(res['masses_eV'])==3
    assert res['scalar_eigenstates_integrated_out']==n_real('A_-4')
    for name in ('C5_smeft_at_MS.npz','C5_smeft_low.npz',
                 'scalar_mass_basis.npz','observables.json','manifest.json'):
        assert (tmp_path/'results'/name).is_file()
    mani=json.loads((tmp_path/'results'/'manifest.json').read_text())
    assert mani['rge_bookkeeping']['old_LLSS_to_C5_direct_leading_log_added_again'] is False
    assert mani['SMEFT_running']['performed']
    assert mani['benchmark_approximation']['other_independently_unavailable_finite_terms_assumed_zero']


@pytest.mark.parametrize('model',ALL_MODELS)
def test_all_16_synthetic_scalar_mass_basis_and_pair_bubble(tmp_path,model):
    a,b=make_fixture(model,tmp_path/'source')
    pair=prepare_pair(model,a,b,tmp_path/'bridge',allow_synthetic=True)
    report=preflight_pair(pair.wilson,pair.background,permit_synthetic=True,expected_model=model)
    assert report['status']=='ready_conditional',report['blockers']
    matching=report['checks']['single_MS_finite_bubble_match_candidate']['details']
    assert matching['eigenstates_integrated_out_at_single_MS']==n_real(model)
    assert matching['unordered_eigenstate_pairs']==n_real(model)*(n_real(model)+1)//2


@pytest.mark.parametrize('model', ALL_MODELS)
def test_missing_labels_recover_only_from_declared_model_basis(tmp_path,model):
    """Producer fields are actually in the fixture's declared canonical/native frame."""
    a,b=make_fixture(model,tmp_path/'fixture')
    am=json.loads(a.read_text());bm=json.loads(b.read_text())
    expected_a=am.pop('scalar_source_real_labels')
    expected_b=bm.pop('scalar_background_real_labels')
    _write_json(a,am);_write_json(b,bm)
    pair=prepare_pair(model,a,b,tmp_path/'adapted',allow_synthetic=True)
    wa=json.loads(pair.wilson.read_text());ba=json.loads(pair.background.read_text())
    assert wa['scalar_source_real_labels']==expected_a
    assert ba['scalar_background_real_labels']==expected_b
    assert wa['scalar_labels_origin']=='declared_standard_A_basis_id'
    assert ba['scalar_labels_origin']=='declared_standard_B_basis_id'
    assert wa['model']==ba['model']==model
    pf=preflight_pair(pair.wilson,pair.background,permit_synthetic=True)
    assert pf['status']=='ready_conditional',pf['blockers']


def test_missing_A_labels_without_documented_basis_fails_closed(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'fixture')
    meta=json.loads(a.read_text())
    meta.pop('scalar_source_real_labels')
    meta['scalar_basis_id']='unmapped_native_S_field_order'
    other=json.loads(b.read_text())
    other['wilson_scalar_basis_id']=meta['scalar_basis_id']
    _write_json(a,meta);_write_json(b,other)
    with pytest.raises(AdapterBlocked,match='scalar_source_real_labels'):
        prepare_pair('A_-4',a,b,tmp_path/'adapted',allow_synthetic=True)
    assert not (tmp_path/'adapted'/'wilson.json').exists()


def test_missing_B_labels_without_documented_basis_fails_closed(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'fixture')
    meta=json.loads(b.read_text())
    meta.pop('scalar_background_real_labels')
    meta['background_scalar_basis_id']='unknown_B_native_components'
    _write_json(b,meta)
    with pytest.raises(AdapterBlocked,match='scalar_background_real_labels'):
        prepare_pair('A_-4',a,b,tmp_path/'adapted',allow_synthetic=True)


def test_reference_Weyl45_requires_correct_scalar_and_lepton_axes(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'fixture')
    meta=json.loads(a.read_text())
    meta.pop('scalar_source_real_labels')
    meta['scalar_basis_id']='T3_RGBeta_reference_Weyl45_real_interleaved_v1'
    other=json.loads(b.read_text())
    other['wilson_scalar_basis_id']=meta['scalar_basis_id']
    _write_json(a,meta);_write_json(b,other)
    # Fixture has only 8 Weyl components; it may not claim the reference 45-Weyl frame.
    with pytest.raises(AdapterBlocked,match='Weyl45 basis ID conflicts'):
        prepare_pair('A_-4',a,b,tmp_path/'adapted',allow_synthetic=True)


def test_B_permuted_full_quartic_is_reindexed_numerically(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'fixture')
    bm=json.loads(b.read_text())
    with np.load(b.parent/'background.npz',allow_pickle=False) as z:
        arrays={k:z[k].copy() for k in z.files}
    original=arrays['quartic_real_full_ms'].copy()
    n=original.shape[0]
    perm=[4,0,5,1,2,3]+list(range(6,n))
    for ax in range(4):
        arrays['quartic_real_full_ms']=np.take(arrays['quartic_real_full_ms'],perm,axis=ax)
    if 'quartic_real_channels_ms' in arrays:
        for ax in range(1,5):
            arrays['quartic_real_channels_ms']=np.take(arrays['quartic_real_channels_ms'],perm,axis=ax)
    np.savez_compressed(b.parent/'background.npz',**arrays)
    bm.pop('format')
    bm.pop('tensor_npz_sha256')
    bm['schema']='t3-parallel-B-background-v1'
    bm['scalar_full_real_labels']=[(['H[1,R]','H[1,I]','H[2,R]','H[2,I]']+
                                   bm['scalar_background_real_labels'])[i] for i in perm]
    _write_json(b,bm)
    pair=prepare_pair('A_-4',a,b,tmp_path/'adapted',allow_synthetic=True)
    with np.load(pair.background.parent/'background.npz',allow_pickle=False) as z:
        np.testing.assert_array_equal(z['quartic_real_full_ms'],original)
    pf=preflight_pair(pair.wilson,pair.background,permit_synthetic=True)
    assert pf['status']=='ready_conditional',pf['blockers']


def test_reconverted_canonical_pair_preserves_provenance(tmp_path):
    a,b=make_fixture('A_-4',tmp_path/'fixture')
    first=prepare_pair('A_-4',a,b,tmp_path/'first',allow_synthetic=True)
    second=prepare_pair('A_-4',first.wilson,first.background,tmp_path/'second',allow_synthetic=True)
    for file in ('wilson.json','background.json'):
        x=json.loads((tmp_path/'first'/file).read_text())
        y=json.loads((tmp_path/'second'/file).read_text())
        for key in ('model','matching_context_id','matched_seed_id',
                    'lagrangian_source_id','mu_f_gev','mu_s_gev','scheme','units'):
            assert x[key]==y[key]
    pf=preflight_pair(second.wilson,second.background,permit_synthetic=True)
    assert pf['status']=='ready_conditional',pf['blockers']
