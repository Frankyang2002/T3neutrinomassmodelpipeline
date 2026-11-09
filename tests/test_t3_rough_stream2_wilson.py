"""Focused adapter, C contract, no-double-count and small ADM tests.

These tests contain synthetic arrays ONLY in pytest's temporary directory.
They must never be reported as genuine T3 saved-model results.
"""
from __future__ import annotations

from itertools import permutations
import json
from pathlib import Path
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from Numerical.parallel_t3.wilson import rough_stream2 as s
from Numerical.parallel_t3.wilson.basis import ALL_MODELS, build_basis
from Numerical.parallel_t3.wilson.coordinates import embed_c5, project_c5
from Numerical.parallel_t3.wilson.adm import beta, total_beta, SYMMETRIC_TERMS


def _test_input(tmp_path:Path, *, c5=False):
    root=tmp_path
    basis=build_basis('A_-4')
    folder=root/s.INPUT_ROOT/'A_m4'
    folder.mkdir(parents=True)
    f=basis.leptons[0][0]
    C=np.zeros((basis.nf,basis.nf,basis.ns,basis.ns),complex)
    C[f,f,4,4]=3e-10+1e-10j
    Y=np.zeros((2,basis.nf,basis.nf,basis.ns),complex)
    lam=np.zeros((2,)+(basis.ns,)*4,float)
    for pp in set(permutations((2,2,4,4))):
        lam[(slice(None),)+pp] = 0.015
    np.savez_compressed(folder/'background.npz',mu_gev=[1e5,1e3],Y=Y,
                        lambda4=lam,g2=[0.,0.],gY=[0.,0.],g3=[0.,0.])
    arrays={'C_psi2phi2_tree_mf':C}
    if c5:arrays['C5_EFT1_mf']=np.diag([1e-16,2e-16,3e-16]).astype(complex)
    np.savez_compressed(folder/'boundary.npz',**arrays)
    m={'model':'A_-4', 'basis_id':s.BASIS_ID,
       'fermion_labels':list(basis.fermion_labels),
       'scalar_labels':list(basis.scalar_labels),
       'full_SM_three_generation_Weyl_basis':True,
       'matching_context_id':'real-context-from-input-111',
       'matched_seed_id':'real-seed-from-input-222',
       'lagrangian_source_id':'real-source-from-input-333',
       'scalar_basis_id':'original-canonical-scalar-coordinate-id',
       'mu_f_gev':1e5,'mu_s_gev':1e3,
       'scheme':'MSbar',
       'operator_basis':'psi2phi2_no_factorial',
       'C_psi2phi2_tree_mf_is_canonical':True,
       'quartic_potential_convention':s.Q_CONVENTION}
    (folder/'background.json').write_text(json.dumps(m))
    (folder/'boundary.json').write_text(json.dumps(m))
    return root,folder,basis


def test_all16_slug_and_full_shapes():
    assert len(ALL_MODELS)==16
    assert s.slug('A_-4')=='A_m4'
    assert s.slug('E_+2')=='E_p2'
    for model in ALL_MODELS:
        b=build_basis(model)
        assert b.nf==45
        assert b.ns==4+2*b.dims[0]+2*b.dims[1]


def test_explicit_zero_boundary_and_full_running_background(tmp_path):
    root,folder,basis=_test_input(tmp_path)
    with pytest.raises(s.BlockedInput,match='No C5 boundary supplied'):
        s.load_canonical(root,folder,'A_-4')
    loaded=s.load_canonical(root,folder,'A_-4',c5_boundary='zero')
    assert loaded.background.y.shape==(2,45,45,basis.ns)
    assert loaded.background.quartic.shape==(2,)+(basis.ns,)*4
    assert loaded.source['boundary']['initial_C5_source']=='explicit_zero_boundary_requested'
    np.testing.assert_allclose(project_c5(loaded.initial,basis)[0],0)
    assert loaded.source['boundary']['LLSS_at_MF_frobenius']>0
    assert loaded.context['matching_context_id']=='real-context-from-input-111'


def test_hard_c5_polynomial_embedding_preserved(tmp_path):
    root,folder,basis=_test_input(tmp_path,c5=True)
    case=s.load_canonical(root,folder,'A_-4')
    c5,_=project_c5(case.initial,basis)
    np.testing.assert_allclose(c5,np.diag([1e-16,2e-16,3e-16]),rtol=1e-12,atol=1e-24)
    assert np.count_nonzero(case.initial[:,:,:4,:4])>3
    with pytest.raises(s.BlockedInput,match='zero C5 conflicts'):
        s.load_canonical(root,folder,'A_-4',c5_boundary='zero')


def test_reject_missing_or_incompatible_context(tmp_path):
    root,folder,_=_test_input(tmp_path)
    p=folder/'boundary.json'
    meta=json.loads(p.read_text())
    meta['matched_seed_id']='different-seed'
    p.write_text(json.dumps(meta))
    with pytest.raises(s.BlockedInput,match='Incompatible Stream1 matched_seed_id'):
        s.load_canonical(root,folder,'A_-4',c5_boundary='zero')
    meta['matched_seed_id']='real-seed-from-input-222'
    meta.pop('matching_context_id')
    p.write_text(json.dumps(meta))
    # Matching context may be sourced from the other Stream1 manifest.
    assert s.load_canonical(root,folder,'A_-4',c5_boundary='zero').context[
        'matching_context_id']=='real-context-from-input-111'
    bg=folder/'background.json'; bmeta=json.loads(bg.read_text())
    bmeta.pop('matching_context_id');bg.write_text(json.dumps(bmeta))
    with pytest.raises(s.BlockedInput,match='Missing Stream1 matching_context_id'):
        s.load_canonical(root,folder,'A_-4',c5_boundary='zero')


def test_reject_native_b_without_canonical_dictionary(tmp_path):
    root,folder,basis=_test_input(tmp_path)
    with np.load(folder/'background.npz') as v:
        y=np.array(v['Y'])
    np.savez_compressed(folder/'background.npz',mu_eft1_gev=[1e5,1e3],
                        eft1__ye=np.zeros((2,3,3),complex))
    with pytest.raises(s.BlockedInput,match='Native Workstream B RGBeta'):
        s.load_canonical(root,folder,'A_-4',c5_boundary='zero')


def test_source_digest_fail_closed(tmp_path):
    root,folder,_=_test_input(tmp_path)
    m=json.loads((folder/'background.json').read_text())
    m['npz_sha256']='0'*64
    (folder/'background.json').write_text(json.dumps(m))
    with pytest.raises(s.BlockedInput,match='npz_sha256 mismatch'):
        s.load_canonical(root,folder,'A_-4',c5_boundary='zero')


def test_channel_mixing_with_existing_six_term_engine():
    """Small 2-Weyl, 2-real-scalar candidate ADM test, NOT physical T3."""
    c=np.zeros((2,2,2,2),complex)
    c[0,0,0,0]=3e-10+1e-10j
    y=np.zeros((2,2,2),complex)
    quartic=np.zeros((2,2,2,2),float)
    for pp in set(permutations((0,0,1,1))):
        quartic[pp]=0.02
    terms,tag=beta(c,y,quartic,variant='fonseca_eq22_six_term')
    assert set(terms)==set(SYMMETRIC_TERMS)
    assert abs(total_beta(terms)[0,0,1,1])>0
    assert abs(terms['scalar_pair'][0,0,1,1])>0

    def f(t,v):
        mat=v.reshape(2,2,2,2)
        z,_=beta(mat,y,quartic,variant='fonseca_eq22_six_term')
        return (s.LOOP*total_beta(z)).reshape(-1)
    sol=solve_ivp(f,(math.log(1e5),math.log(1e3)),c.reshape(-1),
                  method='DOP853',rtol=1e-8,atol=1e-22)
    assert sol.success
    generated=sol.y[:,-1].reshape(2,2,2,2)[0,0,1,1]
    assert abs(generated)>0
    assert not np.allclose(generated,0,atol=0)


def test_export_exact_c_contract_and_no_double_log(tmp_path):
    root,folder,basis=_test_input(tmp_path)
    case=s.load_canonical(root,folder,'A_-4',c5_boundary='zero')
    initial=case.initial
    # Mock solver endpoint checks the Stream2 adapter only: the GENERATED
    # value below is not physics and remains solely in tmp_path.
    end=initial+embed_c5(np.array([[2e-16+1e-17j,1e-16,0],
                                   [1e-16,3e-16,0],[0,0,4e-16]],complex),basis)
    zero=np.zeros_like(initial)
    chs={name:zero for name in SYMMETRIC_TERMS}
    result={'C_initial':initial,'C_final':end,
            'C5_final':project_c5(end,basis)[0],
            'beta_at_MF':zero,'beta_at_MS':zero,
            'beta_channels_at_MF':chs,'beta_channels_at_MS':chs,
            'metadata':{'rhs_evaluations':9,'execution_backends':['test_mock_not_physics'],
                'projections':[{'mu_gev':1e5,'C5_frobenius_gev_inverse':0.,
                                'C_full_frobenius_gev_inverse':float(np.linalg.norm(initial))},
                               {'mu_gev':1e3,'C5_frobenius_gev_inverse':float(
                                    np.linalg.norm(project_c5(end,basis)[0])),
                                'C_full_frobenius_gev_inverse':float(np.linalg.norm(end))}]}}
    out=root/s.OUTPUT_ROOT/'A_m4'
    m=s.write_c_endpoint(case,result,out,adm_variant='fonseca_eq22_six_term',
                         rtol=1e-6,atol=1e-21,steps=1)
    assert m['format']==s.W_FORMAT
    assert m['matching_context_id']==case.context['matching_context_id']
    assert m['legacy_direct_llss_to_c5_log_added'] is False
    assert m['strict_1loop_endpoint_at_ms'] is False
    assert m['C5_generated_from_zero_boundary'] is True
    assert m['lepton_indices']==[list(p) for p in basis.leptons]
    assert m['higgs_source_indices']==[0,1,2,3]
    assert m['scalar_basis_id']=='original-canonical-scalar-coordinate-id'
    with np.load(out/'wilson_at_MS.npz',allow_pickle=False) as arrays:
        assert arrays['C_psi2phi2_ms'].shape==(45,45,basis.ns,basis.ns)
        assert arrays['C_psi2phi2_tree_mf'].shape==(45,45,basis.ns,basis.ns)
        c5=arrays['C5_EFT1_ms']
        np.testing.assert_allclose(c5,project_c5(arrays['C_psi2phi2_ms'],basis)[0])
        assert 'C_psi2phi2_ms_strict_1loop' not in arrays
    assert s.sha256(out/'wilson_at_MS.npz')==m['tensor_npz_sha256']


def test_saved_B_native_probe_reports_noncanonical(tmp_path):
    root=tmp_path
    folder=root/'output/parallel_t3/background/B_all16_qualitative_v1/A_m4'
    folder.mkdir(parents=True)
    uv='a'*64
    (folder/'background_manifest.json').write_text(json.dumps({
        'model':'A_-4','UV_endpoint_sha256':uv,
        'scalar_basis':'native RGBeta invariant coordinates'}))
    boundary=np.ones((3,3,2,2,8,8),complex)*1e-10
    np.savez_compressed(folder/'fermion_boundary_formal.npz',C_LLSS_formal=boundary)
    (folder/'fermion_boundary_manifest.json').write_text(json.dumps({
        'model':'A_-4','same_run_uv_state_sha256':uv,
        'fermion_matching_source_basis':'Matchete GammaCC source',
        'formal_LLsS_tree_available':True,'hard_C5_status':'hard_C5_unavailable',
        'fermion_boundary_npz_sha256':s.sha256(folder/'fermion_boundary_formal.npz')}))
    probe=s.probe_saved_b(root,Path('output/parallel_t3/background/B_all16_qualitative_v1'),
                          'A_-4')
    assert probe['nonzero_saved_LLSS'] is True
    assert probe['canonical_Y_lambda4_available'] is False
    assert probe['status']=='saved_B_native_requires_Stream1_canonical_mapping'


@pytest.mark.skipif(__import__('os').environ.get('T3_STREAM2_RUN_FULL45')!='1',
                    reason='Opt-in full 45-Weyl DOP853 test; can be slow with native engine')
def test_full45_DOP853_same_pipeline_synthetic_smoke(tmp_path):
    """Full-dimensional numerical test, but synthetic input, NOT a real saved model."""
    root,folder,basis=_test_input(tmp_path)
    result=s.run_one(root,'A_-4',c5_boundary='zero',steps=2,
                     rtol=1e-6,atol=1e-21)
    assert result['solver']['method']=='DOP853'
    assert result['C5_EFT1_mf_norm_gev_inverse']==0
    assert result['LLSS_mf_norm_gev_inverse']>0
    assert result['C5_generated_from_LLSS_only_boundary'] is True
    assert result['C5_EFT1_ms_norm_gev_inverse']>0
