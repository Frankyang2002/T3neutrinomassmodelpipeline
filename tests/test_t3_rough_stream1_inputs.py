"""Focused algebra, symmetry and source-fail-closed regression tests.

Tests use algebraic unit inputs only, not fake T3 production trajectories.
"""
from __future__ import annotations

import numpy as np
import pytest

from Numerical.parallel_t3.rough_inputs.builder import (
    assemble_quartic,quartic_unit_tensors,sm_yukawa_tensor,InputBlocked,
    _product,_bilinear,_scalar_mass
)
from Numerical.parallel_t3.wilson.basis import build_basis


def test_basic_quartic_derivative_and_normalisation():
    # V=lambdaH/2*(Hdag H)^2 -> V = lambdaH/8*(sum h_a^2)^2.
    units=quartic_unit_tensors('A_-4')
    assert units['lambdaH'].shape==(12,)*4
    assert np.isclose(units['lambdaH'][0,0,0,0],3)
    assert np.isclose(units['lambdaH'][0,0,1,1],1)
    assert np.isclose(units['lambdaH1'][0,0,4,4],1)
    assert np.max(abs(units['lambdaT3']))>0


def test_all16_quartic_geometry_and_pair_symmetries():
    from Numerical.parallel_t3.background.contracts import MODELS
    for model in MODELS:
        units=quartic_unit_tensors(model)
        state={key:np.asarray(0.01+0j) for key in units}
        if 'lambdaT3' in state:state['lambdaT3']=np.asarray(.02+.03j)
        arr,lost=assemble_quartic(units,state,reduced=False)
        assert not lost
        assert np.isfinite(arr).all()
        assert np.allclose(arr,arr.transpose(1,0,2,3))
        assert np.allclose(arr,arr.transpose(0,2,1,3))
        assert np.allclose(arr,arr.transpose(0,1,3,2))


def test_nonzero_unmapped_quartic_is_never_silent():
    units=quartic_unit_tensors('A_-4')
    state={'lambdaH':np.asarray(.5+0j),
           'lambdaT3GeneralQ999':np.asarray(.06+.02j)}
    with pytest.raises(InputBlocked,match='Unmapped nonzero'):
        assemble_quartic(units,state,reduced=False)
    arr,omitted=assemble_quartic(units,state,reduced=True)
    assert arr.shape==(12,)*4
    assert omitted['lambdaT3GeneralQ999']=={'real':.06,'imag':.02}


def test_full_45_weyl_sm_yukawa_embedding():
    b=build_basis('B_-1')
    yu=np.array([[.8,.1,.2],[.2,.5,.01],[.03,.05,1]],complex)
    yd=np.eye(3,dtype=complex)*.2
    ye=np.eye(3,dtype=complex)*.02
    Y=sm_yukawa_tensor(b,yu,yd,ye)
    assert Y.shape==(45,45,b.ns)
    assert np.allclose(Y,Y.swapaxes(0,1))
    f={x:i for i,x in enumerate(b.fermion_labels)}
    assert np.isclose(Y[f['Q[1,1,1]'],f['uc[2,1]'],2],yu[0,1]/np.sqrt(2))
    assert np.isclose(Y[f['L[1,1]'],f['ec[1]'],1],-1j*.02/np.sqrt(2))


def test_special_hessians_retain_real_imag_mixing():
    for model,key in [('A_-2','muT3GeneralS2Sq'),
                      ('D_+0','muT3GeneralS1Sq'),
                      ('B_-1','muT3GeneralS1S2'),
                      ('C_-1','muT3GeneralS1S2')]:
        obj=_scalar_mass(model,{'mS1Sq':np.asarray(4e6+0j),
                                'mS2Sq':np.asarray(3e6+0j),
                                key:np.asarray(12000+8000j)})
        M=obj.mass2_canonical_gev2
        assert M.shape[0]==M.shape[1]
        assert np.allclose(M,M.T)
        assert np.min(np.linalg.eigvalsh(M))>0
        assert np.max(abs(M-np.diag(np.diag(M))))>0


def test_run_one_model_exports_actual_contracts_from_saved_b_api(monkeypatch,tmp_path):
    """I/O probe only: constructed arrays are test fixtures, never model results."""
    from pathlib import Path
    from Numerical.parallel_t3.rough_inputs import builder as code
    from Numerical.parallel_t3.background.contracts import hash_file, tag

    model='A_-4'
    basis=build_basis(model)
    folder=tmp_path/'output/parallel_t3/background'/'B_all16_qualitative_v1'/tag(model)
    folder.mkdir(parents=True)
    for fname in ('background.npz','background_manifest.json',
                  'fermion_boundary_formal.npz','fermion_boundary_manifest.json'):
        (folder/fname).write_bytes(b'test-fixture-only')
    source=tmp_path/'Lagrangian/interactions/ScalarPotential.wl'
    source.parent.mkdir(parents=True)
    source.write_text('test-fixture-only',encoding='utf8')
    common=('gY','g2','g3','yu','yd','ye','mS1Sq','mS2Sq',
            'lambdaH','lambdaS1','lambdaS2','lambdaH1','lambdaH2',
            'lambda12','lambdaT3','lambdaH2Adj')
    bmeta={'model':model,'input_scales_gev':{'F':1e5,'S':1e3},
           'UV_endpoint_sha256':'test-uv-state',
           'beta_symbol_order_eft1':list(common),
           'beta_symbol_order_uv':list(common)+['MF','y1','y2'],
           'input_provenance':{'local_sources':{'Lagrangian/interactions/ScalarPotential.wl':{
                 'sha256':hash_file(source),'required_for_B_reuse':True}}}}
    base={'gY':.36,'g2':.65,'g3':1.1,'mS1Sq':4e6,'mS2Sq':3e6,
          'lambdaH':.3,'lambdaS1':.1,'lambdaS2':.1,
          'lambdaH1':.01,'lambdaH2':.02,'lambda12':.03,
          'lambdaT3':.01+.002j,'lambdaH2Adj':.001}
    arrays={'mu_uv_gev':np.array([1e6,2e5,1e5]),
            'mu_eft1_gev':np.array([1e5,1e4,1e3]),
            't_uv':np.log([1e6,2e5,1e5]),
            't_eft1':np.log([1e5,1e4,1e3])}
    for key in common:
        x=np.eye(3,dtype=np.complex128)*.1 if key in ('yu','yd','ye') else np.asarray(base[key],complex)
        arr=np.stack((x,x,x)).astype(np.complex128)
        arrays['eft1__'+key]=arr
        arrays['uv__'+key]=arr
    for key in ('MF','y1','y2'):
        arrays['uv__'+key]=np.repeat(np.eye(3)[None,:,:],3,axis=0).astype(complex)
    class Stub: pass
    Stub.manifest=bmeta
    Stub.arrays=arrays
    seed=tmp_path/'output/general_observable/fermion_first_hermitian_portal_v1/A_m4/data/eft1_after_F_wilson_seed.json'
    seed.parent.mkdir(parents=True,exist_ok=True)
    seed.write_text('test-only formal seed',encoding='utf8')
    fmeta={'model':model,'all_three_heavy_fermions_decoupled_at_single_MF':True,
           'same_run_uv_state_sha256':'test-uv-state',
           'source_signature_sha256':'test-formal-signature',
           'source_matchete_seed_sha256':hash_file(seed),
           'source_matchete_seed':'output/general_observable/fermion_first_hermitian_portal_v1/A_m4/data/eft1_after_F_wilson_seed.json',
           'mu_F_gev':1e5}
    formal=np.zeros((3,3,2,2,basis.ns-4,basis.ns-4),dtype=np.complex128)
    formal[0,0,0,0,0,0]=1e-5
    monkeypatch.setattr(code.Background,'load',lambda _folder:Stub())
    monkeypatch.setattr(code,'load_formal_boundary',lambda _folder:({'C_LLSS_formal':formal},fmeta))
    # The source used in this I/O test is a fixture, not a real A(-4) Matchete seed.
    opts=code.Options(root=tmp_path,assume_formal_identity=True,
        assume_sm_weyl_yukawa=True,hard_c5_zero_benchmark=True,
        reduced_potential=True)
    output=code.build_model(model,opts)
    assert output['status']=='written_provisional_A_C_inputs'
    cdir=tmp_path/'output/parallel_t3/rough_inputs/rough_fermion_first_v1/c/A_m4'
    adir=tmp_path/'output/parallel_t3/rough_inputs/rough_fermion_first_v1/a/A_m4'
    import json
    cm=json.loads((cdir/'background.json').read_text())
    assert cm['format']=='t3-parallel-c-background-input-v1'
    assert cm['single_common_scalar_threshold'] is True
    assert cm['quartic_inventory_status']=='full_tensor'
    assert cm['mu_f_gev']==1e5 and cm['mu_s_gev']==1e3
    with np.load(cdir/'background.npz') as cdata:
        assert set(('mass_squared_gev2_ms','R_wilson_from_background',
                    'quartic_real_full_ms','ye','yu','yd'))<=set(cdata.files)
        assert cdata['mass_squared_gev2_ms'].shape==(8,8)
        assert cdata['quartic_real_full_ms'].shape==(12,)*4
    with np.load(adir/'background.npz') as data:
        assert data['Y'].shape==(3,45,45,12)
        assert data['lambda4'].shape==(3,12,12,12,12)
    with np.load(adir/'boundary.npz') as data:
        assert data['C_full_provisional'].shape==(45,45,12,12)
    assert (adir/'saved_native_trajectory.npz').is_file()


def test_unknown_hard_c5_is_not_automatically_zero(tmp_path):
    from Numerical.parallel_t3.rough_inputs.builder import _c5_hard,Options
    meta={'model':'A_-4','UV_endpoint_sha256':'bench-uv'}
    with pytest.raises(InputBlocked,match='Hard C5 is unknown'):
        _c5_hard(tmp_path,Options(root=tmp_path),{},meta)
    hard,origin,_=_c5_hard(tmp_path,Options(root=tmp_path,hard_c5_zero_benchmark=True),{},meta)
    assert np.count_nonzero(hard)==0
    assert origin=='explicit_zero_benchmark_not_physical_vanishing'


def test_nonzero_hard_c5_external_same_run(tmp_path):
    import json
    from Numerical.parallel_t3.rough_inputs.builder import _c5_hard,Options
    from Numerical.parallel_t3.background.contracts import hash_file
    p=tmp_path/'input_hard/A_m4'
    p.mkdir(parents=True)
    c5=np.diag([1e-10,2e-10,3e-10]).astype(complex)
    np.savez_compressed(p/'hard_C5.npz',C5_hard=c5)
    (p/'hard_C5_manifest.json').write_text(json.dumps({
        'model':'A_-4','mu_F_gev':1e5,
        'same_run_uv_state_sha256':'same-uv',
        'hard_c5_file_sha256':hash_file(p/'hard_C5.npz')}),encoding='utf-8')
    got,origin,evidence=_c5_hard(tmp_path,Options(root=tmp_path,hard_c5_root=__import__('pathlib').Path('input_hard')),
                                  {},{'model':'A_-4','UV_endpoint_sha256':'same-uv'})
    assert np.array_equal(got,c5)
    assert origin=='supplied_independent_file'
    assert evidence['sha256']==hash_file(p/'hard_C5.npz')
