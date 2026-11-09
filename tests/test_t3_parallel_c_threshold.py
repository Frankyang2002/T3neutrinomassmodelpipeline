"""Numerical smoke/contracts; synthetic fixtures never certify physical matching."""
from __future__ import annotations
import json
import math
import numpy as np
import pytest
from Numerical.parallel_t3.threshold.models import (ALL_MODELS,from_parameters,
                    MIXED,SELF_S1,SELF_S2,dimensions,n_real,scalar_labels)
from Numerical.parallel_t3.threshold.matching import (master_weinberg_template,
                    b0bar,project_master,bubble_amplitude,mass_eigenbasis)
from Numerical.parallel_t3.threshold.fixtures import make_fixture
from Numerical.parallel_t3.threshold.pipeline import run_case


def test_all16_registry():
    assert len(ALL_MODELS)==16
    assert len(set(ALL_MODELS))==16
    for model in ALL_MODELS:
        d1,d2,df=dimensions(model)
        assert n_real(model)==2*(d1+d2)
        assert len(scalar_labels(model))==n_real(model)
    assert len(MIXED)==2 and len(SELF_S1|SELF_S2)==6


def test_weinberg_master_pivot_quarter_and_flavour_complex():
    t=master_weinberg_template()
    assert t.shape==(2,2,4,4)
    assert t[0,0,2,2]==pytest.approx(.25)
    c5=np.array([[1+2j,.2j,0],[.2j,2-.7j,0],[0,0,.3+1j]])*1e-14
    amp=np.einsum('fg,ijAB->fgijAB',c5,t)
    val,err=project_master(amp)
    assert err<1e-13
    np.testing.assert_allclose(val,c5,rtol=1e-13,atol=1e-28)


def test_b0_degenerate_and_nearly_degenerate():
    x=1.69e6
    assert b0bar(x,x,1000.)==pytest.approx(-math.log(1.69))
    near=b0bar(x*(1+1e-8),x*(1-1e-8),1000.)
    assert near==pytest.approx(-math.log(1.69),rel=1e-13)
    with pytest.raises(ValueError):b0bar(-1,x,1000)


def test_mixed_complex_phase_potential_map():
    for model in MIXED:
        ob=from_parameters(model,1.9e6,2.6e6,1.4e5+1.1e5j)
        R=ob.R_canonical_from_native
        assert R.shape==(8,8)
        np.testing.assert_allclose(R.T@R,np.eye(8),atol=1e-13)
        eig,U=mass_eigenbasis(ob.mass2_canonical_gev2)
        assert abs(eig[0]-eig[3])<1e-6
        assert abs(eig[4]-eig[7])<1e-6
        # Physically equivalent quadratic energy in the native mass basis.
        rng=np.random.default_rng(19)
        for _ in range(3):
            z=rng.normal(size=8)
            np.testing.assert_allclose(z@(R.T@ob.mass2_canonical_gev2@R)@z,
                z@ob.mass2_native_gev2@z,rtol=1e-13)


def test_self_complex_phase_offdiagonal_and_split():
    for model in SELF_S1|SELF_S2:
        d1,d2,_=dimensions(model)
        mu=110000.+83000j
        ob=from_parameters(model,1.9e6,2.5e6,mu)
        s=0 if model in SELF_S1 else 2*d1
        n=d1 if model in SELF_S1 else d2
        for j in range(n):
            block=ob.mass2_canonical_gev2[s+2*j:s+2*j+2,s+2*j:s+2*j+2]
            assert block[0,1]==pytest.approx(-2*mu.imag)
            np.testing.assert_allclose(np.sort(np.linalg.eigvalsh(block)),
                                       [float(np.trace(block)/2-2*abs(mu)),
                                        float(np.trace(block)/2+2*abs(mu))],rtol=1e-12)
        np.testing.assert_allclose(ob.R_canonical_from_native.T@
                                  ob.R_canonical_from_native,np.eye(n_real(model)))


def test_missing_special_mu_not_faked():
    with pytest.raises(ValueError,match='explicitly supplied'):
        from_parameters('B_-1',1.9e6,2.0e6)
    with pytest.raises(ValueError,match='non-positive'):
        from_parameters('E_+0',1e5,2e6,1e5+1e5j)


def test_basis_invariance_for_arbitrary_nondegenerate_mass_and_cross_tensor():
    rng=np.random.default_rng(17)
    n=6
    A=rng.normal(size=(n,n));m=A.T@A+np.diag(np.arange(1,n+1))*1e6
    R=np.linalg.qr(rng.normal(size=(n,n)))[0]
    c=rng.normal(size=(3,3,2,2,n,n))+1j*rng.normal(size=(3,3,2,2,n,n))
    c=.5*(c+c.swapaxes(-2,-1))
    q=rng.normal(size=(4,4,n,n))
    q=.5*(q+q.swapaxes(0,1));q=.5*(q+q.swapaxes(-2,-1))
    before=bubble_amplitude(c,q,m,1200.,detailed=False)['amplitude']
    cm=np.einsum('fgijab,ak,bl->fgijkl',c,R,R,optimize=True)
    qm=np.einsum('ABab,ak,bl->ABkl',q,R,R,optimize=True)
    after=bubble_amplitude(cm,qm,R.T@m@R,1200.,detailed=False)['amplitude']
    np.testing.assert_allclose(before,after,rtol=2e-11,atol=1e-13)


@pytest.mark.parametrize('model',ALL_MODELS)
def test_all16_synthetic_strict_end_to_end_matching(tmp_path,model):
    w,b=make_fixture(model,tmp_path/'input')
    out=tmp_path/'results'
    doc=run_case(wilson_manifest=w,background_manifest=b,output_dir=out,
                 mode='strict_one_loop',match_only=True)
    assert doc['scalar_eigenstates_integrated_out']==n_real(model)
    assert doc['single_scalar_threshold'] is True
    assert doc['rge_bookkeeping']['old_LLSS_to_C5_direct_leading_log_added_again'] is False
    assert doc['normalization_and_status']['synthetic_fixture'] is True
    assert doc['normalization_and_status']['physical_neutrino_mass_certified'] is False
    assert doc['all_unordered_eigenstate_pairs_included']==n_real(model)*(n_real(model)+1)//2
    with np.load(out/'C5_smeft_at_MS.npz',allow_pickle=False) as z:
        np.testing.assert_allclose(z['C5_EFT1_MS']+z['C5_scalar_bubble_MS']+
                                   z['C5_independent_finite_MS'],z['C5_SMEFT_MS'])
        np.testing.assert_allclose(z['C5_SMEFT_MS'],z['C5_SMEFT_MS'].T)
    assert len(doc['quartic_channel_names'])==2
    assert doc['matching_residuals']['finite_scalar_bubble'] < 1e-8


def test_strict_versus_rge_improved_contributions(tmp_path):
    w,b=make_fixture('A_-4',tmp_path/'inputs')
    strict=run_case(wilson_manifest=w,background_manifest=b,
                    output_dir=tmp_path/'strict',mode='strict_one_loop',match_only=True)
    improved=run_case(wilson_manifest=w,background_manifest=b,
                      output_dir=tmp_path/'improved',mode='rge_improved',match_only=True)
    assert strict['rge_bookkeeping']['bubble_LLSS_source']=='tree_order0_at_MF'
    assert improved['rge_bookkeeping']['bubble_LLSS_source']=='evolved_full_at_MS_rge_improved'
    s=np.load(tmp_path/'strict'/'C5_smeft_at_MS.npz')
    r=np.load(tmp_path/'improved'/'C5_smeft_at_MS.npz')
    np.testing.assert_array_equal(s['C5_EFT1_MS'],r['C5_EFT1_MS'])
    assert np.linalg.norm(s['C5_scalar_bubble_MS']-r['C5_scalar_bubble_MS'])>1e-20


def test_stale_seed_and_duplicate_log_rejected(tmp_path):
    w,b=make_fixture('A_-4',tmp_path/'inputs')
    bm=json.loads(b.read_text());bm['lagrangian_source_id']='OLD_STALE_PORTAL'
    b.write_text(json.dumps(bm))
    with pytest.raises(ValueError,match='lagrangian_source_id'):
        run_case(wilson_manifest=w,background_manifest=b,output_dir=tmp_path/'bad',match_only=True)
    bm['lagrangian_source_id']='SYNTHETIC_NOT_ACTUAL_LAGRANGIAN'
    b.write_text(json.dumps(bm))
    wm=json.loads(w.read_text());wm['legacy_direct_llss_to_c5_log_added']=True
    w.write_text(json.dumps(wm))
    with pytest.raises(ValueError,match='duplicate|log'):
        run_case(wilson_manifest=w,background_manifest=b,output_dir=tmp_path/'bad',match_only=True)


def test_strict_needs_tree_A_tensor(tmp_path):
    w,b=make_fixture('B_-1',tmp_path/'inputs')
    path=w.parent/'wilson_at_MS.npz'
    with np.load(path,allow_pickle=False) as z:
        arrays={k:z[k] for k in z.files if k!='C_psi2phi2_tree_mf'}
    np.savez_compressed(path,**arrays)
    wm=json.loads(w.read_text())
    from hashlib import sha256
    wm['tensor_npz_sha256']=sha256(path.read_bytes()).hexdigest()
    w.write_text(json.dumps(wm))
    with pytest.raises(ValueError,match='order-zero'):
        run_case(wilson_manifest=w,background_manifest=b,
                 output_dir=tmp_path/'bad',mode='strict_one_loop',match_only=True)


def test_reuses_existing_smeft_runner_when_repository_is_present(tmp_path):
    """Runs genuine legacy module in user's repo; skipped in isolated ZIP tests."""
    try:
        from Numerical.running.SMWeinbergEvolution import SMWeinbergInitialConditions
        assert SMWeinbergInitialConditions is not None
    except ImportError:
        pytest.skip('Repository SMWeinbergEvolution is not installed in isolated ZIP test tree')
    w,b=make_fixture('A_-4',tmp_path/'input')
    doc=run_case(wilson_manifest=w,background_manifest=b,
                 output_dir=tmp_path/'results',mu_low=173.,mode='rge_improved')
    assert doc['SMEFT_running']['performed'] is True
    assert doc['SMEFT_running']['success'] is True
    assert doc['SMEFT_running']['nfev']>0
    assert np.isfinite(doc['mass_singular_values_eV']).all()


def test_optional_source_supplied_gauge_ward_mass_check():
    from Numerical.parallel_t3.threshold.matching import scalar_gauge_ward_residuals
    n=4
    q=np.zeros((4+n,)*4)
    M=np.diag([2.0e6,2.0e6,3.0e6,3.0e6])
    generator=np.zeros((1,4+n,4+n),complex)
    # Generator G=-i T in real coordinates. T real antisymmetric.
    generator[0,4,5]=-1j
    generator[0,5,4]=1j
    ward=scalar_gauge_ward_residuals(q,M,generator)
    assert ward['passed'] is True
    M[1,1]=2.1e6
    assert scalar_gauge_ward_residuals(q,M,generator)['passed'] is False
