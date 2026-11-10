"""Batch 6 scalar-pole/CT algebra, full source coverage and scan arithmetic.

Designed to run against the installed repository. Does not require re-running
Matchete, UV RGBeta, or any previous numerical batch during pytest.
"""
from __future__ import annotations

import math
import numpy as np
import pytest

from Numerical.parallel_t3.wilson.msbar_counterterm_batch6 import (
    LOOP, ADM_CHANNELS, SCENARIOS, SCALE_FACTORS,
    analytic_ct_from_soft, candidate_adm_counterterms,
    matching_scale_scenarios, native_pole_by_quartic,
    relative, source_scalar_counterterm, symmetric3)


def sym3(x):
    return (x+x.T)/2


def synthetic_counterterm_inputs():
    rng=np.random.default_rng(18)
    n=8  # both heavy multiplets, not one scalar threshold
    c=(rng.normal(size=(3,3,3,2,2,n,n)) +
       1j*rng.normal(size=(3,3,3,2,2,n,n)))*1e-12
    c=(c+c.transpose(0,2,1,4,3,5,6))/2
    c=(c+c.swapaxes(-1,-2))/2
    q=rng.normal(size=(4,4,4,n,n))*0.002
    q=(q+q.swapaxes(1,2))/2
    q=(q+q.swapaxes(-1,-2))/2
    return c,q,('quarticH1','quarticT3','quarticX','quarticY')


def test_explicit_eft1_ct_cancels_calculated_soft_pole():
    from Numerical.parallel_t3.threshold.matching import master_weinberg_template
    c,q,names=synthetic_counterterm_inputs()
    d=source_scalar_counterterm(c,q,names)
    assert d['scalar_components_active_at_MF']==8
    assert d['fermions_removed_together_at_MF']==3
    assert d['CT_per_channel_per_heavy'].shape==(4,3,3,3)
    assert np.linalg.norm(d['EFT1_soft_bare_pole'])>0
    assert np.allclose(d['EFT1_soft_bare_pole']+
                       d['EFT1_MSbar_counterterm_pole'],0,atol=1e-23)
    assert relative(d['CT_per_channel'].sum(axis=0),d['CT_scalar_LLHH'])<1e-12
    assert relative(d['CT_per_heavy'].sum(axis=0),d['CT_scalar_LLHH'])<1e-12
    # Direct contraction independent of the implementation's channel sum.
    amplitude=np.einsum('kfgijab,xABab->fgijAB',c,q,optimize=True)
    T=master_weinberg_template()
    reference=np.einsum('fgijAB,ijAB->fg',amplitude,T.conj(),optimize=True)/np.vdot(T,T)/(2*LOOP)
    assert relative(d['CT_scalar_LLHH'],reference)<1e-12


def test_counterterm_rejects_missing_scalar_components_and_channels():
    c,q,names=synthetic_counterterm_inputs()
    with pytest.raises(ValueError):source_scalar_counterterm(c[:2],q,names)
    with pytest.raises(ValueError):source_scalar_counterterm(c,q[:2],names)
    with pytest.raises(ValueError):source_scalar_counterterm(c,q,(names[0],)*4)
    with pytest.raises(ValueError):source_scalar_counterterm(c,q[:,:,:,:,:6],names)


def test_beta_to_msbar_pole_conversion_uses_full_seven_channels():
    rng=np.random.default_rng(31)
    strict={}
    for name in ADM_CHANNELS:
        strict['C5_EFT1_strict_running_'+name]=sym3(
            rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))*1e-12
    strict['C5_EFT1_strict_running_delta']=sum(strict.values())
    L=math.log(1e3/1e5)
    z=candidate_adm_counterterms(strict,L)
    assert len(z['CT_per_ADM_channel'])==7
    assert np.allclose(z['CT_candidate_total_7channel'],
                       strict['C5_EFT1_strict_running_delta']/(2*L),atol=1e-25)
    assert np.array_equal(z['CT_candidate_scalar_pair'],
                          z['CT_per_ADM_channel']['scalar_pair'])


def test_adm_source_missing_or_tampered_is_not_silently_completed():
    strict={'C5_EFT1_strict_running_'+n:np.eye(3)*1e-12 for n in ADM_CHANNELS}
    strict['C5_EFT1_strict_running_delta']=sum(strict.values())
    with pytest.raises(ValueError):candidate_adm_counterterms(strict,0.1)
    del strict['C5_EFT1_strict_running_mixed_yukawa_gauge']
    with pytest.raises(ValueError):candidate_adm_counterterms(strict,-4.6)


def test_native_raw_pole_retains_original_quartic_names():
    k=np.zeros((4,3,3,3),complex)
    for i in range(4):
        k[i]=np.eye(3)[None,:,:]*(i+1)*1e-12
    names=['lambdaA','lambdaB','lambdaA','lambdaC']
    ans=native_pole_by_quartic(k,names)
    assert set(ans)=={'lambdaA','lambdaB','lambdaC'}
    assert np.allclose(ans['lambdaA'],4*3e-12*np.eye(3))
    assert np.allclose(sum(ans.values()),k.sum(axis=(0,1)))
    with pytest.raises(ValueError):native_pole_by_quartic(k,['other']*4)


def test_analytic_b0_matching_scale_derivative_returns_own_ct():
    ct=np.array([[4e-13+0.1e-13j,2e-13,0],
                 [2e-13,-1e-13,0],[0,0,3e-13]],complex)
    origin=np.eye(3,dtype=complex)*(-1e-11)
    soft={f:origin-2*math.log(f)*ct for f in SCALE_FACTORS}
    calculated=analytic_ct_from_soft(soft,SCALE_FACTORS)
    assert relative(calculated,ct)<1e-13
    soft[2.0]=soft[2.0]+np.eye(3)*0.5e-12
    with pytest.raises(ValueError):analytic_ct_from_soft(soft,SCALE_FACTORS)


def test_matching_scale_nominal_sources_and_three_scenarios():
    m=np.array([0.9,1.,1.1])*1e5
    K=np.zeros((2,3,3,3),complex)
    for j in range(2):
        for k in range(3):K[j,k]=(1+j+k)*1e-13*np.eye(3)
    nh=np.einsum('kfij,f->ij',K,1+2*np.log(1e5/m))
    ct=np.eye(3,dtype=complex)*1.4e-12
    soft={f:np.eye(3,dtype=complex)*2e-11-2*math.log(f)*ct
          for f in SCALE_FACTORS}
    analytic=np.eye(3,dtype=complex)*1.0e-12
    running=np.eye(3,dtype=complex)*2.5e-12
    bubble=np.eye(3,dtype=complex)*3.0e-12
    z=matching_scale_scenarios(
        native_constant_by_kernel_fermion=K,masses_gev=m,
        c5_analytic_hard_mf=analytic,c5_native_formal_mf=nh,
        c5_running_delta_mf_to_ms=running,c5_scalar_bubble_ms=bubble,
        soft_finite=soft)
    assert len(z['scales'])==3
    for f in SCALE_FACTORS:
        assert set(z['scales'][f])==set(SCENARIOS)
    assert np.allclose(z['scales'][1.]['analytic_triangle']['C5_SMEFT_MS'],
                       analytic+running+bubble)
    assert np.allclose(z['scales'][1.]['native_formal']['C5_SMEFT_MS'],
                       nh+running+bubble)
    assert np.allclose(z['scales'][1.]['native_as_UV_minus_soft_hypothesis']['C5_SMEFT_MS'],
                       nh-soft[1.]+running+bubble)
    assert relative(z['scales'][0.5]['analytic_triangle']['hard_boundary_MF_renormalization_scale'],
                    analytic+2*math.log(0.5)*ct)<1e-13


def test_no_extra_physical_threshold_or_native_subtraction_is_inferred():
    m=np.array([1e5,1.1e5,1.2e5]); k=np.zeros((1,3,3,3),complex)
    soft={0.5:np.eye(3)*1e-11, 1.:np.eye(3)*1e-11, 2.:np.eye(3)*1e-11}
    params=dict(native_constant_by_kernel_fermion=k,masses_gev=m,
                c5_analytic_hard_mf=np.zeros((3,3)),
                c5_native_formal_mf=np.zeros((3,3)),
                c5_running_delta_mf_to_ms=np.zeros((3,3)),
                c5_scalar_bubble_ms=np.zeros((3,3)),soft_finite=soft)
    with pytest.raises(ValueError):matching_scale_scenarios(**params,factors=(0.005,1.0))
    with pytest.raises(ValueError):matching_scale_scenarios(**params,factors=(0.5,0.5,1.0))
