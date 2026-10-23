"""Run inside the full T3 repository with: pytest -q tests/test_t3b_common_scale_threshold_logs.py"""
import numpy as np
import pytest
from scipy.linalg import expm

from Numerical.diagnostics.T3BHeavyLogAlignment import (
    flavor_kernel_by_heavy, matrix_to_json,
)
from Numerical.diagnostics.T3BExploratoryCoupledBenchmark import (
    FWD, REV, LOOP, calculate, common_scale_log_boundary,
    coupled_heavy_resolved, fixed_order_log_from_common_boundary,
)


def samples(*, lam=0.3+0.05j):
    zero={'real':np.zeros((3,3)).tolist(), 'imag':np.zeros((3,3)).tolist()}
    return [dict(mu_gev=m, g2=0.60, lambdaH=0.15,lambda12=0.08,
                 lambda12Adj=0.03,lambdaT3=[lam.real,lam.imag],
                 ye=zero,yu=zero,yd=zero)
            for m in (1e5,1e4,1e3)]


def inputs():
    y1=np.diag([0.1,0.2,0.3]).astype(complex)
    y2=np.diag([0.4,0.5,0.6]).astype(complex)
    return y1,y2,np.array([9e4,1.1e5,1.3e5])


def test_fixed_order_mass_resolved_log_identity():
    s=samples();y1,y2,m=inputs()
    observed=fixed_order_log_from_common_boundary(s,y1,y2,m)
    kernel=flavor_kernel_by_heavy(y1,y2,m)
    lam=complex(*s[-1]['lambdaT3'])
    expected=2*FWD*lam/LOOP * np.einsum('r,rij->ij',np.log(1e3/m),kernel)
    np.testing.assert_allclose(observed,expected,rtol=5e-14,atol=1e-24)


def test_exact_coupled_constant_background_with_heavy_log_boundary():
    s=samples(); y1,y2,m=inputs()
    a,b,details=common_scale_log_boundary(s,y1,y2,m,diagonal=False)
    kernel=flavor_kernel_by_heavy(y1,y2,m)
    lam=complex(*s[-1]['lambdaT3'])
    a0=kernel.sum(axis=0)
    b0=FWD*lam/LOOP * np.einsum('r,rij->ij', np.log(1e5/m),kernel)
    t=np.log(1e3/1e5)
    op=np.array([[0,REV*np.conjugate(lam)],[FWD*lam,0]],complex)
    propagator=expm(op*t/LOOP)
    np.testing.assert_allclose(a,propagator[0,0]*a0+propagator[0,1]*b0,
                               rtol=3e-12,atol=1e-20)
    np.testing.assert_allclose(b,propagator[1,0]*a0+propagator[1,1]*b0,
                               rtol=3e-12,atol=1e-20)
    assert details['sequential_partial_decoupling_implemented'] is False
    np.testing.assert_allclose(details['mass_resolved_log_mf_over_mr'],np.log(1e5/m))
    np.testing.assert_allclose(details['evolution_interval_gev'],[1e5,1e3])


def test_degenerate_heavy_masses_recover_previous_when_no_diagonal_running():
    s=samples(lam=0.1+0j)
    y1=np.eye(3)*0.1;y2=np.eye(3)*0.2
    m=np.full(3,1e5)
    a,b,d=common_scale_log_boundary(s,y1,y2,m,diagonal=False)
    old_a,old_b=coupled_heavy_resolved(s,y1,y2,m,diagonal=False)
    np.testing.assert_allclose(a,old_a,rtol=1e-10,atol=1e-19)
    np.testing.assert_allclose(b,old_b,rtol=1e-10,atol=1e-19)
    assert np.linalg.norm(d['heavy_boundary_c5_ordered']['real'])==0


def test_input_validation_for_masses_below_ms():
    s=samples();y1,y2,m=inputs()
    with pytest.raises(ValueError,match='exceed MS'):
        common_scale_log_boundary(s,y1,y2,[999,1e5,1.1e5])


def test_fixed_order_anchor_and_optional_scheme_comparison():
    s=samples();y1,y2,m=inputs()
    direct=fixed_order_log_from_common_boundary(s,y1,y2,m)
    baseline=direct+np.eye(3)*1e-10
    control={'status':'FixedOrderControlComplete',
             'intermediate_running':{'samples':s},
             'mu_f_gev':1e5,'mu_s_gev':1e3,
             'fixed_order_c5_at_ms':matrix_to_json(baseline)}
    normalization={'status':'UnfittedNormalizationPassed',
                   'production_direct_c5':matrix_to_json(direct),
                   'predicted_direct_c5':matrix_to_json(direct)}
    out=calculate(control,normalization,y1,y2,m,
                  compare_legacy_threshold=True)
    assert out['status']=='ExploratoryCoupledMatrixComplete'
    assert out['threshold_scheme']=='common_scale_log_boundary'
    assert out['common_boundary_fixed_order_direct_residual']<1e-12
    assert out['heavy_background_above_mf']=='not used'
    assert out['legacy_threshold_comparison'] is not None
    assert not out['validated_physical_prediction_available']
    norm2=dict(normalization)
    norm2['production_direct_c5']=matrix_to_json(direct*1.03)
    norm2['predicted_direct_c5']=matrix_to_json(direct*1.03)
    with pytest.raises(ValueError,match='does not reproduce production direct'):
        calculate(control,norm2,y1,y2,m)
