"""Run these tests in the full T3 project root."""
import numpy as np
import pytest

from RGE.running.intermediate.T3BC12FlavorRGE import (
    REVERSE_COEFFICIENT, c12_beta_numerator,
    verify_lepton_flavour_action_against_master,
)


def test_one_generation_limit():
    a=np.zeros((3,3),complex);b=a.copy();ye=a.copy()
    a[0,0]=2.+1j;ye[0,0]=.2+.1j
    result=c12_beta_numerator(a,b,ye=ye,g2=.6,lambda12=.04,
                              lambda12Adj=.02,lambdaT3=1.)
    expected=(-3*.6**2+2*.04-.02/2+abs(ye[0,0])**2)*a[0,0]
    np.testing.assert_allclose(result[0,0],expected,rtol=1e-14,atol=1e-14)


def test_non_diagonal_complex_yukawa_and_symmetry():
    rng=np.random.default_rng(42)
    ye=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    a=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));a=(a+a.T)/2
    b=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));b=(b+b.T)/2
    e=ye@ye.conj().T
    expected=(-3*.65**2+2*.1-.02/2)*a +.5*(e@a+a@e.T)
    expected+=REVERSE_COEFFICIENT*np.conj(.1+.02j)*b
    actual=c12_beta_numerator(a,b,ye=ye,g2=.65,lambda12=.1,
                              lambda12Adj=.02,lambdaT3=.1+.02j)
    np.testing.assert_allclose(actual,expected,atol=5e-14,rtol=5e-14)
    np.testing.assert_allclose(actual,actual.T,atol=5e-14)


def test_lepton_basis_covariance():
    rng=np.random.default_rng(123)
    q,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
    ye=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    a=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));a=(a+a.T)/2
    b=np.zeros((3,3),complex)
    opts=dict(g2=.7,lambda12=.03,lambda12Adj=.01,lambdaT3=.04)
    original=c12_beta_numerator(a,b,ye=ye,**opts)
    transformed=c12_beta_numerator(q@a@q.T,b,ye=q@ye,**opts)
    np.testing.assert_allclose(transformed,q@original@q.T,rtol=1e-12,atol=1e-12)


def test_master_tensor_two_generations():
    result=verify_lepton_flavour_action_against_master()
    assert result['status']=='C12LeptonFlavorSourceVerified', result
    assert all(r=='0' for r in result['master_minus_matrix_residuals'].values())
    assert not result['unexpected_terms']


def test_rhs_uses_c12_source_equation():
    from Numerical.diagnostics.T3BExploratoryCoupledBenchmark import _rhs
    rng=np.random.default_rng(5)
    a=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));a=(a+a.T)/2
    c={'g2':.66,'lambdaH':.15,'lambda12':.06,'lambda12Adj':.02,
       'lambdaT3':0.,'ye':np.diag([.01,.02,.03]).astype(complex),
       'yu':np.diag([.002,.01,.7]).astype(complex),
       'yd':np.diag([.001,.003,.03]).astype(complex)}
    state=np.r_[a.ravel(),np.zeros(9,complex)]
    beta=_rhs(0.,state,lambda _:c,diagonal=True,reverse=False)[:9].reshape(3,3)
    expected=c12_beta_numerator(a,np.zeros((3,3),complex),ye=c['ye'],g2=c['g2'],
        lambda12=c['lambda12'],lambda12Adj=c['lambda12Adj'],lambdaT3=0.)
    np.testing.assert_allclose(beta*(16*np.pi**2),expected,rtol=1e-12,atol=1e-14)
