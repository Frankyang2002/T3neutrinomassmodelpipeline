"""Pure numerical tests: source sign, H.c., propagator phase, independent loops."""
import numpy as np
import pytest

from scripts.VerifyA2ExtraQuarticC5 import (
    matching_channels, direct_spectral_propagators, scalar_triangle,
    lift_flavour, normalized_error,
)


def test_complex_conjugation_and_rgbeta_sign():
    p=np.diag([2,3,4]).astype(complex)
    m=np.diag([1,2,3]).astype(complex)
    mu=3+4j
    normal,anomalous=matching_channels(p,m,lambda_t3=0.1+0.2j,
                          lambda_rgbeta=0.4+0.3j,mu_rgbeta=mu,tensor_ratio=1)
    assert np.allclose(normal,0.5*(0.1+0.2j)*(p+m))
    assert np.allclose(anomalous,-0.5*(0.4-0.3j)*(3+4j)/5*(p-m))


def test_zero_and_phase_dependence():
    p=np.eye(3);m=0.8*np.eye(3)
    _,a=matching_channels(p,m,lambda_t3=1,lambda_rgbeta=0.3j,mu_rgbeta=0,tensor_ratio=1)
    assert not np.any(a)
    _,b=matching_channels(p,m,lambda_t3=1,lambda_rgbeta=0,mu_rgbeta=2+3j,tensor_ratio=1)
    assert not np.any(b)
    _,g=direct_spectral_propagators(5,90,70,3+4j)
    expected=0.5*(3+4j)/5*(1/95-1/75)
    assert np.allclose(g,expected,rtol=1e-14)


def test_full_three_flavour_symmetry():
    rng=np.random.default_rng(2512)
    y1=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    y2=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    c=lift_flavour([1e-8,1.1e-8,1.3e-8],y1,y2)
    assert c.shape==(3,3)
    assert np.allclose(c,c.T)
    assert normalized_error(c,c*(1+1e-12))<2e-12


@pytest.mark.parametrize('mf,ms1sq,ms2sq',[(1e5,1.3e6,0.98e6),(1.1e5,1.3e6,1.02e6)])
def test_positive_convergent_scalar_triangle(mf,ms1sq,ms2sq):
    a=float(scalar_triangle(mf,ms1sq,ms2sq))
    b=float(scalar_triangle(mf,ms1sq,ms2sq*1.001))
    assert np.isfinite(a) and a>0
    assert np.isfinite(b) and 0<b<a


def test_fail_closed_invalid_matrices():
    with pytest.raises(ValueError):
        lift_flavour([1,2],np.eye(3),np.eye(3))
    with pytest.raises(ValueError):
        matching_channels(np.eye(2),np.eye(2),lambda_t3=1,
                          lambda_rgbeta=1,mu_rgbeta=1,tensor_ratio=1)
