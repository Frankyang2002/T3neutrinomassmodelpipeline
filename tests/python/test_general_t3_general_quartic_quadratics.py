"""Canonical spherical-basis quadratic mass-mixing prototype tests."""
import json
from pathlib import Path
import numpy as np
import pytest

from Numerical.experimental.GeneralT3ScalarQuadraticMasses import quadratic_mass_hessian,quadratic_structure

ROOT=Path(__file__).resolve().parents[2]
MAN=json.loads((ROOT/'docs/general_scalar_trial_operator_manifest.json').read_text())


def test_eight_extra_quadratics_total():
    assert sum(len(quadratic_structure(MAN,k)) for k in MAN['models'])==8


def test_zero_coupling_recovers_diagonal_masses():
    for name,body in MAN['models'].items():
        m=quadratic_mass_hessian(MAN,name,1350**2,1780**2)
        d1,d2=m['dimensions']
        got=np.sort(m['mass_squared_eigenvalues_gev2'])
        want=np.sort([1350**2]*(2*d1)+[1780**2]*(2*d2))
        assert np.array_equal(got,want)
        assert not m['matchete_mass_basis_verified']
        assert not m['included_in_threshold_matching']


def test_doublet_mixing_gives_known_two_by_two_eigenvalues():
    m1,m2=1350**2,1780**2
    mu=complex(40000,30000)
    expected=[(m1+m2)/2+s*np.sqrt(((m1-m2)/2)**2+abs(mu)**2) for s in (-1,1)]
    for model in ('B_-1','C_-1'):
        eig=np.asarray(quadratic_mass_hessian(MAN,model,m1,m2,mu)['mass_squared_eigenvalues_gev2'])
        assert np.allclose(np.sort(eig),np.sort([expected[0]]*4+[expected[1]]*4),rtol=1e-12)


def test_singlet_holomorphic_mass_splits_two_real_degrees():
    m=quadratic_mass_hessian(MAN,'D_-2',2.0e6,1.4e6,10000+20000j)
    vals=m['mass_squared_eigenvalues_gev2']
    expected=(1.4e6-2*abs(10000+20000j),1.4e6+2*abs(10000+20000j))
    assert any(np.isclose(v,expected[0]) for v in vals)
    assert any(np.isclose(v,expected[1]) for v in vals)


def test_real_hessian_symmetric_and_correct_dimensions():
    for model in ('A_-2','A_+0','B_-1','C_-1','D_-2','D_+0','E_-2','E_+0'):
        result=quadratic_mass_hessian(MAN,model,1.5e6,2.2e6,complex(5e4,2e4))
        h=np.asarray(result['real_scalar_hessian_gev2'])
        d1,d2=result['dimensions']
        assert h.shape==(2*(d1+d2),2*(d1+d2))
        assert np.allclose(h,h.T,rtol=0,atol=1e-12)
        assert np.all(np.isfinite(h))


def test_no_extra_operator_rejects_nonzero_mass_mix():
    with pytest.raises(ValueError,match='no allowed'):
        quadratic_mass_hessian(MAN,'A_-4',1e6,2e6,1+1j)


def test_non_finite_input_rejected():
    with pytest.raises(ValueError,match='Nonfinite'):
        quadratic_mass_hessian(MAN,'B_-1',np.inf,1e6)
    with pytest.raises(ValueError,match='Nonfinite'):
        quadratic_mass_hessian(MAN,'B_-1',1e6,1e6,complex(np.nan,0))


def test_metric_invariant_symmetry_for_j1():
    from Numerical.experimental.GeneralT3ScalarQuadraticMasses import _metric
    g=_metric(3)
    assert np.array_equal(g,g.T)
    assert np.array_equal(g@g,np.eye(3))
    assert np.array_equal(_metric(2).T,-_metric(2))
