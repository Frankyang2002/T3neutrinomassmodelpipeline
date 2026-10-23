import numpy as np
import pytest
from Numerical.diagnostics.T3BHeavyLogAlignment import analyze, flavor_kernel_by_heavy, fit_shape


def test_heavy_resolved_reconstructs_exact_target():
    y1 = np.array([[1.,2.],[2.,1.],[1.,3.]])
    y2 = np.array([[3.,2.],[1.,4.],[2.,1.]])
    masses = np.array([100.,200.])
    kernel = flavor_kernel_by_heavy(y1,y2,masses)
    target = 2.5 * np.einsum('r,rij->ij',np.log(10/masses),kernel)
    report = analyze(target,y1,y2,masses,10,150)
    assert report['heavy_resolved_log']['alignment_residual'] < 1e-12
    assert report['heavy_logs_explain_mismatch']
    assert report['common_log']['alignment_residual'] > 1e-4


def test_symmetric_complex_boundary():
    y1 = np.array([[1+1j],[2-1j],[0+2j]])
    y2 = np.array([[1-2j],[0+1j],[2+0j]])
    k = flavor_kernel_by_heavy(y1,y2,[100])[0]
    np.testing.assert_allclose(k,k.T,atol=0)


def test_nonproportional_matrix_flagged():
    a=np.eye(3,dtype=complex)
    b=np.diag([1.,2.,3.])
    _,r=fit_shape(a,b)
    assert r>0.1


def test_rejects_bad_masses():
    with pytest.raises(ValueError):
        flavor_kernel_by_heavy(np.ones((3,1)),np.ones((3,1)),[-1])
