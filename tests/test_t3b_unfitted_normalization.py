import numpy as np
import pytest
from Numerical.diagnostics.T3BUnfittedDirectNormalization import evaluate


def inputs():
    y1 = np.array([[1., 0.2], [0.3, 1.], [0.1, 0.4]], complex)
    y2 = np.array([[0.1, 0.4], [1., 0.2], [0.3, 0.5]], complex)
    m = np.array([1e5, 1.2e5])
    ms = 1e3
    return y1, y2, m, ms


def test_unfitted_normalization():
    from Numerical.diagnostics.T3BHeavyLogAlignment import flavor_kernel_by_heavy
    y1, y2, m, ms = inputs()
    lam = 1.1
    source = np.einsum('r,rij->ij', np.log(ms/m), flavor_kernel_by_heavy(y1,y2,m))
    target = (2 * 4*np.sqrt(3)/3 * lam/(16*np.pi**2))*source
    res = evaluate(target, y1,y2,m,ms,lam)
    assert res['status'] == 'UnfittedNormalizationPassed'
    assert res['relative_frobenius_residual'] < 1e-14


def test_mismatch_is_rejected():
    y1,y2,m,ms = inputs()
    target = np.eye(3, dtype=complex)
    assert evaluate(target,y1,y2,m,ms,1.1)['status'] == 'UnfittedNormalizationMismatch'


def test_invalid_mass_order():
    y1,y2,m,ms = inputs()
    with pytest.raises(ValueError):
        evaluate(np.eye(3),y1,y2,m,2e5,1)
