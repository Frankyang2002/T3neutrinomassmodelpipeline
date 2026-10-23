import numpy as np
from Numerical.diagnostics.T3BFlavorMixingBridge import mixed_boundary,offdiagonal_evolution


def test_boundary_symmetric_and_one_generation_component():
    a=np.diag([1+1j,2,3]).astype(complex)
    b=np.diag([2,3-1j,4]).astype(complex)
    masses=np.array([10,20,30])
    got=mixed_boundary(a,b,masses)
    assert np.allclose(got,got.T)
    assert np.allclose(np.diag(got),np.conjugate(np.diag(a))*np.conjugate(np.diag(b))/masses)


def test_forward_one_loop_and_reverse_feedback():
    rows=[{'mu_gev':1e5,'lambdaT3':[1,0]},{'mu_gev':1e3,'lambdaT3':[1,0]}]
    c12=np.diag([1e-5,2e-5,3e-5]).astype(complex)
    data=offdiagonal_evolution(rows,c12)
    first=data['c5_at_ms_forward_only_running_lambda']
    ll=data['c5_at_ms_frozen_mf_leading_log']
    a=np.asarray(first['real'])+1j*np.asarray(first['imag'])
    b=np.asarray(ll['real'])+1j*np.asarray(ll['imag'])
    assert np.allclose(a,b,atol=1e-16,rtol=1e-8)
    assert 0<data['relative_feedback_on_offdiagonal_c5']<0.05
    assert data['c12_transpose_residual']<1e-15
    assert data['c5_transpose_residual']<1e-15
