import numpy as np
from Numerical.diagnostics.T3BExploratoryCoupledBenchmark import _gamma, _rhs


def sample():
    return {'g2': 0.6, 'lambdaH': 0.2, 'lambda12': 0.01,
            'lambda12Adj': 0.02, 'lambdaT3': 0j,
            'ye': np.diag([.01, .02, .03]).astype(complex),
            'yu': np.diag([.001, .01, .7]).astype(complex),
            'yd': np.diag([.001, .003, .02]).astype(complex)}


def test_weinberg_matches_repository_formula():
    from RGE.running.weinberg.WeinbergRGE import sm_weinberg_beta
    c = sample()
    b = np.array([[2, 1j, 0], [1j, 3, .3], [0, .3, 4]], complex) * 1e-13
    g12, g5, a = _gamma(c)
    expected = g5*b-1.5*(a@b+b@a.T)
    actual = sm_weinberg_beta(gY=.3, g2_value=c['g2'], g3=.9,
        lambdaH_value=c['lambdaH'], ye=c['ye'], yu=c['yu'],
        yd=c['yd'], K=b)[-1]
    assert np.allclose(expected, actual, rtol=1e-12, atol=1e-28)
    y=np.concatenate([np.zeros(9,complex),b.reshape(-1)])
    derived=_rhs(0.,y,lambda t:c,diagonal=True,reverse=False)[9:].reshape(3,3)
    assert np.allclose(derived*16*np.pi**2, actual, rtol=1e-12, atol=1e-28)


def test_one_generation_c12_coefficient_is_preserved():
    c=sample();c['ye']=np.zeros((3,3),complex);c['ye'][0,0]=.02
    a=np.zeros((3,3),complex);a[0,0]=1.
    y=np.concatenate([a.ravel(),np.zeros(9,complex)])
    beta=_rhs(0.,y,lambda t:c,diagonal=True,reverse=False)[:9].reshape(3,3)*16*np.pi**2
    expected=(-3*.6**2+2*.01-.02/2+.02**2)
    assert abs(beta[0,0]-expected)<1e-12
