import numpy as np
import pytest
from types import SimpleNamespace
from Numerical.running.T3BWilsonBackground import background_from_trajectory


def _state(yukawa_shape=(1,1)):
    rep=SimpleNamespace(shared_scalar=False,d_s1=2,d_s2=2,d_f=1,alpha=-1)
    sm=SimpleNamespace(g2=0.6,lambdaH=0.25, yu=np.ones(yukawa_shape),
                       yd=np.ones(yukawa_shape),ye=np.ones(yukawa_shape))
    return SimpleNamespace(representation=rep,sm=sm,lambda12=0.1,
                           lambda12Adj=0.02,lambdaT3=0.03+0.04j)


def test_exports_scalar_trajectory():
    t=SimpleNamespace(solver_success=True,n_points=2,mu_gev=[1e5,1e3],
                      state_at_index=lambda i:_state())
    result=background_from_trajectory(t)
    assert len(result['samples']) == 2
    assert result['samples'][0]['lambdaT3'] == [0.03, 0.04]


def test_refuses_fake_three_generation_reduction():
    t=SimpleNamespace(solver_success=True,n_points=2,mu_gev=[1e5,1e3],
                      state_at_index=lambda i:_state((3,3)))
    with pytest.raises(ValueError, match='one-generation'):
        background_from_trajectory(t)
