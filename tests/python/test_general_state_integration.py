"""Integration tests for extended canonical T3 state/ODE interfaces."""
import sys
import types
import numpy as np
import pytest

# In the complete project these functions come from common.T3Model.
# This fallback is only for standalone execution of the delivery ZIP.
try:
    from common.T3Model import valid_t3_dimensions
except ImportError:
    common=types.ModuleType('common')
    model=types.ModuleType('common.T3Model')
    model.valid_t3_dimensions=lambda a,b,c:(a,b,c) in ((1,3,2),(2,2,1),(2,2,3),(3,1,2),(3,3,2))
    model.valid_shared_scalar_dimensions=lambda a,c:(a,c) in ((2,1),(2,3))
    sys.modules['common']=common
    sys.modules['common.T3Model']=model

from Numerical.core.State import SMNumericalState,T3Representation,T3UVState
from Numerical.core.StateVector import pack_uv_state,unpack_uv_state
try:
    from Numerical.core.BetaVector import validate_beta_keys
except ImportError:
    beta=types.ModuleType('Numerical.core.BetaVector')
    def validate_beta_keys(values,layout,*,allow_extra=False):
        assert set(values)==set(layout.names)
    beta.validate_beta_keys=validate_beta_keys
    beta.beta_values_to_derivative=lambda *args,**kwargs: None
    sys.modules['Numerical.core.BetaVector']=beta
from Numerical.core.RGBetaEvaluator import state_environment,evaluate_inputform_expression


def _state(extra=None):
    z=np.zeros((3,3),dtype=complex)
    sm=SMNumericalState(.35,.65,1.2,.13,z,z,z)
    return T3UVState(1e7,T3Representation(1,3,2,-2),sm,z,z,np.eye(3),
                     1e6,4e6,.1,.2,.3,.4,.5,.01+.02j,
                     lambdaH2Adj=.02,lambdaS2Adj=.03,
                     general_couplings=extra or {})


def test_original_inventory_unchanged():
    base=_state().validated()
    vector,layout=pack_uv_state(base)
    assert not any(n.startswith('lambdaT3GeneralQ') for n in layout.names)
    assert len(layout.names)==len(set(layout.names))
    again=unpack_uv_state(vector,layout)
    assert again.general_couplings=={}


def test_extra_coupling_roundtrip_and_expression():
    extra={'lambdaT3GeneralQ000004I1':.02+.03j,'muT3GeneralS2Sq':1.5+2.5j,
           'y2Conj':np.eye(3,dtype=complex)*(.01+.04j)}
    initial=_state(extra)
    vector,layout=pack_uv_state(initial)
    final=unpack_uv_state(vector,layout)
    for name,value in extra.items():
        np.testing.assert_allclose(final.general_couplings[name],value)
    environment=state_environment(final)
    assert evaluate_inputform_expression('Bar[lambdaT3GeneralQ000004I1]*muT3GeneralS2Sq',environment)==pytest.approx((.02-.03j)*(1.5+2.5j))
    assert layout.block('y2Conj').shape==(3,3)


def test_invalid_unknown_and_collision_fail():
    with pytest.raises(ValueError): _state({'unsupportedSymbol':0}).validated()
    with pytest.raises(ValueError): _state({'lambdaT3':.1}).validated()
    with pytest.raises(ValueError): _state({'y2Conj':np.zeros((2,2))}).validated()
