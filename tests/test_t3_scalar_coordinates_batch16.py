"""Pure numerical tests for Batch16 coordinate-span logic.

These do not independently validate Matchete/RGBeta physics or run the heavy
repository benchmark; run the CLI on saved Batch13 data for that.
"""
from __future__ import annotations
import numpy as np
import pytest

from Numerical.parallel_t3.complete.two_model_scalar_coordinate_batch16 import (
    fit_sector,rebuilt_quartic_delta,relative,_group,transform_tensor)


def _e(i):
    a=np.zeros((3,3),float)
    a.reshape(-1)[i]=1.
    return a


def test_nontrivial_two_by_two_basis_and_sign():
    # One native invariant is -R1+2R2; the second is 3R1+R2.
    x,y=_e(0),_e(1)
    out=fit_sector({'r1':x,'r2':y},
                   {'n1':-x+2*y,'n2':3*x+y})
    assert out['unique_bijective_coordinate_map']
    np.testing.assert_allclose(out['native_to_reference_coefficients'],
                               [[-1,3],[2,1]],atol=1e-12)
    inv=np.asarray(out['reference_to_native_coefficients'])
    np.testing.assert_allclose(np.array([[-1,3],[2,1]])@inv,np.eye(2),atol=1e-12)
    assert out['reference_direction_fit_relative']<1e-12


def test_partial_native_channel_cannot_reconstruct_two_reference_directions():
    x,y=_e(0),_e(1)
    out=fit_sector({'r1':x,'r2':y},{'n':x+y})
    assert out['native_in_reference_span']
    assert not out['reference_in_native_span']
    assert out['archived_native_columns_missing']==1
    assert not out['unique_bijective_coordinate_map']
    correction,audit=rebuilt_quartic_delta({
        'portal':{'fit':out,'reference':{'r1':x,'r2':y},'native':{'n':x+y}}},
        {'r1':1.,'r2':2.})
    assert not np.any(correction)
    assert audit['portal']['mode']=='RetainRGBetaReference_BasisIncomplete'


def test_unique_native_reconstruction_preserves_quartic_tensor():
    x=np.zeros((12,)*4);y=np.zeros_like(x)
    x[0,0,0,0]=2.;y[1,1,1,1]=4.
    native={'n1':-x+2*y,'n2':3*x+y}
    reference={'r1':x,'r2':y}
    fit=fit_sector(reference,native)
    corr,audit=rebuilt_quartic_delta({
        'self':{'fit':fit,'reference':reference,'native':native}},
        {'r1':0.3,'r2':-0.9})
    assert np.linalg.norm(corr)<1e-12
    assert audit['self']['mode']=='SourceNativeReconstructionOfRGBetaReference'


def test_empty_native_sector_is_explicitly_incomplete():
    out=fit_sector({'r':_e(0)}, {})
    assert out['native_names']==[]
    assert out['reference_to_native_coefficients'] is None
    assert not out['reference_in_native_span']


def test_symmetry_rotation_identity():
    x=np.zeros((12,)*4)
    x[2,2,2,2]=3.
    np.testing.assert_allclose(transform_tensor(x,np.eye(12)),x)


def test_group_inference_is_conservative():
    assert _group('λ')=='H4'
    assert _group('lambdaH2Inv2')=='H2S2_2'
    assert _group('lambdaS1Inv1')=='S1_4'
    assert _group('lambdaT3') is None
    assert _group('lambdaT3GeneralQ000004I1') is None


def test_bad_complex_quartic_rejected():
    with pytest.raises(ValueError,match='expected a real'):
        fit_sector({'r':_e(0)+1j*_e(1)}, {'n':_e(0)})


def test_relative_zero_and_conventional_positive_norm():
    assert relative(np.zeros(2),np.zeros(2))==0
    assert relative(np.array([2,0]),np.array([2,0]))==0
