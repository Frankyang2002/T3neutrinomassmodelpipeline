"""A(-2) native scalar mass-mapping checks; no Wolfram required for unit tests."""
import math
import numpy as np
import pytest

from scripts.VerifyA2BilinearNormalization import (
    integral_triangle, parse_native_metric, verify_mass_map,
)


def example(k='1/Sqrt[3]'):
    return {
        'format':'t3-a-minus2-native-bilinear-tensor-v1',
        'model':'A_-2','status':'NativeMatcheteBilinearMetricExtracted',
        'native_metric_is_factor_times_identity':True,
        'mass_operator_in_full_lagrangian':True,
        'species':['S2dag','S2dag'],
        'native_metric_factor_inputform':k,
        'native_metric_inputform':[[k,'0','0'],['0',k,'0'],['0','0',k]],
    }


def test_native_matchete_tensor_factor_is_read_not_assumed():
    k=parse_native_metric(example())
    assert np.isclose(k,1/math.sqrt(3))
    assert parse_native_metric(example('1'))==1


def test_offdiagonal_tensor_rejected():
    obj=example()
    obj['native_metric_inputform'][0][1]='1'
    with pytest.raises(ValueError,match='CG normalization'):
        parse_native_metric(obj)


@pytest.mark.parametrize('k',[1,1/math.sqrt(3),-1j/math.sqrt(3)])
@pytest.mark.parametrize('phase',[0,0.45,2.31])
def test_lagrangian_vs_potential_real_triplet_mass_map(k,phase):
    m2=1e6
    mu=5000*complex(math.cos(phase),math.sin(phase))
    result=verify_mass_map(k,m2,mu)
    assert result['potential_max_relative_error']<1e-12
    assert np.allclose(sorted(result['masses_squared_gev2']),[m2-10000,m2+10000])
    assert np.isclose(result['mass_splitting_squared_gev2'],20000)


def test_tachyonic_bilinear_rejected():
    with pytest.raises(ValueError,match='tachyonic'):
        verify_mass_map(1,100,60j)


def test_euclidean_triangle_symmetric_and_positive():
    a=integral_triangle(100000,1200,1050)
    b=integral_triangle(100000,1050,1200)
    assert float(a)>0
    assert abs(a/b-1)<1e-12
    assert integral_triangle(100000,1200,1100)<integral_triangle(100000,1200,1000)


def test_one_scale_kernel_against_saved_matchete_if_available():
    from pathlib import Path
    from scripts.VerifyA2BilinearNormalization import ROOT, check_one_scale_loop_shape
    path=ROOT/'output/hypercharge/T3_A_alpha_m2/data/c5_coefficient.txt'
    if not path.is_file():
        pytest.skip('Saved one-scale A(-2) Matchete source not included')
    result=check_one_scale_loop_shape(path,1e5,1200,1.02e6,0.98e6)
    assert result['maximum_relative_ratio_error']<2e-9
