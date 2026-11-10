"""Fast algebra and convention tests; integration tests use local source exports."""
from __future__ import annotations

import math
from pathlib import Path
import numpy as np
import pytest

from Numerical.parallel_t3.complete.two_model_conditional_uvfit_batch12 import (
    NativeSoftFast, flatten_uv, unflatten_uv, candidate_c5_as_smeft,
    PARAMS, B11, FORMAT,
)
from Numerical.parallel_t3.complete.two_model_native_laurent_batch11 import (
    SOFT_PATHS, parse_laurent, evaluate_soft,
)

W='Coupling[T3EFT1W12,{Index[d$$1,Flavor],Index[d$$2,Flavor]},0]'
EX='Coupling[T3EFT1WExtra,{Index[d$$1,Flavor],Index[d$$2,Flavor]},0]'
L='Coupling[lambdaT3,{},0]'
LX='Coupling[lambdaT3ExtraQ020101I1,{},0]'
M1='Coupling[MS1,{},0]'
M2='Coupling[MS2,{},0]'
MM='Coupling[T32MinusMass,{},0]'
MP='Coupling[T32PlusMass,{},0]'


def test_real_18_uv_parameterization_exact_roundtrip():
    a=np.arange(1.,10.).reshape(3,3)*1e-4
    b=(-np.arange(10.,19.).reshape(3,3))*1e-5
    vector=flatten_uv(a,b)
    assert vector.shape==(18,)
    x,y=unflatten_uv(vector)
    np.testing.assert_array_equal(x,a)
    np.testing.assert_array_equal(y,b)
    assert vector[9]==b[0,0]
    with pytest.raises(ValueError,match='18'):
        unflatten_uv(np.zeros(17))
    with pytest.raises(ValueError,match='REAL'):
        flatten_uv(a.astype(complex),b)


def test_gamma_factor_is_explicit_and_majorana_transpose_not_dagger():
    a=np.array([[1+2j,2+3j,3],[2+3j,4+5j,6],[3,6,7j]],complex)*1e-15
    np.testing.assert_array_equal(candidate_c5_as_smeft(a,-2),-2*a)
    with pytest.raises(ValueError,match='signed|Signed'):
        candidate_c5_as_smeft(a,0)
    with pytest.raises(ValueError,match='symmetric'):
        candidate_c5_as_smeft(a+np.eye(3,1)*1e-13,1)


def _compare_fast_with_sympy(text,model,*,w1,x,lam,ex,m1=1200.,m2=1000.,
                             mminus=980.,mplus=1020.):
    fast=NativeSoftFast(text,model)
    reg,pole=parse_laurent(text,model)
    kw=dict(w12=w1,wextra=x,lam=lam,extra=ex,m1=m1,m2=m2,
            mminus=mminus,mplus=mplus,mu_bar=1000.)
    slow= evaluate_soft(reg,pole,**kw)
    new=fast.evaluate(**kw)
    for got,want in zip(new,slow):
        np.testing.assert_allclose(got,want,rtol=2e-12,atol=1e-24)


def test_exact_A_source_vectorized_retains_complex_phases():
    delta=f'({M1}^2-{M2}^2)'
    text=(f'Sqrt[3]*hbar*Bar[{W}]*{L}*'
          f'(1/ϵ+Log[μbar2/{M1}^2]+{M1}^2/{delta})'
          f'+Sqrt[3]*hbar*Bar[{W}]*Bar[{L}]')
    a=np.array([[1+2j,2j,3],[4,5+1j,6],[7,8,9j]],complex)*1e-9
    _compare_fast_with_sympy(text,'A_-4',w1=a,x=np.zeros((3,3)),
                             lam=.11+.04j,ex=0)


def test_D_split_real_and_extra_Wilson_vectorized():
    text=(f'Sqrt[3]*hbar*Bar[{LX}]*Bar[{EX}]*'
          f'({M1}^2+{MM}^2+{MP}^2)/(ϵ*{M1}^2)'
          f'+hbar*{L}*Bar[{W}]*Log[μbar2/{MP}^2]')
    a=np.array([[1+1j,2,3],[4j,5,6],[7,8,9]],complex)*1e-9
    b=np.array([[1j,3,1],[2,1j,3],[4,2,1]],complex)*1e-9
    _compare_fast_with_sympy(text,'D_-2',w1=a,x=b,lam=.13-.03j,ex=.027+.012j)


def test_safeguards_unrecognized_native_coupling():
    text=(f'hbar*Bar[{W}]*{L}*(1/ϵ+Log[μbar2/{M2}^2])'
          '+Coupling[IncorrectExtra,{},0]')
    with pytest.raises(ValueError,match='unrecognised'):
        NativeSoftFast(text,'A_-4')


@pytest.mark.parametrize('model',['A_-4','D_-2'])
def test_saved_native_multiplet_sources_fast_match_symbolic_if_present(model):
    root=Path(__file__).resolve().parents[1]
    path=root/SOFT_PATHS[model]
    if not path.is_file():
        pytest.skip('Actual native saved source available only in full checkout')
    text=path.read_text(encoding='utf-8-sig')
    rng=np.random.default_rng(77)
    a=(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))*1e-10
    b=(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))*1e-10
    _compare_fast_with_sympy(text,model,w1=a,x=b,lam=.1,ex=.017)


def test_packaged_source_paths_are_platform_agnostic():
    assert B11.as_posix().endswith('thesis_native_laurent_batch11_v1')
    assert FORMAT.endswith('batch12-v1')
    assert len(PARAMS)==14
