"""Regression tests for source-tensor and anomalous spectral matching.

No Wolfram or new RG integration. The final source tensor test uses the
repository's existing Stage22 and SU2 CG JSON, so it also catches stale data.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
import pytest

from Numerical.experimental.T3SelfBilinearAnomalousKernel import (
    MODELS, EXTRA, model_alpha, native_self_metric, parse_matchete_sparse,
    verify_native_self_metric, topology_tensor, tensor_ratio,
    anomalous_phase, anomalous_from_unit_kernels, unit_topology_kernels,
)

_A_SPARSE=('SparseArray[Automatic, {2, 2, 3}, 0, {1, '
 '{{0, 3, 6}, {{1, 1}, {1, 3}, {2, 2}, {1, 2}, {2, 1}, {2, 3}}}, '
 '{3^(-1/4), I/3^(1/4), -3^(-1/4), -3^(-1/4), -3^(-1/4), I/3^(1/4)}}]')
_E_SPARSE=('SparseArray[Automatic, {2, 2, 3, 3}, 0, {1, '
 '{{0, 6, 12}, {{1, 1, 2}, {1, 2, 1}, {1, 2, 3}, {1, 3, 2}, {2, 1, 3}, {2, 3, 1}, '
 '{1, 1, 3}, {1, 3, 1}, {2, 1, 2}, {2, 2, 1}, {2, 2, 3}, {2, 3, 2}}}, '
 '{1/Sqrt[2], -(1/Sqrt[2]), (-I)/Sqrt[2], I/Sqrt[2], (-I)/Sqrt[2], I/Sqrt[2], '
 '(-I)/Sqrt[2], I/Sqrt[2], 1/Sqrt[2], -(1/Sqrt[2]), I/Sqrt[2], (-I)/Sqrt[2]}}]')


@pytest.mark.parametrize('model',MODELS)
def test_model_scope(model):
    alpha=model_alpha(model)
    assert alpha in (0,-2)
    assert EXTRA[alpha].startswith('lambdaT3GeneralQ')


def test_sparse_A_native_source():
    t=parse_matchete_sparse(_A_SPARSE)
    assert t.shape==(2,2,3)
    a=3**(-.25)
    assert t[0,0,0]==pytest.approx(a)
    assert t[0,0,2]==pytest.approx(1j*a)
    assert t[1,1,2]==pytest.approx(1j*a)
    assert np.count_nonzero(t)==6


def test_sparse_E_native_source():
    t=parse_matchete_sparse(_E_SPARSE)
    assert t.shape==(2,2,3,3)
    assert t[0,0,0,1]==pytest.approx(2**(-.5))
    assert t[0,1,2,0]==pytest.approx(1j*2**(-.5))
    assert np.count_nonzero(t)==12


@pytest.mark.parametrize('bad',[
    'SparseArray[bad]',
    _A_SPARSE.replace('{0, 3, 6}','{1, 3, 6}'),
    _A_SPARSE.replace('3^(-1/4)', 'x',1),
    _A_SPARSE.replace('{1, 1}, {1, 3}', '{1, 1}, {1, 1}'),
])
def test_sparse_bad_rejected(bad):
    with pytest.raises(ValueError):parse_matchete_sparse(bad)


@pytest.mark.parametrize('dim', [1,3])
def test_native_metric(dim):
    m=native_self_metric(dim)
    assert np.allclose(m,np.eye(dim))


_A_Q=["3^(-1/4)","0","(-I)/3^(1/4)","0","-3^(-1/4)","0",
      "0","-3^(-1/4)","0","-3^(-1/4)","0","(-I)/3^(1/4)"]
_E_Q=["0","1/Sqrt[2]","0","-(1/Sqrt[2])","0","I/Sqrt[2]",
      "0","(-I)/Sqrt[2]","0","0","0","I/Sqrt[2]",
      "0","0","0","(-I)/Sqrt[2]","0","0",
      "0","0","I/Sqrt[2]","0","0","0",
      "(-I)/Sqrt[2]","0","0","0","1/Sqrt[2]","0",
      "-(1/Sqrt[2])","0","(-I)/Sqrt[2]","0","I/Sqrt[2]","0"]


def fixture_native_catalog(model):
    ds1,ds2=(1,3) if model[0]=='A' else (3,1) if model[0]=='D' else (3,3)
    dim=ds1 if model.endswith('+0') else ds2
    metric=np.eye(dim)
    sym=(np.einsum('ab,cd->abcd',metric,metric)+
         np.einsum('ac,bd->abcd',metric,metric)+
         np.einsum('ad,bc->abcd',metric,metric))/3
    def comp(x):
        if abs(x)<1e-12:return '0'
        if abs(x-1)<1e-12:return '1'
        if abs(x-1/3)<1e-12:return '1/3'
        raise AssertionError(f'Unexpected self metric source component {x}')
    alpha=model_alpha(model)
    self_symbol='lambdaT3ExtraQ000400I1' if alpha==0 else 'lambdaT3ExtraQ000004I1'
    hh_symbol='lambdaT3ExtraQ021010I1' if alpha==0 else 'lambdaT3ExtraQ020101I1'
    species=['Hdag','Hdag','S1','S2'] if alpha==0 else ['Hdag','Hdag','S1dag','S2dag']
    comps=_E_Q if model[0]=='E' else _A_Q
    stage={'format':'stage22-exact-tensors-v1','source':'MatcheteInvariantTensors',
          'model_count':16,'quartic_slot_count':50,'models':{model:{'tensors':[
              {'symbol':self_symbol,'dimensions':[dim]*4,
               'components':[comp(v) for v in sym.reshape(-1)]},
              {'symbol':hh_symbol,'species':species,
               'dimensions':[2,2,ds1,ds2],'potential_prefactor':'1',
               'components':comps.copy()}
          ]}}}
    group={'status':'Success','CGRegistry':[
        {'Name':'T3MixCG','TensorInputForm':_E_SPARSE if model[0]=='E' else _A_SPARSE}]}
    return stage,group,ds1,ds2


@pytest.mark.parametrize('model',MODELS)
def test_six_explicit_native_tensor_contractions(model):
    stage,group,ds1,ds2=fixture_native_catalog(model)
    info=tensor_ratio(stage,group,model,ds1,ds2)
    assert info['ratio']==pytest.approx(1.0)
    assert info['relative_tensor_residual']<1e-12
    assert info['native_self_metric_residual']<1e-12


def test_tensor_change_is_not_silently_accepted():
    stage,group,d1,d2=fixture_native_catalog('E_+0')
    stage['models']['E_+0']['tensors'][1]['components'][1]='2'
    with pytest.raises(ValueError):tensor_ratio(stage,group,'E_+0',d1,d2)


@pytest.mark.parametrize('alpha,expected',[
    (0,1),(-2,1),(0,-1),(-2,-1),
])
def test_anomalous_phase_and_sign(alpha,expected):
    mu=complex(0,2*expected)
    phase=anomalous_phase(mu,alpha)
    assert phase==pytest.approx((-1j if alpha==0 else 1j)*expected)


def test_zero_mu_zero_anomalous():
    arr=np.eye(3,dtype=complex)
    z=anomalous_from_unit_kernels(arr,2*arr,mu_sq=0j,alpha=0,matchete_extra=1+2j)
    assert not np.any(z)


def test_expected_spectral_difference_and_quartic_conjugation():
    a=np.eye(3,dtype=complex)*3
    b=np.eye(3,dtype=complex)
    # mu>0: G_anomalous proportional (G_+-G_-)/2; extra HH is Q^dag.
    z=anomalous_from_unit_kernels(a,b,mu_sq=1+0j,alpha=-2,matchete_extra=2+3j)
    assert np.allclose(z,np.eye(3)*(2-3j))
    w=anomalous_from_unit_kernels(a,b,mu_sq=1j,alpha=0,matchete_extra=2+0j)
    assert np.allclose(w,-2j*np.eye(3))


def test_no_accidental_even_mass_dependence():
    a=np.eye(3,dtype=complex)*.6
    b=np.eye(3,dtype=complex)*.6
    assert not np.any(anomalous_from_unit_kernels(a,b,mu_sq=3,alpha=0,matchete_extra=2))


def test_anomalous_bad_shape_and_asymmetry():
    with pytest.raises(ValueError):
        anomalous_from_unit_kernels(np.eye(2),np.eye(2),mu_sq=1,alpha=0,matchete_extra=1)
    a=np.eye(3,dtype=complex);a[0,1]=1
    with pytest.raises(ValueError):
        anomalous_from_unit_kernels(a,np.eye(3),mu_sq=1,alpha=0,matchete_extra=1)


@dataclass(frozen=True)
class EFT:
    mS1Sq:float
    mS2Sq:float
    lambdaT3:complex
    def validated(self):return self

class Kernel:
    def evaluate(self,_uv,eft):
        return np.eye(3)*eft.lambdaT3/(1+eft.mS1Sq/100)

class Basis:
    species='S1'
    eigen_masses_squared=(80,120)


def test_unit_kernel_recovers_original_normal():
    eft=EFT(100,200,3+4j)
    plus=Kernel().evaluate(None,EFT(80,200,eft.lambdaT3))
    minus=Kernel().evaluate(None,EFT(120,200,eft.lambdaT3))
    source=.5*(plus+minus)
    u,v,err=unit_topology_kernels(Kernel(),None,eft,Basis(),source)
    assert err<1e-14
    assert np.allclose(.5*(u+v)*eft.lambdaT3,source)


def test_unit_kernel_rejects_non_linear_source():
    class Nonlinear(Kernel):
        def evaluate(self,uv,eft):
            return np.eye(3)*eft.lambdaT3**2/(1+eft.mS1Sq/100)
    eft=EFT(100,200,2+0j)
    source=.5*(Nonlinear().evaluate(None,EFT(80,200,2))+
               Nonlinear().evaluate(None,EFT(120,200,2)))
    with pytest.raises(ValueError):unit_topology_kernels(Nonlinear(),None,eft,Basis(),source)


def test_repository_native_six_model_tensor_ratios():
    root=Path(__file__).resolve().parents[1]
    saved=root/'output/general_t3_stage22/matchete_tensors.json'
    if not saved.is_file():
        pytest.skip('run with complete local T3 repository source tensor artifacts')
    stage=json.loads(saved.read_text(encoding='utf-8-sig'))
    dims={'A':(1,3),'D':(3,1),'E':(3,3)}
    for model in MODELS:
        file=root/f'output/group_factors/T3_{model[0]}_uv_cg_registry.json'
        assert file.is_file(),file
        group=json.loads(file.read_text(encoding='utf-8-sig'))
        info=tensor_ratio(stage,group,model,*dims[model[0]])
        assert info['ratio']==pytest.approx(1)
        assert info['relative_tensor_residual']<2e-12
