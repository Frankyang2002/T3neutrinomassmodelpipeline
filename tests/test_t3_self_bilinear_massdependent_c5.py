"""Mass-dependent six-model anomalous matching regression tests.

These tests validate arithmetic and source contracts only. The complete
numerical run uses the user's existing repository and RGBeta trajectories.
"""
from dataclasses import dataclass
import numpy as np
import pytest

from Numerical.experimental.T3SelfBilinearAnomalousKernel import MODELS
from Numerical.experimental.T3SelfBilinearMassDependentKernel import (
    source_path,verify_legacy_source,split_unit_kernels,spectral_c5,
)


@dataclass(frozen=True)
class EFT:
    mS1Sq:float
    mS2Sq:float
    lambdaT3:complex
    def validated(self): return self

class Basis:
    species='S1'
    eigen_masses_squared=(90.0,110.0)

class Kernel:
    def evaluate(self, _uv, eft):
        # A symmetric mass-dependent one-loop surrogate (for arithmetic only).
        x=eft.lambdaT3*(1/(eft.mS1Sq+50)+0.5/(eft.mS2Sq+20))
        return x*np.array([[1,0.1j,0],[0.1j,2,0.2],[0,0.2,3]],dtype=complex)

class BadKernel:
    def evaluate(self,_uv,eft):
        return eft.lambdaT3*np.eye(3)/(eft.mS2Sq+20)

@pytest.mark.parametrize('model',MODELS)
def test_six_legacy_paths(model):
    p=source_path('/project',model).as_posix()
    alpha=model.split('_')[1]
    tag=('m' if alpha.startswith('-') else 'p')+alpha[1:]
    assert p.endswith(f'T3_{model[0]}_alpha_{tag}/data/c5_coefficient.txt')


def test_massdependent_source_contract(tmp_path):
    p=tmp_path/'c5_coefficient.txt'
    tokens=['Coupling[MF,','Coupling[MS1,','Coupling[MS2,',
            'Coupling[lambdaT3,','Coupling[y1,','Coupling[y2,']
    p.write_text('('+' + '.join(tokens)+')'*60,encoding='utf-8')
    assert verify_legacy_source(p)['status']=='ExplicitMassDependentLegacyOneScaleSource'
    p.write_text(p.read_text().replace('Coupling[MS1,','NO_MS1') ,encoding='utf-8')
    with pytest.raises(ValueError,match='both scalar masses'):
        verify_legacy_source(p)


def test_hierarchical_mass_insensitive_source_rejected(tmp_path):
    p=tmp_path/'c5_coefficient.txt'
    p.write_text('lambdaT3 / MF'*80,encoding='utf-8')
    with pytest.raises(ValueError,match='both scalar masses'):
        verify_legacy_source(p)


def test_unit_kernel_is_mass_sensitive():
    p,m,response=split_unit_kernels(Kernel(),None,EFT(100,200,3+4j),Basis())
    assert response>1e-4
    assert np.linalg.norm(p-m)>0
    assert p.shape==(3,3)


def test_mass_insensitive_kernel_fails_closed():
    with pytest.raises(ValueError,match='insensitive'):
        split_unit_kernels(BadKernel(),None,EFT(100,200,3+4j),Basis())


def test_normal_and_anomalous_spectral_channels():
    plus,minus,_=split_unit_kernels(Kernel(),None,EFT(100,200,3+4j),Basis())
    d=spectral_c5(plus,minus,lam_t3=3+4j,mu_sq=5+0j,
                  alpha=0,extra_matchete=-0.02+0j,tensor_factor=1+0j)
    assert np.allclose(d['normal'],0.5*(3+4j)*(plus+minus))
    assert np.allclose(d['anomalous'],-0.01*(plus-minus))
    assert d['anomalous_over_normal']>0
    assert np.allclose(d['combined'],d['normal']+d['anomalous'])


def test_conjugation_and_zero_coupling_limits():
    plus,minus,_=split_unit_kernels(Kernel(),None,EFT(100,200,3+4j),Basis())
    d=spectral_c5(plus,minus,lam_t3=2,mu_sq=1j,
                  alpha=0,extra_matchete=1+2j,tensor_factor=1)
    assert np.allclose(d['anomalous'],(1-2j)*(-1j)/2*(plus-minus))
    q=spectral_c5(plus,minus,lam_t3=2,mu_sq=0,
                  alpha=0,extra_matchete=1+2j,tensor_factor=1)
    assert not np.any(q['anomalous'])
    q=spectral_c5(plus,minus,lam_t3=2,mu_sq=1j,
                  alpha=0,extra_matchete=0,tensor_factor=1)
    assert not np.any(q['anomalous'])


def test_nonlinear_topology_source_rejected():
    class Nonlinear(Kernel):
        def evaluate(self,uv,eft):
            return super().evaluate(uv,eft)*eft.lambdaT3
    with pytest.raises(ValueError,match='not linear'):
        split_unit_kernels(Nonlinear(),None,EFT(100,200,3+4j),Basis())


def test_repository_six_massdependent_sources_if_available():
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    candidate=source_path(root,'A_+0')
    if not candidate.is_file():
        pytest.skip('full T3 repository source artifacts absent in isolated harness')
    for model in MODELS:
        assert verify_legacy_source(source_path(root,model))['source_length']>120
