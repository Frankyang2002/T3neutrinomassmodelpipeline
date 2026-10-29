"""Integration contracts for the six-model normal-channel one-loop continuation."""
from dataclasses import dataclass
from types import SimpleNamespace
import numpy as np
import pytest

from scripts import RunT3SelfBilinearNeutrino as stage
from Numerical.experimental.T3SelfBilinearMassBasis import SelfBilinearMassBasis


@dataclass(frozen=True)
class FakeEFT:
    mS1Sq: float
    mS2Sq: float

    def validated(self):
        if min(self.mS1Sq,self.mS2Sq)<=0:
            raise ValueError('negative scalar mass')
        return self


class SingleScalarLineKernel:
    def evaluate(self,_uv,eft):
        # A deliberately non-affine one-propagator mass dependence.
        scale=1.0/(eft.mS1Sq+430)**2 + 0.07/(eft.mS2Sq+200)
        M=np.array([[1,0.2j,0.3], [0.2j,2,0.4],[0.3,0.4,3]],dtype=complex)
        return scale*M


def test_normal_channel_is_half_sum_not_mass_average():
    basis=SelfBilinearMassBasis.build('S1',3,1200,100+55j)
    eft=FakeEFT(1200,1900)
    base,mixed,parts=stage.normal_propagator_c5(SingleScalarLineKernel(),None,eft,basis)
    assert np.allclose(mixed,.5*(parts[0]+parts[1]))
    assert np.linalg.norm(mixed-base)>0
    assert np.allclose(mixed,mixed.T)


def test_nonzero_mass_phase_does_not_affect_normal_channel():
    eft=FakeEFT(1600,2100)
    kernel=SingleScalarLineKernel()
    a=SelfBilinearMassBasis.build('S1',1,1600,50+60j)
    b=SelfBilinearMassBasis.build('S1',1,1600,50-60j)
    _,ca,_=stage.normal_propagator_c5(kernel,None,eft,a)
    _,cb,_=stage.normal_propagator_c5(kernel,None,eft,b)
    assert np.allclose(ca,cb)


def test_zero_self_mass_reproduces_original_kernel():
    eft=FakeEFT(1600,2100)
    basis=SelfBilinearMassBasis.build('S2',3,2100,0)
    original,corrected,parts=stage.normal_propagator_c5(SingleScalarLineKernel(),None,eft,basis)
    assert np.allclose(original,corrected,rtol=1e-14)
    assert np.allclose(parts[0],parts[1],rtol=1e-14)


@pytest.mark.parametrize('model,dim',[
    ('A_+0',1),('A_-2',3),('D_+0',3),('D_-2',1),('E_+0',3),('E_-2',3)])
def test_model_dependent_self_splitting_checks_saved_hessian(model,dim):
    species=stage.split_species(model)
    dims={'A':(1,3),'D':(3,1),'E':(3,3)}[model[0]]
    diag={'mS1Sq':1500,'mS2Sq':2000}
    rep={'d_s1':dims[0],'d_s2':dims[1]}
    mu=40+20j
    mass_sq=diag['mS1Sq' if species=='S1' else 'mS2Sq']
    other_dim=dims[1 if species=='S1' else 0]
    other_sq=diag['mS2Sq' if species=='S1' else 'mS1Sq']
    b=SelfBilinearMassBasis.build(species,dim,mass_sq,mu)
    spectrum=np.sort(np.r_[np.repeat(b.eigen_masses_squared[0],dim),
                            np.repeat(b.eigen_masses_squared[1],dim),
                            np.repeat(other_sq,2*other_dim)]).tolist()
    item={'quadratic_kind':'self_s1' if species=='S1' else 'self_s2',
          'eft1_scalar_endpoint_state':{'representation':rep,**diag},
          'quadratic_mu2_at_scalar_endpoint_gev2':{'real':mu.real,'imag':mu.imag},
          'scalar_endpoint_spectrum':{'mass_squared_eigenvalues_gev2':spectrum}}
    check=stage.candidate_basis(model,item)
    assert check.species==species and check.dim==dim
    item['scalar_endpoint_spectrum']['mass_squared_eigenvalues_gev2'][0]+=3
    with pytest.raises(ValueError,match='Hessian'):
        stage.candidate_basis(model,item)
