"""Tests for Stage18 numerical extra scalar-pair source transport."""
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from Numerical.experimental.GeneralT3ExtraWilsonTransport import (
    compile_stage17_sources, frozen_seed_flavor_contractions,
    integrate_source_over_eft1, run_stage18_model,
)


def certified(model='B_m1'):
    # These are synthetic, not the user's actual Stage17 Matchete coefficients.
    suffixes=((1,1),(1,2),(2,2),(3,3))
    expr='(lambdaHHdagS2S2+conjugate(lambdaHHdagS2S2))*conjugate(y2)**2/MF + (lambdaHHdagS1barS1bar+conjugate(lambdaHHdagS1barS1bar))*conjugate(y1)**2/(2*MF)'
    return {'status':'ExtraEFT1ScalarPairMixingDerived','model':model,
            'source_convention':'lambda_Matchete=-lambda_RGBeta, same named scalar field basis',
            'full_eft1_wilson_running_implemented':False,
            'extra_low_energy_neutrino_mass':None,
            'master_rge_term':'delta_beta_ijab^(1) = (1/2) sum_swap(a,b),c,d delta_lambda_abcd C_ijcd',
            'nonzero_weinberg_real_components':4,
            'components':[{'higgs_real_indices':list(k),'beta_16pi2_expression':expr} for k in suffixes]}


def test_compile_certified_sources_and_mapping():
    cs=compile_stage17_sources('B_m1',certified())
    assert len(cs)==4
    assert cs[0].higgs_pair==(1,1)
    yy=np.ones((3,3),dtype=complex)
    matrix=cs[0].evaluate(2+1j,3+2j,yy,yy)
    assert np.allclose(matrix,7*yy)


def test_full_flavor_lift_symmetric_and_heavy_weighted():
    m=np.array([2.,4.,8.])
    y1=np.array([[1,0,0],[2,0,0],[3,0,0]],dtype=complex)
    y2=np.zeros((3,3),dtype=complex)
    x1,x2=frozen_seed_flavor_contractions(m,y1,y2)
    expected=np.outer([1,2,3],[1,2,3])/2.
    assert np.allclose(x1,expected)
    assert np.allclose(x2,0)
    assert np.allclose(x1,x1.T)


def test_complex_yukawas_are_conjugated_without_hermitian_transpose():
    m=[100.,200.,300.]
    y1=np.diag([1+2j,2-1j,3j]);y2=np.eye(3)
    got,_=frozen_seed_flavor_contractions(m,y1,y2)
    assert abs(got[0,0]-(1-2j)**2/100)<1e-15
    assert got[0,1]==0


def test_unmatched_model_rejected():
    with pytest.raises(ValueError,match='Stage17'):
        compile_stage17_sources('C_m1',certified('B_m1'))


def test_mixed_yukawa_monomials_refused():
    p=certified()
    p['components'][0]['beta_16pi2_expression']='conjugate(y1)*conjugate(y2)*lambdaHHdagS2S2/MF'
    with pytest.raises(ValueError,match='Yukawa monomial'):
        compile_stage17_sources('B_m1',p)


def test_nonlinear_quartic_rejected():
    p=certified()
    p['components'][0]['beta_16pi2_expression']='lambdaHHdagS2S2**2*conjugate(y2)**2/MF'
    with pytest.raises(ValueError,match='linear'):
        compile_stage17_sources('B_m1',p)


def test_unresolved_symbol_rejected():
    p=certified()
    p['components'][0]['beta_16pi2_expression']='unknown*lambdaHHdagS2S2*conjugate(y2)**2/MF'
    with pytest.raises(ValueError,match='Unsupported symbols'):
        compile_stage17_sources('B_m1',p)


def test_tree_seed_not_constant_in_extra_coupling():
    p=certified()
    p['components'][0]['beta_16pi2_expression']='conjugate(y2)**2/MF'
    with pytest.raises(ValueError,match='zero at zero'):
        compile_stage17_sources('B_m1',p)


class FakeEFT:
    def __init__(self,mus,l2,l1):
        self.mu_gev=np.asarray(mus,dtype=float)
        self.states=[SimpleNamespace(
            mu_gev=float(mu), lambdaHHdagS2S2=complex(l2[i]),
            lambdaHHdagS1barS1bar=complex(l1[i]),
            representation=SimpleNamespace(shared_scalar=False,d_s1=2,d_s2=2,alpha=-1))
                     for i,mu in enumerate(mus)]
    def state_at_index(self,k):return self.states[k]


def test_integrated_constant_source_sign_and_units():
    mu=np.geomspace(4000,1000,17)
    t=FakeEFT(mu,np.full(17,2.+0j),np.full(17,3.+0j))
    cs=compile_stage17_sources('B_m1',certified())
    y1=np.zeros((3,3),complex);y2=np.zeros((3,3),complex)
    y1[0,0]=1j;y2[0,0]=2
    out=integrate_source_over_eft1(cs,t,np.array([100.,200.,300.]),y1,y2)
    # f1=(lambda1+conj)/2=3, f2=(lambda2+conj)=4.
    # conj(i)^2/100=-.01 ; conj(2)^2/100=.04 -> beta=.13 GeV^-1.
    expected=.13*np.log(1000/4000)/(16*np.pi**2)
    got=out['higgs_real_components'][0]['integrated_delta_C_ijab_gev_inverse']
    assert got['real'][0][0]==pytest.approx(expected,rel=1e-12)
    assert got['imag'][0][0]==pytest.approx(0.,abs=1e-17)
    assert out['extra_physical_c5_at_mu_s'] is None
    assert out['extra_neutrino_mass_low_ev'] is None


def test_time_linear_source_simpson_exact():
    mus=np.geomspace(4000,1000,17); t=np.log(mus/mus[0])
    l2=2.0+3.0*t; l1=np.zeros(17)
    fake=FakeEFT(mus,l2,l1)
    cs=compile_stage17_sources('B_m1',certified())
    y=np.zeros((3,3),complex);y[0,0]=1
    out=integrate_source_over_eft1(cs,fake,np.array([100.,200.,300.]),np.zeros_like(y),y)
    lo=np.log(1000/4000)
    # beta = (4+6*t)/100, integral=(4*lo+3*lo^2)/100/16pi^2
    expected=(4*lo+3*lo**2)/(100*16*np.pi**2)
    got=out['higgs_real_components'][0]['integrated_delta_C_ijab_gev_inverse']['real'][0][0]
    assert got==pytest.approx(expected,rel=1e-12)


def test_insufficient_or_bad_scales_fail_closed():
    cs=compile_stage17_sources('C_m1',certified('C_m1'))
    y=np.eye(3)
    for mus in [np.geomspace(4000,1000,4),np.geomspace(4000,1000,6),np.array([4000,3000,3100,2000,1000])]:
        fake=FakeEFT(mus,np.ones(len(mus)),np.ones(len(mus)))
        with pytest.raises(ValueError,match='odd'):
            integrate_source_over_eft1(cs,fake,np.array([1.,2.,3.]),y,y)


def test_bad_nsave_rejected_early():
    with pytest.raises(ValueError,match='n_save'):
        run_stage18_model('B_m1',{}, {}, {},certified(), n_save=12)


def test_invalid_matrix_and_mass_rejected():
    with pytest.raises(ValueError,match='positive'):
        frozen_seed_flavor_contractions(np.array([0.,1.,2.]),np.eye(3),np.eye(3))
    with pytest.raises(ValueError,match='3x3'):
        frozen_seed_flavor_contractions(np.array([1.,2.,3.]),np.eye(2),np.eye(3))


def test_real_stage17_reports_if_available():
    import json
    for model in ('B_m1','C_m1'):
        path=ROOT/'output/general_t3_stage17'/f'{model}_extra_eft1_scalar_pair_beta.json'
        if not path.exists():
            pytest.skip('Actual Stage17 reports are only in user local repository')
        parsed=compile_stage17_sources(model,json.loads(path.read_text(encoding='utf-8-sig')))
        assert len(parsed)==4
        for c in parsed:
            x=c.evaluate(1e-12+2e-12j,3e-12-1e-12j,np.eye(3)*1e-15,np.eye(3)*1e-15)
            assert x.shape==(3,3)
