"""Offline tests: no GitHub, Wolfram, fitted files or live RGBeta execution."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import json
import math
from types import SimpleNamespace

import numpy as np
import pytest

from scripts import RunT3ZeroExtraBoundaryBenchmark as b


def rep():
    return SimpleNamespace(d_s1=2,d_s2=2,d_f=1,alpha=-1)


def payload():
    names=['mS1Sq','mS2Sq','lambdaT3','lambdaHHdagS2S2',
           'lambdaHHdagS1barS1bar','lambdaHHdagS1barS2barCross',
           'lambdaS1bar2S2bar2','lambdaS1barS2S2bar2','lambdaS1S1bar2S2bar',
           'lambdaT3GeneralQ110101I1',b.MASS_SYMBOL]
    return {
        'format':b.BETA_FORMAT,'status':b.BETA_STATUS,
        'physical_mapping_certified':False,
        'mass_eigenstate_threshold_matching_performed':False,
        'beta_symbols':names,'report_betas':{k:'0' for k in names},
        'metadata':{
            'GeneralScalarModel':'B_-1','GeneralScalarStage':'UV',
            'GeneralScalarBasisStatus':'CandidateNativeRGBetaUnmappedToMatchete',
            'GeneralScalarQuadraticMassMixing':'IncludedCandidate',
            'QuadraticKind':'mixed','QuadraticMassSymbol':b.MASS_SYMBOL,
            'QuadraticPrefactorConvention':'V2 = mu2 * invariant + h.c.; no 1/2',
            'ReportBetaConvention':b.CONVENTION,
            'SharedScalar':False,'dS1':2,'dS2':2,'dF':1,'alpha':-1,
            'PhysicalMatcheteRGBetaMassBasisCertified':False,
            'GeneralScalarNewBetaKeys':['lambdaT3GeneralQ110101I1'],
            'GeneralScalarAdditionalQuartics':7,
            'GeneralScalarExistingStructuralCount':6,
        },
    }


def params():
    return dict(model='B_-1',stage='UV',rep=rep(),
                base_names=tuple(k for k in payload()['beta_symbols']
                                 if k!=b.MASS_SYMBOL and not k.startswith('lambdaT3GeneralQ')),
                new_names=('lambdaT3GeneralQ110101I1',),
                reused_names=('lambdaHHdagS2S2','lambdaHHdagS1barS1bar',
                              'lambdaHHdagS1barS2barCross','lambdaS1bar2S2bar2',
                              'lambdaS1barS2S2bar2','lambdaS1S1bar2S2bar'))


def snap(q=0j,mu=0j):
    reused={k:b.complex_json(0j) for k in params()['reused_names']}
    reused['lambdaHHdagS2S2']=b.complex_json(q)
    return {'new_quartics':{'lambdaT3GeneralQ110101I1':b.complex_json(0j)},
            'reused_quartics':reused,'quadratic_mu2_gev2':b.complex_json(mu),
            'mS1Sq_gev2':1.0e6,'mS2Sq_gev2':1.5e6,
            'lambdaT3':b.complex_json(.1)}


def test_beta_payload_complete():
    assert len(b.check_payload(payload(),**params()))==11


@pytest.mark.parametrize('key,value',[
    ('physical_mapping_certified',True),
    ('mass_eigenstate_threshold_matching_performed',True),
    ('status','Success'),
    ('format','nope'),
])
def test_reject_false_certification_and_wrong_format(key,value):
    v=payload();v[key]=value
    with pytest.raises(ValueError):b.check_payload(v,**params())


def test_reject_beta_missing():
    v=payload();v['report_betas'].pop(b.MASS_SYMBOL)
    with pytest.raises(ValueError,match='missing/extra'):b.check_payload(v,**params())


def test_reject_beta_head():
    v=payload();v['report_betas'][b.MASS_SYMBOL]='BetaTerm[foo]'
    with pytest.raises(ValueError,match='unresolved'):b.check_payload(v,**params())


def test_reject_wrong_stage():
    v=payload();v['metadata']['GeneralScalarStage']='EFT1'
    with pytest.raises(ValueError,match='metadata'):b.check_payload(v,**params())


def test_reject_wrong_representation():
    v=payload();v['metadata']['dF']=3
    with pytest.raises(ValueError,match='metadata'):b.check_payload(v,**params())


def test_reject_incomplete_registry():
    v=payload();v['metadata']['GeneralScalarAdditionalQuartics']=6
    with pytest.raises(ValueError,match='registry'):b.check_payload(v,**params())


def test_zero_boundary_exact():
    s=snap();b.check_zero_boundary(s,new=params()['new_names'],reused=params()['reused_names'])
    s['quadratic_mu2_gev2']=b.complex_json(1e-30)
    with pytest.raises(ValueError,match='not exactly zero'):
        b.check_zero_boundary(s,new=params()['new_names'],reused=params()['reused_names'])


def test_generation_and_ratios():
    initial=snap();final_uv=snap(q=0.001)
    final_eft=snap(q=0.003,mu=0j)
    final_eft['reused_quartics']['lambdaHHdagS1barS1bar']=b.complex_json(-0.001)
    d=b.inspect_generation(initial,final_uv,final_eft,
                           new=params()['new_names'],reused=params()['reused_names'])
    assert math.isclose(d['direct_C5_candidate_quartics']['lambdaHHdagS2S2']['absolute_ratio_to_lambdaT3'],0.03)
    assert d['direct_C5_candidate_quartics']['lambdaHHdagS1barS1bar']['map_certified'] is False
    assert d['quadratic_fraction_of_min_diagonal_m2']==0.0


def test_invalid_mass_basis():
    initial=snap();finish=snap()
    finish['mS1Sq_gev2']=-1
    with pytest.raises(ValueError,match='non-positive'):
        b.inspect_generation(initial,initial,finish,
                             new=params()['new_names'],reused=params()['reused_names'])


def test_c5_ratio():
    x=np.eye(3,dtype=complex)
    y=1.02*x
    assert math.isclose(b.relative_c5(x,y),0.02,rel_tol=1e-13)
    assert math.isclose(b.relative_c5(x,x),0.0)
    with pytest.raises(ValueError,match='vanishes'):
        b.relative_c5(np.zeros((3,3),complex),y)


def test_c5_matrix_validation():
    ok={'symmetrized_C5_gev_inv':[[b.complex_json(0j) for _ in range(3)] for _ in range(3)]}
    assert b.c5_matrix(ok,'test').shape==(3,3)
    ok['symmetrized_C5_gev_inv'][1][0]=b.complex_json(1j)
    with pytest.raises(ValueError,match='asymmetric'):b.c5_matrix(ok,'test')


def test_integration_without_physics_sources():
    # The RK solver helper is independently tested with an analytically
    # soluble affine 2-coupling toy system, not an invented RGBeta export.
    y0=np.array([0.,0.])
    rhs=lambda t,y:np.array([1.,y[0]])
    sol=b.integrate_segment(y0,rhs,100.,1.,n_points=8)
    logratio=math.log(.01)
    assert np.allclose(sol.y[:,-1],[logratio,0.5*logratio**2],rtol=1e-8,atol=1e-8)


def test_segment_guard():
    with pytest.raises(ValueError):b.integrate_segment(np.zeros(2),lambda t,y:y,1,100)


def test_read_write_and_hash(tmp_path):
    path=tmp_path/'x.json';b.save_json(path,{'x':3})
    assert b.read_json(path)=={'x':3}
    assert len(b.digest(path))==64
    assert b.beta_path('B_-1','UV',tmp_path).name=='B_m1_uv_quadratic_live_rgbeta.json'
    assert b.beta_path('C_-1','EFT1',tmp_path).name=='C_m1_eft1_quadratic_live_rgbeta.json'


def test_summary_no_physical_claim():
    r={'model':'B_-1','generation':{
        'direct_C5_candidate_quartics':{
           'lambdaHHdagS2S2':{'rgbeta_at_muS':b.complex_json(0.001)},
           'lambdaHHdagS1barS1bar':{'rgbeta_at_muS':b.complex_json(0.002)}},
        'quadratic_fraction_of_min_diagonal_m2':0.0},'C5_preview_performed':False}
    s=b.summary_md([r])
    assert 'not certified' in s and '0.001' in s and '0.002' in s
    assert 'NOT a neutrino-mass error bar' in s
