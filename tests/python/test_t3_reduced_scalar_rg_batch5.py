"""Offline Batch 5 tests: exact source algebra, RG generation orders, validation."""
from __future__ import annotations

from copy import deepcopy
import pytest

from scripts import AuditT3ReducedScalarRG as a


def test_expression_cancellation():
    assert a.expand_beta('2*lambdaT3 - 2*lambdaT3') == {}
    assert a.expand_beta('(x+y)^2-x^2-2*x*y-y^2') == {}


def test_rational_and_bar():
    p=a.expand_beta('Bar[lambdaT3]*Tr[y1 . Trans[Bar[y2]]] / 3')
    assert list(p.values()) == [a.Fraction(1,3)]
    assert any('BAR_lambdaT3' in m for m in p)
    assert any('__TR_' in ' '.join(m) for m in p)


def test_positive_and_negative_sources():
    p=a.expand_beta('3*q*x+lambdaT3*Tr[y1 . Trans[Bar[y2]]]')
    assert len(a.source_terms(p,{'q'}))==1
    assert a.source_terms(a.expand_beta('3*q*x + q^2'), {'q'})=={}


def test_barred_extra_zero():
    p=a.expand_beta('4*Bar[q]*mS1Sq+2*q*mS2Sq')
    assert a.source_terms(p,{'q'})=={}


def test_generation_chain_two_loops():
    betas={'q':a.expand_beta('lambdaT3*Tr[y1 . Trans[Bar[y2]]]'),
           'mu':a.expand_beta('6*q*mS1Sq+6*Bar[q]*mS2Sq')}
    assert a.minimum_orders(betas,{'q','mu'})=={'q':1,'mu':2}


def test_generation_chain_three_loops():
    betas={'q1':a.expand_beta('x'),
           'q2':a.expand_beta('q1*y'),
           'mu':a.expand_beta('q2*mS1Sq')}
    assert a.minimum_orders(betas,set(betas))=={'q1':1,'q2':2,'mu':3}


def test_no_source_at_any_finite_order():
    betas={'q':a.expand_beta('q*x'), 'mu':a.expand_beta('q*mS1Sq+mu')}
    assert a.minimum_orders(betas,set(betas))=={'q':None,'mu':None}


def test_nonlin_perturbation_order():
    betas={'q':a.expand_beta('x'),'r':a.expand_beta('q^2'),'mu':a.expand_beta('r*y')}
    assert a.minimum_orders(betas,set(betas))=={'q':1,'r':3,'mu':4}


def test_rational_division_reject_variable():
    with pytest.raises(ValueError,match='divides by a coupling'):
        a.expand_beta('q/x')


def test_unsupported_wolfram_rejected():
    with pytest.raises(ValueError):
        a.expand_beta('Sqrt[muT3GeneralS1Sq]')


def test_no_physical_mapping_claim():
    spec={'new_rgbeta_symbols':[], 'preexisting_rgbeta_symbols':[],
          'total_additional_quartics':0}
    x={'format':'general-t3-quadratic-rgbeta-experimental-v1',
       'status':'SuccessCandidateMassBeta','physical_mapping_certified':False,
       'metadata':{
           'GeneralScalarModel':'A_-4','GeneralScalarStage':'UV',
           'GeneralScalarBasisStatus':'CandidateNativeRGBetaUnmappedToMatchete',
           'QuadraticKind':'none','QuadraticMassSymbol':None,
           'PhysicalMatcheteRGBetaMassBasisCertified':False,
           'ReportBetaConvention':'16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g',
           'GeneralScalarNewBetaKeys':[], 'GeneralScalarAdditionalQuartics':0},
       'report_betas':{'lambdaT3':'2*lambdaT3'}}
    betas,fields=a.check_export('A_-4','UV',x,spec,False)
    assert betas=={} and fields['omitted']==[]
    y=deepcopy(x)
    y['physical_mapping_certified']=True
    with pytest.raises(ValueError,match='certificate invalid'):
        a.check_export('A_-4','UV',y,spec,False)


def test_unexpected_quadratic_source_fails():
    spec={'new_rgbeta_symbols':[], 'preexisting_rgbeta_symbols':[],
          'total_additional_quartics':0}
    x={'format':'general-t3-quadratic-rgbeta-experimental-v1',
       'status':'SuccessCandidateMassBeta','physical_mapping_certified':False,
       'metadata':{
           'GeneralScalarModel':'A_+0','GeneralScalarStage':'EFT1',
           'GeneralScalarBasisStatus':'CandidateNativeRGBetaUnmappedToMatchete',
           'QuadraticKind':'self_s1','QuadraticMassSymbol':'muT3GeneralS1Sq',
           'PhysicalMatcheteRGBetaMassBasisCertified':False,
           'ReportBetaConvention':'16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g',
           'GeneralScalarNewBetaKeys':[], 'GeneralScalarAdditionalQuartics':0},
       'report_betas':{'muT3GeneralS1Sq':'mS1Sq'}}
    with pytest.raises(ValueError,match='directly sourced'):
        a.check_export('A_+0','EFT1',x,spec,True)


def test_direct_registry_match():
    spec={'new_rgbeta_symbols':[], 'preexisting_rgbeta_symbols':[
        {'existing_name':'a','counts':[0,2,0,0,2,0]},
        {'existing_name':'b','counts':[0,2,0,2,0,0]}]}
    assert a._expected_direct_symbols('B_-1',spec)=={'a','b'}
    assert a._expected_direct_symbols('A_+0',spec)==set()


def test_missing_direct_registry_fails():
    with pytest.raises(ValueError,match='not represented'):
        a._expected_direct_symbols('B_-1',{'new_rgbeta_symbols':[],
                                            'preexisting_rgbeta_symbols':[]})


def test_frozen_log_factor():
    assert abs(a.math.log(100)/(16*a.math.pi*a.math.pi)-0.02916)<0.00001


def test_input_sha_changes(tmp_path):
    p=tmp_path/'input.json'
    p.write_text('a')
    h=a.fingerprint(p)
    p.write_text('b')
    assert a.fingerprint(p)!=h


def test_actual_saved_b_c_source_structure():
    # Representative exact terms from the checked B(-1) UV beta expressions.
    mu='muT3GeneralS1S2'
    q1='lambdaS1S1bar2S2bar'
    q2='lambdaS1barS2S2bar2'
    direct='lambdaHHdagS2S2'
    beta_direct=(f'-9*g2^2*{direct} + Bar[lambdaT3]*Tr[y1 . Trans[Bar[y2]]]'
                 f' + 4*{direct}*lambdaH2')
    beta_mu=(f'6*{q1}*mS1Sq + 6*{q2}*mS2Sq + '
             f'12*lambdaS1bar2S2bar2*Bar[{mu}]')
    assert len(a.source_terms(a.expand_beta(beta_direct), {direct,q1,q2})) == 1
    assert a.source_terms(a.expand_beta(beta_mu),{mu,q1,q2,'lambdaS1bar2S2bar2'}) == {}


def test_full_c5_matching_manifest_rejects_missing(tmp_path):
    with pytest.raises(ValueError, match='missing/duplicate'):
        a._validate_saved_matching('B_-1', {'results':[]})
