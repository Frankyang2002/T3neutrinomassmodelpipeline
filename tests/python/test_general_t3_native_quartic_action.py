"""Stage25 native assembly audit: synthetic fixtures plus optional live data."""
from copy import deepcopy
import json
from pathlib import Path

import pytest
from scripts.NativeQuarticActionAudit import FORMAT, EXPECTED_STATUS, verify, run

ROOT = Path(__file__).resolve().parents[2]


def synthetic():
    """Small mock invariant forms, 16 model rows, 50 mock slots, 46 blocks.

    This is test input only; NO live Wolfram or Matchete evidence is fabricated.
    """
    models = [f'm{i}' for i in range(16)]
    dm, lm, blocks = {}, {}, []
    symbols=[]
    for i, m in enumerate(models):
        n = 4 if i < 2 else 3  # 2*4 + 14*3 = 50
        ss=[f'lambda_{i}_{j}' for j in range(n)]
        symbols += [(m,s) for s in ss]
        dm[m]={'declarations':[
            {'symbol':s,'unique_arrangements_inputform':'{{1,2,3,4},{2,1,3,4}}'}
            for s in ss]}
        lm[m]=[
            {'symbol':s,'native_permutations':2,'native_assembly_rules_checked':2,
             'assembly_equals_independent_permutation_reconstruction':True,
             'complex_conjugate_metric_two_verified':True,
             'field_species_positions':[1,2,3,4],
             'native_coupling_value_inputform':f'6*{s}'}
            for s in ss]
    # 4 two-slot blocks -> 50 slots in 46 blocks
    pairs = {(0,0),(0,2),(1,0),(1,2)}
    for m in models:
        i=int(m[1:]); n=4 if i<2 else 3
        j=0
        while j < n:
            take=2 if (i,j) in pairs else 1
            ss=[f'lambda_{i}_{k}' for k in range(j,j+take)]
            blocks.append({'model':m,'rgbeta_couplings':ss,
                           'native_rgbeta_permutation_counts':{s:2 for s in ss}})
            j+=take
    l={'format':FORMAT,'status':EXPECTED_STATUS,'models_checked':16,'slots_checked':50,
       'all_native_permutation_reconstructions_exact':True,
       'all_complex_conjugate_metric_two_verified':True,
       'native_complex_species_to_cartesian_real_intertwiner_verified':False,
       'matchete_relative_potential_sign_verified':False,
       'physical_coupling_mapping_certified':False,'production_modified':False,
       'models':lm}
    d={'models':dm}
    t={'source_prefactors_verified':True,'signed_physical_coupling_map_certified':False,
       'totals':{'slots':50,'blocks':46},'blocks':blocks}
    return l,d,t


def test_mock_complete_assembly_coverage_but_never_certifies():
    report=verify(*synthetic())
    assert report['slots_checked']==50
    assert report['blocks_covered']==46
    assert report['native_permutation_rules_checked']==100
    assert report['physical_coupling_mapping_certified'] is False


def test_mock_missing_model_rejected():
    l,d,t=synthetic(); l['models'].pop('m0')
    with pytest.raises(ValueError,match='model sets differ'):
        verify(l,d,t)


def test_mock_incorrect_permutation_factor_rejected():
    l,d,t=synthetic(); l['models']['m1'][0]['native_assembly_rules_checked']=1
    with pytest.raises(ValueError,match='permutation count'):
        verify(l,d,t)


def test_mock_wrong_complex_metric_rejected():
    l,d,t=synthetic(); l['models']['m1'][0]['complex_conjugate_metric_two_verified']=False
    with pytest.raises(ValueError,match='unverified native'):
        verify(l,d,t)


def test_mock_wrong_native_factor_rejected():
    l,d,t=synthetic(); l['models']['m1'][0]['native_coupling_value_inputform']='12*Matrix[lambda_1_0][]'
    with pytest.raises(ValueError,match='coupling factor'):
        verify(l,d,t)


def test_unevaluated_matrix_form_is_not_the_runtime_coefficient():
    # RGBeta GroupsAndIndices.m: Matrix[m_][] = m.
    l,d,t=synthetic()
    l['models']['m1'][0]['native_coupling_value_inputform']='6*Matrix[lambda_1_0][]'
    with pytest.raises(ValueError, match='coupling factor'):
        verify(l,d,t)


def test_mock_unexpected_certification_rejected():
    for field in ('physical_coupling_mapping_certified','matchete_relative_potential_sign_verified',
                  'native_complex_species_to_cartesian_real_intertwiner_verified'):
        l,d,t=synthetic();l[field]=True
        with pytest.raises(ValueError,match='must not assert'):
            verify(l,d,t)


def test_mock_duplicate_symbol_rejected():
    l,d,t=synthetic();l['models']['m2'][1]['symbol']=l['models']['m2'][0]['symbol']
    with pytest.raises(ValueError,match='duplicate'):
        verify(l,d,t)


def test_mock_stage23_inconsistency_rejected():
    l,d,t=synthetic(); t['blocks'][0]['native_rgbeta_permutation_counts'][t['blocks'][0]['rgbeta_couplings'][0]]=3
    with pytest.raises(ValueError,match='Stage23/native'):
        verify(l,d,t)


def test_optional_live_integration():
    folder=ROOT/'output/general_t3_stage22'
    if not (folder/'stage25_live_native_quartic_assembly.json').is_file():
        pytest.skip('Requires live Stage25 Mathematica native-assembly export')
    result=run(folder,None)
    assert result['slots_checked']==50
    assert result['physical_coupling_mapping_certified'] is False
