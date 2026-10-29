"""Exact source-based SU2 gauge covariance and counterexample regressions."""
from pathlib import Path
import copy
import json

import pytest
import sympy as sp

from scripts.RGBeta110101Covariance import (
    COUNTS, FIELDS, MODELS, EXPECTED, ORIGINAL_FIRST, ORIGINAL_REUSED,
    PROPOSED, _SIGMA, analyze, components, covariance, generator_residual,
)
from scripts.RGBetaTensorAlgebra import exact

INPUT = Path('output/general_t3_stage22')


def load():
    paths=[INPUT/n for n in ('matchete_tensors.json','rgbeta_native_invariants.json','rgbeta_native_declarations.json')]
    if not all(p.is_file() for p in paths):
        pytest.skip('Native Stage22 export evidence must be present for integration tests')
    return [json.loads(p.read_text(encoding='utf-8-sig')) for p in paths]


def test_native_first_is_su2_covariant_under_all_generators():
    result=covariance(EXPECTED[ORIGINAL_FIRST]); assert result['su2_covariant']
    assert [x['nonzero_residual_components'] for x in result['generator_checks'].values()]==[0,0,0]


def test_native_reused_breaks_su2_covariance():
    result=covariance(EXPECTED[ORIGINAL_REUSED]); assert not result['su2_covariant']
    assert [x['nonzero_residual_components'] for x in result['generator_checks'].values()]==[4,0,4]
    assert result['generator_checks']['T1']['first_nonzero'][0]['indices_1_based']==[1,1,1,2]


def test_proposed_replacement_is_covariant_and_not_zero():
    result=covariance(PROPOSED)
    assert result['su2_covariant']
    assert any(v != 0 for v in components(PROPOSED).values())
    assert [x['nonzero_residual_components'] for x in result['generator_checks'].values()]==[0,0,0]


def test_conjugation_rep_is_necessary_for_epsilon_diagnostic():
    native=components(EXPECTED[ORIGINAL_FIRST])
    q=_SIGMA[2]/2
    assert generator_residual(native,q)=={}
    # Alter the generator convention deliberately: four fundamental fields
    # is not the field content of this holomorphic structure.
    idx=(0,0,1,0)
    naive=sum(q[b,idx[ax]]*native[tuple(idx[:ax]+(b,)+idx[ax+1:])]
              for ax in range(4) for b in range(2))
    assert naive!=0


def test_both_models_exact_proposed_matrix_and_original_obstruction():
    result=analyze(*load())
    phase=result['native_source_phase']
    assert phase in ('PreFix','PostFix')
    assert result['counts']=={'models_checked':2,'native_symbols_checked':4,
                              'native_pass':4 if phase=='PostFix' else 2,
                              'native_fail':0 if phase=='PostFix' else 2,
                              'proposed_pass':2,'proposed_full_rank_geometry_blocks':2}
    assert not result['production_modified']
    assert not result['physical_coupling_mapping_certified']
    for m in MODELS:
        r=result['models'][m]
        assert r['original_rgbeta_columns_in_matchete_span']==([True,True] if phase=='PostFix' else [True,False])
        assert r['proposed_exact_geometry']['matchete_rank']==2
        assert r['proposed_exact_geometry']['rgbeta_rank']==2
        assert r['proposed_exact_geometry']['maximum_exact_component_residual']=='0'
        assert r['proposed_exact_geometry']['r_in_m']==[['1/2','1'],['-sqrt(3)/2','0']]
        assert sp.simplify(exact(r['proposed_matrix_determinant'])-sp.sqrt(3)/2)==0


def test_tampered_native_reused_fails_closed():
    mat,inv,decl=load()
    inv=copy.deepcopy(inv)
    bad=next(x for x in inv['models']['B_-1']['invariants'] if x['symbol']==ORIGINAL_REUSED)
    bad['native_invariant_expression_inputform']='eps[SU2L[fund],t3i1,t3i4]*eps[SU2L[fund],t3i2,t3i3]'
    with pytest.raises(ValueError, match='source invariant changed'):
        analyze(mat,inv,decl)


def test_tampered_field_order_fails_closed():
    mat,inv,decl=load()
    inv=copy.deepcopy(inv)
    r=next(x for x in inv['models']['C_-1']['invariants'] if x['symbol']==ORIGINAL_FIRST)
    r['native_fields_inputform']='{Bar[H], H, Bar[S1], Bar[S2]}'
    with pytest.raises(ValueError,match='incoherent fields_inputform'):
        analyze(mat,inv,decl)


def synthetic_postfix(mat, inv, dec):
    for model in MODELS:
        ir=next(x for x in inv['models'][model]['invariants'] if x['symbol']==ORIGINAL_REUSED)
        dr=next(x for x in dec['models'][model]['declarations'] if x['symbol']==ORIGINAL_REUSED)
        ir['native_invariant_expression_inputform']=PROPOSED
        corrected_function='del[SU2L[fund], #1, #3]*eps[SU2L[fund], #2, #4] & '
        ir['native_invariant_function_inputform']=corrected_function
        dr['invariant_inputform']=corrected_function
    return mat, inv, dec


def test_synthetic_postfix_native_pair_is_covariant_and_exact():
    mat,inv,dec=load()
    mat,inv,dec=synthetic_postfix(mat,copy.deepcopy(inv),copy.deepcopy(dec))
    report=analyze(mat,inv,dec)
    assert report['native_source_phase']=='PostFix'
    assert report['counts']['native_pass']==4
    assert report['counts']['native_fail']==0
    for model in MODELS:
        r=report['models'][model]
        assert r['native_reused']['su2_covariant']
        assert r['original_rgbeta_columns_in_matchete_span']==[True,True]
        assert r['proposed_exact_geometry']['r_in_m']==[['1/2','1'],['-sqrt(3)/2','0']]


def test_full_model_source_contains_only_reviewed_corrected_invariant():
    source=Path('RGE/running/rgbeta/T3RGBetaModel.wl')
    if not source.is_file():
        pytest.skip('Full RGBeta model file not available in this isolated test directory')
    text=source.read_text(encoding='utf-8')
    marker='lambdaHHdagS1barS2barCross,'
    import re
    matches=re.findall(r'AddQuartic\[\s*lambdaHHdagS1barS2barCross,([\s\S]*?)SelfConjugate',text)
    assert len(matches)==1
    block=matches[0]
    assert 'del[SU2L @ fund, #1, #3]' in block
    assert 'eps[SU2L @ fund, #2, #4]' in block
    assert 'eps[SU2L @ fund, #1, #3]' not in block
