"""Stage 23 exact source-prefactor identity, including post-Stage24 source phase."""
import copy
import hashlib
import json
from pathlib import Path

import pytest
import sympy as sp

from scripts.Stage23QuarticPrefactorAudit import (
    CORRECTED, CORRECTED_FUNCTION, INPUTS, calculate, exact_matrix,
    parse_permutations, run, scalar_component_prefactor,
)


def test_rgbeta_complex_quartic_factor_cancelled_by_symmetric_real_tensor():
    assert sp.factorial(4) / sp.sqrt(2)**4 == 6
    for n in (1,2,3,4,6,8,12,24):
        assert scalar_component_prefactor(6,n) == sp.Rational(1,4)
    with pytest.raises(ValueError):
        scalar_component_prefactor(24,12)


def test_native_arrangements_parser_and_unique_permutations():
    assert parse_permutations('{{1, 2, 3, 4}, {2, 1, 3, 4}}') == [(1,2,3,4),(2,1,3,4)]
    assert parse_permutations('{{1,2,3,4}}') == [(1,2,3,4)]
    for invalid in ('{{1,1,3,4}}','{{1,2,3,4},{1,2,3,4}}',
                    '{{1,2,3,4}};Exit[]','{}'):
        with pytest.raises(ValueError):
            parse_permutations(invalid)


def test_signed_basis_map_exact_symbolically():
    A = exact_matrix([['1/2','1'],['-sqrt(3)/2','0']])
    S = -A.inv()
    assert [[str(sp.simplify(S[i,j])) for j in range(2)] for i in range(2)] == [
        ['0','2*sqrt(3)/3'],['-1','-sqrt(3)/3']]
    assert A*S == -sp.eye(2)
    with pytest.raises(ValueError):
        exact_matrix([['1','0'],['0','0']])
    with pytest.raises(ValueError):
        exact_matrix([['__import__("os")']])


@pytest.fixture(scope='module')
def corrected_local_fixture(tmp_path_factory):
    """Synthetic 110101 correction on saved native exports; not a live re-export."""
    source = Path('output/general_t3_stage22')
    if not all((source/f).is_file() for f in INPUTS[:3]):
        pytest.skip('Native Stage22 evidence not present for integration fixture')
    from scripts.Stage22RGBetaAlignedGeometry import analyze
    folder = tmp_path_factory.mktemp('stage23_synthetic_postfix')
    data = {f:json.loads((source/f).read_text(encoding='utf-8-sig')) for f in INPUTS[:3]}
    for model in ('B_-1','C_-1'):
        for item in data[INPUTS[1]]['models'][model]['declarations']:
            if item['symbol'] == 'lambdaHHdagS1barS2barCross':
                item['invariant_inputform'] = CORRECTED_FUNCTION.replace('&',' & ')
        for item in data[INPUTS[2]]['models'][model]['invariants']:
            if item['symbol'] == 'lambdaHHdagS1barS2barCross':
                item['native_invariant_expression_inputform'] = CORRECTED
                item['native_invariant_function_inputform'] = CORRECTED_FUNCTION.replace('&',' & ')
    for f,content in data.items():
        (folder/f).write_text(json.dumps(content,indent=2)+'\n',encoding='utf-8')
    geometries = analyze(data[INPUTS[0]],data[INPUTS[2]],data[INPUTS[1]])
    geometries['input_sha256'] = {
        f:hashlib.sha256((folder/f).read_bytes()).hexdigest() for f in INPUTS[:3]}
    (folder/INPUTS[3]).write_text(json.dumps(geometries,indent=2)+'\n',encoding='utf-8')
    return folder


def test_conditional_50_couplings_and_46_blocks(corrected_local_fixture):
    d = run(corrected_local_fixture)
    fixed = {'models':16, 'slots':50, 'blocks':46,
             'DirectGeometryCandidate':32, 'CalibratedTripletCGCandidate':6,
             'CorrectedReusedInvariantCandidate':4}
    historical = {**fixed, 'PseudorealReindexCandidate':8}
    corrected_stage24 = {**fixed, 'CorrectedBCReusedInvariantCandidate':8}
    assert d['totals'] in (historical, corrected_stage24)
    assert len(d['blocks']) == 46
    assert d['source_prefactors_verified'] is True
    assert d['signed_physical_coupling_map_certified'] is False
    for b in d['blocks']:
        assert b['physical_mapping_certified'] is False
    bc = [x for x in d['blocks'] if x['model'] in ('B_-1','C_-1')
          and x['block'] == '110101']
    assert len(bc) == 2
    for x in bc:
        assert x['lambda_rgbeta_from_lambda_matchete_if_assumptions_hold'] == [
            ['0','2*sqrt(3)/3'],['-1','-sqrt(3)/3']]


def test_original_broken_invariant_fails_closed(corrected_local_fixture):
    data = {f:json.loads((corrected_local_fixture/f).read_text(encoding='utf-8'))
            for f in INPUTS}
    for x in data[INPUTS[2]]['models']['B_-1']['invariants']:
        if x['symbol'] == 'lambdaHHdagS1barS2barCross':
            x['native_invariant_expression_inputform'] = (
                'eps[SU2L[fund],t3i1,t3i3]*eps[SU2L[fund],t3i2,t3i4]')
    with pytest.raises(ValueError,match='pre-fix or stale'):
        calculate(data)


def test_wrong_native_factor_and_symmetry_average_fail_closed(corrected_local_fixture):
    data = {f:json.loads((corrected_local_fixture/f).read_text(encoding='utf-8'))
            for f in INPUTS}
    x = data[INPUTS[1]]['models']['A_+0']['declarations'][0]
    x['coupling_factor_inputform'] = x['coupling_factor_inputform'].replace('6*','24*')
    with pytest.raises(ValueError,match='exactly 6'):
        calculate(data)
    data = {f:json.loads((corrected_local_fixture/f).read_text(encoding='utf-8'))
            for f in INPUTS}
    data[INPUTS[1]]['models']['A_+0']['declarations'][0][
        'unique_arrangements_inputform'] = '{{1,1,3,4}}'
    data[INPUTS[2]]['models']['A_+0']['invariants'][0][
        'native_unique_arrangements_inputform'] = '{{1,1,3,4}}'
    with pytest.raises(ValueError,match='Non-unique or non-permutation'):
        calculate(data)


def test_geometry_source_hash_must_be_current(corrected_local_fixture):
    data = json.loads((corrected_local_fixture/INPUTS[3]).read_text(encoding='utf-8'))
    data['input_sha256'][INPUTS[1]] = '0'*64
    (corrected_local_fixture/INPUTS[3]).write_text(
        json.dumps(data)+'\n',encoding='utf-8')
    with pytest.raises(ValueError,match='hashes stale'):
        run(corrected_local_fixture)
