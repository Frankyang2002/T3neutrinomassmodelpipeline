"""Stage 29 exact convention-dependent boundary tests; no synthetic live claims."""
from copy import deepcopy
import json
from pathlib import Path

import pytest
import sympy as sp

from scripts.Stage29ConventionUVHandoff import (
    EXACT_MODELS, build, exact_input, source_symbols,
)


def fixtures():
    models = sorted(EXACT_MODELS)
    inventory = {m: {'total_additional_quartics': 0,
                     'new_rgbeta_symbols': [], 'preexisting_rgbeta_symbols': []}
                 for m in models}
    blocks = []
    seq = 0
    # 4 two-dimensional blocks + 42 one-dimensional blocks = 46 blocks / 50 slots.
    for n in range(46):
        model = models[n % 16]
        width = 2 if n < 4 else 1
        matchete, rgbeta = [], []
        for k in range(width):
            name = ('lambdaExisting' if seq < 12 else 'lambdaT3GeneralQ') + f'{seq:06d}I1'
            m_name = f'lambdaT3ExtraQ{seq:06d}I1'
            if seq < 12:
                inventory[model]['preexisting_rgbeta_symbols'].append({'existing_name':name})
            else:
                inventory[model]['new_rgbeta_symbols'].append({'name':name})
            inventory[model]['total_additional_quartics'] += 1
            rgbeta.append(name); matchete.append(m_name); seq += 1
        if width == 2:
            A = sp.Matrix([[sp.Rational(1,2), 1], [-sp.sqrt(3)/2, 0]])
        else:
            A = sp.Matrix([[1]])
        B = -A.inv()
        cls = 'CalibratedTripletCGCandidate' if n < 6 else 'DirectGeometryCandidate'
        blocks.append({
            'model': model, 'block': f'test-block-{n:02d}',
            'geometry_class': cls,
            'matchete_couplings': matchete, 'rgbeta_couplings': rgbeta,
            'r_in_m': [[str(A[i,j]) for j in range(width)] for i in range(width)],
            'lambda_rgbeta_from_lambda_matchete_if_assumptions_hold':
                [[str(sp.simplify(B[i,j])) for j in range(width)] for i in range(width)],
            'physical_mapping_certified': False,
            'assumed_relative_potential_sign':'-1',
        })
    assert seq == 50
    reg = {
        'totals': {'matchete_extra_quartics':50,'rgbeta_new_vertices':38,
                   'rgbeta_preexisting_structural':12},
        'models':inventory,
    }
    s23 = {'totals':{'slots':50,'blocks':46},'blocks':blocks,
           'signed_physical_coupling_map_certified':False}
    s26 = {'models_checked':16,'slots_checked':50,'blocks_checked':46,
           'signed_polynomial_identity_conditional_verified':True,
           'physical_coupling_mapping_certified':False,
           'blocks':[{'model':b['model'],'block':b['block'],
                      'exact_signed_action_identity_in_provisional_basis':True}
                     for b in blocks]}
    s27 = {'source_native_matchete_tensor_count':6,
           'su2_commutators_casimir_hermiticity_intertwiner_exact':True,
           'physical_coupling_mapping_certified':False}
    s28 = {'models_checked':6, 'rgbeta_delS2_orthogonality_verified':True,
           'matchete_triplet_invariant_metric_and_generators_exact':True,
           'native_rgbeta_numerical_triplet_generators_exported':False,
           'physical_coupling_mapping_certified':False}
    reused={v['existing_name'] for p in inventory.values() for v in p['preexisting_rgbeta_symbols']}
    source={'state':reused|{'lambdaH1','lambdaT3'},
            'vector':reused|{'lambdaH1','lambdaT3'},
            'evaluator':reused|{'lambdaH1','lambdaT3'}}
    return s23,s26,s27,s28,reg,source


def compute(*args):
    data = fixtures() if not args else args
    return build(*data, sp.Rational(1,100))


def test_16_models_50_slots_46_blocks_exact():
    r=compute()
    assert r['models_checked']==16 and r['slots_checked']==50 and r['blocks_checked']==46
    assert r['new_quartic_slots']==38 and r['reused_quartic_slots']==12
    assert sum(len(x['rgbeta_quartics']) for x in r['models'].values())==50
    assert r['block_dimension_counts']=={'1':42,'2':4}
    assert r['conditional_block_reconstruction_verified']


def test_nontrivial_signed_two_by_two_is_not_naive_sign_flip():
    r=compute()
    first=r['models'][sorted(EXACT_MODELS)[0]]['blocks'][0]
    assert len(first['signed_map_exact'])==2
    block=r['models'][sorted(EXACT_MODELS)[0]]
    rr=[sp.sympify(block['rgbeta_quartics'][s]['exact']) for s in first['rgbeta_symbols']]
    a=sp.Matrix([[sp.Rational(1,2),1],[-sp.sqrt(3)/2,0]])
    assert all(sp.simplify(z)==0 for z in a*sp.Matrix(rr)+sp.ones(2,1)/100)
    assert rr != [-sp.Rational(1,100)]*2


def test_numerical_values_and_no_false_rge_claims():
    r=compute()
    assert r['max_abs_conditional_rgbeta_quartic_uv']>0
    assert r['physical_coupling_mapping_certified'] is False
    assert r['one_loop_beta_functions_for_new_quartics_evaluated'] is False
    assert r['rge_evolution_performed'] is False
    assert r['class_E_topology_lambda_t3_factor_applied_to_additional_quartics'] is False
    assert r['production_ode_symbol_coverage']['state']['new_quartic_slots_recognised']==0
    assert r['production_ode_symbol_coverage']['state']['reused_quartic_slots_recognised']==12


@pytest.mark.parametrize('literal', ['0', '-0.01', '1.1','nan','inf','sqrt(2)','1e9999'])
def test_reject_invalid_quartic_input(literal):
    with pytest.raises(ValueError,match='Stage29 fail-closed'):
        exact_input(literal)


def test_decimal_exact_preserved():
    assert exact_input('0.01')==sp.Rational(1,100)
    assert exact_input('2.5e-2')==sp.Rational(1,40)


@pytest.mark.parametrize('where,value',[
    ('23',True),('26',True),('27',True),('28',True),
])
def test_refuse_premature_certificate(where,value):
    parts=list(fixtures())
    obj={'23':(0,'signed_physical_coupling_map_certified'),
         '26':(1,'physical_coupling_mapping_certified'),
         '27':(2,'physical_coupling_mapping_certified'),
         '28':(3,'physical_coupling_mapping_certified')}[where]
    parts[obj[0]][obj[1]]=value
    with pytest.raises(ValueError,match='Stage29 fail-closed'):
        compute(*parts)


def test_reject_modified_2x2_sign():
    parts=list(fixtures())
    parts[0]['blocks'][0]['lambda_rgbeta_from_lambda_matchete_if_assumptions_hold'][0][0]='1'
    with pytest.raises(ValueError,match='signed coupling matrix'):
        compute(*parts)


def test_reject_unregistered_or_missing_quartic_slot():
    parts=list(fixtures())
    parts[4]['models'][sorted(EXACT_MODELS)[0]]['new_rgbeta_symbols'].append({'name':'lambdaT3GeneralQ999999I1'})
    with pytest.raises(ValueError,match='Registry model total'):
        compute(*parts)


def test_reject_missing_26_block_identity():
    parts=list(fixtures())
    parts[1]['blocks'][0]['exact_signed_action_identity_in_provisional_basis']=False
    with pytest.raises(ValueError,match='Stage26 block'):
        compute(*parts)


def test_inspect_production_state_from_ast_without_execution(tmp_path:Path):
    paths={p:tmp_path/(p+'.py') for p in ('state','vector','evaluator')}
    paths['state'].write_text('class T3UVState:\n    lambdaS1: float\n    lambdaT3: complex\n    lambdaH1: float\n    lambda12: float\n    lambdaExisting: float\n')
    paths['vector'].write_text('v="lambdaExisting"\n')
    paths['evaluator'].write_text('v="lambdaExisting"\n')
    s=source_symbols(paths)
    assert 'lambdaExisting' in s['state']
    assert 'lambdaExisting' in s['vector']
    assert 'lambdaExisting' in s['evaluator']


def test_json_report_serializable():
    r=compute()
    assert json.loads(json.dumps(r))['models_checked']==16
