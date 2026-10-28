"""Stage24 exact SU2 covariance, source integrity, and source-phase-aware checks."""
import copy
from pathlib import Path

import pytest
import sympy as sp

from scripts.Stage24BCRemainderCovariance import (
    OLD, PROPOSED, FIELDS, MODELS, covariance, components, calculate,
    group_blocks, global_doublet_scan, run,
)
from scripts.VerifyGeneralT3Stage24LiveBC import verify


@pytest.mark.parametrize('name',sorted(PROPOSED))
def test_native_quartic_is_not_su2_covariant(name):
    audit = covariance(OLD[name],FIELDS[name])
    assert not audit['su2_covariant']
    assert audit['nonzero_residuals_T1_T2_T3'][0] > 0
    assert audit['nonzero_residuals_T1_T2_T3'][2] > 0


@pytest.mark.parametrize('name',sorted(PROPOSED))
def test_proposed_tensor_is_su2_covariant_and_nonzero(name):
    audit = covariance(PROPOSED[name],FIELDS[name])
    assert audit['su2_covariant']
    assert audit['nonzero_residuals_T1_T2_T3'] == [0,0,0]
    assert any(components(PROPOSED[name]).values())


def test_tensors_are_different_and_source_has_four_exact_replacements():
    assert len(PROPOSED) == 4 and set(PROPOSED) == set(OLD) == set(FIELDS)
    for sym in OLD:
        assert any(sp.simplify(components(OLD[sym])[k]-components(PROPOSED[sym])[k]) != 0
                   for k in components(OLD[sym]))


def test_live_verifier_rejects_incomplete_evidence():
    with pytest.raises(ValueError, match='Invalid Wolfram'):
        verify({'totals':{'old_non_covariant_slots':8},
                'physical_coupling_mapping_certified':False}, {})


def test_live_verifier_rejects_incorrect_native_component():
    wl = {'format':'stage24-live-rgbeta-bc-other-quartic-covariance-v1',
          'physical_coupling_mapping_certified':False,'models':{
          m:{s:{'native':{'su2_covariant':False,
                          'nonzero_residual_components_by_generator':[1,0,1],
                          'component_tensor_flat_inputform':['0']*16},
                'proposed':{'su2_covariant':True,
                            'nonzero_residual_components_by_generator':[0,0,0],
                            'component_tensor_flat_inputform':['0']*16}}
             for s in PROPOSED} for m in MODELS}}
    with pytest.raises(ValueError, match='component tensor mismatch'):
        verify({'totals':{'old_non_covariant_slots':8},
                'physical_coupling_mapping_certified':False}, wl)


def test_current_input_integration_if_present():
    folder = Path('output/general_t3_stage22')
    required = ['matchete_tensors.json','rgbeta_native_declarations.json',
                'rgbeta_native_invariants.json','stage22_aligned_geometry_diagnostics.json']
    if not all((folder/name).exists() for name in required):
        pytest.skip('Actual corrected Stage22 evidence is not available')
    try:
        result = run(folder)
    except ValueError as exc:
        if 'stale or unexpectedly changed source expression' in str(exc):
            pytest.skip('Historical Stage24 audit requires pre-four-invariant source')
        raise
    assert result['totals']['old_non_covariant_slots'] == 8
    assert all(not x['global_basis_scan']['common_global_basis_exists']
               for x in result['models'].values())


def test_synthetic_post_110101_evidence_has_eight_exact_candidate_matches(tmp_path):
    """Only a pre-Stage24 native export can exercise this historical witness."""
    import hashlib
    import json
    from scripts.Stage22RGBetaAlignedGeometry import analyze
    from scripts.Stage23QuarticPrefactorAudit import INPUTS, CORRECTED, CORRECTED_FUNCTION
    from scripts.Stage23QuarticPrefactorAudit import compact
    src = Path('output/general_t3_stage22')
    if not all((src/f).is_file() for f in INPUTS[:3]):
        pytest.skip('No saved Stage22 native evidence for synthetic regression fixture')
    data = {f:json.loads((src/f).read_text(encoding='utf-8-sig')) for f in INPUTS[:3]}
    # This test deliberately inspects the *pre*-Stage24 quartet.  After the
    # genuine native correction there is no historical OLD witness to audit.
    current = [row['native_invariant_expression_inputform']
               for m in MODELS
               for row in data[INPUTS[2]]['models'][m]['invariants']
               if row['symbol'] in OLD]
    if len(current) != 8:
        raise AssertionError('Expected eight B/C reused-invariant slots')
    if any(compact(expr) != OLD[name]
           for m in MODELS
           for name in OLD
           for expr in [next(x['native_invariant_expression_inputform']
                             for x in data[INPUTS[2]]['models'][m]['invariants']
                             if x['symbol'] == name)]):
        pytest.skip('Native input is post-Stage24: pre-correction synthetic fixture is inapplicable')
    for model in MODELS:
        for row in data[INPUTS[1]]['models'][model]['declarations']:
            if row['symbol'] == 'lambdaHHdagS1barS2barCross':
                row['invariant_inputform'] = CORRECTED_FUNCTION.replace('&',' & ')
        for row in data[INPUTS[2]]['models'][model]['invariants']:
            if row['symbol'] == 'lambdaHHdagS1barS2barCross':
                row['native_invariant_expression_inputform'] = CORRECTED
                row['native_invariant_function_inputform'] = CORRECTED_FUNCTION.replace('&',' & ')
    for f,v in data.items():
        (tmp_path/f).write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')
    geometry = analyze(data[INPUTS[0]],data[INPUTS[2]],data[INPUTS[1]])
    geometry['input_sha256'] = {
        f:hashlib.sha256((tmp_path/f).read_bytes()).hexdigest() for f in INPUTS[:3]}
    (tmp_path/INPUTS[3]).write_text(json.dumps(geometry,indent=2)+'\n',encoding='utf-8')
    result = run(tmp_path)
    assert result['totals'] == {
        'models':2,'old_non_covariant_slots':8,'proposed_covariant_slots':8,
        'proposed_exact_matchete_1d_matches':8,
        'restricted_global_I_or_J_basis_exists':False}
    assert all(g['global_basis_scan']['choices'] == 8 for g in result['models'].values())
    assert all(g['global_basis_scan']['maximum_aligned_blocks'] == 5
               for g in result['models'].values())
    altered = copy.deepcopy(data[INPUTS[2]])
    for x in altered['models']['B_-1']['invariants']:
        if x['symbol'] == 'lambdaS1bar2S2bar2':
            x['native_invariant_expression_inputform'] = PROPOSED['lambdaS1bar2S2bar2']
    with pytest.raises(ValueError, match='stale or unexpectedly changed'):
        calculate(data[INPUTS[0]],altered,
                  {'totals':{'slots':50,'blocks':46},
                   'signed_physical_coupling_map_certified':False})
