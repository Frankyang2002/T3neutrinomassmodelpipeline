"""Tests for evidence schema and Mathematica probe source (not a Wolfram run)."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from scripts.VerifyT3BCMatchete110101Action import (
    EXPECTED, FORMAT, SOURCE, verify, verify_source,
)

ROOT = Path(__file__).resolve().parents[1]


def fixture() -> dict:
    models = {}
    for name in ('B_-1', 'C_-1'):
        models[name] = {
            'model': name, 'status': 'LiveMatchete110101ActionTermsValidated',
            'matchete_action_validated': True, 'native_model_built': True,
            'scalar_extra_quartics_count': 7,
            'selected_110101_couplings_count': 2,
            'combined_check_lagrangian_passed': True,
            'combined_check_lagrangian_return_inputform': 'True',
            'source_operator_uses_plus_hc': 'verified separately from GeneralScalarOperators.wl',
            'cross_tool_action_equality_verified': False,
            'physical_coupling_mapping_certified': False,
            'physical_matching_certified': False,
            'selected_terms': [{
                'name': term, 'coupling': coupling,
                'coupling_present_in_native_expression': True,
                'check_lagrangian_passed': True,
                'check_lagrangian_return_inputform': 'True',
                'native_expression_inputform': f'PlusHc[{coupling}[] H[1] ...]',
            } for term, coupling in EXPECTED],
        }
    return {'format': FORMAT,
            'status': 'NativeMatcheteInteractionConstructionDiagnostic',
            'source': SOURCE, 'models': models,
            'matchete_scalar_off_diagonal_mass_included': False,
            'cross_tool_action_equality_verified': False,
            'physical_coupling_mapping_certified': False,
            'physical_matching_certified': False}


def test_fixture_is_only_schema_not_runtime():
    data = verify(fixture())
    assert data['matchete_native_terms_validated'] == 4
    assert data['matchete_combined_checks_passed'] == 2
    assert data['cross_tool_action_equality_verified'] is False


@pytest.mark.parametrize('change', [
    lambda d: d['models']['B_-1'].update({'combined_check_lagrangian_passed': False}),
    lambda d: d['models']['C_-1']['selected_terms'][0].update({'check_lagrangian_passed': False}),
    lambda d: d['models']['B_-1'].update({'scalar_extra_quartics_count': 6}),
    lambda d: d['models']['C_-1'].update({'physical_coupling_mapping_certified': True}),
    lambda d: d.update({'cross_tool_action_equality_verified': True}),
    lambda d: d['models']['B_-1']['selected_terms'][1].update({'coupling': 'incorrect'}),
    lambda d: d['models']['B_-1']['selected_terms'][0].update({'native_expression_inputform': '0'}),
    lambda d: d['models']['C_-1'].update({'status': 'Matchete110101ActionBlocked'}),
])
def test_bad_evidence_rejected(change):
    data = deepcopy(fixture())
    change(data)
    with pytest.raises(ValueError):
        verify(data)


def test_current_matchete_source_convention():
    verify_source(ROOT / 'Lagrangian/interactions/GeneralScalarOperators.wl')


def test_probe_source_contains_actual_native_calls():
    src = (ROOT/'scripts/ProbeT3BCMatchete110101Action.wl').read_text(encoding='utf-8')
    for token in ('UsingFrontEnd[Needs["Matchete`"]',
                  'T3ExtendedBuildExtraQuarticCandidates[model]',
                  'CheckLagrangian[',
                  'T3ModelFromClass[cls,-1]',
                  '{k,1,2}',
                  '"cross_tool_action_equality_verified"->False'):
        assert token in src


def test_m_all_native_operations_run_inside_front_end():
    """Matchete FrontEnd lifetime must cover runtime, not only package load."""
    import re
    source = (ROOT / 'scripts/ProbeT3BCMatchete110101Action.wl').read_text(encoding='utf8')
    assert re.search(
        r'result\s*=\s*CheckAbort\s*\[\s*Check\s*\[\s*UsingFrontEnd\s*\[',
        source, re.S,
    )
    scope = source.split('result=CheckAbort[Check[', 1)[1]
    assert scope.index('UsingFrontEnd[') < scope.index('T3BCMOne["B_-1","B"]')
    assert scope.index('UsingFrontEnd[') < scope.index('T3BCMOne["C_-1","C"]')
    assert '"models"-><|' in scope


def test_no_frontend_body_regression():
    """Removing v2's outer scope must make the regression gate fail."""
    import re
    source = (ROOT / 'scripts/ProbeT3BCMatchete110101Action.wl').read_text(encoding='utf8')
    hypothetical = source.replace('  UsingFrontEnd[\n  <|"format"', '  <|"format"')
    required = r'result\s*=\s*CheckAbort\s*\[\s*Check\s*\[\s*UsingFrontEnd\s*\['
    assert re.search(required,source,re.S)
    assert not re.search(required,hypothetical,re.S)
