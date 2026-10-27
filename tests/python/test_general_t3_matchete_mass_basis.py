"""Offline regression for the installed Matchete DetailedOutput baseline adapter.

The real Wolfram kernel is NOT mocked as running; Python-side report contracts,
failure handling and source semantics are tested. Synthetic native values mirror
the *shape*, not the contents, of the user-observed opaque Matchete values.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
SOURCE = ROOT / 'scripts/RunGeneralT3QuadraticMatcheteAudit.py'
SPEC = importlib.util.spec_from_file_location('t3_mass_diag_cli', SOURCE)
assert SPEC and SPEC.loader
cli = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cli)

ALL = ('A_+0', 'A_-2', 'A_-4', 'B_+1', 'B_-1', 'B_-3',
       'C_+1', 'C_-1', 'C_-3', 'D_+0', 'D_+2', 'D_-2',
       'E_+0', 'E_+2', 'E_-2', 'E_-4')
MIXED = {'B_-1', 'C_-1'}
MASS = {'A_+0', 'A_-2', 'B_-1', 'C_-1', 'D_+0', 'D_-2', 'E_+0', 'E_-2'}
QUARTICS = {'A_+0': 5, 'A_-2': 5, 'A_-4': 0, 'B_+1': 1, 'B_-1': 7, 'B_-3': 1,
            'C_+1': 1, 'C_-1': 7, 'C_-3': 1, 'D_+0': 5, 'D_+2': 0, 'D_-2': 5,
            'E_+0': 6, 'E_+2': 0, 'E_-2': 6, 'E_-4': 0}
OPAQUE = {'Charges': 'MatcheteNativeCharges[value]',
          'GaugeAnomalies': 'MatcheteNativeAnomalies[value]',
          'ZeroModeMixing': 'MatcheteNativeZeroModes[value]'}


def diag(*, mixed=False, opaque=None):
    opaque = dict(OPAQUE if opaque is None else opaque)
    return {
        'diagnostic_returned': True,
        'standard_check_passed': not mixed,
        'standard_check_boolean': True,
        'non_mass_basis_checks_passed': True,
        'heavy_mass_basis_reported': True,
        'heavy_mass_basis_passed': not mixed,
        'opaque_metadata_baseline_match': True,
        'opaque_metadata_values_inputform': opaque,
        'failed_checks': [],
        'comparison_method': 'ExplicitCoreChecksPlusOrdinaryBaselineAnchoredOpaqueMetadata',
        'diagnostics': {
            'Open/Complex/InconsistentSpinChains': 'PassEmpty',
            'Hermiticity': 'PassTrue',
            'UncontractedIndices': 'PassEmpty',
            'CanonicallyNormalized': 'PassTrue',
            'HeavyMassBasis': 'FailFalse' if mixed else 'PassTrue',
            'HeavyTadpoles': 'PassEmpty',
            'GaugeInvariant': 'PassTrue',
            'AllObjectsDefined': 'PassTrue',
            **{name: 'OpaqueMatchesValidBaseline' for name in opaque},
        },
    }


def fake_report(model):
    cls, alpha = model.split('_')
    mixed = model in MIXED
    quartics = [f'T3ExtraQ_{model}_{i}' for i in range(QUARTICS[model])]
    quadratics = [f'T3ExtraM_{model}'] if model in MASS else []
    return {
        'status': cli.NON_MASS_PASS if mixed else cli.HEAVY_PASS,
        'class': cls, 'alpha': int(alpha),
        'quadratic_slots': len(quadratics),
        'additional_quartics': len(quartics),
        'baseline_heavy_lagrangian_valid': True,
        'combined_non_mass_basis_checks_passed': True,
        'baseline_calibrated_diagnostic_parser': True,
        'opaque_diagnostics_uninterpreted_but_baseline_identical': True,
        'combined_heavy_mass_basis_passed': not mixed,
        'checked_in_original_heavy_field_definitions': True,
        'operator_consistency_validated_except_mass_basis': mixed,
        'matchete_full_lagrangian_valid': not mixed,
        'mass_basis_diagonalized': False,
        'heavy_scalar_matching_mass_basis_usable': False,
        'mixed_heavy_mass_term_requires_diagonalization': mixed,
        'operator_normalization_certified': False,
        'rgbeta_physical_mapping_certified': False,
        'mass_eigenbasis_and_matching_done': False,
        'matching_executed': False,
        'quartic_operators': quartics,
        'quadratic_operators': quadratics,
        'quartic_individual_checks': [
            {'name': x, 'valid': True, 'diagnostic': diag()} for x in quartics],
        'quadratic_individual_checks': [
            {'name': x, 'valid': True, 'diagnostic': diag(mixed=mixed)} for x in quadratics],
        'baseline_diagnostic': diag(),
        'combined_diagnostic': diag(mixed=mixed),
    }


@pytest.mark.parametrize('model', ALL)
def test_all_16_models_anchor_validated(model):
    data = cli._validate_report(model, fake_report(model), int(model in MASS), QUARTICS[model])
    assert data['status'] == (cli.NON_MASS_PASS if model in MIXED else cli.HEAVY_PASS)


@pytest.mark.parametrize('model', ['B_-1', 'C_-1'])
@pytest.mark.parametrize('flag', [
    'status', 'combined_heavy_mass_basis_passed',
    'matchete_full_lagrangian_valid', 'mass_basis_diagonalized',
    'operator_normalization_certified', 'matching_executed',
    'rgbeta_physical_mapping_certified', 'heavy_scalar_matching_mass_basis_usable',
    'baseline_calibrated_diagnostic_parser',
    'opaque_diagnostics_uninterpreted_but_baseline_identical',
])
def test_mixed_unjustified_certificate_rejected(model, flag):
    result = fake_report(model)
    if flag == 'status':
        result[flag] = cli.HEAVY_PASS
    else:
        result[flag] = not result[flag]
    with pytest.raises(ValueError):
        cli._validate_report(model, result, 1, 7)


@pytest.mark.parametrize('model', ['A_+0', 'B_-1', 'C_-1'])
@pytest.mark.parametrize('field', [
    'Hermiticity', 'GaugeInvariant', 'HeavyMassBasis',
    'CanonicallyNormalized', 'UncontractedIndices',
    'Open/Complex/InconsistentSpinChains', 'HeavyTadpoles',
    'AllObjectsDefined', 'GaugeAnomalies',
])
def test_core_failures_and_missing_native_metadata_rejected(model, field):
    result = fake_report(model)
    obj = result['combined_diagnostic']
    if field == 'GaugeAnomalies':
        del obj['opaque_metadata_values_inputform'][field]
        del obj['diagnostics'][field]
    else:
        obj['diagnostics'][field] = 'UnexpectedValue'
    with pytest.raises(ValueError):
        cli._validate_report(model, result, int(model in MASS), QUARTICS[model])


@pytest.mark.parametrize('field', ['Charges', 'GaugeAnomalies', 'ZeroModeMixing'])
def test_changed_opaque_native_values_are_rejected_even_when_marked_as_passing(field):
    data = fake_report('B_-1')
    data['combined_diagnostic']['opaque_metadata_values_inputform'][field] = 'CHANGED'
    with pytest.raises(ValueError, match='native opaque diagnostic values differ'):
        cli._validate_report('B_-1', data, 1, 7)


def test_non_boolean_baseline_opaque_values_are_legitimate():
    result = fake_report('A_+0')
    assert cli._validate_report('A_+0', result, 1, 5)
    assert result['baseline_diagnostic']['diagnostics']['Charges'] == 'OpaqueMatchesValidBaseline'


def test_missing_native_opaque_dictionary_fails():
    result = fake_report('B_-1')
    del result['baseline_diagnostic']['opaque_metadata_values_inputform']
    with pytest.raises(ValueError):
        cli._validate_report('B_-1', result, 1, 7)


def test_mixed_plain_checker_must_return_false():
    result = fake_report('B_-1')
    result['quadratic_individual_checks'][0]['diagnostic']['standard_check_passed'] = True
    with pytest.raises(ValueError, match='standard_check_passed'):
        cli._validate_report('B_-1', result, 1, 7)


def test_mass_diagonal_basis_fails_when_no_quartic_list():
    row = fake_report('A_+0')
    row['quartic_operators'] = []
    with pytest.raises(ValueError):
        cli._validate_report('A_+0', row, 1, 5)


def test_alpha_encoding_and_report_validation(monkeypatch, tmp_path):
    basis = {'models': {x: {'independent_complex_mass_parameters_in_declared_basis': int(x in MASS)}
                        for x in ALL}}
    monkeypatch.setattr(cli, 'audit', lambda: basis)
    def mock_run(argv, **kwargs):
        assert argv[-2:] == ['B', 'm1']
        Path(argv[3]).write_text(json.dumps(fake_report('B_-1')))
        return SimpleNamespace(stdout='ok', stderr='', returncode=0)
    monkeypatch.setattr(cli.subprocess, 'run', mock_run)
    result = cli.run('B_-1', tmp_path, 'wolframscript')
    assert result['status'] == cli.NON_MASS_PASS
    assert (tmp_path / 'B_m1_matchete.log').is_file()


def test_matchete_executable_failure_does_not_pass(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, 'audit', lambda: {'models': {x: {} for x in ALL}})
    monkeypatch.setattr(cli.subprocess, 'run', lambda *a, **k:
                        SimpleNamespace(stdout='failed', stderr='', returncode=5))
    with pytest.raises(RuntimeError, match='exit=5'):
        cli.run('B_-1', tmp_path, 'wolframscript')


def test_wolfram_baseline_anchor_covers_live_changed_keys():
    src = (ROOT / 'Lagrangian/RunGeneralQuadraticMatcheteAudit.wl').read_text(encoding='utf-8')
    assert 'If[StringQ[key], Return[key]]' in src
    assert 'CheckLagrangian[lag, DetailedOutput -> True]' in src
    assert 'CheckLagrangian[lag]' in src
    assert 'opaqueBaseline = baseline["opaque_metadata_values_inputform"]' in src
    assert 'OpaqueMetadataDiffersFromValidBaseline' in src
    assert 'opaqueStrings, opaqueBaseline' in src
    assert 'ChargeNeutralOrChargesMissing' in src
    assert '"standard_check_passed" -> TrueQ[ordinary]' in src


def test_full_summary_discriminates_14_and_2(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, 'audit', lambda: {'models': {x: {} for x in ALL}})
    monkeypatch.setattr(cli, 'run', lambda m, _o, _e: fake_report(m))
    monkeypatch.setattr(sys, 'argv', ['module', '--all', '--output-dir', str(tmp_path)])
    cli.main()
    summary = json.loads((tmp_path / 'matchete_16_model_summary.json').read_text())
    assert (summary['models'], summary['quartics'], summary['quadratic_slots']) == (16, 50, 8)
    assert (summary['full_matchete_checks'], summary['non_mass_basis_only_checks']) == (14, 2)
    assert summary['mixed_scalar_threshold_matching_done'] is False
