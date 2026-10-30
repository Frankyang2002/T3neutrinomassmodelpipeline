"""Fast source-level regression tests; no Wolfram runtime required."""
from __future__ import annotations

from pathlib import Path

from scripts.AuditGeneralT3QuadraticBasis import predict, audit
from scripts.RunT3CompleteGeneralScalarSymbolic import (
    arguments, validate_betas, validate_build,
)

ROOT = Path(__file__).resolve().parents[1]


def test_complete_mass_inventory():
    result = audit()
    kinds = {k: v['holomorphic_mass_operators']
             for k, v in result['models'].items()}
    assert len(kinds) == 16
    assert sum(map(len, kinds.values())) == 8
    assert kinds['A_-2'] == ['self_s2']
    assert kinds['A_+0'] == ['self_s1']
    assert kinds['B_-1'] == ['mixed']
    assert kinds['C_-1'] == ['mixed']
    assert kinds['E_+2'] == []
    assert predict(3, 3, 0) == ['self_s1']


def test_lagrangian_source_is_opt_in_and_fail_closed():
    source = (ROOT / 'Lagrangian/interactions/T3Topology.wl').read_text()
    builder = (ROOT / 'Lagrangian/model/LagrangianBuilder.wl').read_text()
    assert 'T3_EXPERIMENTAL_GENERAL_QUADRATICS' in source
    assert 'T3ExtendedBuildCandidates[model, 2, structures]' in source
    assert 'T3ExtendedBuildExtraQuarticCandidates' in source
    assert 'RequiresMixedScalarMassDiagonalization' in builder
    assert 'RequiredGeneralScalarOperatorRejected' in builder
    assert 'AdditionalQuadraticsIncluded' in builder


def test_cli_representation_encoding():
    assert arguments('A_-2') == ('A', 'm2', ['1', '3', '2'])
    assert arguments('C_-1') == ('C', 'm1', ['2', '2', '3'])
    assert arguments('E_+2') == ('E', '2', ['3', '3', '2'])


def test_fails_closed_on_missing_bilinear():
    data = {
        'format': 't3-full-general-scalar-lagrangian-candidate-v1',
        'model': 'A_-2', 'status': 'Success',
        'quadratic_expected': 1, 'quadratic_included': 0,
        'quadratic_mass_basis_blocked': 0,
        'matching_mass_basis_ready': True, 'matching_performed': False,
        'physical_matchete_rgbeta_mapping_certified': False,
    }
    try:
        validate_build('A_-2', data, 1)
    except ValueError:
        pass
    else:
        raise AssertionError('Missing bilinear wrongly accepted')


def test_fails_closed_on_unresolved_mass_beta():
    payload = {
        'format': 'general-t3-quadratic-rgbeta-experimental-v1',
        'status': 'SuccessCandidateMassBeta',
        'physical_mapping_certified': False,
        'mass_eigenstate_threshold_matching_performed': False,
        'metadata': {
            'GeneralScalarModel': 'A_-2',
            'GeneralScalarStage': 'UV',
            'GeneralScalarQuadraticMassMixing': 'IncludedCandidate',
            'QuadraticKind': 'self_s2',
            'QuadraticMassSymbol': 'muT3GeneralS2Sq',
            'GeneralScalarAdditionalQuartics': 1,
            'GeneralScalarQuarticBetaBasis': 'RGBetaFullQuarticProjection',
            'PhysicalMatcheteRGBetaMassBasisCertified': False,
            'GeneralScalarNewBetaKeys': ['lambdaExtra'],
        },
        'beta_symbols': ['lambdaExtra', 'mS1Sq', 'mS2Sq', 'muT3GeneralS2Sq'],
        'report_betas': {
            'lambdaExtra': 'x', 'mS1Sq': 'x', 'mS2Sq': 'x',
            'muT3GeneralS2Sq': 'BetaTerm[muT3GeneralS2Sq,1]'
        },
    }
    reg = {'new_rgbeta_symbols': [{'name': 'lambdaExtra'}],
           'total_additional_quartics': 1}
    try:
        validate_betas('A_-2', 'UV', payload, reg, 'self_s2')
    except ValueError:
        pass
    else:
        raise AssertionError('Unresolved mass beta wrongly accepted')
