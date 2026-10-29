"""Static and synthetic-data checks for the new Stage 6 audit.

The local installation of Matchete and RGBeta is NOT simulated or claimed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.AssessGeneralT3C5Response import (  # noqa: E402
    ACTIVE, ACTIVE_MODELS, EXPECTED_MODELS, assess, classify_atom,
)


def _prepare_proof(root: Path, *, missing: str | None = None,
                   nonlinear: bool = False, zero_not_proven: bool = False) -> None:
    root.mkdir(parents=True, exist_ok=True)
    rows = []
    for model in sorted(EXPECTED_MODELS):
        if model == missing:
            continue
        active = model in ACTIVE_MODELS
        rows.append({
            'model': model,
            'status': 'ZeroExtraLimitProven' if active else 'TextIdentical',
            'zero_extra_limit_proven': not (zero_not_proven and model == 'B_m1'),
            'coefficient_independent_of_other_extra_couplings': not nonlinear,
            'linear_in_extra_coupling_atoms': not nonlinear,
            'response_reconstruction_residual_proven_zero': not nonlinear,
            'response_at_zero_proven_zero': True,
        })
        if active:
            artifact = root / model / 'general_c5_algebraic' / 'delta_c5_atoms.json'
            artifact.parent.mkdir(parents=True, exist_ok=True)
            artifact.write_text(json.dumps({
                'model': model,
                'response_reconstruction_residual_proven_zero': not nonlinear,
                'response_at_zero_proven_zero': True,
                'polynomial_in_extra_atoms': not nonlinear,
                'atoms': [
                    {
                        'coupling_atom': f'Coupling[{coupling}, {{}}, 0]',
                        'coefficient_formally_zero': False,
                        'linear_coefficient_inputform':
                            'hbar*Tr[y1 . Trans[Bar[y1]]]' + (
                                '+Coupling[lambdaT3ExtraQ020020I1, {}, 0]' if nonlinear
                                else ''),
                    }
                    for coupling in sorted(ACTIVE)
                ],
            }))
    (root / 'c5_algebraic_report.json').write_text(json.dumps({'results': rows}))


def test_active_models_and_names_are_explicit():
    assert len(EXPECTED_MODELS) == 16
    assert ACTIVE_MODELS == {'B_m1', 'C_m1'}
    assert ACTIVE == {'lambdaT3ExtraQ020020I1', 'lambdaT3ExtraQ020200I1'}


def test_formal_response_verified_for_complete_synthetic_report(tmp_path):
    _prepare_proof(tmp_path)
    report = assess(tmp_path)
    assert report['status'] == 'VerifiedFormalResponse'
    assert report['errors'] == []
    assert report['numeric_matching_implemented'] is False
    active = [m for m in report['results'] if m['extra_c5_derivatives']]
    assert len(active) == 2
    assert all(len(m['extra_c5_derivatives']) == 2 for m in active)


def test_missing_model_cannot_pass(tmp_path):
    _prepare_proof(tmp_path, missing='C_m1')
    report = assess(tmp_path)
    assert report['status'] == 'ReviewRequired'
    assert any('missing' in err.lower() for err in report['errors'])


def test_unproven_zero_limit_cannot_pass(tmp_path):
    _prepare_proof(tmp_path, zero_not_proven=True)
    assert assess(tmp_path)['status'] == 'ReviewRequired'


def test_nonlinear_derivative_cannot_pass(tmp_path):
    _prepare_proof(tmp_path, nonlinear=True)
    report = assess(tmp_path)
    assert report['status'] == 'ReviewRequired'
    assert any('nonlinear' in err for err in report['errors'])


def test_missing_coefficients_cannot_pass(tmp_path):
    _prepare_proof(tmp_path)
    (tmp_path / 'B_m1' / 'general_c5_algebraic' / 'delta_c5_atoms.json').unlink()
    assert assess(tmp_path)['status'] == 'ReviewRequired'


def test_classification_labels_pole_and_matching_log():
    obj = classify_atom({
        'coupling_atom': 'Bar[Coupling[lambdaT3ExtraQ020020I1, {}, 0]]',
        'linear_coefficient_inputform':
            'hbar*Bar[Coupling[y1, {}, 0]]*Log[μbar2/MF^2]/ϵ',
    })
    assert obj['holomorphic_or_conjugate'] == 'conjugate'
    assert obj['yukawa_species_in_coefficient'] == ['y1']
    assert obj['contains_regulator']
    assert obj['contains_matching_log']
    assert not obj['depends_on_other_extra_couplings']


def test_scripts_keep_basis_mapping_unverified_until_direct_check():
    tensor = (ROOT / 'scripts' / 'ExportGeneralT3TensorConvention.wl').read_text()
    rge = (ROOT / 'scripts' / 'AnalyzeGeneralT3RGBetaSources.wl').read_text()
    ps1 = (ROOT / 'scripts' / 'RunGeneralT3ScalarQuarticResponseAudit.ps1').read_text()
    assert '"rgbeta_mapping_verified" -> False' in tensor
    assert '"matchete_rgbeta_coupling_map_verified" -> False' in rge
    assert 'RunGeneralT3Analysis.ps1 -RequireProof' in ps1
    assert 'StringFreeQ' in rge
    assert 'MatcheteTensorConventionVerified' in tensor


def test_zero_derivative_is_not_a_discovered_contribution(tmp_path):
    _prepare_proof(tmp_path)
    path = tmp_path / 'B_m1' / 'general_c5_algebraic' / 'delta_c5_atoms.json'
    data = json.loads(path.read_text())
    data['atoms'].append({
        'coupling_atom': 'Coupling[lambdaT3ExtraQ020020I1, {}, 0]',
        'linear_coefficient_inputform': '0',
        'coefficient_formally_zero': True,
    })
    path.write_text(json.dumps(data))
    assert 'zero-response atoms' in ' '.join(assess(tmp_path)['errors'])


def test_incomplete_new_formal_proof_does_not_pass(tmp_path):
    _prepare_proof(tmp_path)
    path = tmp_path / 'B_m1' / 'general_c5_algebraic' / 'delta_c5_atoms.json'
    data = json.loads(path.read_text())
    data['response_reconstruction_residual_proven_zero'] = False
    path.write_text(json.dumps(data))
    assert assess(tmp_path)['status'] == 'ReviewRequired'


def test_stage6_refreshes_formal_artifacts_by_default():
    ps1 = (ROOT / 'scripts' / 'RunGeneralT3ScalarQuarticResponseAudit.ps1').read_text()
    wl = (ROOT / 'scripts' / 'CompareGeneralT3C5.wl').read_text()
    assert 'SkipFormalRefresh' in ps1
    assert 'Refresh formal C5 difference and response atoms' in ps1
    assert 'newExpr - GeneralC5ZeroExtra[newExpr]' in wl
    assert 'GeneralC5Describe[physicalDifference]' in wl
    assert 'response_reconstruction_residual_proven_zero' in wl


def test_mathematica_atom_checker_has_real_selftest():
    wl = (ROOT / 'scripts' / 'CompareGeneralT3C5.wl').read_text()
    ps = (ROOT / 'scripts' / 'RunGeneralT3ScalarQuarticResponseAudit.ps1').read_text()
    assert 'GeneralC5SelfTest[] :=' in wl
    assert 'args === {"SELFTEST"}' in wl
    assert '"nonlinear_rejected"' in wl
    assert 'Wolfram extra-coupling atom self-test' in ps


def test_wolframscript_file_mode_passes_selftest_as_positional_argument():
    """A dash-leading token can be consumed by wolframscript as an option."""
    runner = (ROOT / 'scripts' / 'RunGeneralT3ScalarQuarticResponseAudit.ps1').read_text()
    wl = (ROOT / 'scripts' / 'CompareGeneralT3C5.wl').read_text()
    import re
    call = re.search(
        r"Invoke-Required 'Wolfram extra-coupling atom self-test'"
        r"(?P<args>.*?)\n\)",
        runner,
        flags=re.DOTALL,
    )
    assert call is not None
    assert "'SELFTEST'" in call.group('args')
    assert "'--selftest'" not in call.group('args')
    assert 'args === {"SELFTEST"}' in wl
    assert 'Received script arguments:' in wl


def test_tensor_export_initialises_sm_gauge_registry_after_reset():
    """Regression for DefineRepresentation::unkwngroup (SU2L missing)."""
    wl = (ROOT / 'scripts' / 'ExportGeneralT3TensorConvention.wl').read_text()
    start = wl.index('T3TensorRun[] :=')
    block = wl[start:]
    reset = block.index('ResetAll[];')
    loaded = block.index('LoadModel["SM"]', reset)
    generated = block.index('T3TensorDescribe[model,', reset)
    assert reset < loaded < generated
    assert 'smBaseline === $Failed || smBaseline === $Aborted' in block
    assert 'SMModelLoadFailed' in block
    assert 'UsingFrontEnd[T3TensorRun[]]' in wl


def test_stage6_runner_supports_fast_tensor_only_retry():
    runner = (ROOT / 'scripts' / 'RunGeneralT3ScalarQuarticResponseAudit.ps1').read_text()
    assert '[switch]$TensorOnly' in runner
    assert 'if ($TensorOnly)' in runner
    assert "Invoke-Required 'Matchete doublet tensor convention'" in runner
    # Standalone tensor retry must exit before formal-recheck invocation.
    assert runner.index('if ($TensorOnly)') < runner.index("Invoke-Required 'Wolfram extra-coupling atom self-test'")
