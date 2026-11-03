"""Synthetic evidence-chain and fail-closed tests; the genuine local run is separate."""
from __future__ import annotations

import copy
import csv
import hashlib
import json
from pathlib import Path

import pytest

from Numerical.lo_patterns_16.final_synthesis import (
    BASE, AUDIT, REPORT, SOURCES, MODELS, REFS, SCENARIOS, CG_TOTAL,
    CG_WEIGHT, CG_RATIO, FORMATS, SynthesisError, generate, sha, verify_inputs,
)


def write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def synthetic_tree(root: Path) -> Path:
    """Generate a complete *synthetic* 8-stage source tree, never UV data."""
    actual_template = Path(__file__).resolve().parents[1] / REPORT
    path = root / REPORT
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(actual_template.read_bytes())
    ref = [m for m in MODELS if m in REFS]
    fake_audit = {'status': 'charge_resolved_group_coefficients', 'count': 16,
                  'models': []}
    for m in MODELS:
        fam = m[0]
        n = 2 if fam == 'E' else 1
        fake_audit['models'].append({
            'class': fam,
            'model': {'alpha': int(m[2:])},
            'channels': [{'nonzero': True, 'complex_weight': CG_WEIGHT[fam],
                          'charge': str(i)} for i in range(n)],
            'total_group_coefficient': CG_TOTAL[fam]})
    write_json(root / AUDIT, fake_audit)
    aud_sha = sha(root / AUDIT)
    models = []
    for m in MODELS:
        values = {'m1_ev': .001, 'm2_ev': .01, 'm3_ev': .05}
        models.append({'model': m, 'family': m[0], 'status': 'conditional_uv_evaluated',
                       'physical_normalisation_certified': False,
                       'candidate_topology_charge_channel_count': 2 if m[0] == 'E' else 1,
                       **values})
        write_json(root / BASE / 'provisional_01/models' / f'{m}.json',
                   {'status': 'conditional_uv_evaluated', 'physical_normalisation_certified': False,
                    'source': {'audit_hashes': {
                        'charge_resolved_loop_weights.json': aud_sha}}})
    s = {'format': FORMATS['batch01'], 'physical_normalisation_certified': False,
         'counts': {'conditional_evaluated': 16, 'blocked': 0, 'physical_source_certified': 0},
         'models': models}
    write_json(root / SOURCES['batch01'], s)
    sh1 = sha(root / SOURCES['batch01'])
    write_json(root / SOURCES['batch02'],
               {'format': FORMATS['batch02'], 'physical_normalisation_certified': False,
                'counts': {'models': 16, 'source_scan_rows': 288},
                'models': [{'model': m} for m in MODELS], 'input_run_sha256': sh1})
    write_json(root / SOURCES['batch03'],
               {'format': FORMATS['batch03'], 'physical_normalisation_certified': False,
                'model_count': 16, 'models': [{'model': m} for m in MODELS],
                'source_summary_sha256': sh1})
    sh3 = sha(root / SOURCES['batch03'])
    m4 = [{'model': m, 'm1_ev': .001, 'm2_ev': .01, 'm3_ev': .05,
           'fermion_type': 'Dirac', 'physical_normalisation_certified': False,
           'lightest_mode': {'geometric_projection_fraction': .2,
                             'coherence_fraction': .25}} for m in MODELS]
    write_json(root / SOURCES['batch04'],
               {'format': FORMATS['batch04'], 'physical_normalisation_certified': False,
                'model_count': 16, 'models': m4,
                'source_sha256': {'source_run_sha256': sh1,
                                  'source_batch03_sha256': sh3}})
    sh4 = sha(root / SOURCES['batch04'])
    write_json(root / SOURCES['batch05'],
               {'format': FORMATS['batch05'], 'physical_normalisation_certified': False,
                'model_count': 16, 'models': [{'model': m} for m in MODELS],
                'source_manifest_sha256': {'study': sh1, 'batch03': sh3, 'batch04': sh4}})
    write_json(root / BASE / 'analysis_batch05_01/SOURCE_MANIFEST.json',
               {'format': 't3-chatA-batch05-source-manifest-v1',
                'physical_normalisation_certified': False})
    sh5 = sha(root / SOURCES['batch05'])
    mh5 = sha(root / BASE / 'analysis_batch05_01/SOURCE_MANIFEST.json')
    write_json(root / SOURCES['batch06'],
               {'format': FORMATS['batch06'], 'physical_normalisation_certified': False,
                'safe_for_physical_source_certification': False,
                'model_count': 16, 'models': MODELS,
                'source_sha256': {'study': sh1, 'batch05': sh5,
                                  'batch05_manifest': mh5}})
    sh6 = sha(root / SOURCES['batch06'])
    write_json(root / BASE / 'analysis_batch06_01/SOURCE_MANIFEST.json',
               {'format': 't3-chatA-batch06-manifest-v1',
                'physical_normalisation_certified': False})
    mh6 = sha(root / BASE / 'analysis_batch06_01/SOURCE_MANIFEST.json')
    write_json(root / BASE / 'analysis_batch06_01/HANDOFF.json',
               {'physical_normalisation_certified': False,
                'source_manifest_sha256': mh6})
    write_json(root / SOURCES['batch07'],
               {'format': FORMATS['batch07'], 'physical_normalisation_certified': False,
                'model_count': 16, 'models': MODELS,
                'source_sha256': {'batch06_factorization': sh6, 'batch06_manifest': mh6,
                                  'study': sh1, 'cg_weight_audit': aud_sha}})
    sh7 = sha(root / SOURCES['batch07'])
    write_json(root / BASE / 'analysis_batch07_01/SOURCE_MANIFEST.json',
               {'format': 't3-chatA-batch07-manifest-v1',
                'physical_normalisation_certified': False})
    mh7 = sha(root / BASE / 'analysis_batch07_01/SOURCE_MANIFEST.json')
    write_json(root / BASE / 'analysis_batch07_01/HANDOFF.json',
               {'physical_normalisation_certified': False, 'manifest_sha256': mh7})
    write_json(root / SOURCES['batch08'],
               {'format': FORMATS['batch08'], 'physical_normalisation_certified': False,
                'model_count': 16, 'models': MODELS,
                'source_sha256': {'batch06': sh6, 'batch06_manifest': mh6,
                                  'batch07': sh7, 'batch07_manifest': mh7, 'study': sh1}})
    write_json(root / BASE / 'analysis_batch08_01/SOURCE_MANIFEST.json',
               {'format': 't3-chatA-batch08-manifest-v1',
                'physical_normalisation_certified': False})
    write_json(root / BASE / 'analysis_batch08_01/HANDOFF.json',
               {'physical_normalisation_certified': False,
                'source_manifest_sha256': sha(root / BASE / 'analysis_batch08_01/SOURCE_MANIFEST.json')})

    base6 = root / BASE / 'analysis_batch06_01/batch06_coefficient_tests.csv'
    write_csv(base6,
              ['scenario', 'model', 'coefficient_residual', 'absolute_scale'],
              [{'scenario': scenario, 'model': model,
                'coefficient_residual': (.0002 if scenario == 'charge_split_diagnostic' else 1e-10),
                'absolute_scale': CG_RATIO[model[0]]}
               for scenario in SCENARIOS for model in MODELS])
    pairs = [('A_-4', 'D_+2'), ('A_-2', 'D_+0'), ('A_+0', 'D_-2'),
             ('B_-3', 'B_+1'), ('C_-3', 'C_+1'),
             ('E_-4', 'E_+2'), ('E_-2', 'E_+0')]
    cols = (['model_a', 'model_b', 'frozen_matrix_proportional_residual'] +
            [f'common_{ref}_matrix_residual' for ref in REFS])
    write_csv(root / BASE / 'analysis_batch05_01/batch05_matched_pairs.csv', cols,
              [dict(model_a=a, model_b=b,
                    frozen_matrix_proportional_residual=.1,
                    **{f'common_{ref}_matrix_residual': 1e-15 for ref in REFS})
               for a, b in pairs])
    write_csv(root / BASE / 'analysis_batch05_01/batch05_common_texture_controls.csv',
              ['model', 'texture_reference', 'm1_over_m3'],
              [dict(model=m, texture_reference=ref,
                    m1_over_m3=.0076484 if ref == 'C_+1' else .60203 if ref == 'C_-1' else .1487)
               for m in MODELS for ref in REFS])
    write_csv(root / BASE / 'analysis_batch07_01/batch07_first_corrections.csv',
              ['model', 'large_mix_cubic_relative_error', 'large_mix_linear_relative_error'],
              [{'model': m, 'large_mix_cubic_relative_error': 1e-8,
                'large_mix_linear_relative_error': 1e-5} for m in MODELS])
    write_csv(root / BASE / 'analysis_batch08_01/batch08_model_summary.csv',
              ['model', 'exact_remainder_order', 'transverse_signal_order',
               'endpoint_nonproportionality'],
              [{'model': m, 'exact_remainder_order': 2.01,
                'transverse_signal_order': '' if m[0] == 'A' else 1.002,
                'endpoint_nonproportionality': 1e-15 if m[0] == 'A' else .0002}
               for m in MODELS])
    write_csv(root / BASE / 'analysis_batch08_01/batch08_derivative_validation.csv',
              ['model', 'central_difference_slope_relative_error'],
              [{'model': m, 'central_difference_slope_relative_error': 1e-5}
               for m in MODELS])
    return root


@pytest.fixture
def root(tmp_path):
    return synthetic_tree(tmp_path)


def tweak_json(path: Path, func):
    data = json.loads(path.read_text(encoding='utf-8'))
    func(data)
    write_json(path, data)


def test_model_catalogue_and_weights():
    assert len(MODELS) == 16 and len(set(MODELS)) == 16
    assert set(CG_RATIO) == set('ABCDE')
    assert CG_RATIO['C'] == pytest.approx(1 / 3**.5)
    assert CG_RATIO['E'] == pytest.approx(2**.5)


def test_synthetic_all16_8_stage_evidence_is_accepted(root):
    d = verify_inputs(root)
    assert len(d['model_table']) == 16
    assert len(d['scenario_maxima']) == 5
    assert len(d['matched_pairs']) == 7
    assert d['physical_normalisation_certified'] is False


def test_generation_creates_only_fresh_files(root):
    before = {str(p): sha(p) for p in root.rglob('*') if p.is_file()}
    d = generate(root)
    output = root / BASE / 'analysis_final_synthesis_01'
    assert len(d['model_table']) == 16
    for name in ('FINAL_VERIFIED_TECHNICAL_REPORT.md', 'consolidated_evidence.json',
                 'SOURCE_MANIFEST.json', 'HANDOFF.md', 'HANDOFF.json',
                 'all16_verified.csv', 'scenario_maxima.csv', 'matched_pairs.csv'):
        assert (output / name).is_file()
    report = (output / 'FINAL_VERIFIED_TECHNICAL_REPORT.md').read_text(encoding='utf-8')
    assert 'C_+1' in report and '0/16' in report and 'factorization' in report
    assert all(sha(Path(k)) == v for k, v in before.items())
    handoff = json.loads((output / 'HANDOFF.json').read_text())
    assert handoff['physically_certified_count'] == 0


def test_output_overwrite_is_refused(root):
    generate(root)
    with pytest.raises(SynthesisError, match='populated'):
        generate(root)


def test_outside_tree_is_refused(root):
    with pytest.raises(SynthesisError, match='inside'):
        generate(root, output_dir=root / 'wrong')


def test_false_physical_promotion_blocked(root):
    tweak_json(root / SOURCES['batch04'],
               lambda d: d.update(physical_normalisation_certified=True))
    with pytest.raises(SynthesisError, match='certification'):
        generate(root)
    assert not (root / BASE / 'analysis_final_synthesis_01').exists()


def test_changed_original_study_blocked(root):
    tweak_json(root / SOURCES['batch01'],
               lambda d: d['models'][0].update(m3_ev=.04))
    with pytest.raises(SynthesisError, match='SHA'):
        generate(root)


def test_duplicate_model_blocked(root):
    tweak_json(root / SOURCES['batch01'],
               lambda d: d['models'][0].update(model=d['models'][1]['model']))
    with pytest.raises(SynthesisError):
        verify_inputs(root)


def test_missing_cg_channel_blocked(root):
    tweak_json(root / AUDIT, lambda d: d['models'][0].update(channels=[]))
    with pytest.raises(SynthesisError):
        verify_inputs(root)


def test_corrupted_charge_count_blocked(root):
    # Change model detail while leaving core JSON untouched.
    model = root / BASE / 'provisional_01/models/A_-4.json'
    tweak_json(model, lambda d: d['source']['audit_hashes'].update(
        {'charge_resolved_loop_weights.json': 'not-the-saved-sha'}))
    with pytest.raises(SynthesisError, match='CG audit SHA'):
        verify_inputs(root)


def test_missing_scenario_blocked(root):
    p = root / BASE / 'analysis_batch06_01/batch06_coefficient_tests.csv'
    data = list(csv.DictReader(p.open()))[:-1]
    write_csv(p, list(data[0]), data)
    with pytest.raises(SynthesisError, match='scenario'):
        verify_inputs(root)


def test_multiplet_factorization_lost_blocked(root):
    p = root / BASE / 'analysis_batch06_01/batch06_coefficient_tests.csv'
    data = list(csv.DictReader(p.open()))
    data[0]['coefficient_residual'] = '.1'
    write_csv(p, list(data[0]), data)
    with pytest.raises(SynthesisError, match='baseline'):
        verify_inputs(root)


def test_wrong_matched_texture_shape_blocked(root):
    p = root / BASE / 'analysis_batch05_01/batch05_matched_pairs.csv'
    data = list(csv.DictReader(p.open()))
    data[0]['common_C_+1_matrix_residual'] = '.01'
    write_csv(p, list(data[0]), data)
    with pytest.raises(SynthesisError, match='proportional'):
        verify_inputs(root)


def test_lost_cubic_improvement_blocked(root):
    p = root / BASE / 'analysis_batch07_01/batch07_first_corrections.csv'
    data = list(csv.DictReader(p.open()))
    data[0]['large_mix_cubic_relative_error'] = '.1'
    write_csv(p, list(data[0]), data)
    with pytest.raises(SynthesisError, match='Cubic'):
        verify_inputs(root)


def test_wrong_quadratic_order_blocked(root):
    p = root / BASE / 'analysis_batch08_01/batch08_model_summary.csv'
    data = list(csv.DictReader(p.open()))
    data[0]['exact_remainder_order'] = '1.01'
    write_csv(p, list(data[0]), data)
    with pytest.raises(SynthesisError, match='quadratic'):
        verify_inputs(root)


def test_wrong_linear_onset_blocked(root):
    p = root / BASE / 'analysis_batch08_01/batch08_model_summary.csv'
    data = list(csv.DictReader(p.open()))
    data[4]['transverse_signal_order'] = '2.5'
    write_csv(p, list(data[0]), data)
    with pytest.raises(SynthesisError, match='linear'):
        verify_inputs(root)


def test_derivative_validation_blocked(root):
    p = root / BASE / 'analysis_batch08_01/batch08_derivative_validation.csv'
    data = list(csv.DictReader(p.open()))
    data[0]['central_difference_slope_relative_error'] = '.01'
    write_csv(p, list(data[0]), data)
    with pytest.raises(SynthesisError, match='derivative'):
        verify_inputs(root)


def test_broken_batch07_manifest_chain(root):
    tweak_json(root / BASE / 'analysis_batch07_01/SOURCE_MANIFEST.json',
               lambda d: d.update(note='tampered'))
    with pytest.raises(SynthesisError, match='SHA'):
        verify_inputs(root)


def test_bad_csv_numeric_nonfinite_blocked(root):
    p = root / BASE / 'analysis_batch08_01/batch08_model_summary.csv'
    data = list(csv.DictReader(p.open()))
    data[1]['exact_remainder_order'] = 'NaN'
    write_csv(p, list(data[0]), data)
    with pytest.raises(SynthesisError, match='invalid number'):
        verify_inputs(root)


def test_malformed_certified_study_count_blocked(root):
    tweak_json(root / SOURCES['batch01'],
               lambda d: d['counts'].update(physical_source_certified=16))
    with pytest.raises(SynthesisError):
        verify_inputs(root)


@pytest.mark.skipif(not __import__('os').environ.get('RUN_LO_FINAL_SYNTHESIS_INTEGRATION'),
                    reason='Opt-in: needs full original Windows Chat A outputs')
def test_opt_in_real_checkout():
    repo = Path(__file__).resolve().parents[1]
    verify_inputs(repo)
