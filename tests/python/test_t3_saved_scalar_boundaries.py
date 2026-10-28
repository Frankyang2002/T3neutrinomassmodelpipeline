"""Offline synthetic regression for the 48-case saved-study scalar audit."""
import csv
import json
from pathlib import Path

import pytest

from scripts import AuditT3SavedStudyScalarBoundaries as a

# Model inventory is synthetic, but matches the 50/38/12/8 registry totals.
COUNTS = {
    'A_+0': 5, 'A_-2': 5, 'A_-4': 0,
    'B_+1': 1, 'B_-1': 7, 'B_-3': 1,
    'C_+1': 1, 'C_-1': 7, 'C_-3': 1,
    'D_+0': 5, 'D_+2': 0, 'D_-2': 5,
    'E_+0': 6, 'E_+2': 0, 'E_-2': 6, 'E_-4': 0,
}
DIMENSIONS = {'A': [1, 3, 2], 'B': [2, 2, 1], 'C': [2, 2, 3],
              'D': [3, 1, 2], 'E': [3, 3, 2]}
BILINEAR = {'A_+0', 'A_-2', 'B_-1', 'C_-1', 'D_+0', 'D_-2', 'E_+0', 'E_-2'}


def write(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding='utf-8')


def synthetic_repository(root: Path):
    registry = {
        'status': 'GeneratedCandidateBasisNotWolframVerified',
        'totals': dict(a.REQUIRED_TOTALS), 'models': {},
    }
    for model in a.MODEL_NAMES:
        d = DIMENSIONS[model[0]]
        alpha = int(model.split('_')[1])
        reused_count = 6 if model in ('B_-1', 'C_-1') else 0
        registry['models'][model] = {
            'model': model, 'dimensions': d, 'alpha': alpha,
            'new_rgbeta_symbols': [
                {'name': f'lambdaT3GeneralQTest{j}_{model[0]}{abs(alpha)}'}
                for j in range(COUNTS[model] - reused_count)],
            'preexisting_rgbeta_symbols': [
                {'existing_name': f'lambdaTestReused_{model[0]}_{j}'}
                for j in range(reused_count)],
            'total_additional_quartics': COUNTS[model],
            'quadratics_not_implemented': int(model in BILINEAR),
        }
    write(root / a.REGISTRY, registry)
    for case, spec in a.STUDIES.items():
        comparisons = []
        for model in a.MODEL_NAMES:
            entry = registry['models'][model]
            filename = a.model_file(entry)
            stem = filename[:-5]
            expected = spec['lambda_ad'] * (0.5 if model.startswith('E') else 1.0)
            ordinary = {'lambdaT3': {'real': expected, 'imag': 0},
                        'lambdaS1': 0.1, 'lambdaH1': 0.02}
            for item in entry['preexisting_rgbeta_symbols']:
                ordinary[item['existing_name']] = {'real': 0.0, 'imag': 0.0}
            config = {
                'representation': {'d_s1': entry['dimensions'][0],
                                   'd_s2': entry['dimensions'][1],
                                   'd_f': entry['dimensions'][2],
                                   'alpha': entry['alpha'], 'shared_scalar': False},
                'base_state': {'mu_gev': 1e7, 'ordinary': ordinary},
                'scales': {'mu_fermion_threshold_gev': 1e5,
                           'mu_scalar_threshold_gev': 1e3,
                           'mu_low_gev': 1e2},
                'lambda_5e9_study': {},
            }
            fit = {'model': stem, 'status': 'Success', 'lambdaT3': expected,
                   'chi2_fit': 0.01}
            if spec['case']:
                config['fermion_start_study'] = {
                    'case': spec['case'], 'lambdaT3': expected,
                    'boundary_definition': 'UV initial parameter values imposed directly at mu_F; UV RGEs disabled',
                }
                comparisons.append({'model': stem, 'status': 'Success',
                                    'lambdaT3': expected, 'chi2': 0.01})
            write(root / spec['root'] / 'fitted' / filename, config)
            write(root / spec['root'] / 'fits' / filename, fit)
        if spec['case']:
            write(root / spec['root'] / 'comparison/optimised/summary.json', {
                'mode': 'optimised', 'case': spec['case'], 'uv_running': False,
                'start_scale': 'mu_F', 'models': comparisons,
            })
    return registry


@pytest.fixture
def repo(tmp_path):
    synthetic_repository(tmp_path)
    return tmp_path


def test_registry_coverage(repo):
    inv = a.registry_inventory(json.loads((repo / a.REGISTRY).read_text()))
    assert len(inv) == 16
    assert sum(len(v[2]) for v in inv.values()) == 38
    assert sum(len(v[3]) for v in inv.values()) == 12
    assert sum(v[4] for v in inv.values()) == 8


def test_complete_audit(repo):
    r = a.build_audit(repo)
    assert r['status'] == 'ExplicitReusedQuarticsZeroAcrossAllCases'
    assert set(r['cases']) == set(a.STUDIES)
    assert len(r['source_sha256']) == 99
    for case in r['cases'].values():
        assert case['models_checked'] == 16
        assert case['passes_chi2_screen'] == 16
        assert case['explicit_zero_reused_slots'] == 12
        assert case['new_slots_unrepresented'] == 38
        assert case['quadratic_slots_unrepresented'] == 8
        assert case['nonzero_reused_slots'] == 0


def test_no_false_explicit_zero_for_new_fields(repo):
    cases = a.build_audit(repo)['cases']
    for case in cases.values():
        b = next(x for x in case['models'] if x['model'] == 'B_-1')
        assert b['explicitly_zero_reused_count'] == 6
        assert b['new_quartics_unrepresented_count'] == 1
        assert b['quadratic_slots_unrepresented'] == 1
    assert all(case['new_slots_unrepresented'] == 38 for case in cases.values())


def test_read_only_historical_inputs(repo):
    before = (repo / a.STUDIES['large_uv_frozen_0p5']['root'] / 'fits' /
              'T3_dS1_2_dS2_2_dF_1_alpha_m1.json').read_bytes()
    r = a.build_audit(repo)
    a.write_audit(repo, r)
    after = (repo / a.STUDIES['large_uv_frozen_0p5']['root'] / 'fits' /
             'T3_dS1_2_dS2_2_dF_1_alpha_m1.json').read_bytes()
    assert after == before
    assert (repo / a.OUTPUT / 'restricted_input_audit.md').is_file()
    assert not (repo / 'output/studies/general_scalar_extended').exists()


def test_csv_48_rows_and_hash_keys_portable(repo):
    r = a.build_audit(repo)
    a.write_audit(repo, r)
    with (repo / a.OUTPUT / 'restricted_input_audit.csv').open(newline='') as f:
        records = list(csv.DictReader(f))
    assert len(records) == 48
    assert all('\\' not in p for p in r['source_sha256'])
    assert all(len(h) == 64 for h in r['source_sha256'].values())


def test_status_poorfit_preserved(repo):
    name = 'T3_dS1_2_dS2_2_dF_1_alpha_m1.json'
    file = repo / a.STUDIES['large_uv_frozen_0p5']['root'] / 'fits' / name
    fit = json.loads(file.read_text())
    fit.update(status='PoorFit', chi2_fit=32.29)
    write(file, fit)
    comp = repo / a.STUDIES['large_uv_frozen_0p5']['root'] / 'comparison/optimised/summary.json'
    data = json.loads(comp.read_text())
    next(x for x in data['models'] if x['model'] == name[:-5])['chi2'] = 32.29
    write(comp, data)
    r = a.build_audit(repo)
    row = next(x for x in r['cases']['large_uv_frozen_0p5']['models'] if x['model'] == 'B_-1')
    assert row['fit_status'] == 'PoorFit' and not row['passes_chi2_screen']
    assert r['cases']['large_uv_frozen_0p5']['passes_chi2_screen'] == 15


def test_detect_nonzero_structural_extra(repo):
    name = 'T3_dS1_2_dS2_2_dF_1_alpha_m1.json'
    f = repo / a.STUDIES['small_uv_running_5e-9']['root'] / 'fitted' / name
    d = json.loads(f.read_text())
    d['base_state']['ordinary']['lambdaTestReused_B_0']['real'] = 0.002
    write(f, d)
    r = a.build_audit(repo)
    assert r['status'] == 'NonzeroReusedQuarticsFound'
    assert r['cases']['small_uv_running_5e-9']['nonzero_reused_slots'] == 1


def test_reject_unsupported_new_quartic_in_old_state(repo):
    name = 'T3_dS1_2_dS2_2_dF_1_alpha_m1.json'
    f = repo / a.STUDIES['small_uv_running_5e-9']['root'] / 'fitted' / name
    d = json.loads(f.read_text())
    d['base_state']['ordinary']['lambdaT3GeneralQTest0_B1'] = {'real': 0, 'imag': 0}
    write(f, d)
    with pytest.raises(ValueError, match='new general quartic appears'):
        a.build_audit(repo)


def test_reject_unexpected_quadratic_in_old_state(repo):
    name = 'T3_dS1_2_dS2_2_dF_1_alpha_m1.json'
    f = repo / a.STUDIES['small_uv_running_5e-9']['root'] / 'fitted' / name
    d = json.loads(f.read_text())
    d['base_state']['ordinary']['muT3GeneralS1S2'] = {'real': 0, 'imag': 0}
    write(f, d)
    with pytest.raises(ValueError, match='candidate general quadratic'):
        a.build_audit(repo)


def test_reject_wrong_lambda(repo):
    f = repo / a.STUDIES['large_uv_frozen_0p5']['root'] / 'fitted' / 'T3_dS1_2_dS2_2_dF_1_alpha_m1.json'
    d = json.loads(f.read_text())
    d['base_state']['ordinary']['lambdaT3']['real'] = 5e-9
    write(f, d)
    with pytest.raises(ValueError, match='lambdaT3'):
        a.build_audit(repo)


def test_reject_mismatched_fit_chi2(repo):
    f = (repo / a.STUDIES['small_uv_frozen_5e-9']['root'] /
         'fits/T3_dS1_2_dS2_2_dF_1_alpha_m1.json')
    d = json.loads(f.read_text())
    d['chi2_fit'] = 23.0
    write(f, d)
    with pytest.raises(ValueError, match='comparison.chi2'):
        a.build_audit(repo)


def test_reject_wrong_frozen_boundary_marker(repo):
    f = (repo / a.STUDIES['small_uv_frozen_5e-9']['root'] /
         'fitted/T3_dS1_2_dS2_2_dF_1_alpha_m1.json')
    d = json.loads(f.read_text())
    d['fermion_start_study']['case'] = '0p5'
    write(f, d)
    with pytest.raises(ValueError, match='frozen fermion-start provenance'):
        a.build_audit(repo)


def test_reject_modified_registry_totals(repo):
    f = repo / a.REGISTRY
    d = json.loads(f.read_text())
    d['totals']['matchete_extra_quartics'] = 51
    write(f, d)
    with pytest.raises(ValueError, match='registry/version/coverage'):
        a.build_audit(repo)


def test_reject_missing_fitted_source(repo):
    f = (repo / a.STUDIES['small_uv_frozen_5e-9']['root'] /
         'fitted/T3_dS1_2_dS2_2_dF_1_alpha_m1.json')
    f.unlink()
    with pytest.raises(ValueError, match='Missing saved source'):
        a.build_audit(repo)


def test_output_summary_makes_boundary_qualification_explicit(repo):
    text = a.summary_md(a.build_audit(repo))
    assert 'not** an independently checked zero-valued coefficient' in text
    assert 'UV frozen' in text
    assert 'general_scalar_extended' in text
    assert 'No fits' in text
