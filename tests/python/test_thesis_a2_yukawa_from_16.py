"""Offline tests of A(-2) reuse from completed 16-model general campaign."""
import json
from hashlib import sha256

import numpy as np
import pytest

from scripts import RunT3Aminus2YukawaBenchmark as bench


def make_campaign(tmp_path, monkeypatch, *, multi=True):
    root = tmp_path / 'campaign'
    run = root / 'thesis16_general_c5_full_v1'
    a2 = run / 'A_m2'
    a2.mkdir(parents=True)
    kernel = tmp_path / 'c5_coefficient.txt'
    kernel.write_text('valid mock general coefficient', encoding='utf-8')
    monkeypatch.setattr(bench, 'CAMPAIGN_OUT', root)
    monkeypatch.setattr(bench.campaign, 'general_coefficient_path', lambda model: kernel)
    ref = {'representation': {'d_s1': 1, 'd_s2': 3, 'd_f': 2, 'alpha': -2}}
    params = {name: (i + 1)*1e-6 for i, name in enumerate(bench.campaign.NAMES)}
    fitted = {
        'representation': ref['representation'],
        'base_state': {'ordinary': {
            n: {'real': [[params[f'{n}_{i}{j}'] for j in range(3)] for i in range(3)],
                'imag': [[0.0]*3 for _ in range(3)]}
            for n in ('y1', 'y2')
        }}}
    settings = {'matching_scheme': 'full-general-Matchete-Yukawa-channel-massbasis-v2',
                'models': ['A_-4', 'A_-2', 'B_-1'] if multi else ['A_-2'],
                'mu2_fraction': 0.01, 'quartic_fraction': 0.1}
    result = {'model': 'A_-2', 'status': 'ExploratoryPoorFit',
              'chi2_final': 10.2, 'best_parameters': params,
              'best_source': 'local_refinement',
              'input_sha256': {'c5_kernel': sha256(kernel.read_bytes()).hexdigest()},
              'fixed_general_scalar_boundary': {
                  'bilinear_fraction_min_mass_squared': 0.01,
                  'quartic_fraction_of_lambdaT3': 0.1}}
    for file, data in [('settings.json', settings),
                       ('A_m2/fitted_config.json', fitted),
                       ('A_m2/result.json', result)]:
        (run/file).write_text(json.dumps(data), encoding='utf-8')
    return run, ref, params


@pytest.mark.parametrize('multi', [True, False])
def test_loads_fitted_yukawas_from_campaign(tmp_path, monkeypatch, multi):
    _, ref, params = make_campaign(tmp_path, monkeypatch, multi=multi)
    best, settings, source = bench.load_fit('thesis16_general_c5_full_v1', ref)
    assert best == params
    assert source['source_campaign_model_count'] == (3 if multi else 1)
    assert source['fit_chi2_original'] == 10.2
    assert source['fit_result_sha256']
    assert source['fit_config_sha256']


def test_rejects_mismatched_fitted_config(tmp_path, monkeypatch):
    run, ref, _ = make_campaign(tmp_path, monkeypatch)
    file = run / 'A_m2/fitted_config.json'
    data = json.loads(file.read_text())
    data['base_state']['ordinary']['y1']['real'][0][0] *= 2
    file.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='differs from optimizer'):
        bench.load_fit('thesis16_general_c5_full_v1', ref)


def test_rejects_missing_A2_result(tmp_path, monkeypatch):
    run, ref, _ = make_campaign(tmp_path, monkeypatch)
    (run / 'A_m2/result.json').unlink()
    with pytest.raises(ValueError, match='A\\(-2\\) must have completed'):
        bench.load_fit('thesis16_general_c5_full_v1', ref)


def test_rejects_wrong_model_campaign(tmp_path, monkeypatch):
    run, ref, _ = make_campaign(tmp_path, monkeypatch)
    file = run / 'settings.json'
    data = json.loads(file.read_text())
    data['models'] = ['A_-4', 'B_-1']
    file.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='does not include A'):
        bench.load_fit('thesis16_general_c5_full_v1', ref)


def test_rejects_old_non_mass_basis_matching_scheme(tmp_path, monkeypatch):
    run, ref, _ = make_campaign(tmp_path, monkeypatch)
    file = run / 'settings.json'
    settings = json.loads(file.read_text())
    settings['matching_scheme'] = 'full-general-Matchete-Yukawa-channel-v1'
    file.write_text(json.dumps(settings))
    with pytest.raises(ValueError, match='mass-basis-v2 campaign'):
        bench.load_fit('thesis16_general_c5_full_v1', ref)
