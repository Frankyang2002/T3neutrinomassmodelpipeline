"""Regression for the two distinct upstream physical-certification flag schemas."""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys

import pytest

from scripts import RunT3SelfBilinearNeutrino as stage


def _write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding='utf-8')


@pytest.fixture
def saved_outputs(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(stage, 'ROOT', tmp_path)
    source = tmp_path / 'output' / 'continuous_full_scalar' / 'continuous_v1'
    out = source / 'self_bilinear_v1'
    for model in stage.continuous.MODELS:
        tag = stage.continuous.model_tag(model)
        trajectory = source / tag / 'continuous_trajectory.json'
        trajectory.parent.mkdir(parents=True, exist_ok=True)
        trajectory.write_text('source '+model+'\n', encoding='utf-8')
        common = {
            'model': model,
            'source_continuous_sha256': sha256(trajectory.read_bytes()).hexdigest(),
            'takagi_masses_ev_descending': [0.05, 0.01, 0.002],
            'oscillation_observables_normal_ordering': {
                'delta_m21_sq_ev2': 7.5e-5,
                'delta_m3l_sq_ev2': 2.5e-3,
                'sin2_theta12': 0.30,
                'sin2_theta13': 0.022,
                'sin2_theta23': 0.55,
            },
        }
        if model in stage.SELF:
            report = {
                **common,
                'format': stage.FORMAT,
                'status': stage.STATUS,
                'full_general_scalar_c5_certified': False,
                'normal_channel_scalar_bilinear_included': True,
                'anomalous_prop_matching_included': False,
                'additional_independent_quartic_direct_matching_included': False,
            }
            path = out / f'{tag}_self_bilinear_c5.json'
        else:
            report = {
                **common,
                'format': stage.neutrino.FORMAT,
                'status': ('CandidateMassEigenBCEvolvedYukawasSMWeinberg'
                           if model in stage.neutrino.MIXED else
                           'CandidateNoAdditionalQuadraticMatchedC5'),
                'full_general_scalar_physical_c5_certified': False,
                'full_sequential_threshold_matching_performed': False,
                'additional_quadratic_c5_omitted': False,
            }
            path = source / 'neutrino_v1' / f'{tag}_neutrino.json'
        _write(path, report)
    return source, out


def test_existing_ten_model_flag_and_six_self_model_flag_are_merged(saved_outputs):
    source, out = saved_outputs
    report = stage.combined_16_summary(source_root=source, output_root=out,
                                       control_output_id='neutrino_v1')
    assert report['models_processed'] == 16
    assert report['models_self_bilinear_normal_channel_resummed'] == 6
    assert report['models_existing_candidate_c5'] == 10
    rows = {r['model']: r for r in report['rows']}
    assert rows['A_-4']['status'] == 'CandidateNoAdditionalQuadraticMatchedC5'
    assert rows['E_-2']['status'] == stage.STATUS
    assert report['physical_full_general_scalar_matching_complete'] is False
    assert (out / 'combined_16_summary.json').exists()
    assert len((out / 'combined_16_summary.csv').read_text().splitlines()) == 17


def test_missing_required_old_certification_flag_fails_closed(saved_outputs):
    source, out = saved_outputs
    model = 'A_-4'
    file = source / 'neutrino_v1' / f'{stage.continuous.model_tag(model)}_neutrino.json'
    payload = json.loads(file.read_text())
    del payload['full_general_scalar_physical_c5_certified']
    _write(file, payload)
    with pytest.raises(ValueError, match='A_-4: incompatible source'):
        stage.combined_16_summary(source_root=source, output_root=out,
                                  control_output_id='neutrino_v1')


def test_overclaim_in_self_bilinear_flag_fails_closed(saved_outputs):
    source, out = saved_outputs
    model = 'A_+0'
    file = out / f'{stage.continuous.model_tag(model)}_self_bilinear_c5.json'
    payload = json.loads(file.read_text())
    payload['anomalous_prop_matching_included'] = True
    _write(file, payload)
    with pytest.raises(ValueError, match='A_\\+0: incompatible source'):
        stage.combined_16_summary(source_root=source, output_root=out,
                                  control_output_id='neutrino_v1')


def test_combine_only_reuses_existing_results_without_recomputation(saved_outputs, monkeypatch):
    source, out = saved_outputs
    monkeypatch.setattr(sys, 'argv', [
        'RunT3SelfBilinearNeutrino', '--all', '--run-id', 'continuous_v1',
        '--source-root', str(source.parent), '--combine-only',
    ])
    def no_recompute(*args, **kwargs):
        raise AssertionError('must not rerun self-bilinear C5')
    monkeypatch.setattr(stage, 'calculate', no_recompute)
    assert stage.main() == 0
    assert (out / 'combined_16_summary.csv').exists()
