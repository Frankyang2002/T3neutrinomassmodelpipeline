"""Offline regression checks for the read-only UV scalar hierarchy diagnostic.

No fake fixture here represents a successful live RGBeta integration.
"""
from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest

from scripts import DiagnoseT3BCUVScalarHierarchy as mod


def row(mu, light, heavy, fmin=8.e4, fmax=1.2e5):
    return {'mu_gev':mu,'m_scalar_light_gev':light,'m_scalar_heavy_gev':heavy,
            'min_fermion_mass_gev':fmin,'max_fermion_mass_gev':fmax}


def test_scalar_eigenmasses_zero_mixing():
    assert mod.scalar_masses(1.e6,4.e6,0j)==(1000.,2000.)


def test_scalar_eigenmasses_complex_mixing():
    # Mathematically equivalent to [[2e6,3e5*i],[-3e5*i,2e6]].
    low,high=mod.scalar_masses(2.e6,2.e6,3e5j)
    assert np.isclose(low**2,1.7e6)
    assert np.isclose(high**2,2.3e6)


@pytest.mark.parametrize('m1,m2,mu',[(0,0,0j),(-1,4,0j),(1,1,2j),
                                  (np.nan,5,0j),(1,1,complex(np.inf,0))])
def test_mass_nonpositive_fails_closed(m1,m2,mu):
    with pytest.raises(ValueError):mod.scalar_masses(m1,m2,mu)


def test_crossing_bracket_lognormal():
    samples=[row(1e6,1e3,1e4),row(1e5,1e4,2e5)]
    z=mod.crossing(samples)
    assert z is not None
    assert 1e5<z['mu_geV_interpolated']<1e6
    assert z['bracket_mu_gev']==[1e6,1e5]
    assert z['numerical_crossing_only']


def test_no_crossing():
    assert mod.crossing([row(1e6,1e3,1e4),row(1e5,1e3,1e4)]) is None


def test_multiple_crossings():
    data=[row(1e6,1e3,1e4),row(1e5,1e4,2e5),row(1e4,1e3,2e3)]
    assert len(mod.crossings(data,threshold='mu_gev',mass='m_scalar_heavy_gev'))==2


def test_crossing_heavy_fermion_reference():
    z=mod.crossing([row(1e7,1e3,1e4),row(1e5,9e4,1.7e5)],
                   threshold='min_fermion_mass_gev',mass='m_scalar_heavy_gev')
    assert z is not None and 1e5<z['mu_geV_interpolated']<1e7


def test_path_guards(tmp_path):
    (tmp_path/'okay.json').write_text('{}',encoding='utf-8')
    assert mod.safe_input(tmp_path,'okay.json').is_file()
    for bad in ('../secret','/abs/file','missing.json'):
        with pytest.raises(ValueError):mod.safe_input(tmp_path,bad)


def test_exact_input_fingerprint(tmp_path):
    p=tmp_path/'original.json';p.write_text('{}',encoding='utf-8')
    h=sha256(p.read_bytes()).hexdigest()
    assert mod.verify_input_hashes(tmp_path,{'prior_input_file_sha256':{'original.json':h}})==1
    p.write_text('{"edited":true}')
    with pytest.raises(ValueError,match='SHA-256 changed'):
        mod.verify_input_hashes(tmp_path,{'prior_input_file_sha256':{'original.json':h}})


def test_empty_input_evidence_fails(tmp_path):
    with pytest.raises(ValueError,match='fingerprints'):
        mod.verify_input_hashes(tmp_path,{})


def test_math_diagnostics_are_separate_not_pole_claim():
    first=row(1e5,3e4,1.4e5)
    last=row(1e3,8e4,1.5e5)
    for r in (first,last):
        r.update({'mS1Sq_gev2':r['m_scalar_light_gev']**2,
                  'mS2Sq_gev2':r['m_scalar_heavy_gev']**2})
    restricted=[dict(first),dict(last)]
    restricted[0]['mS1Sq_gev2']=first['mS1Sq_gev2']-1e3
    restricted[1]['mS2Sq_gev2']=last['mS2Sq_gev2']-5e4
    g={'UV':[first],'EFT1':[last]}
    r={'UV':[restricted[0]],'EFT1':[restricted[1]]}
    info=mod.diagnostics(g,r,{'chi2_final':17000.},'B_-1')
    assert not info['fermion_first_ordering_at_muF']
    assert not info['scalar_decoupling_scale_close_to_muS']
    assert not info['passes_chi2_5p86']
    assert 'Running-parameter' in info['conclusion']


def test_summary_reports_scope():
    sample={'m_scalar_heavy_gev':1.1e5,'m_scalar_light_gev':9e4,
            'min_fermion_mass_gev':8e4}
    diag={'fitted_chi2':123.,'general_at_muF':sample,'general_at_muS':sample,
          'delta_mass_sq_at_muS_over_general':{'mS1Sq_gev2':0.1,'mS2Sq_gev2':0.2},
          'fermion_first_ordering_at_muF':False}
    s=mod.make_summary({'models':{'B_-1':{'diagnostic':diag,
        'general':{'crossings_heavy_scalar_equals_mu':[]}}}})
    assert 'Running mass' in s or 'running mass' in s
    assert 'NO' in s
    assert 'not physical pole masses' in s.lower()


def test_run_id_collision_fails_before_loading_files(tmp_path):
    d=tmp_path/'output/studies/general_scalar_extended/diagnostics/existing'
    d.mkdir(parents=True)
    (d/'marker.txt').write_text('keep',encoding='utf-8')
    with pytest.raises(ValueError,match='output already exists'):
        mod.execute(models=('B_-1',),run_id='existing',root=tmp_path)
    assert (d/'marker.txt').read_text()=='keep'


def test_reject_duplicate_model_before_files(tmp_path):
    with pytest.raises(ValueError,match='models must be unique'):
        mod.execute(models=('B_-1','B_-1'),run_id='fresh',root=tmp_path)


def test_reject_unsupported_model_before_files(tmp_path):
    with pytest.raises(ValueError,match='models must be unique'):
        mod.execute(models=('A_+0',),run_id='fresh',root=tmp_path)


def test_reject_bad_grid_before_files(tmp_path):
    with pytest.raises(ValueError,match='grid must'):
        mod.execute(models=('B_-1',),run_id='fresh',grid=2,root=tmp_path)


def test_signed_boundary_scenarios_not_relabelled():
    assert mod.CASE=='small_uv_running_5e-9'
    assert mod.DEFAULT_RUNS['B_-1'] != mod.DEFAULT_RUNS['C_-1']
    assert mod.FORMAT.endswith('v1')


def _saved_fit_fixture(tmp_path: Path, model: str = 'B_-1') -> tuple[dict, dict]:
    """Contract fixture for the Batch 11/12 fitted-config interface."""
    tag = mod.b9.safe_tag(model)
    folder = (tmp_path / 'output/studies/general_scalar_extended' /
              mod.CASE / mod.DEFAULT_RUNS[model])
    folder.mkdir(parents=True)
    vals = {key: (i + 1) / 100.0 for i, key in enumerate(mod.b11.NAMES)}
    def matrix(offset):
        return {'real': [[vals[mod.b11.NAMES[offset + 3*i + j]]
                          for j in range(3)] for i in range(3)],
                'imag': [[0.0]*3 for _ in range(3)]}
    cfg = {
        'representation': {'d_s1': 2, 'd_s2': 2,
                           'd_f': 1 if model == 'B_-1' else 3,
                           'alpha': -1, 'shared_scalar': False},
        'base_state': {'ordinary': {'y1': matrix(0), 'y2': matrix(9)}},
    }
    prior_inputs = {}
    for i in range(4):
        p = tmp_path / f'prerequisite_{i}.json'
        p.write_text(f'{{"prerequisite": {i}}}', encoding='utf-8')
        prior_inputs[p.name] = sha256(p.read_bytes()).hexdigest()
    fit = {
        'model': model, 'case': mod.CASE,
        'run_id': mod.DEFAULT_RUNS[model],
        'UV_running_enabled': True,
        'one_scale_matchete_expression_reused': True,
        'physical_matchete_rgbeta_map_certified': False,
        'best_parameters': vals,
        'chi2_final': 120.0,
        'source_manifest': {'prior_input_file_sha256': prior_inputs,
                            'quartic_fraction': 0.1, 'mu2_fraction': 0},
    }
    (folder / f'{tag}_exploratory_fit.json').write_text(
        json.dumps(fit), encoding='utf-8')
    (folder / f'{tag}_exploratory_fitted_config.json').write_text(
        json.dumps(cfg), encoding='utf-8')
    return cfg, fit


@pytest.mark.parametrize('model', mod.MODELS)
def test_saved_fit_loading_passes_wrapped_config_to_seed_vector(tmp_path, monkeypatch, model):
    cfg, fit = _saved_fit_fixture(tmp_path, model)
    historic = {'cfg': {'historic': True}, 'new': ('unused',),
                'reused': (), 'uv_on': True}
    monkeypatch.setattr(mod.b9, 'study_source', lambda *args, **kwargs: historic)
    seen = []
    def seed_from_saved_wrapper(source):
        # Fail on the previous bug: a bare config has no `cfg` key.
        assert source['cfg'] == cfg
        seen.append(True)
        ordinary = source['cfg']['base_state']['ordinary']
        return np.concatenate([
            np.asarray(ordinary[key]['real'], dtype=float).ravel()
            for key in ('y1', 'y2')])
    monkeypatch.setattr(mod.b11, 'seed_vector', seed_from_saved_wrapper)
    saved, src, settings = mod.fit_source(
        tmp_path, model, mod.DEFAULT_RUNS[model], {})
    assert seen == [True]
    assert src is historic
    assert saved == fit
    assert settings['verified_source_files'] == 4
    assert settings['quartic_fraction'] == 0.1


def test_saved_fit_loading_rejects_best_vector_mismatch(tmp_path, monkeypatch):
    cfg, _ = _saved_fit_fixture(tmp_path)
    monkeypatch.setattr(mod.b9, 'study_source', lambda *args, **kwargs: {'cfg': cfg})
    def valid_seed(source):
        return np.concatenate([
            np.asarray(source['cfg']['base_state']['ordinary'][key]['real'], dtype=float).ravel()
            for key in ('y1', 'y2')])
    monkeypatch.setattr(mod.b11, 'seed_vector', valid_seed)
    path = (tmp_path / 'output/studies/general_scalar_extended' / mod.CASE /
            mod.DEFAULT_RUNS['B_-1'] / 'B_m1_exploratory_fit.json')
    fit = json.loads(path.read_text(encoding='utf-8'))
    fit['best_parameters']['y1_00'] += 0.125
    path.write_text(json.dumps(fit), encoding='utf-8')
    with pytest.raises(ValueError, match='best parameter vector disagree'):
        mod.fit_source(tmp_path, 'B_-1', mod.DEFAULT_RUNS['B_-1'], {})
