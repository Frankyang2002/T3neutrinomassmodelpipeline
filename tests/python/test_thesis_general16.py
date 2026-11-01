"""Self-contained checks for the new campaign coordinator (no physical solver)."""
from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / 'scripts' / 'RunT3GeneralAll16FromReference.py'
spec = importlib.util.spec_from_file_location('thesis_general16', SCRIPT)
assert spec and spec.loader
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_all_16_unique_and_keys():
    assert len(mod.MODELS) == 16
    assert len(set(mod.MODELS)) == 16
    assert len(mod.NAMES) == 18
    assert mod.ref_path('B_-1').name == 'T3_dS1_2_dS2_2_dF_1_alpha_m1.json'
    assert mod.ref_path('A_+0').name.endswith('alpha_p0.json')
    assert mod.tag('C_-1') == 'C_m1'


def test_matrix_round_trip_and_optimizer_binding():
    a = np.arange(9, dtype=float).reshape(3,3) / 19
    a[0,0] = -0.02
    encoded = mod.matrix_json(a)
    np.testing.assert_allclose(mod.matrix(encoded), a)
    values = {f'y1_{i}{j}':float(a[i,j]) for i in range(3) for j in range(3)}
    values.update({f'y2_{i}{j}':float(-a[i,j]) for i in range(3) for j in range(3)})
    config = {'base_state':{'ordinary':{'y1':encoded,'y2':encoded}},'benchmark_search':{'enabled':True},'sensitivity':{'enabled':True}}
    fitted = mod.config_with_yukawas(config, values)
    np.testing.assert_allclose(mod.matrix(fitted['base_state']['ordinary']['y2']), -a)
    assert config['benchmark_search']['enabled'] is True
    assert not fitted['benchmark_search']['enabled']
    assert not fitted['sensitivity']['enabled']


def test_invalid_args_rejected():
    good = dict(run_id='pilot_01',sobol_points=128,max_nfev=60,seed_count=3,
                target_chi2=1e-8,bound=3,sobol_floor=1e-5,sobol_fraction=.5,
                quartic_fraction=.1,mu2_fraction=.01,extra_yukawa_fraction=.05)
    mod.validate_args(SimpleNamespace(**good))
    for override in ({'sobol_points':127},{'run_id':'bad/path'}, {'mu2_fraction':-.2}):
        with pytest.raises(ValueError):
            mod.validate_args(SimpleNamespace(**dict(good,**override)))


def test_general_c5_requires_matched_full_general_provenance(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, 'ROOT', tmp_path)
    p = mod.general_coefficient_path('E_+2')
    assert p.name == 'c5_coefficient.txt'
    assert 'general_matchete_matching' in str(p)
    with pytest.raises(FileNotFoundError, match='full-GENERAL Matchete C5'):
        mod.ensure_general_matching('E_+2', generate=False, wolframscript='wolframscript')
    p.parent.mkdir(parents=True)
    p.write_text('Coupling[y1] Coupling[y2]')
    mod.save(p.parent/'general_matching_status.json', {
        'model':'E_+2', 'status':'Matched', 'model_selection_verified':True,
        'matching_scheme':'FullGeneralUnsplitMatcheteOneScaleCandidate'})
    assert mod.ensure_general_matching('E_+2', generate=False,
                                       wolframscript='wolframscript') == p


def test_summary_without_physics(tmp_path):
    mod.summarise(tmp_path, [
        {'model':'A_-2','status':'ExploratoryAcceptedFit','chi2_final':0.02,
         'prediction':{'sin2_theta13':0.022},'method':'test','output_folder':'a'},
        {'model':'C_-1','status':'Blocked','message':'missing beta'},
    ])
    import csv,json
    rows = list(csv.DictReader((tmp_path/'summary.csv').open()))
    assert len(rows)==2 and rows[0]['sin2_theta13']=='0.022'
    assert json.loads((tmp_path/'summary.json').read_text())['successful_fit_count']==1
