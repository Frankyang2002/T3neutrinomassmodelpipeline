from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
from types import ModuleType
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]

def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj

prep = load('pilot_prep', 'scripts/PrepareThesisNumericPilot.py')
html = load('pilot_html', 'scripts/BuildThesisNumericHTML.py')


def minimal():
    return {'kind':'t3_pipeline_numerical_results_v1',
            'base_state':{'ordinary':{'lambdaT3':{'real':.1, 'imag':0},
                                      'y1':{'real':[[.0001]*3 for _ in range(3)], 'imag':[[0]*3 for _ in range(3)]},
                                      'y2':{'real':[[.0002]*3 for _ in range(3)], 'imag':[[0]*3 for _ in range(3)]}}}}


def test_models_and_starting_values():
    assert prep.model_key('A_-2') == 'T3_dS1_1_dS2_3_dF_2_alpha_m2'
    assert prep.model_key('E_+2') == 'T3_dS1_3_dS2_3_dF_2_alpha_p2'
    b = minimal()
    out = prep.make_config(b, lam=.2, yukawa_factor=2.)
    assert out['base_state']['ordinary']['lambdaT3']['real'] == .2
    assert out['base_state']['ordinary']['y1']['real'][0][0] == .0002
    assert b['base_state']['ordinary']['y1']['real'][0][0] == .0001
    assert out['benchmark_search']['enabled'] is False
    assert out['sensitivity']['n_points'] == 11
    assert Path(out['sensitivity']['oscillation_target']).is_absolute()


def test_refit_config():
    o = prep.make_config(minimal(), refit=True, sobol_points=64)
    assert o['benchmark_search']['enabled'] is True
    assert len(o['benchmark_search']['parameters']) == 19
    assert o['benchmark_search']['local_refinement']['enabled'] is True


def test_observables_are_read_without_recalculation():
    d={'status':'Success','ordering':'NO','low_energy':{'masses_ev':[.01,.012,.05], 'takagi_residual':1e-12},
       'final_running':{k:[.1,.2] for k in html.OBS},
       'sensitivity':{'status':'Success','directions':{'T1':{'points':[{'gamma':1,'status':'Success','chi2':0.1}]}}}}
    r=html.observable_row(d,'Akey')
    assert r['m3_ev']==.05
    assert r['delta_m21_sq_ev2']==.2
    assert r['chi2_at_gamma1']==.1


def test_html_stages_existing_renderers(tmp_path):
    input_dir=tmp_path/'input'/'T3_dS1_1_dS2_3_dF_2_alpha_m2'/'data'
    input_dir.mkdir(parents=True)
    config=tmp_path/'config.json'
    config.write_text(json.dumps({'representation':{'alpha':-2,'d_s1':1,'d_s2':3,'d_f':2}}))
    d={'status':'Success','ordering':'NO','numerical_config':str(config),
       'low_energy':{'masses_ev':[0,.01,.05]},
       'final_running':{k:[.1] for k in html.OBS},
       'sensitivity':{'status':'Success','target':{'ordering':'NO'},
          'directions':{kind:{'points':[{'gamma':g,'status':'Success',
                 'prediction':{k:.1 for k in html.OBS}} for g in (.8,1.,1.2)]}
                        for kind in ('T1','T2','T3')}}}
    (input_dir/'running_diagnostics.json').write_text(json.dumps(d))
    plotting=ModuleType('Numerical.plotting')
    rge=ModuleType('Numerical.plotting.InteractiveModelComparison')
    ovals=ModuleType('Numerical.plotting.BuildFittedModelEnvelopes')
    def dash(source,target,title):
        Path(target).write_text('<html>original RGE</html>')
    def oval(source,target):
        assert (source/'comparison_summary.json').is_file()
        assert (source/'T3_dS1_1_dS2_3_dF_2_alpha_m2'/'sensitivity_diagnostics.json').is_file()
        Path(target).write_text('<html>original oval</html>')
    rge.write_dashboard=dash
    ovals.build=oval
    with patch.dict(sys.modules,{'Numerical':ModuleType('Numerical'),
                'Numerical.plotting':plotting,
                'Numerical.plotting.InteractiveModelComparison':rge,
                'Numerical.plotting.BuildFittedModelEnvelopes':ovals}):
        rows, path=html.write_outputs(tmp_path/'input', tmp_path/'out')
    assert len(rows)==1 and path is not None and path.is_file()
    assert (tmp_path/'out'/'interactive_RGE.html').is_file()
    assert (tmp_path/'out'/'observables.csv').is_file()
