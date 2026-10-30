"""No numerical-physics dependencies required for the original RGE canvas adapter.

Run in a repository checkout with scripts/T3Frozen64Sources.py installed.
"""
from __future__ import annotations
import ast
import importlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))


def reports():
    return importlib.import_module('Numerical.plotting.BuildT3FullScalarFrozen64HTML')


def test_uses_original_interactive_report_no_scalar_or_custom_canvas():
    report=reports()
    source=Path(report.__file__).read_text(encoding='utf-8')
    assert 'Numerical/plotting/InteractiveModelComparison.py' in source
    assert "node.name == 'write_dashboard'" in source
    assert "t.id == 'template'" in source
    assert 'BuildGeneralScalarLambdaReports' not in source
    assert 'DASHBOARD = ' not in source
    assert 'ctx.lineTo(' not in source
    assert 'not computed' in source


def test_original_rge_ast_template_is_used(tmp_path, monkeypatch):
    report=reports()
    module=tmp_path/'Numerical'/'plotting'/'InteractiveModelComparison.py'
    module.parent.mkdir(parents=True)
    template='''<html>__TITLE__<select id="quantity">__OPTIONS__</select>
<canvas id="plot"></canvas><button id="zoomFinal">final</button>
<select id="yScaleMode"></select><button id="resetZoom">Reset</button>
<script>const DATA=__DATA__; const META=__META__;
function defaultLimits(){return [s.mu_uv||1e7,s.mu_low||100];}</script>
<div class="layout"><aside></aside></div>
Numerical integration direction: UV → IR. Saved benchmark trajectories; no matching or RGE recalculation.</html>'''
    module.write_text('def write_dashboard():\n    template = '+repr(template)+'\n    return template\n',encoding='utf-8')
    monkeypatch.setattr(report,'ROOT',tmp_path)
    assert report.original_interactive_rge_template() == template


def sample_record(module):
    obs=dict(zip(module.policy.FIELDS,(7.5e-5,2.5e-3,.3088,.02248,.47)))
    pt=lambda mu:{'mu_gev':mu,'c5_norm':1e-15,'c5_components':{key:1e-16 for key in
                ('11','12','13','22','23','33')},
                'masses_ev_descending':[0.05,0.009,0.001], 'observables':obs}
    cell={'model':'A_+0','status':'NEUTRINO_COMPLETE'}
    neutrino={'status':'CandidateNoAdditionalQuadraticMatchedC5',
              'takagi_masses_ev_descending':[0.05,0.009,0.001],
              'oscillation_observables_normal_ordering':obs}
    running={'points':[pt(1000),pt(316.227766),pt(100)]}
    return cell,{'neutrino':neutrino},running


def test_real_saved_observables_are_adapted_to_legacy_schema():
    module=reports()
    record=sample_record(module)
    data,meta,options=module.rge_payload([record])
    assert data['quantities']['c5_norm']['A_+0']['x'] == [1000,316.227766,100]
    assert data['quantities']['dm3l']['A_+0']['y'] == [0.0025]*3
    assert data['quantities']['s13']['A_+0']['y'] == [0.02248]*3
    assert data['quantities']['m3']['A_+0']['y'] == [0.05]*3
    assert 'Scalar masses' not in options and 'lambda' not in options
    assert all(not v['weinberg'] for v in meta.values())
    assert data['matched']=={} and data['components']=={}


def test_original_html_receives_real_data_and_preserved_controls():
    module=reports()
    record=sample_record(module)
    template='''<header>__TITLE__
Numerical integration direction: UV → IR. Saved benchmark trajectories; no matching or RGE recalculation.
</header><select id="quantity">__OPTIONS__</select><canvas id="plot"></canvas>
<button id="zoomFinal"></button><select id="yScaleMode"></select>
<button id="resetZoom"></button><div class="layout"><aside></aside></div>
<script>const DATA=__DATA__;const META=__META__;function defaultLimits(){return [s.mu_uv||1e7,s.mu_low||100];}</script>'''
    page=module.render_original_interactive(template,'optimised_low',[record])
    assert '0.0025' in page
    assert '0.02248' in page
    assert 'id="zoomFinal"' in page and 'id="yScaleMode"' in page
    assert '[s.mu_scalar_threshold||1e3,s.mu_low||100]' in page
    assert 'not computed' in page
    assert 'mS1Sq' not in page
    assert 'ctx.lineTo(' not in page
    assert '__DATA__' not in page


def test_mass_outliers_are_not_classified_as_fits():
    module=reports()
    assert module.endpoint_screen(None)=='NOT_CALCULATED'
    assert module.endpoint_screen({'takagi_masses_ev_descending':[382,98,75],
        'oscillation_observables_normal_ordering':{'delta_m3l_sq_ev2':140605}})=='MASS_SCALE_OUTLIER'


def test_trajectory_generator_is_existing_solver():
    from scripts import GenerateT3Frozen64WeinbergRunning as generator
    text=Path(generator.__file__).read_text(encoding='utf-8')
    assert 'run_weinberg_trajectory(' in text
    assert 'prediction_from_mass_matrix(' in text
    assert 'c5_at_1tev_gev_inverse' in text
    assert 'c5_at_100gev_gev_inverse' in text
    assert 'C12 not solved' in text
