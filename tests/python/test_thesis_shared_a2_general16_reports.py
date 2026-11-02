"""Offline tests for shared-A(-2) gamma transformations and report adapter."""
from __future__ import annotations
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest

from scripts import RunT3SharedA2General16Reports as m


def sample_params():
    return {k: float(i+1)/1e4 for i,k in enumerate(m.campaign.NAMES)}


def test_gamma_grid_center_and_log_spacing():
    g = m.gammas(0.7, 1.4, 7)
    assert len(g) == 7
    assert g[3] == 1.0
    assert np.all(np.diff(g)>0)
    assert g[0] == pytest.approx(.7)
    assert g[-1] == pytest.approx(1.4)


@pytest.mark.parametrize('limits', [(1, 2, 5), (.1, 1, 5), (0, 2, 5), (.7, 1.4, 4)])
def test_invalid_gamma_grid(limits):
    with pytest.raises(ValueError):
        m.gammas(*limits)


@pytest.mark.parametrize('direction',m.DIRECTIONS)
def test_gamma_center_identical(direction):
    params = sample_params()
    p, lam = m.perturb(params, direction, 1.0)
    assert p == params
    assert lam['real'] == 1.


@pytest.mark.parametrize('direction,expected',[
    ('T1', {'y1': 1/np.sqrt(4),'y2':1/np.sqrt(4),'lam':4}),
    ('T2', {'y1': 1,'y2':1/4,'lam':4}),
    ('T3', {'y1': 4,'y2':1/4,'lam':1}),
])
def test_gamma_directions(direction,expected):
    params=sample_params()
    p, lam =m.perturb(params,direction,4)
    assert lam['real']==pytest.approx(expected['lam'])
    for name in ('y1','y2'):
        for key in params:
            if key.startswith(name+'_'):
                assert p[key] == pytest.approx(params[key]*expected[name])
    assert params==sample_params()  # no mutation


def test_reference_lambda_copy_only():
    ref={'base_state': {'ordinary':{'lambdaT3':{'real':.1,'imag':0},'y1':{}}}}
    val=m.reference_with_lambda(ref,1.4)
    assert val['base_state']['ordinary']['lambdaT3']['real']==pytest.approx(.14)
    assert ref['base_state']['ordinary']['lambdaT3']['real']==.1


def test_fixed_conjugate_yukawas_from_a2(monkeypatch):
    shared = {'base_state':{'ordinary':{'y1': m.campaign.matrix_json(np.eye(3)),
                                      'y2':m.campaign.matrix_json(2*np.eye(3))}}}
    reference={'base_state':{'ordinary':{'lambdaT3':{'real':.1,'imag':0}}}}
    def stub_fixed(_reference,_uv,**_kw):
        return {'y1Conj':None,'y2Conj':None,'lambdaT3':{'real':.1,'imag':0}}
    monkeypatch.setattr(m.campaign,'fixed_overrides',stub_fixed)
    fixed=m.fixed_for_model('D_+0',reference,{},shared,
        {'quartic_fraction':.1,'mu2_fraction':.01,'extra_yukawa_fraction':.05},
        lambda_scale=1.,conj_source='a2-reference')
    from scripts.RunT3Aminus2YukawaBenchmark import complex_matrix
    assert np.allclose(complex_matrix(fixed['y1Conj']), .05*np.eye(3))
    assert np.allclose(complex_matrix(fixed['y2Conj']), .1*np.eye(3))


def test_run_one_writes_existing_renderer_schema(tmp_path, monkeypatch):
    params=sample_params()
    obs=m.OBS
    target=SimpleNamespace(ordering='NO',observable_names=obs,
        central_values=np.ones(len(obs)),covariance=np.eye(len(obs)),source='test')
    calls=[]
    def mock_evaluator(_model,*_args):
        def evaluate(point, lam_scale):
            calls.append((point,lam_scale))
            prediction=[1e-4,2e-3,.3,.02,.5]
            fit=SimpleNamespace(observable_names=obs,prediction=prediction,chi2=37.5)
            return fit,{'actual_endpoint':'saved'}
        return evaluate, {'matching_scheme':'test conditional', 'matching_source':'test'}
    monkeypatch.setattr(m,'load_evaluator',mock_evaluator)
    grid=m.gammas(.7,1.4,5)
    result=m.run_one('A_-2',params,{}, {},target,grid,'a2-reference',tmp_path/'A_m2')
    assert result['status']=='Success'
    assert len(calls)==15
    scan=json.loads((tmp_path/'A_m2'/'sensitivity_diagnostics.json').read_text())
    assert set(scan['directions'])=={'T1','T2','T3'}
    for d in m.DIRECTIONS:
        pts=scan['directions'][d]['points']
        assert len(pts)==5
        assert pts[2]['gamma']==1
        assert pts[2]['status']=='Success'
        assert set(pts[2]['prediction'])==set(obs)
    cen=json.loads((tmp_path/'A_m2'/'centre.json').read_text())
    assert cen['chi2']==37.5
    assert cen['numerical_details']=={'actual_endpoint':'saved'}


def test_cross_model_report_labels_not_individual_fits(tmp_path):
    m.report_html([{'model':'B_-1','status':'Blocked'}], tmp_path,
                  donor_run='thesis16_general_c5_full_v1', gamma_grid=[.7,1,1.4], oval=True)
    text=(tmp_path/'index.html').read_text()
    assert 'All models use the same fitted' in text
    assert 'B(-1) and C(-1) use conditional' in text
    assert 'sensitivity_oval.html' in text


def test_oval_adapter_preserves_model_names_and_tags(tmp_path, monkeypatch):
    import sys
    from types import ModuleType
    import numpy as np
    out = tmp_path/'campaign'
    out.mkdir()
    rows = [
        {'model':'A_-2','status':'Success','verified_chi2':34.5},
        {'model':'D_+0','status':'Success','verified_chi2':40.},
        {'model':'C_-3','status':'Blocked','message':'missing matching'},
    ]
    m.campaign.save(out/'comparison_summary.json', rows)
    for model in ('A_-2','D_+0'):
        folder=out/m.campaign.tag(model)
        folder.mkdir()
        m.campaign.save(folder/'sensitivity_diagnostics.json',
                        {'model':model,'status':'Success','directions':{'T1':{},'T2':{},'T3':{}}})

    calls=[]
    def pretend_renderer(source, html_path):
        calls.append((source,html_path))
        staged=m.saved_rows(source/'comparison_summary.json')
        assert staged == rows
        for model in ('A_-2','D_+0'):
            assert (source/model/'sensitivity_diagnostics.json').is_file()
        assert not (source/'C_-3'/'sensitivity_diagnostics.json').exists()
        html_path.write_text('Individually fitted T3 model sensitivity — Fitted centres — '
                             'Each model uses its own fitted λT3, Yukawa matrices and other fitted parameters.',
                             encoding='utf-8')
        return {'models':[{'model':'A_-2'},{'model':'D_+0'}]}

    # Replace only the presentation renderer, not the physics source or saver.
    module=ModuleType('Numerical.plotting.BuildFittedModelEnvelopes')
    module.build=pretend_renderer
    monkeypatch.setitem(sys.modules,'Numerical',ModuleType('Numerical'))
    monkeypatch.setitem(sys.modules,'Numerical.plotting',ModuleType('Numerical.plotting'))
    monkeypatch.setitem(sys.modules,'Numerical.plotting.BuildFittedModelEnvelopes',module)
    assert m.build_oval(out) is True
    assert len(calls)==1
    text=(out/'sensitivity_oval.html').read_text(encoding='utf-8')
    assert 'Shared A(-2) UV Yukawa T3 sensitivity' in text
    assert 'Shared-Yukawa centres' in text
    assert 'All models use the same A(-2) fitted UV Y1/Y2' in text


def test_render_saved_only_uses_no_numerical_calls(tmp_path,monkeypatch):
    run_id='test_shared'
    monkeypatch.setattr(m,'OUTPUT',tmp_path)
    out=tmp_path/run_id
    out.mkdir()
    m.campaign.save(out/'settings.json',{'format':m.FORMAT,'donor_run_id':'donor',
                                        'gammas':[0.7,1.,1.4]})
    rows=[{'model':'A_-2','status':'Success','verified_chi2':34.5}]
    m.campaign.save(out/'comparison_summary.json',rows)
    m.campaign.save(out/'manifest.json',{'format':m.FORMAT})
    monkeypatch.setattr(m,'build_oval',lambda folder: True)
    assert m.render_saved_only(run_id)==0
    assert 'sensitivity_oval.html' in (out/'index.html').read_text()
    payload=m.campaign.read(out/'manifest.json')
    assert payload['interactive_oval_generated'] is True
    assert payload['successful_scan_count']==1
