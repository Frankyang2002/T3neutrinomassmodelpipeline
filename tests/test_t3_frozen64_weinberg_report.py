"""Contract tests for the user's exact optimised_rge.html and narrow gamma scan.

Does not run the 64-case physics campaign; tests the display adapter and
strict gamma provenance. Intended project path: tests/test_t3_frozen64_weinberg_report.py.
"""
from __future__ import annotations

import csv
import importlib
import json
import re
from pathlib import Path

import pytest


r = importlib.import_module('Numerical.plotting.BuildT3FullScalarFrozen64HTML')


def sample():
    obs = dict(zip(r.policy.FIELDS, (7.5e-5, 2.5e-3, .3088, .02248, .47)))
    def point(mu):
        return {'mu_gev': mu, 'c5_norm': 1e-15,
                'c5_components': {key: 1e-16 for key in
                                  ('11','12','13','22','23','33')},
                'masses_ev_descending': [0.05,0.009,0.001],
                'observables': obs}
    cell = {'model': 'A_+0', 'benchmark':'optimised_low',
            'status':'NEUTRINO_COMPLETE', 'trajectory_sha256':'t',
            'neutrino_sha256':'n'}
    neutrino = {'status':'CandidateNoAdditionalQuadraticMatchedC5',
                'takagi_masses_ev_descending': [0.05,0.009,0.001],
                'oscillation_observables_normal_ordering':obs}
    return cell, {'neutrino':neutrino}, {'points':[point(1000),point(316.227766),point(100)]}


def test_template_is_the_actual_uploaded_report_without_historical_data():
    source = r.original_interactive_rge_template()
    assert '<title>__TITLE__</title>' in source
    assert source.count('__TITLE__') == 2
    assert 'const DATA=__DATA__;' in source
    assert 'const META=__META__;' in source
    assert 'T3 Optimised - Interactive RGE Running' not in source
    for name in ('quantity','all','none','models','zoomAll','zoomUV','zoomInt',
                 'zoomFinal','yScaleMode','plot','showComponents','caption','status'):
        assert f'id="{name}"' in source
    for fragment in ('canvas{width:100%;height:min(70vh,720px)',
                     'const COLORS=', 'function draw(){', 'canvas.addEventListener',
                     '<div class="layout">', '<div class="legend" id="regionsLegend"'):
        assert fragment in source
    assert len(source) < 25000  # historical 4 MB data is NOT secretly included
    assert 'BuildGeneralScalarLambdaReports' not in Path(r.__file__).read_text(encoding='utf-8')


def test_exact_template_css_js_survives_replacement():
    source = r.original_interactive_rge_template()
    page = r.render_original_interactive(source, 'optimised_low', [sample()])
    assert re.search(r'<style>(.*?)</style>', source, re.S).group(1) == \
           re.search(r'<style>(.*?)</style>', page, re.S).group(1)
    assert page.count('<canvas id="plot"') == 1
    assert page.count('function draw(){') == 1
    assert 'No model data are fabricated' not in page
    assert 'const DATA=__DATA__' not in page
    assert 'Scalar masses' not in page and 'mS1Sq' not in page
    assert '0.0025' in page and '0.02248' in page
    assert 'Intermediate C₁₂ running NOT computed' in page
    assert 'LLSS one-loop kernel η' in page
    assert '[s.mu_scalar_threshold||1e3,s.mu_low||100]' in page
    assert '<div style="margin:12px 20px' not in page  # no added banner/layout
    assert page.count('<optgroup') == 3


def test_real_schema_includes_final_weinberg_and_observables_only():
    data,meta,menu = r.rge_payload([sample()])
    assert data['quantities']['dm3l']['T3-A-p0']['y'] == [0.0025]*3
    assert data['quantities']['s13']['T3-A-p0']['y'] == [.02248]*3
    assert data['quantities']['m1']['T3-A-p0']['y'] == [.001]*3
    assert data['quantities']['m3']['T3-A-p0']['y'] == [.05]*3
    assert data['quantities']['c5_23']['T3-A-p0']['x'] == [1000,316.227766,100]
    assert data['models'][0]=='T3-A-m2'
    assert data['models'][-1]=='T3-E-p2'
    assert data['matched']=={} and data['components']=={}
    assert meta['llss_kernel']['block']=='analytic_kernel'
    assert all(not x['weinberg'] for x in meta.values())
    assert 'lambdaT3' not in menu and 'mS' not in menu


def make_scans(tmp_path, gamma):
    b, m, scan = 'optimised_low','A_+0','gamma_095_105_v1'
    root = tmp_path / scan
    target = root / b / 'A_p0.csv'
    target.parent.mkdir(parents=True,exist_ok=True)
    source_file = tmp_path / 'scripts' / 'GenerateT3FrozenFullScalarGamma.py'
    source_file.parent.mkdir(parents=True,exist_ok=True)
    source_file.write_text('# numerical gamma generator source\n',encoding='utf-8')
    data = {'model':m,'benchmark':b,'trajectory_sha256':'t',
            'neutrino_sha256':'n', 'generator_sha256':r.policy.digest(source_file),
            'gamma':list(gamma)}
    columns = ['model','benchmark','direction','gamma','status',*r.policy.FIELDS]
    with target.open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader()
        for d in r.policy.DIRECTIONS:
            for g in gamma:
                writer.writerow(dict(model=m,benchmark=b,direction=d,gamma=g,
                                     status='NeutrinoEvaluated',
                                     **dict(zip(r.policy.FIELDS, (7.5e-5,2.5e-3,.3088,.02248,.47)))))
    target.with_suffix('.json').write_text(json.dumps({
        'source':data,'csv_sha256':r.policy.digest(target)}),encoding='utf-8')
    return root


def test_narrow_gamma_actual_33_points_and_old_scan_rejected(tmp_path,monkeypatch):
    monkeypatch.setattr(r,'ROOT',tmp_path)
    cell,files,_ = sample()
    grid = r.NARROW_GAMMA
    assert grid == (.95,.96,.97,.98,.99,1.,1.01,1.02,1.03,1.04,1.05)
    old = make_scans(tmp_path, [round(.5+.1*i,1) for i in range(11)])
    with pytest.raises(ValueError,match='0.95'):
        r.scans_for(cell,files,tmp_path,'gamma_095_105_v1')
    narrow = make_scans(tmp_path,grid)
    result = r.scans_for(cell,files,tmp_path,'gamma_095_105_v1')
    assert result['centre']['gamma']==1.
    assert all(len(result['directions'][d])==11 for d in r.policy.DIRECTIONS)
    assert all(result['directions'][d][0]['gamma']==.95 for d in r.policy.DIRECTIONS)
    assert all(result['directions'][d][-1]['gamma']==1.05 for d in r.policy.DIRECTIONS)


def test_mass_outlier_is_not_reported_as_fit():
    assert r.endpoint_screen({'takagi_masses_ev_descending':[382.,98.,75.],
        'oscillation_observables_normal_ordering':{'delta_m3l_sq_ev2':140605.}}) == 'MASS_SCALE_OUTLIER'
