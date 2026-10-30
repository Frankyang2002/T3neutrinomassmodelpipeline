"""Isolated contracts; actual SciPy neutrino integrations run in local checkout."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from scripts import T3Frozen64Sources as policy
from scripts import T3ZeroExtraReportSources as zero
from Numerical.plotting import BuildT3ZeroExtraReferenceReports as report


OBS=dict(zip(policy.FIELDS, [7.537e-5,.002511,.3088,.02248,.47]))


def test_gamma_range_and_two_modes():
    assert zero.MODES == ('held_zero','boundary_zero')
    assert zero.NARROW_GAMMA == (.95,.96,.97,.98,.99,1.,1.01,1.02,1.03,1.04,1.05)
    assert len(zero.NARROW_GAMMA)*len(policy.DIRECTIONS)*len(policy.MODELS)*2 == 1056


def test_summary_requires_all_32_real_results(tmp_path,monkeypatch):
    monkeypatch.setattr(zero,'BASE',tmp_path)
    root=tmp_path/'some_run';root.mkdir()
    settings={'benchmark':'optimised_low','models':list(policy.MODELS),
              'modes':list(zero.MODES)}
    (root/'control_settings.json').write_text(json.dumps(settings))
    rows=[{'benchmark':'optimised_low','mode':mode,'model':m}
          for m in policy.MODELS for mode in zero.MODES]
    (root/'summary.json').write_text(json.dumps({'format':zero.FORMAT,
                                                'settings':settings,'rows':rows}))
    x=zero.load_summary('some_run')
    assert len(x['rows'])==32
    assert x['rows']['held_zero','A_-2']['model']=='A_-2'
    (root/'summary.json').write_text(json.dumps({'format':zero.FORMAT,
                                                'settings':settings,'rows':rows[:-1]}))
    with pytest.raises(ValueError,match='32'):
        zero.load_summary('some_run')


def test_zero_boundary_and_central_hashes_are_checked(tmp_path,monkeypatch):
    monkeypatch.setattr(zero,'ROOT',tmp_path)
    monkeypatch.setattr(zero,'BASE',tmp_path/'studies')
    mode,model='held_zero','A_+0'
    base=zero.BASE/'test_run'/'optimised_low'/mode
    folder=base/'A_p0';folder.mkdir(parents=True)
    tp=folder/'continuous_trajectory.json'
    nroot=base/'neutrino_v1';nroot.mkdir()
    np=nroot/'A_p0_neutrino.json'
    t={'model':model,'boundary_mode':'FROZEN_UV','zero_extra_control_mode':mode,
       'settings':{'benchmark':'optimised_low','zero_extra_control_mode':mode,
                   'frozen_uv':True},'uv_ode_evaluations':0,
       'scalar_endpoint_physically_admissible':True,
       'extra_quartic_initial_magnitudes':{'new':0.,'reused':0.},
       'quadratic_mu2_at_fermion_threshold_gev2':{'real':0.,'imag':0.}}
    tp.write_text(json.dumps(t));np.write_text(json.dumps({'model':model,
        'source_continuous_sha256':policy.digest(tp),
        'full_general_scalar_physical_c5_certified':False,
        'sm_weinberg_rge_performed':True}))
    meta={'benchmark':'optimised_low','rows':{(mode,model):{
        'lambdaT3_input_abs':5e-9,'max_extra_quartic_endpoint':0,
        'abs_mu2_endpoint_gev2':0}}}
    assert zero.sources('test_run',meta,mode,model)['cell']['neutrino_sha256']==policy.digest(np)
    t['extra_quartic_initial_magnitudes']['reused']=.01
    tp.write_text(json.dumps(t))
    # Provenance check correctly takes precedence after an altered source.
    with pytest.raises(ValueError,match='provenance incompatible'):
        zero.sources('test_run',meta,mode,model)
    n=json.loads(np.read_text())
    n['source_continuous_sha256']=policy.digest(tp)
    np.write_text(json.dumps(n))
    with pytest.raises(ValueError,match='initial value not zero'):
        zero.sources('test_run',meta,mode,model)


def test_original_reference_rge_is_used_without_scalar_plots():
    old=report.old
    template=old.original_interactive_rge_template()
    assert 'function draw(){' in template
    assert 'id="zoomFinal"' in template
    points=[]
    for mu in (1000.,316.227766,100.):
        points.append({'mu_gev':mu,'c5_norm':1e-15,
            'c5_components':{k:1e-16 for k in ('11','12','13','22','23','33')},
            'masses_ev_descending':[.05,.009,.001], 'observables':OBS})
    cell={'model':'A_+0','benchmark':'optimised_low'}
    n={'status':'CandidateNoAdditionalQuadraticMatchedC5',
       'takagi_masses_ev_descending':[.05,.009,.001],
       'oscillation_observables_normal_ordering':OBS}
    html=old.operator_html('held_zero optimised_low',[(cell,{'neutrino':n},{'points':points})])
    assert '<canvas id="plot"' in html
    assert 'function draw(){' in html
    assert 'T3-A-p0' in html
    assert 'Saved SM+Weinberg running' in html
    assert 'mS1Sq' not in html and 'lambdaH1Adj' not in html
    assert 'Original' not in html or 'llss_kernel' in html


def test_narrow_scan_validation_refuses_other_modes_and_grid(tmp_path,monkeypatch):
    monkeypatch.setattr(zero,'ROOT',tmp_path)
    monkeypatch.setattr(zero,'BASE',tmp_path/'studies')
    # Fake the source script files purely for digest validation.
    scripts=tmp_path/'scripts';scripts.mkdir()
    (scripts/'GenerateT3FrozenFullScalarGamma.py').write_text('# original gamma\n')
    monkeypatch.setattr(report.gamma_runner,'__file__',str(scripts/'GenerateT3FrozenFullScalarGamma.py'))
    mode,model='held_zero','A_+0'
    run_id='test_run';scan_id='gamma_095_105_v1'
    path=zero.BASE/run_id/scan_id/mode/'A_p0.csv';path.parent.mkdir(parents=True)
    src={'cell':{'model':model,'benchmark':'optimised_low',
                 'trajectory_sha256':'th','neutrino_sha256':'nh'},
         'neutrino':{'oscillation_observables_normal_ordering':OBS,'status':'Candidate'}}
    summary={'benchmark':'optimised_low','summary_sha256':'ssh','settings_sha256':'csh'}
    identity={'model':model,'mode':mode,'benchmark':'optimised_low',
        'gamma':list(zero.NARROW_GAMMA),'trajectory_sha256':'th',
        'neutrino_sha256':'nh','source_summary_sha256':'ssh',
        'source_settings_sha256':'csh',
        'generator_sha256':policy.digest(Path(report.gamma_runner.__file__)),
        'frozen_gamma_calculator_sha256':policy.digest(
            scripts/'GenerateT3FrozenFullScalarGamma.py')}
    fields=['model','benchmark','direction','gamma','status',*policy.FIELDS]
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for direction in policy.DIRECTIONS:
            for g in zero.NARROW_GAMMA:
                writer.writerow({'model':model,'benchmark':'optimised_low',
                    'direction':direction,'gamma':g,'status':'NeutrinoEvaluated',**OBS})
    meta=path.with_suffix('.json')
    meta.write_text(json.dumps({'source':identity,'csv_sha256':policy.digest(path),
                                'rows':33,'failures':0}))
    assert len(report.validated_scan(run_id,mode,src,scan_id,summary)['directions']['T1'])==11
    identity['mode']='boundary_zero'
    meta.write_text(json.dumps({'source':identity,'csv_sha256':policy.digest(path),
                                'rows':33,'failures':0}))
    with pytest.raises(ValueError,match='hash/grid mismatch'):
        report.validated_scan(run_id,mode,src,scan_id,summary)


def test_boundary_only_gamma_passes_model_specific_matching_directory(tmp_path, monkeypatch):
    """B/C previous failure: endpoint requires the individual model folder."""
    import sys
    import types
    from scripts import GenerateT3ZeroExtraGamma as gamma

    original = tmp_path / 'GenerateT3FrozenFullScalarGamma.py'
    original.write_text('# simulated existing physical gamma calculator\n')
    fake = types.ModuleType('scripts.GenerateT3FrozenFullScalarGamma')
    fake.__file__ = str(original)
    calls = []

    def compute(cell, grid, folder, path):
        calls.append((cell['model'], folder, tuple(grid)))
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('gamma,status\n1,NeutrinoEvaluated\n')
        return 33, 0

    fake.compute = compute
    monkeypatch.setitem(sys.modules, 'scripts.GenerateT3FrozenFullScalarGamma', fake)
    monkeypatch.setattr(gamma.policy, 'MODELS', ('B_-1',))
    monkeypatch.setattr(gamma.zero, 'run_root', lambda _: tmp_path)
    monkeypatch.setattr(gamma.zero, 'load_summary', lambda _: {
        'benchmark': 'optimised_low', 'summary_sha256': 's',
        'settings_sha256': 'c'})
    mode_root = tmp_path / 'optimised_low' / 'boundary_zero'
    monkeypatch.setattr(gamma.zero, 'sources', lambda run_id, meta, mode, model: {
        'cell': {'model':model, 'trajectory_sha256':'t','neutrino_sha256':'n'},
        'source_root': mode_root})

    result = gamma.run('zero_extra_all16_low_v1', 'gamma_095_105_v1', mode='boundary_zero')
    assert result['expected'] == result['complete'] == 1
    assert result['mode_selection'] == ['boundary_zero']
    assert calls == [('B_-1', mode_root / 'B_m1', tuple(zero.NARROW_GAMMA))]
    assert not (tmp_path / 'gamma_095_105_v1' / 'held_zero').exists()


def test_report_build_ignores_partial_held_zero_scans(tmp_path, monkeypatch):
    """Boundary-only HTML generation must never inspect the failed held_zero CSV."""
    from Numerical.plotting import BuildT3ZeroExtraReferenceReports as r
    import types

    monkeypatch.setattr(r.policy, 'MODELS', ('B_-1',))
    monkeypatch.setattr(r.zero, 'load_summary', lambda _: {
        'benchmark': 'optimised_low', 'summary_sha256': 's',
        'settings_sha256': 'c'})
    monkeypatch.setattr(r.zero, 'sources', lambda run_id, summary, mode, model: {
        'cell': {'model':model, 'benchmark':'optimised_low'},
        'neutrino': {'status': 'Candidate',
                     'takagi_masses_ev_descending':[.05, .009, .001],
                     'oscillation_observables_normal_ordering':OBS}})
    monkeypatch.setattr(r, 'validated_weinberg', lambda *args: {'points':[]})
    selected=[]
    def validate(*args):
        selected.append(args[1]);
        assert args[1] == 'boundary_zero'
        return {'centre': OBS}
    monkeypatch.setattr(r, 'validated_scan', validate)
    monkeypatch.setattr(r, 'sensitivity_html', lambda *args: ('<html>sensitivity</html>', 1))
    monkeypatch.setattr(r.old, 'operator_html', lambda *args: '<html>reference canvas</html>')
    monkeypatch.setattr(r.old, 'table_rows', lambda records: ('<tr></tr>', '<tr></tr>', ['model'], [['B_-1']]))
    monkeypatch.setattr(r.old, 'make_page', lambda title, body: '<html>'+body+'</html>')
    result = r.build(dest=tmp_path / 'reports', mode='boundary_zero')
    assert result['mode_selection'] == ['boundary_zero']
    assert set(result['counts']) == {'boundary_zero'}
    assert selected == ['boundary_zero']
    assert (tmp_path/'reports'/'boundary_zero'/'interactive_rge.html').exists()
    assert not (tmp_path/'reports'/'held_zero').exists()
