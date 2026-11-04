"""Study-local contract tests; no production physics modules modified.

The driver tests use a deliberately synthetic evaluator; they do NOT validate
physical beta functions or Matchete expressions, which must run in the repo.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from scripts import RunT3GeneralFermionThresholdScan as study


def base_config(model="A_-4"):
    family, alpha = model.split("_")
    d1,d2,df=study.REPS[family]
    def encode(a):
        a = np.asarray(a, complex)
        return {"real": a.real.tolist(), "imag": a.imag.tolist()}
    return {
        "representation": {"d_s1":d1,"d_s2":d2,"d_f":df,
                           "alpha":int(alpha),"shared_scalar":False},
        "ordering":"NO",
        "base_state": {"mu_gev":1e7, "ordinary": {
            "MF":encode(np.diag([1e5,1.1e5,1.2e5])),
            "y1":encode(np.full((3,3),1e-4)),
            "y2":encode(np.full((3,3),2e-4)),
            "lambdaT3":{"real":0.1,"imag":0.0}}},
        "scales":{"mu_fermion_threshold_gev":1e5,
                  "mu_scalar_threshold_gev":1e3,"mu_low_gev":1e2},
        "numerical":{"vev_gev":246.22},
    }


def reference_prediction():
    return {"delta_m21_sq_ev2":7.537e-5,"delta_m3l_sq_ev2":2.511e-3,
            "sin2_theta12":0.3088,"sin2_theta13":0.02248,
            "sin2_theta23":0.47}


def test_rescale_changes_only_physical_mf_and_muf():
    cfg = base_config()
    old_json = json.dumps(cfg, sort_keys=True)
    for factor in study.FACTORS:
        out = study.scale_configuration(cfg, factor)
        assert study.scales(out) == (1e7, 1e5*factor, 1e3, 1e2)
        assert np.allclose(study.matrix(out['base_state']['ordinary']['MF']),
                           factor * study.matrix(cfg['base_state']['ordinary']['MF']))
        for name in ('y1','y2'):
            assert out['base_state']['ordinary'][name] == cfg['base_state']['ordinary'][name]
        restored = study.scale_configuration(out, 1.0)
        # Scaling by 1 retains this point; the original is never overwritten.
        assert study.matrix(restored['base_state']['ordinary']['MF'])[0,0] == 1e5*factor
    assert json.dumps(cfg, sort_keys=True) == old_json


def test_scale_factor_outside_small_scan_is_rejected():
    with pytest.raises(ValueError, match="0.5,1,2"):
        study.scale_configuration(base_config(), 10.)


def test_invalid_hierarchy_rejected():
    cfg = base_config()
    cfg['scales']['mu_fermion_threshold_gev'] = 5e2
    with pytest.raises(ValueError, match="hierarchy"):
        study.scales(cfg)


def test_no_simultaneous_decoupling_by_hierarchy_guard():
    with pytest.raises(ValueError, match="hierarchy"):
        study.ensure_mass_hierarchy(5e4,1e3,np.array([500.,800.,900.]),
                                    np.array([1000.,1200.]))


def test_oscillation_baseline_gate():
    obs = reference_prediction()
    assert study.baseline_agreement(obs, obs, rtol=0.03)['pass']
    shifted = dict(obs, sin2_theta23=0.6)
    assert not study.baseline_agreement(shifted, obs, rtol=0.03)['pass']


def test_relative_direction_masses_and_mixing():
    a = {'masses_ev':[0.001,0.01,0.05], 'sum_masses_ev':0.061,
         'observables':reference_prediction()}
    b = {'masses_ev':[0.002,0.005,0.025], 'sum_masses_ev':0.032,
         'observables':dict(reference_prediction(),sin2_theta23=0.48)}
    t = study.trend(a,b)
    assert t['m1_ev']['direction']=='increase'
    assert t['m3_ev']['direction']=='decrease'
    assert t['m3_ev']['relative_change']==pytest.approx(-0.5)
    assert t['sin2_theta23']['relative_change']==pytest.approx(0.48/0.47-1)


def test_reject_unsafe_manifest_path(tmp_path, monkeypatch):
    monkeypatch.setattr(study,'ROOT',tmp_path)
    for name in ('../../x', 'output\\..\\x', 'C:\\not_local\\file.json','/abs/path'):
        with pytest.raises(ValueError):
            study.repo_path(name)
    assert study.repo_path('output\\valid\\result.json') == tmp_path/'output/valid/result.json'


def test_manifest_must_contain_exact_sixteen():
    rows = [{'model':model, 'status':'ExploratoryAcceptedFit',
             'prediction':reference_prediction(), 'result_path':f'output/{study.tag(model)}/result.json'}
            for model in study.MODELS]
    saved = {'accepted_count':16,'blocked_count':0,'rows':rows}
    assert len(study.select_rows(saved,study.PILOT))==16
    with pytest.raises(ValueError,match='exactly 16'):
        study.select_rows(dict(saved,rows=rows[:-1]),study.PILOT)


def write_fake_study_root(tmp_path, monkeypatch, model='A_-4'):
    monkeypatch.setattr(study,'ROOT',tmp_path)
    monkeypatch.setattr(study,'OUTPUT',tmp_path/'output/studies/thesis_qualitative/fermion_threshold')
    monkeypatch.setattr(study,'MANIFEST',tmp_path/'saved_manifest.json')
    rows=[]
    for m in study.MODELS:
        cfg = base_config(m)
        model_dir=tmp_path/'output'/'fits'/study.tag(m)
        model_dir.mkdir(parents=True)
        chosen = {f'{name}_{i}{j}':float(study.matrix(cfg['base_state']['ordinary'][name])[i,j].real)
                  for name in ('y1','y2') for i in range(3) for j in range(3)}
        result={'model':m,'chi2_final':1e-5,'prediction':reference_prediction(),
                'best_parameters':chosen,
                'full_general_sequential_matching_performed':False}
        (model_dir/'result.json').write_text(json.dumps(result),encoding='utf8')
        (model_dir/'fitted_config.json').write_text(json.dumps(cfg),encoding='utf8')
        if m not in study.MIXED:
            (model_dir/'general_uv_overrides.json').write_text('{}',encoding='utf8')
        rows.append({'model':m,'status':'ExploratoryAcceptedFit',
                     'chi2_final':1e-5,'prediction':reference_prediction(),
                     'result_path':f'output/fits/{study.tag(m)}/result.json'})
    study.MANIFEST.write_text(json.dumps({'accepted_count':16,'blocked_count':0,'rows':rows}),
                              encoding='utf8')
    monkeypatch.setattr(study,'checkout_head', lambda:None)
    monkeypatch.setattr(study,'preflight', lambda *a,**kw:{'model':model,'status':'Ready',
                                                         'inputs':{},'problems':[]})


def test_study_uses_baseline_first_and_varies_physical_threshold(tmp_path, monkeypatch):
    write_fake_study_root(tmp_path,monkeypatch)
    calls=[]
    def synthetic(source, cfg, **kwargs):
        mu_f=study.scales(cfg)[1]
        r=mu_f/1e5
        calls.append((r,float(study.matrix(cfg['base_state']['ordinary']['MF'])[0,0].real)))
        pred={key:value*(1+0.05*(r-1)) for key,value in reference_prediction().items()}
        return {'masses_ev':[0.001/r,0.01/r,0.05/r],
                'sum_masses_ev':0.061/r,'observables':pred,
                'chi2_against_fixed_target':1e-5,
                'uv_endpoint_gev':mu_f,'eft1_initial_gev':mu_f,'eft1_endpoint_gev':1e3}
    monkeypatch.setattr(study,'_nonmixed',synthetic)
    assert study.run(('A_-4',),run_id='test_run',preflight_only=False,
                     baseline_rtol=.03,canonical_run_id='full16_general_v1',
                     wolframscript='unused',rtol=1e-7,atol=1e-10,
                     allow_different_head=False)==0
    assert calls==[(1.,1e5),(.5,5e4),(2.,2e5)]
    folder=study.OUTPUT/'test_run'
    rows=study.read(folder/'summary.json')['rows']
    assert [row['status'] for row in rows]==['ComputedConditional']*3
    assert rows[1]['relative']['m3_ev']['relative_change']==pytest.approx(1.)
    assert rows[2]['relative']['m3_ev']['relative_change']==pytest.approx(-0.5)
    assert (folder/'A_m4'/'factor_0p5.json').exists()
    assert (folder/'A_m4'/'factor_2p0.json').exists()


def test_baseline_disagreement_stops_mass_variations(tmp_path, monkeypatch):
    write_fake_study_root(tmp_path,monkeypatch)
    calls=[]
    def bad(source,cfg,**kwargs):
        calls.append(study.scales(cfg)[1])
        return {'masses_ev':[0.001,0.01,0.05], 'sum_masses_ev':0.061,
                'observables':dict(reference_prediction(),sin2_theta23=0.6),
                'uv_endpoint_gev':1e5,'eft1_initial_gev':1e5,'eft1_endpoint_gev':1e3}
    monkeypatch.setattr(study,'_nonmixed',bad)
    assert study.run(('A_-4',),run_id='bad_ref',preflight_only=False,
                     baseline_rtol=.03,canonical_run_id='full16_general_v1',
                     wolframscript='unused',rtol=1e-7,atol=1e-10,
                     allow_different_head=False)==2
    assert calls==[1e5]
    states=[r['status'] for r in study.read(study.OUTPUT/'bad_ref'/'summary.json')['rows']]
    assert states==['BaselineMismatch','Blocked']
    assert not (study.OUTPUT/'bad_ref'/'A_m4'/'factor_0p5.json').exists()


def test_no_reuse_of_historical_fitted_scale_or_restricted_c5():
    import inspect
    source = inspect.getsource(study)
    assert 'RunFittedFermionScaleScan' not in source
    body = inspect.getsource(study._mixed)
    assert 'bc.CandidateEvaluator(' not in body
    assert 'bc.integrate_point(' not in body
    assert 'mu_uv, mu_f, n_points=3' in body
    assert 'mu_f, mu_s, n_points=3' in body
    assert 'bc.c5_at_endpoint(model, snap, heavy, compiled)' in body
    nonmixed = inspect.getsource(study._nonmixed)
    assert 'GeneralT3FlavourC5' in nonmixed
    assert 'ensure_general_matching(model, generate=False' in nonmixed
