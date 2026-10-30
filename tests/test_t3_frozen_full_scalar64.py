"""Pure contract tests. Full numerical RGEs must be tested in the user checkout."""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import pytest
from scripts import T3Frozen64Sources as p


def complete_manifest():
    return {'format':p.FORMAT,'run_id':'fixture','cells':[
        {'model':m,'benchmark':b,'yukawa_mode':mode,'case':case,'status':'NOT_RUN'}
        for b,(mode,case) in p.BENCHMARKS.items() for m in p.MODELS
    ]}


def test_64_exact_ids_and_stable_schema():
    assert len(p.MODELS)==16
    assert len(p.BENCHMARKS)==4
    assert len(p.expected_keys())==64
    assert len(set(p.MODELS))==16
    assert p.validate_manifest(complete_manifest())


def test_strict_source_separation():
    assert {case for _,case in p.BENCHMARKS.values()} == {
        'small_uv_frozen_5e-9','large_uv_frozen_0p5'}
    assert 'small_uv_running_5e-9' not in str(p.BENCHMARKS)
    assert p.BENCHMARKS['shared_a_m2_low'][0]=='shared_a_m2'
    assert p.tag('E_-4')=='E_m4'


def test_missing_duplicate_and_wrong_strategy_rejected():
    x=complete_manifest();x['cells'].pop()
    with pytest.raises(ValueError):p.validate_manifest(x)
    x=complete_manifest();x['cells'][1]=deepcopy(x['cells'][0])
    with pytest.raises(ValueError):p.validate_manifest(x)
    x=complete_manifest();x['cells'][0]['yukawa_mode']='shared_a_m2'
    with pytest.raises(ValueError):p.validate_manifest(x)


def test_source_paths_are_explicit_and_sha_verified(tmp_path):
    source=tmp_path/'output'/'source.json'
    source.parent.mkdir()
    source.write_text('{"ok":true}',encoding='utf-8')
    assert p.validate_source(tmp_path, 'output/source.json',p.digest(source))==source
    with pytest.raises(ValueError):p.validate_source(tmp_path,'output/source.json','0'*64)
    with pytest.raises(ValueError):p.safe_manifest_path(tmp_path,'../../private.json')
    with pytest.raises(ValueError):p.safe_manifest_path(tmp_path,'C:\\temp\\old.json')


def test_gamma_grid_contract():
    # Import here: the real matching and scalar numerics are intentionally not imported.
    import ast
    script=Path(__file__).resolve().parents[1]/'scripts'/'GenerateT3FrozenFullScalarGamma.py'
    tree=ast.parse(script.read_text(encoding='utf-8'))
    assert any(getattr(n,'name',None)=='perturb' for n in tree.body)
    assert any(getattr(n,'name',None)=='compute' for n in tree.body)


def test_html_report_is_presentation_only():
    script=Path(__file__).resolve().parents[1]/'Numerical/plotting/BuildT3FullScalarFrozen64HTML.py'
    source=script.read_text(encoding='utf-8')
    assert 'full_scalar_sequential_matching_certified' in source
    assert 'validate_source(' in source
    assert 'sensitivity_ovals.html' in source
    assert 'continuous_trajectory.json' in source
    assert 'bulk_48' not in source


def test_empty_64_slot_report_cannot_invent_physics(tmp_path,monkeypatch):
    # Standalone presentation contract without importing experimental physics.
    import sys
    import types
    import importlib
    fake=types.ModuleType('Numerical.plotting.BuildFittedModelEnvelopes')
    fake.HTML='<html><div class="layout"></div><script>const D=__DATA__;</script></html>'
    fake.LABELS={key:key for key in p.FIELDS}
    monkeypatch.setitem(sys.modules,'Numerical.plotting.BuildFittedModelEnvelopes',fake)
    render=importlib.import_module('Numerical.plotting.BuildT3FullScalarFrozen64HTML')
    monkeypatch.setattr(render,'ROOT',tmp_path)
    monkeypatch.setattr(render,'DATA',tmp_path/'runs')
    manifest=complete_manifest()
    (tmp_path/'runs'/'fixture').mkdir(parents=True)
    (tmp_path/'runs'/'fixture'/'manifest.json').write_text(json.dumps(manifest))
    (tmp_path/'docs').mkdir()
    (tmp_path/'docs'/'FROZEN_LAMBDA_TARGET_REFERENCE.json').write_text(json.dumps({
        'target':{'observable_names':p.FIELDS,'central_values':[1]*5,
                  'covariance':[[1 if i==j else 0 for j in range(5)] for i in range(5)]},
        'future_sigma':{},'future_source':'fixture'}))
    out=tmp_path/'reports'
    coverage=render.write_report('fixture','gamma_scans_v1',out)
    assert len(coverage)==64
    assert not any(x['scalar_rge'] or x['candidate_neutrino'] or x['complete_gamma_scan'] for x in coverage)
    assert (out/'index.html').is_file()
    assert len(list(out.glob('*/sensitivity_ovals.html')))==4
    assert len(list(out.glob('*/interactive_rge.html')))==4
    assert 'No verified saved scalar trajectories' in (out/'optimised_low'/'interactive_rge.html').read_text()
