"""Tests for the reduced-scalar adapter to the original lambda-study renderer."""
from __future__ import annotations
import csv
import json
from pathlib import Path
import pytest
from scripts import BuildReducedScalarLambdaStyleReports as report


def write_csv(path,rows):
    with path.open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def synthetic(path):
    path.mkdir(parents=True)
    fits=[];pulls=[];running=[]
    shared_good=('A_-2','A_+0','B_+1','D_+0','D_-2')
    for case in report.CASES:
        for m in report.MODELS:
            ref=(.001 if m=='A_-2' else .002) if m in shared_good else 10.
            fits.append({'case':case,'model':m,'chi2_individual':.001,
                'chi2_shared_A_m2':ref,'individual_accepted':True,
                'shared_A_m2_accepted':ref<=report.CUTOFF,
                'lambdaT3_input':(5e-9 if case==report.CASES[0] else .5)*
                                  (.5 if m.startswith('E_') else 1),
                'c5_source_kind':'hierarchical'})
            for mode,score in (('individual',.001),('shared_A_-2',ref)):
                for obs in report.OBS:
                    pulls.append({'case':case,'model':m,'yukawa_choice':mode,
                        'observable':obs,'prediction':1.,'target':1.,'sigma':.2,
                        'pull':0.,'chi2_five_observables':score})
    for m in report.MODELS:
        running.append({'model':m,'screen_pass':False,
            'mS1_formal_muS_tev':40.,'mS2_formal_muS_tev':60.,'rhoF_at_muF':.42})
    write_csv(path/'frozen_comparison.csv',fits)
    write_csv(path/'oscillation_observable_pulls.csv',pulls)
    write_csv(path/'running_uv_formal_endpoints.csv',running)
    (path/'summary.md').write_text('Synthetic values, only for offline testing.\n')
    (path/'manifest.json').write_text(json.dumps({
        'format':'t3-reduced-final-report-v1',
        'primary_frozen_run':'reduced_frozen_overnight_v2',
        'diagnostic_running_run':'reduced_48_v1',
        'number_of_individual_accepted':32,
        'number_of_shared_accepted':10,
        'number_of_running_historical_scalar_screen_pass':0}))
    return path


def test_complete(tmp_path):
    d=report.load(synthetic(tmp_path/'input'))
    assert len(d['fits'])==32 and len(d['pulls'])==320 and len(d['running'])==16


def test_missing_one_frozen_case(tmp_path):
    source=synthetic(tmp_path/'input')
    x=report.read_csv(source/'frozen_comparison.csv')
    write_csv(source/'frozen_comparison.csv',x[:-1])
    with pytest.raises(ValueError,match='32 frozen'):report.load(source)


def test_reject_source_mismatch(tmp_path):
    source=synthetic(tmp_path/'input')
    x=json.loads((source/'manifest.json').read_text());x['number_of_shared_accepted']=11
    (source/'manifest.json').write_text(json.dumps(x))
    with pytest.raises(ValueError,match='manifest count'):report.load(source)


def test_bad_coupling(tmp_path):
    source=synthetic(tmp_path/'input')
    x=report.read_csv(source/'frozen_comparison.csv');x[0]['lambdaT3_input']=1.
    write_csv(source/'frozen_comparison.csv',x)
    with pytest.raises(ValueError,match='topology coupling'):report.load(source)


def test_bad_observable(tmp_path):
    source=synthetic(tmp_path/'input')
    x=report.read_csv(source/'oscillation_observable_pulls.csv');x[0]['prediction']=2.
    write_csv(source/'oscillation_observable_pulls.csv',x)
    with pytest.raises(ValueError,match='inconsistent observable pull'):report.load(source)


def test_no_overwrite(tmp_path):
    d=report.load(synthetic(tmp_path/'input'))
    out=tmp_path/'output';out.mkdir()
    with pytest.raises(ValueError,match='already exists'):report.build(d,out,html_only=True)


def test_html_and_latex(tmp_path):
    d=report.load(synthetic(tmp_path/'input'))
    out=report.build(d,tmp_path/'output',html_only=True)
    assert (out/'index.html').exists()
    assert (out/'optimised.html').exists()
    assert (out/'shared_a_m2.html').exists()
    assert (out/'running_uv.html').exists()
    assert len(list(out.glob('model_*.html')))==16
    assert (out/'benchmark_report.tex').exists()
    assert 'BuildLambda5e9Displays' in (out/'provenance.json').read_text()
    assert '32/32' in (out/'index.html').read_text()
    assert 'general scalar' not in (out/'index.html').read_text().lower()


def test_latex_escape():
    assert report.tex_escape('A_-2')==r'A\_-2'
    assert report.tex_escape('x & y')==r'x \& y'


def test_no_general_scalar_import():
    import inspect
    text=inspect.getsource(report)
    assert 'BuildGeneralScalar48FinalReports' not in text
