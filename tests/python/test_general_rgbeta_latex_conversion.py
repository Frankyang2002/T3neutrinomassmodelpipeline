"""Functional tests: saved RGBeta expressions can be typeset without rerunning RGBeta."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
from unittest.mock import patch


def load_writer():
    root=Path(__file__).resolve().parents[2]
    source=root/'scripts'/'WriteT3GeneralSymbolicPDFs.py'
    spec=importlib.util.spec_from_file_location('general_report_writer_under_test',source)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_existing_latex_passes_without_wolfram(tmp_path):
    w=load_writer()
    src=tmp_path/'uv_rgbeta.json'
    src.write_text('{}',encoding='utf-8')
    payload={'report_betas':{'g2':'-3*g2^3'},'report_beta_latex':{'g2':r'-3g_2^3'}}
    with patch.object(w.subprocess,'run',side_effect=AssertionError('No conversion needed')):
        assert w.obtain_beta_latex(payload,src,tmp_path/'pdfs','wolframscript')['g2']==r'-3g_2^3'


def test_missing_latex_uses_saved_inputform(tmp_path):
    w=load_writer()
    src=tmp_path/'uv_rgbeta.json'
    payload={'report_betas':{'g2':'-3*g2^3','lambdaH':'lambdaH^2+g2^4'}}
    src.write_text(json.dumps(payload),encoding='utf-8')
    dest=tmp_path/'UV'/'A_m2'
    def fake_run(args,**kw):
        assert str(src) == args[-2]
        assert 'format_saved_rgbeta.wl' in args[2]
        Path(args[-1]).write_text(json.dumps({'report_beta_latex':{
            'g2':r'-3g_2^3','lambdaH':r'\\lambda_H^2+g_2^4'}}),encoding='utf-8')
        return type('Result',(),{'returncode':0,'stdout':'','stderr':''})()
    with patch.object(w.shutil,'which',return_value='wolframscript'),patch.object(w.subprocess,'run',side_effect=fake_run):
        result=w.obtain_beta_latex(payload,src,dest,'wolframscript')
    assert set(result)==set(payload['report_betas'])
    assert 'ToExpression[expression,InputForm]' in (dest/'format_saved_rgbeta.wl').read_text()
    with patch.object(w.subprocess,'run',side_effect=AssertionError('Must reuse cache')):
        assert w.obtain_beta_latex(payload,src,dest,'wolframscript')==result


def test_native_beta_latex_compiles(tmp_path):
    w=load_writer()
    contents=w.preamble('T3 A_-2 UV','Saved one-loop beta functions')
    contents += [r'\begin{align*}',r'16\pi^2\frac{dg_2}{d\ln\mu}&=-3g_2^{3}',r'\end{align*}']
    pdf=w.compile_report(tmp_path,'UV',contents,'pdflatex')
    assert pdf.is_file() and pdf.stat().st_size>1000
