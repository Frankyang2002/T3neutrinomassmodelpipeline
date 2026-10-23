"""Regression tests for PDF replacement of the existing GF_EFT_1_after_F report."""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace, ModuleType
import sys

import pytest
import sympy as sp

from Reports import C12GroupFactorPDF as pdf
from RGE.group_factors.validation.C12RepresentationGroupFactors import project_c12_group_factors


def sample_result(*, g=3, rep=None):
    y = sp.conjugate(sp.Symbol('y1')) * sp.conjugate(sp.Symbol('y2'))
    mf = sp.Symbol('MF')
    entries = {(1,1,5,11): y / (4*mf), (1,1,6,12): -y / (4*mf)}
    beta = {
        '1,1,5,11': {'total': str((g*sp.Symbol('g2')**2)*y/(4*mf)),
                      'scalar_pair': str(-sp.sqrt(3)*sp.Symbol('lambda12Inv2')*y/(3*mf))},
        '1,1,6,12': {'total': str(-(g*sp.Symbol('g2')**2)*y/(4*mf)),
                      'scalar_pair': str(sp.sqrt(3)*sp.Symbol('lambda12Inv2')*y/(3*mf))},
    }
    value = project_c12_group_factors(entries, beta)
    value['representation'] = rep or {'dS1':2,'dS2':2,'dF':1,'alpha':-1,'shared_scalar':False}
    value['sources'] = {}
    return value


def make_inputs(tmp_path, name='source'):
    d=tmp_path/name
    d.mkdir()
    for filename in pdf.INPUT_FILES:
        (d/filename).write_text('{}',encoding='utf-8')
    return d


def make_tex(tmp_path):
    p = tmp_path/'GroupFactors'/'GF_EFT_1_after_F.tex'
    p.parent.mkdir(parents=True)
    p.write_text(r'\documentclass{article}'+'\n'+r'\usepackage{booktabs,longtable}'+'\n'
                 +r'\begin{document}'+'\n'+r'Existing group-factors data.'+'\n'
                 +r'\end{document}'+'\n',encoding='utf-8')
    return p


def test_updates_existing_tex_not_new_report_and_is_idempotent(tmp_path,monkeypatch):
    tex = make_tex(tmp_path)
    d=make_inputs(tmp_path)
    result=sample_result()
    monkeypatch.setattr(pdf,'derive_from_saved_data', lambda x:result)
    first=pdf.update_existing_group_factor_pdf(tex_path=tex,model_data=[('T3-B alpha=-1',d)],compile_pdf=False)
    rendered=tex.read_text()
    assert first['changed']
    assert rendered.count(pdf.BEGIN)==rendered.count(pdf.END)==1
    assert 'Existing group-factors data.' in rendered
    assert 'T3-B alpha=-1' in rendered
    assert 'g_2^{2}' in rendered
    assert 'No' not in first['models'][0]['status']
    assert not tex.with_suffix('.html').exists()
    second=pdf.update_existing_group_factor_pdf(tex_path=tex,model_data=[('T3-B alpha=-1',d)],compile_pdf=False)
    assert not second['changed']
    assert tex.read_text()==rendered


def test_two_representations_not_reusing_factors(tmp_path,monkeypatch):
    tex=make_tex(tmp_path)
    a=make_inputs(tmp_path,'A')
    b=make_inputs(tmp_path,'B')
    monkeypatch.setattr(pdf,'derive_from_saved_data', lambda path: sample_result(
        g=3 if Path(path).name=='A' else 7,
        rep={'dS1':2,'dS2':2,'dF':1,'alpha':-1,'shared_scalar':False}
        if Path(path).name=='A' else
        {'dS1':1,'dS2':3,'dF':2,'alpha':-2,'shared_scalar':False}))
    result=pdf.update_existing_group_factor_pdf(tex_path=tex,
        model_data=[('B model',a),('A model',b)],compile_pdf=False)
    text=tex.read_text()
    assert 'B model' in text and 'A model' in text
    assert '$g_2^{2}\\,C_{12}$' in text and ' & $3$ ' in text and ' & $7$ ' in text
    assert len(result['models'])==2


def test_missing_saved_tensor_input_reported_without_fabricating(tmp_path):
    tex=make_tex(tmp_path)
    d=tmp_path/'no_data'
    d.mkdir()
    with pytest.raises(ValueError, match='saved tensor data missing'):
        pdf.update_existing_group_factor_pdf(tex_path=tex,
            model_data=[('T3-B alpha=-1',d)],compile_pdf=False)


def test_rejects_missing_tex_and_bad_markers(tmp_path):
    with pytest.raises(FileNotFoundError):
        pdf.update_existing_group_factor_pdf(tex_path=tmp_path/'absent.tex',
                                             model_data=[],compile_pdf=False)
    with pytest.raises(ValueError, match='markers'):
        pdf.replace_pdf_tex_section('hello '+pdf.BEGIN+'\\end{document}', 'foo')
    with pytest.raises(ValueError, match='end-of-document'):
        pdf.replace_pdf_tex_section('hello world', 'foo')


def test_compile_invoked_for_changed_tex_only(tmp_path,monkeypatch):
    tex=make_tex(tmp_path)
    d=make_inputs(tmp_path)
    monkeypatch.setattr(pdf,'derive_from_saved_data', lambda x:sample_result())
    called=[]
    m=ModuleType('Reports.ReportGeneration')
    def compile_latex_document(p):
        called.append(p)
        Path(p).with_suffix('.pdf').write_bytes(b'%PDF-1.4\n')
    m.compile_latex_document=compile_latex_document
    monkeypatch.setitem(sys.modules,'Reports.ReportGeneration',m)
    pdf.update_existing_group_factor_pdf(tex_path=tex,
        model_data=[('T3-B',d)],compile_pdf=True)
    pdf.update_existing_group_factor_pdf(tex_path=tex,
        model_data=[('T3-B',d)],compile_pdf=True)
    assert called==[tex]
    assert tex.with_suffix('.pdf').exists()


def test_stage_refresh_only_real_f_first_reports(tmp_path,monkeypatch):
    f=make_tex(tmp_path)
    uv=f.with_name('GF_UV.tex');uv.write_text('untouched')
    wrong=f.with_name('GF_EFT_2_after_S1_S2.tex');wrong.write_text('untouched')
    datadir=tmp_path/'run'/'data';datadir.mkdir(parents=True)
    records=[SimpleNamespace(name='T3_B_alpha_m1',shared_scalar=False,output_dir=datadir.parent,
             summary={'EFTStages':[{'Label':'EFT_1_after_F',
                                     'IntegratedFields':['F'],
                                     'ActiveHeavyFields':['S1','S2']}]})]
    fake=ModuleType('Reports.GroupFactorReports')
    fake._records_at_stage=lambda recs,label: [r for r in recs if any(
        item['Label']==label for item in r.summary['EFTStages'])]
    fake._stage_is_f_first_eft1=lambda r,label: label=='EFT_1_after_F'
    monkeypatch.setitem(sys.modules,'Reports.GroupFactorReports',fake)
    called=[]
    monkeypatch.setattr(pdf,'update_existing_group_factor_pdf',lambda **k:called.append(k) or {'pdf':'expected.pdf'})
    out=pdf.refresh_stage_group_factor_pdfs(records,[uv,f,wrong])
    assert out==[{'pdf':'expected.pdf'}]
    assert called[0]['tex_path']==f
    assert called[0]['model_data']==[('T3_B_alpha_m1',datadir)]


def test_source_projection_flags_noncollinear():
    y=sp.conjugate(sp.Symbol('y1'))*sp.conjugate(sp.Symbol('y2'))
    mf=sp.Symbol('MF')
    entries={(1,1,5,11):y/mf,(1,1,6,12):-y/mf}
    beta={'1,1,5,11':{'total':str(y/mf)},'1,1,6,12':{'total':str(2*y/mf)}}
    result=project_c12_group_factors(entries,beta)
    assert result['sectors']['total']['coefficient'] is None
    assert not result['full_operator_closure_verified']
