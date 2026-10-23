from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
import sys

import pytest
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from RGE.group_factors.validation.C12RepresentationGroupFactors import (
    derive_from_saved_data,
    project_c12_group_factors,
)
from Reports.C12GroupFactorHTML import (
    BEGIN, END, insert_c12_into_existing_html, render_c12_section,
    replace_existing_section,
)


def sample():
    y = sp.conjugate(sp.Symbol("y1")) * sp.conjugate(sp.Symbol("y2"))
    mf = sp.Symbol("MF")
    c = y / (4 * mf)
    r = -3*sp.Symbol('g2')**2 - 4*sp.sqrt(3)*sp.Symbol('lambda12Inv2')/3 + sp.Symbol('ye')*sp.conjugate(sp.Symbol('ye'))
    entries = {(1,1,5,11): c, (1,1,6,12): -c}
    beta = {
        '1,1,5,11': {'total': str(sp.expand(c*r)), 'scalar_pair': str(c*(-4*sp.sqrt(3)*sp.Symbol('lambda12Inv2')/3))},
        '1,1,6,12': {'total': str(sp.expand(-c*r)), 'scalar_pair': str(-c*(-4*sp.sqrt(3)*sp.Symbol('lambda12Inv2')/3))},
        '1,1,3,3': {'total': '2*sqrt(3)*lambdaT3*conjugate(y1)*conjugate(y2)/(3*MF)'},
    }
    return entries, beta


def test_exact_T3B_C12_factor_derived_from_source_form():
    entries, beta = sample()
    result = project_c12_group_factors(entries, beta)
    assert result['seed_support_component_count'] == 2
    assert result['sectors']['total']['collinear_on_c12_seed_support']
    expected = -3*sp.Symbol('g2')**2 - 4*sp.sqrt(3)*sp.Symbol('lambda12Inv2')/3 + sp.Symbol('ye')*sp.conjugate(sp.Symbol('ye'))
    assert sp.simplify(sp.sympify(result['sectors']['total']['coefficient']) - expected) == 0
    assert [1,1,3,3] in result['generated_mixed_components_outside_c12_seed_support']
    assert not result['full_operator_closure_verified']


def test_noncollinear_beta_is_flagged_not_reported_as_unique():
    entries, beta = sample()
    beta['1,1,6,12']['total'] = str(sp.sympify(beta['1,1,6,12']['total']) + sp.conjugate(sp.Symbol('y1'))*sp.conjugate(sp.Symbol('y2')) / sp.Symbol('MF'))
    result = project_c12_group_factors(entries, beta)
    assert result['sectors']['total']['coefficient'] is None
    assert len(result['sectors']['total']['failures']) == 1


def test_zero_or_unmatched_seed_rejected():
    with pytest.raises(ValueError, match='no mixed'):
        project_c12_group_factors({(1,1,5,11):sp.Symbol('y1')**2}, {})


def _derive_mocked_real_adapter_contract_and_paper_conversion(tmp_path, monkeypatch):
    entries, beta = sample()
    d = tmp_path / 'data'
    d.mkdir()
    (d / 'eft1_after_F_wilson_seed.json').write_text('{}')
    (d / 'eft1_rgbeta_rge.json').write_text(json.dumps({'metadata':{
        'dS1':2, 'dS2':2, 'dF':1, 'alpha':-1, 'SharedScalar':False}}))
    (d / 'eft1_wilson_rge.json').write_text(json.dumps({
        'status':'Success','scan_mode':'L/eC one-loop closure',
        'metadata':{'dS1':2,'dS2':2},'betas':beta}))
    from types import ModuleType
    fake_adapter = ModuleType('RGE.running.intermediate.ScalarOnlyTensorAdapters')
    fake_adapter.load_and_build_eft1_wilson_tensor = lambda *args, **kwargs: SimpleNamespace(entries=entries)
    monkeypatch.setitem(sys.modules, fake_adapter.__name__, fake_adapter)
    derived = derive_from_saved_data(d)
    assert derived['t3b_paper_basis']['matches_previous_two_operator_diagonal_coefficient']
    assert derived['representation']['alpha'] == -1
    return derived


def test_existing_html_is_updated_in_place_and_idempotent(tmp_path, monkeypatch):
    derived = _derive_mocked_real_adapter_contract_and_paper_conversion(tmp_path, monkeypatch)
    from Reports import C12GroupFactorHTML as module
    monkeypatch.setattr(module, 'derive_from_saved_data', lambda _dir: derived)
    file = tmp_path / 'GF_EFT_1_after_F.html'
    original = '<!DOCTYPE html><html><body><main><h1>Existing group factor report</h1><p>keep everything</p></main><script>window.x=4</script></body></html>'
    file.write_text(original)
    one = insert_c12_into_existing_html(html_path=file, data_dir=tmp_path/'data', label='T3-B')
    text = file.read_text()
    assert one['total_beta_collinear'] is True
    assert text.count(BEGIN) == text.count(END) == 1
    assert '<p>keep everything</p>' in text and '<script>window.x=4</script>' in text
    insert_c12_into_existing_html(html_path=file, data_dir=tmp_path/'data', label='T3-B')
    assert file.read_text() == text
    assert '<math' in text


def test_insert_rejects_missing_existing_html(tmp_path):
    with pytest.raises(FileNotFoundError):
        insert_c12_into_existing_html(html_path=tmp_path/'missing.html', data_dir=tmp_path)


def test_replace_rejects_broken_markers():
    with pytest.raises(ValueError, match='marker'):
        replace_existing_section('<html><body>' + BEGIN + '</body></html>', BEGIN + END)


def test_html_writer_uses_same_existing_group_factor_file(tmp_path, monkeypatch):
    import Reports.GroupFactorCoupledHTML as writer
    from Reports import C12GroupFactorHTML as inserter
    record_path = tmp_path/'output'
    data_dir = record_path / 'data'
    data_dir.mkdir(parents=True)
    for name in ('eft1_wilson_rge.json', 'eft1_after_F_wilson_seed.json', 'eft1_rgbeta_rge.json'):
        (data_dir/name).write_text('{}')
    tex = tmp_path/'GroupFactors'/'GF_EFT_1_after_F.tex'
    tex.parent.mkdir()
    tex.write_text('existing TeX')
    html = tex.with_suffix('.html')
    html.write_text('<html><body><main>PRE-EXISTING HTML</main></body></html>')
    entries,beta=sample()
    projection=project_c12_group_factors(entries,beta)
    projection['representation']={'dS1':2,'dS2':2,'dF':1,'alpha':-1,'shared_scalar':False}
    projection['sources']={}
    monkeypatch.setattr(inserter,'derive_from_saved_data', lambda x: projection)
    record = SimpleNamespace(d_s1=2,d_s2=2,d_f=1,alpha=-1,shared_scalar=False,name='T3-B',output_dir=record_path,
                             summary={'EFTStages':[{'Label':'EFT_1_after_F','IntegratedFields':['F'],'ActiveHeavyFields':['S1','S2']}]})
    result=writer.write_group_factor_html_companions([record],[tex],preserve_existing_html=True)
    assert result == [html]
    text=html.read_text()
    assert 'PRE-EXISTING HTML' in text
    assert 'Representation-specific' in text
    assert text.count(BEGIN)==1


def test_exact_saved_contract_and_paper_conversion(tmp_path, monkeypatch):
    result = _derive_mocked_real_adapter_contract_and_paper_conversion(tmp_path, monkeypatch)
    assert result['t3b_paper_basis']['matches_previous_two_operator_diagonal_coefficient']


def test_multi_representation_sections_retain_distinct_source_factors(tmp_path, monkeypatch):
    from Reports import C12GroupFactorHTML as module
    entries, beta = sample()
    first = project_c12_group_factors(entries, beta)
    first['representation'] = {'dS1':2,'dS2':2,'dF':1,'alpha':-1,'shared_scalar':False}
    first['sources'] = {}
    second = project_c12_group_factors(entries, beta)
    second['representation'] = {'dS1':1,'dS2':3,'dF':2,'alpha':-2,'shared_scalar':False}
    second['sources'] = {}
    # These data are intentionally illustrative: verify layout keeps models
    # distinct rather than copying the first model's representation metadata.
    monkeypatch.setattr(module, 'derive_from_saved_data',
                        lambda directory: first if Path(directory).name == 'A' else second)
    html = tmp_path/'GF_EFT_1_after_F.html'
    html.write_text('<html><body><main><h1>Already existing</h1></main></body></html>')
    sources = [('Model B',tmp_path/'A'),('Model A',tmp_path/'B')]
    result = module.insert_many_c12_into_existing_html(html_path=html, model_data=sources)
    assert result[0]['representation']['alpha'] == -1
    assert result[1]['representation']['alpha'] == -2
    text = html.read_text()
    assert text.count(BEGIN) == text.count(END) == 1
    assert text.count('class="c12-representation-group-factors"') == 2
    assert 'Model B' in text and 'Model A' in text
    module.insert_many_c12_into_existing_html(html_path=html, model_data=sources)
    assert html.read_text() == text
