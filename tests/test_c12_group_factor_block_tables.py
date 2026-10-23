"""Structure-by-representation PDF table regressions."""
from pathlib import Path
from types import SimpleNamespace
import sys
import types

import pytest
import sympy as sp

from Reports import C12GroupFactorPDF as gf


def fake_projection(*, coeff: str, alpha: int, cls: str) -> dict:
    dims = {"A": (1, 3, 2), "B": (2, 2, 1), "C": (2, 2, 3), "D": (3, 1, 2), "E": (3, 3, 2)}
    s1, s2, df = dims[cls]
    return {
        "representation": {"dS1": s1, "dS2": s2, "dF": df, "alpha": alpha},
        "seed_support_component_count": 20,
        "generated_mixed_components_outside_c12_seed_support": [],
        "sectors": {"total": {"collinear_on_c12_seed_support": True, "coefficient": coeff}},
    }


def test_exact_coefficients_split_into_numeric_group_factors():
    assert gf._factor_terms('-3*g2**2 + 2*lambda12 - lambda12Adj/2 + ye*conjugate(ye)') == {
        'g2**2': -3,
        'lambda12': 2,
        'lambda12Adj': -sp.Rational(1, 2),
        'ye*conjugate(ye)': 1,
    }
    assert gf._factor_terms('-4*sqrt(3)*lambda12Inv2/3')['lambda12Inv2'] == -4*sp.sqrt(3)/3


def test_group_factors_comparison_blocks_three_columns(tmp_path, monkeypatch):
    paths = []
    mapping = {}
    for alpha, value in [(-3, '3*g2**2'), (-1, '-3*g2**2-4*sqrt(3)*lambda12Inv2/3+ye*conjugate(ye)'), (1, '4*g2**2')]:
        label = f'T3-B-{"m"+str(-alpha) if alpha<0 else "p"+str(alpha)}'
        d = tmp_path/label
        d.mkdir()
        for f in gf.INPUT_FILES:
            (d/f).write_text('{}', encoding='utf-8')
        mapping[d] = fake_projection(coeff=value, alpha=alpha, cls='B')
        paths.append((label, d))
    monkeypatch.setattr(gf, 'derive_from_saved_data', lambda d: mapping[d])
    section, statuses = gf.render_c12_pdf_section(paths)
    assert 'Run block 1: T3-B-m3 -- T3-B-p1' in section
    assert section.count(r'\begin{longtable}') == 1
    assert r'\lambda_{12,\mathrm{Inv2}}' in section
    assert r'C_{12}' in section
    assert r'- \frac{4 \sqrt{3}}{3}' in section or r'- \frac{4\sqrt{3}}{3}' in section or r'- \frac{4 \sqrt{3}}{3}' in section
    assert [s['representation']['alpha'] for s in statuses] == [-3, -1, 1]


def test_inserts_in_blank_report_and_preserves_existing_c5(tmp_path, monkeypatch):
    d = tmp_path/'T3-B-m1'; d.mkdir()
    for f in gf.INPUT_FILES:
        (d/f).write_text('{}',encoding='utf-8')
    monkeypatch.setattr(gf, 'derive_from_saved_data', lambda _: fake_projection(coeff='-3*g2**2',alpha=-1,cls='B'))
    tex = tmp_path/'GF_EFT_1_after_F.tex'
    tex.write_text(r'\documentclass{article}'+'\n'+r'\begin{document}'+'\n'+r'\section*{C5 group factors} Existing C5.'+'\n'+r'\end{document}',encoding='utf-8')
    first = gf.update_existing_group_factor_pdf(tex_path=tex,model_data=[('T3-B-m1',d)], compile_pdf=False)
    second = gf.update_existing_group_factor_pdf(tex_path=tex,model_data=[('T3-B-m1',d)], compile_pdf=False)
    source = tex.read_text()
    assert first['changed'] and not second['changed']
    assert 'Existing C5.' in source
    assert source.count(gf.BEGIN) == source.count(gf.END) == 1
    assert r'\begin{longtable}' in source


def test_no_preexisting_c12_text_required(tmp_path, monkeypatch):
    # New generator must never parse old per-model C12 tables.
    d = tmp_path/'T3-A-p0'; d.mkdir()
    for f in gf.INPUT_FILES:
        (d/f).write_text('{}',encoding='utf-8')
    monkeypatch.setattr(gf, 'derive_from_saved_data', lambda _: fake_projection(coeff='-2*lambda12Inv1',alpha=0,cls='A'))
    fragment, _ = gf.render_c12_pdf_section([('T3-A-p0',d)])
    assert r'\lambda_{12,\mathrm{Inv1}}' in fragment
    assert '\\subsection*{T3-A-p0}' not in fragment


def test_auto_discovery_can_use_saved_records(monkeypatch, tmp_path):
    fake = types.ModuleType('scripts.BuildCombinedAnalyticalReports')
    fake.reconstruct_records=lambda path: [SimpleNamespace(name='T3-B-m1', output_dir=tmp_path/'model')]
    monkeypatch.setitem(sys.modules, 'scripts.BuildCombinedAnalyticalReports', fake)
    assert gf._discover_models(tmp_path) == [('T3-B-m1', tmp_path/'model'/'data')]
