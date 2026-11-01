"""Presentation regression: general coupling names stay readable and uniquely indexed."""
import importlib.util
from pathlib import Path
from unittest.mock import patch


def adapter():
    path=Path(__file__).resolve().parents[2]/'Reports/GeneralSymbolicAdapter.py'
    spec=importlib.util.spec_from_file_location('general_adapter_notation_test',path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_generated_quartic_symbols_and_conjugates():
    a=adapter()
    assert a.generated_symbol('lambdaT3GeneralQ110101I1') == r'\lambda_{\mathrm{Q}110101}^{(1)}'
    assert a.generated_symbol('lambdaT3GeneralQ110101I2') != a.generated_symbol('lambdaT3GeneralQ110101I1')
    assert a.generated_symbol('y2Conj') == r'y_{2}^{\mathrm{conj}}'
    assert a.generated_symbol('muT3GeneralS2Sq').startswith(r'\mu_')


def test_math_token_normalisation():
    a=adapter()
    tex=r'4lambdaT3GeneralQ110101I1 (lambdaT3GeneralQ110101I1)^*+\text{lambdaT3GeneralQ110101I2}'
    out=a.format_general_latex(tex)
    assert 'lambdaT3GeneralQ' not in out
    assert r'\lambda_{\mathrm{Q}110101}^{(1)}' in out
    assert r'\lambda_{\mathrm{Q}110101}^{(2)}' in out


def test_coupling_index_source_registry(monkeypatch):
    a=adapter()
    sample={'new_rgbeta_symbols':[{'name':'lambdaT3GeneralQ110101I1',
         'counts':[1,1,0,1,0,1], 'invariant':'del[SU2L, #1, #2]'}],
         'preexisting_rgbeta_symbols':[]}
    import sys, types
    gen=types.ModuleType('Reports.ReportGeneration')
    gen.latex_escape_text=a.esc
    gen.paper_symbol_latex=lambda x:x
    gen.paper_notation_key_lines=lambda **kwargs:[r'\section*{Coupling and interaction notation}']
    with patch.dict(sys.modules,{'Reports':types.ModuleType('Reports'),'Reports.ReportGeneration':gen}):
        with patch.object(a,'registry_entries',lambda _:sample):
            lines=a.notation_index('B_+1',{'interactions':['y2Conj']},['lambdaT3GeneralQ110101I1'])
    all_tex='\n'.join(lines)
    assert r'\lambda_{\mathrm{Q}110101}^{(1)}' in all_tex
    assert r'S_1^\dagger' in all_tex and r'S_2^\dagger' in all_tex
    assert r'y_{2}^{\mathrm{conj}}' in all_tex
    assert r'\begin{longtable}' in all_tex
