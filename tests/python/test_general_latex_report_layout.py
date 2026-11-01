"""Operational tests for the per-model mathematical PDF writer."""
import importlib.util
from pathlib import Path
import shutil
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / 'scripts' / 'WriteT3GeneralSymbolicPDFs.py'
spec = importlib.util.spec_from_file_location('general_math_reports', SCRIPT)
writer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(writer)


def test_model_tag_and_latex_escaping():
    assert writer.tag('A_-2') == 'A_m2'
    assert writer.tag('E_+0') == 'E_p0'
    assert writer.esc('lambda_1') == r'lambda\_1'


def test_wolfram_converter_produces_termwise_texform():
    assert 'TeXForm' in writer.WOLFRAM_CONVERTER
    assert 'InputForm' in writer.WOLFRAM_CONVERTER
    assert 'Export' in writer.WOLFRAM_CONVERTER


@pytest.mark.skipif(shutil.which('pdflatex') is None, reason='pdflatex unavailable')
def test_real_math_pdf_compilation(tmp_path):
    content = writer.preamble('T3 A_-2 UV','Mathematical report')
    content += [r'\section*{Complete BSM UV Lagrangian}']
    content += writer.lagrangian_latex([r'\lambda_1(H^\dagger H)(S_1^\dagger S_1)',
                                        r'\lambda_2(S_1^\dagger S_1)^2'])
    content += [r'\section*{One-loop renormalisable beta functions}',
                r'\begin{align*}',r'16\pi^2\frac{d\lambda_1}{d\ln\mu} &= 3\lambda_1^2 \\',
                r'&\quad - 9g_2^2\lambda_1 \\',r'\end{align*}']
    pdf = writer.compile_report(tmp_path / 'UV' / 'A_m2', 'UV', content, 'pdflatex')
    assert pdf.is_file() and pdf.stat().st_size > 1500
    tex = pdf.with_suffix('.tex').read_text(encoding='utf-8')
    assert r'\begin{adjustbox}' in tex
    assert r'\lambda_1' in tex
    assert 'lstinputlisting' not in tex
