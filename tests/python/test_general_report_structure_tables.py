"""Presentation-only regression tests. No physics verification."""
from pathlib import Path
import importlib.util

SRC = Path(__file__).resolve().parents[2] / 'scripts' / 'WriteT3GeneralSymbolicPDFs.py'
spec = importlib.util.spec_from_file_location('t3_general_reports_under_test', SRC)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.paper_symbol = lambda s: s


def fixture():
    def term(structure, latex, coefficient):
        return dict(structure=structure, structure_latex=latex,
                    coefficient=str(coefficient), coefficient_latex=str(coefficient))
    data = {
        'A_-2': {'UV': {'lambdaX': [term('g2^2 lambdaX', r'g_2^2\lambda_X', -3),
                                             term('lambdaX^2', r'\lambda_X^2', 8)]},
                 'EFT1': {'lambdaX': [term('lambdaX^2', r'\lambda_X^2', 5)]}},
        'A_+0': {'UV': {'lambdaX': [term('g2^2 lambdaX', r'g_2^2\lambda_X', -5)]},
                 'EFT1': {'lambdaX': [term('lambdaX^2', r'\lambda_X^2', 7)]}},
    }
    fields = [dict(model=model, field=f, d_SU2='2', Y='1/2')
              for model in data for f in ('S1', 'S2', 'F')]
    return data, fields


def test_prefactor_union_and_zero_in_comparison():
    data, fields = fixture()
    tex = '\n'.join(module.factor_tex('UV', data, fields))
    assert r'\lambda_X^2' in tex
    assert r'$0$' in tex
    assert r'$-3$' in tex and r'$-5$' in tex
    assert r'\mathcal S' in tex
    assert 'A\\_-2' in tex


def test_eft1_independent_stage():
    data, fields = fixture()
    tex = '\n'.join(module.factor_tex('EFT1', data, fields))
    assert r'\lambda_X^2' in tex
    assert r'$5$' in tex and r'$7$' in tex
    assert r'g_2^2\lambda_X' not in tex


def test_native_rational_aggregation():
    data, fields = fixture()
    data['A_-2']['UV']['lambdaX'].append(dict(
        structure='g2^2 lambdaX', structure_latex=r'g_2^2\lambda_X',
        coefficient='1/2', coefficient_latex=r'\frac{1}{2}'))
    table = dict(module.build_factor_rows(data, 'UV'))['lambdaX']
    assert table['g2^2 lambdaX']['values']['A_-2'] == '-5/2'
    assert module.numeric_display('-5/2') == r'\frac{-5}{2}'


def test_compiles_typeset_math(tmp_path):
    import shutil
    import pytest
    if not shutil.which('pdflatex'):
        pytest.skip('pdflatex unavailable')
    data, fields = fixture()
    tex = module.factor_tex('UV', data, fields)
    pdf = module.compile_report(tmp_path, 'GF_UV', tex, 'pdflatex')
    assert pdf.is_file() and pdf.stat().st_size > 1500
