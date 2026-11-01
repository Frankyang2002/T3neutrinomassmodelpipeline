"""Formatting-only regression for Matchete Extra vs RGBeta General names."""
import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PRESENTATION = ROOT / 'Reports' / 'GeneralSymbolicPresentation.py'


def test_extra_quartics_are_not_rejected():
    module = ast.parse(PRESENTATION.read_text(encoding='utf-8'))
    names = {node.name for node in module.body if isinstance(node, ast.FunctionDef)}
    assert {'generated_symbol', '_codes', 'format_general_latex'} <= names
    source = PRESENTATION.read_text(encoding='utf-8')
    assert 'Unregistered quartic' not in source
    assert "matchete_basis_name" in source
    assert "lambdaT3(?:General|Extra)" in source


def test_matching_extras_have_identical_short_symbol():
    # Exercise the real formatting function without loading Wolfram or the pipeline.
    source = PRESENTATION.read_text(encoding='utf-8')
    start = source.index('def generated_symbol(')
    end = source.index('\n\ndef format_general_latex(', start)
    import re
    scope = {'QUARTIC': re.compile(r'lambdaT3(?:General|Extra)Q([0-9]{6})I([0-9]+)'),
             'MASS': re.compile(r'(?:mu|m)T3[A-Za-z0-9_]+'),
             '_codes': lambda: {'000202': (3, False), '110101': (8, True)}}
    exec(source[start:end], scope)
    symbol = scope['generated_symbol']
    assert symbol('lambdaT3ExtraQ000202I1') == symbol('lambdaT3GeneralQ000202I1')
    assert symbol('lambdaT3ExtraQ110101I2') == r'\lambda_{X_{8}}^{(2)}'
    assert symbol('lambdaT3ExtraQ999999I1').startswith(r'\lambda_{X_')
