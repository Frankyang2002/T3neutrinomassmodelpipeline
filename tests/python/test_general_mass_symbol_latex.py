"""Regression for PDF compilation: generated mass symbols must balance braces."""
from pathlib import Path
import importlib.util


def load_adapter():
    path = Path(__file__).resolve().parents[2] / 'Reports' / 'GeneralSymbolicAdapter.py'
    spec = importlib.util.spec_from_file_location('t3_mass_latex_adapter', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_mass_symbols_have_balanced_braces():
    adapter = load_adapter()
    for name in ('muT3GeneralS1Sq', 'muT3GeneralS2Sq', 'mT3GeneralS12'):
        result = adapter.generated_symbol(name)
        assert result.count('{') == result.count('}'), (name, result)
        assert result.endswith('}}')
    assert adapter.generated_symbol('muT3GeneralS1Sq') == r'\mu_{\mathrm{T3GeneralS1Sq}}'
    assert adapter.generated_symbol('mT3GeneralS12') == r'm_{\mathrm{T3GeneralS12}}'


def test_full_notation_rewrite_mass_symbols():
    adapter = load_adapter()
    source = r'\text{muT3GeneralS1Sq}+\mathrm{muT3GeneralS2Sq}+mT3GeneralS12'
    result = adapter.format_general_latex(source)
    assert 'muT3General' not in result
    assert result.count('{') == result.count('}')
    assert r'\mu_{\mathrm{T3GeneralS1Sq}}' in result


def test_existing_quartic_format_unchanged():
    adapter = load_adapter()
    assert adapter.generated_symbol('lambdaT3GeneralQ110101I1') == r'\lambda_{\mathrm{Q}110101}^{(1)}'
