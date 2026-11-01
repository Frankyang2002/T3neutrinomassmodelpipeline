"""Regression: model-independent, readable symbols for general quartics."""
import importlib.util
import json
from pathlib import Path
from unittest.mock import patch


def adapter():
    path = Path(__file__).resolve().parents[2] / 'Reports' / 'GeneralSymbolicAdapter.py'
    spec = importlib.util.spec_from_file_location('general_adapter_xlabels', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def registry_fixture(tmp_path):
    p = tmp_path / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
    p.parent.mkdir(parents=True)
    p.write_text(json.dumps({'models': {
        'A_+0': {'new_rgbeta_symbols': [
            {'name': 'lambdaT3GeneralQ000211I1'},
            {'name': 'lambdaT3GeneralQ110101I1'},
            {'name': 'lambdaT3GeneralQ110101I2'},
        ]},
        'B_+1': {'new_rgbeta_symbols': [
            {'name': 'lambdaT3GeneralQ110101I1'},
        ]},
    }}), encoding='utf-8')


def test_compact_and_unique_labels(tmp_path):
    a = adapter()
    registry_fixture(tmp_path)
    a.ROOT = tmp_path
    assert a.generated_symbol('lambdaT3GeneralQ000211I1') == r'\lambda_{X_{1}}'
    assert a.generated_symbol('lambdaT3GeneralQ110101I1') == r'\lambda_{X_{2}}^{(1)}'
    assert a.generated_symbol('lambdaT3GeneralQ110101I2') == r'\lambda_{X_{2}}^{(2)}'
    assert a.generated_symbol('y2Conj') == r'y_{2}^{\mathrm{conj}}'


def test_powers_and_conjugates_use_braces(tmp_path):
    a = adapter()
    registry_fixture(tmp_path)
    a.ROOT = tmp_path
    tex = r'\text{lambdaT3GeneralQ110101I1}^2+lambdaT3GeneralQ110101I2^*'
    result = a.format_general_latex(tex)
    assert r'{\lambda_{X_{2}}^{(1)}}^2' in result
    assert r'{\lambda_{X_{2}}^{(2)}}^*' in result
    assert 'Q110101' not in result


def test_registry_order_does_not_affect_labels(tmp_path):
    a = adapter()
    registry_fixture(tmp_path)
    a.ROOT = tmp_path
    first = a._general_quartic_labels()
    p = tmp_path / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
    data = json.loads(p.read_text())
    data['models']['A_+0']['new_rgbeta_symbols'].reverse()
    p.write_text(json.dumps(data))
    assert a._general_quartic_labels() == first


def test_unregistered_name_fails_closed(tmp_path):
    import pytest
    a = adapter()
    registry_fixture(tmp_path)
    a.ROOT = tmp_path
    with pytest.raises(ValueError, match='Unregistered'):
        a.generated_symbol('lambdaT3GeneralQ999999I1')
