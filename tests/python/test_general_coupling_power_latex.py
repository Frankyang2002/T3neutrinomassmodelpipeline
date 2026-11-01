"""Regression for TeXForm coupling powers with short readable labels."""
from test_general_notation_index import adapter, registry_fixture


def test_general_coupling_power_rendering(tmp_path):
    a = adapter()
    registry_fixture(tmp_path)
    a.ROOT = tmp_path
    value = a.format_general_latex(r'\text{lambdaT3GeneralQ000211I1}^2 + 4 lambdaT3GeneralQ110101I2^*')
    assert r'{\lambda_{X_{1}}}^2' in value
    assert r'{\lambda_{X_{2}}^{(2)}}^*' in value


def test_mass_and_yukawa_labels_preserved(tmp_path):
    a = adapter()
    registry_fixture(tmp_path)
    a.ROOT = tmp_path
    value = a.format_general_latex(r'y2Conj^2 + \mathrm{muT3GeneralS1Sq}^2')
    assert r'{y_{2}^{\mathrm{conj}}}^2' in value
    assert r'{\mu_{\mathrm{T3GeneralS1Sq}}}^2' in value
