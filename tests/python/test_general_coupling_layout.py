"""Source-independent tests for the in-pipeline general coupling registry."""
import numpy as np
import pytest

from Numerical.core.GeneralCouplingLayout import (
    fermion_first_projection, layout_from_rgbeta, validate_general_values,
)


def test_general_inventory_preserves_complex_and_matrix_values():
    payload = {"report_betas": {
        "g2": "-g2^3", "y1Conj": "y1Conj", "MF": "MF",
        "lambdaT3GeneralQ000400I1": "lambdaT3GeneralQ000400I1",
        "muT3GeneralS1Sq": "muT3GeneralS1Sq",
    }}
    a = np.eye(3, dtype=complex) * (0.1 + 0.2j)
    m = np.eye(3, dtype=complex) * 1000
    values = {
        "g2": 0.65, "y1Conj": a, "MF": m,
        "lambdaT3GeneralQ000400I1": 0.01 + 0.03j,
        "muT3GeneralS1Sq": 5 + 7j,
    }
    layout = layout_from_rgbeta(payload)
    assert {x.name for x in layout} == set(values)
    out = validate_general_values(payload, values, majorana=True)
    assert out["muT3GeneralS1Sq"].imag == 7
    assert np.array_equal(out["y1Conj"], a)
    with pytest.raises(ValueError, match="mismatch"):
        validate_general_values(payload, {**values, "unregistered": 1})
    bad = dict(values, MF=np.ones((3, 3)) + 1j * np.triu(np.ones((3, 3)), 1))
    with pytest.raises(ValueError, match="symmetric"):
        validate_general_values(payload, bad, majorana=True)


def test_eft1_projection_keeps_added_scalar_couplings():
    eft1 = {"report_betas": {"g2": "g2", "muT3GeneralS1S2": "muT3GeneralS1S2"}}
    result = fermion_first_projection(eft1, {
        "g2": 0.63, "muT3GeneralS1S2": 2 + 3j, "MF": np.eye(3)
    })
    assert set(result) == {"g2", "muT3GeneralS1S2"}
    assert result["muT3GeneralS1S2"] == 2 + 3j


def test_unknown_source_key_is_not_silently_ignored():
    with pytest.raises(ValueError, match="Unrecognised"):
        layout_from_rgbeta({"report_betas": {"mystery": "mystery"}})
