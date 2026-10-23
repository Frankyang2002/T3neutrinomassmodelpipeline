import json
from pathlib import Path
import sympy as sp
import pytest

from validation.ValidateC12C5Mixing import extract_matrix, validate

@pytest.fixture
def audit():
    # Small exact regression derived from the real T3-B alpha=-1 audit.
    def block(ratio):
        return {"proportionality_factor": ratio,
                "proportionality_residual_count": 0,
                "nonzero_target_components": 1}
    empty = {"proportionality_factor": None,
             "proportionality_residual_count": 0,
             "nonzero_target_components": 0}
    return {"mixing": {
      "C12_only": {
        "LLSS": block("-3*g2**2-4*sqrt(3)*lambda12Inv2/3+ye*conjugate(ye)"),
        "C5_HH": block("4*sqrt(3)*lambdaT3*conjugate(y1)*conjugate(y2)/(3*MF*kappa)"),
        "mixed_HS": empty},
      "C5_only": {
        "LLSS": block("8*sqrt(3)*MF*kappa*conjugate(lambdaT3)/(3*conjugate(y1)*conjugate(y2))"),
        "C5_HH": block("-3*g2**2+6*yd*conjugate(yd)-ye*conjugate(ye)+6*yu*conjugate(yu)+2*λ"),
        "mixed_HS": empty}
    }}

def test_matrix_extraction(audit):
    result = validate(audit)
    assert result["status"] == "ChecksPassed"
    assert result["higgs_quartic_mapping"]["status"] == "requires_explicit_mapping"
    assert all(result["checks"].values())

def test_production_consistency(audit):
    prod = {"weinberg_subspace_validation":
        {"beta_kappa_16pi2": "4*sqrt(3)*lambdaT3*conjugate(y1)*conjugate(y2)/(3*MF)"}}
    assert validate(audit, prod)["checks"]["production_forward_matches_audit"]

def test_detects_bad_forward(audit):
    audit["mixing"]["C12_only"]["C5_HH"]["proportionality_factor"] = "0"
    assert validate(audit)["status"] == "ChecksFailed"

def test_rejects_nonclosure(audit):
    audit["mixing"]["C5_only"]["LLSS"]["proportionality_residual_count"] = 2
    with pytest.raises(ValueError, match="Non-closed"):
        extract_matrix(audit)
