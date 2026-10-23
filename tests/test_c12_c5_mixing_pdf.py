"""Explicit regression for the two off-diagonal Wilson RGEs in the GF PDF."""
from __future__ import annotations

from pathlib import Path
import sympy as sp

from Reports import C12GroupFactorPDF as gf


def projected(*, alpha: int, label: str, collinear: bool = True) -> dict:
    return {
        "model": label,
        "status": "Projected",
        "representation": {
            "dS1": 2, "dS2": 2, "dF": 1, "alpha": alpha,
            "shared_scalar": False,
        },
        "terms": {"g2**2": sp.Integer(-3)},
        "total_collinear": collinear,
    }


def test_verified_both_t3b_m1_off_diagonal_group_factors():
    item = projected(alpha=-1, label="T3-B-m1")
    assert sp.simplify(gf._forward_mixing_factor(item) - 4*sp.sqrt(3)/3) == 0
    assert sp.simplify(gf._reverse_mixing_factor(item) - 8*sp.sqrt(3)/3) == 0


def test_other_hypercharge_reverse_factor_unknown_not_zero():
    item = projected(alpha=1, label="T3-B-p1")
    assert gf._reverse_mixing_factor(item) is None
    tex = "\n".join(gf._comparison_table([item], 1))
    assert r"$\lambda_5^* C_5\to\beta_{C_{12}}$" in tex
    assert r"\textit{not verified}" in tex
    assert r"\beta_{\perp}" not in tex


def test_noncollinear_diagnostic_explained_outside_table():
    item = projected(alpha=-1, label="T3-B-m1", collinear=False)
    tex = "\n".join(gf._comparison_table([item], 1))
    assert "Projection status:" in tex
    assert "not collinear" in tex
    assert r"\beta_{\perp}" not in tex


def test_missing_reverse_input_does_not_silently_claim_zero():
    item = projected(alpha=-1, label="T3-B-m1")
    item["status"] = "MissingSavedTensorRGE"
    assert gf._forward_mixing_factor(item) is None
    assert gf._reverse_mixing_factor(item) is None
