"""Regression: only reverse C5 -> C12 mixing belongs in C12 PDF table."""
from __future__ import annotations
import sympy as sp
from Reports import C12GroupFactorPDF as gf


def projected(alpha=-1, label="T3-B-m1"):
    return {"model": label, "status": "Projected",
            "representation": {"dS1": 2, "dS2": 2, "dF": 1,
                               "alpha": alpha, "shared_scalar": False},
            "terms": {"g2**2": sp.Integer(-3)}, "total_collinear": True}


def test_c5_forward_row_removed_from_c12_table():
    item = projected()
    tex = "\n".join(gf._comparison_table([item], 1))
    assert r"\lambda_5^* C_5\to\beta_{C_{12}}" in tex
    assert r"\lambda_5 C_{12}\to\beta_{C_5}" not in tex
    assert r"\beta_{\perp}" not in tex
    assert tex.count("Reverse") == 1


def test_previous_t3_b_m1_reverse_coefficient_preserved():
    assert sp.simplify(gf._reverse_mixing_factor(projected()) - 8*sp.sqrt(3)/3) == 0


def test_other_reverse_factor_stays_unknown_until_audit():
    item = projected(alpha=1, label="T3-B-p1")
    assert gf._reverse_mixing_factor(item) is None
    assert r"\textit{not verified}" in "\n".join(gf._comparison_table([item], 1))
