"""Regression guards for multiplicity-only full-study benchmark normalization."""
from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / "studies" / "FullT3Study.py"


def _study_constants():
    tree = ast.parse(STUDY.read_text(encoding="utf-8"))
    result = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in {"LAMBDA_T3_NORMALISATION", "MODELS", "COMPARISON_SCENARIOS"}:
                    result[target.id] = ast.literal_eval(node.value)
    return result


def test_only_class_e_has_multiplicity_compensation():
    assert _study_constants()["LAMBDA_T3_NORMALISATION"] == {
        "A": 1.0, "B": 1.0, "C": 1.0, "D": 1.0, "E": 0.5,
    }


def test_four_scenarios_and_sixteen_models():
    constants = _study_constants()
    assert len(constants["MODELS"]) == 16
    assert len(constants["COMPARISON_SCENARIOS"]) == 4
    assert len(set(constants["MODELS"])) == 16
    assert len(constants["MODELS"]) * len(constants["COMPARISON_SCENARIOS"]) == 64


def test_summary_exposes_multiplicity_basis():
    source = STUDY.read_text(encoding="utf-8")
    assert "charge-channel multiplicity only" in source
    assert '"physical_group_factors_modified": False' in source
    assert '"normalisation_basis": LAMBDA_T3_NORMALISATION_BASIS' in source
