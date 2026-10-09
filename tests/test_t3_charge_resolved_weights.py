"""Regression checks for charge-basis T3 loop coefficients.

Generated Matchete audits are optional integration fixtures, not source files.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest
import sympy as sp

from Numerical.diagnostics.ResolveT3ChargeWeights import run, transformed, rotations

FIXTURES = Path(os.environ.get(
    "T3_MULTIPLICITY_INPUT_DIR",
    str(Path(__file__).resolve().parents[1] / "Reports/output/full/multiplicity"),
))


def load(name: str):
    path = FIXTURES / name
    if not path.is_file():
        pytest.skip(f"Generated charge-resolved multiplicity fixture unavailable: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def audit():
    return run(load("component_support_audit.json"),
               load("matchete_generator_analysis.json"),
               load("cg_covariance_audit.json"),
               load("charge_channels.json"))


def test_16_models_weights():
    models = audit()["models"]
    assert len(models) == 16
    expected = {
        "A": (-2 * sp.sqrt(3) / 3, 1),
        "B": (-2 * sp.sqrt(3) / 3, 1),
        "C": (sp.Rational(2, 3), 1),
        "D": (-2 * sp.sqrt(3) / 3, 1),
        "E": (2 * sp.sqrt(6) / 3, 2),
    }
    for rec in models:
        total, count = expected[rec["class"]]
        assert sp.simplify(sp.sympify(rec["total_group_coefficient"]) - total) == 0
        assert rec["nonzero_weight_count"] == count
        assert rec["equal_nonzero_weights"]
        assert rec["verified_physical_loop_multiplicity"] is None
        assert rec["compensation_factor"] is None


def test_E_charge_channel_weights():
    for rec in audit()["models"]:
        if rec["class"] != "E":
            continue
        assert len(rec["channels"]) == 2
        for channel in rec["channels"]:
            assert sp.simplify(sp.sympify(channel["complex_weight"]) - sp.sqrt(6) / 3) == 0


def test_transforms_are_unitary_and_preserve_norm():
    us = rotations(load("matchete_generator_analysis.json"))
    case = load("component_support_audit.json")["cases"][0]
    for name, flags in [
        ("T3Y1CG", (True, True, False)),
        ("T3Y2CG", (True, False, True)),
        ("T3MixCG", (False, False, False, True)),
    ]:
        vertex = case["vertices"][name]
        rotated = transformed(vertex, flags, us)
        old_norm = sum(sp.sympify(c["cg_abs_squared"]) for c in vertex["components"])
        new_norm = sum(sp.conjugate(x) * x for x in rotated.values())
        assert sp.simplify(old_norm - new_norm) == 0
