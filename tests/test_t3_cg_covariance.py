"""CG covariance audit tests; generated Matchete fixtures are optional."""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from Numerical.diagnostics.VerifyT3CGCovariance import audit, generators, tensor_residual

DATA = Path(os.environ.get(
    "T3_MULTIPLICITY_INPUT_DIR",
    str(Path(__file__).resolve().parents[1] / "Reports/output/full/multiplicity"),
))


def load(name: str):
    path = DATA / name
    if not path.is_file():
        pytest.skip(f"Generated Matchete multiplicity fixture unavailable: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def test_generator_algebra():
    gs = generators(load("matchete_generator_probe.json"))
    assert set(gs) == {1, 2, 3}


def test_all_exported_cg_tensors_are_covariant():
    component = load("component_support_audit.json")
    probe = load("matchete_generator_probe.json")
    result = audit(component, probe)
    assert result["count"] == 16
    assert all(x["all_covariant"] for x in result["models"])
    assert all(x["verified_loop_multiplicity"] is None for x in result["models"])


def test_corrupt_tensor_is_not_covariant():
    component = load("component_support_audit.json")
    probe = load("matchete_generator_probe.json")
    t = dict(component["cases"][0]["vertices"]["T3MixCG"])
    t["components"] = [dict(x) for x in t["components"]]
    t["components"][0]["cg"] = "99"
    assert tensor_residual(t, (False, False, False, True), generators(probe))
