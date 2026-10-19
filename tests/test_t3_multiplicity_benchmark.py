"""Frozen-reference provenance and class-E multiplicity regression checks.

The current 16-point study reuses validated, frozen inputs; it does not build
four scenarios or reapply class-E scalar multiplicity compensation.
"""
from __future__ import annotations

import json
from pathlib import Path

from studies import FullT3Study

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "Numerical" / "benchmarks" / "fitted_reference" / "benchmark_manifest.json"
CLASS_DIMENSIONS = {
    "A": (1, 3, 2), "B": (2, 2, 1), "C": (2, 2, 3),
    "D": (3, 1, 2), "E": (3, 3, 2),
}


def test_fitted_reference_manifest_has_16_distinct_ordinary_models() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    rows = manifest["models"]
    assert len(rows) == 16
    assert len({row["model"] for row in rows}) == 16
    represented = {
        (row["representation"]["d_s1"],
         row["representation"]["d_s2"],
         row["representation"]["d_f"])
        for row in rows
    }
    assert represented == set(CLASS_DIMENSIONS.values())
    assert all(row["representation"].get("shared_scalar") is False for row in rows)


def test_fitted_reference_configs_are_present_and_validated() -> None:
    # Uses the authoritative study loader, including its SHA-256, representation,
    # disabled-fit and disabled-sensitivity checks.
    models = FullT3Study._load_models()
    assert len(models) == 16
    assert len({name for name, _, _ in models}) == 16
    assert all(config.is_file() for _, _, config in models)


def test_multiplicity_is_not_reapplied_to_frozen_class_e_inputs() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    note = manifest.get("class_e_note", "")
    assert "No extra multiplicity factor" in note
    assert "representation-specific pipeline" in note
    e_rows = [row for row in manifest["models"] if (
        row["representation"]["d_s1"],
        row["representation"]["d_s2"],
        row["representation"]["d_f"],
    ) == CLASS_DIMENSIONS["E"]]
    assert e_rows
    # This tests the full-study layer's behaviour, not an independent derivation
    # of the group-theory multiplicity appearing in the matching calculation.
    source = (ROOT / "studies" / "FullT3Study.py").read_text(encoding="utf-8")
    assert "_fixed_comparison_config" not in source
    assert "No extra multiplicity factor" in note
