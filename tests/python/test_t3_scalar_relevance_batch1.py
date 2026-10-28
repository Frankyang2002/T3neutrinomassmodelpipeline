"""Offline tests for the Batch-1 relevance audit (not Matchete/RGBeta tests)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.AuditT3ScalarRelevance import (
    AuditError,
    classify_sources,
    counts_key,
    load_sources,
    tag,
    write_reports,
)


def sample_sources() -> dict:
    # Two-model synthetic fixture, including exact B(-1) manifest counts.
    models = {
        "A_-4": {
            "extra_quartic_complex_couplings": 0,
            "quadratic_complex_couplings_not_enabled": 0,
            "extra_quartic_structures": [],
            "quadratic_structures_not_enabled": [],
        },
        "B_-1": {
            "extra_quartic_complex_couplings": 7,
            "quadratic_complex_couplings_not_enabled": 1,
            "extra_quartic_structures": [
                {"counts": [0, 0, 0, 1, 1, 2], "field_structure": "S1dag S2 S2dag^2", "invariants": 1},
                {"counts": [0, 0, 0, 2, 0, 2], "field_structure": "S1dag^2 S2dag^2", "invariants": 1},
                {"counts": [0, 0, 1, 2, 0, 1], "field_structure": "S1 S1dag^2 S2dag", "invariants": 1},
                {"counts": [0, 2, 0, 0, 2, 0], "field_structure": "Hdag^2 S2^2", "invariants": 1},
                {"counts": [0, 2, 0, 2, 0, 0], "field_structure": "Hdag^2 S1dag^2", "invariants": 1},
                {"counts": [1, 1, 0, 1, 0, 1], "field_structure": "H Hdag S1dag S2dag", "invariants": 2},
            ],
            "quadratic_structures_not_enabled": [
                {"counts": [0, 0, 0, 1, 0, 1], "field_structure": "S1dag S2dag", "invariants": 1},
            ],
        },
    }
    observed = ["lambdaT3ExtraQ020020I1", "lambdaT3ExtraQ020200I1"]
    physics = []
    matches = []
    audits = {}
    quadratics = {}
    for model, row in models.items():
        is_b = model == "B_-1"
        code = tag(model)
        physics.append({"model": code, "zero_extra_limit_proven": True,
                        "extra_c5_couplings": observed if is_b else []})
        matches.append({
            "model": code,
            "accepted_extra_quartics_in_extended_log": row["extra_quartic_complex_couplings"],
            "extended_log_operator_count_matches_manifest": True,
            "variants": {
                "baseline": {"matched_c5_artifact_available": True},
                "extended": {"matched_c5_artifact_available": True,
                             "extra_couplings_observed": observed if is_b else []},
            },
        })
        audits[model] = {"success": True, "baseline_valid": True,
                         "combined_valid": True,
                         "accepted_extra_quartics": row["extra_quartic_complex_couplings"]}
        quadratics[model] = {
            "independent_complex_mass_parameters_in_declared_basis":
                row["quadratic_complex_couplings_not_enabled"]}
    return {
        "manifest": {"models": models},
        "matching": {"results": matches},
        "physics": {"symbolically_consistent": True,
                    "all_zero_extra_limits_formally_proven": True, "models": physics},
        "audit": {"models": audits},
        "quadratic": {"models": quadratics},
    }


def test_coupling_naming_and_model_tags():
    assert counts_key([0, 2, 0, 0, 2, 0]) == "020020"
    assert tag("A_+0") == "A_p0"
    assert tag("B_-1") == "B_m1"


def test_two_model_classification_no_false_higher_order_claim():
    report = classify_sources(sample_sources(), strict=False)
    summary = report["summary"]
    assert summary["extra_quartic_slots"] == 7
    assert summary["extra_quadratic_slots"] == 1
    assert summary["observed_models"] == ["B_-1"]
    assert summary["observed_additional_quartic_C5_slots"] == 2
    extras = [r for r in report["operators"] if r["kind"] == "additional_quartic"]
    assert len(extras) == 7
    assert len({r["coupling"] for r in extras}) == 7
    assert len([r for r in extras if r["classification"] == "observed_in_saved_one_loop_C5"]) == 2
    assert len([r for r in extras if r["classification"] == "higgs_vertex_not_seen_in_saved_C5"]) == 2
    assert len([r for r in extras if r["classification"] == "no_direct_single_vertex_LLSS_to_LLHH"]) == 3
    quadratic = next(r for r in report["operators"] if r["kind"] == "additional_quadratic")
    assert quadratic["classification"] == "unsuppressed_tree_scalar_mass_effect"
    assert "not suppressed" in quadratic["physics_reason"]


def test_missing_source_refused(tmp_path: Path):
    with pytest.raises(AuditError, match="Missing required source"):
        load_sources(tmp_path)


def test_incomplete_symbolic_c5_refused():
    src = sample_sources()
    next(x for x in src["matching"]["results"] if x["model"] == "B_m1")["variants"]["extended"]["matched_c5_artifact_available"] = False
    with pytest.raises(AuditError, match="Missing extended C5 artifact"):
        classify_sources(src, strict=False)


def test_unmapped_symbolic_coupling_refused():
    src = sample_sources()
    row = next(x for x in src["matching"]["results"] if x["model"] == "B_m1")
    row["variants"]["extended"]["extra_couplings_observed"].append("lambdaT3ExtraQ999999I1")
    src["physics"]["models"][1]["extra_c5_couplings"].append("lambdaT3ExtraQ999999I1")
    with pytest.raises(AuditError, match="absent from the inventory"):
        classify_sources(src, strict=False)


def test_zero_extra_unproven_refused():
    src = sample_sources()
    src["physics"]["models"][0]["zero_extra_limit_proven"] = False
    with pytest.raises(AuditError, match="Zero-extra limit unproven"):
        classify_sources(src, strict=False)


def test_wrong_quadratic_slot_count_refused():
    src = sample_sources()
    src["quadratic"]["models"]["B_-1"]["independent_complex_mass_parameters_in_declared_basis"] = 0
    with pytest.raises(AuditError, match="Quadratic counts differ"):
        classify_sources(src, strict=False)


def test_report_formats_and_provenance(tmp_path: Path):
    report = classify_sources(sample_sources(), strict=False)
    hashes = {"docs/example.json": "a" * 64}
    write_reports(report, hashes, tmp_path / "reports")
    root = tmp_path / "reports"
    full = json.loads((root / "operator_classification.json").read_text(encoding="utf-8"))
    assert full["input_sha256"] == hashes
    assert (root / "operator_classification.csv").read_text(encoding="utf-8").count("B_-1") > 8
    assert "Input hashes" in (root / "summary.md").read_text(encoding="utf-8")


def test_strict_16_model_requirement():
    with pytest.raises(AuditError, match="Model set mismatch"):
        classify_sources(sample_sources(), strict=True)
