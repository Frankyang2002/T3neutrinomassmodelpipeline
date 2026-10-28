"""Offline field-flow and source-artifact consistency tests for Batch 2."""
from __future__ import annotations

import copy
import csv
import hashlib
import json
from pathlib import Path

import pytest

from scripts.AuditT3OneLoopScalarVertices import (
    AuditError,
    analyze_batch1,
    canonical_vertex_test,
    conjugate_counts,
    validate_batch1,
    verify_source_hashes,
    write_result,
)

PROJECT = Path(__file__).resolve().parents[2]
REPORT = PROJECT / "output/scalar_relevance/operator_classification.json"


@pytest.fixture
def source_report():
    if not REPORT.is_file():
        pytest.skip("Run Batch 1 first to create output/scalar_relevance/operator_classification.json")
    return json.loads(REPORT.read_text(encoding="utf-8"))


def test_hermitian_conjugation_is_involution():
    value = (2, 0, 1, 0, 0, 1)
    assert conjugate_counts(conjugate_counts(value)) == value
    assert conjugate_counts(value) == (0, 2, 0, 1, 1, 0)


def test_original_T3_y1y2_vertex_survives():
    result = canonical_vertex_test((2, 0, 1, 0, 0, 1))
    assert result["status"] == "single_vertex_wick_eligible"
    assert result["orientation"] == "hermitian_conjugate"
    assert result["channel"] == "y1_y2"


def test_BC_two_y1_scalar_vertex_survives():
    result = canonical_vertex_test((0, 2, 0, 2, 0, 0))
    assert result["status"] == "single_vertex_wick_eligible"
    assert result["channel"] == "y1_y1"


def test_BC_two_y2_scalar_vertex_survives():
    result = canonical_vertex_test((0, 2, 0, 0, 2, 0))
    assert result["status"] == "single_vertex_wick_eligible"
    assert result["channel"] == "y2_y2"


def test_mixed_H_Hdag_is_not_holomorphic_Weinberg_vertex():
    result = canonical_vertex_test((1, 1, 0, 1, 0, 1))
    assert result["status"] == "higgs_chirality_mismatch"


@pytest.mark.parametrize("counts", [
    (0, 2, 1, 0, 1, 0),    # Hdag^2 S1 S2
    (0, 2, 0, 1, 0, 1),    # Hdag^2 S1dag S2dag
    (0, 2, 2, 0, 0, 0),    # Hdag^2 S1^2
    (0, 2, 0, 0, 0, 2),    # Hdag^2 S2dag^2
])
def test_scalar_orientation_exclusions(counts):
    assert canonical_vertex_test(counts)["status"] == "scalar_wick_orientation_mismatch"


def test_wrong_count_rejected():
    with pytest.raises(AuditError):
        canonical_vertex_test((0, 2, 0, 1, 0, 0))


def test_full_16_model_report_matches_exact_saved_C5(source_report):
    result = analyze_batch1(source_report, source_hashes_verified=False)
    counts = result["counts"]
    assert (counts["models_checked"], counts["new_two_H_vertices"],
            counts["new_direct_wick_eligible"], counts["new_no_H_vertex"]) == (16, 24, 4, 26)
    assert counts["excluded_H_Hdag_by_higgs_orientation"] == 10
    assert counts["excluded_Hdag2_by_scalar_orientation"] == 10
    assert {o["model"] for o in result["operators"]
            if o["status"] == "single_vertex_wick_eligible_and_saved_C5_observed"} == {"B_-1", "C_-1"}
    assert len(result["operators"]) == 186


def test_nonzero_quadratic_mass_is_never_excluded(source_report):
    result = analyze_batch1(source_report, source_hashes_verified=False)
    masses = [o for o in result["operators"] if o["kind"] == "additional_quadratic"]
    assert len(masses) == 8
    assert all(o["status"] == "unsuppressed_tree_level_mass_effect" for o in masses)


def test_falsified_C5_claim_rejected(source_report):
    changed = copy.deepcopy(source_report)
    entry = next(o for o in changed["operators"] if o["kind"] == "additional_quartic"
                 and o["model"] == "B_-1" and o["coupling"].endswith("Q110101I1"))
    entry["observed_in_saved_C5"] = True
    with pytest.raises(AuditError, match="contradicts Wick"):
        analyze_batch1(changed, source_hashes_verified=False)


def test_bad_coupling_count_rejected(source_report):
    changed = copy.deepcopy(source_report)
    entry = next(o for o in changed["operators"] if o["kind"] == "additional_quartic")
    entry["counts_H_Hdag_S1_S1dag_S2_S2dag"] = [0, 2, 0, 2, 0, 0]
    with pytest.raises(AuditError, match="coupling/field-count mismatch"):
        validate_batch1(changed)


def test_source_hash_verification_fails_closed_on_absent_evidence(source_report, tmp_path):
    with pytest.raises(AuditError, match="missing Batch 1 evidence file"):
        verify_source_hashes(source_report, tmp_path)


def test_output_files_and_summary(source_report, tmp_path):
    result = analyze_batch1(source_report, source_hashes_verified=False)
    result["batch1_report_sha256"] = hashlib.sha256(REPORT.read_bytes()).hexdigest()
    paths = write_result(result, tmp_path)
    assert len(paths) == 3
    from_disk = json.loads((tmp_path / "leading_c5_vertex_audit.json").read_text())
    assert from_disk["counts"]["new_direct_wick_eligible"] == 4
    with (tmp_path / "leading_c5_vertex_audit.csv").open(newline="", encoding="utf-8") as f:
        assert len(list(csv.DictReader(f))) == 186
    assert "DIAGONAL" in (tmp_path / "leading_c5_vertex_summary.md").read_text()
