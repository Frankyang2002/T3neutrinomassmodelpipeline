"""Fast offline algebra/CLI tests. No Matchete and no expensive matching."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from Numerical.running.T3SequentialBatchNumerics import (
    DEFAULT_MODELS, MODELS, PhysicsInputError, assemble_llss,
    dminus2_existing_reference, from_archived_v6,
    from_dynamic_case, mass_class, model_tag, normalize_model,
)


def encoded(a):
    x = np.asarray(a, dtype=float)
    return {"real": x.tolist(), "imag": np.zeros_like(x).tolist()}


def baseline(*, model="D_+2"):
    log = np.diag([-9e-16, -8e-16, -7e-16])
    bubble = np.diag([1e-16, 1e-16, 1e-16])
    return model, encoded(log), encoded(bubble), encoded(log + bubble)


def test_catalog_all16_three_mass_classes():
    assert len(MODELS) == 16
    assert set(DEFAULT_MODELS) == {"D_+2", "A_-4", "B_-1"}
    assert sum(mass_class(m) == "ordinary" for m in MODELS) == 8
    assert sum(mass_class(m) == "self_bilinear" for m in MODELS) == 6
    assert sum(mass_class(m) == "mixed" for m in MODELS) == 2
    assert normalize_model("D_p2") == "D_+2"
    assert model_tag("B_-1") == "B_m1"
    with pytest.raises(PhysicsInputError):
        normalize_model("C_+2")


def test_llss_numerical_closure_takagi_and_vev():
    m, running, bubble, total = baseline()
    result = assemble_llss(model=m, leading_raw=running, bubble_raw=bubble,
                           total_raw=total, projection_residual=0,
                           scale_residual=1e-16)
    expected = np.sort(np.array([8, 7, 6]) * 1e-16 * 246.22**2 / 2 * 1e9)
    assert np.allclose(result["masses_eV_per_abs_unknown_NGamma"], expected)
    assert result["C5_relative_closure_residual"] < 1e-14
    assert result["complete_C5"] is False
    assert result["hard_one_loop_matching_independently_computed"] is False


def test_invalid_closure_and_symmetry_are_blocked():
    m, running, bubble, total = baseline()
    bad = encoded(np.diag([-9e-16, -8e-16, -7e-16]))
    with pytest.raises(PhysicsInputError, match="closure"):
        assemble_llss(model=m, leading_raw=running, bubble_raw=bubble,
                      total_raw=bad, projection_residual=0, scale_residual=0)
    asymmetric = np.diag([-8e-16, -7e-16, -6e-16])
    asymmetric[0, 1] += 1e-16
    with pytest.raises(PhysicsInputError):
        assemble_llss(model=m, leading_raw=running, bubble_raw=bubble,
                      total_raw=encoded(asymmetric), projection_residual=0,
                      scale_residual=0)


def test_archived_fixture_and_provenance(tmp_path):
    m, running, bubble, total = baseline()
    path = tmp_path / "D_p2.json"
    masses = sorted((np.linalg.svd(np.array(total["real"]),
                                   compute_uv=False) * 246.22**2 / 2 * 1e9).tolist())
    fixture = {"model": m,
               "format": "t3-conditional-c5-per-unknown-gammacc-mixed-split-v6",
               "complete_one_loop_C5_certified": False,
               "physical_C5_certified": False,
               "response": {
                   "status": "ConditionalOneLoopLLSSResponse_PerUnitUncertifiedSpinorNormalization",
                   "model_source_boundary_sha256": "test-source-hash",
                   "mu_f_gev": 1e5, "mu_s_gev": 1e3,
                   "C5_EFT1_induced_leading_log_per_N_gev_inverse": running,
                   "C5_joint_scalar_finite_bubble_per_N_gev_inverse": bubble,
                   "C5_LLSS_total_one_loop_per_N_gev_inverse": total,
                   "weinberg_beta_projection_residual": 0,
                   "weinberg_bubble_projection_residual": 0,
                   "relative_matching_scale_cancellation_residual": 0},
               "numerical_metrics": {
                   "conditional_neutrino_mass_singular_values_eV_per_abs_N": masses}}
    path.write_text(json.dumps(fixture), encoding="utf-8")
    row = from_archived_v6(path, expected_model="D_+2",
                           vev_gev=246.22, normalization_factor=1)
    assert row["model"] == m
    assert row["source_json_sha256"]
    fixture["model"] = "B_-1"
    path.write_text(json.dumps(fixture), encoding="utf-8")
    with pytest.raises(PhysicsInputError):
        from_archived_v6(path, expected_model="D_+2", vev_gev=246.22,
                         normalization_factor=1)


def test_dynamic_fixture_uses_fixed_order_case_110():
    m, running, bubble, total = baseline()
    case = {"c5_eft1_leading_log": running,
            "c5_scalar_finite_bubble": bubble,
            "c5_at_scalar_threshold": total,
            "weinberg_projection_residual": 0,
            "matching_scale_cancellation_residual": 0}
    raw = {"model": m, "status": "ExploratoryGeneralLLSSComputed",
           "cases": {"110": case}, "source_boundary_sha256": "abcd"}
    row = from_dynamic_case(raw, model=m, vev_gev=246.22,
                            normalization_factor=1, with_smeft=False)
    assert row["source_boundary_sha256"] == "abcd"
    assert row["masses_eV_per_abs_unknown_NGamma"][0] > 0
    with pytest.raises(PhysicsInputError, match="SMEFT-on"):
        from_dynamic_case(raw, model=m, vev_gev=246.22,
                          normalization_factor=1, with_smeft=True)


def test_dminus2_manifest_read_only(tmp_path):
    manifest = tmp_path / "manifest.json"
    assert dminus2_existing_reference(manifest)["status"] == "NotPresentOnThisMachine"
    manifest.write_text(json.dumps({
        "model": "D_-2",
        "status": "SequentialC5BenchmarkReproduced_ProvisionalNormalization",
        "scalar_one_loop_matched": True,
        "fixed_order_c5_candidate_assembled": True,
        "nf1_simple_pole_cancellation_verified": True,
        "boundary_hashes_unchanged": True,
        "neutrino_masses_eV": [0, 0.000563839, 0.001937568],
    }))
    assert dminus2_existing_reference(manifest)["status"] == "PassedLocalDMinus2Reference"
