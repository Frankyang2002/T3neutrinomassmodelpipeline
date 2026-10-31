"""Independent unit tests for Workstream B SU2 generator and phase audit."""
from __future__ import annotations

import copy

import numpy as np
import pytest

from scripts import VerifyT3BCSU2Basis as s


def rotated():
    theta=0.37
    p=np.array([[np.cos(theta),-np.sin(theta),0],
                [np.sin(theta), np.cos(theta),0],
                [0,0,1.]])
    return np.einsum("ab,bij->aij",p,s.PAULI)


def test_canonical_generators_are_normalized_su2():
    r=s.generators_check(s.PAULI)
    assert r["max_su2_commutator_residual"] < 1e-13
    assert r["max_generator_reconstruction_residual"] < 1e-13
    assert np.allclose(r["adjoint_rotation_matchete_to_canonical"],np.eye(3))


def test_nontrivial_native_adjoint_rotation():
    r=s.generators_check(rotated())
    assert not np.allclose(r["adjoint_rotation_matchete_to_canonical"],np.eye(3))
    assert r["max_adjoint_orthogonality_residual"] < 1e-12
    assert r["adjoint_rotation_determinant"] == pytest.approx(1)


def test_native_generator_bad_trace_rejected():
    g=s.PAULI.copy()
    g[0]+=0.1*s.IDENTITY
    with pytest.raises(ValueError,match="non-traceless|Non-Hermitian"):
        s.generators_check(g)


def test_native_generator_wrong_factor_two_rejected():
    with pytest.raises(ValueError,match="normalized"):
        s.generators_check(2*s.PAULI)


def test_native_generator_wrong_phase_hermiticity_rejected():
    with pytest.raises(ValueError,match="Non-Hermitian"):
        s.generators_check(1j*s.PAULI)


def test_component_arrays_dimension_guard():
    with pytest.raises(ValueError,match="dimensions"):
        s.matrix_payload(np.ones((2,2,2)),np.zeros((2,2,2)))


def test_component_arrays_nonfinite_guard():
    a=np.zeros((3,2,2))
    a[1,1,1]=float("nan")
    with pytest.raises(ValueError,match="Nonfinite"):
        s.matrix_payload(a,np.zeros((3,2,2)))


def test_registered_cg_invariance_and_random_su2_transformations():
    r=s.invariant_check(rotated(),s.EPS)
    assert r["max_finite_su2_covariance_residual"] < 1e-12
    assert r["max_gauge_rotated_bilinear_residual"] < 1e-12
    assert r["max_rephased_cross_tool_action_residual"] < 1e-12
    assert r["unrephased_coupling_negative_control"] > 0.01
    assert not r["independent_scalar_u1_phase_fixed_by_generators"]


def test_native_cg_epsilon_wrong_normalization_detected():
    c=s.EPS.copy()
    c[1,0]*=2
    with pytest.raises(ValueError,match="epsilon"):
        s.invariant_check(s.PAULI,c)


def test_native_cg_wrong_phase_is_conditional_not_physical():
    c=1j*s.EPS
    r=s.invariant_check(s.PAULI,c)
    assert r["canonical_ratio_if_same_component_phase"] == pytest.approx([0,-1])
    assert not r["independent_scalar_u1_phase_fixed_by_generators"]


def test_negative_control_fails_under_independent_rephase():
    r=s.invariant_check(s.PAULI,s.EPS)
    assert r["unrephased_coupling_negative_control"] > 0.01


def test_reversed_generator_order_detected_as_adjoint_orientation():
    g=s.PAULI.copy()
    g[[0,1]]=g[[1,0]]
    r=s.generators_check(g)
    assert r["adjoint_rotation_determinant"] == pytest.approx(-1)
    assert r["max_su2_commutator_residual"] < 1e-12


def test_native_report_success_and_uncertified_flags_required():
    row={"model":"B_-1","status":"LiveMatcheteSU2GeneratorsCaptured"}
    report={"format":s.FORMAT,"physical_mapping_certified":False,
        "relative_field_phase_certified":False,
        "rgbeta_numeric_generators_extracted":False,
        "models":{"B_-1":row,"C_-1":{"model":"C_-1"}}}
    assert set(s.check_generator_document(report)) == set(s.MODELS)
    bad=copy.deepcopy(report)
    bad["relative_field_phase_certified"]=True
    with pytest.raises(ValueError,match="format or flags"):
        s.check_generator_document(bad)


def test_malformed_generator_report_rejected():
    with pytest.raises(ValueError,match="Missing native B/C"):
        s.check_generator_document({"format":s.FORMAT,
            "relative_field_phase_certified":False,
            "physical_mapping_certified":False,
            "rgbeta_numeric_generators_extracted":False,
            "models":{"B_-1":{}}})


def test_full_synthetic_audit_requires_prior_certificate_and_preserves_phase_flag(monkeypatch):
    from scripts import VerifyT3BCMixedBilinearMap as prior
    monkeypatch.setattr(prior, "check_source", lambda _: {
        "Lagrangian/interactions/GeneralScalarOperators.wl": "sha_m",
        "RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl": "sha_r",
        "RGE/running/rgbeta/T3RGBetaModel.wl": "sha_g",
    })
    mock_mat = {"format":prior.MAT_FORMAT, "physical_mapping_certified":False,
        "source_sha256":{"Lagrangian/interactions/GeneralScalarOperators.wl":"sha_m"},"models":{}}
    mock_rgb = {"format":prior.RGB_FORMAT,"physical_mapping_certified":False,
        "source_sha256":{
            "RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl":"sha_r",
            "RGE/running/rgbeta/T3RGBetaModel.wl":"sha_g"},"models":{}}
    for model in s.MODELS:
        mock_mat["models"][model] = {
            "model":model,"status":"LiveMatcheteMixedBilinearCGCaptured",
            "counts":[0,0,0,1,0,1],"matchete_coupling":prior.COUPLING_M,
            "matchete_cg":"T3ExtraCGM000101I1",
            "tensor_accessor":"Matchete`GetCGTensor",
            "cg_present_in_actual_plus_hc_expression":True,
            "cg_components_inputform":[["0","1"],["-1","0"]],
            "physical_mapping_certified":False}
        mock_rgb["models"][model] = {
            "model":model,"status":"LiveRGBetaMixedBilinearRegisteredAndLamDetected",
            "rgbeta_mass_symbol":prior.COUPLING_R,
            "rgbeta_fields":["Bar[S1]","Bar[S2]"],
            "rgbeta_self_conjugate":False,
            "rgbeta_unique_permutations":12,
            "native_lam_symbol":"RGBeta`PackageScope`Lam",
            "native_massive_lam_evaluated":True,
            "native_lam_mass_sparse_array_rules_count":24,
            "native_rgbeta_epsilon_component_phase_verified":False,
            "rgbeta_mass_coupling_function_inputform":"12*muT3GeneralS1S2",
            "rgbeta_invariant_inputform":"eps[SU2L[fund],i,j]",
            "rgbeta_potential_convention":"+muR * eps * Bar[S1] Bar[S2] + h.c.; - sign in Lagrangian",
            "physical_mapping_certified":False}
    old=prior.audit(mock_mat,mock_rgb)
    generator={"format":s.FORMAT,"relative_field_phase_certified":False,
        "physical_mapping_certified":False,"rgbeta_numeric_generators_extracted":False,
        "models":{}}
    for model in s.MODELS:
        generator["models"][model]={
            "model":model,"status":"LiveMatcheteSU2GeneratorsCaptured",
            "accessor":"Matchete`Generators[SU[2], {1}]",
            "generator_dimensions":[3,2,2],
            "generator_real":s.PAULI.real.tolist(),
            "generator_imag":s.PAULI.imag.tolist(),
            "relative_field_phase_certified":False,
            "physical_mapping_certified":False}
    report=s.audit(generator,mock_mat,mock_rgb,old)
    assert report["physical_mapping_certified"] is False
    assert len(report["models"]) == 2
    assert all(row["conditional_mu_rgbeta_over_mu_matchete"] == "-1"
               for row in report["models"])
    bad=copy.deepcopy(generator)
    bad["models"]["C_-1"]["status"]="Blocked"
    with pytest.raises(ValueError,match="not successful"):
        s.audit(bad,mock_mat,mock_rgb,old)
