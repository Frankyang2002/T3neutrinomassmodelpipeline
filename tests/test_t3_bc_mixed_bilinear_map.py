"""Independent algebra and fail-closed synthetic-input tests, no Wolfram dependency."""
from __future__ import annotations

from copy import deepcopy
import numpy as np
import pytest
import sympy as sp

from scripts import VerifyT3BCMixedBilinearMap as audit


def synthetic_live(k: str = "1/Sqrt[2]"):
    mat={"format":audit.MAT_FORMAT,"physical_mapping_certified":False,"models":{},
         "source_sha256":{"Lagrangian/interactions/GeneralScalarOperators.wl":"sha_m"}}
    rgb={"format":audit.RGB_FORMAT,"physical_mapping_certified":False,"models":{},
         "source_sha256":{
            "RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl":"sha_r",
            "RGE/running/rgbeta/T3RGBetaModel.wl":"sha_g"}}
    for model in audit.MODELS:
        mat["models"][model]={
            "model":model,"status":"LiveMatcheteMixedBilinearCGCaptured",
            "counts":[0,0,0,1,0,1],"matchete_coupling":audit.COUPLING_M,
            "matchete_cg":"T3ExtraCGM000101I1",
            "tensor_accessor":"Matchete`GetCGTensor",
            "cg_present_in_actual_plus_hc_expression":True,
            "cg_components_inputform":[["0",k],["-("+k+")","0"]],
            "physical_mapping_certified":False}
        rgb["models"][model]={
            "model":model,"status":"LiveRGBetaMixedBilinearRegisteredAndLamDetected",
            "rgbeta_mass_symbol":audit.COUPLING_R,
            "rgbeta_fields":["Bar[S1]","Bar[S2]"],
            "rgbeta_self_conjugate":False,
            "rgbeta_unique_permutations":12,
            "native_lam_symbol":"RGBeta`PackageScope`Lam",
            "native_massive_lam_evaluated":True,
            "native_lam_mass_sparse_array_rules_count":12,
            "native_rgbeta_epsilon_component_phase_verified":False,
            "rgbeta_mass_coupling_function_inputform":"12*muT3GeneralS1S2",
            "rgbeta_invariant_inputform":"eps[SU2L[fund],i,j]",
            "rgbeta_potential_convention":"+muR * eps * Bar[S1] Bar[S2] + h.c.; - sign in Lagrangian",
            "physical_mapping_certified":False}
    return mat,rgb


def test_matchete_cg_factor_is_not_assumed():
    for k in ("1/Sqrt[2]", "-Sqrt[2]/3", "I/2"):
        mat,_=synthetic_live(k)
        C=audit.check_matchete(mat,"B_-1")
        factor=sp.simplify(C[0,1])
        assert C == factor*audit.EPS
        diagnostics=audit.check_equivalence(C,factor)
        assert diagnostics["max_action_relative_residual"]<1e-12
        assert diagnostics["wrong_complex_conjugation_negative_control"]>1e-3


def test_reject_non_epsilon_tensor():
    mat,_=synthetic_live()
    mat["models"]["B_-1"]["cg_components_inputform"][1][0]="1"
    C=audit.check_matchete(mat,"B_-1")
    k=C[0,1]
    assert (C-k*audit.EPS)!=sp.zeros(2)


def test_rgbeta_factor_and_field_order_fail_closed():
    _,rgb=synthetic_live()
    assert audit.check_rgbeta(rgb,"C_-1")["internal_real_tensor_factor"]==12
    wrong=deepcopy(rgb)
    wrong["models"]["C_-1"]["rgbeta_mass_coupling_function_inputform"]="24*muT3GeneralS1S2"
    with pytest.raises(ValueError,match="normalization"):
        audit.check_rgbeta(wrong,"C_-1")
    wrong=deepcopy(rgb)
    wrong["models"]["C_-1"]["rgbeta_fields"]=["Bar[S2]","Bar[S1]"]
    with pytest.raises(ValueError,match="contract"):
        audit.check_rgbeta(wrong,"C_-1")


@pytest.mark.parametrize("form", (
    "12*muT3GeneralS1S2",
    "12*Global`muT3GeneralS1S2",
))
def test_evaluated_rgbeta_scalar_mass_factor_twelve(form):
    _, rgb = synthetic_live()
    rgb["models"]["B_-1"]["rgbeta_mass_coupling_function_inputform"] = form
    result = audit.check_rgbeta(rgb,"B_-1")
    assert result["internal_real_tensor_factor"] == 12
    assert result["evaluated_mass_coupling_function_inputform"] == form


@pytest.mark.parametrize("form", (
    "12*Matrix[muT3GeneralS1S2][]",  # unevaluated: RGBeta collapses Matrix[m_][]
    "12*Matrix[muT3GeneralS1S2]",
    "6*muT3GeneralS1S2", "24*muT3GeneralS1S2",
    "-12*muT3GeneralS1S2", "12*Bar[muT3GeneralS1S2]",
    "12*muT3GeneralS1S2+0*otherSymbol",
    "12*muT3GeneralS2S1",
))
def test_reject_non_native_factor_or_symbol(form):
    _, rgb = synthetic_live()
    rgb["models"]["B_-1"]["rgbeta_mass_coupling_function_inputform"] = form
    with pytest.raises(ValueError,match="normalization/symbol"):
        audit.check_rgbeta(rgb,"B_-1")


def test_no_runtime_evidence_never_accepted():
    mat,rgb=synthetic_live()
    mat["models"]["B_-1"]["status"]="Blocked"
    with pytest.raises(ValueError):audit.check_matchete(mat,"B_-1")
    _,rgb=synthetic_live()
    rgb["models"]["B_-1"]["native_massive_lam_evaluated"]=False
    with pytest.raises(ValueError):audit.check_rgbeta(rgb,"B_-1")


def test_wrong_conjugation_distinguishable():
    C=sp.Matrix([[0,sp.sqrt(2)/2],[-sp.sqrt(2)/2,0]])
    report=audit.check_equivalence(C,sp.sqrt(2)/2)
    assert report["max_action_relative_residual"]<1e-12
    assert report["max_real_hessian_difference"]<1e-8
    assert report["max_mass_eigenvalue_residual"]<1e-8
    assert report["wrong_complex_conjugation_negative_control"]>1e-3
    assert report["wrong_factor_two_negative_control"]>1e-3


def test_component_parser_rejects_unknown_symbols():
    assert audit.parse_component("Sqrt[2]/2")==sp.sqrt(2)/2
    for inp in ("xyz", "Sqrt[2] + f[3]", "1/0"):
        with pytest.raises((ValueError,TypeError,ZeroDivisionError,sp.SympifyError)):
            audit.parse_component(inp)


def test_full_audit_returns_uncertified(monkeypatch):
    monkeypatch.setattr(audit,"check_source",lambda root:{
        "Lagrangian/interactions/GeneralScalarOperators.wl":"sha_m",
        "RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl":"sha_r",
        "RGE/running/rgbeta/T3RGBetaModel.wl":"sha_g"})
    mat,rgb=synthetic_live()
    result=audit.audit(mat,rgb)
    assert result["physical_mapping_certified"] is False
    assert len(result["models"])==2
    for model in result["models"]:
        assert model["conditional_mu_rgbeta_over_mu_matchete_exact"] in ("-sqrt(2)/2","-1/sqrt(2)")
        assert model["requires_conjugation_under_same_holomorphic_field_word"] is False


def test_stale_probe_rejected(monkeypatch):
    monkeypatch.setattr(audit,"check_source",lambda root:{
        "Lagrangian/interactions/GeneralScalarOperators.wl":"new_sha",
        "RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl":"sha_r",
        "RGE/running/rgbeta/T3RGBetaModel.wl":"sha_g"})
    mat,rgb=synthetic_live()
    with pytest.raises(ValueError,match="source changed"):
        audit.audit(mat,rgb)


def test_su2_bilinear_covariance():
    rng=np.random.default_rng(613)
    eps=audit.EPS_NUM
    for i in range(20):
        a,b=rng.normal(size=2)+1j*rng.normal(size=2)
        norm=np.sqrt(abs(a)**2+abs(b)**2)
        a,b=a/norm,b/norm
        g=np.array([[a,b],[-b.conjugate(),a.conjugate()]])
        assert np.allclose(g.conjugate().T@eps@g.conjugate(),eps,rtol=0,atol=5e-15)
        x=rng.normal(size=2)+1j*rng.normal(size=2)
        y=rng.normal(size=2)+1j*rng.normal(size=2)
        assert np.allclose(audit.bilinear(g@x,g@y,eps),
                           audit.bilinear(x,y,eps),rtol=0,atol=5e-14)
