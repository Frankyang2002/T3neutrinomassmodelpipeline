"""Source-independent synthetic end-to-end batch and coupling-ledger regression."""
from __future__ import annotations

import json
from pathlib import Path
from fractions import Fraction

from scripts import RunT3GeneralBatchVerification as batch
from scripts.T3GeneralGroupFactorProjection import (
    expected_gauge_betas, expected_scalar_gauge_factors,
    expected_yukawa_gauge_factors,
)


def _fixtures(folder: Path, corrupt: bool = False):
    model, dims, alpha = "A_-2", (1, 3, 2), -2
    entry = {
        "new_rgbeta_symbols": [{"name": "lambdaT3GeneralQ000004I1",
                               "counts": [0, 0, 0, 0, 0, 4]}],
        "preexisting_rgbeta_symbols": [],
    }
    for stage in ("UV", "EFT1"):
        betas = {g: f"({b})*{g}^3" for g, b in
                 expected_gauge_betas(dims, alpha, stage).items()}
        if corrupt and stage == "UV":
            betas["g3"] = "-6*g3^3"
        for item in batch.scalar_interactions(model, entry, "self_s2"):
            factors = expected_scalar_gauge_factors(item["counts"], dims, alpha)
            terms = [f"({x})*{g}^2*{item['name']}" for g, x in factors.items() if x]
            betas[item["name"]] = "+".join(terms) if terms else "lambdaH^2"
        if stage == "UV":
            fac = expected_yukawa_gauge_factors(dims, alpha)
            for key in ("y1", "y2", "y2Conj"):
                betas[key] = "+".join(f"({x})*{g}^2*Matrix[{key}][$i,$j]"
                                       for g, x in fac.items())
        (folder / f"{stage.lower()}_rgbeta.json").write_text(
            json.dumps({"report_betas": betas}), encoding="utf-8")
    return entry


def test_source_independent_a2_native_beta_projection(tmp_path):
    entry = _fixtures(tmp_path)
    checks, scalars, terms = batch.project_model("A_-2", entry, "self_s2", tmp_path)
    assert len(checks) > 50
    assert all(r["status"] == batch.STATUS_PASS for r in checks)
    assert {x["stage"] for x in terms} == {"UV", "EFT1"}
    assert {x["sector"] for x in scalars} >= {"baseline_quartic", "additional_bilinear"}


def test_wrong_exact_gauge_group_factor_is_detected(tmp_path):
    entry = _fixtures(tmp_path, corrupt=True)
    checks, _, _ = batch.project_model("A_-2", entry, "self_s2", tmp_path)
    mismatch = [r for r in checks if r["status"] == batch.STATUS_MISMATCH]
    assert len(mismatch) == 1
    assert mismatch[0]["stage"] == "UV" and mismatch[0]["beta"] == "g3"
    assert mismatch[0]["residual"] == "1"


def test_multicontraction_quartic_dictionary_is_not_falsely_identified():
    model = "B_-1"
    entry = {
        "new_rgbeta_symbols": [
            {"name": "lambdaT3GeneralQ110101I1", "counts": [1, 1, 0, 1, 0, 1]}],
        "preexisting_rgbeta_symbols": [
            {"existing_name": "lambdaHHdagS1barS2barCross", "counts": [1, 1, 0, 1, 0, 1]}],
    }
    manifest = {"models": {model: {
        "extra_quartic_structures": [
            {"counts": [1, 1, 0, 1, 0, 1], "invariants": 2}],
        "quadratic_structures_not_enabled": [
            {"counts": [0, 0, 0, 1, 0, 1]}],
    }}}
    rows = batch.coupling_ledger(model, entry, manifest, "mixed")
    quartics = [r for r in rows if r["sector"] == "extra_quartic"]
    assert len(quartics) == 2
    assert all(r["correspondence"] == "BasisBlockUndetermined" for r in quartics)
    assert all(r["normalization"] == batch.STATUS_UNCERTIFIED for r in quartics)


def test_report_model_accepts_normalized_structural_status(tmp_path):
    batch.report_model(
        "A_-2", tmp_path, [], [], [],
        {"model": "A_-2", "structural_status": "VerifiedStructural"},
        tmp_path, "full16_general_symbolic",
    )
    report = (tmp_path / "models" / "A_m2.md").read_text(encoding="utf-8")
    assert "Source-matched structural inventory:** VerifiedStructural" in report
