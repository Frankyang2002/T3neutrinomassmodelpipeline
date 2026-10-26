"""Independent SU(2) tensor check for the C5-active B/C(-1) quartics.

RGBeta source: aethomsen/RGBeta Package/GroupsAndIndices.m, SU(n)
branch, where delS2 delS2 is the projector onto Sym^2(fund). The
corresponding source-model contractions are in
RGE/running/rgbeta/T3RGBetaModel.wl.

This checks the tensor *shape and factor*. It is NOT a proof that Matchete's
and RGBeta's complex Lagrangian coefficients have the same overall sign/phase
or Hermitian-conjugation convention.
"""
from __future__ import annotations

from pathlib import Path
import json
import re
from typing import Any

import numpy as np

MAPPING = {
    "lambdaT3ExtraQ020020I1": "lambdaHHdagS2S2",
    "lambdaT3ExtraQ020200I1": "lambdaHHdagS1barS1bar",
}
MODELS = ("B_m1", "C_m1")


def symmetric_doublet_projector() -> np.ndarray:
    """P_ab,cd = (delta_ac delta_bd + delta_ad delta_bc)/2."""
    ident = np.eye(2, dtype=int)
    return (np.einsum("ac,bd->abcd", ident, ident) +
            np.einsum("ad,bc->abcd", ident, ident)) / 2.0


def epsilon_doublet() -> np.ndarray:
    """Choose epsilon_12=+1 as in the canonical Matchete tensor audit."""
    return np.array([[0, 1], [-1, 0]], dtype=int)


def all_bar_tensor_from_projector() -> np.ndarray:
    """Apply SU(2) pseudo-reality on two indices to obtain eps-eps."""
    projector = symmetric_doublet_projector()
    eps = epsilon_doublet()
    return np.einsum("abuv,uc,vd->abcd", projector, eps, eps)


def _proportionality(source: np.ndarray, target: np.ndarray) -> tuple[float, float]:
    den = float(np.vdot(target, target).real)
    if den <= 0:
        raise ValueError("Cannot compare to zero tensor")
    factor = float(np.vdot(target, source).real / den)
    residual = float(np.max(np.abs(source - factor * target)))
    return factor, residual


def _matchete_entries(document: dict[str, Any]) -> dict[tuple[str, str], dict]:
    if document.get("status") != "MatcheteTensorConventionVerified":
        raise ValueError("Matchete tensor validation was not successful")
    entries = {}
    for block in document.get("results", []):
        model = str(block.get("model", ""))
        for tensor in block.get("tensors", []):
            key = (model, tensor.get("coupling"))
            if key in entries:
                raise ValueError("Duplicate Matchete tensor result")
            if tensor.get("status") != "CanonicalTensorProportional":
                raise ValueError(f"Invalid Matchete tensor result: {key}")
            if str(tensor.get("matchete_over_canonical_factor_inputform")) != "1":
                raise ValueError(f"Nonunit Matchete tensor ratio: {key}")
            if tensor.get("all_component_residuals_zero") is not True:
                raise ValueError(f"Matchete tensor residual was not zero: {key}")
            entries[key] = tensor
    if set(entries) != {(m, c) for m in MODELS for c in MAPPING}:
        raise ValueError("Expected exactly the four certified Matchete tensors")
    return entries


def _verify_project_rgbeta_source(path: Path) -> None:
    """Fail closed if project RGBeta coupling definitions have changed."""
    source = Path(path).read_text(encoding="utf-8-sig")
    for coupling, fields in (
        ("lambdaHHdagS2S2", r"\{Bar\s*@\s*H\s*,\s*Bar\s*@\s*H\s*,\s*S2\s*,\s*S2\}"),
        ("lambdaHHdagS1barS1bar", r"\{Bar\s*@\s*H\s*,\s*Bar\s*@\s*H\s*,\s*Bar\s*@\s*S1\s*,\s*Bar\s*@\s*S1\}"),
    ):
        pattern = re.compile(r"AddQuartic\s*\[\s*"+coupling+r"\s*,\s*"+fields+r"\s*,\s*GroupInvariant\s*->\s*\(\s*delS2\s*\[", re.S)
        if not pattern.search(source):
            raise ValueError("Current RGBeta source no longer has the certified delS2 structure for " + coupling)


def build_tensor_mapping(matchete_report: Path, rgbeta_model_source: Path) -> dict[str, Any]:
    """Report tensor-level factors and explicit remaining coupling ambiguity."""
    obj = json.loads(Path(matchete_report).read_text(encoding="utf-8-sig"))
    _matchete_entries(obj)
    _verify_project_rgbeta_source(rgbeta_model_source)
    proj = symmetric_doublet_projector()
    eps = epsilon_doublet()
    ref_eps = (np.einsum("ac,bd->abcd", eps, eps) +
               np.einsum("ad,bc->abcd", eps, eps)) / 2.0
    transformed = all_bar_tensor_from_projector()
    delta_factor, delta_residual = _proportionality(proj, proj)
    eps_factor, eps_residual = _proportionality(transformed, ref_eps)
    if delta_residual > 1e-14 or eps_residual > 1e-14:
        raise ValueError("RGBeta SU(2) projector does not match canonical tensors")
    data = []
    for model in MODELS:
        for mat, rg in MAPPING.items():
            factor = delta_factor if mat.endswith("020020I1") else eps_factor
            data.append({
                "model": model, "matchete_coupling": mat,
                "rgbeta_coupling": rg,
                "rgbeta_over_canonical_tensor_factor": factor,
                "matchete_over_canonical_tensor_factor": 1,
                "rgbeta_over_matchete_tensor_factor": factor,
                "relative_phase_convention": "epsilon_12 = +1; two epsilon conversions applied for the all-bar structure",
                "full_lagrangian_coupling_identification_verified": False,
            })
    return {
        "status": "DoubletGroupTensorRatiosVerified",
        "scope": "ordinary split-scalar B/C(-1) only",
        "source_rgbeta_project_model": str(rgbeta_model_source),
        "matchete_reference": str(matchete_report),
        "rgbeta_projector_identity": "delS2(a,i,j) delS2(a,k,l) = (delta_ik delta_jl+delta_il delta_jk)/2, from aethomsen/RGBeta Package/GroupsAndIndices.m",
        "delta_projector_residual": delta_residual,
        "epsilon_projector_residual": eps_residual,
        "tensor_ratio_only": True,
        "full_lagrangian_coupling_identification_verified": False,
        "caveat": "The RGBeta index contraction factor is verified from its documented/source SU(2) projector. Relative phases or signs and scalar potential versus interaction Lagrangian conventions still need a vertex-level check before inserting RGBeta values into Matchete C5.",
        "mappings": data,
    }
