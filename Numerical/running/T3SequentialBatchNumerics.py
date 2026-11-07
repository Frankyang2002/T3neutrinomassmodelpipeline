"""Source-traceable, *conditional* all-16 T3 sequential LLSS numerical handoff.

Physics retained here: UV -> all F -> EFT1(S1,S2), direct LLSS-to-Weinberg
one-loop logarithm, simultaneous all-scalar finite bubble, optional SMEFT RGE.
The independent one-loop hard matching is NOT calculated in this module.

This module deliberately does not infer an unknown Matchete GammaCC/Weyl factor,
convert conditional coefficients into physically certified masses, or supply
missing operators. It validates both freshly computed and archived outputs.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

ALPHAS: dict[str, tuple[int, ...]] = {
    "A": (-4, -2, 0), "B": (-3, -1, 1), "C": (-3, -1, 1),
    "D": (-2, 0, 2), "E": (-4, -2, 0, 2),
}
MODEL_DIMENSIONS = {"A": (1, 3, 2), "B": (2, 2, 1),
                    "C": (2, 2, 3), "D": (3, 1, 2), "E": (3, 3, 2)}
MODELS = tuple(f"{family}_{alpha:+d}" for family, numbers in ALPHAS.items()
               for alpha in numbers)
DEFAULT_MODELS = ("D_+2", "A_-4", "B_-1")
MIXED = frozenset(("B_-1", "C_-1"))
SELF = frozenset(("A_-2", "A_+0", "D_-2", "D_+0", "E_-2", "E_+0"))
ARCHIVED_REL = Path("output/general_observable/lo_joint_scalar_matching/") / (
    "exploratory_all16_v6/C5_per_unknown_NGamma_exploratory")
EXPECTED_ARCHIVED_FORMAT = "t3-conditional-c5-per-unknown-gammacc-mixed-split-v6"
FORMAT = "t3-first-batch-conditional-sequential-llss-v1"


class PhysicsInputError(ValueError):
    """Missing/invalid source evidence; never substitute a synthetic result."""


def normalize_model(model: str) -> str:
    text = model.strip()
    if "_" not in text:
        raise PhysicsInputError(f"Invalid model name: {model}")
    family, suffix = text.split("_", 1)
    if family not in ALPHAS:
        raise PhysicsInputError(f"Invalid model family: {family}")
    if suffix.startswith("m"):
        suffix = "-" + suffix[1:]
    elif suffix.startswith("p"):
        suffix = "+" + suffix[1:]
    try:
        result = f"{family}_{int(suffix):+d}"
    except ValueError as exc:
        raise PhysicsInputError(f"Invalid alpha: {model}") from exc
    if result not in MODELS:
        raise PhysicsInputError(f"Not one of the 16 T3 models: {result}")
    return result


def model_tag(model: str) -> str:
    return normalize_model(model).replace("+", "p").replace("-", "m")


def mass_class(model: str) -> str:
    model = normalize_model(model)
    return "mixed" if model in MIXED else "self_bilinear" if model in SELF else "ordinary"


def read_json(path: Path) -> dict[str, Any]:
    path = Path(path)
    if not path.is_file():
        raise PhysicsInputError(f"Required source missing: {path}")
    try:
        obj = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, ValueError) as exc:
        raise PhysicsInputError(f"Invalid source JSON {path}: {exc}") from exc
    if not isinstance(obj, dict):
        raise PhysicsInputError(f"JSON object expected: {path}")
    return obj


def fingerprint(path: Path) -> str:
    h = sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def matrix(raw: Any, label: str) -> np.ndarray:
    if not isinstance(raw, Mapping) or not {"real", "imag"} <= raw.keys():
        raise PhysicsInputError(f"{label}: missing complex real/imag representation")
    try:
        a = np.asarray(raw["real"], dtype=float)
        b = np.asarray(raw["imag"], dtype=float)
    except (TypeError, ValueError) as exc:
        raise PhysicsInputError(f"{label}: invalid complex entries") from exc
    if a.shape != (3, 3) or b.shape != (3, 3):
        raise PhysicsInputError(f"{label}: C5 must have 3x3 flavour axes")
    c = a + 1j * b
    if not np.isfinite(c).all():
        raise PhysicsInputError(f"{label}: nonfinite C5")
    return c


def encoded(z: np.ndarray) -> dict[str, Any]:
    return {"real": z.real.tolist(), "imag": z.imag.tolist()}


def _real_quantity(raw: Any, label: str) -> float:
    if isinstance(raw, bool) or raw is None:
        raise PhysicsInputError(f"{label}: numerical value required")
    try:
        x = float(raw)
    except (TypeError, ValueError) as exc:
        raise PhysicsInputError(f"{label}: invalid number") from exc
    if not math.isfinite(x):
        raise PhysicsInputError(f"{label}: nonfinite number")
    return x


def assemble_llss(*, model: str, leading_raw: dict, bubble_raw: dict,
                  total_raw: dict, projection_residual: Any,
                  scale_residual: Any, vev_gev: float = 246.22,
                  normalization_factor: float = 1.0) -> dict[str, Any]:
    """Reassemble and check LLSS sector; no independent hard C5 is inferred.

    ``normalization_factor`` is an *explicit choice*. By default the stored
    convention is per unit unknown signed N_GammaCC.
    """
    model = normalize_model(model)
    vev = _real_quantity(vev_gev, "vev_gev")
    factor = _real_quantity(normalization_factor, "normalization_factor")
    if vev <= 0 or factor == 0:
        raise PhysicsInputError("vev must be positive and factor nonzero")
    leading = matrix(leading_raw, "LLSS direct logarithm")
    bubble = matrix(bubble_raw, "LLSS scalar bubble")
    total = matrix(total_raw, "LLSS total")
    den = max(float(np.linalg.norm(total)), 1e-40)
    closure = float(np.linalg.norm(total - leading - bubble) / den)
    symmetry = float(np.linalg.norm(total - total.T) / den)
    if closure > 1e-8 or symmetry > 1e-8:
        raise PhysicsInputError(f"{model}: C5 invalid (closure={closure}, symmetry={symmetry})")
    projection = _real_quantity(projection_residual, "projection_residual")
    scale = _real_quantity(scale_residual, "scale_residual")
    if projection < 0 or scale < 0 or projection > 1e-8 or scale > 1e-7:
        raise PhysicsInputError(f"{model}: Ward/scale check failed: {projection}, {scale}")
    z = total * factor
    # Singular values are Takagi magnitudes for symmetric complex matrices.
    masses = np.sort(np.linalg.svd(z, compute_uv=False) * vev**2 / 2 * 1e9)
    return {
        "model": model, "mass_class": mass_class(model),
        "status": "ConditionalLLSSContributionComputed",
        "C5_LLSS_total_GeV_inverse": encoded(z),
        "C5_LLSS_direct_log_GeV_inverse": encoded(leading * factor),
        "C5_LLSS_joint_scalar_bubble_GeV_inverse": encoded(bubble * factor),
        "C5_LLSS_frobenius_GeV_inverse": float(np.linalg.norm(z)),
        "masses_eV_per_abs_unknown_NGamma": masses.tolist(),
        "vev_GeV": vev,
        "explicit_assumed_spinor_factor": factor,
        "C5_relative_closure_residual": closure,
        "C5_relative_symmetry_residual": symmetry,
        "SU2_Weinberg_projection_residual": projection,
        "matching_scale_cancellation_residual": scale,
        "all_heavy_fermions_integrated_before_scalars": True,
        "both_scalar_multiplets_integrated_together_in_LLSS_bubble": True,
        "hard_one_loop_matching_independently_computed": False,
        "complete_C5": False,
        "model_specific_absolute_SU2_and_spinor_normalization_certified": False,
        "physical_neutrino_masses_certified": False,
    }


def from_archived_v6(path: Path, *, expected_model: str,
                     vev_gev: float, normalization_factor: float) -> dict[str, Any]:
    """Re-evaluate existing source-backed v6 matrices; never rerun matching."""
    obj = read_json(path)
    model = normalize_model(expected_model)
    if obj.get("model") != model or obj.get("format") != EXPECTED_ARCHIVED_FORMAT:
        raise PhysicsInputError(f"{model}: incorrect archived model/format")
    if obj.get("complete_one_loop_C5_certified") is not False or (
            obj.get("physical_C5_certified") is not False):
        raise PhysicsInputError(f"{model}: unexpected archive certification state")
    response = obj.get("response")
    if not isinstance(response, dict) or response.get("model_source_boundary_sha256") is None:
        raise PhysicsInputError(f"{model}: missing source provenance")
    if response.get("status") != "ConditionalOneLoopLLSSResponse_PerUnitUncertifiedSpinorNormalization":
        raise PhysicsInputError(f"{model}: wrong LLSS response status")
    for key in ("mu_f_gev", "mu_s_gev"):
        _real_quantity(response.get(key), key)
    if not float(response["mu_f_gev"]) > float(response["mu_s_gev"]) > 0:
        raise PhysicsInputError(f"{model}: inconsistent matching scale ordering")
    row = assemble_llss(
        model=model,
        leading_raw=response["C5_EFT1_induced_leading_log_per_N_gev_inverse"],
        bubble_raw=response["C5_joint_scalar_finite_bubble_per_N_gev_inverse"],
        total_raw=response["C5_LLSS_total_one_loop_per_N_gev_inverse"],
        projection_residual=max(response["weinberg_beta_projection_residual"],
                                response["weinberg_bubble_projection_residual"]),
        scale_residual=response["relative_matching_scale_cancellation_residual"],
        vev_gev=vev_gev, normalization_factor=normalization_factor,
    )
    row.update({
        "source_mode": "archived_v6_provisional_LLSS_revalidated_no_new_matching",
        "model_source_boundary_sha256": response["model_source_boundary_sha256"],
        "source_json_sha256": fingerprint(path),
        "mu_F_GeV": float(response["mu_f_gev"]),
        "mu_S_GeV": float(response["mu_s_gev"]),
        "scalar_basis_certified": False,
        "matchete_to_RGBeta_quartic_map_certified": False,
    })
    reference = obj.get("numerical_metrics", {}).get(
        "conditional_neutrino_mass_singular_values_eV_per_abs_N")
    if reference is not None and vev_gev == 246.22 and normalization_factor == 1:
        a = np.sort(np.asarray(reference, dtype=float))
        b = np.asarray(row["masses_eV_per_abs_unknown_NGamma"])
        if a.shape != (3,) or not np.allclose(a, b, rtol=1e-10, atol=1e-15):
            raise PhysicsInputError(f"{model}: archived numerical metrics differ from recomputed Takagi magnitudes")
    return row


def from_dynamic_case(result: dict[str, Any], *, model: str,
                      vev_gev: float, normalization_factor: float,
                      with_smeft: bool) -> dict[str, Any]:
    """Consume a new evaluation from RunT3GeneralSequentialLeadingLog.process.

    The EFT1-on, Wilson-on, SMEFT-off case is always ``110``. SMEFT running
    is optional and is reported separately from the C5 scalar threshold.
    """
    if result.get("model") != normalize_model(model) or (
            result.get("status") != "ExploratoryGeneralLLSSComputed"):
        raise PhysicsInputError(f"{model}: invalid source-backed dynamic LLSS result")
    cases = result.get("cases")
    if not isinstance(cases, dict) or "110" not in cases:
        raise PhysicsInputError(f"{model}: missing EFT1-on/Wilson-on case 110")
    case = cases["110"]
    c5 = case.get("c5_at_scalar_threshold")
    row = assemble_llss(
        model=model,
        leading_raw=case["c5_eft1_leading_log"],
        bubble_raw=case["c5_scalar_finite_bubble"],
        total_raw=c5,
        projection_residual=case["weinberg_projection_residual"],
        scale_residual=case["matching_scale_cancellation_residual"],
        vev_gev=vev_gev, normalization_factor=normalization_factor,
    )
    row.update({
        "source_mode": "fresh_python_full_general_LLSS_source_matched_to_reference",
        "source_boundary_sha256": result.get("source_boundary_sha256"),
        "running_state_sha256": result.get("running_state_sha256"),
        "scalar_basis_certified": False,
        "matchete_to_RGBeta_quartic_map_certified": False,
    })
    if with_smeft:
        low = cases.get("111")
        if not isinstance(low, dict) or not low.get("smeft_running_on"):
            raise PhysicsInputError(f"{model}: SMEFT-on case 111 missing")
        low_masses = low.get("masses_ev_ascending")
        if low_masses is None:
            raise PhysicsInputError(f"{model}: low scale masses missing")
        vals = np.asarray(low_masses, dtype=float)
        if vals.shape != (3,) or not np.isfinite(vals).all() or (vals < 0).any():
            raise PhysicsInputError(f"{model}: invalid low-scale masses")
        row["smeft_evolved_masses_eV_per_abs_unknown_NGamma"] = (
            np.sort(vals * abs(normalization_factor)).tolist())
    return row


def dminus2_existing_reference(path: Path) -> dict[str, Any]:
    """Read-only status of the LOCAL completed reference; never changes it."""
    if not path.is_file():
        return {"status": "NotPresentOnThisMachine", "manifest": str(path)}
    doc = read_json(path)
    fields = ("scalar_one_loop_matched", "fixed_order_c5_candidate_assembled",
              "nf1_simple_pole_cancellation_verified", "boundary_hashes_unchanged")
    good = doc.get("model") == "D_-2" and (
        doc.get("status") == "SequentialC5BenchmarkReproduced_ProvisionalNormalization") and (
        all(doc.get(name) is True for name in fields))
    masses = doc.get("neutrino_masses_eV")
    if good:
        good = (isinstance(masses, list) and len(masses) == 3 and
                all(isinstance(x, (int, float)) and math.isfinite(x) and x >= 0 for x in masses))
    return {"status": "PassedLocalDMinus2Reference" if good else "NotVerified",
            "manifest": str(path), "sha256": fingerprint(path),
            "neutrino_masses_eV": masses if good else None}
