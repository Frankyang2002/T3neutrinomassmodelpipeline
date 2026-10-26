"""Conditional RGBeta -> Matchete C5 response for ordinary B/C(-1).

NO complete UV-to-EFT matching is performed. The RGBeta-to-Matchete
complex coupling map remains unverified: both independent couplings are
scanned over {z, -z, conjugate(z), -conjugate(z)}. The coupled UV endpoint
is ABOVE the heavy thresholds while the Matchete hard kernels use their
individual MF matching scales. Toy MF, Yukawas, and scalar masses remain
frozen. The output is a sensitivity estimate, never a physical prediction.
"""
from __future__ import annotations

import itertools
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from Numerical.experimental.GeneralT3C5Pilot import ACTIVE, MODELS, ResponseKernel, evaluate_response

RGBETA_BY_MAT = {
    "lambdaT3ExtraQ020020I1": "lambdaHHdagS2S2",
    "lambdaT3ExtraQ020200I1": "lambdaHHdagS1barS1bar",
}
MAP_CHOICES = ("same", "negative", "conjugated", "negative_conjugated")


def _parse_complex(value: Any, name: str) -> complex:
    if not isinstance(value, Mapping) or set(value) != {"real", "imag"}:
        raise ValueError(f"{name}: expected only real/imag keys")
    z = complex(float(value["real"]), float(value["imag"]))
    if not np.isfinite([z.real, z.imag]).all():
        raise ValueError(f"{name}: nonfinite value")
    return z


def _encode(z: complex) -> dict[str, float]:
    return {"real": float(complex(z).real), "imag": float(complex(z).imag)}


def _transform(z: complex, convention: str) -> complex:
    if convention not in MAP_CHOICES:
        raise ValueError("Unknown complex mapping scenario")
    x = z.conjugate() if "conjugated" in convention else z
    return -x if "negative" in convention else x


def _validate_tensor_mapping(document: Mapping[str, Any], model: str) -> None:
    if document.get("status") != "DoubletGroupTensorRatiosVerified":
        raise ValueError("Stage 13 tensor verification must have succeeded")
    if document.get("full_lagrangian_coupling_identification_verified") is not False:
        raise ValueError("Expected a tensor-only, not a certified full coupling map")
    matches = [r for r in document.get("mappings", []) if r.get("model") == model]
    found = {}
    for row in matches:
        name = row.get("matchete_coupling")
        if name not in RGBETA_BY_MAT or name in found:
            raise ValueError("Duplicate or unknown Matchete tensor mapping")
        if row.get("rgbeta_coupling") != RGBETA_BY_MAT[name]:
            raise ValueError("RGBeta field-structure mapping differs from certified convention")
        if float(row.get("rgbeta_over_matchete_tensor_factor", float("nan"))) != 1.0:
            raise ValueError("Nonunit or invalid RGBeta-to-Matchete tensor factor")
        if row.get("full_lagrangian_coupling_identification_verified") is not False:
            raise ValueError("Expected tensor-only mapping row")
        found[name] = True
    if set(found) != set(RGBETA_BY_MAT):
        raise ValueError("Missing B/C(-1) tensor mapping")


def _induced_rgbeta_couplings(document: Mapping[str, Any], model: str,
                              parameters: Mapping[str, Any]) -> dict[str, complex]:
    if document.get("status") != "CoupledUVOneLoopRGBetaOnly" or document.get("model") != model:
        raise ValueError("Expected verified model-specific Stage 13 UV result")
    if document.get("solver_success") is not True:
        raise ValueError("Stage 13 ODE integration not successful")
    if document.get("full_matchete_rgbeta_coupling_mapping_verified") is not False:
        raise ValueError("Stage 13 must not claim a full coupling identification")
    if document.get("includes_mass_threshold_matching") is not False:
        raise ValueError("Stage 13 UV-only threshold boundary changed")
    initial_mf = np.asarray(parameters["MF"], dtype=float)
    if initial_mf.ndim != 1 or initial_mf.size != 3 or np.any(initial_mf <= 0):
        raise ValueError("Expected three initial positive heavy masses")
    start = float(document.get("uv_scale_gev", float("nan")))
    endpoint = float(document.get("endpoint_gev", float("nan")))
    if not np.isfinite([start, endpoint]).all() or not (start > endpoint > max(initial_mf)):
        raise ValueError("Expected UV-only endpoint above all initial heavy masses")
    if not np.isclose(start, float(parameters["scale_uv_gev"]), rtol=1e-12):
        raise ValueError("Stage 13 and toy config UV scales differ")
    records = document.get("active_couplings")
    if not isinstance(records, list) or len(records) != len(RGBETA_BY_MAT):
        raise ValueError("Expected two Stage 13 active quartics")
    found: dict[str, complex] = {}
    for row in records:
        name = row.get("coupling")
        if name not in RGBETA_BY_MAT.values() or name in found:
            raise ValueError("Duplicate or unknown RGBeta quartic")
        initial = _parse_complex(row.get("uv_initial_rgbeta"), f"{name} initial")
        end = _parse_complex(row.get("uv_endpoint_rgbeta"), f"{name} endpoint")
        delta = _parse_complex(row.get("change_coupled"), f"{name} change")
        if not np.isclose(end - initial, delta, rtol=1e-9, atol=1e-25):
            raise ValueError("Stage 13 reported delta inconsistent with endpoint")
        if abs(initial) > 1e-25:
            raise ValueError("Induced-from-zero analysis requires zero UV initial quartic")
        found[name] = delta
    if set(found) != set(RGBETA_BY_MAT.values()):
        raise ValueError("Incomplete Stage 13 active quartics")
    return found


def _validate_reference(reference: Mapping[str, Any], model: str) -> float:
    if reference.get("model") != model or reference.get("status") != "ExperimentalHardBoundaryC5":
        raise ValueError("Expected Stage 11 model-specific Matchete pilot")
    if reference.get("units") != "GeV^-1":
        raise ValueError("Reference hard-C5 has wrong units")
    number = float(reference["delta_c5"]["frobenius_norm"])
    if not np.isfinite(number) or number <= 0:
        raise ValueError("Reference Matchete C5 norm must be positive")
    return number


def conditional_uv_to_c5(
    model: str,
    parameters: Mapping[str, Any],
    stage13: Mapping[str, Any],
    tensor_map: Mapping[str, Any],
    kernels: tuple[ResponseKernel, ...],
    reference: Mapping[str, Any],
    *,
    vev_gev: float = 246.22,
) -> dict[str, Any]:
    """Evaluate all 16 *hypothetical* discrete complex-coupling maps.

    Does not insert anything into production numerical/RGE code. In particular
    no threshold matching, UV-to-IR Lagrangian sign, or Yukawa running is proved.
    """
    if model not in MODELS:
        raise ValueError("Only ordinary B/C(-1) supported")
    _validate_tensor_mapping(tensor_map, model)
    induced = _induced_rgbeta_couplings(stage13, model, parameters)
    reference_norm = _validate_reference(reference, model)
    vev = float(vev_gev)
    if not np.isfinite(vev) or vev <= 0:
        raise ValueError("Invalid Higgs vev")
    if {k.coupling for k in kernels} != set(RGBETA_BY_MAT) or len(kernels) != 2:
        raise ValueError("Exactly two proven Matchete derivative kernels required")
    names = tuple(sorted(RGBETA_BY_MAT))
    estimates = []
    for conventions in itertools.product(MAP_CHOICES, repeat=2):
        copied = dict(parameters)
        copied["matchete_extra_couplings"] = {
            name: _encode(_transform(induced[RGBETA_BY_MAT[name]], convention))
            for name, convention in zip(names, conventions)
        }
        outcome = evaluate_response(kernels, copied)
        c5 = float(outcome["delta_c5"]["frobenius_norm"])
        singular = np.linalg.svd(
            -(vev ** 2 / 2.0) * 1e9 *
            (np.asarray(outcome["delta_c5"]["real"]) +
             1j * np.asarray(outcome["delta_c5"]["imag"])),
            compute_uv=False,
        )
        estimates.append({
            "mapping_choices": dict(zip(names, conventions)),
            "c5_norm_gev_inverse": c5,
            "extra_only_mnu_norm_ev": (vev ** 2 / 2.0) * 1e9 * c5,
            "extra_only_takagi_singular_values_ev": sorted(float(s) for s in singular),
            "relative_to_stage11_matchete_toy_c5_norm": c5 / reference_norm,
        })
    summary = {
        "minimum_extra_only_mnu_norm_ev": min(e["extra_only_mnu_norm_ev"] for e in estimates),
        "maximum_extra_only_mnu_norm_ev": max(e["extra_only_mnu_norm_ev"] for e in estimates),
        "minimum_relative_to_stage11_reference": min(e["relative_to_stage11_matchete_toy_c5_norm"] for e in estimates),
        "maximum_relative_to_stage11_reference": max(e["relative_to_stage11_matchete_toy_c5_norm"] for e in estimates),
    }
    return {
        "status": "ConditionalUVInducedC5SensitivityOnly",
        "model": model,
        "induced_from_zero_rgbeta_quartics": {k: _encode(v) for k, v in induced.items()},
        "stage13_uv_endpoint_gev": float(stage13["endpoint_gev"]),
        "number_of_hypothetical_complex_maps": len(estimates),
        "possible_discrete_mappings": list(MAP_CHOICES),
        "mapping_cases": estimates,
        "range": summary,
        "reference_matchete_toy_c5_norm_gev_inverse": reference_norm,
        "full_lagrangian_mapping_verified": False,
        "threshold_consistent_matching": False,
        "rg_improved_neutrino_mass": False,
        "fixed_toy_mass_and_yukawa_inputs": True,
        "disclaimer": (
            "Illustrative sensitivity only: RGBeta UV quartics evolved to an above-threshold "
            "endpoint; hard Matchete kernels at individual MF with frozen toy MF, Yukawas "
            "and scalar masses. Four discrete sign/conjugation possibilities per quartic "
            "are sampled, NOT an exhaustive space of complex convention phases. "
            "No matched running, decoupling, or full neutrino spectrum is calculated."
        ),
    }
