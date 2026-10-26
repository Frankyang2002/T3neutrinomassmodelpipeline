"""Extra-only neutrino mass diagnostics for the experimentally matched B/C(-1) C5.

No baseline C5, RGBeta-to-Matchete map, RGE evolution, scalar mixing,
or mass-hierarchy corrections are inferred. For the project's convention:
    m_nu [eV] = -(vev_gev**2 / 2) * C5 [GeV^-1] * 1e9.
The Takagi values are the singular values of the symmetric Majorana matrix;
no PMNS angles or physical oscillation fit are inferred from extra-only terms.
"""
from __future__ import annotations

import json
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np

EXPECTED_MODELS = frozenset({"B_m1", "C_m1"})
EXPECTED_COUPLINGS = frozenset({
    "lambdaT3ExtraQ020020I1", "lambdaT3ExtraQ020200I1"
})
PHASES_DEG = (0, 90, 180, 270)


def _as_matrix(payload: dict[str, Any], *, context: str) -> np.ndarray:
    if not isinstance(payload, dict) or not {"real", "imag", "frobenius_norm"} <= payload.keys():
        raise ValueError(f"{context}: expected numerical complex-matrix block")
    r = np.asarray(payload["real"], dtype=float)
    i = np.asarray(payload["imag"], dtype=float)
    if r.shape != (3, 3) or i.shape != r.shape:
        raise ValueError(f"{context}: expected 3x3 real and imaginary components")
    matrix = np.asarray(r + 1j * i, dtype=np.complex128)
    if not np.isfinite(matrix).all():
        raise ValueError(f"{context}: nonfinite matrix")
    if not np.allclose(matrix, matrix.T, rtol=1e-12, atol=1e-26):
        raise ValueError(f"{context}: matrix not symmetric")
    recorded = float(payload["frobenius_norm"])
    measured = float(np.linalg.norm(matrix, ord="fro"))
    if not np.isfinite(recorded) or not np.isclose(recorded, measured, rtol=1e-9, atol=1e-26):
        raise ValueError(f"{context}: Frobenius norm does not match components")
    return matrix


def _matrix_payload(matrix: np.ndarray) -> dict[str, Any]:
    return {
        "real": np.asarray(matrix.real, dtype=float).tolist(),
        "imag": np.asarray(matrix.imag, dtype=float).tolist(),
        "frobenius_norm": float(np.linalg.norm(matrix, ord="fro")),
    }


def _eigen_data(c5: np.ndarray, vev_gev: float) -> dict[str, Any]:
    # Same convention as physics.NeutrinoMass.numerical_neutrino_mass_matrix.
    mnu_ev = -(vev_gev ** 2 / 2.0) * c5 * 1.0e9
    # For a complex symmetric Majorana matrix these are the physical
    # Takagi singular values, sorted; we do NOT infer a PMNS rotation.
    singular_values = np.sort(np.linalg.svd(mnu_ev, compute_uv=False))
    return {
        "c5": _matrix_payload(c5),
        "extra_neutrino_mass_matrix_ev": _matrix_payload(mnu_ev),
        "extra_only_takagi_singular_values_ev": [float(x) for x in singular_values],
        "sum_singular_values_ev": float(np.sum(singular_values)),
    }


def analyze_extra_mass_response(payload: dict[str, Any], *, vev_gev: float = 246.22) -> dict[str, Any]:
    """Compute scenario responses solely from an authenticated Stage-11 C5 payload."""
    if not isinstance(payload, dict) or payload.get("status") != "ExperimentalHardBoundaryC5":
        raise ValueError("Stage-11 hard C5 pilot status required")
    model = payload.get("model")
    if model not in EXPECTED_MODELS:
        raise ValueError("Only B_m1 and C_m1 are supported")
    if payload.get("units") != "GeV^-1" or payload.get("finite_mass_ratio_logs_preserved") is not True:
        raise ValueError("Incompatible Matchete hard C5 units or matching convention")
    vev = float(vev_gev)
    if not np.isfinite(vev) or vev <= 0:
        raise ValueError("vev_gev must be finite and positive")
    contributions = payload.get("contributions")
    if not isinstance(contributions, dict) or set(contributions) != EXPECTED_COUPLINGS:
        raise ValueError("Expected exactly two verified C5-active quartic contributions")
    matrices = {k: _as_matrix(v, context=k) for k, v in contributions.items()}
    delta = _as_matrix(payload.get("delta_c5"), context="delta_c5")
    if not np.allclose(sum(matrices.values()), delta, rtol=1e-9, atol=1e-26):
        raise ValueError("Individual quartic contributions do not sum to delta_c5")

    kernels = payload.get("kernels", [])
    if not isinstance(kernels, list) or len(kernels) != 2:
        raise ValueError("Expected two Matchete coefficient metadata records")
    orientation = {}
    for kernel in kernels:
        if not isinstance(kernel, dict):
            raise ValueError("Invalid Matchete kernel metadata")
        name = kernel.get("coupling")
        if name not in EXPECTED_COUPLINGS or name in orientation or not isinstance(kernel.get("conjugated"), bool):
            raise ValueError("Incomplete or repeated Matchete coupling orientation")
        orientation[name] = kernel["conjugated"]
    if set(orientation) != EXPECTED_COUPLINGS:
        raise ValueError("Missing coupling orientation")

    names = sorted(EXPECTED_COUPLINGS)
    a, b = (matrices[name] for name in names)
    scenarios = {}
    for label, factors in {
        "both_as_input": {names[0]: 1.0, names[1]: 1.0},
        "only_first": {names[0]: 1.0, names[1]: 0.0},
        "only_second": {names[0]: 0.0, names[1]: 1.0},
        "extra_couplings_zero": {names[0]: 0.0, names[1]: 0.0},
    }.items():
        total = sum((factors[k] * matrices[k] for k in names), np.zeros((3, 3), complex))
        scenarios[label] = _eigen_data(total, vev)
    # Rotate the *physical quartic coupling* relative to its input value;
    # the C5 atom rotates with the conjugate phase for a Bar[Coupling] atom.
    scan = []
    for p0, p1 in product(PHASES_DEG, repeat=2):
        angles = (p0, p1)
        terms = []
        for name, deg in zip(names, angles):
            multiplier = np.exp(1j * np.deg2rad(deg))
            terms.append(np.conjugate(multiplier) * matrices[name]
                         if orientation[name] else multiplier * matrices[name])
        total = sum(terms, np.zeros((3, 3), complex))
        mass_spectrum = _eigen_data(total, vev)
        scan.append({
            "physical_coupling_phase_shifts_deg": dict(zip(names, angles)),
            "c5_norm_gev_inverse": mass_spectrum["c5"]["frobenius_norm"],
            "mnu_norm_ev": mass_spectrum["extra_neutrino_mass_matrix_ev"]["frobenius_norm"],
            "takagi_singular_values_ev": mass_spectrum["extra_only_takagi_singular_values_ev"],
        })

    combined_norm_sq = float(np.linalg.norm(delta, "fro") ** 2)
    individual_sq = float(sum(np.linalg.norm(mat, "fro") ** 2 for mat in matrices.values()))
    cross = float(2 * np.vdot(a, b).real)
    if not np.isclose(combined_norm_sq, individual_sq + cross, rtol=1e-9, atol=1e-40):
        raise RuntimeError("C5 interference identity failed")
    min_s = min(scan, key=lambda x: x["mnu_norm_ev"])
    max_s = max(scan, key=lambda x: x["mnu_norm_ev"])
    return {
        "status": "ExtraOnlyNeutrinoMassImpact",
        "model": model,
        "vev_gev": vev,
        "mass_convention": "m_nu[eV]=-(v^2/2)*C5[GeV^-1]*1e9",
        "couplings_in_order": names,
        "conjugated_coupling_atoms": orientation,
        "scenarios": scenarios,
        "c5_norm_squared_decomposition_gev_minus2": {
            "sum_individual_norm_squared": individual_sq,
            "interference_2ReTrAstarB": cross,
            "combined_norm_squared": combined_norm_sq,
        },
        "phase_grid_degrees": list(PHASES_DEG),
        "phase_scan": scan,
        "minimum_mass_norm_phase": min_s,
        "maximum_mass_norm_phase": max_s,
        "limitations": [
            "Only extra hard C5 evaluated, not original C5 or full neutrino spectrum.",
            "Phase scan varies toy Matchete couplings; it is not a fitted oscillation scan.",
            "Takagi values are singular values of extra-only m_nu; not physical mass predictions.",
            "No Matchete/RGBeta conversion, coupled RGE, or quadratic mixing is assumed.",
        ],
    }


def load_and_analyze(path: Path, *, vev_gev: float = 246.22) -> dict[str, Any]:
    return analyze_extra_mass_response(json.loads(Path(path).read_text(encoding="utf-8-sig")), vev_gev=vev_gev)
