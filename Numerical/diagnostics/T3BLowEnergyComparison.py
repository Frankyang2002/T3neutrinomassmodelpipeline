"""Compare fitted T3-B fixed-order and exploratory coupled C5 at low energy.

Use the *same* production SM+Weinberg RGE and scalar-threshold SM
background for both C5 inputs. The upstream coupled result still uses a
common-scale logarithmic boundary and is not a validated sequential EFT.

This module never changes production C5, UV RGEs, or matching outputs.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np


_TINY = 1.0e-30


def _matrix_from_json(value: Any, *, name: str) -> np.ndarray:
    """Load complex matrices in the format emitted by the benchmark/control."""
    if not isinstance(value, dict) or "real" not in value or "imag" not in value:
        raise ValueError(f"{name}: expected {{'real':..., 'imag':...}}")
    real = np.asarray(value["real"], dtype=float)
    imag = np.asarray(value["imag"], dtype=float)
    if real.shape != (3, 3) or imag.shape != (3, 3):
        raise ValueError(f"{name}: expected 3x3 real and imaginary arrays")
    result = real + 1j * imag
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name}: matrix contains non-finite entries")
    return result


def _matrix_to_json(matrix: np.ndarray) -> dict:
    matrix = np.asarray(matrix, dtype=complex)
    if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
        raise ValueError("Expected finite complex 3x3 matrix")
    return {"real": matrix.real.tolist(), "imag": matrix.imag.tolist(),
            "frobenius_norm": float(np.linalg.norm(matrix))}


def _relative_residual(actual: np.ndarray, target: np.ndarray) -> float:
    actual = np.asarray(actual, dtype=complex)
    target = np.asarray(target, dtype=complex)
    return float(np.linalg.norm(actual - target) /
                 max(np.linalg.norm(target), _TINY))


def _require_majorana(matrix: np.ndarray, name: str, *, tolerance: float = 1e-9) -> None:
    matrix = np.asarray(matrix)
    if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
        raise ValueError(f"{name} must be a finite 3x3 complex matrix")
    error = _relative_residual(matrix, matrix.T)
    if error > tolerance:
        raise ValueError(f"{name} is not symmetric: relative residual {error:.3e}")


def _extract_threshold_matrices(control: dict, coupled: dict) -> tuple[np.ndarray, np.ndarray, dict]:
    """Validate status, common scheme and the non-double-counting substitution."""
    if control.get("status") != "FixedOrderControlComplete":
        raise ValueError("Control file is not a successful FixedOrderControlComplete output")
    if coupled.get("status") != "ExploratoryCoupledMatrixComplete":
        raise ValueError("Coupled file is not a successful ExploratoryCoupledMatrixComplete output")
    if coupled.get("threshold_scheme") != "common_scale_log_boundary":
        raise ValueError("Require the common_scale_log_boundary coupled output")
    if float(coupled.get("common_boundary_fixed_order_direct_residual", 1.0)) > 1e-8:
        raise ValueError("Upstream matching-log check did not pass")

    fixed = _matrix_from_json(control["fixed_order_c5_at_ms"], name="fixed C5(MS)")
    exploratory = _matrix_from_json(coupled["exploratory_total_c5_at_ms"],
                                    name="exploratory C5(MS)")
    direct = _matrix_from_json(coupled["production_direct_c5_at_ms"],
                               name="production direct C5(MS)")
    hard = _matrix_from_json(coupled["production_hard_c5_at_ms"],
                             name="production hard C5(MS)")
    evolved = _matrix_from_json(coupled["coupled_eft1_c5_at_ms"],
                                name="evolved intermediate C5(MS)")
    for name, matrix in (("fixed", fixed), ("exploratory", exploratory),
                         ("hard", hard), ("direct", direct), ("evolved", evolved)):
        _require_majorana(matrix, name)
    fixed_decomposition_error = _relative_residual(hard + direct, fixed)
    coupled_decomposition_error = _relative_residual(hard + evolved, exploratory)
    if max(fixed_decomposition_error, coupled_decomposition_error) > 1e-8:
        raise ValueError("Inconsistent hard/direct/coupled C5 decomposition at MS")
    relative_change = _relative_residual(exploratory, fixed)
    recorded = float(coupled["relative_difference_from_fixed_order"])
    if not np.isclose(relative_change, recorded, rtol=1e-7, atol=1e-10):
        raise ValueError("C5 comparison differs from recorded coupled output")
    return fixed, exploratory, {
        "fixed_order_decomposition_residual": fixed_decomposition_error,
        "coupled_decomposition_residual": coupled_decomposition_error,
        "relative_c5_change_at_ms": relative_change,
    }


def _check_saved_ms_background(control: dict, trajectory: Any,
                               *, relative_tolerance: float = 2e-6) -> dict:
    """Catch mismatched fitted configs/RGBeta artifacts using the saved MS sample."""
    samples = control["intermediate_running"]["samples"]
    if not samples:
        raise ValueError("Control JSON contains no intermediate running samples")
    saved = samples[-1]
    state = trajectory.eft1_threshold_state
    ms_saved = float(saved["mu_gev"])
    ms_rebuilt = float(state.mu_gev)
    if not np.isclose(ms_saved, ms_rebuilt, rtol=1e-10):
        raise ValueError("Control scalar scale and rebuilt trajectory disagree")

    def complex_scalar(data: Any) -> complex:
        if isinstance(data, list):
            if len(data) != 2:
                raise ValueError("Complex scalar must have two components")
            return complex(float(data[0]), float(data[1]))
        return complex(data)

    comparisons = {
        "g2": (complex_scalar(saved["g2"]), state.sm.g2),
        "lambdaH": (complex_scalar(saved["lambdaH"]), state.sm.lambdaH),
        "lambda12": (complex_scalar(saved["lambda12"]), state.lambda12),
        "lambda12Adj": (complex_scalar(saved["lambda12Adj"]), state.lambda12Adj),
        "lambdaT3": (complex_scalar(saved["lambdaT3"]), state.lambdaT3),
        "ye": (_matrix_from_json(saved["ye"], name="saved Ye"), state.sm.ye),
        "yu": (_matrix_from_json(saved["yu"], name="saved Yu"), state.sm.yu),
        "yd": (_matrix_from_json(saved["yd"], name="saved Yd"), state.sm.yd),
    }
    residuals = {}
    for name, (old, new) in comparisons.items():
        old = np.asarray(old, dtype=complex)
        new = np.asarray(new, dtype=complex)
        residual = float(np.linalg.norm(old - new) /
                         max(np.linalg.norm(old), np.linalg.norm(new), 1e-15))
        residuals[name] = residual
        if not np.isfinite(residual) or residual > relative_tolerance:
            raise ValueError(f"Rebuilt fitted MS background differs in {name}: {residual:.3e}")
    return residuals


def _mixing_angles(pmns_abs: np.ndarray) -> dict:
    """PDG mixing angles extracted from |U_PMNS|, independent of Majorana phases."""
    u = np.asarray(pmns_abs, dtype=float)
    if u.shape != (3, 3) or not np.all(np.isfinite(u)):
        raise ValueError("|U_PMNS| must be finite 3x3")
    s13_sq = float(np.clip(u[0, 2] ** 2, 0.0, 1.0))
    c13_sq = 1.0 - s13_sq
    if c13_sq < 1e-12:
        raise ValueError("PMNS angles undefined for theta13 ~ 90 degrees")
    s12_sq = float(np.clip(u[0, 1] ** 2 / c13_sq, 0.0, 1.0))
    s23_sq = float(np.clip(u[1, 2] ** 2 / c13_sq, 0.0, 1.0))
    return {
        "sin2_theta12": s12_sq, "sin2_theta13": s13_sq,
        "sin2_theta23": s23_sq,
        "theta12_deg": float(np.degrees(np.arcsin(np.sqrt(s12_sq)))),
        "theta13_deg": float(np.degrees(np.arcsin(np.sqrt(s13_sq)))),
        "theta23_deg": float(np.degrees(np.arcsin(np.sqrt(s23_sq)))),
    }


def _observables_for_evolved_state(state: Any, *, ordering: str, vev_gev: float) -> dict:
    """Use the repository's C5->mnu, charged-lepton basis and PMNS machinery."""
    from physics.NeutrinoTrajectory import (
        neutrino_mass_from_c5, charged_lepton_mass_basis_matrix,
    )
    from physics.NeutrinoObservables import calculate_neutrino_observables

    mass_gev = neutrino_mass_from_c5(state.K, vev_gev=vev_gev)
    charged_basis_mass = charged_lepton_mass_basis_matrix(mass_gev, state.ye)
    obs = calculate_neutrino_observables(charged_basis_mass, ordering=ordering)
    if not np.isfinite(obs.takagi_residual) or obs.takagi_residual > 1e-7:
        raise RuntimeError(f"Takagi diagonalisation failed: {obs.takagi_residual}")
    return {
        "ordering": obs.ordering,
        "masses_ev": np.asarray(obs.masses_ev, dtype=float).tolist(),
        "delta_m21_sq_ev2": float(obs.delta_m21_sq_ev2),
        "delta_m31_sq_ev2": float(obs.delta_m31_sq_ev2),
        "delta_m32_sq_ev2": float(obs.delta_m32_sq_ev2),
        "pmns_abs": np.asarray(obs.pmns_abs, dtype=float).tolist(),
        "mixing_angles": _mixing_angles(obs.pmns_abs),
        "takagi_residual": float(obs.takagi_residual),
    }


def _comparisons(a: dict, b: dict) -> dict:
    m_a = np.asarray(a["masses_ev"], dtype=float)
    m_b = np.asarray(b["masses_ev"], dtype=float)
    keys = ("delta_m21_sq_ev2", "delta_m31_sq_ev2", "delta_m32_sq_ev2")
    angle_keys = ("theta12_deg", "theta13_deg", "theta23_deg")
    return {
        "mass_differences_ev": (m_b - m_a).tolist(),
        "relative_mass_differences": [
            None if abs(float(x)) < _TINY else float((y-x)/x)
            for x, y in zip(m_a, m_b)
        ],
        "mass_squared_splitting_differences_ev2": {
            k: float(b[k] - a[k]) for k in keys
        },
        "mixing_angle_differences_deg": {
            k: float(b["mixing_angles"][k] - a["mixing_angles"][k])
            for k in angle_keys
        },
    }


def run(*, config: Path, uv_rgbeta: Path, eft1_rgbeta: Path,
        control: Path, coupled: Path, mu_low_gev: float | None = None,
        ordering: str = "NO", vev_gev: float | None = None) -> dict:
    """Run two independent SM+Weinberg evolutions with identical SM boundary."""
    from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
    from Numerical.fitting.ScanCLI import build_uv_state_from_config
    from Numerical.running.T3Trajectory import run_renormalisable_t3_trajectory
    from Numerical.running.ScalarThresholdBoundary import build_sm_weinberg_initial_conditions
    from Numerical.running.SMWeinbergEvolution import evolve_sm_weinberg

    config_data = json.loads(Path(config).read_text(encoding="utf-8"))
    control_data = json.loads(Path(control).read_text(encoding="utf-8"))
    coupled_data = json.loads(Path(coupled).read_text(encoding="utf-8"))
    fixed, candidate, consistency = _extract_threshold_matrices(control_data, coupled_data)
    scales = config_data["scales"]
    mu_f = float(scales["mu_fermion_threshold_gev"])
    mu_s = float(scales["mu_scalar_threshold_gev"])
    mu_low = float(mu_low_gev if mu_low_gev is not None else scales["mu_low_gev"])
    vev = float(vev_gev if vev_gev is not None else
                config_data.get("numerical", {}).get("vev_gev", 246.22))
    if not (np.isfinite(mu_f) and np.isfinite(mu_s) and np.isfinite(mu_low)
            and mu_f > mu_s > mu_low > 0):
        raise ValueError("Need finite mu_F > mu_S > mu_low > 0")
    if not np.isfinite(vev) or vev <= 0:
        raise ValueError("vev_gev must be positive finite")
    if ordering.upper() not in ("NO", "IO", "AUTO"):
        raise ValueError("ordering must be NO, IO, or AUTO")
    if not np.isclose(mu_f, float(control_data["mu_f_gev"]), rtol=1e-12):
        raise ValueError("Config mu_F disagrees with control")
    if not np.isclose(mu_s, float(control_data["mu_s_gev"]), rtol=1e-12):
        raise ValueError("Config mu_S disagrees with control")
    coupled_interval = coupled_data.get("threshold_scheme_details", {}).get("evolution_interval_gev")
    if coupled_interval is not None and not np.allclose(coupled_interval,
            [mu_f, mu_s], rtol=1e-12, atol=0):
        raise ValueError("Coupled threshold interval does not match configuration")

    raw = dict(config_data)
    raw.setdefault("scan", {"bindings": {}})
    initial = build_uv_state_from_config(SimpleNamespace(raw=raw), {}).validated()
    rep = initial.representation
    if (rep.d_s1, rep.d_s2, rep.d_f, rep.alpha, rep.shared_scalar) != (2, 2, 1, -1, False):
        raise ValueError("This comparison is restricted to ordinary T3-B alpha=-1")
    numerics = config_data.get("numerical", {})
    ren = run_renormalisable_t3_trajectory(
        initial, load_rgbeta_payload(uv_rgbeta), load_rgbeta_payload(eft1_rgbeta),
        mu_fermion_threshold_gev=mu_f, mu_scalar_threshold_gev=mu_s,
        rtol=float(numerics.get("rtol", 1e-8)),
        atol=float(numerics.get("atol", 1e-11)),
    )
    consistency["ms_background_residuals"] = _check_saved_ms_background(control_data, ren)

    results = {}
    c5_low = {}
    for label, at_ms in (("fixed_order", fixed), ("exploratory_coupled", candidate)):
        boundary = build_sm_weinberg_initial_conditions(ren.final_sm_boundary, at_ms)
        evolved = evolve_sm_weinberg(boundary, mu_s, mu_low)
        if not evolved.solver_success:
            raise RuntimeError(f"SM+Weinberg solver failed for {label}")
        _require_majorana(evolved.K, f"{label} low-scale C5")
        c5_low[label] = np.asarray(evolved.K, dtype=complex)
        results[label] = {
            "observables": _observables_for_evolved_state(
                evolved, ordering=ordering.upper(), vev_gev=vev),
            "solver_nfev": int(evolved.nfev),
            "c5_at_low_scale": _matrix_to_json(evolved.K),
        }

    fixed_obs = results["fixed_order"]["observables"]
    coupled_obs = results["exploratory_coupled"]["observables"]
    low_c5_ratio = _relative_residual(c5_low["exploratory_coupled"], c5_low["fixed_order"])
    return {
        "status": "ExploratoryLowEnergyComparisonComplete",
        "model": control_data.get("model", "T3-B alpha=-1"),
        "mu_fermion_threshold_gev": mu_f,
        "mu_scalar_threshold_gev": mu_s,
        "mu_low_gev": mu_low,
        "vev_gev": vev,
        "ordering_requested": ordering.upper(),
        "source_inputs": {"config": str(config), "control": str(control),
                          "coupled": str(coupled)},
        "threshold_scheme": coupled_data["threshold_scheme"],
        "ms_input_checks": consistency,
        "relative_c5_difference_at_ms": consistency["relative_c5_change_at_ms"],
        "relative_c5_difference_at_low_scale": low_c5_ratio,
        "fixed_order": results["fixed_order"],
        "exploratory_coupled": results["exploratory_coupled"],
        "differences": _comparisons(fixed_obs, coupled_obs),
        "validated_physical_prediction_available": False,
        "limitations": [
            "Coupled intermediate-EFT result is exploratory; full three-generation operator closure not established",
            "Common-scale heavy matching logs are not complete sequential heavy-fermion decoupling",
            "No complete higher-order scalar-threshold subtraction/matching proof",
            "SM+Weinberg unbroken-phase running is extended to the configured low scale; electroweak threshold matching is not added",
        ],
    }


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--uv-rgbeta", required=True, type=Path)
    p.add_argument("--eft1-rgbeta", required=True, type=Path)
    p.add_argument("--control", required=True, type=Path)
    p.add_argument("--coupled", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--mu-low", type=float, default=None,
                   help="Default: scales.mu_low_gev in fitted config (100 GeV)")
    p.add_argument("--ordering", choices=("NO", "IO", "AUTO"), default="NO")
    p.add_argument("--vev-gev", type=float, default=None)
    args = p.parse_args(argv)
    payload = run(config=args.config, uv_rgbeta=args.uv_rgbeta,
                  eft1_rgbeta=args.eft1_rgbeta, control=args.control,
                  coupled=args.coupled, mu_low_gev=args.mu_low,
                  ordering=args.ordering, vev_gev=args.vev_gev)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    summary = {
        "status": payload["status"], "mu_low_gev": payload["mu_low_gev"],
        "relative_c5_difference_at_ms": payload["relative_c5_difference_at_ms"],
        "relative_c5_difference_at_low_scale": payload["relative_c5_difference_at_low_scale"],
        "fixed_order": payload["fixed_order"]["observables"],
        "exploratory_coupled": payload["exploratory_coupled"]["observables"],
        "differences": payload["differences"],
        "validated_physical_prediction_available": False,
    }
    print(json.dumps(summary, indent=2))
    print(json.dumps({"output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
