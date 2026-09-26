"""Pipeline-integrated numerical running and figure generation.

A config with

    "kind": "t3_pipeline_numerical_results_v1"

describes a concrete numerical benchmark point. The matching/RGE artifact paths
are taken from the current RunRecord, so the default config is not tied to a
particular output directory.

This module is numerical orchestration only. It reuses the existing UV,
intermediate, final-C5, final SM+Weinberg, neutrino-observable and figure
implementations.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from common.RunRecords import RunRecord
from Numerical.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.RGBetaEvaluator import load_rgbeta_payload
from Numerical.RunningDiagnostics import (
    _complex_matrix_payload,
    _intermediate_running_payload,
    _log_save_scales,
    _uv_running_payload,
    mixing_angles_from_pmns_abs,
)
from Numerical.RunningResultFigures import generate_running_result_figures
from Numerical.IntermediateWeinbergDiagnostics import sample_intermediate_direct_weinberg
from Numerical.FinalC5ContributionDiagnostics import evaluate_final_c5_contributions
from Numerical.BenchmarkSensitivity import run_t3_sensitivity_scan
from Numerical.ScalarThresholdBoundary import build_sm_weinberg_initial_conditions
from Numerical.ScanCLI import build_uv_state_from_config
from Numerical.T3Trajectory import run_renormalisable_t3_trajectory
from Numerical.WeinbergTrajectory import run_weinberg_trajectory
from physics.NeutrinoObservables import calculate_neutrino_observables
from physics.NeutrinoTrajectory import (
    charged_lepton_mass_basis_matrix,
    neutrino_mass_from_c5,
    scale_dependent_neutrino_observables,
)


PIPELINE_NUMERICAL_CONFIG_KIND = "t3_pipeline_numerical_results_v1"


def _load_json_object(path: Path) -> dict[str, Any]:
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object.")
    return payload


def is_pipeline_numerical_results_config(path: Path) -> bool:
    """Return whether a JSON file uses the pipeline-results numerical schema."""

    try:
        payload = _load_json_object(Path(path))
    except Exception:
        return False

    return payload.get("kind") == PIPELINE_NUMERICAL_CONFIG_KIND


def _record_artifact(
    record: RunRecord,
    summary_key: str,
    *,
    fallback_name: str | None = None,
) -> Path:
    relative = record.summary.get(summary_key)

    if relative:
        path = record.output_dir / str(relative)
    elif fallback_name is not None:
        path = record.output_dir / "data" / fallback_name
    else:
        raise KeyError(
            f"{record.name} summary does not contain {summary_key!r}."
        )

    if not path.is_file():
        raise FileNotFoundError(
            f"Required numerical artifact does not exist: {path}"
        )

    return path


def _validate_record_against_config(
    record: RunRecord,
    payload: dict[str, Any],
) -> None:
    representation = payload.get("representation")
    if not isinstance(representation, dict):
        raise ValueError("Numerical config requires a representation object.")

    expected = (
        int(representation["d_s1"]),
        int(representation["d_s2"]),
        int(representation["d_f"]),
        int(representation["alpha"]),
        bool(representation.get("shared_scalar", False)),
    )
    actual = (
        int(record.d_s1),
        int(record.d_s2),
        int(record.d_f),
        int(record.alpha),
        bool(record.shared_scalar),
    )

    if actual != expected:
        raise ValueError(
            "Numerical config representation does not match the pipeline model: "
            f"config={expected}, run={actual}."
        )

    if record.shared_scalar:
        raise NotImplementedError(
            "Pipeline-integrated running figures currently support the ordinary "
            "split-scalar branch only."
        )


def _build_default_uv_state(payload: dict[str, Any]):
    """Build the config's fixed UV point using the canonical ScanCLI state parser."""

    raw = dict(payload)
    raw.setdefault("scan", {"bindings": {}})
    config_view = SimpleNamespace(raw=raw)
    return build_uv_state_from_config(config_view, {}).validated()


def _run_and_build_payload(
    *,
    record: RunRecord,
    payload: dict[str, Any],
    config_path: Path,
) -> dict[str, Any]:
    _validate_record_against_config(record, payload)

    scales = payload.get("scales")
    if not isinstance(scales, dict):
        raise ValueError("Numerical config requires a scales object.")

    numerical = payload.get("numerical", {})
    if not isinstance(numerical, dict):
        raise ValueError("numerical must be an object.")

    mu_f = float(scales["mu_fermion_threshold_gev"])
    mu_s = float(scales["mu_scalar_threshold_gev"])
    mu_low = float(scales["mu_low_gev"])
    vev = float(numerical.get("vev_gev", 246.22))
    rtol = float(numerical.get("rtol", 1.0e-8))
    atol = float(numerical.get("atol", 1.0e-11))
    n_scale_points = int(numerical.get("n_scale_points", 120))
    ordering = str(payload.get("ordering", "NO")).upper()

    if ordering not in {"NO", "IO", "AUTO"}:
        raise ValueError("ordering must be NO, IO or AUTO.")

    uv_rgbeta_path = _record_artifact(
        record,
        "UVRGEFile",
        fallback_name="uv_rgbeta_rge.json",
    )
    intermediate_rgbeta_path = _record_artifact(
        record,
        "EFT1RenormalisableRGEFile",
        fallback_name="eft1_rgbeta_rge.json",
    )

    final_weinberg_relative = record.summary.get(
        "FinalWeinbergCoefficientFile"
    )
    if final_weinberg_relative:
        final_weinberg_path = record.output_dir / str(
            final_weinberg_relative
        )
        if not final_weinberg_path.is_file():
            raise FileNotFoundError(
                f"Required final C5 artifact does not exist: {final_weinberg_path}"
            )
    else:
        final_weinberg_path = _record_artifact(
            record,
            "WeinbergCoefficientFile",
            fallback_name="final_weinberg_coefficient.json",
        )

    uv_payload = load_rgbeta_payload(uv_rgbeta_path)
    intermediate_payload = load_rgbeta_payload(intermediate_rgbeta_path)
    uv_state = _build_default_uv_state(payload)

    renormalisable = run_renormalisable_t3_trajectory(
        uv_state,
        uv_payload,
        intermediate_payload,
        mu_fermion_threshold_gev=mu_f,
        mu_scalar_threshold_gev=mu_s,
        uv_save_scales_gev=_log_save_scales(
            uv_state.mu_gev,
            mu_f,
            n_scale_points,
        ),
        eft1_save_scales_gev=_log_save_scales(
            mu_f,
            mu_s,
            n_scale_points,
        ),
        rtol=rtol,
        atol=atol,
    )

    intermediate_weinberg = sample_intermediate_direct_weinberg(
        final_weinberg_path,
        renormalisable,
        _log_save_scales(
            mu_f,
            mu_s,
            n_scale_points,
        ),
    )

    c5_contributions = evaluate_final_c5_contributions(
        final_weinberg_path,
        renormalisable,
    )

    c5_evaluator = FinalC5TrajectoryEvaluator(final_weinberg_path)
    c5_threshold = np.asarray(
        c5_evaluator({}, renormalisable),
        dtype=complex,
    )

    initial = build_sm_weinberg_initial_conditions(
        renormalisable.final_sm_boundary,
        c5_threshold,
    )

    final_trajectory = run_weinberg_trajectory(
        initial,
        mu_s,
        mu_low,
        save_scales_gev=_log_save_scales(
            mu_s,
            mu_low,
            n_scale_points,
        ),
        rtol=rtol,
        atol=atol,
    )

    scale_points = scale_dependent_neutrino_observables(
        final_trajectory,
        vev_gev=vev,
        ordering=ordering,
    )

    running_scales: list[float] = []
    masses: list[list[float]] = []
    c5_abs: list[list[list[float]]] = []
    mnu_flavor_abs_ev: list[list[list[float]]] = []
    dm21: list[float] = []
    dm3l: list[float] = []
    s12: list[float] = []
    s13: list[float] = []
    s23: list[float] = []

    for index, point in enumerate(scale_points):
        obs = point.observables
        running = final_trajectory.point_at_index(index)

        running_scales.append(float(point.mu_gev))
        masses.append(np.asarray(obs.masses_ev, dtype=float).tolist())
        c5_abs.append(np.abs(point.c5).tolist())
        dm21.append(float(obs.delta_m21_sq_ev2))
        dm3l.append(
            float(
                obs.delta_m31_sq_ev2
                if obs.ordering.upper() == "NO"
                else obs.delta_m32_sq_ev2
            )
        )

        angle12, angle13, angle23 = mixing_angles_from_pmns_abs(
            obs.pmns_abs
        )
        s12.append(angle12)
        s13.append(angle13)
        s23.append(angle23)

        mass_flavor = charged_lepton_mass_basis_matrix(
            point.mass_matrix_gev,
            running.ye,
        )
        mnu_flavor_abs_ev.append(
            (np.abs(mass_flavor) * 1.0e9).tolist()
        )

    low = final_trajectory.final_point
    mnu_low = charged_lepton_mass_basis_matrix(
        neutrino_mass_from_c5(
            low.c5,
            vev_gev=vev,
        ),
        low.ye,
    )
    low_observables = calculate_neutrino_observables(
        mnu_low,
        ordering=ordering,
    )

    sensitivity_payload = None
    sensitivity_cfg = payload.get("sensitivity", {})
    if isinstance(sensitivity_cfg, dict) and sensitivity_cfg.get("enabled", False):
        target_value = sensitivity_cfg.get(
            "oscillation_target",
            "configs/nufit_6_1_ic24_no.json",
        )
        target_path = Path(target_value)
        if not target_path.is_absolute():
            target_path = Path(__file__).resolve().parents[1] / target_path

        sensitivity_payload = run_t3_sensitivity_scan(
            base_payload=payload,
            uv_rgbeta_path=uv_rgbeta_path,
            intermediate_rgbeta_path=intermediate_rgbeta_path,
            final_weinberg_path=final_weinberg_path,
            oscillation_target_path=target_path,
        )

    return {
        "status": "Success",
        "source": "PipelineNumericalResults",
        "numerical_config": str(Path(config_path).resolve()),
        "ordering": low_observables.ordering,
        "vev_gev": vev,
        "scales_gev": {
            "mu_uv": float(uv_state.mu_gev),
            "mu_fermion_threshold": mu_f,
            "mu_scalar_threshold": mu_s,
            "mu_low": mu_low,
        },
        "uv_running": _uv_running_payload(renormalisable.uv),
        "intermediate_running": _intermediate_running_payload(
            renormalisable.intermediate
        ),
        "intermediate_direct_weinberg": {
            "mu_gev": intermediate_weinberg.mu_gev.tolist(),
            "delta_c5_abs": intermediate_weinberg.delta_c5_abs.tolist(),
            "scheme": "fixed_order_one_loop_direct_LLSS_to_Weinberg",
            "note": (
                "This is the direct O(hbar) LLSS -> Weinberg contribution. "
                "LLSS self-running is not inserted because that would first "
                "affect C5 at O(hbar^2)."
            ),
        },
        "c5_threshold_contributions": c5_contributions.as_json_dict(),
        "scalar_threshold": {
            "c5": _complex_matrix_payload(c5_threshold),
        },
        "final_running": {
            "mu_gev": running_scales,
            "masses_ev": masses,
            "c5_abs": c5_abs,
            "mnu_charged_lepton_abs_ev": mnu_flavor_abs_ev,
            "delta_m21_sq_ev2": dm21,
            "delta_m3l_sq_ev2": dm3l,
            "sin2_theta12": s12,
            "sin2_theta13": s13,
            "sin2_theta23": s23,
        },
        "low_energy": {
            "masses_ev": np.asarray(
                low_observables.masses_ev,
                dtype=float,
            ).tolist(),
            "delta_m21_sq_ev2": float(
                low_observables.delta_m21_sq_ev2
            ),
            "delta_m31_sq_ev2": float(
                low_observables.delta_m31_sq_ev2
            ),
            "delta_m32_sq_ev2": float(
                low_observables.delta_m32_sq_ev2
            ),
            "pmns_abs": np.asarray(
                low_observables.pmns_abs,
                dtype=float,
            ).tolist(),
            "mnu_charged_lepton_basis_gev": _complex_matrix_payload(
                mnu_low
            ),
            "c5": _complex_matrix_payload(low.c5),
            "takagi_residual": float(
                low_observables.takagi_residual
            ),
        },
        "sensitivity": sensitivity_payload,
        "notes": {
            "benchmark_status": (
                "The default config is a reproducible numerical benchmark point, "
                "not an oscillation-data best fit."
            ),
            "intermediate_dimension_five_running": (
                "Not numerically sampled here. The scalar-only numerical state "
                "contains renormalisable parameters only; LLSS Wilson running "
                "remains in the RGE/matching pipeline."
            ),
        },
    }


def run_pipeline_numerical_results(
    record: RunRecord,
    numerical_config: Path,
) -> bool:
    """Generate running diagnostics and figures for one pipeline record."""

    summary = record.summary
    print(
        f"  {record.name}: starting integrated numerical results...",
        flush=True,
    )

    try:
        payload = _load_json_object(Path(numerical_config))
        if payload.get("kind") != PIPELINE_NUMERICAL_CONFIG_KIND:
            raise ValueError(
                "Config is not a pipeline numerical-results config."
            )

        diagnostics_payload = _run_and_build_payload(
            record=record,
            payload=payload,
            config_path=Path(numerical_config),
        )

        data_dir = record.output_dir / "data"
        figure_dir = record.output_dir / "figures"
        data_dir.mkdir(parents=True, exist_ok=True)
        figure_dir.mkdir(parents=True, exist_ok=True)

        diagnostics_path = data_dir / "running_diagnostics.json"
        temporary = diagnostics_path.with_name(
            diagnostics_path.name + ".tmp"
        )
        temporary.write_text(
            json.dumps(diagnostics_payload, indent=2),
            encoding="utf-8",
        )
        temporary.replace(diagnostics_path)

        generate_running_result_figures(
            diagnostics_path,
            figure_dir,
        )

    except Exception as exc:
        summary["PipelineNumericalResultsStatus"] = "Failed"
        summary["PipelineNumericalResultsError"] = str(exc)
        print(
            f"  {record.name}: integrated numerical results failed: {exc}"
        )
        return False

    summary["PipelineNumericalResultsStatus"] = "Success"
    summary["PipelineNumericalResultsConfig"] = str(
        Path(numerical_config).resolve()
    )
    summary["RunningDiagnosticsFile"] = (
        diagnostics_path.relative_to(record.output_dir).as_posix()
    )
    summary["RunningFiguresDirectory"] = (
        figure_dir.relative_to(record.output_dir).as_posix()
    )

    print(
        f"  {record.name}: integrated numerical results=Success"
        f" -> {figure_dir}",
        flush=True,
    )
    return True
