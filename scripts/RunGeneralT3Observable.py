"""All-16 general-T3 numerical observable orchestration (conditional one-scale C5).

Run from project root:

    python -X utf8 -m scripts.RunGeneralT3Observable --all \
      --canonical-run-id full16_general_v1 \
      --generate-symbolic --generate-matching \
      --output-root output/general_observable/all16_v1

The general model inventory and beta functions come from the existing canonical
full-general symbolic exports. The UV and EFT1 RGEs use the existing numerical
runners, and the Matchete one-scale C5 uses GeneralT3FlavourC5. This module
adds no new beta functions or matching formulas.

IMPORTANT: The C5 is *conditional direct one-scale* matching. It is NOT
fermion-first EFT1 Wilson matching, Wilson evolution, a full renormalized
sequential scalar threshold, or a certified Matchete/RGBeta coupling map.
B_-1 and C_-1 reuse their archived one-scale mass-eigenstate Matchete
coefficient with a CONDITIONAL signed quartic dictionary.  Their physical
Matchete/RGBeta normalization is NOT certified.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import sys
import traceback
from typing import Any, Mapping

import numpy as np

from scripts.RunT3GeneralMatcheteMatching import ALL as MODELS, MIXED, tag

ROOT = Path(__file__).resolve().parents[1]
FORMAT = "t3-general-all16-numerical-observables-v1"
C5_STATUS = "ProvisionalDirectMatching_NotFullEFT1"
SM_NAMES = ("gY", "g2", "g3", "lambdaH")


def json_load(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return obj


def json_save(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def source_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def encode(obj: Any) -> Any:
    if isinstance(obj, np.ndarray):
        return encode(obj.tolist())
    if isinstance(obj, np.generic):
        return encode(obj.item())
    if isinstance(obj, complex):
        return {"real": float(obj.real), "imag": float(obj.imag)}
    if isinstance(obj, (tuple, list)):
        return [encode(x) for x in obj]
    if isinstance(obj, Mapping):
        return {str(k): encode(v) for k, v in obj.items()}
    return obj


def complex_matrix(matrix: np.ndarray) -> list[list[dict[str, float]]]:
    array = np.asarray(matrix, dtype=complex)
    if array.shape != (3, 3) or not np.isfinite(array).all():
        raise ValueError("Expected finite complex 3x3 matrix")
    return encode(array)


def flatten_parameter(name: str, value: Any):
    """Long-form CSV rows: parameter key, real value, imaginary value."""
    if isinstance(value, Mapping) and set(value) == {"real", "imag"}:
        re, im = value["real"], value["imag"]
        if not isinstance(re, (int, float)) or not isinstance(im, (int, float)):
            raise ValueError(f"Non-scalar complex value in {name}")
        yield name, float(re), float(im)
    elif isinstance(value, (list, tuple)):
        for i, entry in enumerate(value, start=1):
            yield from flatten_parameter(f"{name}[{i}]", entry)
    elif isinstance(value, (float, int)):
        yield name, float(value), 0.0
    else:
        raise ValueError(f"Unexpected serialized parameter {name}: {type(value)}")


def save_samples(path: Path, result: Any, state_values: Any) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow(["mu_gev", "parameter", "real", "imag"])
        for i, mu in enumerate(result.mu_gev):
            values = state_values(result.state_at_index(i))
            for name, value in sorted(values.items()):
                for entry, re, im in flatten_parameter(name, value):
                    if not np.isfinite((re, im)).all():
                        raise ValueError(f"Nonfinite running parameter {entry} at mu={mu}")
                    writer.writerow([float(mu), entry, re, im])
                    count += 1
    return count


def diagonal(values: tuple[float, float, float]) -> list[list[float]]:
    return [[float(values[i]) if i == j else 0.0 for j in range(3)]
            for i in range(3)]


def illustrative_overrides(uv_payload: Mapping[str, Any], mu_f: float,
                           mu_s: float) -> dict[str, Any]:
    """Controlled CP-real benchmark; all allowed UV parameter keys still exist.

    Only explicit starting values are chosen here. No interaction is removed
    from the symbolic inventory or from the RGEs. The CP-real subspace does
    not validate complex self-bilinear matching.
    """
    from Numerical.core.GeneralCouplingLayout import spec_for

    betas = uv_payload.get("report_betas")
    if not isinstance(betas, Mapping):
        raise ValueError("Missing UV coupling registry")
    ans: dict[str, Any] = {
        "g3": 0.765492,
        "yu": diagonal((1e-5, 0.003, 0.5)),
        "y1": diagonal((0.0004, 0.0004, 0.0004)),
        "y2": diagonal((0.0004, 0.0004, 0.0004)),
        "MF": diagonal((1.5 * mu_f, 1.8 * mu_f, 2.0 * mu_f)),
    }
    for name in betas:
        if name in ("y1Conj", "y2Conj"):
            ans[name] = diagonal((0.0004, 0.0004, 0.0004))
        elif name.startswith("muT3General"):
            ans[name] = 0.02 * mu_s ** 2
        elif (spec_for(name).kind == "complex"
              and name not in ans):
            ans[name] = 0.01
    return {key: val for key, val in ans.items() if key in betas}


def user_overrides(path: Path | None, model: str) -> dict[str, Any]:
    if path is None:
        return {}
    obj = json_load(path)
    if "shared" not in obj and "models" not in obj:
        return obj
    if set(obj) - {"shared", "models"}:
        raise ValueError("Override file with shared/models has unknown top-level keys")
    shared = obj.get("shared", {})
    models = obj.get("models", {})
    if not isinstance(shared, dict) or not isinstance(models, dict):
        raise ValueError("shared and models must be JSON objects")
    if set(models) - set(MODELS):
        raise ValueError("Unknown model override keys: " + str(set(models) - set(MODELS)))
    local = models.get(model, {})
    if not isinstance(local, dict):
        raise ValueError(f"Overrides for {model} must be an object")
    return {**shared, **local}


def hierarchy_diagnostics(uv_at_f: Any, eft_at_s: Any,
                          mu_f: float, mu_s: float) -> dict[str, Any]:
    masses = np.sort(np.linalg.svd(np.asarray(uv_at_f.MF, dtype=complex),
                                   compute_uv=False))
    scalar_sq = (float(eft_at_s.mS1Sq), float(eft_at_s.mS2Sq))
    scalar_mass = [float(np.sqrt(v)) if v > 0 else None for v in scalar_sq]
    return {
        "grouped_fermion_matching_scale_gev": mu_f,
        "heavy_fermion_singular_masses_at_F_gev": masses.tolist(),
        "all_fermion_masses_above_mu_F": bool(np.all(masses > mu_f)),
        "scalar_matching_scale_gev": mu_s,
        "unbroken_scalar_mass_parameters_at_S_gev": scalar_mass,
        "all_unsplit_scalar_masses_above_mu_S": bool(all(
            v is not None and v > mu_s for v in scalar_mass)),
        "scalar_mixed_or_self_bilinear_mass_basis_checked": False,
        "note": "Scalar masses are unbroken unsplit parameters, not all physical eigenmasses",
    }


def coefficient_path(model: str, matching_dir: Path, generate: bool,
                     wolframscript: str) -> Path:
    folder = matching_dir / tag(model)
    candidate = folder / "c5_coefficient.txt"
    status_path = folder / "general_matching_status.json"

    def ready() -> bool:
        if not status_path.is_file() or not candidate.is_file() or candidate.stat().st_size == 0:
            return False
        state = json_load(status_path)
        return (state.get("model") == model and state.get("status") == "Matched"
                and state.get("matching_scheme") in
                ("FullGeneralUnsplitMatcheteOneScaleCandidate",
                 "FullGeneralSelfBilinearRealMassEigenfieldsOneScaleCandidate"))

    if not ready() and generate:
        from scripts.RunT3GeneralMatcheteMatching import run_model
        run_model(model, matching_dir, wolframscript, reuse=False)
    if not ready():
        raise FileNotFoundError(f"{model}: no validated general one-scale Matchete result at {candidate}; "
                                "use --generate-matching or --matching-dir")
    return candidate


def run_one(model: str, args: argparse.Namespace) -> dict[str, Any]:
    """Execute UV, EFT1, conditional direct C5, and SM running in one pass."""
    from Numerical.orchestration.GeneralRunningBenchmark import (
        _state_values, initial_state, load_canonical_general_exports,
    )
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
    from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state

    folder = args.output_root / tag(model)
    folder.mkdir(parents=True, exist_ok=True)
    row: dict[str, Any] = {"model": model, "output_folder": str(folder),
                           "rge_status": "NotRun", "c5_status": "NotRun"}
    stage = "load_general_symbolic_RGE"
    try:
        uv_data, eft_data, provenance = load_canonical_general_exports(
            model, args.canonical_run_id, generate=args.generate_symbolic,
            wolframscript=args.wolframscript,
        )
        benchmark = (illustrative_overrides(uv_data, args.mu_f, args.mu_s)
                     if args.benchmark == "illustrative-real" else {})
        benchmark.update(user_overrides(args.initial_json, model))
        stage = "UV_initial_state"
        _, initial = initial_state(uv_data, mu_uv=args.mu_uv,
                                   mu_f=args.mu_f, mu_s=args.mu_s,
                                   overrides=benchmark)
        stage = "UV_numerical_RGE"
        uv = run_uv_segment(initial, uv_data, args.mu_f,
                            save_scales_gev=np.geomspace(args.mu_uv, args.mu_f, args.uv_points),
                            rtol=args.rtol, atol=args.atol)
        stage = "fermion_boundary_renormalisable_only"
        boundary = build_intermediate_scalar_boundary(uv.final_state)
        _, layout = pack_intermediate_scalar_state(boundary)
        beta_names = set(eft_data.get("report_betas", {}))
        if beta_names != set(layout.names):
            raise ValueError("EFT1 RGE/state mismatch: missing="
                             + repr(sorted(beta_names - set(layout.names)))
                             + ", extra=" + repr(sorted(set(layout.names) - beta_names)))
        stage = "EFT1_numerical_renormalisable_RGE"
        eft = run_intermediate_scalar_segment(
            boundary, eft_data, args.mu_s,
            save_scales_gev=np.geomspace(args.mu_f, args.mu_s, args.eft1_points),
            rtol=args.rtol, atol=args.atol)
        diag = hierarchy_diagnostics(uv.final_state, eft.final_state,
                                     args.mu_f, args.mu_s)
        running = {
            "model": model,
            "status": "NumericalRenormalisableRunningOnly",
            "general_matching_performed": False,
            "scales_gev": {"UV": args.mu_uv, "F": args.mu_f, "S": args.mu_s},
            "uv_solver": {"success": uv.solver_success, "nfev": uv.nfev},
            "eft1_solver": {"success": eft.solver_success, "nfev": eft.nfev},
            "uv_initial": _state_values(uv.initial_state),
            "uv_at_fermion_threshold": _state_values(uv.final_state),
            "eft1_at_fermion_threshold": _state_values(boundary),
            "eft1_at_scalar_threshold": _state_values(eft.final_state),
            "symbolic_provenance": provenance,
            "physical_mass_parameters_available": False,
            "basis_mapping_required_for_physical_interpretation": True,
            "illustrative_benchmark": args.benchmark == "illustrative-real",
            "input_override_keys": sorted(benchmark),
            "hierarchy_diagnostics": diag,
            "eft1_wilson_running_performed": False,
        }
        json_save(folder / "running.json", running)
        row["uv_trajectory_entries"] = save_samples(folder / "uv_running.csv", uv, _state_values)
        row["eft1_trajectory_entries"] = save_samples(folder / "eft1_running.csv", eft, _state_values)
        row["rge_status"] = "Success"
        row["rge_output"] = str(folder / "running.json")
        row["uv_coupling_names"] = len(uv_data["report_betas"])
        row["eft1_coupling_names"] = len(eft_data["report_betas"])
        row["hierarchy"] = diag
        json_save(folder / "inventory.json", {
            "model": model, "format": FORMAT,
            "uv_coupling_keys": sorted(uv_data["report_betas"]),
            "eft1_coupling_keys": sorted(eft_data["report_betas"]),
            "source_provenance": provenance,
            "benchmark": args.benchmark,
            "override_keys": sorted(benchmark),
            "physical_Matchete_RGBeta_parameter_map_certified": False,
            "general_symmetry_basis_certification": False,
        })
    except Exception as exc:
        row["rge_status"] = "Failed"
        row["c5_status"] = "NotRun_RGEFailed"
        row["failure_stage"] = stage
        row["error"] = f"{type(exc).__name__}: {exc}"
        row["traceback"] = traceback.format_exc(limit=7)
        json_save(folder / "status.json", row)
        return row

    if args.no_c5:
        row["c5_status"] = "NotRequested"
        json_save(folder / "status.json", row)
        return row
    stage = "Matchete_general_one_scale_matching"
    try:
        from Numerical.matching.GeneralT3FlavourC5 import GeneralT3FlavourC5
        from scripts.RunT3GeneralSymbolicStudy import expected_yukawa
        from Numerical.running.ScalarThresholdBoundary import (
            build_final_sm_boundary, build_sm_weinberg_initial_conditions,
        )
        from Numerical.running.SMWeinbergEvolution import evolve_sm_weinberg
        from physics.NeutrinoMass import numerical_neutrino_mass_matrix
        from physics.NeutrinoTrajectory import charged_lepton_mass_basis_matrix
        from physics.NeutrinoObservables import calculate_neutrino_observables

        if model in MIXED:
            from Numerical.matching.GeneralT3MixedFlavourC5 import GeneralT3MixedFlavourC5
            matcher = GeneralT3MixedFlavourC5(model, source_dir=args.mixed_source_dir)
            coeff = matcher.coefficient_path
        else:
            coeff = coefficient_path(model, args.matching_dir,
                                     args.generate_matching, args.wolframscript)
            matcher = GeneralT3FlavourC5(model, coeff, wolframscript=args.wolframscript)
            expected = expected_yukawa(model)
            if expected and expected not in matcher.required_yukawas:
                raise RuntimeError(f"{model}: symbolic C5 omits allowed conjugate Yukawa {expected}")
        c5 = np.asarray(matcher.evaluate(uv.final_state, eft.final_state), dtype=complex)
        if c5.shape != (3, 3) or not np.isfinite(c5).all():
            raise ValueError("Matchete result is not a finite 3x3 C5")
        if not np.allclose(c5, c5.T, rtol=1e-10, atol=1e-25):
            raise ValueError("Matchete C5 fails Majorana symmetry")
        sm_boundary = build_final_sm_boundary(eft.final_state)
        initial_sm = build_sm_weinberg_initial_conditions(sm_boundary, c5)
        matched: dict[str, Any] = {
            "model": model, "status": C5_STATUS,
            "matching_kind": matcher.last_diagnostics.get(
                "matching_kind", "general-Matchete-one-scale-flavour-C5"),
            "running_input": str(folder / "running.json"),
            "coefficient": str(coeff),
            "coefficient_sha256": source_hash(coeff),
            "hbar_value": float(matcher.hbar_value),
            "C5_at_scalar_threshold": complex_matrix(c5),
            "c5_norm_at_scalar_gev_inv": float(np.linalg.norm(c5)),
            "general_c5_diagnostics": matcher.last_diagnostics,
            "coupling_mapping_assumptions": matcher.last_diagnostics.get(
                "quartic_basis_mapping_assumptions", []),
            "physical_general_c5_certified": False,
            "sequential_eft_wilson_matching_performed": False,
            "eft1_wilson_running_performed": False,
            "scalar_threshold_sequential_matching_performed": False,
            "low_energy_rge_status": "NotRun",
            "low_energy_rge_errors": [],
        }
        json_save(folder / "c5_conditional.json", matched)
        row["c5_status"] = "ThresholdC5Available"
        row["c5_output"] = str(folder / "c5_conditional.json")

        stage = "SM_Weinberg_low_energy_RGE"
        try:
            evolved = evolve_sm_weinberg(initial_sm, args.mu_s, args.mu_low,
                                         rtol=args.sm_rtol, atol=args.sm_atol)
            mass_matrix = numerical_neutrino_mass_matrix(evolved.K, vev_gev=args.vev)
            lepton_basis_mass = charged_lepton_mass_basis_matrix(mass_matrix, evolved.ye)
            obs = calculate_neutrino_observables(lepton_basis_mass,
                                                 ordering=args.ordering)
            if obs.takagi_residual > 1e-7:
                raise RuntimeError("Unacceptable Takagi residual: " + str(obs.takagi_residual))
            matched.update({
                "low_energy_rge_status": "Success",
                "low_energy_rge_errors": [],
                "observables_scale_gev": float(args.mu_low),
                "C5_at_low_scale": complex_matrix(evolved.K),
                "neutrino_masses_ev": obs.masses_ev.tolist(),
                "delta_m21_sq_ev2": obs.delta_m21_sq_ev2,
                "delta_m31_sq_ev2": obs.delta_m31_sq_ev2,
                "delta_m32_sq_ev2": obs.delta_m32_sq_ev2,
                "pmns_abs": obs.pmns_abs.tolist(),
                "takagi_residual": obs.takagi_residual,
                "charged_lepton_basis_used": True,
            })
            json_save(folder / "c5_conditional.json", matched)
            row["c5_status"] = "ConditionalC5AndSMRunningSuccess"
            row["neutrino_masses_ev"] = obs.masses_ev.tolist()
        except (ValueError, RuntimeError, FloatingPointError, OverflowError) as exc:
            matched.update({"low_energy_rge_status": "LowEnergyRGEFailed_ThresholdOnly",
                            "low_energy_rge_errors": [f"{type(exc).__name__}: {exc}"],
                            "observables_scale_gev": float(args.mu_s),
                            "C5_at_low_scale": None})
            json_save(folder / "c5_conditional.json", matched)
            row["c5_status"] = "ThresholdC5Available_LowRunningFailed"
            row["low_energy_error"] = str(exc)

        if row["c5_status"] == "ConditionalC5AndSMRunningSuccess" and not args.no_plots:
            stage = "sample_C5_trajectory_and_plots"
            from scripts.PlotGeneralT3C5Running import run as plot_running
            plot_args = argparse.Namespace(
                running_json=folder / "running.json",
                c5_json=folder / "c5_conditional.json",
                output_dir=folder / "c5_trajectory",
                points=args.sm_points, mu_low=args.mu_low,
                vev=args.vev, rtol=args.sm_rtol, atol=args.sm_atol,
            )
            plot_summary = plot_running(plot_args)
            row["c5_trajectory_csv"] = str(folder / "c5_trajectory" / "c5_running.csv")
            row["c5_endpoint_relative_check"] = plot_summary.get(
                "sampled_vs_public_endpoint_rel")
    except Exception as exc:
        row["failure_stage"] = stage
        row["c5_error"] = f"{type(exc).__name__}: {exc}"
        row["traceback"] = traceback.format_exc(limit=9)
        if row["c5_status"] in ("NotRun", "ThresholdC5Available"):
            row["c5_status"] = "C5StageFailed"
        elif stage == "sample_C5_trajectory_and_plots":
            row["c5_status"] = "C5AndSMRunningSuccess_PlotFailed"
    json_save(folder / "status.json", row)
    return row


def select_models(args: argparse.Namespace) -> tuple[str, ...]:
    if args.model:
        result = (args.model,)
    elif args.models:
        result = tuple(dict.fromkeys(x.strip() for x in args.models.split(",") if x.strip()))
    else:
        result = tuple(MODELS)
    bad = set(result) - set(MODELS)
    if not result or bad:
        raise ValueError("Invalid model selection; unknown model(s): " + repr(sorted(bad)))
    return result


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    select = p.add_mutually_exclusive_group()
    select.add_argument("--model", choices=MODELS)
    select.add_argument("--models", help="Comma-separated model names (default all 16)")
    select.add_argument("--all", action="store_true", help="Explicitly select all 16")
    p.add_argument("--list-models", action="store_true")
    p.add_argument("--canonical-run-id", default="full16_general_v1")
    p.add_argument("--generate-symbolic", action="store_true")
    p.add_argument("--generate-matching", action="store_true")
    p.add_argument("--matching-dir", type=Path, default=ROOT / "output/general_matchete_matching")
    p.add_argument("--mixed-source-dir", type=Path,
                   default=ROOT / "output/general_t3_mass_eigen_matchete/one_loop_probe",
                   help="Archived B(-1)/C(-1) mass-eigenstate Matchete C5 sources")
    p.add_argument("--output-root", type=Path,
                   default=ROOT / "output/general_observable/all16_v1")
    p.add_argument("--initial-json", type=Path,
                   help="Plain shared overrides, or {shared:{},models:{MODEL:{}}}; values override benchmark")
    p.add_argument("--benchmark", choices=("illustrative-real", "raw"),
                   default="illustrative-real",
                   help="Explicit CP-real illustrative inputs, or original runner baseline")
    p.add_argument("--mu-uv", type=float, default=1e7)
    p.add_argument("--mu-f", type=float, default=1e5)
    p.add_argument("--mu-s", type=float, default=1e3)
    p.add_argument("--mu-low", type=float, default=91.1876)
    p.add_argument("--uv-points", type=int, default=13)
    p.add_argument("--eft1-points", type=int, default=13)
    p.add_argument("--sm-points", type=int, default=41)
    p.add_argument("--rtol", type=float, default=1e-7)
    p.add_argument("--atol", type=float, default=1e-10)
    p.add_argument("--sm-rtol", type=float, default=1e-9)
    p.add_argument("--sm-atol", type=float, default=1e-21)
    p.add_argument("--vev", type=float, default=246.22)
    p.add_argument("--ordering", choices=("AUTO", "NO", "IO"), default="AUTO")
    p.add_argument("--wolframscript", default="wolframscript")
    p.add_argument("--no-c5", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    p = build_parser()
    args = p.parse_args(argv)
    if args.list_models:
        print(" ".join(MODELS))
        return 0
    if not (np.isfinite([args.mu_uv, args.mu_f, args.mu_s, args.mu_low]).all()
            and args.mu_uv > args.mu_f > args.mu_s > args.mu_low > 0):
        p.error("Require finite mu_UV > mu_F > mu_S > mu_low > 0")
    if min(args.uv_points, args.eft1_points, args.sm_points) < 3:
        p.error("Each trajectory must have at least 3 sample points")
    if min(args.rtol, args.atol, args.sm_rtol, args.sm_atol) <= 0:
        p.error("ODE tolerances must be positive")
    try:
        models = select_models(args)
    except ValueError as exc:
        p.error(str(exc))
    args.output_root.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, model in enumerate(models, start=1):
        print(f"[{i}/{len(models)}] {model}: running general UV/EFT1", flush=True)
        result = run_one(model, args)
        rows.append(result)
        print(f"    RGE: {result['rge_status']}; C5: {result['c5_status']}", flush=True)
        if result.get("error") or result.get("c5_error"):
            print("    " + str(result.get("error", result.get("c5_error"))), flush=True)
        json_save(args.output_root / "summary.json", {
            "format": FORMAT,
            "canonical_run_id": args.canonical_run_id,
            "benchmark": args.benchmark,
            "scales_gev": {"UV": args.mu_uv, "F": args.mu_f,
                           "S": args.mu_s, "low": args.mu_low},
            "selected_models": models,
            "rge_success_count": sum(x["rge_status"] == "Success" for x in rows),
            "conditional_c5_low_success_count": sum(
                x["c5_status"] in ("ConditionalC5AndSMRunningSuccess",
                                    "C5AndSMRunningSuccess_PlotFailed") for x in rows),
            "mixed_model_conditional_c5_count": sum(
                x["model"] in MIXED and x["c5_status"] in
                ("ConditionalC5AndSMRunningSuccess", "C5AndSMRunningSuccess_PlotFailed")
                for x in rows),
            "physically_certified_sequential_eft_count": 0,
            "matchete_rgbeta_coupling_dictionary_certified": False,
            "results": rows,
        })
    count = sum(x["rge_status"] == "Success" for x in rows)
    c5count = sum(x["c5_status"] in ("ConditionalC5AndSMRunningSuccess",
                                         "C5AndSMRunningSuccess_PlotFailed") for x in rows)
    print(f"General UV/EFT1: {count}/{len(models)}. Conditional C5+SM: "
          f"{c5count}/{len(models)}. Sequential EFT: 0 certified.", flush=True)
    # B/C(-1) missing generic C5 adapter is a known scoped limitation, not an RGE failure.
    failures = [x for x in rows if (x["rge_status"] != "Success" or
        (not args.no_c5 and x["model"] not in MIXED and
         x["c5_status"] != "ConditionalC5AndSMRunningSuccess"))]
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
