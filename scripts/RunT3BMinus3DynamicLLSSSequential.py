"""B(-3): dynamic F-tree LLSS -> leading-log C5 -> SMEFT numerical benchmark.

Stage 3. Each parameter point is freshly run through the existing full-general
UV and scalar-only EFT1 RGBeta ODEs. The actual evolved UV fermion masses and
Yukawas generate a NEW numerical Matchete-tree LLSS tensor at muF. The saved
source-backed HHSS scalar channels are evaluated with the NEW EFT1(muS) state,
and their one-loop direct LLSS -> Weinberg leading log is computed. This
REPLACES, and is never added on top of, the analytic LLSS direct log in the
existing B(-3) sequential hard+running kernel. The saved MSbar scalar-threshold
hard remainder is retained. A separate scalar bubble is exported for diagnosis
but NOT added a second time.

This is a conditional fixed-one-loop process, NOT a certified all-operator
Wilson evolution: the Matchete GammaCC-to-numerical-Weyl normalization is a
user-declared factor, the native extra-quartic basis map is provisional, and
the direct leading log uses muS scalar couplings (no C_LLSS self-transport).
Nothing is changed in the existing one-scale fitter or symbolic matching.

PowerShell from repository root:
  python -X utf8 -m scripts.RunT3BMinus3DynamicLLSSSequential --check-only --beta-source archived
  python -X utf8 -m scripts.RunT3BMinus3DynamicLLSSSequential --run-id llss_baseline_v1 --beta-source archived
  python -X utf8 -m scripts.RunT3BMinus3DynamicLLSSSequential --run-id llss_y2_v1 --beta-source archived --yukawa-scale 2

Inputs: existing B(-3) Stage 1/2 scripts, saved Matchete sequential matching,
original F-first Wilson exports, quartic CG cache, complete scalar registry.
Outputs: output/general_sequential_pilot/B_m3_dynamic_llss_sequential/<run-id>/
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Any

import numpy as np

from scripts import RunT3BMinus3DynamicSequential as stage2
from scripts import RunT3BMinus3SequentialSMEFT as stage1

ROOT = Path(__file__).resolve().parents[1]
MODEL = "B_-3"
DEFAULT_OUT = ROOT / "output/general_sequential_pilot/B_m3_dynamic_llss_sequential"
FORMAT = "t3-b-minus3-dynamic-f-tree-llss-one-loop-sequential-smeft-v1"


def _path(x: Path) -> Path:
    return stage2.local(Path(x))


def _source_files() -> dict[str, Path]:
    from scripts.RunT3GeneralSequentialLeadingLog import (
        DEFAULT_RELEASE, DEFAULT_V4, DEFAULT_SEEDS, DEFAULT_REGISTRY,
        source_paths,
    )
    return source_paths(
        dict(release=DEFAULT_RELEASE, v4=DEFAULT_V4, seeds=DEFAULT_SEEDS,
             running=ROOT / "output/general_observable/all16_v1",
             registry=DEFAULT_REGISTRY), MODEL,
    )


def preflight_scalar_sources() -> tuple[dict, Any, Any, dict]:
    """Load original, hash-bound Matchete CG quartics; NOT a frozen LLSS input."""
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import read_v4_cached_channels
    paths = _source_files()
    required = ("boundary", "channels", "quartic_seed", "registry")
    missing = {k: str(paths[k]) for k in required if not paths[k].is_file()}
    if missing:
        raise FileNotFoundError(f"B(-3) LLSS/quartic source files absent: {missing}")
    reference = load_source(paths["boundary"], expected_model=MODEL)
    quartics = read_v4_cached_channels(paths["channels"],
                                       paths["quartic_seed"], reference)
    registry = stage1.read_object(paths["registry"])
    if MODEL not in registry.get("models", {}):
        raise ValueError("Full-general scalar registry missing B(-3)")
    if reference.model != MODEL or quartics.model != MODEL:
        raise ValueError("B(-3) scalar/fermion source model mismatch")
    return paths, reference, quartics, registry


def dynamic_llss_response(running: dict, config: dict, reference: Any,
                          quartics: Any, registry: dict, *, factor: float,
                          check_scale_cancellation: bool) -> tuple[Any, dict, dict]:
    """Fitted muF Matchete tree -> muS HHSS tensors -> direct one-loop C5.

    Returns the dynamic SourceBoundary, a decoded numerical LLSS result, and a
    provenance record. The original static reference is used to verify the
    symbolic quartic cache only, NOT as the numerical fermion boundary.
    """
    from Numerical.running.GeneralT3DynamicFermionWilson import nonmixed_boundary
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import RunningSnapshot
    from Numerical.running.GeneralT3SequentialLeadingLog import scalar_tensors, response, decode

    scales = config["scales"]
    mu_f = float(scales["mu_fermion_threshold_gev"])
    mu_s = float(scales["mu_scalar_threshold_gev"])
    uv = running["uv_at_fermion_threshold"]
    eft = running["eft1_at_scalar_threshold"]
    heavy = stage1.fitted_fermion_state(uv, config, mu_f)
    dynamic, details = nonmixed_boundary(MODEL, heavy, repo_root=ROOT, mu_f=mu_f)
    if dynamic.model != reference.model or tuple(dynamic.labels) != tuple(reference.labels):
        raise ValueError("Dynamic Matchete F-tree and cached quartic scalar bases differ")
    if dynamic.raw.shape != (3, 3, 2, 2, 8, 8):
        raise ValueError(f"B(-3) dynamic LLSS dimensions are wrong: {dynamic.raw.shape}")
    if not np.isfinite(dynamic.raw).all():
        raise ValueError("Dynamic F-tree LLSS tensor has nonfinite entries")

    snapshot = RunningSnapshot(
        model=MODEL,
        source=Path("dynamic_full_general_eft1_endpoint"),
        sha256=sha256(json.dumps(eft, sort_keys=True).encode()).hexdigest(),
        mu_f=mu_f, mu_s=mu_s,
        scalar_mass1_sq=stage1.finite_real(eft["mS1Sq"], "mS1Sq"),
        scalar_mass2_sq=stage1.finite_real(eft["mS2Sq"], "mS2Sq"),
        couplings=eft,
        evidence={"fresh_full_general_EFT1": True},
    )
    bundle = scalar_tensors(
        MODEL, dynamic, quartics, snapshot, registry,
        include_general_extra=True, include_portal_completion=True,
    )
    raw = response(
        dynamic, bundle, mu_f=mu_f, mu_s=mu_s,
        spinor_factor=complex(factor),
        check_scale_cancellation=check_scale_cancellation,
    )
    direct = decode(raw["c5_eft1_leading_log"])
    bubble = decode(raw["c5_scalar_finite_bubble"])
    if direct.shape != (3, 3) or bubble.shape != (3, 3):
        raise ValueError("LLSS scalar matching did not produce a 3x3 C5 matrix")
    details_out = {
        "source": "dynamic_current_uv_matchete_tree",
        "computed_from_current_uv_at_muF": True,
        "F_tree_source_diagnostics": details,
        "F_tree_source_hash": dynamic.hash,
        "F_tree_shape": list(dynamic.raw.shape),
        "F_tree_frobenius_norm_gev_inverse": float(np.linalg.norm(dynamic.raw)),
        "HHSS_model": bundle.model,
        "HHSS_channels": list(bundle.names),
        "HHSS_channel_count": len(bundle.names),
        "HHSS_provenance": bundle.provenance,
        "HHSS_quartic_and_mass_parameters_from_new_eft1_at_muS": True,
        "fermion_first_tree_passed_to_scalar_response": True,
        "wilson_normalization_factor_assumed": float(factor),
        "direct_C5_leading_log_gev_inverse": stage1.encode_matrix(direct),
        "separately_calculated_scalar_bubble_gev_inverse": stage1.encode_matrix(bubble),
        "scalar_bubble_added_to_final_C5": False,
        "scalar_bubble_nonaddition_reason": (
            "Saved renormalised stage-2 Matchete hard remainder is retained; "
            "its scalar finite contribution is not independently subtracted, "
            "so adding the new bubble would risk double counting"
        ),
        "matching_scale_cancellation_residual": raw["matching_scale_cancellation_residual"],
        "weinberg_beta_projection_residual": raw["weinberg_beta_projection_residual"],
        "weinberg_bubble_projection_residual": raw["weinberg_bubble_projection_residual"],
        "additional_scalar_quartics_from_current_endpoint": True,
        "full_Wilson_self_transport_performed": False,
        "intermediate_direct_Wilson_mixing_computed": True,
        "one_loop_fixed_order_approximation": True,
    }
    return dynamic, raw, details_out


def _numeric_sympy_3x3(symbolic: Any, subs: dict, name: str) -> np.ndarray:
    import sympy as sp
    if symbolic.shape != (3, 3):
        raise ValueError(f"{name} must be 3x3")
    result = np.empty((3, 3), dtype=complex)
    for i in range(3):
        for j in range(3):
            expr = sp.N(symbolic[i, j].subs(subs), 18)
            if getattr(expr, "free_symbols", None):
                raise ValueError(f"{name}[{i},{j}] has unbound parameters: "
                                 f"{sorted(map(str, expr.free_symbols))}")
            result[i, j] = complex(expr)
    if not np.isfinite(result).all():
        raise ValueError(f"{name} contains nonfinite values")
    if not np.allclose(result, result.T, rtol=1e-9, atol=1e-25):
        raise ValueError(f"{name} is not a symmetric Majorana coefficient")
    return (result + result.T) / 2


def saved_threshold_components(c5_path: Path, running: dict,
                               config: dict) -> tuple[np.ndarray, np.ndarray, dict]:
    """Evaluate saved stage-2 renormalised remainder and original direct log.

    The two terms are separable in the existing FinalWeinbergAdapter. We retain
    the first and REPLACE the second with the dynamic LLSS response.
    """
    from Numerical.running.FinalC5TrajectoryAdapter import (
        FinalC5EvaluationInputs, _heavy_mass_basis_data,
    )
    from Numerical.running.SMWeinbergStage import (
        _sympy_substitutions_from_config, _hierarchical_scale_substitutions,
    )
    from RGE.matching.FinalWeinbergAdapter import load_hierarchical_majorana_c5

    mu_f = float(config["scales"]["mu_fermion_threshold_gev"])
    mu_s = float(config["scales"]["mu_scalar_threshold_gev"])
    eft = running["eft1_at_scalar_threshold"]
    heavy_state = stage1.fitted_fermion_state(
        running["uv_at_fermion_threshold"], config, mu_f
    )
    heavy = _heavy_mass_basis_data(
        heavy_state,
        offdiagonal_atol=1e-7, imaginary_atol=1e-9,
        takagi_residual_rtol=1e-8, takagi_residual_atol=1e-7,
    )
    mass1 = stage1.finite_real(eft["mS1Sq"], "mS1Sq")
    mass2 = stage1.finite_real(eft["mS2Sq"], "mS2Sq")
    if mass1 <= 0 or mass2 <= 0:
        raise ValueError("Cannot integrate out tachyonic scalar multiplet")
    values = FinalC5EvaluationInputs(
        heavy_masses_gev=np.asarray(heavy.masses_gev, float),
        y1=np.asarray(heavy.y1, complex),
        y2=np.asarray(heavy.y2, complex),
        scalar_mass1_gev=math.sqrt(mass1),
        scalar_mass2_gev=math.sqrt(mass2),
        grouped_scalar_matching_scale_gev=mu_s,
        lambda_t3=complex(stage1.decode(eft["lambdaT3"])),
    )
    model_config = values.as_config()
    symbolic = load_hierarchical_majorana_c5(
        c5_path, n_lepton=3,
        n_heavy=len(heavy.masses_gev), split_heavy_masses=True,
    )
    subs = _sympy_substitutions_from_config(model_config)
    subs.update(_hierarchical_scale_substitutions(model_config))
    hard = _numeric_sympy_3x3(symbolic["K_hard"], subs,
                              "saved_scalar_threshold_hard")
    old_log = _numeric_sympy_3x3(symbolic["K_running"], subs,
                                 "saved_symbolic_LLSS_direct_log")
    return hard, old_log, {
        "heavy_mass_basis_residual": float(max(heavy.takagi_residual, heavy.basis_residual)),
        "fermion_masses_gev": np.asarray(heavy.masses_gev, float).tolist(),
        "scalar_masses_gev": [math.sqrt(mass1), math.sqrt(mass2)],
        "lambdaT3_at_muS": {"real": float(values.lambda_t3.real),
                           "imag": float(values.lambda_t3.imag)},
        "hard_remainder_is_preexisting_symbolic_Matchete_stage2": True,
        "hard_remainder_recomputed_from_current_threshold_parameters": True,
    }


def smeft_observables(c5_at_s: np.ndarray, running: dict, config: dict,
                      target_path: Path) -> dict:
    """Run existing SMEFT RGEs and oscillation observable stack for new C5."""
    from Numerical.running.SMWeinbergEvolution import (
        SMWeinbergInitialConditions, evolve_sm_weinberg,
    )
    from Numerical.fitting.OscillationFit import (
        OscillationFitTarget, evaluate_mass_matrix_fit,
    )
    from physics.NeutrinoTrajectory import (
        charged_lepton_mass_basis_matrix, neutrino_mass_from_c5,
    )
    from physics.NeutrinoObservables import calculate_neutrino_observables

    eft = running["eft1_at_scalar_threshold"]
    numerical = config.get("numerical", {})
    mu_s = float(config["scales"]["mu_scalar_threshold_gev"])
    mu_low = float(config["scales"]["mu_low_gev"])
    vev = float(numerical.get("vev_gev", 246.22))
    initial = SMWeinbergInitialConditions(
        gY=stage1.finite_real(eft["gY"], "gY"),
        g2=stage1.finite_real(eft["g2"], "g2"),
        g3=stage1.finite_real(eft["g3"], "g3"),
        lambdaH=stage1.finite_real(eft["lambdaH"], "lambdaH"),
        yu=stage1.finite_matrix(eft["yu"], "yu"),
        yd=stage1.finite_matrix(eft["yd"], "yd"),
        ye=stage1.finite_matrix(eft["ye"], "ye"),
        K=c5_at_s,
    ).validated()
    end = evolve_sm_weinberg(
        initial, mu_s, mu_low,
        rtol=float(numerical.get("rtol", 1e-8)),
        atol=float(numerical.get("atol", 1e-11)),
    )
    if not end.solver_success:
        raise RuntimeError("SMEFT integration reported failure")
    mnu = neutrino_mass_from_c5(end.K, vev_gev=vev)
    lepton_mass = charged_lepton_mass_basis_matrix(mnu, end.ye)
    target = OscillationFitTarget.from_json(target_path)
    fit = evaluate_mass_matrix_fit(lepton_mass, target)
    obs = calculate_neutrino_observables(lepton_mass, ordering=target.ordering)
    return {
        "smeft": {
            "RGEs_evolved": True,
            "solver_success": True,
            "solver_nfev": int(end.nfev),
            "C5_at_muS_gev_inverse": stage1.encode_matrix(c5_at_s),
            "C5_at_muLow_gev_inverse": stage1.encode_matrix(end.K),
        },
        "neutrinos": {
            "ordering": obs.ordering,
            "vev_gev": vev,
            "masses_ev_m1_m2_m3": [float(x) for x in obs.masses_ev],
            "mass_matrix_gev": stage1.encode_matrix(mnu),
            "mass_matrix_charged_lepton_basis_gev": stage1.encode_matrix(lepton_mass),
            "PMNS_absolute": np.asarray(obs.pmns_abs).tolist(),
            "takagi_residual": float(obs.takagi_residual),
            "oscillation_observables": {
                name: float(value) for name, value
                in zip(fit.observable_names, fit.prediction)
            },
            "oscillation_chi2": float(fit.chi2),
            "target_ordering": str(target.ordering),
        },
    }


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", default="llss_baseline_v1")
    ap.add_argument("--check-only", action="store_true")
    ap.add_argument("--canonical-run-id", default="full16_general_v1")
    ap.add_argument("--beta-source", choices=("auto", "canonical", "archived"),
                    default="archived")
    ap.add_argument("--matching-json", type=Path, default=stage2.DEFAULT_C5)
    ap.add_argument("--fit-root", type=Path, default=stage2.DEFAULT_FIT)
    ap.add_argument("--param-overrides", type=Path)
    ap.add_argument("--yukawa-scale", type=float, default=1.0)
    ap.add_argument("--target", type=Path)
    ap.add_argument("--wilson-factor", type=float, default=1.0,
                    help="Declared Matchete GammaCC / numerical Weyl normalization, default 1")
    ap.add_argument("--check-scale-cancellation", action="store_true",
                    help="Optional conditional one-loop matching scale cancellation check")
    ap.add_argument("--output-root", type=Path, default=DEFAULT_OUT)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", args.run_id):
        raise ValueError("--run-id must be a single folder name")
    if not math.isfinite(args.wilson_factor) or args.wilson_factor == 0:
        raise ValueError("--wilson-factor must be finite and nonzero")

    c5_path = _path(args.matching_json)
    fit_root = _path(args.fit_root)
    fit_path = fit_root / "result.json"
    config_path = fit_root / "fitted_config.json"
    overrides_path = fit_root / "general_uv_overrides.json"
    c5_obj, fit, config, pilot = stage1.source_preflight(
        c5_path, fit_path, config_path
    )
    base = stage1.read_object(overrides_path)
    changes = stage1.read_object(_path(args.param_overrides)) if args.param_overrides else None
    overrides = stage2.build_uv_overrides(base, changes, args.yukawa_scale)
    uv_beta, eft_beta, provenance = stage2.load_beta_pair(
        args.canonical_run_id, args.beta_source
    )
    if set(overrides).difference(uv_beta.get("report_betas", {})):
        raise ValueError("UV parameters not all present in RGBeta beta inventory")
    paths, reference, quartics, registry = preflight_scalar_sources()
    if args.check_only:
        print("B(-3) Matchete sequential hard kernel: OK")
        print("B(-3) source-verified F-tree reference and HHSS channels: OK")
        print(f"Source channel names: {len(quartics.by_name)}")
        print("RGBeta UV/EFT1 beta pair: OK (" + str(provenance.get("source")) + ")")
        print("NOTICE: full LLSS gamma and spinor normalization are not certified")
        print("No integration or files created by --check-only")
        return 0

    output = _path(args.output_root) / args.run_id
    if output.exists():
        raise FileExistsError(f"Output folder already exists: {output}")
    if args.target is None:
        from Numerical.orchestration.RunSmallLambdaRefit import TARGET
        target = Path(TARGET)
    else:
        target = _path(args.target)
    if not target.is_file():
        raise FileNotFoundError(f"Missing target oscillation data: {target}")

    running = stage2.run_dynamic_rges(uv_beta, eft_beta, config, overrides)
    boundary, response_data, llss_detail = dynamic_llss_response(
        running, config, reference, quartics, registry,
        factor=args.wilson_factor,
        check_scale_cancellation=args.check_scale_cancellation,
    )
    hard, old_log, hard_detail = saved_threshold_components(
        c5_path, running, config
    )
    from Numerical.running.GeneralT3SequentialLeadingLog import decode
    new_log = decode(response_data["c5_eft1_leading_log"])
    c5_at_s = (hard + new_log)
    c5_at_s = 0.5 * (c5_at_s + c5_at_s.T)
    if not np.isfinite(c5_at_s).all():
        raise ValueError("Nonfinite dynamically matched C5 at scalar threshold")
    physical = smeft_observables(c5_at_s, running, config, target)
    baseline = changes is None and args.yukawa_scale == 1.0
    result = {
        "format": FORMAT,
        "model": MODEL,
        "status": "DynamicFTreeLLSSLeadingLogAndSavedScalarHardToSMEFTComputed",
        "scales_gev": {"UV": config["base_state"]["mu_gev"],
                       **running["scales_gev"],
                       "low": config["scales"]["mu_low_gev"]},
        "renormalisable_running": {
            "fresh_uv_solver": running["uv_solver"],
            "fresh_eft1_solver": running["eft1_solver"],
            "uv_initial": running["uv_initial"],
            "UV_at_fermion_threshold": running["uv_at_fermion_threshold"],
            "EFT1_at_fermion_threshold": running["eft1_at_fermion_threshold"],
            "EFT1_at_scalar_threshold": running["eft1_at_scalar_threshold"],
        },
        "fermion_first_LLSS": llss_detail,
        "scalar_threshold_saved_Matchete": hard_detail,
        "matching_decomposition_gev_inverse": {
            "hard_remainder_from_saved_Matchete_scalar_threshold": stage1.encode_matrix(hard),
            "original_symbolic_direct_log_replaced_not_added": stage1.encode_matrix(old_log),
            "dynamic_LLSS_direct_log_used": stage1.encode_matrix(new_log),
            "independent_dynamic_scalar_bubble_not_added_twice": (
                response_data["c5_scalar_finite_bubble"]
            ),
            "C5_current_final_at_muS": stage1.encode_matrix(c5_at_s),
            "C5_original_symbolic_at_current_parameters": (
                stage1.encode_matrix(hard + old_log)
            ),
            "change_in_direct_log_frobenius": float(np.linalg.norm(new_log - old_log)),
            "direct_log_original_frobenius": float(np.linalg.norm(old_log)),
            "direct_log_dynamic_frobenius": float(np.linalg.norm(new_log)),
        },
        **physical,
        "sources": {
            "beta_source": provenance,
            "symbolic_final_C5_file": str(c5_path.resolve()),
            "original_fitted_parameter_file": str(overrides_path.resolve()),
            "parameter_overrides_file": (
                str(_path(args.param_overrides).resolve()) if args.param_overrides else None
            ),
            "original_F_tree_reference_file": str(paths["boundary"].resolve()),
            "native_HHSS_cache": str(paths["channels"].resolve()),
            "full_scalar_quartic_registry": str(paths["registry"].resolve()),
            "F_tree_operator_source": str(boundary.path),
            "spinor_factor_assumption": args.wilson_factor,
        },
        "flags": {
            "full_general_UV_RGE_rerun": True,
            "fermion_first_Matchete_tree_evaluated_at_current_muF": True,
            "fermion_removed_from_EFT1_renormalisable_running": True,
            "full_general_EFT1_RGE_rerun": True,
            "dynamic_LLSS_direct_one_loop_mixing_fed_to_C5": True,
            "both_scalars_integrated_at_one_matching_scale": True,
            "saved_Matchete_scalar_threshold_hard_reused": True,
            "SMEFT_Weinberg_RGE_rerun": True,
            "full_LLSS_self_transport_performed": False,
            "hard_matching_recomputed_in_Wolfram_for_each_parameter_point": False,
            "Matchete_Weyl_normalization_certified": False,
            "complete_all_operator_one_loop_prediction_certified": False,
        },
        "scope": (
            "B(-3) freshly run UV and EFT1 renormalisable RGEs; Matchete "
            "F-first tree LLSS assembled from current Yukawas and MF; "
            "muS general HHSS quartics contracted with the LLSS tensor to "
            "produce direct fixed-one-loop EFT1 leading log; saved "
            "renormalised Matchete scalar-threshold hard remainder reused "
            "WITHOUT adding separately calculated scalar bubble twice; "
            "SM+Weinberg running and neutrino observables recomputed. "
            "Partial one-loop Wilson mixing, conditional basis and convention."
        ),
    }
    if baseline:
        saved = fit["numerical_details"]
        result["baseline_reproduction"] = {
            "UV_endpoint_relative": stage2.reproduction_diagnostic(
                running["uv_at_fermion_threshold"],
                saved["uv_at_fermion_threshold"], stage2.REPRO_KEYS_UV,
            ),
            "EFT1_endpoint_relative": stage2.reproduction_diagnostic(
                running["eft1_at_scalar_threshold"],
                saved["eft1_at_scalar_threshold"], stage2.REPRO_KEYS_EFT1,
            ),
        }

    # Full dynamic source and the summed-matrix decomposition are saved for
    # independent inspection; none of the original matching artifacts change.
    stage1.write_outputs(output, result)
    np.savez_compressed(
        output / "fermion_tree_LLSS_at_MF.npz",
        C_LLSS_gev_inverse=boundary.raw,
        C_LLSS_by_source_gev_inverse=boundary.by_source,
        C5_direct_from_dynamic_LLSS_gev_inverse=new_log,
        C5_saved_scalar_hard_gev_inverse=hard,
        C5_combined_gev_inverse=c5_at_s,
        mu_F_gev=np.asarray(float(config["scales"]["mu_fermion_threshold_gev"])),
        mu_S_gev=np.asarray(float(config["scales"]["mu_scalar_threshold_gev"])),
        model=np.asarray(MODEL),
    )
    print("Wrote:", output / "result.json")
    print("Fresh UV/EFT1 RGE evaluations:",
          running["uv_solver"]["nfev"], "/", running["eft1_solver"]["nfev"])
    print("Dynamic F-tree LLSS dimensions:", list(boundary.raw.shape))
    print("Scalar HHSS channels:", llss_detail["HHSS_channel_count"])
    print("Analytic saved direct-log norm:", f"{np.linalg.norm(old_log):.9g}")
    print("Dynamic tensor direct-log norm:", f"{np.linalg.norm(new_log):.9g}")
    print("chi2 =", f"{result['neutrinos']['oscillation_chi2']:.9g}")
    print("masses (eV) =", ", ".join(
        f"{x:.9g}" for x in result["neutrinos"]["masses_ev_m1_m2_m3"]
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
