"""B(-3): dynamic fermion-first LLSS + BOTH-scalar matching + SMEFT.

Stage 4 of the opt-in exploratory B(-3) numerical process.  This driver
REUSES the installed Stage 1, 2, and 3 scripts without changing them.

For each UV parameter point:
  1. run the full-general UV RGBeta RGEs to mu_F;
  2. integrate out all fermions at tree level to numerical C_LLSS(mu_F),
     using the saved source-derived Matchete operator geometry;
  3. remove the heavy fermions from the renormalisable EFT1 state, and run
     the scalar-only EFT1 RGBeta RGEs to mu_S;
  4. obtain the complete B(-3) candidate HHSS tensor at mu_S and calculate
     direct LLSS -> C5 one-loop leading-log mixing;
  5. integrate out ALL eight real scalar components (both scalar doublets)
     in a single one-loop finite scalar bubble at mu_S;
  6. evolve SM + the resulting Weinberg coefficient to mu_low, and evaluate
     neutrino observables/chi-square.

The LO scalar/Wilson contribution is formed ONCE:
    C5_LLSS(mu_S) = Delta C5_direct_LLSS(mu_F -> mu_S)
                    + Delta C5_joint_scalar_bubble(mu_S).
This is the existing joint-scalar matching kernel's c5_LLSS_total.
The saved analytic direct log is NOT added again.

IMPORTANT: The existing saved B(-3) Matchete hard term cannot yet be proven
independent of the numerical scalar bubble in these two conventions. Therefore
BOTH conditional interpretations are computed and published:
    omit   = C5_LLSS only (no saved hard term; default primary result),
    retain = C5_LLSS + saved Matchete hard term (candidate; possible overlap).
The --hard-policy argument chooses which is the primary reported result; it
DOES NOT verify the physical overlap. No arbitrary subtraction or numerical
rescaling is applied. This is not a certified complete one-loop matching.

The Matchete GammaCC/Weyl factor is declared with --wilson-factor (default 1).
The scalar Matchete/RGBeta basis crosswalk is conditional. LLSS self-running
would enter the one-loop scalar bubble at O((16 pi^2)^-2), and is omitted at
this fixed-one-loop order. It is not a full seven-term Wilson evolution.

PowerShell, from repository root:
  python -X utf8 -m scripts.RunT3BMinus3DynamicScalarThreshold --check-only --beta-source archived
  python -X utf8 -m scripts.RunT3BMinus3DynamicScalarThreshold --run-id scalar_baseline_v1 --beta-source archived
  python -X utf8 -m scripts.RunT3BMinus3DynamicScalarThreshold --run-id scalar_y2_v1 --beta-source archived --yukawa-scale 2

Outputs (new only):
  output/general_sequential_pilot/B_m3_dynamic_scalar_threshold/<run-id>/
    result.json, observables.csv, observables_alternative.csv,
    dynamic_matching_C5.npz
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
from typing import Any

import numpy as np

from scripts import RunT3BMinus3SequentialSMEFT as stage1
from scripts import RunT3BMinus3DynamicSequential as stage2
from scripts import RunT3BMinus3DynamicLLSSSequential as stage3

ROOT = Path(__file__).resolve().parents[1]
MODEL = "B_-3"
DEFAULT_OUT = ROOT / "output/general_sequential_pilot/B_m3_dynamic_scalar_threshold"
FORMAT = "t3-b-minus-3-dynamic-all-scalar-threshold-smeft-v1"


def _matrix(raw: Any, label: str) -> np.ndarray:
    """Read an existing encoded complex matrix with strict dimension checks."""
    from Numerical.running.GeneralT3SequentialLeadingLog import decode

    value = np.asarray(decode(raw), dtype=complex)
    if value.shape != (3, 3) or not np.isfinite(value).all():
        raise ValueError(f"{label}: expected finite complex symmetric 3x3 C5")
    # Absolute tolerance is deliberately small relative to the O(1e-15) C5.
    if not np.allclose(value, value.T, rtol=1e-8, atol=1e-28):
        raise ValueError(f"{label}: coefficient is not symmetric in lepton flavour")
    return 0.5 * (value + value.T)


def assemble_matching(response: dict, saved_hard: np.ndarray,
                      saved_analytic_log: np.ndarray) -> dict[str, np.ndarray | float]:
    """One dynamic soft threshold; never add the analytic log or bubble twice."""
    direct = _matrix(response["c5_eft1_leading_log"], "LLSS direct leading log")
    bubble = _matrix(response["c5_scalar_finite_bubble"], "joint scalar bubble")
    total_from_source = _matrix(response["c5_LLSS_total"], "source combined LLSS")
    hard = np.asarray(saved_hard, dtype=complex)
    analytic_log = np.asarray(saved_analytic_log, dtype=complex)
    for label, value in (("saved hard", hard), ("old analytic log", analytic_log)):
        if value.shape != (3, 3) or not np.isfinite(value).all() or not np.allclose(
            value, value.T, rtol=1e-8, atol=1e-28
        ):
            raise ValueError(f"Invalid {label} coefficient")
    total = direct + bubble
    discrepancy = float(np.linalg.norm(total - total_from_source))
    scale = max(float(np.linalg.norm(total_from_source)),
                float(np.linalg.norm(direct)),
                float(np.linalg.norm(bubble)), 1e-30)
    if discrepancy > scale * 1e-9 + 1e-28:
        raise ValueError("Joint-scalar source total disagrees with direct log + bubble")
    return {
        "direct": direct,
        "bubble": bubble,
        "dynamic_llss_only": total,
        "saved_hard": hard,
        "saved_analytic_log_diagnostic": analytic_log,
        "dynamic_llss_plus_saved_hard": total + hard,
        "previous_stage3": hard + direct,
        "identity_absolute_residual": discrepancy,
        "identity_relative_residual": discrepancy / scale,
    }


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-id", default="scalar_baseline_v1")
    p.add_argument("--check-only", action="store_true")
    p.add_argument("--beta-source", default="archived",
                   choices=("auto", "canonical", "archived"))
    p.add_argument("--canonical-run-id", default="full16_general_v1")
    p.add_argument("--matching-json", type=Path, default=stage2.DEFAULT_C5)
    p.add_argument("--fit-root", type=Path, default=stage2.DEFAULT_FIT)
    p.add_argument("--param-overrides", type=Path)
    p.add_argument("--yukawa-scale", type=float, default=1.0)
    p.add_argument("--target", type=Path)
    p.add_argument("--wilson-factor", type=float, default=1.0,
                   help="Assumed GammaCC -> Weyl normalization; no implicit fitting")
    p.add_argument("--hard-policy", choices=("omit", "retain"), default="omit",
                   help="Primary interpretation; BOTH are always calculated")
    p.add_argument("--skip-scale-check", action="store_true",
                   help="Do not run conditional fixed-coupling scale-cancellation check")
    p.add_argument("--output-root", type=Path, default=DEFAULT_OUT)
    return p


def _compact_physics(physical: dict) -> dict:
    """The two SMEFT branches have identical schema, not identical physics."""
    return {
        "smeft": physical["smeft"],
        "neutrinos": physical["neutrinos"],
    }


def _alternative_csv(path: Path, physical: dict) -> None:
    """Write explicit alternative observables, without modifying Stage 1 I/O."""
    n = physical["neutrinos"]
    lines = ["quantity,value,unit"]
    for key, value in n["oscillation_observables"].items():
        unit = "dimensionless" if key.startswith("sin2") else "eV^2"
        lines.append(f"{key},{float(value):.17g},{unit}")
    for index, mass in enumerate(n["masses_ev_m1_m2_m3"], 1):
        lines.append(f"m{index},{float(mass):.17g},eV")
    lines.append(f"chi2,{float(n['oscillation_chi2']):.17g},dimensionless")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", args.run_id) is None:
        raise ValueError("--run-id must be a safe single directory name")
    if not math.isfinite(args.wilson_factor) or args.wilson_factor == 0:
        raise ValueError("--wilson-factor must be finite and nonzero")

    c5_path = stage2.local(args.matching_json)
    fit_root = stage2.local(args.fit_root)
    fit_path = fit_root / "result.json"
    cfg_path = fit_root / "fitted_config.json"
    override_path = fit_root / "general_uv_overrides.json"
    _c5, saved_fit, config, _pilot = stage1.source_preflight(
        c5_path, fit_path, cfg_path
    )
    base = stage1.read_object(override_path)
    changes = (stage1.read_object(stage2.local(args.param_overrides))
               if args.param_overrides is not None else None)
    overrides = stage2.build_uv_overrides(base, changes, args.yukawa_scale)
    uv_beta, eft_beta, beta_provenance = stage2.load_beta_pair(
        args.canonical_run_id, args.beta_source
    )
    missing = set(overrides).difference(uv_beta.get("report_betas", {}))
    if missing:
        raise ValueError(f"UV override has couplings without beta functions: {sorted(missing)}")
    paths, reference, quartics, registry = stage3.preflight_scalar_sources()

    if args.check_only:
        print("B(-3): archived symbolic hard and direct C5 terms: OK")
        print("B(-3): dynamic F-tree Matchete and scalar HHSS sources: OK")
        print("B(-3): UV and EFT1 beta inventories: OK")
        print("Scalar native HHSS source channels:", len(quartics.by_name))
        print("Hard treatment: TWO conditional results; primary =", args.hard_policy)
        print("NOTICE: scalar-vs-saved-hard overlap is not established")
        print("No integrations or files written by --check-only")
        return 0

    out_dir = stage2.local(args.output_root) / args.run_id
    if out_dir.exists():
        raise FileExistsError(f"Output already exists: {out_dir}")
    if args.target is None:
        from Numerical.orchestration.RunSmallLambdaRefit import TARGET
        target = Path(TARGET)
    else:
        target = stage2.local(args.target)
    if not target.is_file():
        raise FileNotFoundError(f"Oscillation target absent: {target}")

    running = stage2.run_dynamic_rges(uv_beta, eft_beta, config, overrides)
    f_tree, response, f_detail = stage3.dynamic_llss_response(
        running, config, reference, quartics, registry,
        factor=args.wilson_factor,
        check_scale_cancellation=not args.skip_scale_check,
    )
    hard, analytic_log, hard_detail = stage3.saved_threshold_components(
        c5_path, running, config
    )
    pieces = assemble_matching(response, hard, analytic_log)
    soft = np.asarray(pieces["dynamic_llss_only"], complex)
    with_hard = np.asarray(pieces["dynamic_llss_plus_saved_hard"], complex)
    chosen = soft if args.hard_policy == "omit" else with_hard
    alternative = with_hard if args.hard_policy == "omit" else soft

    # Both branches run the genuine SM+Weinberg solver; results are not
    # rescalings of Stage 3 or copies of the old one-scale C5 prediction.
    primary_physics = stage3.smeft_observables(chosen, running, config, target)
    alternative_physics = stage3.smeft_observables(alternative, running, config, target)
    other_policy = "retain" if args.hard_policy == "omit" else "omit"

    matching = {
        "LLSS_direct_leading_log_gev_inverse": stage1.encode_matrix(pieces["direct"]),
        "finite_joint_scalar_bubble_gev_inverse": stage1.encode_matrix(pieces["bubble"]),
        "LLSS_direct_plus_bubble_gev_inverse": stage1.encode_matrix(soft),
        "saved_renormalized_hard_diagnostic_gev_inverse": stage1.encode_matrix(hard),
        "saved_analytic_log_not_added_gev_inverse": stage1.encode_matrix(analytic_log),
        "dynamic_LLSS_plus_saved_hard_candidate_gev_inverse": stage1.encode_matrix(with_hard),
        "previous_stage3_hard_plus_direct_gev_inverse": stage1.encode_matrix(pieces["previous_stage3"]),
        "sum_identity_absolute_residual": pieces["identity_absolute_residual"],
        "sum_identity_relative_residual": pieces["identity_relative_residual"],
        "direct_log_frobenius": float(np.linalg.norm(pieces["direct"])),
        "bubble_frobenius": float(np.linalg.norm(pieces["bubble"])),
        "llss_total_frobenius": float(np.linalg.norm(soft)),
        "saved_hard_frobenius": float(np.linalg.norm(hard)),
        "primary_policy": args.hard_policy,
        "primary_coefficient_at_muS_gev_inverse": stage1.encode_matrix(chosen),
        "alternative_policy": other_policy,
        "alternative_coefficient_at_muS_gev_inverse": stage1.encode_matrix(alternative),
        "both_scalar_multiplets_removed_together": True,
        "dynamic_scalar_loop_added_exactly_once": True,
        "dynamic_LLSS_direct_log_added_exactly_once": True,
        "saved_analytic_LLSS_direct_log_added_to_result": False,
        "saved_hard_and_scalar_bubble_independence_established": False,
        "hard_subtraction_or_unverified_overlap_correction_applied": False,
        "hard_policy_interpretation": {
            "omit": "Direct LLSS log + finite two-scalar bubble; does not include saved hard.",
            "retain": "Same LLSS log + scalar bubble plus saved hard. May overlap; candidate only.",
        },
        "scalar_Feynman_integral_scheme": "one_loop_MSbar_B0bar_joint_scalar",
        "fixed_order": "one_loop_LLSS_induced_sector",
    }
    result = {
        "format": FORMAT,
        "model": MODEL,
        "status": "DynamicFermionFirstAndDynamicBothScalarsToSMEFTConditional",
        "scales_gev": {
            "UV": config["base_state"]["mu_gev"],
            **running["scales_gev"],
            "low": config["scales"]["mu_low_gev"],
        },
        "renormalisable_running": {
            "fresh_uv_solver": running["uv_solver"],
            "fresh_eft1_solver": running["eft1_solver"],
            "uv_initial": running["uv_initial"],
            "UV_at_fermion_threshold": running["uv_at_fermion_threshold"],
            "EFT1_at_fermion_threshold": running["eft1_at_fermion_threshold"],
            "EFT1_at_scalar_threshold": running["eft1_at_scalar_threshold"],
        },
        "fermion_first_LLSS": f_detail,
        "hard_threshold_reference": hard_detail,
        "matching_decomposition_gev_inverse": matching,
        "smeft": primary_physics["smeft"],
        "neutrinos": primary_physics["neutrinos"],
        "alternative_policy_result": {
            "hard_policy": other_policy,
            **_compact_physics(alternative_physics),
        },
        "sources": {
            "beta_source": beta_provenance,
            "matching_json": str(c5_path.resolve()),
            "fitted_uv_overrides": str(override_path.resolve()),
            "current_parameter_override_file": (
                str(stage2.local(args.param_overrides).resolve())
                if args.param_overrides is not None else None
            ),
            "reference_F_tree": str(paths["boundary"].resolve()),
            "native_HHSS_cache": str(paths["channels"].resolve()),
            "quartic_seed": str(paths["quartic_seed"].resolve()),
            "scalar_registry": str(paths["registry"].resolve()),
            "dynamic_F_tree_seed": str(f_tree.path),
            "spinor_factor_assumed": float(args.wilson_factor),
        },
        "flags": {
            "fresh_UV_RGBeta_RGE": True,
            "integrated_out_all_fermions_at_muF": True,
            "dynamic_matchete_tree_LLSS_from_current_UV_Yukawas": True,
            "fermion_absent_from_EFT1_renormalisable_RGE": True,
            "fresh_full_general_scalar_EFT1_RGBeta_RGE": True,
            "fixed_order_direct_LLSS_to_Weinberg_RG_mixing": True,
            "dynamic_joint_scalar_one_loop_matching": True,
            "both_scalar_multiplets_removed_at_muS": True,
            "low_energy_SM_Weinberg_RGE": True,
            "SMEFT_prediction_both_hard_policies_computed": True,
            "C_LLSS_self_running_computed": False,
            "complete_general_dimension_five_operator_running": False,
            "hard_vs_soft_overlap_certified": False,
            "Matchete_Weyl_spinor_normalization_certified": False,
            "RGBeta_Matchete_general_quartic_basis_certified": False,
            "symbolic_threshold_recomputed_each_point": False,
        },
        "scope": (
            "B(-3) full-general renormalisable UV/EFT1 running with both "
            "heavy fermions excluded after muF; current-point Matchete "
            "fermion tree C_LLSS; LLSS->C5 one-loop leading log; simultaneous "
            "finite scalar bubble for both scalar doublets, added once; "
            "SM+Weinberg evolution and oscillation observables for TWO explicit "
            "candidate hard treatments. Normalization/basis/hard-overlap remain "
            "conditional, not independently certified."
        ),
    }
    if changes is None and args.yukawa_scale == 1.0:
        prior = saved_fit["numerical_details"]
        result["baseline_reproduction"] = {
            "UV_endpoint_relative": stage2.reproduction_diagnostic(
                running["uv_at_fermion_threshold"],
                prior["uv_at_fermion_threshold"], stage2.REPRO_KEYS_UV,
            ),
            "EFT1_endpoint_relative": stage2.reproduction_diagnostic(
                running["eft1_at_scalar_threshold"],
                prior["eft1_at_scalar_threshold"], stage2.REPRO_KEYS_EFT1,
            ),
        }

    stage1.write_outputs(out_dir, result)
    _alternative_csv(out_dir / "observables_alternative.csv", alternative_physics)
    np.savez_compressed(
        out_dir / "dynamic_matching_C5.npz",
        C_LLSS_muF_gev_inverse=np.asarray(f_tree.raw, complex),
        C_LLSS_by_source_muF_gev_inverse=np.asarray(f_tree.by_source, complex),
        C5_direct_log_gev_inverse=np.asarray(pieces["direct"], complex),
        C5_joint_scalar_bubble_gev_inverse=np.asarray(pieces["bubble"], complex),
        C5_LLSS_total_gev_inverse=soft,
        C5_saved_hard_gev_inverse=hard,
        C5_with_saved_hard_candidate_gev_inverse=with_hard,
        C5_primary_muS_gev_inverse=chosen,
        C5_alternative_muS_gev_inverse=alternative,
        muF_gev=np.asarray(float(config["scales"]["mu_fermion_threshold_gev"])),
        muS_gev=np.asarray(float(config["scales"]["mu_scalar_threshold_gev"])),
    )

    print("Wrote:", out_dir / "result.json")
    print("Fresh UV/EFT1 RGE evaluations:",
          running["uv_solver"]["nfev"], "/", running["eft1_solver"]["nfev"])
    print("Dynamic F-tree LLSS tensor:", list(f_tree.raw.shape))
    print("HHSS channels:", f_detail["HHSS_channel_count"])
    print("Direct-log norm:", f"{np.linalg.norm(pieces['direct']):.9g}")
    print("Dynamic scalar-bubble norm:", f"{np.linalg.norm(pieces['bubble']):.9g}")
    print("Dynamic LLSS total norm:", f"{np.linalg.norm(soft):.9g}")
    print("Scale-cancellation residual:", f_detail["matching_scale_cancellation_residual"])
    for name, phys in ((args.hard_policy, primary_physics),
                       (other_policy, alternative_physics)):
        obs = phys["neutrinos"]
        print(f"hard-policy={name}: chi2 = {obs['oscillation_chi2']:.9g}")
        print("  masses (eV):", ", ".join(
            f"{m:.9g}" for m in obs["masses_ev_m1_m2_m3"]
        ))
    print("NOTICE: hard/bubble overlap and GammaCC/Weyl normalization remain provisional")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
