"""Compare UV/EFT1 renormalisable RG running on/off in the NEW 16-T3 EFT chain.

All four cases use the SAME SAVED UV fitted parameters and the SAME thresholds:

    U1_E1  UV and EFT1 renormalisable beta functions evolved
    U1_E0  UV evolved, EFT1 renormalisable couplings held fixed
    U0_E1  UV couplings fixed, EFT1 renormalisable beta functions evolved
    U0_E0  both renormalisable stages held fixed

Every case recalculates its actual fermion-tree LLSS boundary, direct one-loop
LLSS -> C5 logarithm, joint finite two-scalar bubble and SM+Weinberg evolution.
The Wilson direct log and final SMEFT evolution are ON in ALL four cases.
Thus 'frozen' refers ONLY to renormalisable beta evolution in the specified
UV/EFT1 interval. Threshold matching, the fixed-order Wilson logarithm and
final SMEFT evolution are NEVER silently removed.

The B_-1/C_-1 calculation uses native source-exact Matchete CGs with U(muS)
in the muF LLSS source, as in RunT3All16DynamicSequential16. The source-
basis identification, absolute GammaCC normalisation, missing hard one-loop
matching and uncoupled full EFT1 Wilson transport remain unverified.

No fit is performed; all cases use the saved fitted UV input parameters.
Writes new output folders only, never changes prior batch16_v1 or source files.

PowerShell from repository root:
    python -X utf8 -m scripts.CompareT3SequentialRunningFrozen16 --model B_-1 --run-id pilot_bm1
    python -X utf8 -m scripts.CompareT3SequentialRunningFrozen16 --model C_-1 --run-id pilot_cm1
    python -X utf8 -m scripts.CompareT3SequentialRunningFrozen16 --all --run-id running_ablation16_v1
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from scripts import RunT3All16DynamicSequential as base
from scripts import RunT3All16DynamicSequential16 as mixed

OUT = base.ROOT / "output/general_sequential_pilot/all16_running_frozen_ablation"
CASE_KEYS = ((True, True), (True, False), (False, True), (False, False))
OBS_KEYS = ("delta_m21_sq_ev2", "delta_m3l_sq_ev2",
            "sin2_theta12", "sin2_theta13", "sin2_theta23")
SCHEMA = "t3-sequential-soft-16-uv-eft1-renormalisable-ablation-v1"


def case_name(uv_on: bool, eft_on: bool) -> str:
    return f"U{int(uv_on)}_E{int(eft_on)}"


def _assert_finite_real(value: Any, name: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{name}: nonfinite numerical result")
    return number


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _enable_mixed_hooks() -> None:
    """Use the same source-specific hooks as the successful 16-model runner."""
    base.fit_input = mixed._mixed_fit_input
    base.check_model = mixed._check_model
    base.dynamic_soft_sector = mixed._dynamic_soft_sector
    base.FORMAT = mixed.FORMAT


def _initial_state(prepared: Any):
    from Numerical.orchestration.GeneralRunningBenchmark import initial_state
    scales = prepared.config["scales"]
    mu_uv = float(prepared.config["base_state"]["mu_gev"])
    mu_f = float(scales["mu_fermion_threshold_gev"])
    mu_s = float(scales["mu_scalar_threshold_gev"])
    model, state = initial_state(
        prepared.uv_beta, mu_uv=mu_uv, mu_f=mu_f, mu_s=mu_s,
        overrides=prepared.overrides,
    )
    if model != prepared.model:
        raise ValueError(f"Prepared and reconstructed UV model differ: {model}")
    return state


def _reconstruct_running(prepared: Any, uv_on: bool, eft_on: bool,
                         uv_cache: dict[bool, tuple]) -> dict:
    """Reuse the repository's RGE solvers, or transport states unchanged in scale.

    Fermion matching/projecting is performed even in frozen scenarios. No
    UV/fermion state is used below muF, and no heavy fermion beta survives EFT1.
    """
    from Numerical.orchestration.GeneralRunningBenchmark import _state_values
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.FermionThresholdBoundary import (
        build_intermediate_scalar_boundary,
    )
    from Numerical.running.IntermediateScalarRunner import (
        run_intermediate_scalar_segment,
    )
    from Numerical.running.IntermediateScalarStateVector import (
        pack_intermediate_scalar_state,
    )

    scales = prepared.config["scales"]
    mu_uv = float(prepared.config["base_state"]["mu_gev"])
    mu_f = float(scales["mu_fermion_threshold_gev"])
    mu_s = float(scales["mu_scalar_threshold_gev"])
    if not (mu_uv > mu_f > mu_s > 0):
        raise ValueError("Expected UV > fermion threshold > scalar threshold")

    if uv_on not in uv_cache:
        initial = _initial_state(prepared)
        if uv_on:
            trajectory = run_uv_segment(
                initial, prepared.uv_beta, mu_f, rtol=1e-7, atol=1e-10
            )
            if not trajectory.solver_success or trajectory.nfev <= 0:
                raise RuntimeError(f"{prepared.model}: UV RGE solver failed")
            final_uv = trajectory.final_state
            uv_meta = {"success": True, "nfev": int(trajectory.nfev),
                       "evolved": True, "frozen": False}
        else:
            final_uv = replace(initial, mu_gev=mu_f).validated()
            uv_meta = {"success": True, "nfev": 0,
                       "evolved": False, "frozen": True}
        uv_cache[uv_on] = (initial, final_uv, uv_meta)

    initial, final_uv, uv_meta = uv_cache[uv_on]
    eft_start = build_intermediate_scalar_boundary(
        final_uv, matching_scale_gev=mu_f
    )
    _, layout = pack_intermediate_scalar_state(eft_start)
    expected = set(prepared.eft_beta.get("report_betas", {}))
    if set(layout.names) != expected:
        raise ValueError(
            f"{prepared.model}: EFT1 beta/state inventory mismatch: "
            f"missing={sorted(expected-set(layout.names))}, "
            f"extra={sorted(set(layout.names)-expected)}"
        )
    if eft_on:
        trajectory = run_intermediate_scalar_segment(
            eft_start, prepared.eft_beta, mu_s, rtol=1e-7, atol=1e-10
        )
        if not trajectory.solver_success or trajectory.nfev <= 0:
            raise RuntimeError(f"{prepared.model}: EFT1 RGE solver failed")
        eft_end = trajectory.final_state
        eft_meta = {"success": True, "nfev": int(trajectory.nfev),
                    "evolved": True, "frozen": False}
    else:
        eft_end = replace(eft_start, mu_gev=mu_s).validated()
        eft_meta = {"success": True, "nfev": 0,
                    "evolved": False, "frozen": True}

    return {
        "model": prepared.model,
        "status": "NumericalRenormalisableRunningOnly",
        "general_matching_performed": False,
        "scales_gev": {"UV": mu_uv, "F": mu_f, "S": mu_s},
        "uv_solver": dict(uv_meta),
        "eft1_solver": dict(eft_meta),
        "uv_initial": _state_values(initial),
        "uv_at_fermion_threshold": _state_values(final_uv),
        "eft1_at_fermion_threshold": _state_values(eft_start),
        "eft1_at_scalar_threshold": _state_values(eft_end),
    }


def _observables(result: dict) -> dict:
    obs = result["neutrinos"]
    pred = obs["oscillation_observables"]
    if not all(k in pred for k in OBS_KEYS):
        raise ValueError(f"Missing fitted observable keys: {set(OBS_KEYS)-set(pred)}")
    masses = [_assert_finite_real(x, "mass") for x in obs["masses_ev_m1_m2_m3"]]
    if len(masses) != 3:
        raise ValueError("Expected three neutrino mass eigenvalues")
    out = {"chi2": _assert_finite_real(obs["oscillation_chi2"], "chi2"),
           **{k: _assert_finite_real(pred[k], k) for k in OBS_KEYS},
           **{f"m{i}_ev": m for i, m in enumerate(masses, start=1)}}
    for key in ("sin2_theta12", "sin2_theta13", "sin2_theta23"):
        value = out[key]
        if not 0 <= value <= 1:
            raise ValueError(f"{key} outside [0,1]")
        out[key.replace("sin2_", "")+"_deg"] = math.degrees(math.asin(math.sqrt(value)))
    return out


def _shape_alignment(reference: np.ndarray, candidate: np.ndarray) -> dict:
    """Best complex rescaling of candidate onto reference; flavour residual."""
    ref = np.asarray(reference, dtype=complex)
    test = np.asarray(candidate, dtype=complex)
    if ref.shape != (3, 3) or test.shape != (3, 3):
        raise ValueError("C5 alignment requires 3x3 matrices")
    ref_norm = float(np.linalg.norm(ref))
    test_norm = float(np.linalg.norm(test))
    if min(ref_norm, test_norm) <= 1e-35:
        return {"status": "too_small_for_stable_alignment"}
    a = np.vdot(test, ref) / np.vdot(test, test)
    residual = float(np.linalg.norm(ref-a*test) / ref_norm)
    return {"status": "ok",
            "fit_reference_equals_a_times_candidate": {"real": float(a.real),
                                                       "imag": float(a.imag)},
            "residual_after_complex_rescaling": residual,
            "norm_ratio_reference_over_candidate": ref_norm/test_norm}


def _comparisons(records: dict[str, dict]) -> dict:
    """Stage effects are *matched* contrasts, not fit-quality statements."""
    full = records["U1_E1"]["observables"]
    case_observables = {k: v["observables"] for k, v in records.items()}
    keys = tuple(full.keys())
    def contrast(after: dict, before: dict) -> dict:
        output = {}
        for key in keys:
            delta = after[key] - before[key]
            output[key] = {"difference": delta,
                           "percent_of_before": (100*delta/before[key]
                                                 if before[key] != 0 else None)}
        return output
    a = case_observables
    interaction = {
        key: a["U1_E1"][key] - a["U1_E0"][key]
             - a["U0_E1"][key] + a["U0_E0"][key]
        for key in keys
    }
    return {
        "full_minus_all_renormalisable_frozen": contrast(a["U1_E1"], a["U0_E0"]),
        "UV_effect_with_EFT1_running": contrast(a["U1_E1"], a["U0_E1"]),
        "UV_effect_with_EFT1_frozen": contrast(a["U1_E0"], a["U0_E0"]),
        "EFT1_effect_with_UV_running": contrast(a["U1_E1"], a["U1_E0"]),
        "EFT1_effect_with_UV_frozen": contrast(a["U0_E1"], a["U0_E0"]),
        "UV_x_EFT1_nonlinear_interaction": interaction,
    }


def _reference_comparison(prepared: Any, run: dict, reference_id: str) -> dict:
    """Optional exact baseline replay check. Never use old C5 as an input."""
    path = base.OUTPUT_ROOT / reference_id / base.model_slug(prepared.model) / "result.json"
    if not path.is_file():
        return {"status": "not_available", "expected_path": str(path)}
    old = base.stage1.read_object(path)
    if old.get("model") != prepared.model:
        raise ValueError("Baseline reference file belongs to another model")
    old_c5 = base.stage1.finite_matrix(
        old["scalar_matching"]["C5_used_at_S_gev_inverse"], "old C5"
    )
    new_c5 = base.stage1.finite_matrix(
        run["scalar_matching"]["C5_used_at_S_gev_inverse"], "new C5"
    )
    norm = max(float(np.linalg.norm(old_c5)), float(np.linalg.norm(new_c5)), 1e-30)
    delta_c5 = float(np.linalg.norm(old_c5-new_c5) / norm)
    old_chi = float(old["neutrinos"]["oscillation_chi2"])
    new_chi = float(run["neutrinos"]["oscillation_chi2"])
    chi_diff = abs(new_chi-old_chi)/max(abs(new_chi),abs(old_chi),1e-30)
    return {"status": "checked", "reference_path": str(path),
            "C5_relative_difference": delta_c5,
            "chi2_relative_difference": float(chi_diff),
            "reproduces_reference_to_1e_minus_7": bool(delta_c5 < 1e-7 and chi_diff < 1e-7),
            "reference_run_id": reference_id}


def _calculate_case(prepared: Any, args: argparse.Namespace,
                    running: dict, directory: Path,
                    uv_on: bool, eft_on: bool) -> dict:
    """Reuse the production 16-model matching and SMEFT implementation unchanged."""
    original = base.fresh_running
    try:
        def selected(p):
            if p is not prepared:
                raise ValueError("Unexpected model entry for injected RGE snapshot")
            return running
        base.fresh_running = selected
        result = base.calculate_model(prepared, args, directory)
    finally:
        base.fresh_running = original

    # The underlying writer assumes both RGE stages are always on. Correct its
    # descriptive flags in this newly generated file; do not alter any physics.
    result["format"] = SCHEMA
    result["run_switches"] = {
        "UV_renormalisable_RGE": uv_on,
        "EFT1_renormalisable_RGE": eft_on,
        "EFT1_LLSS_to_C5_fixed_one_loop_log": True,
        "both_scalar_matching": True,
        "SMEFT_RGE": True,
        "frozen_stage_semantics": "couplings held constant; matching scales advanced",
        "UV_nfev": running["uv_solver"]["nfev"],
        "EFT1_nfev": running["eft1_solver"]["nfev"],
    }
    result["flags"]["UV_evolved_fresh"] = uv_on
    result["flags"]["EFT1_renormalisable_RGE_evolved"] = eft_on
    result["flags"]["mixed_muS_scalar_rotation_reused_at_F"] = (
        prepared.model in mixed.MIXED
    )
    result["scope"] = (
        "UV/EFT1 renormalisable RGE ablation: frozen stages transport the "
        "same renormalisable state to the endpoint scale; all fermions are "
        "integrated at F, both scalars at S; direct fixed-order LLSS->C5 "
        "log and final SMEFT evolution always on. This is the exploratory "
        "soft sector; hard matching and physical normalization unverified."
    )
    _write_json(directory / "result.json", result)
    return result


def _row_for(model: str, case: str, result: dict) -> dict:
    obs = _observables(result)
    c5 = base.stage1.finite_matrix(
        result["scalar_matching"]["C5_used_at_S_gev_inverse"], "C5 at S"
    )
    return {"model": model, "case": case,
            "UV_RGE": int(result["run_switches"]["UV_renormalisable_RGE"]),
            "EFT1_RGE": int(result["run_switches"]["EFT1_renormalisable_RGE"]),
            "UV_nfev": result["run_switches"]["UV_nfev"],
            "EFT1_nfev": result["run_switches"]["EFT1_nfev"],
            "C5_norm_gev_inv": float(np.linalg.norm(c5)), **obs}


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parser() -> argparse.ArgumentParser:
    p = base.parser()
    p.description = __doc__
    p.set_defaults(output_root=OUT)
    p.add_argument("--reference-run-id", default="batch16_v1",
                   help="Optional prior full-run ID for exact U1_E1 replay check")
    return p


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    _enable_mixed_hooks()
    if not isinstance(args.run_id, str) or not base.re.fullmatch(
        r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", args.run_id
    ):
        raise ValueError("Invalid --run-id")
    if args.all and args.parameter_overrides:
        raise ValueError("--parameter-overrides is only supported for a single model")
    if not math.isfinite(args.yukawa_scale) or args.yukawa_scale <= 0:
        raise ValueError("--yukawa-scale must be finite and positive")
    if not math.isfinite(args.wilson_factor) or args.wilson_factor == 0:
        raise ValueError("--wilson-factor must be finite and nonzero")
    if args.b3_hard_policy != "omit":
        raise ValueError("Only --b3-hard-policy omit is supported for like-for-like soft-sector ablation")
    targets = base.ALL_MODELS if args.all else (args.model,)
    if not args.check_only and not base.absolute(args.target).is_file():
        raise FileNotFoundError(f"Oscillation target unavailable: {args.target}")
    root = base.absolute(args.output_root) / args.run_id
    if not args.check_only and root.exists():
        raise FileExistsError(f"Output root already exists: {root}")

    rows: list[dict] = []
    effect_rows: list[dict] = []
    summary = {"format": SCHEMA, "run_id": args.run_id, "models": {},
               "target": str(base.absolute(args.target)),
               "reference_run_id": args.reference_run_id,
               "case_definitions": {
                   case_name(u, e): {"UV_renormalisable_RGE": u,
                                     "EFT1_renormalisable_RGE": e,
                                     "LLSS_to_C5_log": True,
                                     "SMEFT_RGE": True}
                   for u, e in CASE_KEYS},
               "scope": "No refitting. New dynamic soft-only sequential C5; hard match and absolute basis/normalisation remain open. For mixed B/C, U(muS) also enters the F boundary and changes with EFT1 running."}
    failed = 0
    for model in targets:
        try:
            prepared = mixed._preflight(model, args)
            if args.check_only:
                print(f"{model}: READY, UV/EFT1 beta files and matched scalar sources", flush=True)
                summary["models"][model] = {"status": "ready"}
                continue
            cache: dict[bool, tuple] = {}
            records = {}
            for uv_on, eft_on in CASE_KEYS:
                case = case_name(uv_on, eft_on)
                running = _reconstruct_running(prepared, uv_on, eft_on, cache)
                folder = root / base.model_slug(model) / case
                result = _calculate_case(prepared, args, running, folder, uv_on, eft_on)
                row = _row_for(model, case, result)
                rows.append(row)
                scenario_obs = _observables(result)
                scenario_obs["C5_norm_gev_inv"] = row["C5_norm_gev_inv"]
                records[case] = {"observables": scenario_obs,
                                 "result_json": str(folder / "result.json")}
                print(f"{model} {case}: chi2={row['chi2']:.8g}; "
                      f"dm21={row['delta_m21_sq_ev2']:.8g}; "
                      f"dm3l={row['delta_m3l_sq_ev2']:.8g}; "
                      f"sin2theta12={row['sin2_theta12']:.7g}; "
                      f"sin2theta23={row['sin2_theta23']:.7g}", flush=True)
            contrasts = _comparisons(records)
            matrices = {}
            for c in records:
                record = base.stage1.read_object(Path(records[c]["result_json"]))
                matrices[c] = base.stage1.finite_matrix(
                    record["scalar_matching"]["C5_used_at_S_gev_inverse"],
                    "C5 at scalar threshold"
                )
            contrasts["C5_shape_change_vs_full"] = {
                c: _shape_alignment(matrices["U1_E1"], matrix)
                for c, matrix in matrices.items() if c != "U1_E1"
            }
            full = base.stage1.read_object(root / base.model_slug(model) / "U1_E1" / "result.json")
            ref = _reference_comparison(prepared, full, args.reference_run_id)
            if ref["status"] == "checked" and not ref["reproduces_reference_to_1e_minus_7"]:
                print(f"WARNING {model}: U1_E1 does not reproduce {args.reference_run_id}; "
                      f"relative C5 difference {ref['C5_relative_difference']:.3g}", flush=True)
            summary["models"][model] = {
                "status": "completed", "scenarios": records,
                "stage_effects": contrasts,
                "full_run_reproduction": ref,
            }
            for contrast_name in (
                "full_minus_all_renormalisable_frozen",
                "UV_effect_with_EFT1_running",
                "EFT1_effect_with_UV_running",
                "UV_effect_with_EFT1_frozen",
                "EFT1_effect_with_UV_frozen",
            ):
                values = contrasts[contrast_name]
                effect_rows.append({
                    "model": model, "effect": contrast_name,
                    **{key+"_delta": values[key]["difference"]
                       for key in ("chi2", "delta_m21_sq_ev2", "delta_m3l_sq_ev2",
                                   "theta12_deg", "theta13_deg", "theta23_deg",
                                   "C5_norm_gev_inv")},
                    **{key+"_percent": values[key]["percent_of_before"]
                       for key in ("delta_m21_sq_ev2", "delta_m3l_sq_ev2",
                                   "C5_norm_gev_inv")},
                })
            # A compact per-model numerical table, including angle differences.
            _write_csv(root / base.model_slug(model) / "cases.csv",
                       [r for r in rows if r["model"] == model])
            print(f"{model}: EFT1 effect at UV-running: "
                  f"delta(dm21)={contrasts['EFT1_effect_with_UV_running']['delta_m21_sq_ev2']['difference']:+.5g} eV^2, "
                  f"delta(theta12)={contrasts['EFT1_effect_with_UV_running']['theta12_deg']['difference']:+.5g} deg; "
                  f"UV effect with EFT1-running: "
                  f"delta(dm21)={contrasts['UV_effect_with_EFT1_running']['delta_m21_sq_ev2']['difference']:+.5g} eV^2", flush=True)
        except Exception as exc:
            failed += 1
            summary["models"][model] = {"status": "failed", "error": str(exc),
                                        "type": type(exc).__name__}
            print(f"{model}: FAILED {type(exc).__name__}: {exc}", flush=True)
    summary["counts"] = {"models_requested": len(targets),
                         "models_completed": sum(x["status"] == "completed"
                                                 for x in summary["models"].values()),
                         "models_failed": failed,
                         "cases_completed": len(rows)}
    if args.check_only:
        print("Check only: no RGEs, matching or output generated.")
    else:
        _write_csv(root / "cases.csv", rows)
        _write_csv(root / "effects.csv", effect_rows)
        _write_json(root / "summary.json", summary)
        print(f"Summary: {root / 'summary.json'}")
        print(f"Cases: {root / 'cases.csv'}")
        print(f"Stage differences: {root / 'effects.csv'}")
    print("Counts:", summary["counts"])
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
