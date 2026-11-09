"""Rebuild the D(-2)/A(-4) shared-Yukawa thesis benchmarks with sequential EFT.

Physics: UV -> integrate out ALL fermions -> EFT1(S1,S2) -> integrate out
ALL scalars jointly -> SMEFT. No one-scale C5 is added. The current
conditional fixed-one-loop LLSS log + bubble (hard matching omitted) is used
in every evaluated point, including optimization and frozen comparisons.

Stages (write-once under output/general_sequential_pilot/thesis_sequential_refresh):
    --phase fit       Optimise ONLY D_-2 and A_-4, holding other inputs fixed.
    --phase transfer  Fix each optimised (Y1,Y2) and run all 16 representations,
                      including three one-stage-frozen tests plus U0_E0.
    --phase all       Fit both donors, then run both 16-model transfer campaigns.

Every transfer keeps the *target model's own* scalar/SM/mass inputs. Additional
model-specific Yukawa channels y1Conj/y2Conj are ZERO by default: the thesis
Appendix B sets D(-2)'s extra Y2Conj to zero. Choose --extra-yukawas saved
only to retain historical model-specific inputs (not a pure shared-pair test).

The new 16-model wrapper and running-frozen comparison scripts must already be
installed. No old results, historical figures or GitHub content are modified.

PowerShell:
    python -X utf8 -m scripts.RunT3SequentialThesisRefresh --phase fit --run-id thesis_seq_v1
    python -X utf8 -m scripts.RunT3SequentialThesisRefresh --phase transfer --run-id thesis_seq_v1
    python -X utf8 -m scripts.PlotT3SequentialThesisRefresh --run-id thesis_seq_v1
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import traceback
from typing import Any

import numpy as np
from scipy.optimize import least_squares

from scripts import RunT3All16DynamicSequential as base
from scripts import RunT3All16DynamicSequential16 as mixed
from scripts import CompareT3SequentialRunningFrozen16 as ablation

ROOT = base.ROOT
OUT = ROOT / "output/general_sequential_pilot/thesis_sequential_refresh"
DONORS = ("D_-2", "A_-4")
MODEL_ORDER = base.ALL_MODELS
OBS_KEYS = ("delta_m21_sq_ev2", "delta_m3l_sq_ev2",
            "sin2_theta12", "sin2_theta13", "sin2_theta23")
FORMAT = "t3-thesis-Dm2-Am4-sequential-soft-fixed-Yukawa-v1"
FIT_THRESHOLD = 5.86  # project-specific screening, not an official NuFIT confidence level


def jread(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(str(path))
    v = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(v, dict):
        raise ValueError(f"Expected JSON object at {path}")
    return v


def jwrite(path: Path, v: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(v, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def csvwrite(path: Path, records: list[dict]) -> None:
    if not records:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]))
        w.writeheader()
        w.writerows(records)


def hash_path(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def model_args(model: str, args: argparse.Namespace) -> argparse.Namespace:
    return base.parser().parse_args([
        "--model", model, "--run-id", args.run_id,
        "--campaign-root", str(args.campaign_root),
        "--canonical-run-id", args.canonical_run_id,
        "--beta-source", args.beta_source,
        "--target", str(args.target),
        "--wilson-factor", str(args.wilson_factor),
        "--b3-hard-policy", "omit",
    ])


def prepare(model: str, args: argparse.Namespace):
    p = mixed._preflight(model, model_args(model, args))
    if p.model != model:
        raise ValueError("Prepared model identity mismatch")
    return p


def apply_extra_yukawa_policy(params: dict, policy: str) -> dict:
    new = dict(params)
    if policy == "zero":
        zeros = np.zeros((3, 3), dtype=complex)
        for name in ("y1Conj", "y2Conj"):
            if name in new:
                new[name] = base.stage2._matrix_to_entries(zeros)
    elif policy != "saved":
        raise ValueError(f"Unexpected extra Yukawa policy: {policy}")
    return new


def put_yukawas(params: dict, y1: np.ndarray, y2: np.ndarray) -> dict:
    params = dict(params)
    for name, matrix in (("y1", y1), ("y2", y2)):
        z = np.asarray(matrix, dtype=complex)
        if z.shape != (3, 3) or not np.isfinite(z).all():
            raise ValueError(f"Invalid {name} matrix")
        params[name] = base.stage2._matrix_to_entries(z)
    return params


def y_from_params(params: dict) -> tuple[np.ndarray, np.ndarray]:
    return tuple(base.stage1.finite_matrix(params[name], name)
                 for name in ("y1", "y2"))


def real_y_vector(params: dict) -> np.ndarray:
    y1, y2 = y_from_params(params)
    if max(np.max(np.abs(y1.imag)), np.max(np.abs(y2.imag))) > 1e-13:
        raise ValueError("This thesis fit is restricted to 18 real Yukawa entries")
    return np.concatenate((y1.real.ravel(), y2.real.ravel()))


def pair_from_real_vector(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    v = np.asarray(x, dtype=float)
    if v.shape != (18,) or not np.isfinite(v).all():
        raise ValueError("Expected 18 real finite Yukawa entries")
    return v[:9].reshape((3, 3)).copy(), v[9:].reshape((3, 3)).copy()


def extract_obs(payload: dict) -> dict:
    n = payload["neutrinos"]
    out = {k: float(n["oscillation_observables"][k]) for k in OBS_KEYS}
    out["chi2"] = float(n["oscillation_chi2"])
    masses = [float(x) for x in n["masses_ev_m1_m2_m3"]]
    if len(masses) != 3 or not all(math.isfinite(x) and x >= 0 for x in masses):
        raise ValueError("Invalid neutrino masses")
    out.update({f"m{i}_ev": masses[i-1] for i in (1, 2, 3)})
    if not all(math.isfinite(x) for x in out.values()):
        raise ValueError("Nonfinite neutrino observable")
    return out


def target_data(path: Path) -> tuple[np.ndarray, np.ndarray]:
    t = jread(path)
    central = np.asarray([t["central_values"][k] for k in OBS_KEYS], float)
    sigma = np.asarray([t["one_sigma_errors"][k] for k in OBS_KEYS], float)
    if not (np.isfinite(central).all() and np.isfinite(sigma).all()
            and (sigma > 0).all()):
        raise ValueError("Invalid diagonal Gaussian target")
    return central, sigma


def predict(prepared, running: dict, target: Path, factor: float,
            *, sm_on: bool, check_scale: bool) -> tuple[dict, np.ndarray]:
    """Pure evaluation: no output writes, no archived one-scale C5 substitution."""
    source, bundle, response, source_info = base.dynamic_soft_sector(
        prepared, running, factor=factor, check_scale=check_scale)
    c5 = base.symmetric_c5(response["c5_LLSS_total"], "sequential LLSS soft C5")
    if sm_on:
        obs = base.stage3.smeft_observables(c5, running, prepared.config, target)
    else:
        obs = smeft_frozen(c5, running, prepared.config, target)
    sm_info = obs["smeft"]
    return {
        "neutrinos": obs["neutrinos"],
        "smeft": sm_info,
        "C5_at_muS": base.stage1.encode_matrix(c5),
        "direct_log_norm": float(np.linalg.norm(base.symmetric_c5(
            response["c5_eft1_leading_log"], "log"))),
        "scalar_bubble_norm": float(np.linalg.norm(base.symmetric_c5(
            response["c5_scalar_finite_bubble"], "bubble"))),
        "C5_norm": float(np.linalg.norm(c5)),
        "scalar_mass2_at_S": np.linalg.eigvalsh(bundle.mass2).tolist(),
        # response() deliberately returns None when the optional scale check
        # is disabled during the many optimizer trial evaluations.  Keep this
        # diagnostic absent, rather than casting None to float.  The final
        # fitted point and all transfer points still enable the check.
        "matching_scale_cancellation": (
            None if response["matching_scale_cancellation_residual"] is None
            else float(response["matching_scale_cancellation_residual"])),
        "beta_projection_residual": float(response["weinberg_beta_projection_residual"]),
        "bubble_projection_residual": float(response["weinberg_bubble_projection_residual"]),
        "native_source_hash": source.hash,
        "native_source_count": len(source.source_indices),
    }, c5


def smeft_frozen(c5: np.ndarray, running: dict, config: dict,
                 target: Path) -> dict:
    """Same scalar-threshold initial condition, but zero final SMEFT evolution."""
    from Numerical.fitting.OscillationFit import (
        OscillationFitTarget, evaluate_mass_matrix_fit,
    )
    from physics.NeutrinoTrajectory import (
        charged_lepton_mass_basis_matrix, neutrino_mass_from_c5,
    )
    from physics.NeutrinoObservables import calculate_neutrino_observables

    eft = running["eft1_at_scalar_threshold"]
    vev = float(config.get("numerical", {}).get("vev_gev", 246.22))
    ye = base.stage1.finite_matrix(eft["ye"], "ye")
    mnu = neutrino_mass_from_c5(c5, vev_gev=vev)
    charged = charged_lepton_mass_basis_matrix(mnu, ye)
    target_obj = OscillationFitTarget.from_json(target)
    fit = evaluate_mass_matrix_fit(charged, target_obj)
    obs = calculate_neutrino_observables(charged, ordering=target_obj.ordering)
    return {
        "smeft": {
            "RGEs_evolved": False,
            "solver_success": True,
            "solver_nfev": 0,
            "C5_at_muS_gev_inverse": base.stage1.encode_matrix(c5),
            "C5_at_muLow_gev_inverse": base.stage1.encode_matrix(c5),
            "note": "Held all SM+Weinberg couplings fixed from muS to muLow; no ODE",
        },
        "neutrinos": {
            "ordering": obs.ordering,
            "vev_gev": vev,
            "masses_ev_m1_m2_m3": [float(x) for x in obs.masses_ev],
            "mass_matrix_gev": base.stage1.encode_matrix(mnu),
            "mass_matrix_charged_lepton_basis_gev": base.stage1.encode_matrix(charged),
            "PMNS_absolute": np.asarray(obs.pmns_abs).tolist(),
            "takagi_residual": float(obs.takagi_residual),
            "oscillation_observables": {
                name: float(value) for name, value in zip(
                    fit.observable_names, fit.prediction)},
            "oscillation_chi2": float(fit.chi2),
            "target_ordering": str(target_obj.ordering),
        },
    }


def evaluate_full(prepared, overrides: dict, args: argparse.Namespace,
                  *, check_scale: bool) -> dict:
    p = replace(prepared, overrides=overrides)
    running = base.fresh_running(p)
    predicted, _ = predict(p, running, args.target, args.wilson_factor,
                           sm_on=True, check_scale=check_scale)
    return predicted


def fit_donor(model: str, prepared, args: argparse.Namespace, out: Path) -> dict:
    initial_params = apply_extra_yukawa_policy(prepared.overrides, args.extra_yukawas)
    x_init = real_y_vector(initial_params)
    central, sigma = target_data(args.target)
    scale = max(float(np.max(np.abs(x_init))), 1e-8)
    x0 = x_init / scale
    count, failures = 0, 0
    best_chi2 = float("inf")
    best_x = x0.copy()

    def residual(x: np.ndarray) -> np.ndarray:
        nonlocal count, failures, best_chi2, best_x
        count += 1
        y1, y2 = pair_from_real_vector(scale * x)
        overrides = put_yukawas(initial_params, y1, y2)
        try:
            prediction = evaluate_full(prepared, overrides, args, check_scale=False)
            observables = prediction["neutrinos"]["oscillation_observables"]
            r = (np.asarray([observables[k] for k in OBS_KEYS]) - central) / sigma
            if not np.isfinite(r).all():
                raise ValueError("Nonfinite fit residual")
            chi2 = float(r @ r)
            if chi2 < best_chi2:
                best_chi2, best_x = chi2, x.copy()
                if count % 20 != 0 or chi2 < FIT_THRESHOLD:
                    print(f"{model}: evaluation {count}, best chi2={chi2:.8g}", flush=True)
            return r
        except (ValueError, FloatingPointError, OverflowError, RuntimeError) as exc:
            failures += 1
            if failures < 4:
                print(f"{model}: rejected trial ({type(exc).__name__}): {exc}", flush=True)
            return np.full(5, float(args.invalid_penalty))

    first_residual = residual(x0)
    if not np.isfinite(first_residual).all() or np.max(np.abs(first_residual)) >= args.invalid_penalty:
        raise RuntimeError(f"{model}: initial fitted seed cannot be evaluated")
    initial_chi2 = float(first_residual @ first_residual)
    opt = least_squares(
        residual, x0, method="trf", jac="2-point", diff_step=args.diff_step,
        x_scale="jac", max_nfev=args.max_nfev,
        ftol=1e-9, xtol=1e-9, gtol=1e-9,
        bounds=(-np.full(18, args.max_yukawa / scale),
                 np.full(18, args.max_yukawa / scale)),
    )
    # Keep the actual best finite trial (optimizer may stop at an invalid point).
    candidate_resid = residual(opt.x)
    if float(candidate_resid @ candidate_resid) < best_chi2:
        best_chi2 = float(candidate_resid @ candidate_resid)
        best_x = opt.x.copy()
    y1, y2 = pair_from_real_vector(best_x * scale)
    final_overrides = put_yukawas(initial_params, y1, y2)
    final = evaluate_full(prepared, final_overrides, args, check_scale=True)
    chi2 = float(final["neutrinos"]["oscillation_chi2"])
    status = "accepted_conditional_fit" if chi2 <= FIT_THRESHOLD else "poor_conditional_fit"
    result = {
        "format": FORMAT,
        "kind": "donor_fit",
        "model": model,
        "status": status,
        "chi2_screening_threshold": FIT_THRESHOLD,
        "initial_chi2": initial_chi2,
        "final_chi2": chi2,
        "optimizer": {
            "method": "scipy.optimize.least_squares/trf, 18 real Yukawa parameters",
            "success": bool(opt.success), "message": str(opt.message),
            "scipy_reported_nfev": int(opt.nfev),
            "scipy_reported_njev": None if opt.njev is None else int(opt.njev),
            "actual_physics_evaluations": count,
            "failed_trial_evaluations": failures,
            "max_nfev": args.max_nfev,
            "max_abs_yukawa_bound": args.max_yukawa,
            "finite_difference_step": args.diff_step,
            "best_seen_chi2": best_chi2,
        },
        "y1": base.stage1.encode_matrix(y1),
        "y2": base.stage1.encode_matrix(y2),
        "initial_y1": base.stage1.encode_matrix(pair_from_real_vector(x_init)[0]),
        "initial_y2": base.stage1.encode_matrix(pair_from_real_vector(x_init)[1]),
        "observables": extract_obs(final),
        "C5_at_muS": final["C5_at_muS"],
        "C5_at_muLow": final["smeft"]["C5_at_muLow_gev_inverse"],
        "donor_fit_source": str(prepared.fit_root.resolve()),
        "donor_original_fit_chi2": prepared.fit.get("chi2_final"),
        "source_hashes": {
            "fitted_config": hash_path(prepared.fit_root / "fitted_config.json"),
            "saved_result": hash_path(prepared.fit_root / "result.json"),
            "target": hash_path(args.target),
        },
        "coupling_policy": {
            "extra_yukawas": args.extra_yukawas,
            "other_model_parameters": "saved target-model benchmark inputs, frozen during fit",
        },
        "limitations": "Conditional LLSS log + finite scalar bubble; hard matching and physical normalization not certified",
    }
    jwrite(out / "fits" / (base.model_slug(model) + ".json"), result)
    print(f"{model}: FIT {status}, chi2={chi2:.9g}, actual evaluations={count}", flush=True)
    return result


def case_prediction(prepared, running: dict, args: argparse.Namespace,
                    *, sm_on: bool) -> dict:
    result, c5 = predict(prepared, running, args.target, args.wilson_factor,
                         sm_on=sm_on, check_scale=True)
    output = {
        "observables": extract_obs(result),
        "C5_at_muS": result["C5_at_muS"],
        "C5_at_muLow": result["smeft"]["C5_at_muLow_gev_inverse"],
        "direct_log_norm": result["direct_log_norm"],
        "scalar_bubble_norm": result["scalar_bubble_norm"],
        "C5_norm": result["C5_norm"],
        "scale_cancellation_residual": result["matching_scale_cancellation"],
        "beta_projection_residual": result["beta_projection_residual"],
        "bubble_projection_residual": result["bubble_projection_residual"],
        "scalar_mass2_at_S": result["scalar_mass2_at_S"],
        "native_source_hash": result["native_source_hash"],
    }
    return output


def stage_effects(cases: dict[str, dict]) -> dict:
    """One-stage-off relative to U1_E1_S1, using the thesis definition."""
    full = cases["U1_E1_S1"]["observables"]
    comparisons = {"UV": "U0_E1_S1", "EFT1": "U1_E0_S1", "SMEFT": "U1_E1_S0"}
    result = {}
    for stage, frozen in comparisons.items():
        old = cases[frozen]["observables"]
        result[stage] = {
            k: {"delta": float(full[k] - old[k]),
                "percent": (100.0 * (full[k]-old[k]) / old[k] if old[k] != 0 else None)}
            for k in OBS_KEYS
        }
    return result


def transfer_model(prepared, donor: dict, args: argparse.Namespace,
                   directory: Path) -> dict:
    overrides = apply_extra_yukawa_policy(prepared.overrides, args.extra_yukawas)
    y1 = base.stage1.finite_matrix(donor["y1"], "donor y1")
    y2 = base.stage1.finite_matrix(donor["y2"], "donor y2")
    overrides = put_yukawas(overrides, y1, y2)
    p = replace(prepared, overrides=overrides)
    uv_cache = {}
    cases = {}
    for uv_on, eft_on in ablation.CASE_KEYS:
        running = ablation._reconstruct_running(p, uv_on, eft_on, uv_cache)
        case = f"U{int(uv_on)}_E{int(eft_on)}_S1"
        cases[case] = case_prediction(p, running, args, sm_on=True)
        cases[case]["UV_nfev"] = running["uv_solver"]["nfev"]
        cases[case]["EFT1_nfev"] = running["eft1_solver"]["nfev"]
        if uv_on and eft_on:
            # The scalar-threshold C5 is *identical* with SMEFT on/off. Do not
            # repeat the expensive Matchete/HHSS response just to freeze SMEFT.
            c5_at_s = base.stage1.finite_matrix(cases[case]["C5_at_muS"], "C5 at S")
            frozen_physics = smeft_frozen(c5_at_s, running, p.config, args.target)
            cases["U1_E1_S0"] = {
                **cases[case],
                "observables": extract_obs(frozen_physics),
                "C5_at_muLow": frozen_physics["smeft"]["C5_at_muLow_gev_inverse"],
            }
    effects = stage_effects(cases)
    output = {
        "format": FORMAT, "kind": "shared_Yukawa_transfer",
        "model": p.model, "reference": donor["model"],
        "reference_fit_chi2": donor["final_chi2"],
        "original_model_fit_chi2": p.fit.get("chi2_final"),
        "status": "completed",
        "cases": cases,
        "effects": effects,
        "extra_yukawa_policy": args.extra_yukawas,
        "input_hashes": {
            "fitted_config": hash_path(p.fit_root / "fitted_config.json"),
            "target": hash_path(args.target),
        },
        "limitations": "Conditional soft-only EFT, spinor normalization/hard matching unverified",
    }
    jwrite(directory / "transfer.json", output)
    return output


def fit_phase(out: Path, args: argparse.Namespace):
    failures = []
    for model in DONORS:
        try:
            file = out / "fits" / (base.model_slug(model) + ".json")
            if file.exists():
                previous = jread(file)
                if previous.get("format") != FORMAT or previous.get("model") != model:
                    raise ValueError(f"Existing donor fit does not match run: {file}")
                print(f"{model}: previously fitted, reusing write-once output", flush=True)
                continue
            p = prepare(model, args)
            fit_donor(model, p, args, out)
        except Exception as exc:
            failures.append({"model": model, "error": str(exc),
                             "type": type(exc).__name__,
                             "traceback": traceback.format_exc()})
            print(f"{model}: FIT FAILED {type(exc).__name__}: {exc}", flush=True)
    jwrite(out / "fit_phase_summary.json", {
        "format": FORMAT, "failures": failures,
        "accepted_count": sum(jread(out / "fits" / (base.model_slug(m)+".json"))
                              ["status"] == "accepted_conditional_fit"
                              for m in DONORS
                              if (out/"fits"/(base.model_slug(m)+".json")).is_file()),
    })
    if failures:
        raise RuntimeError(f"Failed {len(failures)} donor fits; see fit_phase_summary.json")


def transfer_phase(out: Path, args: argparse.Namespace):
    donors = {}
    for model in DONORS:
        path = out / "fits" / (base.model_slug(model) + ".json")
        fit = jread(path)
        if fit.get("format") != FORMAT or fit.get("model") != model:
            raise ValueError(f"Invalid donor fit provenance: {path}")
        if fit.get("status") != "accepted_conditional_fit" and not args.include_poor:
            raise RuntimeError(
                f"{model}: chi2={fit.get('final_chi2')} was not accepted; "
                "do not generate thesis comparison as an accepted benchmark. "
                "Use --include-poor for clearly labelled provisional plots.")
        donors[model] = fit
    all_rows, effects_rows, errors = [], [], []
    for donor_model, donor in donors.items():
        donor_folder = out / "transfers" / base.model_slug(donor_model)
        for model in MODEL_ORDER:
            dest = donor_folder / base.model_slug(model)
            if dest.exists():
                print(f"{donor_model}->{model}: previously written; reading existing result", flush=True)
                current = jread(dest / "transfer.json")
            else:
                try:
                    p = prepare(model, args)
                    current = transfer_model(p, donor, args, dest)
                    print(f"{donor_model}->{model}: chi2={current['cases']['U1_E1_S1']['observables']['chi2']:.7g}", flush=True)
                except Exception as exc:
                    errors.append({"reference": donor_model, "model": model,
                                   "error": str(exc), "type": type(exc).__name__})
                    print(f"{donor_model}->{model}: FAILED {type(exc).__name__}: {exc}", flush=True)
                    continue
            if current.get("format") != FORMAT or current.get("model") != model or current.get("reference") != donor_model:
                errors.append({"reference": donor_model, "model": model, "error": "Existing transfer provenance mismatch"})
                continue
            if current.get("extra_yukawa_policy") != args.extra_yukawas:
                errors.append({"reference": donor_model, "model": model, "error": "Extra Yukawa policy differs from this run"})
                continue
            full = current["cases"]["U1_E1_S1"]["observables"]
            all_rows.append({"reference": donor_model, "model": model, **full})
            for stage, stage_info in current["effects"].items():
                effects_rows.append({"reference": donor_model, "model": model, "stage": stage,
                                     **{k+"_percent": stage_info[k]["percent"] for k in OBS_KEYS}})
    csvwrite(out / "transfer_observables.csv", all_rows)
    csvwrite(out / "running_effects.csv", effects_rows)
    jwrite(out / "transfer_phase_summary.json", {
        "format": FORMAT, "run_id": args.run_id,
        "donor_status": {m: f["status"] for m, f in donors.items()},
        "models_requested": len(DONORS)*len(MODEL_ORDER),
        "models_completed": len(all_rows), "models_failed": len(errors),
        "errors": errors, "extra_yukawa_policy": args.extra_yukawas,
        "running_stage_semantics": "UV, EFT1 are renormalisable RGE on/off; fixed-order LLSS log remains on. SMEFT frozen only in S0 case",
    })
    if errors:
        raise RuntimeError(f"{len(errors)} transfers failed; see transfer_phase_summary.json")


def parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--phase", choices=("fit", "transfer", "all"), required=True)
    p.add_argument("--run-id", required=True)
    p.add_argument("--campaign-root", type=Path, default=base.FIT_CAMPAIGN)
    p.add_argument("--canonical-run-id", default="full16_general_v1")
    p.add_argument("--beta-source", choices=("archived", "canonical", "auto"), default="archived")
    p.add_argument("--target", type=Path, default=ROOT/"configs/nufit_6_1_ic24_no.json")
    p.add_argument("--output-root", type=Path, default=OUT)
    p.add_argument("--wilson-factor", type=float, default=1.0)
    p.add_argument("--extra-yukawas", choices=("zero", "saved"), default="zero")
    p.add_argument("--max-nfev", type=int, default=60)
    p.add_argument("--max-yukawa", type=float, default=0.02)
    p.add_argument("--diff-step", type=float, default=1e-3)
    p.add_argument("--invalid-penalty", type=float, default=1e4)
    p.add_argument("--include-poor", action="store_true",
                   help="Allow poorly fitted donors to be transferred (labelled provisional)")
    p.add_argument("--check-only", action="store_true")
    return p


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", args.run_id):
        raise ValueError("Unsafe run-id")
    if not (args.max_nfev >= 1 and 0 < args.max_yukawa < 1 and
            0 < args.diff_step < 1 and args.invalid_penalty >= 10 and
            math.isfinite(args.wilson_factor) and args.wilson_factor != 0):
        raise ValueError("Invalid fitting or Wilson factor controls")
    args.campaign_root = base.absolute(args.campaign_root)
    args.target = base.absolute(args.target)
    args.output_root = base.absolute(args.output_root)
    if not args.target.is_file():
        raise FileNotFoundError(args.target)
    ablation._enable_mixed_hooks()
    if args.check_only:
        for model in MODEL_ORDER:
            p = prepare(model, args)
            print(f"{model}: READY, beta inventories {len(p.uv_beta['report_betas'])}/{len(p.eft_beta['report_betas'])}", flush=True)
        print("No fits, RGEs or files written (--check-only)")
        return 0
    out = args.output_root / args.run_id
    out.mkdir(parents=True, exist_ok=True)
    manifest = out / "manifest.json"
    if manifest.exists():
        existing = jread(manifest)
        if (existing.get("format") != FORMAT or
            existing.get("target_sha256") != hash_path(args.target) or
            existing.get("extra_yukawa_policy") != args.extra_yukawas or
            existing.get("beta_source") != args.beta_source or
            existing.get("wilson_factor") != args.wilson_factor or
            existing.get("canonical_run_id") != args.canonical_run_id or
            existing.get("campaign_root") != str(args.campaign_root)):
            raise RuntimeError("Run configuration differs from existing manifest; choose a new run-id")
    else:
        jwrite(manifest, {
            "format": FORMAT, "run_id": args.run_id,
            "campaign_root": str(args.campaign_root),
            "target": str(args.target), "target_sha256": hash_path(args.target),
            "beta_source": args.beta_source,
            "canonical_run_id": args.canonical_run_id,
            "wilson_factor": args.wilson_factor,
            "extra_yukawa_policy": args.extra_yukawas,
            "donors": list(DONORS), "models": list(MODEL_ORDER),
            "scales": {"UV": 1e7, "F": 1e5, "S": 1e3, "low": 1e2},
            "method": "Dynamic native F-tree + LLSS log + joint scalar bubble + SMEFT",
            "hard_matching_included": False,
            "physical_normalization_certified": False,
            "sensitivity_plot_regenerated": False,
        })
    if args.phase in ("fit", "all"):
        fit_phase(out, args)
    if args.phase in ("transfer", "all"):
        transfer_phase(out, args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
