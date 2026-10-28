"""Original interactive RGE and sensitivity oval reports for REDUCED scalar fits.

The two HTML renderers are imported, unchanged, from the existing pipeline:
  Numerical.plotting.InteractiveModelComparison.write_dashboard
  Numerical.plotting.BuildFittedModelEnvelopes.build

The new inputs are generated from the SAME reduced-scalar evaluator and matching
used by RunT3ReducedScalarSobolComparison; original lambda report inputs are
never substituted. No refitting, new symbolic matching, or UV-running fitting.

    python -m scripts.BuildReducedScalarExistingInteractiveReports \
        --case small_uv_frozen_5e-9 --mode optimised --model A_-2 --kind all
    python -m scripts.BuildReducedScalarExistingInteractiveReports --kind all
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CASES = ("small_uv_frozen_5e-9", "large_uv_frozen_0p5")
MODES = ("optimised", "shared_a_m2")
DIRECTIONS = ("T1", "T2", "T3")
OBS = ("delta_m21_sq_ev2", "delta_m3l_sq_ev2", "sin2_theta12",
       "sin2_theta13", "sin2_theta23")
SOURCE_RUN = "reduced_frozen_overnight_v2"
FORMAT = "t3-reduced-existing-interactive-diagnostics-v1"
DESCRIPTIONS = {
    "T1": "lambdaT3 -> gamma lambdaT3; y1,y2 -> gamma^(-1/2) y1,y2",
    "T2": "lambdaT3 -> gamma lambdaT3; y2 -> gamma^(-1) y2",
    "T3": "y1 -> gamma y1; y2 -> gamma^(-1) y2",
}


def require(test, message):
    if not test:
        raise ValueError("Reduced interactive reports: " + message)


def read(path):
    path = Path(path)
    require(path.is_file(), "missing " + str(path))
    return json.loads(path.read_text(encoding="utf-8-sig"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(content, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)


def validate_run_id(value):
    require(bool(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", value)),
            "invalid run ID")
    return value


def gamma_grid(low, high, points):
    require(math.isfinite(low) and math.isfinite(high) and
            0 < low < 1 < high and points >= 3 and points % 2 == 1,
            "gamma range must straddle 1 with an odd number of >=3 points")
    # Geometric symmetry about 1 ensures the original three gamma=1 centres.
    require(math.isclose(low * high, 1.0, rel_tol=1e-10, abs_tol=1e-13),
            "gamma min/max must be reciprocal")
    a = np.geomspace(low, high, points)
    a[points // 2] = 1.0
    return a.tolist()


def scale_inputs(vector, lam, direction, gamma):
    """The exact T1/T2/T3 directions of the original BenchmarkSensitivity."""
    x = np.asarray(vector, dtype=float).copy()
    require(x.shape == (18,) and np.isfinite(x).all() and gamma > 0,
            "invalid input vector or gamma")
    if direction == "T1":
        x *= gamma ** (-0.5)
        lam *= gamma
    elif direction == "T2":
        x[9:] *= gamma ** (-1.0)
        lam *= gamma
    elif direction == "T3":
        x[:9] *= gamma
        x[9:] /= gamma
    else:
        raise ValueError("Unknown sensitivity direction: " + str(direction))
    return x, float(lam)


def source_dir(run_id):
    return ROOT / "output/studies/reduced_scalar_sobol" / validate_run_id(run_id)


def diagnostics_root(run_id):
    return ROOT / "output/studies/reduced_scalar_interactive_diagnostics" / validate_run_id(run_id)


def report_root(run_id):
    return ROOT / "Reports/output/display_reduced_scalar" / validate_run_id(run_id)


def find_eft_beta(fit, model, reduced, bulk):
    """Reuse exact prior native EFT1 export recorded by this accepted fit."""
    candidates = []
    for relative, hash_value in (fit.get("source_sha256") or {}).items():
        p = ROOT / relative
        if p.name != bulk.safe(model) + ".json" or "/betas/" not in ("/" + str(relative).replace("\\", "/")):
            continue
        if not p.is_file():
            continue
        require(sha(p) == hash_value, "modified EFT1 source: " + str(p))
        candidates.append(p)
    require(len(candidates) == 1,
            model + ": expected one verified saved native EFT1 beta export, got " + str(candidates))
    payload = bulk.read(candidates[0])
    require(payload.get("model") == model and payload.get("status") == "Success",
            model + ": invalid native EFT1 export")
    return payload, candidates[0]


def load_context(case, model, run_id):
    from scripts import RunT3GeneralScalar48Campaign as bulk
    from scripts import RunT3ReducedScalarSobolComparison as reduced
    from Numerical.orchestration import RunSmallLambdaRefit as refit
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    require(case in CASES, "unsupported frozen case " + str(case))
    registry = bulk.read(bulk.REGISTRY)
    require(len(registry["models"]) == 16 and model in registry["models"],
            "unsupported T3 model " + model)
    entry, new, reused = bulk.inventory(registry, model)
    path = source_dir(run_id) / case / bulk.safe(model) / "fit.json"
    fit = read(path)
    require(fit.get("model") == model and fit.get("case") == case and
            fit.get("status") == "Accepted" and
            fit.get("format") == reduced.FORMAT and
            fit.get("extra_quadratic_boundary_exact_zero") is True and
            fit.get("additional_quadratic_evolution_included") is False and
            set(fit.get("new_additional_quartic_symbols", [])) == set(new) and
            set(fit.get("reused_additional_quartic_symbols_zeroed", [])) == set(reused),
            model + ": saved result does not match reduced-scalar convention")
    beta, beta_file = find_eft_beta(fit, model, reduced, bulk)
    hist = bulk.historical(case, model, entry)
    require(hist["uv_on"] is False, "UV-running source passed to frozen dashboard")
    target = OscillationFitTarget.from_json(refit.TARGET)
    kernel, kernel_file, kind = reduced.matching(model)
    require(fit["c5_source_kind"] == kind and
            (ROOT / fit["c5_source"]).resolve() == kernel_file.resolve(),
            "C5 matching source differs from completed fit")
    donor_file = source_dir(run_id) / case / bulk.safe("A_-2") / "fit.json"
    donor = read(donor_file)
    require(donor.get("status") == "Accepted", "A(-2) reference not accepted")
    return {"case": case, "model": model, "entry": entry, "new": new,
            "reused": reused, "fit": fit, "fit_path": path,
            "donor": donor, "donor_path": donor_file,
            "hist": hist, "beta": beta, "beta_file": beta_file,
            "kernel": kernel, "kernel_file": kernel_file, "target": target,
            "bulk": bulk, "reduced": reduced,
            "long_name": bulk.stem(entry).removesuffix(".json")}


def chosen_parameters(ctx, mode):
    require(mode in MODES, "unknown mode")
    block = (ctx["fit"] if mode == "optimised" else ctx["donor"]).get("best_mass_safe") or {}
    names = ctx["bulk"].NAMES
    require(all(k in block.get("parameters", {}) for k in names), "missing saved Yukawa entries")
    x = np.array([block["parameters"][n] for n in names], dtype=float)
    require(x.shape == (18,) and np.isfinite(x).all(), "invalid saved fitted Yukawas")
    target_block = block if mode == "optimised" else ctx["fit"].get("shared_a2", {})
    require(target_block.get("chi2") is not None and target_block.get("prediction"),
            "no saved evaluation for selected Yukawa mode")
    return x, target_block


def is_close(a, b, *, atol=1e-9, rtol=2e-5):
    return math.isclose(float(a), float(b), abs_tol=atol, rel_tol=rtol)


def check_endpoint(result, reference, *, where):
    require(result.get("status") == "NeutrinoEvaluated", where + ": baseline failed screen")
    require(is_close(result["chi2"], reference["chi2"], atol=2e-5),
            where + ": endpoint chi2 differs from stored fit")
    for key in OBS:
        require(is_close(result["prediction"][key], reference["prediction"][key],
                         atol=1e-11, rtol=3e-5), where + ": endpoint " + key + " differs")


def gamma_hist(original, lam):
    h = dict(original)
    h["cfg"] = copy.deepcopy(original["cfg"])
    value = h["cfg"]["base_state"]["ordinary"]["lambdaT3"]
    require(isinstance(value, dict) and "real" in value and "imag" in value,
            "unsupported lambdaT3 serialization")
    value["real"] = float(lam)
    value["imag"] = 0.0
    h["lambda_fixed"] = float(lam)
    return h


# The frozen reduced-scalar evaluator integrates EFT1 in __init__, before
# ReducedEvaluator.evaluate can catch and mark a failed gamma point. For large
# lambdaT3, that boundary RGE can fail (e.g. an ODE step size underflow).
# These are genuine unavailable trial points, not missing source artifacts.
# Never turn unexpected metadata/matching/validation errors into scan failures.
EXPECTED_SCALAR_RGE_FAILURES = (
    "Coupled EFT1 integration failed:",
    "Nonfinite coupled EFT1 beta",
    "must be positive to define a physical scalar threshold",
)


def is_expected_scalar_rge_failure(exc):
    return (isinstance(exc, ValueError) and
            any(marker in str(exc) for marker in EXPECTED_SCALAR_RGE_FAILURES))


def evaluate_gamma(ctx, vector, gamma, direction, cache):
    x, lam = scale_inputs(vector, ctx["hist"]["lambda_fixed"], direction, gamma)
    if lam not in cache:
        h = gamma_hist(ctx["hist"], lam)
        try:
            cache[lam] = ctx["reduced"].ReducedEvaluator(
                ctx["case"], ctx["model"], ctx["entry"], ctx["new"], ctx["reused"],
                h, ctx["beta"], ctx["kernel"], ctx["target"])
        except ValueError as exc:
            if not is_expected_scalar_rge_failure(exc):
                raise
            # The gamma=1 reference is mandatory and cannot be skipped.
            if math.isclose(float(gamma), 1.0, rel_tol=0, abs_tol=1e-12):
                raise
            cache[lam] = {"status": "NumericalFailure", "reason": str(exc)[:500]}
    item = cache[lam]
    if isinstance(item, dict):
        return dict(item)
    return item.evaluate(x, f"sensitivity_{direction}_{gamma:.9g}")


def target_record(target):
    return {"ordering": target.ordering,
            "observable_names": list(target.observable_names),
            "central_values": np.asarray(target.central_values, dtype=float).tolist(),
            "covariance": np.asarray(target.covariance, dtype=float).tolist(),
            "source": target.source,
            "note": "Gaussianised NuFIT fit target, not official joint NuFIT contours"}


def sensitivity_record(ctx, mode, gamma_values):
    vector, reference = chosen_parameters(ctx, mode)
    cache = {}
    directions = {}
    for direction in DIRECTIONS:
        points = []
        for gamma in gamma_values:
            row = evaluate_gamma(ctx, vector, gamma, direction, cache)
            point = {"gamma": float(gamma), "status": row["status"]}
            if row["status"] == "NeutrinoEvaluated":
                point.update(chi2=float(row["chi2"]), prediction=row["prediction"])
                if is_close(gamma, 1., rtol=0, atol=1e-12):
                    check_endpoint(row, reference, where=ctx["model"]+"/"+mode)
            else:
                point["message"] = str(row.get("reason", "screen failed"))
            points.append(point)
        directions[direction] = {"definition": DESCRIPTIONS[direction], "points": points}
    centres = [directions[d]["points"][len(gamma_values)//2] for d in DIRECTIONS]
    require(all(c["status"] == "NeutrinoEvaluated" for c in centres),
            ctx["model"] + ": failed gamma=1 centre")
    for c in centres[1:]:
        for key in OBS:
            require(is_close(c["prediction"][key], centres[0]["prediction"][key],
                             atol=1e-12, rtol=1e-8),
                    ctx["model"] + ": inconsistent T1/T2/T3 gamma=1 centre")
    return {"format": FORMAT, "status": "Success", "case": ctx["case"],
            "mode": mode, "model": ctx["model"], "source_fit_sha256": sha(ctx["fit_path"]),
            "source_A_m2_fit_sha256": sha(ctx["donor_path"]),
            "gamma_min": float(gamma_values[0]), "gamma_max": float(gamma_values[-1]),
            "n_points": len(gamma_values), "spacing": "log",
            "n_failed_points": sum(p["status"] != "NeutrinoEvaluated"
                                   for d in DIRECTIONS for p in directions[d]["points"]),
            "directions": directions, "target": target_record(ctx["target"]),
            "notes": {"uv_running": "Disabled", "eft1": "Candidate extra-quartic 1-loop RGEs",
                      "C5": "Existing fixed-order leading one-loop matching; legacy cases not harmonised",
                      "ovals": "Empirical sensitivity envelopes; NOT confidence regions"}}


def _low_running(ctx, uv, eft, c5, vev, n_points):
    """One continuous SM+Weinberg solve; never interpolate a fitted observable."""
    from scipy.integrate import solve_ivp
    from Numerical.running.ScalarThresholdBoundary import (
        build_final_sm_boundary, build_sm_weinberg_initial_conditions)
    from Numerical.running.SMWeinbergEvolution import _beta, _pack, _unpack
    from Numerical.fitting.ParameterScan import (
        neutrino_mass_from_c5, charged_lepton_mass_basis_matrix)
    from Numerical.fitting.OscillationFit import evaluate_mass_matrix_fit
    from physics.NeutrinoObservables import calculate_neutrino_observables
    mu_s, mu_low = 1e3, 1e2
    boundary = build_final_sm_boundary(eft, matching_scale_gev=mu_s)
    ini = build_sm_weinberg_initial_conditions(boundary, c5)
    mus = np.geomspace(mu_s, mu_low, n_points)
    solution = solve_ivp(_beta, (math.log(mu_s), math.log(mu_low)), _pack(ini),
                         t_eval=np.log(mus), method="DOP853", rtol=1e-8, atol=1e-11)
    require(solution.success and np.isfinite(solution.y).all(),
            "SM+Weinberg sampled evolution failed")
    curves = {key: [] for key in OBS}
    c5_abs = []
    masses = []
    for item in solution.y.T:
        _, _, _, _, ye, _, _, k = _unpack(item)
        k = (k + k.T)/2
        neutrino_mass = neutrino_mass_from_c5(k, vev_gev=vev)
        mflavour = charged_lepton_mass_basis_matrix(neutrino_mass, ye)
        fit = evaluate_mass_matrix_fit(mflavour, ctx["target"])
        prediction = dict(zip(fit.observable_names, map(float, fit.prediction)))
        for name in OBS:
            curves[name].append(float(prediction[name]))
        physics = calculate_neutrino_observables(mflavour,
                                                  ordering=ctx["target"].ordering)
        eigen = np.asarray(physics.masses_ev, dtype=float)
        require(eigen.shape == (3,) and np.isfinite(eigen).all(),
                "invalid neutrino masses in sampled Weinberg evolution")
        masses.append(eigen.tolist())
        c5_abs.append(np.abs(k).tolist())
    return {"mu_gev": [float(m) for m in mus], "masses_ev": masses,
            "c5_abs": c5_abs, **curves}


def rge_record(ctx, mode, n_points):
    from Numerical.core.StateVector import unpack_uv_state
    from scripts import Stage31CoupledUVRunning as s31
    vector, reference = chosen_parameters(ctx, mode)
    bulk = ctx["bulk"]
    cfg = bulk.update_yukawas(ctx["hist"], vector)
    packed, layout = s31.prepare_state(cfg, ctx["entry"], ctx["reused"],
                                       {name: 0j for name in ctx["reused"]})
    uv = unpack_uv_state(packed, layout, mu_gev=1e5)
    eft, masses, qs, _ = bulk.run_eft(uv, tuple(0j for _ in ctx["new"]),
                                      ctx["beta"], ctx["model"], ctx["new"], ctx["reused"])
    screen = bulk.testscreen(uv, eft, masses, tuple(0j for _ in ctx["new"]), qs)
    require(screen["screen"], ctx["model"] + ": saved endpoint fails physical screen")
    c5 = ctx["kernel"].evaluate(uv, eft)
    from scripts.RunT3GeneralScalarBCExploratoryFit import low_energy_fit
    endpoint, _ = low_energy_fit(c5, eft, ctx["target"],
                                  vev=float(cfg.get("numerical", {}).get("vev_gev", 246.22)))
    prediction = dict(zip(endpoint.observable_names, map(float, endpoint.prediction)))
    check_endpoint({"status": "NeutrinoEvaluated", "chi2": float(endpoint.chi2),
                    "prediction": prediction}, reference, where=ctx["model"]+"/"+mode)
    final = _low_running(ctx, uv, eft, c5,
                         float(cfg.get("numerical", {}).get("vev_gev", 246.22)), n_points)
    for name in OBS:
        require(is_close(final[name][-1], prediction[name], atol=1e-11, rtol=3e-5),
                "sampled Weinberg endpoint mismatch for " + name)
    return {"format": FORMAT, "status": "Success", "source": "reduced-scalar-frozen-RGE-sampled",
            "case": ctx["case"], "mode": mode, "model": ctx["model"],
            "source_fit_sha256": sha(ctx["fit_path"]),
            "source_A_m2_fit_sha256": sha(ctx["donor_path"]),
            "scales_gev": {"mu_fermion_threshold": 1e5,
                           "mu_scalar_threshold": 1e3, "mu_low": 1e2},
            "low_energy": {"chi2": float(endpoint.chi2), "prediction": prediction},
            "final_running": final,
            "notes": {"UV": "Frozen at muF: no physical UV running trajectory",
                      "EFT1": "Renormalisable scalar-only running included in endpoint; not plotted above muS",
                      "C5": "Fixed-order leading one-loop coefficient; no fabricated intermediate C5",
                      "scalar_masses_muS_tev": screen["scalar_masses_muS_tev"],
                      "matching": ctx["fit"]["c5_source_kind"]}}


def reuse_compatible(path, ctx, mode, *, n_points=None, gamma_values=None):
    if not Path(path).is_file():
        return False
    try:
        p = read(path)
        if (p.get("format") != FORMAT or p.get("status") != "Success" or
            p.get("mode") != mode or p.get("model") != ctx["model"] or
            p.get("case") != ctx["case"] or
            p.get("source_fit_sha256") != sha(ctx["fit_path"]) or
            p.get("source_A_m2_fit_sha256") != sha(ctx["donor_path"])):
            return False
        if n_points is not None:
            return len(p.get("final_running", {}).get("mu_gev", [])) == n_points
        if gamma_values is not None:
            return (p.get("n_points") == len(gamma_values) and
                    is_close(p.get("gamma_min"), gamma_values[0]) and
                    is_close(p.get("gamma_max"), gamma_values[-1]))
        return True
    except (ValueError, TypeError, KeyError, json.JSONDecodeError):
        return False


def render_existing(source, destination, mode, case, kind):
    from Numerical.plotting.InteractiveModelComparison import write_dashboard
    from Numerical.plotting.BuildFittedModelEnvelopes import build as build_ovals
    output = destination / case / "reused" / mode
    output.mkdir(parents=True, exist_ok=True)
    if kind in ("rge", "all"):
        write_dashboard(source / "models", output / "interactive_rge.html",
                        f"T3 reduced scalar — UV off — {case} — {mode}")
    if kind in ("sensitivity", "all"):
        items = []
        for model_folder in sorted((source / "sensitivity").glob("*")):
            if not model_folder.is_dir():
                continue
            scan_path = model_folder / "sensitivity_diagnostics.json"
            if not scan_path.is_file():
                continue
            scan = read(scan_path)
            centre = scan["directions"]["T1"]["points"][scan["n_points"]//2]
            items.append({"model": model_folder.name, "status": "Success",
                          "verified_chi2": float(centre["chi2"]),
                          "fit_stage": "reduced_scalar_frozen_" + mode})
        require(bool(items), "no saved sensitivity scans to render")
        write_json(source / "sensitivity" / "comparison_summary.json", items)
        build_ovals(source / "sensitivity", output / "sensitivity_ovals.html")
    return output


def run(args):
    from scripts import RunT3GeneralScalar48Campaign as bulk
    chosen_cases = CASES if args.case == "all" else (args.case,)
    chosen_modes = MODES if args.mode == "all" else (args.mode,)
    registry = bulk.read(bulk.REGISTRY)
    models = list(registry["models"])
    if args.model:
        selected = [model for model in models if args.model in (
            model, bulk.safe(model), bulk.stem(registry["models"][model]).removesuffix(".json"))]
        require(len(selected) == 1, "model not found: " + args.model)
        models = selected
    gamma_values = gamma_grid(args.gamma_min, args.gamma_max, args.points)
    source = diagnostics_root(args.source_run)
    dest = report_root(args.source_run)
    completed = []
    for case in chosen_cases:
        for mode in chosen_modes:
            current = (source / case / mode if not args.model else
                       source / "pilots" / bulk.safe(args.model) / case / mode)
            failures = []
            for model in models:
                ctx = load_context(case, model, args.source_run)
                if args.kind in ("rge", "all"):
                    target = current / "models" / ctx["long_name"] / "running_diagnostics.json"
                    if not args.force and reuse_compatible(target, ctx, mode, n_points=args.rge_points):
                        print("REUSE RGE", case, mode, model, flush=True)
                    else:
                        print("CALCULATE RGE", case, mode, model, flush=True)
                        try:
                            write_json(target, rge_record(ctx, mode, args.rge_points))
                        except Exception as exc:
                            failures.append((model, "rge", type(exc).__name__ + ": " + str(exc)))
                            print("FAILED RGE", case, mode, model, str(exc), flush=True)
                if args.kind in ("sensitivity", "all"):
                    target = current / "sensitivity" / ctx["long_name"] / "sensitivity_diagnostics.json"
                    if not args.force and reuse_compatible(target, ctx, mode, gamma_values=gamma_values):
                        print("REUSE SENSITIVITY", case, mode, model, flush=True)
                    else:
                        print("CALCULATE SENSITIVITY", case, mode, model, flush=True)
                        try:
                            scan = sensitivity_record(ctx, mode, gamma_values)
                            write_json(target, scan)
                            print("SENSITIVITY SAVED", case, mode, model,
                                  "failed gamma points:", scan["n_failed_points"], flush=True)
                        except Exception as exc:
                            failures.append((model, "sensitivity", type(exc).__name__ + ": " + str(exc)))
                            print("FAILED SENSITIVITY", case, mode, model, str(exc), flush=True)
            # Never generate a deceptively complete all-model dashboard for a pilot selection.
            # In pilot mode, write into separate pilot subdirectory.
            output_dir = dest if not args.model else dest / "pilots" / bulk.safe(args.model)
            if not failures:
                try:
                    rendered = render_existing(current, output_dir, mode, case, args.kind)
                    completed.append(str(rendered))
                    print("REPORTS", rendered, flush=True)
                except Exception as exc:
                    failures.append(("ALL", "renderer", type(exc).__name__ + ": " + str(exc)))
                    print("FAILED DASHBOARD", str(exc), flush=True)
            if failures:
                print("ERRORS", case, mode, failures, flush=True)
                raise RuntimeError(f"{case}/{mode}: {len(failures)} errors; see log")
    return completed


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-run", default=SOURCE_RUN)
    p.add_argument("--case", choices=(*CASES, "all"), default="all")
    p.add_argument("--mode", choices=(*MODES, "all"), default="all")
    p.add_argument("--kind", choices=("rge", "sensitivity", "all"), default="all")
    p.add_argument("--model", help="Model such as A_-2 or T3_dS1_1_dS2_3_dF_2_alpha_m2")
    p.add_argument("--gamma-min", type=float, default=.1)
    p.add_argument("--gamma-max", type=float, default=10.)
    p.add_argument("--points", type=int, default=11)
    p.add_argument("--rge-points", type=int, default=41)
    p.add_argument("--force", action="store_true")
    args = p.parse_args()
    require(args.rge_points >= 3, "rge-points must be >=3")
    validate_run_id(args.source_run)
    completed = run(args)
    for item in completed:
        print("DONE", item)


if __name__ == "__main__":
    main()
