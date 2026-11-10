#!/usr/bin/env python3
"""Thesis-input pilot for A_-4 and D_-2 using the existing T3 strict-one-loop chain.

This is an adapter/driver, NOT a replacement of EFT formulae. Every objective
call evaluates fresh UV and EFT1 background RGEs, source-resolved F LLSS/hard
triangle, strict EFT1 Wilson, simultaneous scalar matching, and SMEFT running.
Only the analytic_triangle prescription is used for the first two-model refit.

Local-only data. No modification of GitHub, Matchete/RGBeta source, or the
completed original Batch 1--8 results. Candidate scratch outputs are isolated
under output/parallel_t3/ and removed unless --keep-work is selected.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
from dataclasses import dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
import shutil
import sys
import traceback
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

CONFIG = Path("configs/thesis_numerics/t3_thesis_A4_Dm2_pilot.json")
OUTPUT = Path("output/parallel_t3/thesis_A4_Dm2_pilot_v1")
MODELS = ("A_-4", "D_-2")
ORDER = ("delta_m21_sq_ev2", "delta_m3l_sq_ev2", "sin2_theta12", "sin2_theta13", "sin2_theta23")
SUCCESS_B1 = "NUMERICAL_FHARD_TRIANGLE_CANDIDATE_NOT_CERTIFIED"
SUCCESS_B2 = "NUMERICAL_STRICT_ONE_LOOP_INCLUDED_TOPOLOGIES_CANDIDATE_NOT_CERTIFIED"
ADAPTER_REVISION = "v2_source_beta_key_override_fix"
FIT_REVISION = "v3_sobol_screened_bounded_least_squares"


def fail(message: str) -> None:
    raise ValueError("Thesis pilot: " + message)


def read_json(path: Path) -> dict:
    val = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(val, dict):
        fail("Expected JSON object: " + str(path))
    return val


def safe_parent(root: Path, path: Path) -> Path:
    root = root.resolve()
    target = (root / path).resolve()
    if not target.is_relative_to(root / "output" / "parallel_t3"):
        fail("Outputs must be inside the project's output/parallel_t3/ directory")
    return target


def _sha(value: Any) -> str:
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def load_config(root: Path, path: Path) -> dict:
    full = root / path if not path.is_absolute() else path
    cfg = read_json(full)
    if cfg.get("schema") != "t3-thesis-benchmark-two-model-pilot-v1":
        fail("Unexpected benchmark schema")
    if set(cfg.get("models", {})) != set(MODELS):
        fail("This first pilot is restricted to A_-4 and D_-2")
    scales = cfg.get("scales_GeV", {})
    if scales != {"UV": 1e7, "F": 1e5, "S": 1e3, "low": 100.0}:
        fail("Thesis pilot scale hierarchy changed unexpectedly")
    if cfg["lambdaT3_UV"] != 0.1 or cfg["lambda_other_UV"] != 0.01:
        fail("The thesis T3 scalar-coupling benchmark must be preserved")
    if cfg["class_E_extra_half_factor"] is not False:
        fail("Class-E half multiplier must be absent")
    for model in MODELS:
        for name in ("y1_UV_times_1e4", "y2_UV_times_1e4"):
            mat = np.asarray(cfg["models"][model][name], dtype=float)
            if mat.shape != (3, 3) or not np.isfinite(mat).all():
                fail(model + ": invalid appendix Yukawa matrix " + name)
    return cfg


def initial_parameters(cfg: dict, model: str) -> np.ndarray:
    scale = float(cfg["fitting"]["input_Yukawa_scale"])
    return np.concatenate([
        np.asarray(cfg["models"][model]["y1_UV_times_1e4"], float).ravel(),
        np.asarray(cfg["models"][model]["y2_UV_times_1e4"], float).ravel()
    ]) * scale


def _complex_items(m: np.ndarray) -> list[list[dict[str, float]]]:
    arr = np.asarray(m, complex)
    if arr.shape != (3, 3) or not np.isfinite(arr).all():
        fail("Yukawa matrix must be finite complex 3x3")
    return [[{"real": float(x.real), "imag": float(x.imag)} for x in row] for row in arr]


def overrides_for(cfg: dict, model: str, names: list[str], vec: np.ndarray, prior: dict) -> dict:
    """Use thesis Table 4.1 for specified values; preserve existing SM-Y defaults.

    Every beta symbol in 'names' is already validated by the upstream runner.
    Extra quartics are 0.1*lambdaT3=0.01; unspecified extra bilinears are
    explicitly set to zero, rather than accidentally inherited from development.
    """
    v = np.asarray(vec, dtype=float)
    if v.shape != (18,) or not np.isfinite(v).all():
        fail("Expected 18 finite real Yukawa entries")
    if np.max(abs(v)) >= math.sqrt(4 * math.pi):
        fail("BSM Yukawas exceed perturbative search bound")
    # Workstream B's _bench_overrides is intentionally SPARSE: it only supplies
    # a subset of the UV beta symbols. GeneralRunningBenchmark.initial_state
    # fills the remaining source-recognised coordinates from its initializer.
    # The thesis adapter explicitly provides the scalar coordinates needed for
    # the thesis benchmark; do not compare against the sparse *prior* keys.
    known = set(names)
    if len(known) != len(names):
        fail("Duplicate source UV beta keys")
    if not set(prior).issubset(known):
        fail("Workstream B supplied unknown UV beta keys: " + repr(sorted(set(prior) - known)))
    required = {"MF", "gY", "g2", "g3", "lambdaH", "lambdaT3",
                "mS1Sq", "mS2Sq", "y1", "y2"}
    if not required.issubset(known):
        fail("Source UV beta inventory missing thesis benchmark fields: " + repr(sorted(required - known)))
    result = dict(prior)
    f = cfg["heavy_fermion_masses_UV_GeV"]
    result["MF"] = _complex_items(np.diag(np.asarray(f, float)))
    result.update({k: float(val) for k, val in cfg["gauge_UV"].items()})
    result["lambdaH"] = float(cfg["lambdaH_UV"])
    result.update({k: float(val) for k, val in cfg["scalar_masses_squared_UV_GeV2"].items()})
    for name in names:
        if name.startswith("lambda") and name not in ("lambdaH", "lambdaT3"):
            result[name] = float(cfg["lambda_other_UV"])
        elif name.startswith("muT3General"):
            result[name] = {"real": 0., "imag": 0.}
        elif name in ("y1Conj", "y2Conj"):
            result[name] = _complex_items(np.zeros((3, 3)))
    result["lambdaT3"] = float(cfg["lambdaT3_UV"])
    result["y1"] = _complex_items(v[:9].reshape(3, 3))
    result["y2"] = _complex_items(v[9:].reshape(3, 3))
    # B's existing SM flavour input is not invented here: its _bench_overrides
    # provides yu/yd/ye until a source-exact thesis SM input file is supplied.
    unknown = set(result) - known
    if unknown:
        fail("Thesis benchmark supplies unsupported UV beta keys: " + repr(sorted(unknown)))
    return result


@contextmanager
def patched_attributes(patches: list[tuple[object, str, Any]]):
    original = [(module, name, getattr(module, name)) for module, name, _ in patches]
    try:
        for module, name, value in patches:
            setattr(module, name, value)
        yield
    finally:
        for module, name, value in reversed(original):
            setattr(module, name, value)


def _get_prediction(record: dict) -> dict:
    if not record.get("mixing_defined"):
        fail("New EFT calculation has undefined neutrino mixing: " + str(record.get("mixing_blocker")))
    result = {
        "delta_m21_sq_ev2": float(record["dm21_sq_eV2"]),
        "delta_m3l_sq_ev2": float(record["dm31_sq_eV2"]),
        "sin2_theta12": float(record["sin2_theta12"]),
        "sin2_theta13": float(record["sin2_theta13"]),
        "sin2_theta23": float(record["sin2_theta23"]),
    }
    if not all(math.isfinite(x) for x in result.values()):
        fail("Nonfinite oscillation observables")
    return result


def residuals(prediction: dict, target: dict) -> np.ndarray:
    return np.array([(prediction[key]-target["central_values"][key]) /
                     target["one_sigma_errors"][key] for key in ORDER], dtype=float)


def _run_full_chain(root: Path, cfg: dict, model: str, vec: np.ndarray, trial_hash: str,
                    *, keep_work: bool) -> dict:
    """Execute checked existing source modules; DO NOT recompute from archived C5."""
    from Numerical.parallel_t3.background import runner as b
    from Numerical.parallel_t3.rough_inputs import builder
    from Numerical.parallel_t3.wilson import fermion_hard_triangle_batch1 as fermion_hard
    from Numerical.parallel_t3.wilson import quartic_channels_batch2 as qchannels
    from scripts import RunT3FHardTriangleBatch1 as stage1
    from scripts import RunT3StrictOneLoopBatch2 as stage2
    from Numerical.parallel_t3.background.contracts import tag

    token = model.replace("_", "").replace("-", "m").replace("+", "p") + "_" + trial_hash[:12]
    b_run = "ThesisB_" + token
    f_run = "ThesisF_" + token
    scratch = root / "output/parallel_t3/thesis_A4_Dm2_pilot_v1/_scratch" / token
    # Each trial owns its *own* B background, F rough inputs, Batch1 and Batch2.
    b_path = root / "output/parallel_t3/background" / b_run / tag(model)
    f_path = root / "output/parallel_t3/rough_inputs" / f_run
    stage1_out = scratch / "batch1"
    stage2_out = scratch / "batch2"
    paths = (b_path, f_path, scratch)
    if any(p.exists() for p in paths):
        fail("A trial's scratch data already exist. Inspect before retry: " + token)
    options = b.Options(root=root, run_id=b_run,
                        mu_uv=cfg["scales_GeV"]["UV"],
                        mu_f=cfg["scales_GeV"]["F"],
                        mu_s=cfg["scales_GeV"]["S"])
    original_bench = b._bench_overrides

    def candidate_overrides(name: str, keys: list[str], opt: Any) -> dict:
        if name != model:
            fail("Unexpected model in the two-model trial")
        base = original_bench(name, keys, opt)
        return overrides_for(cfg, model, keys, vec, base)

    original_rough = stage1.run_rough_case

    def rough_at_100(*a: Any, **kw: Any) -> Any:
        kw["mu_low"] = 100.
        kw["vev"] = float(cfg["vev_GeV"])
        return original_rough(*a, **kw)

    patches = [
        (b, "_bench_overrides", candidate_overrides),
        (builder, "B_RUN", b_run),
        (fermion_hard, "B_RUN", b_run),
        (qchannels, "B_RUN", b_run),
        (stage1, "run_rough_case", rough_at_100),
        (stage1, "RUN_ID", f_run),
        (stage1, "INPUT", Path("output/parallel_t3/rough_inputs") / f_run),
        (stage2, "BASE_INPUT", Path("output/parallel_t3/rough_inputs") / f_run / "a"),
    ]
    try:
        with patched_attributes(patches):
            br = b.run_one(options, model)
            if br.get("background") != "CandidateRenormalisableBackgroundSuccess" or not br.get("fermion_boundary", "").startswith("FormalSourceExact"):
                fail("Fresh B/background or fermion boundary blocked: " + repr(br))
            a = stage1.execute(root, output=stage1_out, models=(model,))
            one = a.get("results", {}).get(model, {})
            if one.get("status") != SUCCESS_B1:
                fail("Fresh source-resolved hard triangle blocked: " + repr(one.get("blocker", one.get("stage"))))
            c = stage2.run_batch(root, batch1=stage1_out, output=stage2_out,
                                 models=(model,), mu_low=100., vev=float(cfg["vev_GeV"]))
            two = c.get("results", {}).get(model, {})
            if two.get("status") != SUCCESS_B2:
                fail("Strict EFT1/two-scalar/SMEFT stage blocked: " + repr(two.get("blocker", two.get("stage"))))
            obs = read_json(stage2_out / "results" / tag(model) / "observables.json")
            pred = _get_prediction(obs)
            # The saved manifest guarantees matching at the common S threshold.
            return {"prediction": pred,
                    "masses_eV": [float(x) for x in obs["masses_eV"]],
                    "mixing_defined": True,
                    "physics_status": SUCCESS_B2,
                    "EFT_order": "UV -> ALL F -> EFT1(S1,S2) -> BOTH S -> SMEFT",
                    "scenario": "analytic_triangle_strict_one_loop",
                    "scalar_quartic_scope": "existing_reduced_source_tensor; unmapped_extra_quartics omitted from LLHH",
                    "mu_low_GeV": 100.,
                    "scratch_kept": bool(keep_work),
                    "source_records": {"B": str((b_path / 'background_manifest.json').relative_to(root)),
                                       "F_hard": str((stage1_out / 'summary.json').relative_to(root)),
                                       "EFT1_SMEFT": str((stage2_out / 'summary.json').relative_to(root))}}
    finally:
        if not keep_work:
            for path in paths:
                if path.is_dir():
                    shutil.rmtree(path)


@dataclass
class Evaluator:
    root: Path
    config: dict
    model: str
    target: dict
    keep_work: bool
    limit: int
    accepted: list[dict]
    attempts: int = 0

    def output_dir(self) -> Path:
        return self.root / OUTPUT / self.model.replace("-", "m").replace("+", "p")

    def evaluate(self, x: np.ndarray) -> dict:
        x = np.asarray(x, float)
        if x.shape != (18,) or not np.isfinite(x).all():
            fail("Invalid optimiser candidate")
        payload = {"schema": self.config["schema"], "cfg_digest": _sha(self.config),
                   "adapter_revision": ADAPTER_REVISION,
                   "model": self.model, "x": x.tolist(), "mu_low_GeV": 100.}
        digest = _sha(payload)
        path = self.output_dir() / "evaluations" / (digest + ".json")
        if path.is_file():
            data = read_json(path)
            if data.get("candidate_digest") != digest or data.get("candidate") != x.tolist():
                fail("Cached candidate does not match requested input")
            # A blocked candidate may have failed due to a corrected adapter or
            # temporarily missing source data. Never treat old BLOCKED as final.
            if data.get("status") == "COMPLETE_CONDITIONAL":
                return data
            print(self.model, "retrying previously blocked candidate", digest[:12], flush=True)
        if self.attempts >= self.limit:
            raise RuntimeError("Pilot evaluation budget reached")
        self.attempts += 1
        result = {"candidate_digest": digest, "candidate": x.tolist(),
                  "model": self.model, "config_digest": _sha(self.config),
                  "mu_low_GeV": 100., "classification": "CONDITIONAL_NOT_PHYSICAL_MSbar",
                  "status": "BLOCKED"}
        try:
            values = _run_full_chain(self.root, self.config, self.model, x, digest,
                                     keep_work=self.keep_work)
            err = residuals(values["prediction"], self.target)
            result.update(status="COMPLETE_CONDITIONAL", chi2_approx=float(err @ err),
                          residuals=err.tolist(), **values)
        except Exception as exc:
            result.update(error=f"{type(exc).__name__}: {exc}",
                          traceback_tail=traceback.format_exc(limit=5))
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        if result["status"] == "COMPLETE_CONDITIONAL":
            self.accepted.append(result)
        print(self.model, "eval", self.attempts, result["status"],
              result.get("chi2_approx", result.get("error")), flush=True)
        return result


def _trial_box(start: np.ndarray, perturbation: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Finite, symmetric trial bounds around the thesis seed for all 18 real entries.

    The 2e-5 absolute floor allows small/zero seed entries to move. This box
    is a numerical search choice, not a physical restriction on allowed Yukawas.
    """
    cap = math.sqrt(4 * math.pi) * (1 - 1e-8)
    width = perturbation * np.maximum(np.abs(start), 2e-5)
    lower = np.maximum(start-width, -cap)
    upper = np.minimum(start+width, cap)
    if not np.all(upper > lower):
        fail("Degenerate Yukawa optimisation bounds")
    return lower, upper, (upper-lower)/2


def fit_pilot(ev: Evaluator, start: np.ndarray, *, sobol_points: int,
              screen_points: int, active_params: int, starts: int,
              local_max_nfev: int, perturbation: float,
              rng_seed: int) -> dict:
    """Actual Sobol start search + sensitivity-screened SciPy least_squares.

    No frozen Wilson coefficients: EVERY new function/Jacobian point calls
    Evaluator.evaluate, which recomputes the source-resolved EFT pipeline.
    The budget limits NEW physics evaluations, not cached requests. The
    returned best-fit is the best verified conditional trial, not a proof
    of global optimality or a complete physical MSbar prediction.
    """
    from scipy.optimize import least_squares
    from scipy.stats import qmc

    lower, upper, widths = _trial_box(start, perturbation)
    centre = (lower+upper)/2
    candidates: dict[str, dict] = {}
    history: list[dict] = []
    best: dict | None = None

    def register(result: dict, phase: str) -> None:
        nonlocal best
        if result.get("status") != "COMPLETE_CONDITIONAL":
            history.append({"phase": phase, "status": "BLOCKED",
                            "reason": result.get("error")})
            return
        candidates[result["candidate_digest"]] = result
        if best is None or float(result["chi2_approx"]) < float(best["chi2_approx"]):
            best = result
        history.append({"phase": phase, "status": "COMPLETE_CONDITIONAL",
                        "candidate_digest": result["candidate_digest"],
                        "chi2_approx": result["chi2_approx"]})

    seed = ev.evaluate(start)
    register(seed, "thesis_seed")
    if seed["status"] != "COMPLETE_CONDITIONAL":
        return {"status": "SEED_BLOCKED", "seed": seed, "best": None,
                "new_evaluations": ev.attempts, "fit_revision": FIT_REVISION}

    # A scrambled 18D Sobol design searches *all* real Yukawa entries.
    # Each evaluated vector is inside the same predeclared bounded box.
    n_draw = 2 ** math.ceil(math.log2(sobol_points))
    sobol = qmc.Sobol(d=18, scramble=True, seed=rng_seed)
    design = sobol.random_base2(m=int(math.log2(n_draw)))
    for unit in design[:sobol_points]:
        if ev.attempts >= ev.limit:
            break
        register(ev.evaluate(lower+(upper-lower)*unit), "sobol")
    valid_sobol = [v for v in candidates.values() if v["status"] == "COMPLETE_CONDITIONAL"]
    if not valid_sobol:
        fail("No valid candidate remained after Sobol sampling")
    valid_sobol.sort(key=lambda x: float(x["chi2_approx"]))

    # First estimate full 18-coordinate sensitivities from the Sobol design.
    # Then, if budget permits, replace those proxies by direct one-sided
    # finite differences at the best Sobol point. Only a few selected entries
    # are subsequently optimised; the others remain fixed for each start.
    anchor = np.asarray(best["candidate"], float)
    r0 = np.asarray(best["residuals"], float)
    scores = np.zeros(18, float)
    for entry in valid_sobol:
        x = np.asarray(entry["candidate"], float)
        scaled = (x-anchor)/widths
        dr = np.asarray(entry["residuals"], float)-r0
        scores += abs(scaled) * float(np.linalg.norm(dr))
    scores /= max(len(valid_sobol), 1)
    screened: list[dict] = []
    for index in range(min(18, screen_points)):
        if ev.attempts >= ev.limit:
            break
        delta = min(0.05*widths[index], upper[index]-anchor[index])
        if delta <= max(1e-18, 1e-10*widths[index]):
            delta = -min(0.05*widths[index], anchor[index]-lower[index])
        if abs(delta) <= max(1e-18, 1e-10*widths[index]):
            continue
        x = anchor.copy()
        x[index] += delta
        result = ev.evaluate(x)
        register(result, "sensitivity")
        if result["status"] == "COMPLETE_CONDITIONAL":
            dr = np.asarray(result["residuals"], float)-r0
            scores[index] = float(np.linalg.norm(dr)) * widths[index] / abs(delta)
            screened.append({"index": index, "score": float(scores[index]),
                             "candidate_digest": result["candidate_digest"]})

    active = np.argsort(-scores, kind="stable")[:active_params]
    if not np.isfinite(scores).all():
        fail("Nonfinite Yukawa sensitivity screening")
    print(ev.model, "Sobol complete; best chi2", best["chi2_approx"],
          "active Yukawa indices", active.tolist(), flush=True)
    least_squares_runs: list[dict] = []

    # Choose multiple *independent* Sobol starts, with the best current result
    # included first. A local optimum is conditional on the active subspace.
    selected = []
    for c in sorted(candidates.values(), key=lambda x: float(x["chi2_approx"])):
        if len(selected) >= starts:
            break
        if c["candidate_digest"] not in {v["candidate_digest"] for v in selected}:
            selected.append(c)

    for ordinal, chosen in enumerate(selected, start=1):
        if ev.attempts >= ev.limit:
            break
        frozen = np.asarray(chosen["candidate"], float)
        z0 = np.clip((frozen[active]-centre[active])/widths[active], -1., 1.)
        inner_cache: dict[tuple[float, ...], np.ndarray] = {}

        def objective(z: np.ndarray) -> np.ndarray:
            key = tuple(np.asarray(z, float))
            if key in inner_cache:
                return inner_cache[key]
            x = frozen.copy()
            x[active] = centre[active] + widths[active]*np.asarray(z, float)
            entry = ev.evaluate(x)
            register(entry, "least_squares")
            if entry["status"] == "COMPLETE_CONDITIONAL":
                value = np.asarray(entry["residuals"], float)
            else:
                # A mass-unsafe or ill-defined-mixing point has no usable fit.
                # Keep least_squares finite, but never report it as a result.
                value = np.full(len(ORDER), 1e6, dtype=float)
            inner_cache[key] = value
            return value

        def finite_difference_jacobian(z: np.ndarray) -> np.ndarray:
            z = np.asarray(z, float)
            value = objective(z)
            jac = np.empty((len(ORDER), len(active)), dtype=float)
            for j in range(len(active)):
                step = 0.005 if z[j] < 0.99 else -0.005
                probe = z.copy()
                probe[j] += step
                probe[j] = np.clip(probe[j], -1., 1.)
                actual = probe[j]-z[j]
                if abs(actual) < 1e-12:
                    jac[:,j] = 0.
                else:
                    jac[:,j] = (objective(probe)-value)/actual
            return jac

        info = {"start_number": ordinal, "start_digest": chosen["candidate_digest"],
                "start_chi2": chosen["chi2_approx"],
                "active_parameter_indices": active.tolist(),
                "local_solver": "scipy.optimize.least_squares_trf_bounded"}
        try:
            fitted = least_squares(objective, z0,
                                   jac=finite_difference_jacobian,
                                   method="trf", bounds=(-np.ones(len(active)), np.ones(len(active))),
                                   x_scale="jac", max_nfev=local_max_nfev,
                                   ftol=1e-6, xtol=1e-6, gtol=1e-6)
            info.update(status="FINISHED", solver_success=bool(fitted.success),
                        nfev=int(fitted.nfev), njev=int(fitted.njev),
                        message=str(fitted.message),
                        solver_final_cost=float(fitted.cost))
        except RuntimeError as exc:
            if str(exc) != "Pilot evaluation budget reached":
                raise
            info.update(status="BUDGET_EXHAUSTED", message=str(exc))
        least_squares_runs.append(info)
        print(ev.model, "least_squares start", ordinal, info["status"],
              "best chi2", best["chi2_approx"], flush=True)

    refined = any(v["status"] in ("FINISHED", "BUDGET_EXHAUSTED")
                  for v in least_squares_runs)
    result_status = "PILOT_BEST_FOUND" if refined else "SOBOL_ONLY_NO_LOCAL_FIT"
    return {"status": result_status, "seed": seed, "best": best,
            "new_evaluations": ev.attempts, "fit_revision": FIT_REVISION,
            "search_method": "scrambled Sobol + screened bounded scipy least_squares",
            "search_box_lower": lower.tolist(), "search_box_upper": upper.tolist(),
            "sobol_points_requested": sobol_points,
            "direct_sensitivity_evaluations": screened,
            "active_parameter_indices": active.tolist(),
            "least_squares_runs": least_squares_runs,
            "candidate_history": history,
            "global_optimality_certified": False,
            "fit_acceptable_project_cut": best["chi2_approx"] <= 5.86}


def save_summary(ev: Evaluator, out: dict, *, model: str) -> None:
    folder = ev.output_dir()
    folder.mkdir(parents=True, exist_ok=True)
    out = dict(out)
    out.update(model=model, schema="t3-thesis-two-model-conditional-pilot-result-v1",
               target="configs/nufit_6_1_ic24_no.json",
               thesis_benchmark_config=str(CONFIG),
               physical_MSbar_matching_certified=False)
    (folder/"pilot_summary.json").write_text(json.dumps(out, indent=2, sort_keys=True,
                                                         allow_nan=False)+"\n", encoding="utf-8")
    valid = out.get("best") or out.get("seed")
    if valid and valid.get("status") == "COMPLETE_CONDITIONAL":
        vec = np.asarray(valid["candidate"], float)
        matrices = {"Y1_UV": vec[:9].reshape(3,3).tolist(), "Y2_UV": vec[9:].reshape(3,3).tolist()}
        (folder/"best_yukawas.json").write_text(json.dumps({"model":model,**matrices,
            "chi2_approx":valid["chi2_approx"],"status":out["status"],
            "source_candidate_digest":valid["candidate_digest"]}, indent=2)+"\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plan", action="store_true", help="Show input values only; no numerical calculation")
    p.add_argument("--baseline", action="store_true", help="Evaluate thesis Appendix-B Yukawas without optimisation")
    p.add_argument("--fit", action="store_true", help="Sobol sample + sensitivity screening + bounded least_squares refit")
    p.add_argument("--model", choices=MODELS, help="Run only one of the two pilot models")
    p.add_argument("--config", type=Path, default=CONFIG)
    p.add_argument("--max-evals", type=int, default=80,
                   help="Maximum NEW complete EFT evaluations; cached successful points are free")
    p.add_argument("--sobol-points", type=int, default=16,
                   help="Number of scrambled Sobol samples across all 18 real Yukawa entries")
    p.add_argument("--screen-points", type=int, default=18,
                   help="Number of direct coordinate sensitivity probes (0..18)")
    p.add_argument("--active-params", type=int, default=5,
                   help="Highest-sensitivity Yukawa coordinates in the local least-squares fit")
    p.add_argument("--starts", type=int, default=1,
                   help="Number of best Sobol seeds for independent least-squares refinements")
    p.add_argument("--local-max-nfev", type=int, default=12,
                   help="Maximum SciPy least_squares objective calls per start (Jacobian trials also cost EFT evaluations)")
    p.add_argument("--perturbation", type=float, default=1.0,
                   help="Fractional half-width around thesis Yukawa seeds; each entry has >=2e-5 scale")
    p.add_argument("--seed", type=int, default=17)
    p.add_argument("--keep-work", action="store_true", help="Keep all independent full-chain source artifacts")
    args = p.parse_args(argv)
    if sum((args.plan, args.baseline, args.fit)) != 1:
        p.error("Select exactly one of --plan, --baseline, or --fit")
    if (args.max_evals < 1 or args.sobol_points < 0 or
            not 0 <= args.screen_points <= 18 or
            not 1 <= args.active_params <= 18 or args.starts < 1 or
            args.local_max_nfev < 1 or not (0 < args.perturbation <= 20)):
        p.error("Invalid Sobol / least-squares search limits")
    if args.fit:
        # Initial seed, Sobol, sensitivity, and one complete local Jacobian.
        minimum = 1 + args.sobol_points + args.screen_points + args.active_params + 2
        if args.max_evals < minimum:
            p.error(f"Sobol plus least_squares needs --max-evals >= {minimum} "
                    f"for this configuration (got {args.max_evals}); "
                    "e.g. use --max-evals 80")
        if args.sobol_points < 1:
            p.error("--fit requires at least one Sobol point")
    cfg = load_config(ROOT, args.config)
    models = (args.model,) if args.model else MODELS
    target = read_json(ROOT / cfg["fitting"]["target"])
    if target.get("ordering") != "NO" or not all(k in target.get("central_values", {}) and
                                                   k in target.get("one_sigma_errors", {}) for k in ORDER):
        fail("Missing normal-ordering five-observable target")
    if args.plan:
        print("THESIS PILOT PLAN (no numerical execution)")
        print("Models:", ", ".join(models), "scales:", cfg["scales_GeV"])
        print("lambdaT3:", cfg["lambdaT3_UV"], "lambda_other:", cfg["lambda_other_UV"])
        print("Old fitted real Yukawa seeds: thesis Appendix B; no class-E factor")
        print("Matching prescription: analytic_triangle strict one-loop (native alternatives not fitted)")
        print("SM yu/yd/ye: Workstream B defaults (not specified in thesis Table 4.1)")
        print("Results root:", ROOT / OUTPUT)
        return 0
    folder = safe_parent(ROOT, OUTPUT)
    folder.mkdir(parents=True, exist_ok=True)
    all_results = []
    for model in models:
        start = initial_parameters(cfg, model)
        ev = Evaluator(ROOT, cfg, model, target, args.keep_work, args.max_evals, [])
        if args.baseline:
            one = ev.evaluate(start)
            outcome = {"status":"BASELINE_COMPLETE" if one["status"] == "COMPLETE_CONDITIONAL" else "BASELINE_BLOCKED",
                       "seed":one, "best":one if one["status"]=="COMPLETE_CONDITIONAL" else None,
                       "new_evaluations":ev.attempts}
        else:
            outcome = fit_pilot(ev, start,
                                sobol_points=args.sobol_points,
                                screen_points=args.screen_points,
                                active_params=args.active_params,
                                starts=args.starts,
                                local_max_nfev=args.local_max_nfev,
                                perturbation=args.perturbation,
                                rng_seed=args.seed)
        save_summary(ev, outcome, model=model)
        all_results.append({"model":model,"status":outcome["status"],
                            "chi2_approx":(outcome.get("best") or {}).get("chi2_approx"),
                            "new_evaluations":ev.attempts})
        if outcome["status"] not in ("BASELINE_COMPLETE", "PILOT_BEST_FOUND"):
            seed_data = outcome.get("seed") or {}
            print(model, "PILOT BLOCKED:", seed_data.get("error", outcome["status"]),
                  file=sys.stderr, flush=True)
    # Preserve the status of whichever model was run in an earlier invocation.
    # Running A_-4 after D_-2 must not erase D_-2 from the overview.
    for previous in MODELS:
        if previous in models:
            continue
        summary_file = folder / previous.replace("-", "m").replace("+", "p") / "pilot_summary.json"
        if summary_file.is_file():
            historic = read_json(summary_file)
            all_results.append({"model": previous,
                                "status": historic.get("status", "UNKNOWN"),
                                "chi2_approx": (historic.get("best") or {}).get("chi2_approx"),
                                "new_evaluations": historic.get("new_evaluations", 0)})
    all_results.sort(key=lambda x: MODELS.index(x["model"]))
    target_file = folder / "two_model_overview.csv"
    with target_file.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["model","status","chi2_approx","new_evaluations"])
        writer.writeheader(); writer.writerows(all_results)
    print("Saved pilot summary:",target_file)
    # Previous-model records are retained in the overview for convenience,
    # but MUST NOT determine the exit status of the selected model(s).
    selected = [x for x in all_results if x["model"] in models]
    return 0 if all(x["status"] in ("BASELINE_COMPLETE", "PILOT_BEST_FOUND")
                    for x in selected) else 2


if __name__ == "__main__":
    sys.exit(main())
