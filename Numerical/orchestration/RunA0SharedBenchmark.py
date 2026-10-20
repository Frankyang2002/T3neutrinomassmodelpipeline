"""Run one fixed Class-A(alpha=0) Yukawa pair through every T3 model.

Local study only. Never modifies individually fitted references or GitHub.
Requires the pre-existing RGBeta and final-C5 model artifacts.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.stats import chi2 as chi2_dist

from Numerical.orchestration import RunSmallLambdaRefit as base
from Numerical.orchestration.PipelineNumericalResults import _run_and_build_payload
from Numerical.plotting.RunningResultFigures import generate_running_result_figures

ROOT = base.STUDY / "a0_shared_benchmark"
DONOR = "T3_dS1_1_dS2_3_dF_2_alpha_p0"
FIELDS = ("delta_m21_sq_ev2", "delta_m3l_sq_ev2", "sin2_theta12", "sin2_theta13", "sin2_theta23")


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def load_inputs(artifact_map):
    references = base.original_rows()
    source_path = base.STUDY / "fitted" / (DONOR + ".json")
    if not source_path.is_file():
        raise FileNotFoundError(f"Missing A(0) fitted Yukawas: {source_path}")
    source = base.load(source_path)
    donor_ordinary = source["base_state"]["ordinary"]
    shared = {y: base.matr(donor_ordinary[y]) for y in ("y1", "y2")}
    mapping = base.load(artifact_map) if artifact_map else None
    models = []
    for row in references:
        name = row["model"]
        cfg_path = base.STUDY / "fitted" / (name + ".json")
        cfg = base.load(cfg_path)
        ordinary = cfg["base_state"]["ordinary"]
        expected = 5e-8 if base.cls(row) == "E" else 1e-7
        if not np.isclose(float(ordinary["lambdaT3"]["real"]), expected, rtol=1e-10, atol=0):
            raise ValueError(f"{name}: unexpected lambdaT3")
        for key in ("y1", "y2"):
            if np.any(np.abs(base.matr(ordinary[key]).imag) > 1e-12):
                raise ValueError(f"{name}: original imaginary Yukawa not zero; this study expects real 18-parameter fits")
        common = copy.deepcopy(cfg)
        for key in ("y1", "y2"):
            base.setmat(common["base_state"]["ordinary"][key], shared[key])
        common["benchmark_search"] = {"enabled": False}
        common["sensitivity"] = {"enabled": False}
        common["a0_shared_study"] = {"donor": DONOR, "role": "unchanged A(0) shared Yukawa pair", "lambdaT3": expected}
        models.append((row, common, base.locate_artifacts(name, mapping)))
    return shared, models


def make_record(row, output_folder, artifacts):
    r = row["representation"]
    return SimpleNamespace(
        name=row["model"], output_dir=output_folder,
        d_s1=int(r["d_s1"]), d_s2=int(r["d_s2"]), d_f=int(r["d_f"]),
        alpha=int(r["alpha"]), shared_scalar=bool(r.get("shared_scalar", False)),
        summary={"UVRGEFile": str(artifacts[0].resolve()),
                 "EFT1RenormalisableRGEFile": str(artifacts[1].resolve()),
                 "FinalWeinbergCoefficientFile": str(artifacts[2].resolve())},
    )


def chi2_value(row, common, artifacts, shared):
    parameters, evaluator = base.get_evaluator(common, artifacts)
    v = np.r_[shared["y1"].real.ravel(), shared["y2"].real.ravel()]
    fit = base.evaluate_point(evaluator, parameters, v)
    return float(fit.chi2), {k: float(val) for k, val in zip(fit.observable_names, fit.prediction)}


def plot_summary(rows, cutoff):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib unavailable: no summary plots")
        return
    ok = [r for r in rows if r["status"] == "Success"]
    if not ok:
        return
    ok.sort(key=lambda x: x["model"])
    fig, ax = plt.subplots(figsize=(13, 5))
    ax.bar(range(len(ok)), [max(r["chi2"], 1e-17) for r in ok])
    ax.axhline(cutoff, color="firebrick", ls="--", label=f"Nominal 95% cutoff {cutoff:.3f}")
    ax.axhline(5.86, color="darkorange", ls=":", label="Stricter screening cutoff 5.86")
    ax.set_yscale("log")
    ax.set_ylabel(r"Shared A(0) Yukawa $\chi^2$")
    ax.set_xticks(range(len(ok)), [r["model"] for r in ok], rotation=90, fontsize=6)
    ax.legend(); fig.tight_layout()
    fig.savefig(ROOT / "all_models_chi2.png", dpi=180)
    plt.close(fig)
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    target = OscillationFitTarget.from_json(base.TARGET)
    names = list(target.observable_names)
    for field in FIELDS:
        if field not in names:
            continue
        idx = names.index(field)
        values = [r["predictions"][field] for r in ok]
        center = float(target.central_values[idx])
        sigma = float(np.sqrt(target.covariance[idx, idx]))
        fig, ax = plt.subplots(figsize=(13, 4.2))
        ax.scatter(range(len(ok)), values, s=23)
        ax.axhline(center, color="black", lw=1, label="NuFIT central")
        ax.axhspan(center - sigma, center + sigma, alpha=.18, label="Marginal 1-sigma band")
        ax.set_xticks(range(len(ok)), [r["model"] for r in ok], rotation=90, fontsize=6)
        ax.set_ylabel(field); ax.legend(); fig.tight_layout()
        fig.savefig(ROOT / (field + "_comparison.png"), dpi=170)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--artifact-map", type=Path)
    ap.add_argument("--selected-only", action="store_true", help="Run five strict-pass models and the sixth nominal 95% model; default all 16")
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--resume", action="store_true", help="Reuse successful saved diagnostics and result records")
    args = ap.parse_args()
    cutoff = float(chi2_dist.ppf(.95, 5))
    shared, models = load_inputs(args.artifact_map)
    atomic_json(ROOT / "shared_yukawas.json", {
        "source_model": DONOR,
        "y1": {"real": shared["y1"].real.tolist(), "imag": shared["y1"].imag.tolist()},
        "y2": {"real": shared["y2"].real.tolist(), "imag": shared["y2"].imag.tolist()},
        "max_abs_y1": float(np.max(np.abs(shared["y1"]))),
        "max_abs_y2": float(np.max(np.abs(shared["y2"]))),
        "rms_abs_y1": float(np.sqrt(np.mean(np.abs(shared["y1"])**2))),
        "rms_abs_y2": float(np.sqrt(np.mean(np.abs(shared["y2"])**2))),
    })
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    target = OscillationFitTarget.from_json(base.TARGET)
    if len(target.observable_names) != 5:
        raise ValueError("Expected five fit observables")
    approved = {"T3_dS1_1_dS2_3_dF_2_alpha_m2", DONOR,
                "T3_dS1_2_dS2_2_dF_1_alpha_m3", "T3_dS1_2_dS2_2_dF_1_alpha_p1",
                "T3_dS1_3_dS2_1_dF_2_alpha_m2", "T3_dS1_3_dS2_1_dF_2_alpha_p0"}
    if args.selected_only:
        models = [entry for entry in models if entry[0]["model"] in approved]
    rows = []
    ROOT.mkdir(parents=True, exist_ok=True)
    for i, (row, config, artifacts) in enumerate(models, 1):
        name = row["model"]
        folder = ROOT / "models" / name
        diagnostics = folder / "running_diagnostics.json"
        result_file = folder / "evaluation.json"
        config_file = folder / "shared_config.json"
        if args.resume and result_file.is_file() and diagnostics.is_file():
            rec = base.load(result_file)
            if rec.get("status") == "Success":
                rows.append(rec)
                print(f"[{i}/{len(models)}] REUSED {name}: chi2={rec['chi2']:.6g}", flush=True)
                continue
        folder.mkdir(parents=True, exist_ok=True)
        atomic_json(config_file, config)
        record = {"model": name, "class": base.cls(row), "lambdaT3": float(config["base_state"]["ordinary"]["lambdaT3"]["real"]), "donor": DONOR}
        try:
            score, prediction = chi2_value(row, config, artifacts, shared)
            synthetic = make_record(row, folder, artifacts)
            payload = _run_and_build_payload(record=synthetic, payload=config, config_path=config_file)
            # Independent cross-check against the full diagnostic low-energy predictions.
            last = payload["final_running"]
            diagnostics_prediction = {key: float(last[key][-1]) for key in FIELDS}
            residual = np.array([diagnostics_prediction[k] - target.central_values[list(target.observable_names).index(k)] for k in target.observable_names])
            whitened = np.linalg.solve(np.linalg.cholesky(target.covariance), residual)
            diagnostic_chi2 = float(whitened @ whitened)
            if not np.isclose(diagnostic_chi2, score, rtol=1e-4, atol=1e-5):
                raise RuntimeError(f"Cross-evaluator chi2={score:.7g} differs from trajectory chi2={diagnostic_chi2:.7g}")
            atomic_json(diagnostics, payload)
            if not args.no_figures:
                generate_running_result_figures(diagnostics, folder / "figures")
            record.update({"status": "Success", "chi2": score, "pass_strict_5p86": score < 5.86,
                           "pass_nominal_95": score <= cutoff, "predictions": prediction,
                           "diagnostic_chi2": diagnostic_chi2})
            print(f"[{i}/{len(models)}] {name}: chi2={score:.6g}, strict={score < 5.86}", flush=True)
        except Exception as exc:
            record.update({"status": "Failed", "error": f"{type(exc).__name__}: {exc}"})
            print(f"[{i}/{len(models)}] FAILED {name}: {record['error']}", flush=True)
        atomic_json(result_file, record)
        rows.append(record)
        atomic_json(ROOT / "summary.json", {"source": DONOR, "chi2_nominal_95_cutoff": cutoff, "models": rows})
    with (ROOT / "comparison.csv").open("w", newline="", encoding="utf-8") as f:
        fieldnames = ["model", "class", "lambdaT3", "status", "chi2", "pass_strict_5p86", "pass_nominal_95", *FIELDS, "error"]
        writer = csv.DictWriter(f, fieldnames=fieldnames); writer.writeheader()
        for rec in rows:
            writer.writerow({**{key: rec.get(key, "") for key in fieldnames},
                             **{key: rec.get("predictions", {}).get(key, "") for key in FIELDS}})
    if not args.no_figures:
        plot_summary(rows, cutoff)
    print(f"Complete: {sum(r['status']=='Success' for r in rows)}/{len(rows)} evaluated. "
          f"Strict={sum(r.get('pass_strict_5p86',False) for r in rows)}, "
          f"nominal95={sum(r.get('pass_nominal_95',False) for r in rows)}. Output: {ROOT}", flush=True)


if __name__ == "__main__":
    main()
