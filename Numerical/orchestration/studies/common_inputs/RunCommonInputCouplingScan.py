"""Numerical-only coupling-strength scan around saved common-input T3 points.

All matching and beta-function expressions are loaded from PRECOMPUTED JSON
artifacts. No symbolic calculations, fitting, or benchmark search is performed.

Separate one-parameter scans:
  yukawa: y1, y2 -> k*y1, k*y2; lambdaT3 remains fixed
  lambda:  lambdaT3 -> k*lambdaT3; y1, y2 remain fixed

The existing class-E 1/2 lambdaT3 convention is preserved by multiplying
its saved common-input value (0.05), not replacing it with 0.1.

Each point compares the full RGE result with the existing frozen-parameter
control. The control includes one-loop C5 matching but is not a complete
fixed-order EFT calculation.

From the repository root:
    python -m Numerical.orchestration.studies.common_inputs.RunCommonInputCouplingScan
    python -m Numerical.orchestration.studies.common_inputs.RunCommonInputCouplingScan --run --limit 2
    python -m Numerical.orchestration.studies.common_inputs.RunCommonInputCouplingScan --run
    python -m Numerical.orchestration.studies.common_inputs.RunCommonInputCouplingScan --run --all-models
"""
from __future__ import annotations

import argparse
import csv
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from Numerical.orchestration.support.ResumePipelineNumericalResults import (
    PROJECT_ROOT, OUTPUT_ROOT, REQUIRED, discover_cases, parse_model_key,
)
from Numerical.orchestration.studies.common_inputs.RunCommonInputFermionScale import SELECTED
from Numerical.orchestration.studies.scalar_rge.RunCommonScalarComparison import evaluate
from Numerical.orchestration.studies.scalar_rge.RunCommonScalarRGControls import fixed_point
from Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison import OBS
from Numerical.fitting.OscillationFit import OscillationFitTarget

SCENARIO = "smallY_smallL"
DEFAULT_FACTORS = (0.5, 0.75, 1.0, 1.5, 2.0)
MODES = ("yukawa", "lambda")


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def make_signature(config_path: Path, artifacts: tuple[Path, ...], factors: list[float], modes: list[str]) -> str:
    h = hashlib.sha256(config_path.read_bytes())
    h.update(json.dumps({"factors": factors, "modes": modes}, sort_keys=True).encode())
    for path in artifacts:
        stat = path.stat()
        h.update(str(path.resolve()).encode())
        h.update(f"{stat.st_size}:{stat.st_mtime_ns}".encode())
    return h.hexdigest()


def deform(payload: dict, mode: str, factor: float) -> dict:
    """Scale exactly the requested numerical couplings without altering representation."""
    if mode not in MODES or not math.isfinite(factor) or factor <= 0:
        raise ValueError("Invalid coupling mode or factor")
    output = deepcopy(payload)
    ordinary = output["base_state"]["ordinary"]
    if mode == "yukawa":
        for name in ("y1", "y2"):
            spec = ordinary[name]
            real = np.asarray(spec["real"], dtype=float)
            imag = np.asarray(spec.get("imag", np.zeros_like(real)), dtype=float)
            if real.shape != (3, 3) or imag.shape != real.shape:
                raise ValueError(f"Unexpected {name} shape")
            spec["real"] = (factor * real).tolist()
            spec["imag"] = (factor * imag).tolist()
    else:
        spec = ordinary["lambdaT3"]
        if isinstance(spec, dict):
            spec["real"] = factor * float(spec["real"])
            spec["imag"] = factor * float(spec.get("imag", 0.0))
        else:
            ordinary["lambdaT3"] = factor * float(spec)
    output["benchmark_search"] = {**output.get("benchmark_search", {}), "enabled": False}
    output["sensitivity"] = {**output.get("sensitivity", {}), "enabled": False}
    output["scan"] = {"bindings": {}}
    return output


def short_label(key: str) -> str:
    dims = parse_model_key(key)
    mapping = {(1, 3, 2): "A", (2, 2, 1): "B", (2, 2, 3): "C", (3, 1, 2): "D", (3, 3, 2): "E"}
    a = dims[3]
    return f"{mapping.get(dims[:3], '?')}-{'m' if a < 0 else 'p'}{abs(a)}"


def save_plots(rows: list[dict], dest: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib unavailable; numerical JSON and CSV are complete.")
        return
    dest.mkdir(parents=True, exist_ok=True)
    labels = {
        "delta_m21_sq_ev2": r"$\Delta m^2_{21}$ (eV$^2$)",
        "delta_m3l_sq_ev2": r"$\Delta m^2_{3\ell}$ (eV$^2$)",
        "sin2_theta12": r"$\sin^2\theta_{12}$",
        "sin2_theta13": r"$\sin^2\theta_{13}$",
        "sin2_theta23": r"$\sin^2\theta_{23}$",
    }
    models = list(dict.fromkeys(r["model"] for r in rows))
    palette = plt.get_cmap("tab20")
    for mode in MODES:
        if not any(r["mode"] == mode and r["status"] == "Success" for r in rows):
            continue
        for observable in OBS:
            fig, axes = plt.subplots(1, 2, figsize=(13, 5))
            for index, model in enumerate(models):
                points = sorted((r for r in rows if r["model"] == model and r["mode"] == mode),
                                key=lambda r: r["factor"])
                chunks, current = [], []
                for point in points:
                    if point["status"] == "Success":
                        current.append(point)
                    elif current:
                        chunks.append(current)
                        current = []
                if current:
                    chunks.append(current)
                for j, chunk in enumerate(chunks):
                    factors = [r["factor"] for r in chunk]
                    color = palette(index % 20)
                    axes[0].plot(factors, [r["running"][observable] for r in chunk], "-o",
                                 lw=1.3, ms=3, color=color, label=short_label(model) if j == 0 else None)
                    axes[0].plot(factors, [r["fixed"][observable] for r in chunk], "--",
                                 lw=1.1, color=color)
                    axes[1].plot(factors, [r["relative_percent"][observable] for r in chunk], "-o",
                                 lw=1.3, ms=3, color=color)
            axes[0].set(title="Solid: full RGE; dashed: frozen parameters", ylabel=labels[observable])
            axes[1].set(title="Additional RGE effect", ylabel=r"$100(O_{\rm RGE}-O_{\rm frozen})/O_{\rm frozen}$ (%)")
            for ax in axes:
                ax.set_xscale("log")
                ax.axvline(1, lw=0.7, ls=":", color="grey")
                ax.set_xlabel("Common coupling-strength multiplier k")
                ax.grid(alpha=0.25)
            axes[0].legend(fontsize=7, ncol=2)
            fig.suptitle("Common-input T3: " + ("both Yukawas scaled together" if mode == "yukawa" else "lambdaT3 scaled alone"))
            fig.tight_layout()
            fig.savefig(dest / f"{mode}_{observable}_running_vs_frozen.png", dpi=170)
            plt.close(fig)
    print("Wrote comparison figures:", dest)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--scenario", default=SCENARIO)
    p.add_argument("--run", action="store_true", help="Run numerical scans (default: inventory)")
    p.add_argument("--limit", type=int, default=0, help="Limit selected models; 0 means all selected")
    p.add_argument("--all-models", action="store_true", help="Include all available common-input models")
    p.add_argument("--factors", type=float, nargs="+", default=list(DEFAULT_FACTORS))
    p.add_argument("--modes", nargs="+", choices=MODES, default=list(MODES))
    p.add_argument("--force", action="store_true", help="Recalculate completed scans")
    args = p.parse_args(argv)
    if args.limit < 0 or not args.factors or any(not math.isfinite(v) or v <= 0 for v in args.factors):
        p.error("Require nonnegative limit and finite positive factors")
    factors = sorted(set(args.factors) | {1.0})
    modes = list(dict.fromkeys(args.modes))
    common_root = PROJECT_ROOT / "Reports" / "output" / "full" / "common_input_rge" / args.scenario
    destination = PROJECT_ROOT / "Reports" / "output" / "full" / "common_input_coupling_scan" / args.scenario
    target_path = PROJECT_ROOT / "configs" / "nufit_6_1_ic24_no.json"
    if not common_root.is_dir():
        raise FileNotFoundError(f"Missing common-input configurations: {common_root}")
    cases = {}
    for scenario, key, case in discover_cases(OUTPUT_ROOT):
        if scenario == args.scenario:
            if key in cases:
                raise RuntimeError("Duplicate comparison artifacts: " + key)
            cases[key] = case
    all_available = sorted(p.name for p in common_root.iterdir()
                           if p.is_dir() and (p / "common_input_config.json").is_file())
    models = all_available if args.all_models else [name for name in SELECTED if name in all_available]
    models = models[:args.limit or None]
    if not models:
        print("No selected common-input models found")
        return 1
    target = OscillationFitTarget.from_json(target_path) if args.run else None
    print(f"Models: {len(models)}; modes: {modes}; factors: {factors}")
    rows_all = []
    has_errors = False
    for model in models:
        path = common_root / model / "common_input_config.json"
        if model not in cases:
            print("Missing numerical artifacts for", model)
            has_errors = True
            continue
        artifacts = tuple(cases[model] / "data" / name for name in REQUIRED)
        if any(not p.is_file() for p in artifacts):
            print("Missing precomputed RGE/matching JSON for", model)
            has_errors = True
            continue
        base = read_json(path)
        rep = base["representation"]
        if tuple(int(rep[k]) for k in ("d_s1", "d_s2", "d_f", "alpha")) != parse_model_key(model):
            raise ValueError("Representation mismatch for " + model)
        if base.get("benchmark_search", {}).get("enabled", False):
            raise ValueError("Refusing fitted/search-enabled config for " + model)
        print("Ready:", short_label(model))
        if not args.run:
            continue
        sig = make_signature(path, artifacts, factors, modes)
        saved = destination / model / "coupling_scan.json"
        if saved.is_file() and not args.force:
            cached = read_json(saved)
            if cached.get("signature") != sig:
                raise RuntimeError(f"Existing scan settings changed for {model}; rerun with --force")
            rows = cached["points"]
            print("Reused cached:", short_label(model))
        else:
            rows = []
            for mode in modes:
                for factor in factors:
                    row = {"model": model, "mode": mode, "factor": factor, "status": "Failed", "error": ""}
                    try:
                        modified = deform(base, mode, factor)
                        running = evaluate(modified, list(artifacts), target)["prediction"]
                        frozen = fixed_point(modified, artifacts[2], target)
                        row["running"] = {o: float(running[o]) for o in OBS}
                        row["fixed"] = {o: float(frozen[o]) for o in OBS}
                        row["delta"] = {o: row["running"][o] - row["fixed"][o] for o in OBS}
                        row["relative_percent"] = {
                            o: (100 * row["delta"][o] / row["fixed"][o]
                                if row["fixed"][o] != 0 else None) for o in OBS
                        }
                        if any(not math.isfinite(v) for name in ("running", "fixed", "delta")
                               for v in row[name].values()):
                            raise ValueError("Non-finite observable")
                        row["status"] = "Success"
                    except Exception as exc:
                        row["error"] = f"{type(exc).__name__}: {exc}"
                    rows.append(row)
                    print(f"  {short_label(model):7} {mode:6} k={factor:g}: {row['status']} {row['error']}", flush=True)
            saved.parent.mkdir(parents=True, exist_ok=True)
            saved.write_text(json.dumps({"signature": sig, "config": str(path),
                                         "artifacts": [str(p) for p in artifacts],
                                         "factors": factors, "modes": modes, "points": rows},
                                        indent=2, allow_nan=False) + "\n", encoding="utf-8")
        rows_all.extend(rows)
        has_errors |= any(r["status"] != "Success" for r in rows)
    if not args.run:
        print("Inventory only; pass --run to evaluate numerical points.")
        return 1 if has_errors else 0
    if not rows_all:
        print("No numerical results produced")
        return 1
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "scan.json").write_text(json.dumps(rows_all, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    fields = ["model", "mode", "factor", "status", "error"]
    for ob in OBS:
        fields.extend(pre + "_" + ob for pre in ("running", "fixed", "delta", "relative_percent"))
    with (destination / "scan.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows_all:
            flat = {k: row.get(k, "") for k in fields[:5]}
            for ob in OBS:
                for pre in ("running", "fixed", "delta", "relative_percent"):
                    flat[pre + "_" + ob] = (row.get(pre) or {}).get(ob, "")
            writer.writerow(flat)
    save_plots(rows_all, destination / "plots")
    good = sum(r["status"] == "Success" for r in rows_all)
    print(f"Completed {good}/{len(rows_all)} points; output: {destination}")
    return 1 if has_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
