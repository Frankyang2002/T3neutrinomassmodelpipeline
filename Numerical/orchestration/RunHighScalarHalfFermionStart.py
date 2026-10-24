"""lambdaT3 = 0.5: fermion-threshold-start 16-model Sobol/fit/benchmark/reports.

Thin local-only orchestrator of the existing 011 pipeline; no UV beta evaluation.
See docs/HIGH_SCALAR_HALF_FERMION_START.md for PowerShell commands.

The 0p5 case uses lambdaT3=0.5 for T3 classes A--D and 0.25 for E,
following the existing fermion-start study's class-E factor-of-two convention.
"""
from __future__ import annotations

import argparse
import html
import json
import math
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

CASE = "0p5"
COUPLINGS = (0.5, 0.25)
MODES = ("optimised", "shared_a_m2")
OBSERVABLES = (
    "delta_m21_sq_ev2", "delta_m3l_sq_ev2",
    "sin2_theta12", "sin2_theta13", "sin2_theta23",
)


def _modules():
    """Register one extra case in memory, without editing upstream modules."""
    from Numerical.orchestration import RunFermionStartLambdaStudy as study
    existing = study.CASES.get(CASE)
    if existing is not None and tuple(existing) != COUPLINGS:
        raise ValueError(f"Existing {CASE} configuration differs: {existing}")
    study.CASES[CASE] = COUPLINGS
    return study


def _replay_main(main, argv):
    """Call the unchanged CLI inside this process so its case registry is shared."""
    before = sys.argv
    try:
        sys.argv = [main.__module__, *argv]
        main()
    finally:
        sys.argv = before


def _inputs(study, args):
    from Numerical.orchestration import RunSmallLambdaRefit as base
    from Numerical.orchestration import RunLambda5e9Refit as old
    missing = []
    for row in base.original_rows():
        model = row["model"]
        old_fitted = old.OUT / "fitted" / (model + ".json")
        if not old_fitted.is_file():
            missing.append(str(old_fitted))
        if args.check_artifacts:
            try:
                mapping = base.load(args.artifact_map) if args.artifact_map else None
                study.artifacts_for(model, mapping)
            except (FileNotFoundError, KeyError, ValueError) as exc:
                missing.append(f"{model}: {exc}")
    if missing:
        raise FileNotFoundError("Missing required local source/artifacts:\n  " + "\n  ".join(missing))
    print("Preflight: 16 fitted seed configs available. " +
          ("All 16 numerical artifact sets located." if args.check_artifacts else
           "Artifact validation skipped (use --check-artifacts)."), flush=True)


def _prepare(study):
    """Avoid regenerating input JSON after Sobol checkpoints have been created."""
    root = study.ROOT / CASE
    manifest = root / "prepared_manifest.json"
    if manifest.exists():
        if not all(study.path_for(CASE, "prepared", r["model"]).is_file()
                   for r in _original_rows()):
            raise RuntimeError("Incomplete prepared set exists. Inspect manually before overwriting.")
        for row in _original_rows():
            path = study.path_for(CASE, "prepared", row["model"])
            cfg = json.loads(path.read_text(encoding="utf-8"))
            target = study.cfg_lambda(CASE, row)
            found = float(cfg["base_state"]["ordinary"]["lambdaT3"]["real"])
            if not math.isclose(found, target, rel_tol=1e-12, abs_tol=0):
                raise RuntimeError(f"Existing prepared input has wrong scalar coupling: {path}")
        print(f"Reusing prepared inputs without overwriting: {manifest}", flush=True)
    else:
        if (root / "fits").exists() or (root / "sobol").exists():
            raise RuntimeError("A partial fit/checkpoint exists but no prepared manifest; inspect first")
        study.prepare(CASE)
    _save_provenance(study)


def _original_rows():
    from Numerical.orchestration import RunSmallLambdaRefit as base
    return base.original_rows()


def _save_provenance(study):
    folder = study.ROOT / CASE
    folder.mkdir(parents=True, exist_ok=True)
    from Numerical.orchestration import RunSmallLambdaRefit as base
    project = base.ROOT
    try:
        revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=project,
                                  capture_output=True, text=True, timeout=5,
                                  check=True).stdout.strip()
    except (FileNotFoundError, OSError, subprocess.CalledProcessError,
            subprocess.TimeoutExpired):
        revision = None
    payload = {
        "case": CASE,
        "lambdaT3_A_to_D": COUPLINGS[0],
        "lambdaT3_E": COUPLINGS[1],
        "source_fitted_configs": "output/studies/lambda_5e-9_refit/fitted/",
        "boundary": "mu_F",
        "uv_running": False,
        "fermion_matching": True,
        "eft1_running": True,
        "scalar_matching": True,
        "weinberg_running": True,
        "donor": study.DONOR,
        "fit": "Existing Sobol+ranked 5/8/12-parameter least-squares optimiser",
        "repository_revision_local": revision,
        "warning": "Fitted oscillation observables are calibration targets, not predictions. "
                   "The standard fixed-order EFT path is not the exploratory coupled C12-C5 solver.",
    }
    path = folder / "workflow_provenance.json"
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _fit(study, args):
    study.fit(CASE, args)


def _compare(study, args, mode):
    upstream_mode = "shared" if mode == "shared_a_m2" else "optimised"
    study.compare(CASE, SimpleNamespace(mode=upstream_mode,
                                         artifact_map=args.artifact_map))


def _ensure_compared(study, modes):
    missing = [str(study.ROOT / CASE / "comparison" / m / "summary.json")
               for m in modes if not (study.ROOT / CASE / "comparison" / m / "summary.json").is_file()]
    if missing:
        raise FileNotFoundError("Run compare first. Missing:\n" + "\n".join(missing))


def _running(study, args, modes):
    """Use the existing true-EFT1 and Weinberg running renderer unchanged."""
    from Numerical.orchestration import BuildFermionStartRunningReports as running
    running.CASES[CASE] = COUPLINGS
    mode = modes[0] if len(modes) == 1 else "all"
    cmd = ["--case", CASE, "--mode", mode, "--samples", str(args.samples)]
    if args.model:
        cmd += ["--model", args.model]
    if args.artifact_map:
        cmd += ["--artifact-map", str(args.artifact_map)]
    _replay_main(running.main, cmd)


def _reports(study, args, modes):
    """Invoke existing InteractiveModelComparison and fitted-envelope builders."""
    from Numerical.orchestration import BuildFermionStartExistingReports as reports
    reports.CASES = tuple(dict.fromkeys((*reports.CASES, CASE)))
    mode = modes[0] if len(modes) == 1 else "all"
    cmd = ["--case", CASE, "--mode", mode, "--kind", args.kind,
           "--gamma-min", str(args.gamma_min), "--gamma-max", str(args.gamma_max),
           "--points", str(args.points), "--rge-points", str(args.rge_points)]
    if args.model:
        cmd += ["--model", args.model]
    if args.artifact_map:
        cmd += ["--artifact-map", str(args.artifact_map)]
    if args.force_reports:
        cmd += ["--force"]
    _replay_main(reports.main, cmd)


def _read_comparisons(study):
    summaries = {}
    for mode in MODES:
        path = study.ROOT / CASE / "comparison" / mode / "summary.json"
        if path.exists():
            source = json.loads(path.read_text(encoding="utf-8"))
            summaries[mode] = {r["model"]: r for r in source["models"]}
    return summaries


def _plots(study):
    """Build static endpoint figures, alongside existing interactive HTML builders."""
    from Numerical.orchestration.ModelNames import shorten
    from Numerical.orchestration import RunSmallLambdaRefit as base
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        raise ImportError("Static endpoint PNGs require matplotlib and numpy") from exc

    datasets = _read_comparisons(study)
    if len(datasets) != 2:
        raise FileNotFoundError("Both comparisons are required for endpoint figures")
    dest = base.ROOT / "Reports" / "output" / "display_fermion_start" / CASE
    figdir = dest / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    names = [r["model"] for r in base.original_rows()]
    labels = [shorten(m) for m in names]
    x = np.arange(len(names))
    result_files = []

    def values(mode, key):
        out = []
        for model in names:
            r = datasets[mode].get(model, {})
            v = r.get("prediction", {}).get(key) if key != "chi2" else r.get("chi2")
            out.append(float(v) if r.get("status") == "Success" and v is not None else np.nan)
        return np.asarray(out)

    groups = [
        ("endpoint_mixing_angles.png", OBSERVABLES[2:], "Mixing-angle observables"),
        ("endpoint_mass_splittings.png", OBSERVABLES[:2], "Mass-squared differences"),
        ("endpoint_chi2.png", ("chi2",), "Endpoint chi-squared (log10(1+chi2))"),
    ]
    for filename, fields, title in groups:
        fig, axes = plt.subplots(len(fields), 1, figsize=(15, 3.25 * len(fields) + 1),
                                 sharex=True, squeeze=False)
        for j, key in enumerate(fields):
            ax = axes[j, 0]
            for mode in MODES:
                raw = values(mode, key)
                y = np.log10(1.0 + raw) if key == "chi2" else raw
                ax.plot(x, y, "o-", markersize=3, label=mode)
            ax.set_ylabel("log10(1 + χ²)" if key == "chi2" else key)
            ax.grid(True, alpha=0.25)
            if j == 0:
                ax.legend()
        axes[-1, 0].set_xticks(x, labels, rotation=70, ha="right")
        fig.suptitle(f"{title} | λT3=0.5 (class E=0.25) | UV frozen")
        fig.tight_layout()
        fig.savefig(figdir / filename, dpi=155)
        plt.close(fig)
        result_files.append(filename)
        print("Figure:", figdir / filename, flush=True)
    content = ["<h1>High scalar 0.5: saved endpoint graphs</h1>",
               "<p>UV RGEs frozen; fixed-order 011. Class E uses 0.25.</p>",
               '<p><a href="../overview.html">Overview</a></p>']
    content += [f'<figure><img src="{html.escape(name)}" style="max-width:100%" '
                f'alt="{html.escape(name)}"></figure>' for name in result_files]
    (figdir / "index.html").write_text("<!doctype html><html><meta charset='utf-8'><body>" +
                                       "\n".join(content) + "</body></html>", encoding="utf-8")


def _audit_running(study):
    """Simple perturbative magnitude SCREEN, not a full theoretical validation."""
    from Numerical.orchestration import RunSmallLambdaRefit as base
    source = base.ROOT / "Reports" / "output" / "display_fermion_start" / CASE / "running"
    limit = 4.0 * math.pi
    rows = []
    for file in sorted(source.glob("*/models/*/trajectory.json")):
        data = json.loads(file.read_text(encoding="utf-8"))
        scalars = data.get("eft1", {}).get("series", {})
        choices = {key: list(values) for key, values in scalars.items()
                   if key.startswith("lambda") and isinstance(values, list)}
        if not choices:
            rows.append({"model": data["model"], "mode": data["mode"],
                         "status": "NoQuarticTrajectory"})
            continue
        mag = {}
        for k, values in choices.items():
            try:
                items = [float(v) for v in values]
                mag[k] = max(abs(v) for v in items) if items and all(math.isfinite(v) for v in items) else math.inf
            except (TypeError, ValueError):
                mag[k] = math.inf
        maximum = max(mag.values())
        rows.append({"model": data["model"], "mode": data["mode"],
                     "status": "PassMagnitudeScreen" if maximum < limit else "FailMagnitudeScreen",
                     "max_abs_quartic": maximum, "per_quartic": mag})
    path = study.ROOT / CASE / "eft1_perturbativity_screen.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"criterion": "max |quartic| < 4 pi along sampled EFT1 trajectory",
                                "caveat": "Necessary practical screen only. Not a unitarity, vacuum-stability, "
                                          "Landau-pole or high-order perturbativity proof. "
                                          "Does not cover gamma-varied sensitivity trajectories.",
                                "rows": rows}, indent=2) + "\n", encoding="utf-8")
    print("Perturbativity screen:", path, "records=", len(rows), flush=True)
    return rows


def _overview(study):
    from Numerical.orchestration import RunSmallLambdaRefit as base
    display = base.ROOT / "Reports" / "output" / "display_fermion_start" / CASE
    display.mkdir(parents=True, exist_ok=True)
    links = [
        ("Endpoint table", "index.html"),
        ("Endpoint graphs (PNG)", "figures/index.html"),
        ("Optimised EFT1 + Weinberg running", "running/optimised/index.html"),
        ("Shared A(-2) EFT1 + Weinberg running", "running/shared_a_m2/index.html"),
        ("Optimised RGE dashboard (existing generator)", "reused/optimised/interactive_rge.html"),
        ("Shared A(-2) RGE dashboard (existing generator)", "reused/shared_a_m2/interactive_rge.html"),
        ("Optimised sensitivity ovals (existing generator)", "reused/optimised/sensitivity_ovals.html"),
        ("Shared A(-2) sensitivity ovals (existing generator)", "reused/shared_a_m2/sensitivity_ovals.html"),
    ]
    available = [(title, link) for title, link in links if (display / link).is_file()]
    missing = [(title, link) for title, link in links if not (display / link).is_file()]
    body = "".join(f'<li><a href="{html.escape(link)}">{html.escape(title)}</a></li>'
                   for title, link in available)
    unavailable = "".join(f'<li>{html.escape(title)} (not generated)</li>' for title, _ in missing)
    markup = ("<!doctype html><html><head><meta charset='utf-8'><title>High scalar fermion-start</title>"
              "<style>body{font:16px system-ui;max-width:900px;margin:32px auto;padding:0 24px}li{padding:7px}</style>"
              "</head><body><h1>λT3=0.5 fermion-start study</h1><p>UV running frozen. "
              "T3 A–D: 0.5; Class E: 0.25. The study uses fermion matching, "
              "EFT1 running, scalar matching and SM+Weinberg running.</p>"
              "<p><strong>Optimised</strong> means each model's fitted Yukawas; "
              "<strong>Shared A(-2)</strong> means the new high-scalar A(-2) fitted Yukawas "
              "imposed at the fermion threshold for all models.</p>"
              "<p>Standard reports are fixed-order (not coupled C12–C5), and fitted oscillation "
              "observables are not independent predictions.</p><h2>Available reports</h2><ul>" +
              body + "</ul><h2>Not yet generated</h2><ul>" + unavailable + "</ul></body></html>")
    (display / "overview.html").write_text(markup, encoding="utf-8")
    print("Overview:", display / "overview.html", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("preflight", "prepare", "fit", "compare", "display",
                                          "running", "reports", "all"))
    parser.add_argument("--mode", choices=(*MODES, "both"), default="both")
    parser.add_argument("--model", help="Full model key or short T3-A-m2 label (fit/report pilot)")
    parser.add_argument("--artifact-map", type=Path)
    parser.add_argument("--check-artifacts", action="store_true", help="Preflight checks all matching/RGE JSON")
    parser.add_argument("--sobol-points", type=int, default=256)
    parser.add_argument("--sobol-seed", type=int, default=20261022)
    parser.add_argument("--sobol-fraction", type=float, default=0.6)
    parser.add_argument("--sobol-floor", type=float, default=1e-5)
    parser.add_argument("--bound", type=float, default=0.01,
                        help="Absolute Yukawa bound; new 0.5 seeds are O(1e-4)")
    parser.add_argument("--seed-count", type=int, default=5)
    parser.add_argument("--accept-chi2", type=float, default=5.86)
    parser.add_argument("--resume", action="store_true", help="Reuse accepted fits and matching Sobol checkpoints")
    parser.add_argument("--kind", choices=("rge", "sensitivity", "all"), default="all")
    parser.add_argument("--gamma-min", type=float, default=0.5)
    parser.add_argument("--gamma-max", type=float, default=2.0)
    parser.add_argument("--points", type=int, default=11)
    parser.add_argument("--rge-points", type=int, default=41)
    parser.add_argument("--samples", type=int, default=15)
    parser.add_argument("--force-reports", action="store_true")
    args = parser.parse_args()
    if args.sobol_points < 2 or args.sobol_points & (args.sobol_points - 1):
        parser.error("--sobol-points must be a power of two >=2")
    if args.bound <= 0 or args.sobol_floor <= 0 or args.sobol_fraction <= 0:
        parser.error("--bound, --sobol-floor and --sobol-fraction must be positive")
    if args.stage == "all" and args.model:
        parser.error("--model is for pilot stages; omit it for the complete 16-model workflow")

    study = _modules()
    modes = MODES if args.mode == "both" else (args.mode,)
    if args.stage == "preflight":
        _inputs(study, args)
        return
    if args.stage in ("prepare", "all"):
        _inputs(study, args)
        _prepare(study)
    if args.stage in ("fit", "all"):
        _fit(study, args)
    if args.stage in ("compare", "all"):
        for mode in modes:
            _compare(study, args, mode)
    if args.stage in ("display", "all"):
        _ensure_compared(study, MODES)
        study.display(CASE)
        _plots(study)
        _overview(study)
    if args.stage in ("running", "all"):
        _ensure_compared(study, modes)
        _running(study, args, modes)
        _audit_running(study)
        _overview(study)
    if args.stage in ("reports", "all"):
        _ensure_compared(study, modes)
        _reports(study, args, modes)
        _overview(study)


if __name__ == "__main__":
    main()
