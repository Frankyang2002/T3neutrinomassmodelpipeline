"""Numerical-only T1/T2/T3 sensitivity envelopes about the COMMON-INPUT RGE points.

Requires completed RunCommonInputRGE output and existing RGBeta/matching JSON.
No symbolic regeneration, fitting, or modifications of fitted results.

Usage from repository root:
    python -m Numerical.orchestration.RunCommonInputEnvelopes
    python -m Numerical.orchestration.RunCommonInputEnvelopes --run --limit 2 --points 7 --gamma-min 0.5 --gamma-max 2
    python -m Numerical.orchestration.RunCommonInputEnvelopes --run
    python -m Numerical.orchestration.RunCommonInputEnvelopes --display-only
    python -m Numerical.orchestration.RunCommonInputEnvelopes --run --points 25 --force
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from Numerical.fitting.BenchmarkSensitivity import gamma_grid, run_t3_sensitivity_scan
from Numerical.orchestration.ResumePipelineNumericalResults import (
    PROJECT_ROOT, OUTPUT_ROOT, REQUIRED, discover_cases,
)
from Numerical.plotting.BuildFittedModelEnvelopes import build as build_envelopes

SCENARIO = "smallY_smallL"
COMMON_ROOT = PROJECT_ROOT / "Reports" / "output" / "full" / "common_input_rge"
DEST_ROOT = PROJECT_ROOT / "Reports" / "output" / "full" / "common_input_sensitivity"
DISPLAY_ROOT = PROJECT_ROOT / "Reports" / "output" / "display"
TARGET = PROJECT_ROOT / "configs" / "nufit_6_1_ic24_no.json"
DIRECTIONS = ("T1", "T2", "T3")
OBS = ("delta_m21_sq_ev2", "delta_m3l_sq_ev2", "sin2_theta12", "sin2_theta13", "sin2_theta23")


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def save_json(path: Path, obj: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def signature(config: dict, paths: list[Path], target: Path, params: dict) -> str:
    h = hashlib.sha256(json.dumps({"config": config, "scan": params}, sort_keys=True).encode("utf-8"))
    for path in [*paths, target]:
        stat = path.stat()
        h.update(f"{path.resolve()}:{stat.st_size}:{stat.st_mtime_ns}".encode("utf-8"))
    return h.hexdigest()


def verify_centre(scan: dict, points: int) -> dict:
    centres = [scan["directions"][d]["points"][points // 2] for d in DIRECTIONS]
    if any(abs(float(c["gamma"]) - 1.0) > 1e-12 for c in centres):
        raise RuntimeError("The gamma=1 centre is missing")
    if any(c.get("status") != "Success" for c in centres):
        raise RuntimeError("At least one gamma=1 baseline calculation failed")
    ref = centres[0]["prediction"]
    for c in centres[1:]:
        for observable in OBS:
            a, b = float(ref[observable]), float(c["prediction"][observable])
            if not math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-12):
                raise RuntimeError(f"Inconsistent gamma=1 centre: {observable}")
    return centres[0]


def make_common_envelopes(source: Path, output: Path) -> None:
    """Reuse original display geometry, relabel it to avoid false fit claims."""
    data = build_envelopes(source, output, include_warnings=False)
    html = output.read_text(encoding="utf-8")
    html = html.replace("Individually fitted T3 model sensitivity", "Common-input T3 model sensitivity")
    html = html.replace(
        "Each model uses its own fitted λT3, Yukawa matrices and other fitted parameters. Data from saved scans only — no RGE recomputation.",
        "Unrefitted common-input benchmark. T1/T2/T3 numerical variations; each curve is centred on its own common-input prediction. Class E uses half the reference lambdaT3. No symbolic recomputation."
    )
    html = html.replace("Fitted T3 prediction-space sensitivity", "Common-input T3 prediction-space sensitivity")
    html = html.replace("Fitted centres", "Common-input centres")
    html = html.replace("about each fitted centre", "about each common-input centre")
    html = html.replace(" fitted models plotted; ", " common-input models plotted; ")
    # Inject interaction into the SAME script scope as the existing manual/bounds/draw
    # variables. No changes to the shared fitted-envelope plotter.
    pan_script = r"""
// Right-mouse drag to pan the visible plot; preserve the current zoom span.
let rightPan = null;
canvas.style.cursor = 'crosshair';
canvas.addEventListener('contextmenu', event => event.preventDefault());
canvas.addEventListener('pointerdown', event => {
    if (event.button !== 2) return;
    event.preventDefault();
    const rect = canvas.getBoundingClientRect();
    const left = 92, right = 20, top = 28, bottom = 75;
    const px = event.clientX - rect.left, py = event.clientY - rect.top;
    if (px < left || px > rect.width - right || py < top || py > rect.height - bottom) return;
    rightPan = {x: event.clientX, y: event.clientY,
                view: (manual || bounds()).slice(),
                width: Math.max(1, rect.width - left - right),
                height: Math.max(1, rect.height - top - bottom)};
    canvas.setPointerCapture(event.pointerId);
    canvas.style.cursor = 'grabbing';
});
canvas.addEventListener('pointermove', event => {
    if (!rightPan) return;
    const p = rightPan;
    const dx = (event.clientX - p.x) / p.width * (p.view[1] - p.view[0]);
    const dy = (event.clientY - p.y) / p.height * (p.view[3] - p.view[2]);
    manual = [p.view[0] - dx, p.view[1] - dx, p.view[2] + dy, p.view[3] + dy];
    draw();
});
function endRightPan() { rightPan = null; canvas.style.cursor = 'crosshair'; }
canvas.addEventListener('pointerup', endRightPan);
canvas.addEventListener('pointercancel', endRightPan);
canvas.addEventListener('lostpointercapture', endRightPan);
"""
    marker = '</script></body></html>'
    if html.count(marker) != 1:
        raise RuntimeError('Unexpected envelope HTML script structure; cannot install right-drag pan')
    html = html.replace(marker, pan_script + marker)
    output.write_text(html, encoding="utf-8")
    print(f"Envelope HTML: {output} ({len(data['models'])} models)")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", action="store_true", help="Run numerical sensitivity scans (default inventory only)")
    p.add_argument("--display-only", action="store_true", help="Rebuild HTML from saved sensitivity scans; no numerical work")
    p.add_argument("--limit", type=int, default=0, help="Only the first N sorted models (0: all)")
    p.add_argument("--force", action="store_true", help="Recompute existing sensitivity scans")
    p.add_argument("--gamma-min", type=float, default=0.1)
    p.add_argument("--gamma-max", type=float, default=10.0)
    p.add_argument("--points", type=int, default=25)
    p.add_argument("--spacing", choices=("log", "linear"), default="log")
    args = p.parse_args(argv)
    if args.display_only:
        if args.run:
            p.error("--display-only and --run are mutually exclusive")
        make_common_envelopes(DEST_ROOT / SCENARIO,
                              DISPLAY_ROOT / "common_input_model_envelopes_all.html")
        return 0
    if args.limit < 0:
        p.error("--limit must be nonnegative")
    try:
        gamma_grid(args.gamma_min, args.gamma_max, args.points, args.spacing)
    except ValueError as e:
        p.error(str(e))
    if not TARGET.is_file():
        p.error(f"Missing oscillation target for observable extraction: {TARGET}")
    cases = {}
    for scenario, model, case in discover_cases(OUTPUT_ROOT):
        if scenario != SCENARIO:
            continue
        if model in cases:
            p.error(f"Duplicate symbolic artifact directory for {model}")
        cases[model] = case
    selected = sorted(cases.items())[:args.limit or None]
    output = DEST_ROOT / SCENARIO
    rows = []
    options = {"gamma_min": args.gamma_min, "gamma_max": args.gamma_max,
               "n_points": args.points, "spacing": args.spacing}
    for model, case in selected:
        baseline = COMMON_ROOT / SCENARIO / model
        source = baseline / "common_input_config.json"
        diagnostics = baseline / "running_diagnostics.json"
        destination = output / model / "sensitivity_diagnostics.json"
        row = {"model": model, "status": "Ready", "verified_chi2": "",
               "scan_json": str(destination), "n_failed_points": "", "message": "",
               "source_config": str(source)}
        try:
            if not source.is_file() or not diagnostics.is_file():
                raise FileNotFoundError(f"RunCommonInputRGE baseline missing for {model}")
            baseline_diag = read_json(diagnostics)
            if baseline_diag.get("status") != "Success" or "common_input_rge" not in baseline_diag:
                raise ValueError("Not a successful common-input diagnostic")
            artifact_paths = [case / "data" / name for name in REQUIRED]
            missing = [str(v) for v in artifact_paths if not v.is_file()]
            if missing:
                raise FileNotFoundError("Missing existing numerical artifacts: " + ", ".join(missing))
            config = read_json(source)
            config["benchmark_search"] = {**config.get("benchmark_search", {}), "enabled": False}
            config["scan"] = {"bindings": {}}
            config["sensitivity"] = {"enabled": True, **options}
            key = signature(config, artifact_paths, TARGET, options)
            if not args.run:
                row["status"] = "AlreadyComplete" if destination.is_file() else "Ready"
                rows.append(row)
                print(f"{row['status']:18} {model}", flush=True)
                continue
            if destination.is_file() and not args.force:
                scan = read_json(destination)
                if scan.get("common_input_sensitivity", {}).get("signature") != key:
                    raise RuntimeError("Existing scan has changed parameters or artifacts; use --force")
            else:
                scan = run_t3_sensitivity_scan(
                    base_payload=config, uv_rgbeta_path=artifact_paths[0],
                    intermediate_rgbeta_path=artifact_paths[1],
                    final_weinberg_path=artifact_paths[2],
                    oscillation_target_path=TARGET,
                )
                scan["common_input_sensitivity"] = {
                    "signature": key, "source_config": str(source),
                    "source_diagnostics": str(diagnostics),
                    "note": "Unrefitted common-input gamma variations; not confidence regions",
                }
                save_json(destination, scan)
            centre = verify_centre(scan, args.points)
            row["verified_chi2"] = float(centre["chi2"])
            row["n_failed_points"] = int(scan.get("n_failed_points", 0))
            # Common-input baselines are not required to fit the NuFIT central point.
            row["status"] = "Success"
            if row["n_failed_points"]:
                row["message"] = f"{row['n_failed_points']} nonbaseline scan points failed"
        except Exception as exc:
            row["status"] = "Failed"
            row["message"] = f"{type(exc).__name__}: {exc}"
        rows.append(row)
        print(f"{row['status']:18} {model}: {row['message']}", flush=True)
    if not args.run:
        print("Inventory only. Add --run to compute sensitivity scans.")
        return 0
    save_json(output / "comparison_summary.json", rows)
    successful = [r for r in rows if r["status"] == "Success"]
    if successful:
        try:
            make_common_envelopes(output, DISPLAY_ROOT / "common_input_model_envelopes_all.html")
        except Exception as exc:
            print(f"Dashboard error: {type(exc).__name__}: {exc}")
            return 1
    print(f"Success: {len(successful)}/{len(rows)}; results: {output}")
    return int(any(r["status"] == "Failed" for r in rows))


if __name__ == "__main__":
    raise SystemExit(main())
