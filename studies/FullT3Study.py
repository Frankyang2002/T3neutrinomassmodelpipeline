"""Run the 16 ordinary T3 models using their frozen fitted-reference inputs.

The fitted numerical configurations are authoritative. This module does not
construct new Yukawa textures, alter fitted couplings, or perform optimisation.
The normal pipeline owns model construction, matching and numerical evolution.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

# Legacy four-scenario comparison constants retained for report/test compatibility.
# They do not control the 16-model fitted-reference full-study execution.
COMPARISON_SCENARIOS = {
    "smallY_smallL": (0.005, 0.1),
    "smallY_largeL": (0.005, 1.0),
    "largeY_smallL": (0.5, 0.1),
    "largeY_largeL": (0.5, 1.0),
}
LAMBDA_T3_NORMALISATION = {
    "A": 1.0, "B": 1.0, "C": 1.0, "D": 1.0, "E": 0.5,
}

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_SCRIPT = PROJECT_ROOT / "pipeline.py"
BENCHMARK_ROOT = PROJECT_ROOT / "Numerical" / "benchmarks" / "fitted_reference"
MANIFEST_PATH = BENCHMARK_ROOT / "benchmark_manifest.json"
OUTPUT_ROOT = PROJECT_ROOT / "output" / "full" / "fitted_reference"
REPORT_ROOT = PROJECT_ROOT / "Reports" / "output" / "full" / "fitted_reference"


def _read_object(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return obj


def _write_object(path: Path, obj: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_models() -> list[tuple[str, tuple[int, int, int, int], Path]]:
    manifest = _read_object(MANIFEST_PATH)
    rows = manifest.get("models")
    if not isinstance(rows, list) or len(rows) != 16:
        raise ValueError("Expected exactly 16 model entries in fitted-reference manifest")
    models: list[tuple[str, tuple[int, int, int, int], Path]] = []
    seen: set[str] = set()
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("Malformed model entry in fitted-reference manifest")
        name = str(row["model"])
        if name in seen:
            raise ValueError(f"Duplicate model in manifest: {name}")
        seen.add(name)
        rep = row["representation"]
        dims = (int(rep["d_s1"]), int(rep["d_s2"]), int(rep["d_f"]), int(rep["alpha"]))
        if rep.get("shared_scalar", False):
            raise ValueError(f"Unexpected shared-scalar model: {name}")
        config_relative = Path(str(row["config_path"]))
        if config_relative.is_absolute() or ".." in config_relative.parts:
            raise ValueError(f"Unsafe config path in manifest: {config_relative}")
        config_path = PROJECT_ROOT / config_relative
        if not config_path.is_file():
            raise FileNotFoundError(f"Missing frozen config for {name}: {config_path}")
        expected_hash = row.get("sha256")
        if expected_hash and _sha256(config_path) != expected_hash:
            raise ValueError(f"Frozen configuration hash mismatch: {config_path}")
        payload = _read_object(config_path)
        actual = payload.get("representation", {})
        if tuple(actual.get(k) for k in ("d_s1", "d_s2", "d_f", "alpha")) != dims:
            raise ValueError(f"Representation mismatch in {config_path}")
        if payload.get("benchmark_search", {}).get("enabled", False):
            raise ValueError(f"Frozen benchmark unexpectedly enables a search: {config_path}")
        if payload.get("sensitivity", {}).get("enabled", False):
            raise ValueError(f"Frozen benchmark unexpectedly enables sensitivity: {config_path}")
        models.append((name, dims, config_path))
    return models


def _common_args(args: argparse.Namespace) -> list[str]:
    values: list[str] = []
    if getattr(args, "debug_reports", False):
        values.append("--debug-reports")
    if getattr(args, "force", False):
        values.append("--force")
    groups = getattr(args, "threshold", None)
    scales = getattr(args, "threshold_scale", None)
    if groups:
        for group in groups:
            values.extend(("--threshold", *group))
    else:
        # The numerical backend uses the fermion-first intermediate scalar EFT.
        values.extend(("--threshold", "F", "--threshold", "S1", "S2"))
    if scales:
        for scale in scales:
            values.extend(("--threshold-scale", str(scale)))
    return values


def run_full_study(args: argparse.Namespace) -> int:
    """Execute one full pipeline run per frozen, individually fitted model."""
    if getattr(args, "numerical", None) is not None:
        raise ValueError(
            "--full uses the 16 frozen files in Numerical/benchmarks/fitted_reference; "
            "do not also supply --numerical. Run an individual model for a custom config."
        )
    models = _load_models()  # Validate every input before launching a costly run.
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    REPORT_ROOT.mkdir(parents=True, exist_ok=True)
    common_args = _common_args(args)
    started = time.time()
    results: list[dict[str, Any]] = []
    status = 0
    print(f"FULL FITTED-REFERENCE T3 STUDY ({len(models)} models)", flush=True)
    for name, (ds1, ds2, df, alpha), config_path in models:
        command = [
            sys.executable, str(PIPELINE_SCRIPT),
            "--dims", str(ds1), str(ds2), str(df),
            "--alpha", str(alpha),
            "--study", f"full/fitted_reference/{name}",
            "--numerical", str(config_path),
            # Avoid silently reusing an older generated_models config.
            "--reset-numerical-configs",
            *common_args,
        ]
        print("\n> " + " ".join(command), flush=True)
        run_started = time.time()
        return_code = int(subprocess.run(command, cwd=PROJECT_ROOT, check=False).returncode)
        status = max(status, return_code)
        results.append({
            "model": name,
            "config_path": str(config_path.relative_to(PROJECT_ROOT)),
            "config_sha256": _sha256(config_path),
            "return_code": return_code,
            "runtime_seconds": time.time() - run_started,
        })
    summary_path = OUTPUT_ROOT / "full_fitted_reference_summary.json"
    _write_object(summary_path, {
        "status": "Success" if status == 0 else "Failed",
        "study_mode": "fitted_reference",
        "optimizer_in_full_study": False,
        "model_count": len(models),
        "manifest_path": str(MANIFEST_PATH.relative_to(PROJECT_ROOT)),
        "results": results,
        "runtime_seconds": time.time() - started,
    })
    print(f"Fitted reference study summary: {summary_path}", flush=True)
    # Study-level postprocessing: group frozen fitted inputs only after all runs.
    # No fitted file or generated numerical config is overwritten.
    try:
        from Numerical.orchestration.RunYukawaGrouping import run_grouping
        grouping = run_grouping()
        print(f"Yukawa grouping: {len(grouping['groups'])} groups", flush=True)
    except Exception as exc:
        print(f"Yukawa grouping failed: {type(exc).__name__}: {exc}", flush=True)
        status = max(status, 1)
    return status
