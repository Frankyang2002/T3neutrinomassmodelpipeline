"""Complete 16-model numerical T3 study.

The full study has two numerical modes:
  1. optimal: each neutral-compatible ordinary T3 model uses its own persistent
     optimized benchmark config;
  2. comparison: every model is run at the same Yukawa/scalar-coupling scale
     for four small/large combinations.

Shared-scalar numerical models are intentionally not included here yet because
PipelineNumericalResults currently rejects shared_scalar records.
"""
from __future__ import annotations

from copy import deepcopy
import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from Numerical.InteractiveModelComparison import write_dashboard

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_SCRIPT = PROJECT_ROOT / "pipeline.py"
OUTPUT_ROOT = PROJECT_ROOT / "output" / "full"
REPORT_ROOT = PROJECT_ROOT / "Reports" / "output" / "full"
DEFAULT_TEMPLATE = PROJECT_ROOT / "configs" / "t3_numerical_default_B_alpha_m1.json"
ACTIVE_CONFIG_DIR = PROJECT_ROOT / "configs" / "generated_models"
CONFIG_SET_ROOT = PROJECT_ROOT / "configs" / "generated_model_sets"

# Exactly the neutral-compatible ordinary T3 class/alpha points selected by the
# current T3Model neutral-state condition over alpha=-4,...,2.
MODELS = (
    ("A", 1, 3, 2, -4), ("A", 1, 3, 2, -2), ("A", 1, 3, 2, 0),
    ("B", 2, 2, 1, -3), ("B", 2, 2, 1, -1), ("B", 2, 2, 1, 1),
    ("C", 2, 2, 3, -3), ("C", 2, 2, 3, -1), ("C", 2, 2, 3, 1),
    ("D", 3, 1, 2, -2), ("D", 3, 1, 2, 0), ("D", 3, 1, 2, 2),
    ("E", 3, 3, 2, -4), ("E", 3, 3, 2, -2), ("E", 3, 3, 2, 0), ("E", 3, 3, 2, 2),
)

COMPARISON_SCENARIOS = {
    "smallY_smallL": (0.005, 0.01),
    "smallY_largeL": (0.005, 0.25),
    "largeY_smallL": (0.5, 0.01),
    "largeY_largeL": (0.5, 0.25),
}


def _model_key(ds1: int, ds2: int, df: int, alpha: int) -> str:
    a = f"p{alpha}" if alpha >= 0 else f"m{abs(alpha)}"
    return f"T3_dS1_{ds1}_dS2_{ds2}_dF_{df}_alpha_{a}"


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    tmp.replace(path)


def _common_args(args: argparse.Namespace) -> list[str]:
    result: list[str] = []
    if getattr(args, "debug_reports", False):
        result.append("--debug-reports")
    if getattr(args, "force", False):
        result.append("--force")
    if getattr(args, "threshold", None):
        for group in args.threshold:
            result.append("--threshold")
            result.extend(group)
    if getattr(args, "threshold_scale", None):
        for scale in args.threshold_scale:
            result.extend(["--threshold-scale", str(scale)])
    return result


def _run(command: list[str]) -> int:
    print("\n> " + " ".join(command), flush=True)
    return int(subprocess.run(command, cwd=PROJECT_ROOT).returncode)


def _snapshot_active(destination: Path) -> None:
    if destination.exists():
        shutil.rmtree(destination)
    if ACTIVE_CONFIG_DIR.exists():
        shutil.copytree(ACTIVE_CONFIG_DIR, destination)
    else:
        destination.mkdir(parents=True, exist_ok=True)


def _restore_active(source: Path) -> None:
    if ACTIVE_CONFIG_DIR.exists():
        shutil.rmtree(ACTIVE_CONFIG_DIR)
    if source.exists():
        shutil.copytree(source, ACTIVE_CONFIG_DIR)
    else:
        ACTIVE_CONFIG_DIR.mkdir(parents=True, exist_ok=True)


OPTIONAL_REAL_QUARTICS = (
    "lambdaH1Adj",
    "lambdaH2Adj",
    "lambdaS1Adj",
    "lambdaS2Adj",
    "lambda12Adj",
    "lambda12Cross",
)

SPECIAL_DOUBLET_COMPLEX_QUARTICS = (
    "lambdaHHdagS2S2",
    "lambdaHHdagS1barS1bar",
    "lambdaS1bar2S2bar2",
    "lambdaS1barS2S2bar2",
    "lambdaS1S1bar2S2bar",
    "lambdaHHdagS1barS2barCross",
)


def _normalise_representation_quartics(
    ordinary: dict[str, Any],
    ds1: int,
    ds2: int,
    alpha: int,
) -> None:
    """Make the optional quartic key set match Numerical.State exactly.

    Existing values are preserved when the coupling is valid for the target
    representation. Missing required couplings are initialized to zero.
    Couplings that do not exist for the target representation are removed.
    """
    required_real = {
        "lambdaH1Adj": ds1 > 1,
        "lambdaH2Adj": ds2 > 1,
        "lambdaS1Adj": ds1 == 3,
        "lambdaS2Adj": ds2 == 3,
        "lambda12Adj": ds1 > 1 and ds2 > 1,
        "lambda12Cross": ds1 == 3 and ds2 == 3,
    }

    for name in OPTIONAL_REAL_QUARTICS:
        if required_real[name]:
            ordinary.setdefault(name, 0.0)
        else:
            ordinary.pop(name, None)

    special_doublet_case = ds1 == 2 and ds2 == 2 and alpha == -1
    for name in SPECIAL_DOUBLET_COMPLEX_QUARTICS:
        if special_doublet_case:
            ordinary.setdefault(name, {"real": 0.0, "imag": 0.0})
        else:
            ordinary.pop(name, None)


def _retarget_config(
    template: dict[str, Any],
    ds1: int,
    ds2: int,
    df: int,
    alpha: int,
) -> dict[str, Any]:
    """Retarget a payload and normalize its representation-dependent quartics."""
    cfg = deepcopy(template)
    cfg["representation"] = {
        "d_s1": ds1,
        "d_s2": ds2,
        "d_f": df,
        "alpha": alpha,
        "shared_scalar": False,
    }
    base_state = cfg.get("base_state")
    if not isinstance(base_state, dict):
        raise ValueError("Numerical config requires base_state object.")
    ordinary = base_state.get("ordinary")
    if not isinstance(ordinary, dict):
        raise ValueError("Numerical config requires base_state.ordinary object.")
    _normalise_representation_quartics(ordinary, ds1, ds2, alpha)
    return cfg


def _normalise_saved_model_config(
    path: Path,
    ds1: int,
    ds2: int,
    df: int,
    alpha: int,
) -> None:
    """Repair only the representation-dependent key set of a saved config."""
    if not path.is_file():
        return
    cfg = _retarget_config(_load_json(path), ds1, ds2, df, alpha)
    _write_json(path, cfg)

def _fixed_comparison_config(
    template: dict[str, Any],
    ds1: int,
    ds2: int,
    df: int,
    alpha: int,
    yukawa: float,
    scalar: float,
) -> dict[str, Any]:
    cfg = _retarget_config(template, ds1, ds2, df, alpha)
    ordinary = cfg["base_state"]["ordinary"]

    # Common real diagonal flavor texture. This deliberately changes only the
    # common Yukawa magnitude, leaving the same texture in all 16 models.
    for name in ("y1", "y2"):
        ordinary[name]["real"] = [
            [yukawa if i == j else 0.0 for j in range(3)]
            for i in range(3)
        ]
        ordinary[name]["imag"] = [[0.0] * 3 for _ in range(3)]

    # lambdaT3 is the scalar comparison parameter.  Other independent
    # representation-specific quartics are present when required, but fixed
    # to zero so they do not introduce an additional varying comparison axis.
    ordinary["lambdaT3"] = {"real": scalar, "imag": 0.0}

    search = cfg.setdefault("benchmark_search", {})
    search["enabled"] = False
    search["use_current_model_representation"] = False
    cfg.setdefault("sensitivity", {})["enabled"] = False
    cfg["comparison_point"] = {
        "yukawa_diagonal": yukawa,
        "lambdaT3_real": scalar,
        "representation_specific_extra_quartics": "required optional quartics fixed to 0",
        "texture": "y1=y2=y*I3; Im(y1)=Im(y2)=0; Im(lambdaT3)=0",
    }
    return cfg


def _dashboard(raw_dir: Path, report_dir: Path, title: str) -> str | None:
    try:
        path = write_dashboard(
            raw_dir,
            report_dir / "interactive_comparison.html",
            title,
        )
        print(f"Interactive comparison: {path}")
        return str(path)
    except FileNotFoundError as exc:
        print(f"Interactive comparison skipped: {exc}")
        return None


def run_full_study(args: argparse.Namespace) -> int:
    """Run all 16 ordinary neutral-compatible models in optimal + comparison modes."""
    template_path = (
        Path(args.numerical)
        if getattr(args, "numerical", None)
        else DEFAULT_TEMPLATE
    )
    if not template_path.is_absolute():
        template_path = PROJECT_ROOT / template_path
    template = _load_json(template_path)

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    REPORT_ROOT.mkdir(parents=True, exist_ok=True)
    optimal_store = CONFIG_SET_ROOT / "optimal"
    comparison_store = CONFIG_SET_ROOT / "comparison"
    common = _common_args(args)
    started = time.time()
    results: list[dict[str, Any]] = []
    status = 0

    print("=" * 72)
    print("FULL 16-MODEL T3 NUMERICAL STUDY")
    print("=" * 72)
    print("Modes: per-model optimal benchmark + four common-parameter comparisons")
    print("Shared-scalar numerical branch: deferred until its numerical state is implemented")

    # OPTIMAL: run each model with a representation-clean template.
    optimal_started = time.time()
    optimal_rc = 0
    optimal_input_dir = optimal_store / "input"
    if optimal_input_dir.exists():
        shutil.rmtree(optimal_input_dir)
    optimal_input_dir.mkdir(parents=True, exist_ok=True)

    for model_class, ds1, ds2, df, alpha in MODELS:
        key = _model_key(ds1, ds2, df, alpha)
        clean_template = _retarget_config(template, ds1, ds2, df, alpha)
        search = clean_template.setdefault("benchmark_search", {})
        search["enabled"] = True
        search["use_current_model_representation"] = False

        cfg_path = optimal_input_dir / f"{key}.json"
        _write_json(cfg_path, clean_template)

        saved = ACTIVE_CONFIG_DIR / f"{key}.json"
        if saved.is_file() and not getattr(args, "reset_numerical_configs", False):
            _normalise_saved_model_config(saved, ds1, ds2, df, alpha)

        cmd = [
            sys.executable,
            str(PIPELINE_SCRIPT),
            "--dims",
            str(ds1),
            str(ds2),
            str(df),
            "--alpha",
            str(alpha),
            "--study",
            f"full/optimal/{key}",
            "--numerical",
            str(cfg_path),
            *common,
        ]
        if getattr(args, "reset_numerical_configs", False):
            cmd.append("--reset-numerical-configs")

        rc = _run(cmd)
        optimal_rc = max(optimal_rc, rc)
        status = max(status, rc)

    _snapshot_active(optimal_store)
    optimal_dash = _dashboard(
        OUTPUT_ROOT / "optimal",
        REPORT_ROOT / "optimal",
        "T3 optimal-benchmark comparison (16 models)",
    )
    results.append(
        {
            "mode": "optimal",
            "return_code": optimal_rc,
            "runtime_seconds": time.time() - optimal_started,
            "dashboard": optimal_dash,
        }
    )

    # Preserve optimal generated configs while fixed comparison runs use the legacy
    # active generated_models directory expected by pipeline.py.
    for scenario, (yukawa, scalar) in COMPARISON_SCENARIOS.items():
        scenario_started = time.time()
        scenario_rc = 0
        input_dir = comparison_store / scenario / "input"
        resolved_dir = comparison_store / scenario / "resolved"
        if input_dir.exists():
            shutil.rmtree(input_dir)
        if resolved_dir.exists():
            shutil.rmtree(resolved_dir)
        input_dir.mkdir(parents=True, exist_ok=True)
        resolved_dir.mkdir(parents=True, exist_ok=True)

        for model_class, ds1, ds2, df, alpha in MODELS:
            key = _model_key(ds1, ds2, df, alpha)
            cfg_path = input_dir / f"{key}.json"
            _write_json(
                cfg_path,
                _fixed_comparison_config(
                    template,
                    ds1,
                    ds2,
                    df,
                    alpha,
                    yukawa,
                    scalar,
                ),
            )
            study = f"full/comparison/{scenario}/{key}"
            cmd = [
                sys.executable,
                str(PIPELINE_SCRIPT),
                "--dims",
                str(ds1),
                str(ds2),
                str(df),
                "--alpha",
                str(alpha),
                "--study",
                study,
                "--numerical",
                str(cfg_path),
                "--reset-numerical-configs",
                *common,
            ]
            rc = _run(cmd)
            scenario_rc = max(scenario_rc, rc)
            status = max(status, rc)
            generated = ACTIVE_CONFIG_DIR / f"{key}.json"
            if generated.is_file():
                shutil.copy2(generated, resolved_dir / generated.name)

        raw = OUTPUT_ROOT / "comparison" / scenario
        report = REPORT_ROOT / "comparison" / scenario
        dash = _dashboard(
            raw,
            report,
            f"T3 common-parameter comparison: {scenario} "
            f"(Y={yukawa:g}, λT3={scalar:g})",
        )
        results.append(
            {
                "mode": "comparison",
                "scenario": scenario,
                "yukawa": yukawa,
                "lambdaT3": scalar,
                "return_code": scenario_rc,
                "runtime_seconds": time.time() - scenario_started,
                "dashboard": dash,
            }
        )

    # Leave the user's active generated config directory in optimal mode.
    _restore_active(optimal_store)

    summary = {
        "status": "Success" if status == 0 else "Failed",
        "ordinary_model_count": len(MODELS),
        "shared_scalar_models_included": False,
        "comparison_definition": {
            "yukawa_texture": "y1=y2=y*I3, real",
            "scalar_parameter": "Re(lambdaT3), Im(lambdaT3)=0",
            "representation_specific_extra_quartics": "required optional quartics fixed to 0",
            "scenarios": {
                k: {"yukawa": v[0], "lambdaT3": v[1]}
                for k, v in COMPARISON_SCENARIOS.items()
            },
        },
        "results": results,
        "runtime_seconds": time.time() - started,
    }
    _write_json(OUTPUT_ROOT / "full_numerical_summary.json", summary)
    print("\n" + "=" * 72)
    print("FULL NUMERICAL SUMMARY")
    print("=" * 72)
    print(f"Status: {summary['status']}")
    print(f"Models: {len(MODELS)} ordinary neutral-compatible")
    print(f"Summary: {OUTPUT_ROOT / 'full_numerical_summary.json'}")
    print(f"Interactive reports: {REPORT_ROOT}")
    return status
