"""Complete 16-model comparison-only numerical T3 study."""
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

from Numerical.plotting.InteractiveModelComparison import write_dashboard


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_SCRIPT = PROJECT_ROOT / "pipeline.py"
OUTPUT_ROOT = PROJECT_ROOT / "output" / "full"
REPORT_ROOT = PROJECT_ROOT / "Reports" / "output" / "full"
DEFAULT_TEMPLATE = PROJECT_ROOT / "configs" / "t3_numerical_default_B_alpha_m1.json"
ACTIVE_CONFIG_DIR = PROJECT_ROOT / "configs" / "generated_models"
CONFIG_SET_ROOT = PROJECT_ROOT / "configs" / "generated_model_sets"

MODELS = (
    ("A", 1, 3, 2, -4), ("A", 1, 3, 2, -2), ("A", 1, 3, 2, 0),
    ("B", 2, 2, 1, -3), ("B", 2, 2, 1, -1), ("B", 2, 2, 1, 1),
    ("C", 2, 2, 3, -3), ("C", 2, 2, 3, -1), ("C", 2, 2, 3, 1),
    ("D", 3, 1, 2, -2), ("D", 3, 1, 2, 0), ("D", 3, 1, 2, 2),
    ("E", 3, 3, 2, -4), ("E", 3, 3, 2, -2), ("E", 3, 3, 2, 0),
    ("E", 3, 3, 2, 2),
)

COMPARISON_SCENARIOS = {
    "smallY_smallL": (0.005, 0.1),
    "smallY_largeL": (0.005, 1.0),
    "largeY_smallL": (0.5, 0.1),
    "largeY_largeL": (0.5, 1.0),
}

OPTIONAL_REAL_QUARTICS = (
    "lambdaH1Adj", "lambdaH2Adj", "lambdaS1Adj", "lambdaS2Adj",
    "lambda12Adj", "lambda12Cross",
)
SPECIAL_DOUBLET_COMPLEX_QUARTICS = (
    "lambdaHHdagS2S2",
    "lambdaHHdagS1barS1bar",
    "lambdaS1bar2S2bar2",
    "lambdaS1barS2S2bar2",
    "lambdaS1S1bar2S2bar",
    "lambdaHHdagS1barS2barCross",
)
COMPARISON_Y1_TEXTURE = (
    (1, .2, -.1),
    (.2, .8, .3),
    (-.1, .3, .6),
)
COMPARISON_Y2_TEXTURE = (
    (.7, -.3, .2),
    (.4, 1, -.2),
    (.1, .3, .8),
)


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
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _common_args(args: argparse.Namespace) -> list[str]:
    result: list[str] = []
    if getattr(args, "debug_reports", False):
        result.append("--debug-reports")
    if getattr(args, "force", False):
        result.append("--force")
    if getattr(args, "threshold", None):
        for group in args.threshold:
            result += ["--threshold", *group]
    if getattr(args, "threshold_scale", None):
        for value in args.threshold_scale:
            result += ["--threshold-scale", str(value)]
    return result


def _run(command: list[str]) -> int:
    print("\n> " + " ".join(command), flush=True)
    return int(
        subprocess.run(
            command,
            cwd=PROJECT_ROOT,
        ).returncode
    )


def _normalise_representation_quartics(
    ordinary: dict[str, Any],
    ds1: int,
    ds2: int,
    alpha: int,
) -> None:
    required = {
        "lambdaH1Adj": ds1 > 1,
        "lambdaH2Adj": ds2 > 1,
        "lambdaS1Adj": ds1 == 3,
        "lambdaS2Adj": ds2 == 3,
        "lambda12Adj": ds1 > 1 and ds2 > 1,
        "lambda12Cross": ds1 == 3 and ds2 == 3,
    }
    for name in OPTIONAL_REAL_QUARTICS:
        if required[name]:
            ordinary.setdefault(name, 0.0)
        else:
            ordinary.pop(name, None)

    special = ds1 == 2 and ds2 == 2 and alpha == -1
    for name in SPECIAL_DOUBLET_COMPLEX_QUARTICS:
        if special:
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
    config = deepcopy(template)
    config["representation"] = {
        "d_s1": ds1,
        "d_s2": ds2,
        "d_f": df,
        "alpha": alpha,
        "shared_scalar": False,
    }

    base_state = config.get("base_state")
    if not isinstance(base_state, dict):
        raise ValueError("Numerical config requires base_state object.")

    ordinary = base_state.get("ordinary")
    if not isinstance(ordinary, dict):
        raise ValueError("Numerical config requires base_state.ordinary object.")

    _normalise_representation_quartics(
        ordinary,
        ds1,
        ds2,
        alpha,
    )
    return config


def _scaled_real_texture(
    texture: tuple[tuple[float, ...], ...],
    scale: float,
) -> list[list[float]]:
    return [
        [scale * float(value) for value in row]
        for row in texture
    ]


def _fixed_comparison_config(
    template: dict[str, Any],
    ds1: int,
    ds2: int,
    df: int,
    alpha: int,
    yukawa: float,
    scalar: float,
) -> dict[str, Any]:
    config = _retarget_config(
        template,
        ds1,
        ds2,
        df,
        alpha,
    )
    ordinary = config["base_state"]["ordinary"]

    ordinary["y1"]["real"] = _scaled_real_texture(
        COMPARISON_Y1_TEXTURE,
        yukawa,
    )
    ordinary["y2"]["real"] = _scaled_real_texture(
        COMPARISON_Y2_TEXTURE,
        yukawa,
    )
    ordinary["y1"]["imag"] = [[0.0] * 3 for _ in range(3)]
    ordinary["y2"]["imag"] = [[0.0] * 3 for _ in range(3)]
    ordinary["lambdaT3"] = {
        "real": scalar,
        "imag": 0.0,
    }

    # The full study is deliberately comparison-only.  The automatic Sobol +
    # local least-squares benchmark search remains available elsewhere, but it
    # is never invoked by this orchestrator.
    search = config.setdefault("benchmark_search", {})
    search["enabled"] = False
    search["use_current_model_representation"] = False

    config.setdefault("sensitivity", {})["enabled"] = False

    config["comparison_point"] = {
        "yukawa_scale": yukawa,
        "lambdaT3_real": scalar,
        "representation_specific_extra_quartics":
            "required optional quartics fixed to 0",
        "texture": (
            "fixed real non-diagonal y1/y2 textures scaled by common Y; "
            "Im(y1)=Im(y2)=0; Im(lambdaT3)=0"
        ),
        "y1_texture": [list(row) for row in COMPARISON_Y1_TEXTURE],
        "y2_texture": [list(row) for row in COMPARISON_Y2_TEXTURE],
    }
    return config


def _dashboard(
    raw: Path,
    report: Path,
    title: str,
) -> str | None:
    try:
        path = write_dashboard(
            raw,
            report / "interactive_comparison.html",
            title,
        )
        print(f"Interactive comparison: {path}")
        return str(path)
    except FileNotFoundError as exc:
        print(f"Interactive comparison skipped: {exc}")
        return None


def run_full_study(args: argparse.Namespace) -> int:
    """Run the standard full study using only fixed comparison benchmarks.

    Optimisation is intentionally detached from this workflow.  The optimizer
    and Sobol benchmark-search modules are preserved and can still be invoked
    explicitly, but ``pipeline.py --full`` no longer runs them.
    """

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

    comparison_store = CONFIG_SET_ROOT / "comparison"
    common = _common_args(args)
    started = time.time()
    results: list[dict[str, Any]] = []
    status = 0

    print("=" * 72)
    print("FULL 16-MODEL T3 COMPARISON STUDY")
    print("=" * 72)
    print("Mode: four fixed common-parameter comparison benchmarks")
    print("Automatic optimal/Sobol benchmark search: detached from --full")
    print("Comparison scales: Y={0.005,0.5}, lambdaT3={0.1,1.0}")
    print("Mass thresholds: MF=100 TeV, MS=1 TeV; UV=1e7 GeV; low=100 GeV")

    for scenario, (yukawa, scalar) in COMPARISON_SCENARIOS.items():
        scenario_started = time.time()
        scenario_rc = 0

        input_dir = comparison_store / scenario / "input"
        resolved_dir = comparison_store / scenario / "resolved"

        for directory in (input_dir, resolved_dir):
            if directory.exists():
                shutil.rmtree(directory)
            directory.mkdir(parents=True, exist_ok=True)

        for _, ds1, ds2, df, alpha in MODELS:
            key = _model_key(ds1, ds2, df, alpha)
            config_path = input_dir / f"{key}.json"

            _write_json(
                config_path,
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

            command = [
                sys.executable,
                str(PIPELINE_SCRIPT),
                "--dims", str(ds1), str(ds2), str(df),
                "--alpha", str(alpha),
                "--study", f"full/comparison/{scenario}/{key}",
                "--numerical", str(config_path),
                "--reset-numerical-configs",
                *common,
            ]
            rc = _run(command)
            scenario_rc = max(scenario_rc, rc)
            status = max(status, rc)

            generated = ACTIVE_CONFIG_DIR / f"{key}.json"
            if generated.is_file():
                shutil.copy2(
                    generated,
                    resolved_dir / generated.name,
                )

        dashboard = _dashboard(
            OUTPUT_ROOT / "comparison" / scenario,
            REPORT_ROOT / "comparison" / scenario,
            (
                f"T3 common-parameter comparison: {scenario} "
                f"(Y={yukawa:g}, λT3={scalar:g})"
            ),
        )
        results.append(
            {
                "mode": "comparison",
                "scenario": scenario,
                "yukawa": yukawa,
                "lambdaT3": scalar,
                "return_code": scenario_rc,
                "runtime_seconds": time.time() - scenario_started,
                "dashboard": dashboard,
            }
        )

    summary = {
        "status": "Success" if status == 0 else "Failed",
        "study_mode": "comparison_only",
        "optimizer_in_full_study": False,
        "ordinary_model_count": len(MODELS),
        "shared_scalar_models_included": False,
        "comparison_definition": {
            "yukawa_texture":
                "fixed real non-diagonal y1/y2 textures scaled by common Y",
            "scalar_parameter":
                "Re(lambdaT3), Im(lambdaT3)=0",
            "representation_specific_extra_quartics":
                "required optional quartics fixed to 0",
            "scenarios": {
                key: {
                    "yukawa": values[0],
                    "lambdaT3": values[1],
                }
                for key, values in COMPARISON_SCENARIOS.items()
            },
        },
        "results": results,
        "runtime_seconds": time.time() - started,
    }

    summary_path = OUTPUT_ROOT / "full_numerical_summary.json"
    _write_json(summary_path, summary)

    print("\n" + "=" * 72)
    print("FULL NUMERICAL SUMMARY")
    print("=" * 72)
    print(f"Status: {summary['status']}")
    print(f"Models: {len(MODELS)} ordinary neutral-compatible")
    print("Optimal benchmark search: not run by --full")
    print(f"Summary: {summary_path}")
    print(f"Interactive reports: {REPORT_ROOT}")
    return status
