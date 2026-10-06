"""Fast structural numerical-config check for the 16 ordinary T3 models.

This intentionally does not run Matchete, RGBeta, matching, ODE integration,
benchmark optimisation, or report generation.  It checks that the numerical
configuration can construct the canonical UV state and the post-F scalar-only
intermediate state for every neutral-compatible ordinary T3 model.
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from Numerical.running.FermionThresholdBoundary import (
    build_intermediate_scalar_boundary,
)
from Numerical.fitting.ScanCLI import build_uv_state_from_config
from studies.FullT3Study import (
    DEFAULT_TEMPLATE,
    MODELS,
    _fixed_comparison_config,
    _model_key,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ACTIVE_CONFIG_DIR = PROJECT_ROOT / "configs" / "generated_models"

CHECK_YUKAWA = 0.005
CHECK_SCALAR = 0.01


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object.")
    return value


def _state_from_payload(payload: dict[str, Any]):
    raw = dict(payload)
    raw.setdefault("scan", {"bindings": {}})
    return build_uv_state_from_config(SimpleNamespace(raw=raw), {}).validated()


def _check_payload(payload: dict[str, Any]) -> None:
    """Validate both UV and post-F numerical state construction."""
    uv = _state_from_payload(payload)

    intermediate = build_intermediate_scalar_boundary(
        uv,
        matching_scale_gev=float(uv.mu_gev),
    )
    intermediate.validated()


def _label(model_class: str, ds1: int, ds2: int, df: int, alpha: int) -> str:
    return (
        f"T3-{model_class} alpha={alpha:+d} "
        f"(dS1,dS2,dF)=({ds1},{ds2},{df})"
    )


def _check_comparison(template: dict[str, Any]) -> tuple[int, list[str]]:
    print("\nCOMPARISON CONFIG STRUCTURE")
    print("-" * 72)

    passed = 0
    failures: list[str] = []

    for model_class, ds1, ds2, df, alpha in MODELS:
        label = _label(model_class, ds1, ds2, df, alpha)
        try:
            payload = _fixed_comparison_config(
                template,
                ds1,
                ds2,
                df,
                alpha,
                CHECK_YUKAWA,
                CHECK_SCALAR,
            )
            _check_payload(payload)
        except Exception as exc:
            failures.append(f"{label}: {type(exc).__name__}: {exc}")
            print(f"FAIL  {label}")
            print(f"      {type(exc).__name__}: {exc}")
        else:
            passed += 1
            print(f"PASS  {label}")

    return passed, failures


def _check_saved_optimal() -> tuple[int, int, list[str]]:
    print("\nSAVED OPTIMAL CONFIGS")
    print("-" * 72)

    checked = 0
    passed = 0
    failures: list[str] = []

    for model_class, ds1, ds2, df, alpha in MODELS:
        key = _model_key(ds1, ds2, df, alpha)
        path = ACTIVE_CONFIG_DIR / f"{key}.json"
        label = _label(model_class, ds1, ds2, df, alpha)

        if not path.is_file():
            print(f"SKIP  {label}  (no saved config)")
            continue

        checked += 1
        try:
            _check_payload(_load_json(path))
        except Exception as exc:
            failures.append(f"{label}: {type(exc).__name__}: {exc}")
            print(f"FAIL  {label}")
            print(f"      {type(exc).__name__}: {exc}")
        else:
            passed += 1
            print(f"PASS  {label}")

    return checked, passed, failures


def main() -> int:
    template = _load_json(DEFAULT_TEMPLATE)

    print("=" * 72)
    print("FAST T3 NUMERICAL CONFIG CHECK")
    print("=" * 72)
    print("No Matchete / RGBeta / matching / ODE running / reports")
    print(
        "Checks: canonical UV state -> post-F scalar-only intermediate state"
    )

    comparison_passed, comparison_failures = _check_comparison(template)
    optimal_checked, optimal_passed, optimal_failures = _check_saved_optimal()

    print("\n" + "=" * 72)
    print("CHECK SUMMARY")
    print("=" * 72)
    print(
        f"Comparison structure: {comparison_passed}/{len(MODELS)} passed"
    )
    if optimal_checked:
        print(
            f"Saved optimal configs: {optimal_passed}/{optimal_checked} passed "
            f"({len(MODELS) - optimal_checked} not present)"
        )
    else:
        print("Saved optimal configs: none present; structural check skipped")

    failures = comparison_failures + optimal_failures
    if failures:
        print("\nFailures:")
        for item in failures:
            print(f"  - {item}")
        return 1

    print("\nAll checked numerical states are structurally valid.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
