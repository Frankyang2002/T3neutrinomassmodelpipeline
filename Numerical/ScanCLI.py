"""Command-line driver for reproducible ordinary-T3 parameter scans.

Usage
-----

    python -m Numerical.ScanCLI configs/my_scan.json
    python -m Numerical.ScanCLI configs/my_scan.json --dry-run

The configuration file specifies

- the representation and base UV numerical state;
- parameter bindings into that base state;
- grid or random scan generation;
- explicit threshold and low-energy scales;
- the UV/intermediate RGBeta JSON payloads;
- the authoritative final-Weinberg JSON;
- the external oscillation-fit target;
- the output JSON path.

No NuFIT values or matched C5 expressions are hard-coded here.

Current production scope
------------------------
The concrete ``FinalC5TrajectoryAdapter`` currently supports the ordinary
split-scalar branch. Shared-scalar scans are rejected by this CLI rather than
silently routed through an unvalidated matching convention.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from Numerical.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.OscillationFit import OscillationFitTarget
from Numerical.ParameterScan import (
    ParameterScanResult,
    ScanParameter,
    grid_parameter_points,
    make_t3_point_evaluator,
    random_parameter_points,
    run_parameter_scan,
    write_scan_result,
)
from Numerical.RGBetaEvaluator import load_rgbeta_payload
from Numerical.State import (
    SMNumericalState,
    SharedT3UVState,
    T3Representation,
    T3UVState,
)


UVState = T3UVState | SharedT3UVState


@dataclass(frozen=True)
class ResolvedScanConfig:
    """Validated scan configuration with absolute file paths."""

    source_path: Path
    raw: dict[str, Any]
    uv_rgbeta_path: Path
    eft1_rgbeta_path: Path
    final_weinberg_path: Path
    oscillation_target_path: Path
    output_path: Path
    checkpoint_path: Path


def _read_json_object(path: Path) -> dict[str, Any]:
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))

    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object.")

    return payload


def _resolve_path(base_dir: Path, value: Any, name: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty path string.")

    path = Path(value)

    if not path.is_absolute():
        path = base_dir / path

    return path.resolve()


def load_scan_config(path: Path) -> ResolvedScanConfig:
    """Load and validate top-level scan configuration."""

    source_path = Path(path).resolve()
    raw = _read_json_object(source_path)

    files = raw.get("files")
    if not isinstance(files, dict):
        raise ValueError("scan config requires a 'files' object.")

    output = raw.get("output")
    if not isinstance(output, dict):
        raise ValueError("scan config requires an 'output' object.")

    base_dir = source_path.parent

    resolved = ResolvedScanConfig(
        source_path=source_path,
        raw=raw,
        uv_rgbeta_path=_resolve_path(
            base_dir,
            files.get("uv_rgbeta"),
            "files.uv_rgbeta",
        ),
        eft1_rgbeta_path=_resolve_path(
            base_dir,
            files.get("eft1_rgbeta"),
            "files.eft1_rgbeta",
        ),
        final_weinberg_path=_resolve_path(
            base_dir,
            files.get("final_weinberg"),
            "files.final_weinberg",
        ),
        oscillation_target_path=_resolve_path(
            base_dir,
            files.get("oscillation_target"),
            "files.oscillation_target",
        ),
        output_path=_resolve_path(
            base_dir,
            output.get("result_json"),
            "output.result_json",
        ),
        checkpoint_path=_resolve_path(
            base_dir,
            output.get(
                "checkpoint_json",
                str(
                    Path(output.get("result_json")).with_name(
                        Path(output.get("result_json")).stem
                        + ".checkpoint.json"
                    )
                ),
            ),
            "output.checkpoint_json",
        ),
    )

    _validate_static_config(resolved)
    return resolved


def _validate_static_config(config: ResolvedScanConfig) -> None:
    raw = config.raw

    representation = raw.get("representation")
    if not isinstance(representation, dict):
        raise ValueError(
            "scan config requires a 'representation' object."
        )

    if bool(representation.get("shared_scalar", False)):
        raise ValueError(
            "ScanCLI currently supports only the ordinary split-scalar branch "
            "because FinalC5TrajectoryAdapter has no validated shared-scalar adapter."
        )

    scales = raw.get("scales")
    if not isinstance(scales, dict):
        raise ValueError("scan config requires a 'scales' object.")

    required_scales = (
        "mu_fermion_threshold_gev",
        "mu_scalar_threshold_gev",
        "mu_low_gev",
    )

    for name in required_scales:
        value = float(scales[name])
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(
                f"scales.{name} must be finite and positive."
            )

    base_state = raw.get("base_state")
    if not isinstance(base_state, dict):
        raise ValueError(
            "scan config requires a 'base_state' object."
        )

    scan = raw.get("scan")
    if not isinstance(scan, dict):
        raise ValueError("scan config requires a 'scan' object.")

    mode = str(scan.get("mode", "")).lower()
    if mode not in {"grid", "random"}:
        raise ValueError(
            "scan.mode must be either 'grid' or 'random'."
        )


def _path_segments(path: str) -> tuple[str, ...]:
    parts = tuple(
        part.strip()
        for part in str(path).split(".")
        if part.strip()
    )

    if not parts:
        raise ValueError("Parameter binding target cannot be empty.")

    return parts


def _set_bound_value(
    root: dict[str, Any],
    target: str,
    value: float,
) -> None:
    """Set one scalar inside nested JSON-like dict/list state data.

    List indices are expressed as numeric dotted path components, e.g.

        ordinary.y1.real.0.1
    """

    parts = _path_segments(target)
    current: Any = root

    for part in parts[:-1]:
        if isinstance(current, dict):
            if part not in current:
                raise KeyError(
                    f"Parameter binding path does not exist: {target}"
                )
            current = current[part]
            continue

        if isinstance(current, list):
            try:
                index = int(part)
            except ValueError as exc:
                raise KeyError(
                    f"Expected list index in binding path: {target}"
                ) from exc

            if index < 0 or index >= len(current):
                raise IndexError(
                    f"List index out of range in binding path: {target}"
                )
            current = current[index]
            continue

        raise KeyError(
            f"Cannot descend through binding path: {target}"
        )

    last = parts[-1]

    if isinstance(current, dict):
        if last not in current:
            raise KeyError(
                f"Parameter binding path does not exist: {target}"
            )
        current[last] = float(value)
        return

    if isinstance(current, list):
        try:
            index = int(last)
        except ValueError as exc:
            raise KeyError(
                f"Expected list index in binding path: {target}"
            ) from exc

        if index < 0 or index >= len(current):
            raise IndexError(
                f"List index out of range in binding path: {target}"
            )

        current[index] = float(value)
        return

    raise KeyError(
        f"Cannot assign through binding path: {target}"
    )


def _complex_scalar(spec: Any, name: str) -> complex:
    if isinstance(spec, (int, float)):
        value = complex(float(spec), 0.0)
    elif isinstance(spec, dict):
        value = complex(
            float(spec.get("real", 0.0)),
            float(spec.get("imag", 0.0)),
        )
    else:
        raise ValueError(
            f"{name} must be a number or {{real, imag}} object."
        )

    if not (
        np.isfinite(value.real)
        and np.isfinite(value.imag)
    ):
        raise ValueError(f"{name} must be finite.")

    return value


def _complex_matrix(spec: Any, name: str) -> np.ndarray:
    if not isinstance(spec, dict):
        raise ValueError(
            f"{name} must contain 'real' and optional 'imag' arrays."
        )

    real = np.asarray(spec.get("real"), dtype=float)

    if "imag" in spec:
        imag = np.asarray(spec["imag"], dtype=float)
    else:
        imag = np.zeros_like(real)

    if real.shape != (3, 3) or imag.shape != (3, 3):
        raise ValueError(
            f"{name} real/imag arrays must both have shape (3, 3)."
        )

    matrix = real + 1j * imag

    if not np.all(np.isfinite(matrix.real)) or not np.all(
        np.isfinite(matrix.imag)
    ):
        raise ValueError(f"{name} must contain finite entries.")

    return matrix


def _build_representation(raw: Mapping[str, Any]) -> T3Representation:
    return T3Representation(
        d_s1=int(raw["d_s1"]),
        d_s2=int(raw["d_s2"]),
        d_f=int(raw["d_f"]),
        alpha=int(raw["alpha"]),
        shared_scalar=bool(raw.get("shared_scalar", False)),
    ).validated()


def _metadata_int(
    metadata: Mapping[str, Any],
    key: str,
    *,
    payload_name: str,
) -> int:
    if key not in metadata:
        raise ValueError(
            f"{payload_name} metadata is missing required key {key!r}."
        )

    try:
        return int(metadata[key])
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"{payload_name} metadata key {key!r} must be an integer."
        ) from exc


def _validate_rgbeta_payload_compatibility(
    *,
    payload: Mapping[str, Any],
    representation: T3Representation,
    payload_name: str,
    eft1: bool,
) -> None:
    """Cross-check a loaded RGBeta payload against the scan representation.

    ``load_rgbeta_payload`` verifies the generic JSON contract. This function
    adds the model/stage compatibility check required by the scan CLI so a
    dry-run cannot succeed with a representation-mismatched RGBeta export.
    """

    metadata = payload.get("metadata")

    if not isinstance(metadata, Mapping):
        raise ValueError(
            f"{payload_name} requires a metadata object."
        )

    if bool(metadata.get("SharedScalar", False)) != bool(
        representation.shared_scalar
    ):
        raise ValueError(
            f"{payload_name} SharedScalar metadata does not match the scan "
            "representation."
        )

    expected = {
        "dS1": representation.d_s1,
        "dS2": representation.d_s2,
        "dF": representation.d_f,
        "alpha": representation.alpha,
    }

    mismatches: list[str] = []

    for key, expected_value in expected.items():
        actual = _metadata_int(
            metadata,
            key,
            payload_name=payload_name,
        )

        if actual != expected_value:
            mismatches.append(
                f"{key}={actual} (config expects {expected_value})"
            )

    if mismatches:
        raise ValueError(
            f"{payload_name} representation mismatch: "
            + ", ".join(mismatches)
        )

    if eft1:
        integrated = metadata.get("IntegratedField")

        if integrated != "F":
            raise ValueError(
                f"{payload_name} must describe the post-F EFT1 stage; "
                f"IntegratedField={integrated!r}."
            )

        active = metadata.get("ActiveBSMFields")

        if not isinstance(active, list):
            raise ValueError(
                f"{payload_name} metadata requires ActiveBSMFields list."
            )

        if set(active) != {"S1", "S2"}:
            raise ValueError(
                f"{payload_name} must have active BSM fields S1 and S2; "
                f"got {active!r}."
            )


def build_uv_state_from_config(
    config: ResolvedScanConfig,
    parameters: Mapping[str, float],
) -> UVState:
    """Apply scan bindings to the base state and build a validated UV state."""

    raw = config.raw
    mutable = deepcopy(raw["base_state"])

    scan = raw["scan"]
    bindings = scan.get("bindings", {})

    if not isinstance(bindings, dict):
        raise ValueError("scan.bindings must be an object.")

    for parameter_name, value in parameters.items():
        if parameter_name not in bindings:
            raise KeyError(
                f"No scan.bindings entry for parameter {parameter_name!r}."
            )

        _set_bound_value(
            mutable,
            str(bindings[parameter_name]),
            float(value),
        )

    representation = _build_representation(
        raw["representation"]
    )

    if representation.shared_scalar:
        raise ValueError(
            "Shared-scalar state construction is outside current ScanCLI scope."
        )

    sm_raw = mutable.get("sm")
    ordinary = mutable.get("ordinary")

    if not isinstance(sm_raw, dict):
        raise ValueError("base_state.sm must be an object.")

    if not isinstance(ordinary, dict):
        raise ValueError(
            "base_state.ordinary must be an object."
        )

    sm = SMNumericalState(
        gY=float(sm_raw["gY"]),
        g2=float(sm_raw["g2"]),
        g3=float(sm_raw["g3"]),
        lambdaH=float(sm_raw["lambdaH"]),
        yu=_complex_matrix(sm_raw["yu"], "base_state.sm.yu"),
        yd=_complex_matrix(sm_raw["yd"], "base_state.sm.yd"),
        ye=_complex_matrix(sm_raw["ye"], "base_state.sm.ye"),
    ).validated()

    required = {
        "mu_gev": float(mutable["mu_gev"]),
        "representation": representation,
        "sm": sm,
        "y1": _complex_matrix(
            ordinary["y1"],
            "base_state.ordinary.y1",
        ),
        "y2": _complex_matrix(
            ordinary["y2"],
            "base_state.ordinary.y2",
        ),
        "MF": _complex_matrix(
            ordinary["MF"],
            "base_state.ordinary.MF",
        ),
        "mS1Sq": float(ordinary["mS1Sq"]),
        "mS2Sq": float(ordinary["mS2Sq"]),
        "lambdaS1": float(ordinary["lambdaS1"]),
        "lambdaS2": float(ordinary["lambdaS2"]),
        "lambdaH1": float(ordinary["lambdaH1"]),
        "lambdaH2": float(ordinary["lambdaH2"]),
        "lambda12": float(ordinary["lambda12"]),
        "lambdaT3": _complex_scalar(
            ordinary["lambdaT3"],
            "base_state.ordinary.lambdaT3",
        ),
    }

    optional_names = (
        "lambdaH1Adj",
        "lambdaH2Adj",
        "lambdaS1Adj",
        "lambdaS2Adj",
        "lambda12Adj",
        "lambda12Cross",
        "lambdaHHdagS2S2",
        "lambdaHHdagS1barS1bar",
        "lambdaS1bar2S2bar2",
        "lambdaS1barS2S2bar2",
        "lambdaS1S1bar2S2bar",
        "lambdaHHdagS1barS2barCross",
    )

    complex_optional = {
        "lambdaHHdagS2S2",
        "lambdaHHdagS1barS1bar",
        "lambdaS1bar2S2bar2",
        "lambdaS1barS2S2bar2",
        "lambdaS1S1bar2S2bar",
        "lambdaHHdagS1barS2barCross",
    }

    for name in optional_names:
        if name not in ordinary:
            continue

        required[name] = (
            _complex_scalar(
                ordinary[name],
                f"base_state.ordinary.{name}",
            )
            if name in complex_optional
            else float(ordinary[name])
        )

    return T3UVState(**required).validated()


def build_scan_points(
    config: ResolvedScanConfig,
) -> tuple[dict[str, float], ...]:
    """Construct grid or random parameter points from config."""

    scan = config.raw["scan"]
    mode = str(scan["mode"]).lower()
    bindings = scan.get("bindings")

    if not isinstance(bindings, dict) or not bindings:
        raise ValueError(
            "scan.bindings must be a non-empty object."
        )

    if mode == "grid":
        values = scan.get("values")

        if not isinstance(values, dict):
            raise ValueError(
                "grid scan requires scan.values object."
            )

        if set(values) != set(bindings):
            raise ValueError(
                "grid scan.values keys must exactly match scan.bindings keys."
            )

        return grid_parameter_points(values)

    definitions = scan.get("parameters")

    if not isinstance(definitions, list) or not definitions:
        raise ValueError(
            "random scan requires non-empty scan.parameters list."
        )

    parameters: list[ScanParameter] = []

    for item in definitions:
        if not isinstance(item, dict):
            raise ValueError(
                "Each random scan.parameters entry must be an object."
            )

        name = str(item["name"])

        if name not in bindings:
            raise ValueError(
                f"Random scan parameter {name!r} has no scan.bindings entry."
            )

        parameters.append(
            ScanParameter(
                name=name,
                low=float(item["low"]),
                high=float(item["high"]),
                scale=str(item.get("scale", "linear")),
            )
        )

    if {parameter.name for parameter in parameters} != set(bindings):
        raise ValueError(
            "random scan.parameters names must exactly match scan.bindings keys."
        )

    return random_parameter_points(
        parameters,
        int(scan["n_points"]),
        seed=(
            None
            if scan.get("seed") is None
            else int(scan["seed"])
        ),
    )


def validate_scan_inputs(
    config: ResolvedScanConfig,
) -> tuple[dict[str, float], ...]:
    """Validate files, state schema, bindings and point generation."""

    for name, path in (
        ("UV RGBeta payload", config.uv_rgbeta_path),
        ("EFT1 RGBeta payload", config.eft1_rgbeta_path),
        ("final Weinberg coefficient", config.final_weinberg_path),
        ("oscillation target", config.oscillation_target_path),
    ):
        if not path.is_file():
            raise FileNotFoundError(
                f"{name} does not exist: {path}"
            )

    points = build_scan_points(config)

    build_uv_state_from_config(
        config,
        points[0],
    )

    representation = _build_representation(
        config.raw["representation"]
    )

    uv_payload = load_rgbeta_payload(
        config.uv_rgbeta_path
    )
    eft1_payload = load_rgbeta_payload(
        config.eft1_rgbeta_path
    )

    _validate_rgbeta_payload_compatibility(
        payload=uv_payload,
        representation=representation,
        payload_name="UV RGBeta payload",
        eft1=False,
    )
    _validate_rgbeta_payload_compatibility(
        payload=eft1_payload,
        representation=representation,
        payload_name="EFT1 RGBeta payload",
        eft1=True,
    )

    OscillationFitTarget.from_json(
        config.oscillation_target_path
    )

    final_payload = _read_json_object(
        config.final_weinberg_path
    )

    if final_payload.get("status") != "Success":
        raise ValueError(
            "final Weinberg JSON status is not Success."
        )

    combined = final_payload.get("combined", {})
    if not isinstance(combined, dict) or not combined.get(
        "ready_for_physical_majorana_numerics",
        False,
    ):
        raise ValueError(
            "final Weinberg JSON is not marked ready for physical Majorana numerics."
        )

    return points


def execute_scan(
    config: ResolvedScanConfig,
) -> ParameterScanResult:
    """Execute the configured end-to-end ordinary-T3 numerical scan."""

    points = validate_scan_inputs(config)

    uv_payload = load_rgbeta_payload(
        config.uv_rgbeta_path
    )
    eft1_payload = load_rgbeta_payload(
        config.eft1_rgbeta_path
    )
    target = OscillationFitTarget.from_json(
        config.oscillation_target_path
    )

    scales = config.raw["scales"]
    numerical = config.raw.get("numerical", {})
    if not isinstance(numerical, dict):
        raise ValueError("numerical must be an object when supplied.")

    c5_builder = FinalC5TrajectoryEvaluator(
        config.final_weinberg_path
    )

    evaluator = make_t3_point_evaluator(
        uv_state_builder=lambda parameters: build_uv_state_from_config(
            config,
            parameters,
        ),
        c5_builder=c5_builder,
        uv_rgbeta_payload=uv_payload,
        eft1_rgbeta_payload=eft1_payload,
        target=target,
        mu_fermion_threshold_gev=float(
            scales["mu_fermion_threshold_gev"]
        ),
        mu_scalar_threshold_gev=float(
            scales["mu_scalar_threshold_gev"]
        ),
        mu_low_gev=float(scales["mu_low_gev"]),
        vev_gev=float(
            numerical.get("vev_gev", 246.22)
        ),
        rtol=float(numerical.get("rtol", 1.0e-8)),
        atol=float(numerical.get("atol", 1.0e-11)),
    )

    scan = config.raw["scan"]

    def checkpoint_callback(partial: ParameterScanResult) -> None:
        write_scan_result(
            partial,
            config.checkpoint_path,
        )

    def progress_callback(
        completed: int,
        total: int,
        point,
    ) -> None:
        if point.succeeded:
            detail = f"chi2={point.chi2:.6g}"
        else:
            message = (point.error or "unknown error").replace("\n", " ")
            if len(message) > 140:
                message = message[:137] + "..."
            detail = f"FAILED: {message}"

        print(
            f"[{completed}/{total}] point={point.index} "
            f"{point.status} {detail}",
            flush=True,
        )

    result = run_parameter_scan(
        points,
        evaluator,
        accepted_chi2_max=(
            None
            if scan.get("accepted_chi2_max") is None
            else float(scan["accepted_chi2_max"])
        ),
        fail_fast=bool(scan.get("fail_fast", False)),
        progress_callback=progress_callback,
        checkpoint_callback=checkpoint_callback,
    )

    write_scan_result(
        result,
        config.output_path,
    )

    return result


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run the numerical ordinary-T3 parameter scan from one JSON "
            "configuration file."
        )
    )
    parser.add_argument(
        "config",
        type=Path,
        help="Path to scan JSON configuration.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Validate paths, state construction, parameter bindings and "
            "external JSON contracts without running the RGEs."
        ),
    )
    return parser


def main() -> int:
    args = _build_parser().parse_args()

    try:
        config = load_scan_config(args.config)
        points = validate_scan_inputs(config)

        if args.dry_run:
            print(
                json.dumps(
                    {
                        "status": "Validated",
                        "config": str(config.source_path),
                        "point_count": len(points),
                        "output": str(config.output_path),
                        "checkpoint": str(config.checkpoint_path),
                    },
                    indent=2,
                )
            )
            return 0

        result = execute_scan(config)
        best = result.best_point

        print(
            json.dumps(
                {
                    "status": (
                        "Success"
                        if result.successful_points
                        else "Failed"
                    ),
                    "point_count": len(result.points),
                    "successful_point_count": len(
                        result.successful_points
                    ),
                    "failed_point_count": len(
                        result.failed_points
                    ),
                    "accepted_point_count": len(
                        result.accepted_points
                    ),
                    "best_point_index": (
                        None if best is None else best.index
                    ),
                    "best_chi2": (
                        None if best is None else best.chi2
                    ),
                    "output": str(config.output_path),
                    "checkpoint": str(config.checkpoint_path),
                },
                indent=2,
            )
        )

        return 0 if result.successful_points else 1

    except Exception as exc:
        print(
            json.dumps(
                {
                    "status": "Failed",
                    "error": str(exc),
                },
                indent=2,
            )
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
