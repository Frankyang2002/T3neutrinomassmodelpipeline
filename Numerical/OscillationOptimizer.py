"""
Efficient local oscillation-fit optimiser for the numerical T3 pipeline.

Strategy
--------
The expensive object in this project is one complete T3 point evaluation:

    UV state
      -> UV running
      -> F threshold
      -> EFT1 running
      -> scalar threshold
      -> matched C5
      -> SM + Weinberg running
      -> neutrino observables
      -> Gaussian residual vector.

A brute-force random scan spends most evaluations away from promising regions.
This module instead:

1. loads either an existing successful scan point or the best point from a
   previous optimiser run as a seed;
2. maps every bounded scan parameter to a dimensionless unit coordinate z in [0, 1];
3. estimates a local finite-difference sensitivity of the five whitened
   oscillation residuals to every scan parameter;
4. keeps only the most locally influential parameters (default: 5);
5. runs bounded nonlinear least squares in that reduced subspace.

The remaining parameters stay fixed at the seed values.  A later optimisation
can increase ``--active-count`` if the reduced fit plateaus.

The optimiser reuses the same authoritative numerical evaluator as ScanCLI.
It does not introduce a new matching formula or a new neutrino observable
definition.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import numpy as np
from scipy.optimize import least_squares


ResidualEvaluator = Callable[[Mapping[str, float]], Any]


@dataclass(frozen=True)
class OptimizerParameter:
    """One bounded optimisation parameter."""

    name: str
    low: float
    high: float
    scale: str = "linear"

    def validated(self) -> "OptimizerParameter":
        name = str(self.name).strip()
        low = float(self.low)
        high = float(self.high)
        scale = str(self.scale).lower()

        if not name:
            raise ValueError("Optimizer parameter name cannot be empty.")
        if not np.isfinite(low) or not np.isfinite(high):
            raise ValueError(f"Bounds for {name} must be finite.")
        if not low < high:
            raise ValueError(f"Parameter {name} requires low < high.")
        if scale not in {"linear", "log"}:
            raise ValueError(
                f"Parameter {name} scale must be 'linear' or 'log'."
            )
        if scale == "log" and low <= 0.0:
            raise ValueError(
                f"Log-scaled parameter {name} requires low > 0."
            )

        return OptimizerParameter(
            name=name,
            low=low,
            high=high,
            scale=scale,
        )

    def to_unit(self, value: float) -> float:
        """Map a physical value to z in [0, 1]."""

        parameter = self.validated()
        value = float(value)

        if parameter.scale == "linear":
            z = (value - parameter.low) / (
                parameter.high - parameter.low
            )
        else:
            z = (
                math.log(value) - math.log(parameter.low)
            ) / (
                math.log(parameter.high)
                - math.log(parameter.low)
            )

        return float(z)

    def from_unit(self, z: float) -> float:
        """Map z in [0, 1] back to the physical parameter."""

        parameter = self.validated()
        z = float(z)

        if parameter.scale == "linear":
            return float(
                parameter.low
                + z * (parameter.high - parameter.low)
            )

        log_value = (
            math.log(parameter.low)
            + z
            * (
                math.log(parameter.high)
                - math.log(parameter.low)
            )
        )
        return float(math.exp(log_value))


@dataclass(frozen=True)
class SensitivityResult:
    """Local residual sensitivity for one parameter."""

    name: str
    score: float
    derivative: tuple[float, ...]
    step_low: float
    step_high: float
    status: str
    error: str | None = None


class EvaluationRecorder:
    """Evaluate the expensive pipeline while retaining/checkpointing progress."""

    def __init__(
        self,
        *,
        evaluator: ResidualEvaluator,
        parameters: Sequence[OptimizerParameter],
        observable_names: Sequence[str],
        checkpoint_path: Path | None = None,
        failure_residual: float = 1.0e4,
    ) -> None:
        self.evaluator = evaluator
        self.parameters = tuple(
            parameter.validated()
            for parameter in parameters
        )
        self.observable_names = tuple(str(x) for x in observable_names)
        self.checkpoint_path = (
            None
            if checkpoint_path is None
            else Path(checkpoint_path)
        )
        self.failure_residual = float(failure_residual)

        if not np.isfinite(self.failure_residual):
            raise ValueError("failure_residual must be finite.")
        if self.failure_residual <= 0.0:
            raise ValueError("failure_residual must be positive.")

        self.history: list[dict[str, Any]] = []
        self.best: dict[str, Any] | None = None
        self.phase: str = "initialisation"

    def physical_from_unit(
        self,
        z: Sequence[float],
    ) -> dict[str, float]:
        values = np.asarray(z, dtype=float)
        if values.shape != (len(self.parameters),):
            raise ValueError(
                f"Expected unit vector shape ({len(self.parameters)},), "
                f"got {values.shape}."
            )

        return {
            parameter.name: parameter.from_unit(float(value))
            for parameter, value in zip(self.parameters, values)
        }

    def _write_checkpoint(self) -> None:
        if self.checkpoint_path is None:
            return

        payload = {
            "status": "Running",
            "phase": self.phase,
            "total_pipeline_evaluations": len(self.history),
            "best": self.best,
            "history": self.history,
        }

        path = self.checkpoint_path
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_text(
            json.dumps(payload, indent=2),
            encoding="utf-8",
        )
        temporary.replace(path)

    def evaluate_unit(
        self,
        z: Sequence[float],
    ) -> np.ndarray:
        """Return whitened residuals, converting failures into a large penalty."""

        parameters = self.physical_from_unit(z)
        record: dict[str, Any] = {
            "evaluation": len(self.history),
            "phase": self.phase,
            "parameters": parameters,
        }

        try:
            fit = self.evaluator(parameters)
            residual = np.asarray(
                fit.whitened_residual,
                dtype=float,
            )

            if residual.shape != (len(self.observable_names),):
                raise RuntimeError(
                    "Evaluator returned whitened residual with unexpected "
                    f"shape {residual.shape}; expected "
                    f"({len(self.observable_names)},)."
                )
            if not np.all(np.isfinite(residual)):
                raise RuntimeError(
                    "Evaluator returned non-finite whitened residuals."
                )

            chi2 = float(np.dot(residual, residual))
            prediction = {
                name: float(value)
                for name, value in zip(
                    fit.observable_names,
                    fit.prediction,
                )
            }

            record.update(
                {
                    "status": "Success",
                    "chi2": chi2,
                    "ordering": str(fit.ordering),
                    "prediction": prediction,
                    "whitened_residual": residual.tolist(),
                    "error": None,
                }
            )

            if (
                self.best is None
                or chi2 < float(self.best["chi2"])
            ):
                self.best = {
                    "chi2": chi2,
                    "parameters": parameters,
                    "ordering": str(fit.ordering),
                    "prediction": prediction,
                    "observable_names": list(fit.observable_names),
                    "whitened_residual": residual.tolist(),
                }

        except Exception as exc:
            residual = np.full(
                len(self.observable_names),
                self.failure_residual,
                dtype=float,
            )
            record.update(
                {
                    "status": "Failed",
                    "chi2": None,
                    "ordering": None,
                    "prediction": None,
                    "whitened_residual": None,
                    "error": str(exc),
                }
            )

        self.history.append(record)
        self._write_checkpoint()

        if record["status"] == "Success":
            print(
                f"[eval {record['evaluation']}] "
                f"{self.phase} chi2={record['chi2']:.8g}",
                flush=True,
            )
        else:
            message = str(record["error"]).replace("\n", " ")
            if len(message) > 140:
                message = message[:137] + "..."
            print(
                f"[eval {record['evaluation']}] "
                f"{self.phase} FAILED: {message}",
                flush=True,
            )

        return residual


def parameters_from_scan_config(
    raw_config: Mapping[str, Any],
) -> tuple[OptimizerParameter, ...]:
    """Read bounded random-scan parameters in the exact config order."""

    scan = raw_config.get("scan")
    if not isinstance(scan, Mapping):
        raise ValueError("scan config requires a 'scan' object.")

    mode = str(scan.get("mode", "")).lower()
    if mode != "random":
        raise ValueError(
            "OscillationOptimizer currently requires a random-scan config "
            "with explicit bounds."
        )

    definitions = scan.get("parameters")
    bindings = scan.get("bindings")

    if not isinstance(definitions, list) or not definitions:
        raise ValueError(
            "scan.parameters must be a non-empty list."
        )
    if not isinstance(bindings, Mapping) or not bindings:
        raise ValueError(
            "scan.bindings must be a non-empty object."
        )

    parameters = tuple(
        OptimizerParameter(
            name=str(item["name"]),
            low=float(item["low"]),
            high=float(item["high"]),
            scale=str(item.get("scale", "linear")),
        ).validated()
        for item in definitions
    )

    names = tuple(parameter.name for parameter in parameters)

    if len(set(names)) != len(names):
        raise ValueError("Optimizer parameter names must be unique.")
    if set(names) != set(str(name) for name in bindings):
        raise ValueError(
            "scan.parameters names must exactly match scan.bindings keys."
        )

    return parameters


def unit_vector_from_parameters(
    parameters: Sequence[OptimizerParameter],
    values: Mapping[str, float],
    *,
    tolerance: float = 1.0e-10,
) -> np.ndarray:
    """Convert a complete physical parameter mapping to unit coordinates."""

    expected = {parameter.name for parameter in parameters}
    supplied = set(values)

    missing = expected - supplied
    if missing:
        raise ValueError(
            "Seed point is missing parameter(s): "
            + ", ".join(sorted(missing))
        )

    z = np.asarray(
        [
            parameter.to_unit(float(values[parameter.name]))
            for parameter in parameters
        ],
        dtype=float,
    )

    if np.any(z < -tolerance) or np.any(z > 1.0 + tolerance):
        offending = [
            parameter.name
            for parameter, value in zip(parameters, z)
            if value < -tolerance or value > 1.0 + tolerance
        ]
        raise ValueError(
            "Seed lies outside optimisation bounds for: "
            + ", ".join(offending)
        )

    return np.clip(z, 0.0, 1.0)


def load_seed_point(
    scan_result_path: Path,
    *,
    seed_index: int | None = None,
) -> dict[str, Any]:
    """Load one successful seed point from an existing scan result."""

    path = Path(scan_result_path)
    payload = json.loads(path.read_text(encoding="utf-8"))

    if not isinstance(payload, dict):
        raise ValueError("Seed result JSON must contain an object.")

    if seed_index is None:
        raw_index = payload.get("best_point_index")
        if raw_index is None:
            raise ValueError(
                "Seed result has no best_point_index; supply --seed-index."
            )
        seed_index = int(raw_index)

    points = payload.get("points")
    if not isinstance(points, list):
        raise ValueError("Seed result JSON must contain a points list.")

    for point in points:
        if (
            isinstance(point, dict)
            and int(point.get("index", -1)) == int(seed_index)
        ):
            if point.get("status") != "Success":
                raise ValueError(
                    f"Seed point {seed_index} is not successful."
                )
            if not isinstance(point.get("parameters"), dict):
                raise ValueError(
                    f"Seed point {seed_index} has no parameter mapping."
                )
            return {
                "source_kind": "scan_result",
                "source_path": str(path),
                "index": int(seed_index),
                "chi2": (
                    None
                    if point.get("chi2") is None
                    else float(point["chi2"])
                ),
                "parameters": dict(point["parameters"]),
                "ordering": point.get("ordering"),
                "prediction": point.get("prediction"),
            }

    raise ValueError(
        f"Seed point index {seed_index} was not found."
    )


def load_optimizer_best(
    optimizer_result_path: Path,
) -> dict[str, Any]:
    """Load the ``best`` point from a previous OscillationOptimizer result."""

    path = Path(optimizer_result_path)
    payload = json.loads(path.read_text(encoding="utf-8"))

    if not isinstance(payload, dict):
        raise ValueError("Optimizer result JSON must contain an object.")

    if payload.get("status") != "Success":
        raise ValueError(
            "Previous optimizer result does not have status='Success'."
        )

    best = payload.get("best")
    if not isinstance(best, dict):
        raise ValueError(
            "Previous optimizer result does not contain a 'best' object."
        )

    parameters = best.get("parameters")
    if not isinstance(parameters, dict):
        raise ValueError(
            "Previous optimizer best point has no parameter mapping."
        )

    chi2 = best.get("chi2")
    if chi2 is None:
        raise ValueError(
            "Previous optimizer best point has no chi2 value."
        )

    return {
        "source_kind": "optimizer_result",
        "source_path": str(path),
        "index": None,
        "chi2": float(chi2),
        "parameters": dict(parameters),
        "ordering": best.get("ordering"),
        "prediction": best.get("prediction"),
    }


def estimate_local_sensitivities(
    recorder: EvaluationRecorder,
    seed_z: Sequence[float],
    *,
    step: float = 0.02,
) -> tuple[np.ndarray, tuple[SensitivityResult, ...]]:
    """Finite-difference residual sensitivity in dimensionless unit coordinates."""

    step = float(step)
    if not (0.0 < step <= 0.25):
        raise ValueError(
            "sensitivity step must satisfy 0 < step <= 0.25."
        )

    seed = np.asarray(seed_z, dtype=float)
    if seed.shape != (len(recorder.parameters),):
        raise ValueError("seed_z has the wrong shape.")

    recorder.phase = "sensitivity_seed"
    baseline = recorder.evaluate_unit(seed)

    results: list[SensitivityResult] = []

    for index, parameter in enumerate(recorder.parameters):
        low_z = seed.copy()
        high_z = seed.copy()
        low_z[index] = max(0.0, seed[index] - step)
        high_z[index] = min(1.0, seed[index] + step)

        denominator = float(
            high_z[index] - low_z[index]
        )

        if denominator <= 0.0:
            results.append(
                SensitivityResult(
                    name=parameter.name,
                    score=0.0,
                    derivative=tuple(
                        0.0 for _ in recorder.observable_names
                    ),
                    step_low=float(low_z[index]),
                    step_high=float(high_z[index]),
                    status="Degenerate",
                    error="No finite-difference interval is available.",
                )
            )
            continue

        try:
            if low_z[index] == seed[index]:
                residual_low = baseline
            else:
                recorder.phase = f"sensitivity:{parameter.name}:low"
                residual_low = recorder.evaluate_unit(low_z)

            if high_z[index] == seed[index]:
                residual_high = baseline
            else:
                recorder.phase = f"sensitivity:{parameter.name}:high"
                residual_high = recorder.evaluate_unit(high_z)

            derivative = (
                residual_high - residual_low
            ) / denominator
            score = float(np.linalg.norm(derivative))

            results.append(
                SensitivityResult(
                    name=parameter.name,
                    score=score,
                    derivative=tuple(float(x) for x in derivative),
                    step_low=float(low_z[index]),
                    step_high=float(high_z[index]),
                    status="Success",
                    error=None,
                )
            )

        except Exception as exc:
            results.append(
                SensitivityResult(
                    name=parameter.name,
                    score=0.0,
                    derivative=tuple(
                        0.0 for _ in recorder.observable_names
                    ),
                    step_low=float(low_z[index]),
                    step_high=float(high_z[index]),
                    status="Failed",
                    error=str(exc),
                )
            )

    ranked = tuple(
        sorted(
            results,
            key=lambda result: result.score,
            reverse=True,
        )
    )

    return baseline, ranked


def optimise_reduced_least_squares(
    recorder: EvaluationRecorder,
    seed_z: Sequence[float],
    sensitivity: Sequence[SensitivityResult],
    *,
    active_count: int = 5,
    max_nfev: int = 50,
    xtol: float = 1.0e-6,
    ftol: float = 1.0e-6,
    gtol: float = 1.0e-6,
) -> Any:
    """Optimise the most sensitive parameters while holding the rest fixed."""

    count = int(active_count)
    if count <= 0:
        raise ValueError("active_count must be positive.")
    if count > len(recorder.parameters):
        raise ValueError(
            "active_count cannot exceed the number of scan parameters."
        )

    ranked_names = [
        result.name
        for result in sensitivity
        if result.status == "Success"
    ]

    if len(ranked_names) < count:
        raise RuntimeError(
            f"Only {len(ranked_names)} parameters have usable sensitivity "
            f"estimates, fewer than active_count={count}."
        )

    active_names = tuple(ranked_names[:count])
    name_to_index = {
        parameter.name: index
        for index, parameter in enumerate(recorder.parameters)
    }
    active_indices = np.asarray(
        [name_to_index[name] for name in active_names],
        dtype=int,
    )

    seed = np.asarray(seed_z, dtype=float)
    active_seed = seed[active_indices].copy()

    def reduced_residual(active_z: np.ndarray) -> np.ndarray:
        full_z = seed.copy()
        full_z[active_indices] = active_z
        recorder.phase = "least_squares"
        return recorder.evaluate_unit(full_z)

    solution = least_squares(
        reduced_residual,
        active_seed,
        bounds=(
            np.zeros(count, dtype=float),
            np.ones(count, dtype=float),
        ),
        method="trf",
        jac="2-point",
        max_nfev=int(max_nfev),
        xtol=float(xtol),
        ftol=float(ftol),
        gtol=float(gtol),
        verbose=0,
    )

    return solution, active_names


def _build_project_evaluator(
    scan_config_path: Path,
) -> tuple[Any, Any]:
    """Construct exactly the same physical point evaluator used by ScanCLI."""

    from Numerical.FinalC5Bridge import FinalC5TrajectoryBuilder
    from Numerical.OscillationFit import OscillationFitTarget
    from Numerical.ParameterScan import make_t3_point_evaluator
    from Numerical.RGBetaEvaluator import load_rgbeta_payload
    from Numerical.ScanCLI import (
        build_uv_state_from_config,
        load_scan_config,
        validate_scan_inputs,
    )

    config = load_scan_config(scan_config_path)

    # This validates external files, metadata compatibility, target, final C5
    # readiness and at least one fully bound UV state.  Point generation is
    # cheap; it does not execute RG running.
    validate_scan_inputs(config)

    uv_payload = load_rgbeta_payload(config.uv_rgbeta_path)
    eft1_payload = load_rgbeta_payload(config.eft1_rgbeta_path)
    target = OscillationFitTarget.from_json(
        config.oscillation_target_path
    )

    scales = config.raw["scales"]
    numerical = config.raw.get("numerical", {})
    if not isinstance(numerical, dict):
        raise ValueError("numerical must be an object when supplied.")

    c5_builder = FinalC5TrajectoryBuilder(
        config.final_weinberg_path
    )

    evaluator = make_t3_point_evaluator(
        uv_state_builder=lambda values: build_uv_state_from_config(
            config,
            values,
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
        vev_gev=float(numerical.get("vev_gev", 246.22)),
        rtol=float(numerical.get("rtol", 1.0e-8)),
        atol=float(numerical.get("atol", 1.0e-11)),
    )

    return config, evaluator


def _write_final_output(
    *,
    output_path: Path,
    config: Any,
    seed_point: Mapping[str, Any],
    recorder: EvaluationRecorder,
    baseline_residual: np.ndarray,
    sensitivity: Sequence[SensitivityResult],
    solution: Any,
    active_names: Sequence[str],
) -> Path:
    output_path = Path(output_path)

    successful = [
        entry
        for entry in recorder.history
        if entry["status"] == "Success"
    ]
    failed = [
        entry
        for entry in recorder.history
        if entry["status"] != "Success"
    ]

    baseline_chi2 = float(
        np.dot(baseline_residual, baseline_residual)
    )

    payload = {
        "status": (
            "Success"
            if recorder.best is not None
            else "Failed"
        ),
        "scan_config": str(config.source_path),
        "seed_source_kind": str(seed_point.get("source_kind", "unknown")),
        "seed_source_path": seed_point.get("source_path"),
        "seed_point_index": (
            None
            if seed_point.get("index") is None
            else int(seed_point["index"])
        ),
        "seed_stored_chi2": (
            None
            if seed_point.get("chi2") is None
            else float(seed_point["chi2"])
        ),
        "seed_recomputed_chi2": baseline_chi2,
        "active_parameter_count": len(active_names),
        "active_parameters": list(active_names),
        "total_pipeline_evaluations": len(recorder.history),
        "successful_evaluations": len(successful),
        "failed_evaluations": len(failed),
        "optimizer": {
            "method": "scipy.optimize.least_squares",
            "algorithm": "trf",
            "jacobian": "2-point",
            "success": bool(solution.success),
            "status": int(solution.status),
            "message": str(solution.message),
            "nfev": int(solution.nfev),
            "njev": (
                None
                if solution.njev is None
                else int(solution.njev)
            ),
            "cost": float(solution.cost),
            "optimality": float(solution.optimality),
        },
        "sensitivity_ranking": [
            {
                "name": result.name,
                "score": result.score,
                "derivative": list(result.derivative),
                "step_low": result.step_low,
                "step_high": result.step_high,
                "status": result.status,
                "error": result.error,
            }
            for result in sensitivity
        ],
        "best": recorder.best,
        "history": recorder.history,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(
        output_path.name + ".tmp"
    )
    temporary.write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )
    temporary.replace(output_path)

    return output_path


def run_optimizer(
    *,
    scan_config_path: Path,
    seed_result_path: Path | None,
    seed_optimizer_result_path: Path | None,
    output_path: Path,
    seed_index: int | None = None,
    checkpoint_path: Path | None = None,
    active_count: int = 5,
    sensitivity_step: float = 0.02,
    max_nfev: int = 50,
) -> Path:
    """Execute sensitivity screening plus reduced bounded least squares."""

    config, evaluator = _build_project_evaluator(
        Path(scan_config_path)
    )
    parameters = parameters_from_scan_config(
        config.raw
    )
    if (seed_result_path is None) == (seed_optimizer_result_path is None):
        raise ValueError(
            "Supply exactly one of seed_result_path or "
            "seed_optimizer_result_path."
        )

    if seed_optimizer_result_path is not None:
        if seed_index is not None:
            raise ValueError(
                "seed_index is only valid with a scan-result seed."
            )
        seed_point = load_optimizer_best(
            Path(seed_optimizer_result_path)
        )
    else:
        seed_point = load_seed_point(
            Path(seed_result_path),
            seed_index=seed_index,
        )

    seed_z = unit_vector_from_parameters(
        parameters,
        seed_point["parameters"],
    )

    # Build the target once here only to obtain the canonical residual names.
    from Numerical.OscillationFit import OscillationFitTarget

    target = OscillationFitTarget.from_json(
        config.oscillation_target_path
    )

    recorder = EvaluationRecorder(
        evaluator=evaluator,
        parameters=parameters,
        observable_names=target.observable_names,
        checkpoint_path=checkpoint_path,
    )

    baseline, sensitivity = estimate_local_sensitivities(
        recorder,
        seed_z,
        step=sensitivity_step,
    )

    print("\nLocal sensitivity ranking:", flush=True)
    for rank, result in enumerate(sensitivity, start=1):
        print(
            f"  {rank:2d}. {result.name:20s} "
            f"score={result.score:.6g}",
            flush=True,
        )

    solution, active_names = optimise_reduced_least_squares(
        recorder,
        seed_z,
        sensitivity,
        active_count=active_count,
        max_nfev=max_nfev,
    )

    output = _write_final_output(
        output_path=Path(output_path),
        config=config,
        seed_point=seed_point,
        recorder=recorder,
        baseline_residual=baseline,
        sensitivity=sensitivity,
        solution=solution,
        active_names=active_names,
    )

    print(
        f"\nOptimisation complete. Best chi2="
        f"{None if recorder.best is None else recorder.best['chi2']}",
        flush=True,
    )
    print(f"Result written to: {output}", flush=True)

    return output


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Sensitivity-screened bounded least-squares optimisation "
            "of a T3 oscillation fit."
        )
    )
    parser.add_argument(
        "scan_config",
        type=Path,
        help=(
            "Existing random-scan config. Its parameter bounds define the "
            "optimisation box."
        ),
    )
    seed_group = parser.add_mutually_exclusive_group(required=True)
    seed_group.add_argument(
        "--seed-result",
        type=Path,
        default=None,
        help="Existing scan result JSON containing the seed point.",
    )
    seed_group.add_argument(
        "--seed-optimizer-result",
        type=Path,
        default=None,
        help=(
            "Previous OscillationOptimizer result JSON. The stored 'best' "
            "point is used as the new warm start."
        ),
    )
    parser.add_argument(
        "--seed-index",
        type=int,
        default=None,
        help=(
            "Successful point index for --seed-result. Default: seed "
            "result's best_point_index. Not valid with "
            "--seed-optimizer-result."
        ),
    )
    parser.add_argument(
        "--active-count",
        type=int,
        default=5,
        help=(
            "Number of most sensitive parameters to optimise. "
            "Default: 5."
        ),
    )
    parser.add_argument(
        "--sensitivity-step",
        type=float,
        default=0.02,
        help=(
            "Finite-difference step in the dimensionless [0,1] parameter "
            "box. Default: 0.02."
        ),
    )
    parser.add_argument(
        "--max-nfev",
        type=int,
        default=50,
        help=(
            "SciPy least_squares max_nfev after sensitivity screening. "
            "Default: 50."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Final optimisation JSON.",
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=None,
        help=(
            "Optional atomic checkpoint JSON updated after every expensive "
            "pipeline evaluation."
        ),
    )
    return parser


def main() -> int:
    args = _build_parser().parse_args()

    checkpoint = args.checkpoint
    if checkpoint is None:
        output = Path(args.output)
        checkpoint = output.with_name(
            output.stem + ".checkpoint.json"
        )

    run_optimizer(
        scan_config_path=args.scan_config,
        seed_result_path=args.seed_result,
        seed_optimizer_result_path=args.seed_optimizer_result,
        output_path=args.output,
        seed_index=args.seed_index,
        checkpoint_path=checkpoint,
        active_count=args.active_count,
        sensitivity_step=args.sensitivity_step,
        max_nfev=args.max_nfev,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
