"""Parameter-scan driver for the numerical T3 pipeline.

This layer is deliberately split into two parts:

1. generic deterministic scan machinery;
2. a T3-aware evaluator that composes the numerical trajectory, an externally
   supplied matched C5 builder, low-energy Weinberg running, and the oscillation
   chi-square layer.

The matched C5 is *not* invented here. It must be supplied by ``c5_builder``,
which is the bridge to the existing symbolic/Matchete matching pipeline or to
a separately validated numerical matching implementation.

A scan point therefore follows

    parameters
      -> UV state builder
      -> UV running
      -> F threshold
      -> scalar-only intermediate running
      -> scalar threshold
      -> c5_builder(parameters, trajectory)
      -> SM + C5 running
      -> m_nu(mu_low)
      -> oscillation chi^2.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import json
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np

from Numerical.OscillationFit import (
    OscillationFitResult,
    OscillationFitTarget,
    evaluate_mass_matrix_fit,
)
from Numerical.T3Trajectory import (
    T3RenormalisableTrajectory,
    continue_with_weinberg_running,
    run_renormalisable_t3_trajectory,
)
from Numerical.State import SharedT3UVState, T3UVState
from physics.NeutrinoTrajectory import (
    charged_lepton_mass_basis_matrix,
    neutrino_mass_from_c5,
)


UVState = T3UVState | SharedT3UVState
ParameterPoint = dict[str, float]
PointEvaluator = Callable[[Mapping[str, float]], OscillationFitResult]
ProgressCallback = Callable[[int, int, "ScanPointResult"], None]
CheckpointCallback = Callable[["ParameterScanResult"], None]
UVStateBuilder = Callable[[Mapping[str, float]], UVState]
C5Builder = Callable[
    [Mapping[str, float], T3RenormalisableTrajectory],
    np.ndarray,
]


@dataclass(frozen=True)
class ScanParameter:
    """One random-scan parameter."""

    name: str
    low: float
    high: float
    scale: str = "linear"

    def validated(self) -> "ScanParameter":
        name = str(self.name).strip()
        low = float(self.low)
        high = float(self.high)
        scale = str(self.scale).lower()

        if not name:
            raise ValueError("Scan parameter name cannot be empty.")

        if not np.isfinite(low) or not np.isfinite(high):
            raise ValueError(
                f"Scan bounds for {name} must be finite."
            )

        if not low < high:
            raise ValueError(
                f"Scan parameter {name} requires low < high."
            )

        if scale not in {"linear", "log"}:
            raise ValueError(
                f"Scan parameter {name} scale must be 'linear' or 'log'."
            )

        if scale == "log" and low <= 0.0:
            raise ValueError(
                f"Log-scanned parameter {name} requires low > 0."
            )

        return ScanParameter(
            name=name,
            low=low,
            high=high,
            scale=scale,
        )


@dataclass(frozen=True)
class ScanPointResult:
    """Result of evaluating one parameter point."""

    index: int
    parameters: ParameterPoint
    status: str
    chi2: float | None
    error: str | None
    ordering: str | None
    prediction: dict[str, float] | None

    @property
    def succeeded(self) -> bool:
        return self.status == "Success"


@dataclass(frozen=True)
class ParameterScanResult:
    """Complete scan output."""

    points: tuple[ScanPointResult, ...]
    accepted_chi2_max: float | None = None

    @property
    def successful_points(self) -> tuple[ScanPointResult, ...]:
        return tuple(
            point
            for point in self.points
            if point.succeeded
        )

    @property
    def failed_points(self) -> tuple[ScanPointResult, ...]:
        return tuple(
            point
            for point in self.points
            if not point.succeeded
        )

    @property
    def best_point(self) -> ScanPointResult | None:
        successful = self.successful_points
        if not successful:
            return None
        return min(
            successful,
            key=lambda point: float(point.chi2),
        )

    @property
    def accepted_points(self) -> tuple[ScanPointResult, ...]:
        if self.accepted_chi2_max is None:
            return self.successful_points

        return tuple(
            point
            for point in self.successful_points
            if float(point.chi2) <= self.accepted_chi2_max
        )


def grid_parameter_points(
    values: Mapping[str, Sequence[float]],
) -> tuple[ParameterPoint, ...]:
    """Return the Cartesian product of explicit parameter values."""

    if not values:
        raise ValueError(
            "At least one grid parameter is required."
        )

    names = tuple(values)

    sequences: list[tuple[float, ...]] = []

    for name in names:
        raw = tuple(float(value) for value in values[name])

        if not raw:
            raise ValueError(
                f"Grid parameter {name} has no values."
            )

        if not np.all(np.isfinite(raw)):
            raise ValueError(
                f"Grid parameter {name} contains non-finite values."
            )

        sequences.append(raw)

    return tuple(
        {
            name: float(value)
            for name, value in zip(names, combination)
        }
        for combination in product(*sequences)
    )


def random_parameter_points(
    parameters: Sequence[ScanParameter],
    n_points: int,
    *,
    seed: int | None = None,
) -> tuple[ParameterPoint, ...]:
    """Generate reproducible uniform linear/log parameter points."""

    count = int(n_points)

    if count <= 0:
        raise ValueError("n_points must be positive.")

    validated = tuple(
        parameter.validated()
        for parameter in parameters
    )

    if not validated:
        raise ValueError(
            "At least one random-scan parameter is required."
        )

    names = [parameter.name for parameter in validated]

    if len(set(names)) != len(names):
        raise ValueError(
            "Random-scan parameter names must be unique."
        )

    rng = np.random.default_rng(seed)

    points: list[ParameterPoint] = []

    for _ in range(count):
        point: ParameterPoint = {}

        for parameter in validated:
            u = float(rng.random())

            if parameter.scale == "linear":
                value = parameter.low + u * (
                    parameter.high - parameter.low
                )
            else:
                log_low = np.log(parameter.low)
                log_high = np.log(parameter.high)
                value = float(
                    np.exp(
                        log_low
                        + u * (log_high - log_low)
                    )
                )

            point[parameter.name] = float(value)

        points.append(point)

    return tuple(points)


def _prediction_dict(
    fit_result: OscillationFitResult,
) -> dict[str, float]:
    return {
        name: float(value)
        for name, value in zip(
            fit_result.observable_names,
            fit_result.prediction,
        )
    }


def run_parameter_scan(
    points: Iterable[Mapping[str, float]],
    evaluator: PointEvaluator,
    *,
    accepted_chi2_max: float | None = None,
    fail_fast: bool = False,
    progress_callback: ProgressCallback | None = None,
    checkpoint_callback: CheckpointCallback | None = None,
) -> ParameterScanResult:
    """Evaluate points sequentially with optional live progress/checkpoint hooks."""

    cutoff: float | None

    if accepted_chi2_max is None:
        cutoff = None
    else:
        cutoff = float(accepted_chi2_max)

        if not np.isfinite(cutoff) or cutoff < 0.0:
            raise ValueError(
                "accepted_chi2_max must be finite and non-negative."
            )

    point_sequence = tuple(points)
    total_points = len(point_sequence)
    results: list[ScanPointResult] = []

    for index, raw_parameters in enumerate(point_sequence):
        parameters = {
            str(name): float(value)
            for name, value in raw_parameters.items()
        }

        if not parameters:
            raise ValueError(
                f"Scan point {index} contains no parameters."
            )

        if not np.all(
            np.isfinite(list(parameters.values()))
        ):
            raise ValueError(
                f"Scan point {index} contains non-finite parameters."
            )

        try:
            fit = evaluator(parameters)

            chi2 = float(fit.chi2)

            if not np.isfinite(chi2) or chi2 < 0.0:
                raise RuntimeError(
                    f"Evaluator returned invalid chi2={chi2}."
                )

            results.append(
                ScanPointResult(
                    index=index,
                    parameters=parameters,
                    status="Success",
                    chi2=chi2,
                    error=None,
                    ordering=fit.ordering,
                    prediction=_prediction_dict(fit),
                )
            )

        except Exception as exc:
            if fail_fast:
                raise

            results.append(
                ScanPointResult(
                    index=index,
                    parameters=parameters,
                    status="Failed",
                    chi2=None,
                    error=str(exc),
                    ordering=None,
                    prediction=None,
                )
            )

        partial_result = ParameterScanResult(
            points=tuple(results),
            accepted_chi2_max=cutoff,
        )

        if checkpoint_callback is not None:
            checkpoint_callback(partial_result)

        if progress_callback is not None:
            progress_callback(
                index + 1,
                total_points,
                results[-1],
            )

    return ParameterScanResult(
        points=tuple(results),
        accepted_chi2_max=cutoff,
    )


def make_t3_point_evaluator(
    *,
    uv_state_builder: UVStateBuilder,
    c5_builder: C5Builder,
    uv_rgbeta_payload: Mapping[str, Any],
    eft1_rgbeta_payload: Mapping[str, Any],
    target: OscillationFitTarget,
    mu_fermion_threshold_gev: float,
    mu_scalar_threshold_gev: float,
    mu_low_gev: float,
    vev_gev: float = 246.22,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-11,
) -> PointEvaluator:
    """Build a point evaluator for the complete numerical T3 trajectory.

    ``uv_state_builder(parameters)`` maps scan parameters to a validated UV
    numerical state.

    ``c5_builder(parameters, trajectory)`` must return the physical symmetric
    C5 matrix at the scalar threshold. It is intentionally external because
    the authoritative project matching remains in the existing Matchete/
    symbolic pipeline.
    """

    target = target.validated()
    mu_f = float(mu_fermion_threshold_gev)
    mu_s = float(mu_scalar_threshold_gev)
    mu_low = float(mu_low_gev)
    vev = float(vev_gev)

    def evaluator(
        parameters: Mapping[str, float],
    ) -> OscillationFitResult:
        uv_state = uv_state_builder(
            parameters
        ).validated()

        trajectory = run_renormalisable_t3_trajectory(
            uv_state,
            uv_rgbeta_payload,
            eft1_rgbeta_payload,
            mu_fermion_threshold_gev=mu_f,
            mu_scalar_threshold_gev=mu_s,
            rtol=rtol,
            atol=atol,
        )

        c5 = np.asarray(
            c5_builder(parameters, trajectory),
            dtype=complex,
        )

        full = continue_with_weinberg_running(
            trajectory,
            c5,
            mu_low,
            rtol=rtol,
            atol=atol,
        )

        mass_matrix = neutrino_mass_from_c5(
            full.weinberg.K,
            vev_gev=vev,
        )

        mass_matrix_flavor = charged_lepton_mass_basis_matrix(
            mass_matrix,
            full.weinberg.ye,
        )

        return evaluate_mass_matrix_fit(
            mass_matrix_flavor,
            target,
        )

    return evaluator


def scan_result_to_json_dict(
    result: ParameterScanResult,
) -> dict[str, Any]:
    """Return a stable JSON-safe scan payload."""

    best = result.best_point

    return {
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
        "accepted_chi2_max": result.accepted_chi2_max,
        "accepted_point_count": len(
            result.accepted_points
        ),
        "best_point_index": (
            best.index if best is not None else None
        ),
        "best_chi2": (
            best.chi2 if best is not None else None
        ),
        "points": [
            {
                "index": point.index,
                "parameters": point.parameters,
                "status": point.status,
                "chi2": point.chi2,
                "error": point.error,
                "ordering": point.ordering,
                "prediction": point.prediction,
            }
            for point in result.points
        ],
    }


def write_scan_result(
    result: ParameterScanResult,
    path: Path,
) -> Path:
    """Write the complete scan result as JSON."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(
        scan_result_to_json_dict(result),
        indent=2,
    )
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)
    return path
