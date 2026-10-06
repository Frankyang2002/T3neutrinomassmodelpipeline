"""Automatic Sobol benchmark search for ordinary T3 numerical runs.

This module chooses a reproducible UV benchmark when a model does not already
have a compatible cached benchmark.

The search is intentionally simple:

    bounded Sobol sample
        -> full authoritative numerical T3 evaluator
        -> lowest oscillation chi-square
        -> cache and reuse

It does not introduce a new matching formula, RGE, or neutrino observable.
The expensive point evaluator is exactly the same implementation used by the
existing parameter-scan machinery.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Mapping, Sequence

import numpy as np
from scipy.stats import qmc

from common.RunRecords import RunRecord
from Numerical.running.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.fitting.OscillationFit import OscillationFitTarget
from Numerical.fitting.OscillationOptimizer import (
    EvaluationRecorder,
    OptimizerParameter,
    estimate_local_sensitivities,
    optimise_reduced_least_squares,
    unit_vector_from_parameters,
)
from Numerical.fitting.ParameterScan import make_t3_point_evaluator
from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
from Numerical.fitting.ScanCLI import build_uv_state_from_config


PROJECT_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_SCHEMA = "t3_sobol_refined_benchmark_v2"
LEGACY_SOBOL_SCHEMA = "t3_sobol_benchmark_v1"


@dataclass(frozen=True)
class BenchmarkParameter:
    """One bounded real parameter in the automatic benchmark search."""

    name: str
    binding: str
    low: float
    high: float
    scale: str = "linear"

    def validated(self) -> "BenchmarkParameter":
        name = str(self.name).strip()
        binding = str(self.binding).strip()
        low = float(self.low)
        high = float(self.high)
        scale = str(self.scale).lower()

        if not name:
            raise ValueError("Benchmark parameter name cannot be empty.")
        if not binding:
            raise ValueError(f"Benchmark parameter {name} requires a binding.")
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

        return BenchmarkParameter(
            name=name,
            binding=binding,
            low=low,
            high=high,
            scale=scale,
        )

    def from_unit(self, value: float) -> float:
        """Map one Sobol coordinate in [0,1) to the physical parameter."""

        parameter = self.validated()
        z = float(value)

        if not 0.0 <= z < 1.0:
            raise ValueError(
                f"Sobol coordinate for {parameter.name} must satisfy 0 <= z < 1."
            )

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
class BenchmarkResolution:
    """Resolved benchmark payload plus provenance for the numerical pipeline."""

    payload: dict[str, Any]
    metadata: dict[str, Any]
    cache_path: Path | None


def default_ordinary_t3_parameters() -> tuple[BenchmarkParameter, ...]:
    """Return the canonical first-pass real-flavor benchmark search profile."""

    parameters: list[BenchmarkParameter] = [
        BenchmarkParameter(
            name="lambdaT3_real",
            binding="ordinary.lambdaT3.real",
            low=1.0e-3,
            high=3.0e-1,
            scale="log",
        )
    ]

    for yukawa in ("y1", "y2"):
        for i in range(3):
            for j in range(3):
                parameters.append(
                    BenchmarkParameter(
                        name=f"{yukawa}_{i+1}{j+1}",
                        binding=f"ordinary.{yukawa}.real.{i}.{j}",
                        low=-5.0e-1,
                        high=5.0e-1,
                        scale="linear",
                    )
                )

    return tuple(parameter.validated() for parameter in parameters)


def _parameters_from_config(
    payload: Mapping[str, Any],
) -> tuple[BenchmarkParameter, ...]:
    search = payload.get("benchmark_search", {})
    if not isinstance(search, Mapping):
        raise ValueError("benchmark_search must be an object.")

    raw_parameters = search.get("parameters")
    if raw_parameters is None:
        return default_ordinary_t3_parameters()

    if not isinstance(raw_parameters, list) or not raw_parameters:
        raise ValueError(
            "benchmark_search.parameters must be a non-empty list when supplied."
        )

    return tuple(
        BenchmarkParameter(
            name=str(item["name"]),
            binding=str(item["binding"]),
            low=float(item["low"]),
            high=float(item["high"]),
            scale=str(item.get("scale", "linear")),
        ).validated()
        for item in raw_parameters
    )


def _is_power_of_two(value: int) -> bool:
    return value > 0 and value & (value - 1) == 0


def generate_sobol_parameter_points(
    parameters: Sequence[BenchmarkParameter],
    n_points: int,
    *,
    seed: int,
) -> tuple[dict[str, float], ...]:
    """Generate a scrambled reproducible Sobol sample in physical coordinates."""

    validated = tuple(parameter.validated() for parameter in parameters)
    count = int(n_points)

    if not validated:
        raise ValueError("At least one Sobol parameter is required.")
    if not _is_power_of_two(count):
        raise ValueError(
            "Sobol benchmark n_points must be a positive power of two "
            "(for example 64, 128, 256 or 512)."
        )

    sampler = qmc.Sobol(
        d=len(validated),
        scramble=True,
        seed=int(seed),
    )
    unit_points = sampler.random_base2(
        m=int(math.log2(count))
    )

    return tuple(
        {
            parameter.name: parameter.from_unit(float(z))
            for parameter, z in zip(validated, row)
        }
        for row in unit_points
    )


def _representation_payload(record: RunRecord) -> dict[str, Any]:
    return {
        "d_s1": int(record.d_s1),
        "d_s2": int(record.d_s2),
        "d_f": int(record.d_f),
        "alpha": int(record.alpha),
        "shared_scalar": bool(record.shared_scalar),
    }


def model_benchmark_key(record: RunRecord) -> str:
    """Return a stable filesystem-safe identity for one T3 representation."""

    if record.shared_scalar:
        return (
            f"T3_shared_dS_{record.d_s1}_dF_{record.d_f}"
            f"_alpha_{record.alpha:+d}"
        ).replace("+", "p").replace("-", "m")

    return (
        f"T3_dS1_{record.d_s1}_dS2_{record.d_s2}_dF_{record.d_f}"
        f"_alpha_{record.alpha:+d}"
    ).replace("+", "p").replace("-", "m")


def retarget_payload_to_record(
    record: RunRecord,
    payload: Mapping[str, Any],
) -> dict[str, Any]:
    """Retarget a numerical template to the current model when explicitly enabled."""

    result = deepcopy(dict(payload))
    search = result.get("benchmark_search", {})

    if (
        isinstance(search, Mapping)
        and bool(search.get("enabled", False))
        and bool(search.get("use_current_model_representation", False))
    ):
        result["representation"] = _representation_payload(record)

    return result


def _set_nested_value(
    root: dict[str, Any],
    binding: str,
    value: float,
) -> None:
    """Set a dotted dictionary/list path relative to base_state."""

    parts = [part for part in str(binding).split(".") if part]
    if not parts:
        raise ValueError("Empty benchmark parameter binding.")

    current: Any = root

    for part in parts[:-1]:
        if isinstance(current, list):
            current = current[int(part)]
        elif isinstance(current, dict):
            current = current[part]
        else:
            raise TypeError(
                f"Cannot traverse benchmark binding {binding!r} through "
                f"{type(current).__name__}."
            )

    final = parts[-1]
    if isinstance(current, list):
        current[int(final)] = float(value)
    elif isinstance(current, dict):
        current[final] = float(value)
    else:
        raise TypeError(
            f"Cannot assign benchmark binding {binding!r} through "
            f"{type(current).__name__}."
        )


def apply_benchmark_parameters(
    payload: Mapping[str, Any],
    parameters: Mapping[str, float],
    definitions: Sequence[BenchmarkParameter],
) -> dict[str, Any]:
    """Return a copy of a numerical payload with benchmark values applied."""

    result = deepcopy(dict(payload))
    base_state = result.get("base_state")

    if not isinstance(base_state, dict):
        raise ValueError("Numerical payload requires a base_state object.")

    by_name = {
        parameter.name: parameter.validated()
        for parameter in definitions
    }

    missing = set(parameters) - set(by_name)
    if missing:
        raise KeyError(
            "Unknown benchmark parameter(s): "
            + ", ".join(sorted(missing))
        )

    for name, value in parameters.items():
        _set_nested_value(
            base_state,
            by_name[name].binding,
            float(value),
        )

    return result


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_sha256(value: Any) -> str:
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _effective_seed(base_seed: int, model_key: str) -> int:
    model_hash = hashlib.sha256(model_key.encode("utf-8")).digest()
    model_bits = int.from_bytes(model_hash[:4], "big")
    return int((int(base_seed) ^ model_bits) % (2**32))


def _resolve_target_path(
    payload: Mapping[str, Any],
) -> Path:
    search = payload.get("benchmark_search", {})
    if not isinstance(search, Mapping):
        raise ValueError("benchmark_search must be an object.")

    raw = search.get(
        "oscillation_target",
        "configs/nufit_6_1_ic24_no.json",
    )
    path = Path(str(raw))
    if not path.is_absolute():
        path = PROJECT_ROOT / path

    if not path.is_file():
        raise FileNotFoundError(
            f"Benchmark oscillation target does not exist: {path}"
        )

    return path.resolve()


def _compatibility_signature(
    *,
    record: RunRecord,
    payload: Mapping[str, Any],
    parameters: Sequence[BenchmarkParameter],
    target_path: Path,
    uv_rgbeta_path: Path,
    intermediate_rgbeta_path: Path,
    final_weinberg_path: Path,
    n_points: int,
    target_chi2: float,
    effective_seed: int,
) -> dict[str, Any]:
    scales = payload.get("scales")
    if not isinstance(scales, Mapping):
        raise ValueError("Numerical payload requires a scales object.")

    base_state = payload.get("base_state")
    if not isinstance(base_state, Mapping):
        raise ValueError("Numerical payload requires a base_state object.")

    return {
        "model": _representation_payload(record),
        "scales": {
            "mu_fermion_threshold_gev": float(
                scales["mu_fermion_threshold_gev"]
            ),
            "mu_scalar_threshold_gev": float(
                scales["mu_scalar_threshold_gev"]
            ),
            "mu_low_gev": float(scales["mu_low_gev"]),
        },
        "ordering": str(payload.get("ordering", "NO")).upper(),
        "parameter_profile": [
            {
                "name": parameter.name,
                "binding": parameter.binding,
                "low": parameter.low,
                "high": parameter.high,
                "scale": parameter.scale,
            }
            for parameter in parameters
        ],
        "n_points": int(n_points),
        "target_chi2": float(target_chi2),
        "effective_seed": int(effective_seed),
        "base_state_sha256": _json_sha256(base_state),
        "inputs_sha256": {
            "oscillation_target": _file_sha256(target_path),
            "uv_rgbeta": _file_sha256(uv_rgbeta_path),
            "intermediate_rgbeta": _file_sha256(intermediate_rgbeta_path),
            "final_weinberg": _file_sha256(final_weinberg_path),
        },
    }


def _cache_path(
    record: RunRecord,
) -> Path:
    return (
        PROJECT_ROOT
        / "output"
        / "benchmarks"
        / model_benchmark_key(record)
        / "sobol_benchmark.json"
    )


def _read_compatible_cache(
    path: Path,
    signature: Mapping[str, Any],
) -> dict[str, Any] | None:
    if not path.is_file():
        return None

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None

    if not isinstance(payload, dict):
        return None
    if payload.get("kind") != BENCHMARK_SCHEMA:
        return None
    if payload.get("compatibility") != dict(signature):
        return None

    best = payload.get("best")
    if not isinstance(best, dict):
        return None
    if not isinstance(best.get("parameters"), dict):
        return None
    if best.get("chi2") is None:
        return None

    return payload


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )
    temporary.replace(path)


def _build_evaluator(
    *,
    payload: Mapping[str, Any],
    parameters: Sequence[BenchmarkParameter],
    uv_rgbeta_path: Path,
    intermediate_rgbeta_path: Path,
    final_weinberg_path: Path,
    target_path: Path,
) -> Callable[[Mapping[str, float]], Any]:
    scales = payload["scales"]
    numerical = payload.get("numerical", {})

    if not isinstance(scales, Mapping):
        raise ValueError("Numerical payload requires a scales object.")
    if not isinstance(numerical, Mapping):
        raise ValueError("numerical must be an object.")

    bindings = {
        parameter.name: parameter.binding
        for parameter in parameters
    }

    config_raw = deepcopy(dict(payload))
    config_raw["scan"] = {
        "bindings": bindings,
    }
    config = SimpleNamespace(raw=config_raw)

    # Validate the zero-cost state-construction boundary before expensive work.
    first_point = {
        parameter.name: (
            math.sqrt(parameter.low * parameter.high)
            if parameter.scale == "log"
            else 0.5 * (parameter.low + parameter.high)
        )
        for parameter in parameters
    }
    build_uv_state_from_config(config, first_point)

    uv_payload = load_rgbeta_payload(uv_rgbeta_path)
    intermediate_payload = load_rgbeta_payload(intermediate_rgbeta_path)
    target = OscillationFitTarget.from_json(target_path)
    c5_builder = FinalC5TrajectoryEvaluator(final_weinberg_path)

    return make_t3_point_evaluator(
        uv_state_builder=lambda values: build_uv_state_from_config(
            config,
            values,
        ),
        c5_builder=c5_builder,
        uv_rgbeta_payload=uv_payload,
        eft1_rgbeta_payload=intermediate_payload,
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


def _prediction_dict(fit: Any) -> dict[str, float]:
    return {
        str(name): float(value)
        for name, value in zip(
            fit.observable_names,
            fit.prediction,
        )
    }


def _run_search(
    *,
    record: RunRecord,
    payload: Mapping[str, Any],
    definitions: Sequence[BenchmarkParameter],
    evaluator: Callable[[Mapping[str, float]], Any],
    cache_path: Path,
    signature: Mapping[str, Any],
    n_points: int,
    target_chi2: float,
    keep_best: int,
    effective_seed: int,
) -> dict[str, Any]:
    points = generate_sobol_parameter_points(
        definitions,
        n_points,
        seed=effective_seed,
    )

    successful: list[dict[str, Any]] = []
    failed_count = 0
    stopped_early = False

    print(
        f"  {record.name}: no compatible optimized benchmark found; "
        f"starting Sobol search ({n_points} points)...",
        flush=True,
    )

    for index, values in enumerate(points):
        try:
            fit = evaluator(values)
            chi2 = float(fit.chi2)

            if not np.isfinite(chi2):
                raise RuntimeError("Fit returned non-finite chi2.")

            entry = {
                "index": int(index),
                "chi2": chi2,
                "parameters": {
                    name: float(value)
                    for name, value in values.items()
                },
                "ordering": str(fit.ordering),
                "prediction": _prediction_dict(fit),
                "whitened_residual": np.asarray(
                    fit.whitened_residual,
                    dtype=float,
                ).tolist(),
            }
            successful.append(entry)
            successful.sort(key=lambda item: float(item["chi2"]))
            del successful[max(1, int(keep_best)):]

            if index == 0 or successful[0]["index"] == index:
                print(
                    f"  {record.name}: Sobol point {index + 1}/{n_points} "
                    f"new best chi2={chi2:.8g}",
                    flush=True,
                )

            if chi2 <= target_chi2:
                stopped_early = True
                break

        except Exception as exc:
            failed_count += 1
            if failed_count <= 3:
                message = str(exc).replace("\n", " ")
                if len(message) > 140:
                    message = message[:137] + "..."
                print(
                    f"  {record.name}: Sobol point {index + 1}/{n_points} "
                    f"failed: {message}",
                    flush=True,
                )

    if not successful:
        raise RuntimeError(
            "Sobol benchmark search produced no successful numerical points."
        )

    best = successful[0]
    evaluated_count = int(best["index"]) + 1 if stopped_early else n_points

    result = {
        "kind": BENCHMARK_SCHEMA,
        "status": (
            "Accepted"
            if float(best["chi2"]) <= target_chi2
            else "BestAvailable"
        ),
        "method": "scrambled_sobol",
        "compatibility": dict(signature),
        "evaluated_point_count": int(evaluated_count),
        "failed_point_count": int(failed_count),
        "stopped_early": bool(stopped_early),
        "best": best,
        "retained_best": successful,
    }

    _atomic_write_json(cache_path, result)

    print(
        f"  {record.name}: Sobol benchmark "
        f"{result['status']} chi2={float(best['chi2']):.8g} "
        f"-> {cache_path}",
        flush=True,
    )

    return result


def _local_refinement_config(payload: Mapping[str, Any]) -> dict[str, Any]:
    search = payload.get("benchmark_search", {})
    if not isinstance(search, Mapping):
        raise ValueError("benchmark_search must be an object.")

    raw = search.get("local_refinement", {})
    if not isinstance(raw, Mapping):
        raise ValueError("benchmark_search.local_refinement must be an object.")

    enabled = bool(raw.get("enabled", True))
    seed_count = int(raw.get("seed_count", 5))
    active_counts = tuple(int(x) for x in raw.get("active_counts", [5, 8, 12]))
    max_nfev = tuple(int(x) for x in raw.get("max_nfev", [30, 40, 60]))
    sensitivity_step = float(raw.get("sensitivity_step", 0.02))

    if seed_count <= 0:
        raise ValueError("local_refinement.seed_count must be positive.")
    if not active_counts or any(x <= 0 for x in active_counts):
        raise ValueError("local_refinement.active_counts must be positive.")
    if len(max_nfev) != len(active_counts):
        raise ValueError(
            "local_refinement.max_nfev must match active_counts in length."
        )
    if any(x <= 0 for x in max_nfev):
        raise ValueError("local_refinement.max_nfev values must be positive.")
    if not (0.0 < sensitivity_step <= 0.25):
        raise ValueError(
            "local_refinement.sensitivity_step must satisfy 0 < step <= 0.25."
        )

    return {
        "enabled": enabled,
        "seed_count": seed_count,
        "active_counts": list(active_counts),
        "max_nfev": list(max_nfev),
        "sensitivity_step": sensitivity_step,
    }


def _optimizer_parameters(
    definitions: Sequence[BenchmarkParameter],
) -> tuple[OptimizerParameter, ...]:
    return tuple(
        OptimizerParameter(
            name=parameter.name,
            low=parameter.low,
            high=parameter.high,
            scale=parameter.scale,
        ).validated()
        for parameter in definitions
    )


def _read_json_cache(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    return payload if isinstance(payload, dict) else None


def _signature_without_local(signature: Mapping[str, Any]) -> dict[str, Any]:
    result = dict(signature)
    result.pop("local_refinement", None)
    return result


def _completed_cache(
    path: Path,
    signature: Mapping[str, Any],
) -> dict[str, Any] | None:
    payload = _read_json_cache(path)
    if payload is None:
        return None
    if payload.get("kind") != BENCHMARK_SCHEMA:
        return None
    if payload.get("compatibility") != dict(signature):
        return None
    best = payload.get("best")
    if not isinstance(best, dict) or not isinstance(best.get("parameters"), dict):
        return None
    return payload


def _reusable_sobol_phase(
    path: Path,
    base_signature: Mapping[str, Any],
) -> dict[str, Any] | None:
    """Reuse the expensive v1 Sobol sample when only refinement is new."""
    payload = _read_json_cache(path)
    if payload is None:
        return None

    kind = payload.get("kind")
    compatibility = payload.get("compatibility")
    if not isinstance(compatibility, dict):
        return None

    if kind == LEGACY_SOBOL_SCHEMA:
        comparable = compatibility
    elif kind == BENCHMARK_SCHEMA:
        comparable = _signature_without_local(compatibility)
    else:
        return None

    if comparable != dict(base_signature):
        return None

    retained = payload.get("retained_best")
    if not isinstance(retained, list) or not retained:
        sobol = payload.get("sobol")
        if isinstance(sobol, dict):
            retained = sobol.get("retained_best")
    if not isinstance(retained, list) or not retained:
        return None

    return payload


def _optimizer_candidate(best: Mapping[str, Any], source: str) -> dict[str, Any]:
    return {
        "index": None,
        "chi2": float(best["chi2"]),
        "parameters": {
            str(name): float(value)
            for name, value in best["parameters"].items()
        },
        "ordering": best.get("ordering"),
        "prediction": best.get("prediction"),
        "whitened_residual": best.get("whitened_residual"),
        "source": source,
    }


def _run_local_stage(
    *,
    evaluator: Callable[[Mapping[str, float]], Any],
    observable_names: Sequence[str],
    parameters: Sequence[OptimizerParameter],
    seed: Mapping[str, Any],
    active_count: int,
    max_nfev: int,
    sensitivity_step: float,
    label: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    seed_z = unit_vector_from_parameters(parameters, seed["parameters"])
    recorder = EvaluationRecorder(
        evaluator=evaluator,
        parameters=parameters,
        observable_names=observable_names,
        checkpoint_path=None,
    )
    baseline, sensitivity = estimate_local_sensitivities(
        recorder,
        seed_z,
        step=sensitivity_step,
    )
    solution, active_names = optimise_reduced_least_squares(
        recorder,
        seed_z,
        sensitivity,
        active_count=active_count,
        max_nfev=max_nfev,
    )
    if recorder.best is None:
        raise RuntimeError(f"Local stage {label} produced no successful point.")

    candidate = _optimizer_candidate(recorder.best, label)
    stage = {
        "stage": label,
        "seed_chi2": float(seed["chi2"]),
        "baseline_recomputed_chi2": float(np.dot(baseline, baseline)),
        "active_parameter_count": int(active_count),
        "active_parameters": list(active_names),
        "max_nfev": int(max_nfev),
        "optimizer_nfev": int(solution.nfev),
        "optimizer_success": bool(solution.success),
        "optimizer_message": str(solution.message),
        "pipeline_evaluations": len(recorder.history),
        "best_chi2": float(candidate["chi2"]),
        "sensitivity_ranking": [
            {
                "name": item.name,
                "score": float(item.score),
                "status": item.status,
            }
            for item in sensitivity
        ],
    }
    return candidate, stage


def refine_sobol_candidates(
    *,
    evaluator: Callable[[Mapping[str, float]], Any],
    observable_names: Sequence[str],
    definitions: Sequence[BenchmarkParameter],
    candidates: Sequence[Mapping[str, Any]],
    target_chi2: float,
    seed_count: int,
    active_counts: Sequence[int],
    max_nfev: Sequence[int],
    sensitivity_step: float,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Refine several Sobol seeds cheaply, then deepen only the best branch."""
    ordered = sorted(
        (dict(candidate) for candidate in candidates),
        key=lambda item: float(item["chi2"]),
    )
    if not ordered:
        raise ValueError("At least one Sobol candidate is required.")

    selected = ordered[: min(int(seed_count), len(ordered))]
    params = _optimizer_parameters(definitions)
    counts = tuple(int(x) for x in active_counts)
    budgets = tuple(int(x) for x in max_nfev)
    if not counts or len(counts) != len(budgets):
        raise ValueError("active_counts and max_nfev must be non-empty and aligned.")
    if any(count > len(params) for count in counts):
        raise ValueError("active_count cannot exceed the benchmark parameter count.")

    if float(ordered[0]["chi2"]) <= float(target_chi2):
        return ordered[0], {
            "enabled": True,
            "skipped": True,
            "reason": "Sobol stage already reached target chi2.",
            "seed_count": len(selected),
            "stopped_early": True,
            "stages": [],
            "pipeline_evaluations": 0,
        }

    history: list[dict[str, Any]] = []
    all_candidates = list(ordered)
    first_refined: list[dict[str, Any]] = []

    first_count = counts[0]
    first_budget = budgets[0]
    for rank, seed in enumerate(selected, start=1):
        label = f"multistart_seed_{rank}_active_{first_count}"
        candidate, stage = _run_local_stage(
            evaluator=evaluator,
            observable_names=observable_names,
            parameters=params,
            seed=seed,
            active_count=first_count,
            max_nfev=first_budget,
            sensitivity_step=sensitivity_step,
            label=label,
        )
        first_refined.append(candidate)
        all_candidates.append(candidate)
        history.append(stage)
        print(
            f"    {label}: {float(seed['chi2']):.8g} -> "
            f"{float(candidate['chi2']):.8g}",
            flush=True,
        )
        if float(candidate["chi2"]) <= float(target_chi2):
            best = min(all_candidates, key=lambda x: float(x["chi2"]))
            return best, {
                "enabled": True,
                "skipped": False,
                "seed_count": len(selected),
                "stopped_early": True,
                "stages": history,
                "pipeline_evaluations": sum(
                    int(item["pipeline_evaluations"]) for item in history
                ),
            }

    branch = min(first_refined, key=lambda x: float(x["chi2"]))
    for count, budget in zip(counts[1:], budgets[1:]):
        label = f"best_branch_active_{count}"
        candidate, stage = _run_local_stage(
            evaluator=evaluator,
            observable_names=observable_names,
            parameters=params,
            seed=branch,
            active_count=count,
            max_nfev=budget,
            sensitivity_step=sensitivity_step,
            label=label,
        )
        history.append(stage)
        all_candidates.append(candidate)
        print(
            f"    {label}: {float(branch['chi2']):.8g} -> "
            f"{float(candidate['chi2']):.8g}",
            flush=True,
        )
        if float(candidate["chi2"]) < float(branch["chi2"]):
            branch = candidate
        if float(branch["chi2"]) <= float(target_chi2):
            break

    best = min(all_candidates, key=lambda x: float(x["chi2"]))
    return best, {
        "enabled": True,
        "skipped": False,
        "seed_count": len(selected),
        "stopped_early": float(best["chi2"]) <= float(target_chi2),
        "stages": history,
        "pipeline_evaluations": sum(
            int(item["pipeline_evaluations"]) for item in history
        ),
    }


def resolve_model_benchmark(
    *,
    record: RunRecord,
    payload: Mapping[str, Any],
    uv_rgbeta_path: Path,
    intermediate_rgbeta_path: Path,
    final_weinberg_path: Path,
) -> BenchmarkResolution:
    """Load or find one model benchmark using Sobol plus local refinement."""
    prepared = retarget_payload_to_record(record, payload)
    search = prepared.get("benchmark_search", {})
    if not isinstance(search, Mapping):
        raise ValueError("benchmark_search must be an object.")
    if not bool(search.get("enabled", False)):
        return BenchmarkResolution(
            payload=prepared,
            metadata={"status": "Disabled", "method": None},
            cache_path=None,
        )
    if record.shared_scalar:
        raise NotImplementedError(
            "Automatic benchmark search currently supports only the ordinary "
            "split-scalar T3 branch."
        )
    if str(search.get("method", "sobol")).lower() != "sobol":
        raise ValueError("Automatic benchmark search currently supports Sobol only.")

    n_points = int(search.get("n_points", 256))
    target_chi2 = float(search.get("target_chi2", 10.0))
    keep_best = int(search.get("keep_best", 10))
    base_seed = int(search.get("seed", 20260927))
    local = _local_refinement_config(prepared)
    if target_chi2 <= 0.0 or not np.isfinite(target_chi2):
        raise ValueError("benchmark_search.target_chi2 must be positive and finite.")
    if keep_best <= 0:
        raise ValueError("benchmark_search.keep_best must be positive.")
    if int(local["seed_count"]) > keep_best:
        raise ValueError("local_refinement.seed_count cannot exceed keep_best.")

    definitions = _parameters_from_config(prepared)
    target_path = _resolve_target_path(prepared)
    effective_seed = _effective_seed(base_seed, model_benchmark_key(record))
    base_signature = _compatibility_signature(
        record=record,
        payload=prepared,
        parameters=definitions,
        target_path=target_path,
        uv_rgbeta_path=Path(uv_rgbeta_path),
        intermediate_rgbeta_path=Path(intermediate_rgbeta_path),
        final_weinberg_path=Path(final_weinberg_path),
        n_points=n_points,
        target_chi2=target_chi2,
        effective_seed=effective_seed,
    )
    signature = dict(base_signature)
    signature["local_refinement"] = dict(local)
    cache_path = _cache_path(record)

    completed = _completed_cache(cache_path, signature)
    if completed is not None:
        best = completed["best"]
        resolved = apply_benchmark_parameters(
            prepared, best["parameters"], definitions
        )
        print(
            f"  {record.name}: using cached refined benchmark "
            f"chi2={float(best['chi2']):.8g} -> {cache_path}",
            flush=True,
        )
        return BenchmarkResolution(
            payload=resolved,
            metadata={
                "status": "Cached",
                "method": "sobol_plus_local_least_squares",
                "chi2": float(best["chi2"]),
                "benchmark_status": str(completed.get("status")),
                "cache_path": str(cache_path),
            },
            cache_path=cache_path,
        )

    evaluator = _build_evaluator(
        payload=prepared,
        parameters=definitions,
        uv_rgbeta_path=Path(uv_rgbeta_path),
        intermediate_rgbeta_path=Path(intermediate_rgbeta_path),
        final_weinberg_path=Path(final_weinberg_path),
        target_path=target_path,
    )
    target = OscillationFitTarget.from_json(target_path)

    reusable = _reusable_sobol_phase(cache_path, base_signature)
    if reusable is not None:
        retained = reusable.get("retained_best")
        if not isinstance(retained, list):
            retained = reusable["sobol"]["retained_best"]
        print(
            f"  {record.name}: reusing existing Sobol phase "
            f"({len(retained)} retained seeds).",
            flush=True,
        )
        sobol = {
            "evaluated_point_count": int(reusable.get("evaluated_point_count", 0)),
            "failed_point_count": int(reusable.get("failed_point_count", 0)),
            "stopped_early": bool(reusable.get("stopped_early", False)),
            "best": dict(retained[0]),
            "retained_best": list(retained),
            "reused_from_cache": True,
        }
    else:
        raw = _run_search(
            record=record,
            payload=prepared,
            definitions=definitions,
            evaluator=evaluator,
            cache_path=cache_path,
            signature=base_signature,
            n_points=n_points,
            target_chi2=target_chi2,
            keep_best=keep_best,
            effective_seed=effective_seed,
        )
        sobol = {
            "evaluated_point_count": int(raw["evaluated_point_count"]),
            "failed_point_count": int(raw["failed_point_count"]),
            "stopped_early": bool(raw["stopped_early"]),
            "best": dict(raw["best"]),
            "retained_best": list(raw["retained_best"]),
            "reused_from_cache": False,
        }

    best = dict(sobol["best"])
    local_result = {
        "enabled": bool(local["enabled"]),
        "skipped": True,
        "reason": "Local refinement disabled.",
        "stages": [],
        "pipeline_evaluations": 0,
    }
    if bool(local["enabled"]) and float(best["chi2"]) > target_chi2:
        print(
            f"  {record.name}: refining top {local['seed_count']} Sobol seeds "
            f"with active counts {local['active_counts']}...",
            flush=True,
        )
        refined, local_result = refine_sobol_candidates(
            evaluator=evaluator,
            observable_names=target.observable_names,
            definitions=definitions,
            candidates=sobol["retained_best"],
            target_chi2=target_chi2,
            seed_count=int(local["seed_count"]),
            active_counts=tuple(local["active_counts"]),
            max_nfev=tuple(local["max_nfev"]),
            sensitivity_step=float(local["sensitivity_step"]),
        )
        if float(refined["chi2"]) < float(best["chi2"]):
            best = refined

    result = {
        "kind": BENCHMARK_SCHEMA,
        "status": "Accepted" if float(best["chi2"]) <= target_chi2 else "BestAvailable",
        "method": "sobol_plus_local_least_squares",
        "compatibility": signature,
        "sobol": sobol,
        "retained_best": sobol["retained_best"],
        "local_refinement": local_result,
        "best": best,
    }
    _atomic_write_json(cache_path, result)
    print(
        f"  {record.name}: benchmark {result['status']} "
        f"chi2={float(best['chi2']):.8g} -> {cache_path}",
        flush=True,
    )

    resolved = apply_benchmark_parameters(prepared, best["parameters"], definitions)
    return BenchmarkResolution(
        payload=resolved,
        metadata={
            "status": "Searched",
            "method": "sobol_plus_local_least_squares",
            "chi2": float(best["chi2"]),
            "benchmark_status": str(result["status"]),
            "cache_path": str(cache_path),
            "sobol_reused_from_cache": bool(sobol["reused_from_cache"]),
            "evaluated_point_count": int(sobol["evaluated_point_count"]),
            "failed_point_count": int(sobol["failed_point_count"]),
            "local_pipeline_evaluations": int(
                local_result.get("pipeline_evaluations", 0)
            ),
        },
        cache_path=cache_path,
    )
