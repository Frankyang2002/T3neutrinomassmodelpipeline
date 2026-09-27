from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from Numerical.OscillationFit import OscillationFitResult

import numpy as np

from Numerical.SobolBenchmarkSearch import (
    BenchmarkParameter,
    apply_benchmark_parameters,
    default_ordinary_t3_parameters,
    generate_sobol_parameter_points,
    refine_sobol_candidates,
    retarget_payload_to_record,
)


def test_sobol_points_are_reproducible_and_bounded() -> None:
    parameters = (
        BenchmarkParameter(
            name="linear",
            binding="ordinary.x",
            low=-2.0,
            high=3.0,
            scale="linear",
        ),
        BenchmarkParameter(
            name="log",
            binding="ordinary.y",
            low=1.0e-3,
            high=1.0e1,
            scale="log",
        ),
    )

    first = generate_sobol_parameter_points(
        parameters,
        64,
        seed=12345,
    )
    second = generate_sobol_parameter_points(
        parameters,
        64,
        seed=12345,
    )

    assert first == second
    assert len(first) == 64

    for point in first:
        assert -2.0 <= point["linear"] < 3.0
        assert 1.0e-3 <= point["log"] < 1.0e1


def test_default_profile_contains_19_real_parameters() -> None:
    parameters = default_ordinary_t3_parameters()

    assert len(parameters) == 19
    assert parameters[0].name == "lambdaT3_real"
    assert sum(parameter.name.startswith("y1_") for parameter in parameters) == 9
    assert sum(parameter.name.startswith("y2_") for parameter in parameters) == 9


def test_benchmark_parameters_are_applied_to_base_state() -> None:
    payload = {
        "base_state": {
            "ordinary": {
                "lambdaT3": {"real": 0.005, "imag": 0.002},
                "y1": {
                    "real": [[0.0, 0.0], [0.0, 0.0]],
                },
            }
        }
    }
    definitions = (
        BenchmarkParameter(
            name="lambdaT3_real",
            binding="ordinary.lambdaT3.real",
            low=1.0e-3,
            high=3.0e-1,
            scale="log",
        ),
        BenchmarkParameter(
            name="y1_12",
            binding="ordinary.y1.real.0.1",
            low=-0.5,
            high=0.5,
            scale="linear",
        ),
    )

    resolved = apply_benchmark_parameters(
        payload,
        {
            "lambdaT3_real": 0.12,
            "y1_12": -0.23,
        },
        definitions,
    )

    assert resolved["base_state"]["ordinary"]["lambdaT3"]["real"] == 0.12
    assert resolved["base_state"]["ordinary"]["lambdaT3"]["imag"] == 0.002
    assert resolved["base_state"]["ordinary"]["y1"]["real"][0][1] == -0.23


def test_retarget_is_opt_in() -> None:
    record = SimpleNamespace(
        d_s1=3,
        d_s2=4,
        d_f=2,
        alpha=1,
        shared_scalar=False,
    )
    payload = {
        "representation": {
            "d_s1": 2,
            "d_s2": 2,
            "d_f": 1,
            "alpha": -1,
            "shared_scalar": False,
        },
        "benchmark_search": {
            "enabled": True,
            "use_current_model_representation": True,
        },
    }

    resolved = retarget_payload_to_record(record, payload)

    assert resolved["representation"] == {
        "d_s1": 3,
        "d_s2": 4,
        "d_f": 2,
        "alpha": 1,
        "shared_scalar": False,
    }


def _synthetic_fit(parameters) -> OscillationFitResult:
    x = float(parameters["x"])
    y = float(parameters["y"])
    residual = np.array([5.0 * (x - 0.2), 2.0 * (y + 0.35)], dtype=float)
    return OscillationFitResult(
        ordering="NO",
        chi2=float(np.dot(residual, residual)),
        observable_names=("sin2_theta12", "sin2_theta13"),
        prediction=residual.copy(),
        target=np.zeros(2, dtype=float),
        residual=residual.copy(),
        whitened_residual=residual.copy(),
        source="synthetic",
    )


def test_local_refinement_improves_sobol_seed() -> None:
    definitions = (
        BenchmarkParameter("x", "ordinary.x", -1.0, 1.0, "linear"),
        BenchmarkParameter("y", "ordinary.y", -1.0, 1.0, "linear"),
    )
    candidates = [
        {"index": 0, "chi2": _synthetic_fit({"x": 0.8, "y": 0.4}).chi2,
         "parameters": {"x": 0.8, "y": 0.4}},
        {"index": 1, "chi2": _synthetic_fit({"x": -0.4, "y": -0.1}).chi2,
         "parameters": {"x": -0.4, "y": -0.1}},
    ]
    best, details = refine_sobol_candidates(
        evaluator=_synthetic_fit,
        observable_names=("sin2_theta12", "sin2_theta13"),
        definitions=definitions,
        candidates=candidates,
        target_chi2=1.0e-8,
        seed_count=2,
        active_counts=(1, 2),
        max_nfev=(20, 30),
        sensitivity_step=0.02,
    )
    assert best["chi2"] < min(item["chi2"] for item in candidates)
    assert best["chi2"] < 1.0e-6
    assert details["stages"]
