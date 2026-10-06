from __future__ import annotations

import json
from pathlib import Path

from Numerical.orchestration.PipelineNumericalResults import (
    PIPELINE_NUMERICAL_CONFIG_KIND,
    is_pipeline_numerical_results_config,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = (
    PROJECT_ROOT / "configs" / "t3_numerical_default_B_alpha_m1.json"
)


def test_default_config_uses_pipeline_numerical_schema() -> None:
    payload = json.loads(DEFAULT_CONFIG.read_text(encoding="utf-8"))

    assert payload["kind"] == PIPELINE_NUMERICAL_CONFIG_KIND
    assert payload["representation"] == {
        "d_s1": 2,
        "d_s2": 2,
        "d_f": 1,
        "alpha": -1,
        "shared_scalar": False,
    }
    assert payload["base_state"]["mu_gev"] == 1.0e7
    assert payload["scales"]["mu_fermion_threshold_gev"] == 1.0e5
    assert payload["scales"]["mu_scalar_threshold_gev"] == 1.0e3
    assert payload["scales"]["mu_low_gev"] == 1.0e2
    assert payload["numerical"]["n_scale_points"] >= 2


def test_default_config_is_detected_as_integrated_results_config() -> None:
    assert is_pipeline_numerical_results_config(DEFAULT_CONFIG)


def test_legacy_or_unrelated_json_is_not_detected(tmp_path: Path) -> None:
    path = tmp_path / "legacy.json"
    path.write_text(
        json.dumps({"mu_initial_gev": 1.0e5}),
        encoding="utf-8",
    )

    assert not is_pipeline_numerical_results_config(path)


def test_default_config_enables_automatic_sobol_benchmark_search() -> None:
    payload = json.loads(DEFAULT_CONFIG.read_text(encoding="utf-8"))
    search = payload["benchmark_search"]

    assert search["enabled"] is True
    assert search["method"] == "sobol"
    assert search["use_current_model_representation"] is True
    assert search["n_points"] == 256
    assert search["target_chi2"] == 10.0
    assert search["oscillation_target"].endswith(
        "nufit_6_1_ic24_no.json"
    )
    local = search["local_refinement"]
    assert local["enabled"] is True
    assert local["seed_count"] == 5
    assert local["active_counts"] == [5, 8, 12]
    assert local["max_nfev"] == [30, 40, 60]


def test_default_config_enables_controlled_sensitivity_scan() -> None:
    payload = json.loads(DEFAULT_CONFIG.read_text(encoding="utf-8"))
    sensitivity = payload["sensitivity"]

    assert sensitivity["enabled"] is True
    assert sensitivity["gamma_min"] == 0.95
    assert sensitivity["gamma_max"] == 1.05
    assert sensitivity["n_points"] == 11
    assert sensitivity["oscillation_target"].endswith(
        "nufit_6_1_ic24_no.json"
    )
