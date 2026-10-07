"""Tests for the symbolic fixed-order LLSS diagnostic sampler."""

from __future__ import annotations

from pathlib import Path

import pytest

from Numerical.diagnostics import IntermediateWilsonDiagnostics as diagnostics


def _fake_transport(**kwargs):
    mu_high = float(kwargs["mu_high"])
    mu_low = float(kwargs["mu_low"])
    log_value = 0.0 if mu_low == mu_high else -1.0
    running = {} if log_value == 0.0 else {"0,0,0,0": "-beta0"}
    return {
        "tree_boundary_component_count": 1,
        "beta_component_count": 1,
        "running_correction_component_count": len(running),
        "generated_by_running_component_count": 0,
        "tree_boundary_components": {"0,0,0,0": "C0"},
        "one_loop_running_components": running,
        "full_fixed_order_components": {
            "0,0,0,0": (
                "C0" if log_value == 0.0 else "C0-hbar*beta0"
            )
        },
    }


def test_fixed_order_llss_sampling_is_diagnostic_only(monkeypatch) -> None:
    monkeypatch.setattr(
        diagnostics,
        "run_component_wilson_transport",
        _fake_transport,
    )

    result = diagnostics.sample_fixed_order_llss_transport(
        wilson_seed_path=Path("seed.json"),
        wilson_rge_path=Path("rge.json"),
        rgbeta_path=Path("rgbeta.json"),
        mu_high_gev=1.0e5,
        scales_gev=[1.0e5, 1.0e4, 1.0e3],
    )

    assert result["status"] == "Success"
    assert result["diagnostic_only"] is True
    assert result["rg_improved"] is False
    assert result["used_in_authoritative_final_c5"] is False
    assert result["equal_scale_running_vanishes"] is True
    assert len(result["points"]) == 3
    assert result["points"][0]["one_loop_running_components"] == {}
    assert (
        result["points"][-1]["full_fixed_order_components"]["0,0,0,0"]
        == "C0-hbar*beta0"
    )


def test_fixed_order_llss_sampling_rejects_scale_above_threshold(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        diagnostics,
        "run_component_wilson_transport",
        _fake_transport,
    )

    with pytest.raises(ValueError, match="above the upper threshold"):
        diagnostics.sample_fixed_order_llss_transport(
            wilson_seed_path=Path("seed.json"),
            wilson_rge_path=Path("rge.json"),
            rgbeta_path=Path("rgbeta.json"),
            mu_high_gev=1.0e5,
            scales_gev=[1.0e5, 2.0e5],
        )
