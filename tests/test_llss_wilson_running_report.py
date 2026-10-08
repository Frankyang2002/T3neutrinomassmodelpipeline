"""Tests for plotting LLSS transport separately from numerical C5."""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from Numerical.plotting import LLSSWilsonRunningReport as report


def _transport(**kwargs):
    ratio = math.log(float(kwargs["mu_low"]) / float(kwargs["mu_high"]))
    return {
        "status": "Success", "tree_boundary_component_count": 1,
        "beta_component_count": 2, "running_correction_component_count": 2,
        "generated_by_running_component_count": 1,
        "tree_boundary_components": {"0,0,0,0": "2*y**2"},
        "one_loop_running_components": {
            "0,0,0,0": str(3 * ratio) + "*y**2",
            "0,0,1,1": str(-ratio) + "*y**2",
        },
    }


def _source(tmp_path: Path) -> Path:
    data = tmp_path / "comparison" / "smallY_smallL" / "T3_A_alpha_m2" / "data"
    data.mkdir(parents=True)
    for name in report.REQUIRED:
        (data / name).write_text("{}", encoding="utf-8")
    return data


def test_generate_symbolic_report_without_inventing_numeric_values(tmp_path, monkeypatch):
    monkeypatch.setattr(report, "run_component_wilson_transport", _transport)
    data = _source(tmp_path)
    out = tmp_path / "reports" / "T3_A_alpha_m2"
    result = report.report_model(data, out, mu_f=1e5, mu_s=1e3)
    assert result.is_file()
    assert (out / "llss_transport_kernel.png").is_file()
    assert (out / "llss_component_support.png").is_file()
    assert not (out / "llss_component_running.png").exists()
    payload = json.loads(result.read_text())
    assert payload["numerical_assignments"] is None
    assert payload["numeric_component_curves"] == []
    assert "No numerical component curves" in (out / "LLSS_README.md").read_text(encoding="utf-8")


def test_complete_explicit_values_produce_component_plots(tmp_path, monkeypatch):
    monkeypatch.setattr(report, "run_component_wilson_transport", _transport)
    data = _source(tmp_path)
    out = tmp_path / "reports" / "T3_A_alpha_m2"
    result = report.report_model(data, out, mu_f=1e5, mu_s=1e3, values={"y": 0.1})
    assert (out / "llss_component_running.png").is_file()
    payload = json.loads(result.read_text())
    assert len(payload["numeric_component_curves"]) == 2
    assert payload["used_in_authoritative_final_c5"] is False


def test_missing_symbol_rejected_without_defaulting_to_one(tmp_path, monkeypatch):
    monkeypatch.setattr(report, "run_component_wilson_transport", _transport)
    data = _source(tmp_path)
    with pytest.raises(ValueError, match="Missing one-generation assignments: y"):
        report.report_model(data, tmp_path / "reports", mu_f=1e5, mu_s=1e3, values={})


def test_report_caches_matching_sources(tmp_path, monkeypatch):
    calls = []
    def transport(**kwargs):
        calls.append(1)
        return _transport(**kwargs)
    monkeypatch.setattr(report, "run_component_wilson_transport", transport)
    data = _source(tmp_path)
    out = tmp_path / "reports"
    report.report_model(data, out, mu_f=1e5, mu_s=1e3)
    report.report_model(data, out, mu_f=1e5, mu_s=1e3)
    assert len(calls) == 1
    (data / report.REQUIRED[0]).write_text('{"updated":true}', encoding="utf-8")
    report.report_model(data, out, mu_f=1e5, mu_s=1e3)
    assert len(calls) == 2


def test_find_saved_model_without_rerunning_numerical_pipeline(tmp_path, monkeypatch):
    monkeypatch.setattr(report, "run_component_wilson_transport", _transport)
    data = _source(tmp_path)
    results = report.build_reports(data.parents[2], tmp_path / "output", model="T3_A_alpha_m2")
    assert len(results) == 1
    assert results[0].parent.name == "T3_A_alpha_m2"


def test_reports_use_benchmark_key_not_dimension_folder(tmp_path, monkeypatch):
    monkeypatch.setattr(report, "run_component_wilson_transport", _transport)
    data = _source(tmp_path)
    new_data = data.parent.parent / "T3_dS1_1_dS2_3_dF_2_alpha_m2" / "data"
    data.parent.rename(new_data.parent)
    new_data.joinpath("running_diagnostics.json").write_text(json.dumps({
        "status": "Success", "benchmark_search": {"model": {
            "model_key": "T3_A_alpha_m2"}},
    }), encoding="utf-8")
    results = report.build_reports(new_data.parents[2], tmp_path / "reports",
                                   model="T3_A_alpha_m2")
    assert len(results) == 1
    assert results[0].parent.name == "T3_A_alpha_m2"
