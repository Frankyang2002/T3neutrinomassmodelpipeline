"""Regression tests for additive opt-in orchestration (no production edits)."""
import importlib.util
import json
import sys
import types
from pathlib import Path

import numpy as np
import pytest


MOD_PATH = Path(__file__).resolve().parents[1] / "Numerical" / "orchestration" / "T3BIntegratedCoupledPipeline.py"
spec = importlib.util.spec_from_file_location("test_integrated_module", MOD_PATH)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def config():
    return {"kind": "t3_pipeline_numerical_results_v1",
            "representation": {"d_s1": 2, "d_s2": 2, "d_f": 1,
                               "alpha": -1, "shared_scalar": False},
            "benchmark_search": {"enabled": False}}


def test_matrix_serialization_preserves_complex_phases_and_shape():
    a = np.eye(3, dtype=complex) * (1.25 + 3.5j)
    b = mod._matrix_read(mod._matrix_json(a))
    assert np.array_equal(a, b)
    assert mod._relative(a, b) == 0


def test_model_guard_allows_fitted_t3b():
    assert mod._model_check(config()) == (2, 2, 1, -1, False)


def test_model_guard_rejects_other_or_sobol():
    c = config()
    c["representation"]["d_f"] = 3
    with pytest.raises(ValueError, match="only implemented"):
        mod._model_check(c)
    c = config()
    c["benchmark_search"]["enabled"] = True
    with pytest.raises(ValueError, match="fixed"):
        mod._model_check(c)


def test_background_rows_are_ordered_and_complex():
    class Point:
        def __init__(self, scale):
            self.mu_gev = scale
            self.sm = types.SimpleNamespace(g2=0.65, lambdaH=0.2,
                 ye=np.eye(3), yu=np.eye(3), yd=np.eye(3))
            self.lambda12, self.lambda12Adj = 0.4, 0.3
            self.lambdaT3 = 0.12+0.03j
    class Eft:
        n_points=3
        def state_at_index(self, n): return Point((1e5, 1e4, 1e3)[n])
    rows = mod._intermediate_samples(types.SimpleNamespace(eft1=Eft()))
    assert [r["mu_gev"] for r in rows] == [1e5, 1e4, 1e3]
    assert rows[0]["lambdaT3"] == [0.12, 0.03]
    assert np.array_equal(mod._matrix_read(rows[0]["ye"]), np.eye(3))


def install_stub_pipeline(monkeypatch, fake):
    # Import only on invoking run(); no official production module is modified.
    pkg = types.ModuleType("Numerical.orchestration")
    pkg.__path__ = []
    stub = types.ModuleType("Numerical.orchestration.PipelineNumericalResults")
    stub.run_pipeline_numerical_results = fake
    monkeypatch.setitem(sys.modules, "Numerical", types.ModuleType("Numerical"))
    monkeypatch.setitem(sys.modules, "Numerical.orchestration", pkg)
    monkeypatch.setitem(sys.modules, "Numerical.orchestration.PipelineNumericalResults", stub)
    records_pkg = types.ModuleType("common")
    records_pkg.__path__ = []
    record_stub = types.ModuleType("common.RunRecords")
    class RunRecord:
        def __init__(self, **kwargs): self.__dict__.update(kwargs)
    record_stub.RunRecord = RunRecord
    monkeypatch.setitem(sys.modules, "common", records_pkg)
    monkeypatch.setitem(sys.modules, "common.RunRecords", record_stub)


def test_disabled_flag_calls_only_official_pipeline(monkeypatch, tmp_path):
    conf = tmp_path / "input.json"; conf.write_text(json.dumps(config()))
    outdir = tmp_path / "output"
    seen = []
    def production(record, path):
        seen.append((record, path))
        return True
    install_stub_pipeline(monkeypatch, production)
    monkeypatch.setattr(mod, "compute_opt_in_coupled", lambda **k: pytest.fail("Unexpected opt-in calculation"))
    report = mod.run(config=conf, output_dir=outdir, enable_coupled=False)
    assert report["status"] == "IntegratedFixedOrderComplete"
    assert report["coupled_ll_enabled"] is False
    assert len(seen) == 1
    assert not (outdir / "data" / "t3b_coupled_ll_exploratory.json").exists()


def test_enabled_writes_separate_artifact_and_leaves_fixed_unmodified(monkeypatch, tmp_path):
    conf = tmp_path / "input.json"; conf.write_text(json.dumps(config()))
    outdir = tmp_path / "output"
    official_path = outdir / "data" / "running_diagnostics.json"
    content = '{"status":"Success","production":true}\n'
    def production(record, path):
        official_path.parent.mkdir(parents=True, exist_ok=True)
        official_path.write_text(content)
        return True
    install_stub_pipeline(monkeypatch, production)
    seen = {}
    def coupled(**kwargs):
        seen.update(kwargs)
        return {"coupled_at_ms": {"relative_difference_from_fixed_order": .0249,
                                  "relative_reverse_feedback": .000015},
                "low_energy": {"relative_c5_difference_at_low_scale": .0249},
                "perturbative_accuracy": {"bookkeeping": {"terms": {
                    "order_2": {"norm_relative_to_order1": .033},
                    "order_3": {"reverse_contribution_relative_to_order1": .000015}
                }}}}
    monkeypatch.setattr(mod, "compute_opt_in_coupled", coupled)
    report = mod.run(config=conf, output_dir=outdir, enable_coupled=True)
    assert report["status"] == "IntegratedFixedAndExploratoryCoupledComplete"
    assert report["relative_c5_change_at_ms"] == .0249
    assert seen["production"]["production"] is True
    assert official_path.read_text() == content
    extra = json.loads((outdir / "data" / "t3b_coupled_ll_exploratory.json").read_text())
    assert extra["coupled_at_ms"]["relative_difference_from_fixed_order"] == .0249


def test_pipeline_failure_stops_before_coupled(monkeypatch, tmp_path):
    conf = tmp_path / "input.json"; conf.write_text(json.dumps(config()))
    def failure(record, path):
        record.summary["PipelineNumericalResultsError"] = "test failure"
        return False
    install_stub_pipeline(monkeypatch, failure)
    with pytest.raises(RuntimeError, match="test failure"):
        mod.run(config=conf, output_dir=tmp_path / "output", enable_coupled=True)
