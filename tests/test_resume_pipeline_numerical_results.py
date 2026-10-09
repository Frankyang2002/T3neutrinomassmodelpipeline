"""Fast tests of numerical-only resume discovery and safety gates."""
from __future__ import annotations
import json
from pathlib import Path

from Numerical.orchestration.ResumePipelineNumericalResults import (
    REQUIRED, discover_cases, parse_model_key, resume,
)


def _case(tmp_path: Path, scenario: str = "largeY_largeL", key: str = "T3_dS1_3_dS2_3_dF_2_alpha_p2"):
    output = tmp_path / "output" / "comparison"
    configs = tmp_path / "configs" / "comparison"
    path = output / scenario / key / "T3_E_alpha_p2"
    (path / "data").mkdir(parents=True)
    config = configs / scenario / "input" / f"{key}.json"
    config.parent.mkdir(parents=True)
    config.write_text(json.dumps({"kind": "t3_pipeline_numerical_results_v1",
        "representation": {"d_s1": 3, "d_s2": 3, "d_f": 2, "alpha": 2},
        "benchmark_search": {"enabled": False}}), encoding="utf-8")
    return output, configs, path


def test_parse_and_discover(tmp_path):
    root, _, _ = _case(tmp_path)
    assert parse_model_key("T3_dS1_3_dS2_3_dF_2_alpha_m4") == (3, 3, 2, -4)
    assert len(discover_cases(root)) == 1


def test_dry_run_and_missing_artifacts(tmp_path):
    root, configs, case = _case(tmp_path)
    calls = []
    fake = lambda *_args: calls.append(1) or True
    first = resume(root, configs, execute=True, runner=fake)
    assert first["counts"] == {"MissingArtifacts": 1}
    for filename in REQUIRED:
        (case / "data" / filename).write_text("{}", encoding="utf-8")
    planned = resume(root, configs, runner=fake)
    assert planned["counts"] == {"Ready": 1}
    assert calls == []
    done = resume(root, configs, execute=True, limit=1, runner=fake)
    assert done["counts"] == {"Success": 1} and len(calls) == 1


def test_skip_complete_and_reject_search(tmp_path):
    root, configs, case = _case(tmp_path)
    for filename in REQUIRED:
        (case / "data" / filename).write_text("{}", encoding="utf-8")
    (case / "data" / "running_diagnostics.json").write_text('{"status":"Success"}', encoding="utf-8")
    assert resume(root, configs)["counts"] == {"AlreadyComplete": 1}
    cfg = next(configs.rglob("*.json"))
    obj = json.loads(cfg.read_text())
    obj["benchmark_search"]["enabled"] = True
    cfg.write_text(json.dumps(obj), encoding="utf-8")
    assert resume(root, configs, force=True)["counts"] == {"InvalidConfig": 1}


def test_resume_disables_sensitivity_without_changing_original(tmp_path):
    root, configs, case = _case(tmp_path)
    for filename in REQUIRED:
        (case / "data" / filename).write_text("{}", encoding="utf-8")
    original = next(configs.rglob("*.json"))
    payload = json.loads(original.read_text(encoding="utf-8"))
    payload["sensitivity"] = {"enabled": True, "oscillation_target": "configs/nufit_6_1_ic24_no.json"}
    original.write_text(json.dumps(payload), encoding="utf-8")
    seen = []
    def fake_runner(_record, config):
        effective = json.loads(config.read_text(encoding="utf-8"))
        seen.append(config)
        assert effective["sensitivity"]["enabled"] is False
        assert effective["benchmark_search"]["enabled"] is False
        assert config != original
        return True
    result = resume(root, configs, execute=True, limit=1, runner=fake_runner)
    assert result["counts"] == {"Success": 1}
    assert len(seen) == 1 and seen[0].is_file()
    assert json.loads(original.read_text(encoding="utf-8"))["sensitivity"]["enabled"] is True


def test_input_preferred_to_stale_resolved(tmp_path):
    root, configs, case = _case(tmp_path)
    for filename in REQUIRED:
        (case / "data" / filename).write_text("{}", encoding="utf-8")
    input_file = next(configs.rglob("*.json"))
    original = json.loads(input_file.read_text(encoding="utf-8"))
    original["base_state"] = {"ordinary": {"y1": {"real": [[0.005]]},
                                          "lambdaT3": {"real": 0.1}}}
    input_file.write_text(json.dumps(original), encoding="utf-8")
    stale = input_file.parent.parent / "resolved" / input_file.name
    stale.parent.mkdir(parents=True)
    corrupted = json.loads(json.dumps(original))
    corrupted["base_state"]["ordinary"]["y1"]["real"][0][0] = -0.20606645548
    corrupted["base_state"]["ordinary"]["lambdaT3"]["real"] = 0.001
    stale.write_text(json.dumps(corrupted), encoding="utf-8")

    def check(_record, path):
        used = json.loads(path.read_text(encoding="utf-8"))
        assert used["base_state"]["ordinary"]["y1"]["real"][0][0] == 0.005
        assert used["base_state"]["ordinary"]["lambdaT3"]["real"] == 0.1
        return True

    result = resume(root, configs, execute=True, force=True, runner=check)
    assert result["counts"] == {"Success": 1}
