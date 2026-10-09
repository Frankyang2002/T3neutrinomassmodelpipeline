"""Unit tests for the read-only full-comparison observables export."""
from __future__ import annotations
import json
from pathlib import Path
from Numerical.diagnostics.ExtractFullComparisonResults import collect, extract_case, write_results


def make_case(root: Path, scenario="largeY_largeL", model="T3_dS1_3_dS2_3_dF_2_alpha_p2", cls="T3_E_alpha_p2") -> Path:
    path = root / scenario / model / cls / "data" / "running_diagnostics.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"status":"Success","ordering":"NO","low_energy":{
        "masses_ev":[0.01,0.02,0.05],"delta_m21_sq_ev2":0.0003,
        "delta_m31_sq_ev2":0.0024,"delta_m32_sq_ev2":0.0021,
        "pmns_abs":[[0.8,0.5,0.1],[0.4,0.6,0.7],[0.3,0.6,0.7]],
        "takagi_residual":1e-14}}),encoding="utf-8")
    return path


def test_extraction_and_angles(tmp_path):
    path=make_case(tmp_path)
    result=extract_case(path,tmp_path)
    assert result["model_class"] == "E" and result["alpha"] == 2
    assert abs(result["sin2_theta13"]-0.01)<1e-14
    assert abs(result["sin2_theta12"]-0.25/0.99)<1e-14
    assert abs(result["sin2_theta23"]-0.49/0.99)<1e-14


def test_malformed_diagnostic_recorded_without_stopping(tmp_path):
    path=make_case(tmp_path)
    bad=make_case(tmp_path,scenario="smallY_smallL")
    bad.write_text("not json",encoding="utf-8")
    rows,errors=collect(tmp_path)
    assert len(rows)==1 and len(errors)==1
    assert str(bad)==errors[0]["path"]


def test_writes_csv_json_and_incomplete_status(tmp_path):
    input_root=tmp_path/"input"
    make_case(input_root)
    out=tmp_path/"out"
    summary=write_results(input_root,out)
    assert summary["status"]=="Incomplete" and summary["extracted_cases"]==1
    assert (out/"t3_full_comparison_observables.csv").is_file()
    assert len(json.loads((out/"t3_full_comparison_observables.json").read_text())["rows"])==1


def test_missing_or_failed_case_does_not_inflate_count(tmp_path):
    path=make_case(tmp_path)
    data=json.loads(path.read_text());data["status"]="Failed";path.write_text(json.dumps(data))
    rows,errors=collect(tmp_path)
    assert rows==[] and len(errors)==1
