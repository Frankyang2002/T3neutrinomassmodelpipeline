from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
from Numerical.plotting import WithinModelComparison as w

def _payload(model:str,dims=(1,3,2))->dict:
    mu=[1000.0,316.0,100.0]; c=np.ones((3,3,3)).tolist()
    return {"status":"Success","benchmark_search":{"model":{
        "model_key":model,"d_s1":dims[0],"d_s2":dims[1],"d_f":dims[2]}},
        "final_running":{"mu_gev":mu,"c5_abs":c,
        "delta_m21_sq_ev2":[7e-5,7.2e-5,7.4e-5],
        "delta_m3l_sq_ev2":[2.4e-3,2.45e-3,2.5e-3],
        "sin2_theta12":[.30,.31,.32],"sin2_theta13":[.021,.022,.023],
        "sin2_theta23":[.54,.55,.56]},
        "intermediate_direct_weinberg":{"mu_gev":[1e5,1e4,1e3],"delta_c5_abs":c}}

def test_collect_and_build_four_scenarios(tmp_path:Path)->None:
    model="T3_dS1_1_dS2_3_dF_2_alpha_p0"; source=tmp_path/"comparison"
    for scenario in w.SCENARIO_ORDER:
        p=source/scenario/model/"data"; p.mkdir(parents=True)
        (p/"running_diagnostics.json").write_text(json.dumps(_payload(model)),encoding="utf-8")
    collected=w.collect(source); assert set(collected)=={model}; assert tuple(collected[model])==w.SCENARIO_ORDER
    output=tmp_path/"reports"; w.build_within_model_comparisons(source,output); model_dir=output/model
    for name in ("direct_weinberg_norm_comparison.png","c5_norm_comparison.png",
                 "dm21_comparison.png","dm3l_comparison.png","mixing_comparison.png"):
        assert (model_dir/name).is_file()
    text=(model_dir/"README.md").read_text(encoding="utf-8")
    assert "lambda_T3^ref=1" in text
    assert "lambda_T3^eff=1" in text

def test_class_c_label_uses_effective_lambda()->None:
    model="T3_dS1_2_dS2_2_dF_3_alpha_m1"; payload=_payload(model,(2,2,3))
    ref,eff,factor,cls=w._lambda_values(model,"smallY_smallL",payload)
    assert cls=="C" and ref==0.1
    assert math.isclose(factor,math.sqrt(3.0))
    assert math.isclose(eff,0.1*math.sqrt(3.0))
    label=w._label(model,"smallY_smallL",payload)
    assert "ref" in label and "eff" in label

def test_scenario_values_are_current()->None:
    assert w.SCENARIOS=={"smallY_smallL":(0.005,0.1),"smallY_largeL":(0.005,1.0),
                         "largeY_smallL":(0.5,0.1),"largeY_largeL":(0.5,1.0)}


def test_class_encoded_model_names_without_dimension_metadata() -> None:
    # The actual regression: production diagnostics can use T3_A_alpha_m2
    # without recording d_s1/d_s2/d_f in benchmark_search.model.
    for model, expected_class, expected_factor in (
        ("T3_A_alpha_m2", "A", 1.0),
        ("T3_B_alpha_m1", "B", 1.0),
        ("T3_C_alpha_p1", "C", math.sqrt(3.0)),
        ("T3_D_alpha_p0", "D", 1.0),
        ("T3_E_alpha_m2", "E", 1.0 / math.sqrt(2.0)),
    ):
        payload = {"benchmark_search": {"model": {"model_key": model}}}
        ref, effective, factor, cls = w._lambda_values(
            model, "smallY_smallL", payload
        )
        assert cls == expected_class
        assert math.isclose(factor, expected_factor)
        assert math.isclose(effective, ref * expected_factor)


def test_build_report_from_class_encoded_model_key(tmp_path: Path) -> None:
    # Test the real path through collect -> plotting -> README for the
    # no-dimensions diagnostics responsible for the observed traceback.
    model = "T3_A_alpha_m2"
    source = tmp_path / "comparison"
    for scenario in w.SCENARIO_ORDER:
        path = source / scenario / model / "data" / "running_diagnostics.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = _payload(model)
        payload["benchmark_search"]["model"] = {"model_key": model}
        path.write_text(json.dumps(payload), encoding="utf-8")

    output = tmp_path / "reports"
    w.build_within_model_comparisons(source, output)
    report = output / model
    for figure in (
        "direct_weinberg_norm_comparison.png",
        "c5_norm_comparison.png",
        "dm21_comparison.png",
        "dm3l_comparison.png",
        "mixing_comparison.png",
    ):
        assert (report / figure).is_file()
    assert "class=A" in (report / "README.md").read_text(encoding="utf-8")
