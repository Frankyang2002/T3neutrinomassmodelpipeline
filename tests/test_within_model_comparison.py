from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from Numerical.plotting import WithinModelComparison as w


def _payload(model: str) -> dict:
    mu=[1000.0,316.0,100.0]
    c=np.ones((3,3,3)).tolist()
    direct=np.ones((3,3,3)).tolist()
    return {
        "status":"Success",
        "benchmark_search":{"model":{"model_key":model}},
        "final_running":{
            "mu_gev":mu,
            "c5_abs":c,
            "delta_m21_sq_ev2":[7e-5,7.2e-5,7.4e-5],
            "delta_m3l_sq_ev2":[2.4e-3,2.45e-3,2.5e-3],
            "sin2_theta12":[.30,.31,.32],
            "sin2_theta13":[.021,.022,.023],
            "sin2_theta23":[.54,.55,.56],
        },
        "intermediate_direct_weinberg":{
            "mu_gev":[1e5,1e4,1e3],
            "delta_c5_abs":direct,
        },
    }


def test_collect_and_build_four_scenarios(tmp_path: Path) -> None:
    model="T3_dS1_1_dS2_3_dF_2_alpha_p0"
    source=tmp_path/"comparison"
    for scenario in w.SCENARIO_ORDER:
        p=source/scenario/model/"data"
        p.mkdir(parents=True)
        (p/"running_diagnostics.json").write_text(
            json.dumps(_payload(model)),encoding="utf-8"
        )

    collected=w.collect(source)
    assert set(collected)=={model}
    assert tuple(collected[model])==w.SCENARIO_ORDER

    output=tmp_path/"reports"
    w.build_within_model_comparisons(source,output)
    model_dir=output/model
    assert (model_dir/"direct_weinberg_norm_comparison.png").is_file()
    assert (model_dir/"c5_norm_comparison.png").is_file()
    assert (model_dir/"dm21_comparison.png").is_file()
    assert (model_dir/"dm3l_comparison.png").is_file()
    assert (model_dir/"mixing_comparison.png").is_file()
    text=(model_dir/"README.md").read_text(encoding="utf-8")
    assert "No physics was recomputed" in text
    assert "lambda_T3=1" in text


def test_scenario_values_are_current() -> None:
    assert w.SCENARIOS == {
        "smallY_smallL":(0.005,0.1),
        "smallY_largeL":(0.005,1.0),
        "largeY_smallL":(0.5,0.1),
        "largeY_largeL":(0.5,1.0),
    }
