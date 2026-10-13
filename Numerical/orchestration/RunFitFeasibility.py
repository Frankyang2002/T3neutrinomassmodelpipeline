"""One-case NuFIT feasibility test: Yukawa scale calibration + 256-point flavour scan.

Does not modify fixed comparison configs, running diagnostics or benchmark cache.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import qmc
from scipy.optimize import least_squares

from Numerical.fitting.SobolBenchmarkSearch import (
    BenchmarkParameter, _build_evaluator, apply_benchmark_parameters,
)
from Numerical.orchestration.ResumePipelineNumericalResults import (
    PROJECT_ROOT, OUTPUT_ROOT, CONFIG_ROOT, REQUIRED,
    discover_cases, config_for_case, inspect_case,
)


def prediction(fit):
    return {str(k): float(v) for k,v in zip(fit.observable_names, fit.prediction)}


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--scenario", default="smallY_smallL")
    p.add_argument("--model", default="T3_dS1_1_dS2_3_dF_2_alpha_m2")
    p.add_argument("--points", type=int, default=256)
    p.add_argument("--local-evals", type=int, default=30)
    p.add_argument("--run", action="store_true")
    a=p.parse_args(argv)
    if a.points<2 or a.points & (a.points-1):
        p.error("--points must be a power of 2")
    matches=[(s,k,c) for s,k,c in discover_cases(OUTPUT_ROOT)
             if s==a.scenario and k==a.model]
    if len(matches)!=1:
        p.error(f"Expected one matching case, got {len(matches)}")
    sc,key,case=matches[0]
    check=inspect_case(sc,key,case,CONFIG_ROOT)
    if check["missing_artifacts"] or check["config_error"]:
        p.error(str(check["missing_artifacts"] or check["config_error"]))
    source=config_for_case(sc,key,CONFIG_ROOT)
    original=json.loads(source.read_text(encoding="utf-8-sig"))
    ordering=str(original.get("ordering","NO")).upper()
    if ordering not in ("NO","IO"):
        p.error("Explicit NO/IO ordering required")
    target=PROJECT_ROOT/"configs"/f"nufit_6_1_ic24_{ordering.lower()}.json"
    if not target.is_file():
        p.error(f"Missing {target}")
    data=case/"data"
    artifacts=[data/name for name in REQUIRED]
    out=PROJECT_ROOT/"Reports"/"output"/"full"/"fit_feasibility"/sc/key
    print("Case:",case,flush=True)
    print("NuFIT target:",target,flush=True)
    print("Output:",out,flush=True)
    if not a.run:
        print("Read-only dry run; add --run.",flush=True)
        return 0

    def evaluator(payload,definitions=()):
        return _build_evaluator(
            payload=payload, parameters=definitions,
            uv_rgbeta_path=artifacts[0],
            intermediate_rgbeta_path=artifacts[1],
            final_weinberg_path=artifacts[2],
            target_path=target,
        )

    def scale_yukawas(factor):
        payload=deepcopy(original)
        for y in ("y1","y2"):
            spec=payload["base_state"]["ordinary"][y]
            for part in ("real","imag"):
                if part in spec:
                    spec[part]=(np.asarray(spec[part],dtype=float)*factor).tolist()
        return payload

    # Logarithmic mass-scale calibration. Includes unscaled benchmark.
    calibration=[]
    for exponent in np.linspace(-7.0,0.0,15):
        factor=float(10.0**exponent)
        try:
            fit=evaluator(scale_yukawas(factor))({})
            item={"factor":factor,"chi2":float(fit.chi2),
                  "prediction":prediction(fit),"status":"Success"}
        except Exception as exc:
            item={"factor":factor,"status":"Failed","error":str(exc)}
        calibration.append(item)
        print(f"calibration factor={factor:.4g} chi2={item.get('chi2','Failed')}",flush=True)
    valid=[v for v in calibration if v["status"]=="Success" and math.isfinite(v["chi2"])]
    if not valid:
        raise RuntimeError("All calibration evaluations failed")
    best_scale=min(valid,key=lambda v:v["chi2"])
    base=scale_yukawas(best_scale["factor"])

    # Fit all 18 real Yukawa entries, around the calibrated mass scale.
    # Per-entry linear bounds retain sign freedom and allow texture changes.
    defs=[]
    for y in ("y1","y2"):
        mat=np.asarray(base["base_state"]["ordinary"][y]["real"],dtype=float)
        amplitude=max(float(np.max(np.abs(mat))),1e-16)
        for i in range(3):
            for j in range(3):
                defs.append(BenchmarkParameter(
                    name=f"{y}_{i}{j}",
                    binding=f"ordinary.{y}.real.{i}.{j}",
                    low=-3.0*amplitude,high=3.0*amplitude,
                ))
    fn=evaluator(base,defs)
    sampler=qmc.Sobol(d=len(defs),scramble=True,seed=20260927)
    sobol=sampler.random_base2(m=int(math.log2(a.points)))
    candidates=[]
    def assess(values,phase):
        try:
            fit=fn(values)
            chi2=float(fit.chi2)
            if not math.isfinite(chi2):
                raise ValueError("Nonfinite chi2")
            return {"phase":phase,"chi2":chi2,
                    "prediction":prediction(fit),
                    "parameters":{k:float(v) for k,v in values.items()},
                    "status":"Success"}
        except Exception as exc:
            return {"phase":phase,"status":"Failed","error":str(exc)}
    # Preserve the calibrated original flavour texture as a candidate.
    initial={d.name:float(np.asarray(base["base_state"]["ordinary"]
        [d.name.split("_")[0]]["real"])[int(d.name[-2]),int(d.name[-1])]) for d in defs}
    candidates.append(assess(initial,"calibrated_baseline"))
    for idx,unit in enumerate(sobol):
        values={d.name:float(d.low+(d.high-d.low)*z) for d,z in zip(defs,unit)}
        candidates.append(assess(values,"sobol"))
        if (idx+1)%32==0:
            successful=[x for x in candidates if x["status"]=="Success"]
            best=min(successful,key=lambda x:x["chi2"]) if successful else None
            print(f"Sobol {idx+1}/{a.points}, best chi2={best['chi2'] if best else 'none'}",flush=True)
    successful=sorted((x for x in candidates if x["status"]=="Success"),
                      key=lambda x:x["chi2"])
    if not successful:
        raise RuntimeError("No valid flavour evaluations")
    best=successful[0]

    # Limited local improvement on six most sensitive coordinates.
    # Residuals are whitened by the repository's NuFIT covariance.
    local_info={"status":"Skipped"}
    if a.local_evals>0:
        from Numerical.fitting.OscillationOptimizer import estimate_local_sensitivities
        # Avoid private optimizer assumptions: rank coordinates by finite
        # difference of chi2, then fit six coordinates with least_squares.
        x0=np.array([best["parameters"][d.name] for d in defs])
        ranks=[]
        for i,d in enumerate(defs):
            trial=x0.copy()
            trial[i]=min(d.high,max(d.low,trial[i]+0.01*(d.high-d.low)))
            vals={v.name:float(z) for v,z in zip(defs,trial)}
            result=assess(vals,"sensitivity_probe")
            ranks.append((abs(result.get("chi2",best["chi2"])-best["chi2"]),i))
        active=[i for _,i in sorted(ranks,reverse=True)[:6]]
        xbase=x0.copy()
        selected=[defs[i] for i in active]
        low=np.array([d.low for d in selected])
        high=np.array([d.high for d in selected])
        def residual(z):
            full=xbase.copy()
            full[active]=z
            values={d.name:float(v) for d,v in zip(defs,full)}
            try:
                fit=fn(values)
                return np.asarray(fit.whitened_residual,dtype=float)
            except Exception:
                # Failure is penalised, not treated as a valid zero residual.
                return np.full(len(best["prediction"]),1e10)
        try:
            result=least_squares(residual,xbase[active],bounds=(low,high),
                                 max_nfev=a.local_evals)
            xbase[active]=result.x
            improved=assess({d.name:float(v) for d,v in zip(defs,xbase)},"local")
            local_info={"status":improved["status"],"nfev":int(result.nfev),
                        "active_parameters":[defs[i].name for i in active],
                        "chi2":improved.get("chi2")}
            if improved["status"]=="Success" and improved["chi2"]<best["chi2"]:
                best=improved
        except Exception as exc:
            local_info={"status":"Failed","error":str(exc)}
    fitted=apply_benchmark_parameters(base,best["parameters"],defs)
    fitted["benchmark_search"]={"enabled":False}
    fitted["sensitivity"]={"enabled":False}
    summary={"scenario":sc,"model":key,"target":str(target),
             "calibration":calibration,"best_scale":best_scale,
             "flavour_points":a.points,"successful_flavour_points":len(successful),
             "local":local_info,"best":best}
    out.mkdir(parents=True,exist_ok=True)
    (out/"feasibility_summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    (out/"fitted_candidate_config.json").write_text(json.dumps(fitted,indent=2)+"\n",encoding="utf-8")
    print("Best chi2:",best["chi2"],flush=True)
    print("Best prediction:",json.dumps(best["prediction"],indent=2),flush=True)
    print("Saved:",out,flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
