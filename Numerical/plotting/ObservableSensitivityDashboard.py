"""Interactive observable-plane T3 sensitivity dashboards from saved diagnostics.

Reads separate `sensitivity_diagnostics.json` `directions.{T1,T2,T3}.points` and NuFIT 6.1
configuration files. Does not evaluate RGEs, matching or missing scan points.

Run:
  python -m Numerical.plotting.ObservableSensitivityDashboard
  python -m Numerical.plotting.ObservableSensitivityDashboard --help
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[2]
SCENARIOS = ("smallY_smallL", "smallY_largeL", "largeY_smallL", "largeY_largeL")
DIRECTIONS = ("T1", "T2", "T3")
OBS = {
    "sin2_theta12": ("sin²θ₁₂", 1.0),
    "sin2_theta13": ("sin²θ₁₃", 1.0),
    "sin2_theta23": ("sin²θ₂₃", 1.0),
    "delta_m21_sq_ev2": ("Δm²₂₁ [10⁻⁵ eV²]", 1e5),
    "delta_m3l_sq_ev2": ("Δm²₃ℓ [10⁻³ eV²]", 1e3),
}
DEFAULT_PAIRS = (
    ("sin2_theta13", "delta_m3l_sq_ev2"),
    ("sin2_theta23", "delta_m3l_sq_ev2"),
    ("sin2_theta13", "sin2_theta23"),
)
COLORS = {"T1": "#1f77b4", "T2": "#ff7f0e", "T3": "#2ca02c"}
SCENARIO_COLORS = {
    "smallY_smallL": "#1f77b4",
    "smallY_largeL": "#ff7f0e",
    "largeY_smallL": "#2ca02c",
    "largeY_largeL": "#d62728",
}


def load_nufit(ordering: str, root: Path = ROOT) -> dict:
    name = f"nufit_6_1_ic24_{ordering.lower()}.json"
    payload = json.loads((root / "configs" / name).read_text(encoding="utf-8"))
    if payload.get("ordering") != ordering:
        raise ValueError(f"NuFIT ordering mismatch in {name}")
    for key in OBS:
        if key not in payload["central_values"] or key not in payload["one_sigma_errors"]:
            raise ValueError(f"Missing NuFIT observable {key}")
    return payload


def _finite(value, context: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite {context}")
    return result


def _normalise_prediction(pred: dict, ordering: str) -> dict:
    """Map aliases only when mathematically justified by the ordering.

    NuFIT defines dm3l = dm31 for NO, dm32 for IO.
    """
    result = {}
    for key in OBS:
        if key in pred:
            result[key] = _finite(pred[key], key)
    if "delta_m3l_sq_ev2" not in result:
        if ordering == "NO" and "delta_m31_sq_ev2" in pred:
            result["delta_m3l_sq_ev2"] = _finite(pred["delta_m31_sq_ev2"], "dm31")
        elif ordering == "IO" and "delta_m32_sq_ev2" in pred:
            result["delta_m3l_sq_ev2"] = _finite(pred["delta_m32_sq_ev2"], "dm32")
    return result


def load_saved(root: Path) -> tuple[list[dict], list[str]]:
    """Read valid sensitivity trajectories; skip absent or malformed data."""
    records, warnings = [], []
    paths = sorted(root.glob("*/*/T3_*/data/sensitivity_diagnostics.json"))
    for path in paths:
        rel = path.relative_to(root)
        scenario, model, class_key = rel.parts[:3]
        if scenario not in SCENARIOS:
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            sensitivity = payload
            if sensitivity.get("status") != "Success":
                warnings.append(f"{scenario}/{model}: sensitivity status is not Success")
                continue
            ordering = str(payload.get("ordering", sensitivity.get("target", {}).get("ordering", "")))
            if ordering not in {"NO", "IO"}:
                raise ValueError(f"unknown ordering {ordering!r}")
            directions = sensitivity.get("directions", {})
            for direction in DIRECTIONS:
                section = directions.get(direction, {})
                points = section.get("points", []) if isinstance(section, dict) else []
                valid = []
                for point in points:
                    if not isinstance(point, dict) or not isinstance(point.get("prediction"), dict):
                        continue
                    prediction = _normalise_prediction(point["prediction"], ordering)
                    if not prediction:
                        continue
                    valid.append({
                        "gamma": _finite(point["gamma"], "gamma"),
                        "chi2": _finite(point["chi2"], "chi2") if point.get("chi2") is not None else None,
                        "prediction": prediction,
                    })
                valid.sort(key=lambda p: p["gamma"])
                if valid:
                    records.append({
                        "scenario": scenario, "model": model, "class_key": class_key,
                        "ordering": ordering, "direction": direction, "points": valid,
                    })
        except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as exc:
            warnings.append(f"{rel}: {exc}")
    return records, warnings


def build_html(records: list[dict], nufit: dict, mode: str, title: str) -> str:
    """Self-contained page logic; Plotly JS is loaded from its public CDN."""
    if mode not in {"within", "between"}:
        raise ValueError(mode)
    import html
    import json as _json
    dataset = _json.dumps(records, ensure_ascii=False, allow_nan=False).replace("</", "<\\/")
    cfg = _json.dumps(nufit, ensure_ascii=False, allow_nan=False).replace("</", "<\\/")
    obs = _json.dumps(OBS, ensure_ascii=False)
    scenario_colors = _json.dumps(SCENARIO_COLORS)
    direction_colors = _json.dumps(COLORS)
    heading = html.escape(title)
    # The bands are independent 1D Gaussianized intervals, not a 2D NuFIT likelihood.
    return r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
body{font-family:system-ui,Segoe UI,Arial,sans-serif;margin:24px;max-width:1500px;color:#1d2935}
h1{font-size:1.45rem;margin-bottom:4px}.note{color:#52616e;font-size:.86rem;line-height:1.5}
.controls{display:flex;gap:14px;flex-wrap:wrap;margin:20px 0}
label{font-size:.85rem;display:flex;flex-direction:column;gap:5px}
select{padding:7px 9px;border:1px solid #bcc8d2;border-radius:5px;min-width:160px}
#plots{display:grid;grid-template-columns:repeat(auto-fit,minmax(390px,1fr));gap:14px}
.plot{height:410px;border:1px solid #e4e8ed;border-radius:6px}
#status{margin:14px 0;font-size:.9rem}
</style></head><body>
<h1>__TITLE__</h1>
<div class="note">NuFIT 6.1 IC24, with SK-atm. Bands show independent Gaussianized ±1σ / ±3σ
intervals from the repository's NuFIT configuration; they are <b>not</b> official
two-dimensional NuFIT contours. Only saved sensitivity scan points are drawn.
Plotly requires an internet connection for the JavaScript library.</div>
<div class="controls">
<label>Ordering<select id="ordering"></select></label>
<label>__FILTER_LABEL__<select id="filter"></select></label>
<label>__GROUP_LABEL__<select id="group"></select></label>
<label>Direction<select id="direction"></select></label>
<label>Horizontal observable<select id="x"></select></label>
<label>Vertical observable<select id="y"></select></label>
</div>
<div id="status"></div><div id="plots"></div>
<script>
const data=__DATA__, nufit=__NUFIT__, obs=__OBS__, scColors=__SC_COLORS__, dirColors=__DIR_COLORS__;
const modelPalette=["#1f77b4","#d62728","#2ca02c","#9467bd","#ff7f0e","#17becf","#8c564b","#e377c2","#7f7f7f","#bcbd22"];
const modelColors=Object.fromEntries(uniqueModels().map((m,i)=>[m,modelPalette[i%modelPalette.length]]));
function uniqueModels(){return [...new Set(data.map(r=>r.model))].sort()}
const mode="__MODE__", ids=["ordering","filter","group","direction","x","y"];
const el=id=>document.getElementById(id);
const unique=a=>[...new Set(a)].sort();
function options(id,values,preferred){const e=el(id), old=e.value;e.replaceChildren();
  values.forEach(v=>{const o=document.createElement("option");o.value=v;o.textContent=v;e.appendChild(o)});
  e.value=values.includes(old)?old:(values.includes(preferred)?preferred:values[0]);}
const obsKeys=Object.keys(obs);
options("x",obsKeys,"sin2_theta13");options("y",obsKeys,"delta_m3l_sq_ev2");
options("ordering",unique(data.map(r=>r.ordering)),"NO");
function updateControls(){
 const ordering=el("ordering").value;
 const ordered=data.filter(r=>r.ordering===ordering);
 if(mode==="within"){
   options("filter",unique(ordered.map(r=>r.model)));
   options("group",["All scenarios",...unique(ordered.filter(r=>r.model===el("filter").value).map(r=>r.scenario))],"All scenarios");
 }else{
   options("filter",unique(ordered.map(r=>r.scenario)));
   options("group",["All models",...unique(ordered.filter(r=>r.scenario===el("filter").value).map(r=>r.model))],"All models");
 }
 options("direction",["All directions","T1","T2","T3"],"All directions");
}
function bands(x,y){
 const c=nufit.central_values,s=nufit.one_sigma_errors;
 if(!(x in c)||!(y in c))return [];
 const shapes=[];
 for(const k of [3,1]){
  const fill=k===3?"rgba(55,160,240,.13)":"rgba(20,135,230,.26)";
  shapes.push({type:"rect",xref:"x",yref:"paper",x0:(c[x]-k*s[x])*obs[x][1],x1:(c[x]+k*s[x])*obs[x][1],y0:0,y1:1,fillcolor:fill,line:{width:0},layer:"below"});
  shapes.push({type:"rect",xref:"paper",yref:"y",x0:0,x1:1,y0:(c[y]-k*s[y])*obs[y][1],y1:(c[y]+k*s[y])*obs[y][1],fillcolor:fill,line:{width:0},layer:"below"});
 }
 return shapes;
}
function draw(){
 if(typeof Plotly==="undefined"){el("status").textContent="Plotly could not load. Connect to the internet and reopen this HTML.";return;}
 const ordering=el("ordering").value,filter=el("filter").value,group=el("group").value,dir=el("direction").value;
 const x=el("x").value,y=el("y").value;
 let rows=data.filter(r=>r.ordering===ordering && (mode==="within"?r.model===filter:r.scenario===filter));
 if(mode==="within"&&group!=="All scenarios")rows=rows.filter(r=>r.scenario===group);
 if(mode==="between"&&group!=="All models")rows=rows.filter(r=>r.model===group);
 if(dir!=="All directions")rows=rows.filter(r=>r.direction===dir);
 const traces=[];
 rows.forEach(r=>{
  const points=r.points.filter(p=>x in p.prediction && y in p.prediction);
  if(!points.length)return;
  const label=(mode==="within"?r.scenario:r.model)+" · "+r.direction;
  const color=mode==="within"?scColors[r.scenario]:modelColors[r.model];
  traces.push({type:"scatter",mode:"lines+markers",name:label,
    x:points.map(p=>p.prediction[x]*obs[x][1]),y:points.map(p=>p.prediction[y]*obs[y][1]),
    customdata:points.map(p=>[p.gamma,p.chi2,r.model,r.scenario]),
    hovertemplate:"γ=%{customdata[0]:.5f}<br>χ²=%{customdata[1]}<br>Model=%{customdata[2]}<br>Scenario=%{customdata[3]}<br>x=%{x:.6g}<br>y=%{y:.6g}<extra>%{fullData.name}</extra>",
    line:{color:color,width:2,dash:r.direction==="T2"?"dash":r.direction==="T3"?"dot":"solid"},
    marker:{size:6,color:color}});
 });
 const plots=el("plots");plots.replaceChildren();
 const div=document.createElement("div");div.className="plot";plots.appendChild(div);
 const label=obs[x][0]+" vs "+obs[y][0];
 Plotly.newPlot(div,traces,{
  title:{text:label,font:{size:14}},xaxis:{title:obs[x][0],zeroline:false,showgrid:true},
  yaxis:{title:obs[y][0],zeroline:false,showgrid:true},
  shapes:bands(x,y),hovermode:"closest",legend:{orientation:"h",y:-.28},
  margin:{l:78,r:20,t:60,b:110},paper_bgcolor:"#fff",plot_bgcolor:"#fff"
 },{responsive:true,displaylogo:false,toImageButtonOptions:{format:"png",scale:3}});
 el("status").textContent=traces.length+" saved trajectories shown; "+rows.length+" candidate records. "+
 (traces.length?"Zoom or hover to inspect points.":"No saved sensitivity points for this selection.");
}
ids.forEach(id=>el(id).addEventListener("change",()=>{
 if(id==="ordering"||id==="filter")updateControls();
 draw();
}));
updateControls();draw();
</script></body></html>""".replace("__TITLE__",heading).replace("__DATA__",dataset).replace("__NUFIT__",cfg).replace("__OBS__",obs).replace("__SC_COLORS__",scenario_colors).replace("__DIR_COLORS__",direction_colors).replace("__MODE__",mode).replace("__FILTER_LABEL__","Model" if mode=="within" else "Scenario").replace("__GROUP_LABEL__","Scenario" if mode=="within" else "Model")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=ROOT / "output" / "full" / "comparison")
    parser.add_argument("--output-root", type=Path, default=ROOT / "Reports" / "output" / "full" / "observable_sensitivity")
    args = parser.parse_args(argv)
    if not args.input_root.is_dir():
        parser.error(f"Missing input directory: {args.input_root}")
    records, warnings = load_saved(args.input_root)
    args.output_root.mkdir(parents=True, exist_ok=True)
    for ordering in ("NO", "IO"):
        cfg = load_nufit(ordering)
        subset = [r for r in records if r["ordering"] == ordering]
        if not subset:
            continue
        for mode in ("within", "between"):
            dest = args.output_root / f"{mode}_model_{ordering.lower()}.html"
            dest.write_text(build_html(subset, cfg, mode, f"T3 {mode}-model sensitivity ({ordering})"), encoding="utf-8")
            print(f"Created {dest}")
    manifest = {
        "sensitivity_records": len(records),
        "models": sorted({r["model"] for r in records}),
        "scenarios": sorted({r["scenario"] for r in records}),
        "warnings": warnings,
        "note": "No numerical points were generated. Only saved successful sensitivity trajectories are plotted."
    }
    (args.output_root / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n",encoding="utf-8")
    print(f"Found {len(records)} saved sensitivity direction records; {len(warnings)} missing/invalid diagnostics")
    if not records:
        print("No saved sensitivity trajectories found. The pipeline must first produce sensitivity-enabled numerical diagnostics.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
