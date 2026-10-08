"""Interactive cross-model T3 running dashboard from saved numerical diagnostics.

Presentation-only: neither RGEs nor matching are recomputed.

Weinberg views show a threshold-anchored continuation |C5_hard(MS) +
Delta_C5_direct(mu)| above MS, joined continuously to the stored final SM+
Weinberg C5 below MS.  This continuation is a display reconstruction, not an
intermediate-EFT Wilson coefficient.  Optional dotted overlays reveal each
component separately; the hard term is held fixed at its MS value.

Older benchmark outputs save only absolute intermediate matrices. The separate
WeinbergMatchedContinuation module reconstructs their signed real one-loop
trajectory only after stringent checks. Complex ambiguous cases are skipped.
No artificial complex phase or running is introduced.

An explicit LLSS selector draws only the universal fixed-order one-loop kernel
eta(mu)=ln(mu/MF)/(16*pi**2).  It is not a numerical LLSS Wilson coefficient:
the saved benchmark diagnostics do not contain its full-flavor trajectory.
"""
from __future__ import annotations

import argparse
import html
import json
import math
import shutil
from pathlib import Path
from typing import Any

import numpy as np

from Numerical.plotting.WeinbergMatchedContinuation import (
    reconstruct_matched_continuation,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
COMPARISON_SCENARIOS = (
    "smallY_smallL", "smallY_largeL", "largeY_smallL", "largeY_largeL"
)

# Each tuple: display label, diagnostics block, array key, unit, selector.
QUANTITIES: dict[str, tuple[Any, ...]] = {
    "dm21": ("Solar splitting Δm²21", "final_running", "delta_m21_sq_ev2", "eV²", False),
    "dm3l": ("Atmospheric splitting |Δm²3ℓ|", "final_running", "delta_m3l_sq_ev2", "eV²", True),
    "s12": ("sin²θ12", "final_running", "sin2_theta12", "", False),
    "s13": ("sin²θ13", "final_running", "sin2_theta13", "", False),
    "s23": ("sin²θ23", "final_running", "sin2_theta23", "", False),
    "m1": ("m1", "final_running", "masses_ev", "eV", 0),
    "m2": ("m2", "final_running", "masses_ev", "eV", 1),
    "m3": ("m3", "final_running", "masses_ev", "eV", 2),
    "y1norm": ("||y1||F", "uv_running", "y1_frobenius_norm", "", False),
    "y2norm": ("||y2||F", "uv_running", "y2_frobenius_norm", "", False),
    "lambdaT3": ("|λT3|", "uv_running", "lambdaT3_abs", "", False),
}

for i in range(3):
    for j in range(i, 3):
        suffix = f"{i+1}{j+1}"
        QUANTITIES[f"c5_{suffix}"] = (
            f"|C5^{{{suffix}}}| — hard + direct continuation", "final_running",
            "c5_abs", "GeV⁻¹", (i, j),
        )
        QUANTITIES[f"direct_c5_{suffix}"] = (
            f"|ΔC5^{{{suffix}}}| — direct intermediate", "intermediate_direct_weinberg",
            "delta_c5_abs", "GeV⁻¹", (i, j),
        )
QUANTITIES["c5_norm"] = (
    "||C5||F — hard + direct continuation", "final_running", "c5_abs",
    "GeV⁻¹", "frobenius",
)
QUANTITIES["direct_c5_norm"] = (
    "||ΔC5||F — direct intermediate", "intermediate_direct_weinberg",
    "delta_c5_abs", "GeV⁻¹", "frobenius",
)


def _extract(payload: dict[str, Any], spec: tuple[Any, ...]) -> tuple[list[float], list[float]]:
    """Select a physical diagnostic; never infer a missing trajectory."""
    _, block_name, key, _, selector = spec
    block = payload[block_name]
    x = [float(v) for v in block["mu_gev"]]
    raw = block[key]
    if isinstance(selector, tuple):
        i, j = selector
        y = [float(row[i][j]) for row in raw]
    elif selector == "frobenius":
        # The saved matrices contain absolute values.  Their squared entries
        # give the Frobenius norm exactly, including off-diagonal duplicates.
        y = [math.sqrt(sum(float(v) ** 2 for row in matrix for v in row))
             for matrix in raw]
    elif isinstance(selector, int) and not isinstance(selector, bool):
        y = [float(row[selector]) for row in raw]
    else:
        y = [float(v) for v in raw]
        if selector is True:
            y = [abs(v) for v in y]
    if len(x) != len(y) or not x:
        raise ValueError(f"Invalid saved curve shape: {block_name}.{key}")
    if not all(math.isfinite(v) and v > 0 for v in x):
        raise ValueError("Renormalisation scales must be positive and finite")
    if not all(math.isfinite(v) for v in y):
        raise ValueError("Saved curve contains non-finite values")
    return x, y


def _model_name(path: Path, payload: dict[str, Any]) -> str:
    model_dir = path.parent.parent
    search = payload.get("benchmark_search")
    rep = search.get("model") if isinstance(search, dict) else None
    if isinstance(rep, dict) and rep.get("model_key"):
        return str(rep["model_key"])
    return model_dir.name


def _is_canonical_optimal_diagnostic(study_dir: Path, path: Path) -> bool:
    """Ignore short aliases in the older optimal comparison directory."""
    study_dir = Path(study_dir)
    if study_dir.name != "optimal":
        return True
    try:
        relative = path.relative_to(study_dir)
    except ValueError:
        return False
    return bool(relative.parts and relative.parts[0].startswith("T3_dS1_"))


def collect(study_dir: Path) -> dict[str, Any]:
    """Read stored curves and build validated matched-continuation views."""
    result: dict[str, Any] = {
        "quantities": {}, "overlays": {}, "matched": {}, "components": {},
        "continuation_methods": {}, "continuation_issues": {},
        "models": [], "scales": {},
    }
    study_dir = Path(study_dir)
    files = [p for p in sorted(study_dir.rglob("running_diagnostics.json"))
             if _is_canonical_optimal_diagnostic(study_dir, p)]
    for path in files:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if not isinstance(payload, dict) or payload.get("status") != "Success":
            continue
        name = _model_name(path, payload)
        result["models"].append(name)
        scales = payload.get("scales_gev", {})
        if isinstance(scales, dict) and not result["scales"]:
            result["scales"] = {
                str(k): float(v) for k, v in scales.items()
                if isinstance(v, (int, float)) and math.isfinite(v) and v > 0
            }
        for qkey, spec in QUANTITIES.items():
            try:
                x, y = _extract(payload, spec)
            except (KeyError, TypeError, ValueError, IndexError):
                continue
            result["quantities"].setdefault(qkey, {})[name] = {"x": x, "y": y}

        # Overlay the existing, independently sampled direct LLSS -> C5 mixing.
        # This is numerically physical; unlike LLSS self-running, it *does*
        # exist as a full-flavor matrix in the saved benchmark diagnostics.
        for i in range(3):
            for j in range(i, 3):
                suffix = f"{i+1}{j+1}"
                qkey = f"c5_{suffix}"
                if name not in result["quantities"].get(qkey, {}):
                    continue
                try:
                    x, y = _extract(payload, QUANTITIES[f"direct_c5_{suffix}"])
                except (KeyError, TypeError, ValueError, IndexError):
                    continue
                result["overlays"].setdefault(qkey, {})[name] = {"x": x, "y": y}
        if name in result["quantities"].get("c5_norm", {}):
            try:
                x, y = _extract(payload, QUANTITIES["direct_c5_norm"])
            except (KeyError, TypeError, ValueError, IndexError):
                pass
            else:
                result["overlays"].setdefault("c5_norm", {})[name] = {"x": x, "y": y}

        # The matched continuation must add COMPLEX coefficients before
        # applying absolute values; summing |hard| + |direct| would be wrong.
        # A skipped model remains visible below MS, with an explicit warning.
        try:
            matched = reconstruct_matched_continuation(payload)
        except (KeyError, TypeError, ValueError, IndexError) as exc:
            result["continuation_issues"][name] = str(exc)
            continue
        result["continuation_methods"][name] = matched.method
        magnitudes = {
            "total": np.abs(matched.combined),
            "hard": np.broadcast_to(np.abs(matched.hard), matched.direct.shape),
            "direct": np.abs(matched.direct),
        }
        for qkey in ("c5_norm",) + tuple(
            f"c5_{i+1}{j+1}" for i in range(3) for j in range(i, 3)
        ):
            if name not in result["quantities"].get(qkey, {}):
                continue
            selector = QUANTITIES[qkey][4]
            def values(matrix: np.ndarray) -> list[float]:
                if selector == "frobenius":
                    return np.sqrt(np.sum(matrix * matrix, axis=(1, 2))).tolist()
                i, j = selector
                return matrix[:, i, j].tolist()
            x = matched.mu_gev.tolist()
            result["matched"].setdefault(qkey, {})[name] = {
                "x": x, "y": values(magnitudes["total"]),
            }
            result["components"].setdefault(qkey, {})[name] = {
                "hard": {"x": x, "y": values(magnitudes["hard"])},
                "direct": {"x": x, "y": values(magnitudes["direct"])},
            }
    result["models"] = sorted(set(result["models"]))
    return result


def write_dashboard(study_dir: Path, output_path: Path, title: str) -> Path:
    data = collect(study_dir)
    if not data["models"]:
        raise FileNotFoundError(f"No successful running_diagnostics.json files under {study_dir}")

    metadata = {
        key: {"label": spec[0], "unit": spec[3], "block": spec[1],
              "weinberg": key.startswith("c5_")}
        for key, spec in QUANTITIES.items()
    }
    metadata["llss_kernel"] = {
        "label": "LLSS one-loop kernel η", "unit": "dimensionless",
        "block": "analytic_kernel", "weinberg": False,
    }
    # The direct-only arrays remain extractable for compatibility/tests, but
    # the default UI shows one coherent Weinberg view per matrix component.
    c5_keys = ["c5_norm"] + [f"c5_{i+1}{j+1}" for i in range(3) for j in range(i, 3)]
    other_keys = [key for key in QUANTITIES if not key.startswith(("c5_", "direct_c5_"))]
    def group(label: str, keys: list[str]) -> str:
        values = "".join(
            f'<option value="{html.escape(key)}">{html.escape(str(metadata[key]["label"]))}</option>'
            for key in keys
        )
        return f'<optgroup label="{html.escape(label)}">{values}</optgroup>'

    options = (
        group("Weinberg operator — complete matching history", c5_keys)
        + group("Weinberg-like LLSS — analytic diagnostic", ["llss_kernel"])
        + group("Other running observables", other_keys)
    )
    payload_json = json.dumps(data, separators=(",", ":"))
    metadata_json = json.dumps(metadata, separators=(",", ":"))
    title_html = html.escape(title)

    template = r"""<!doctype html>
<html lang="en">
<head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
<style>
:root{color-scheme:light;--ink:#17243b;--muted:#556375;--line:#d7dee8;--paper:#fff;--surface:#f7f9fc;--accent:#225ba8}
*{box-sizing:border-box}
[hidden]{display:none!important}
body{font-family:system-ui,-apple-system,Segoe UI,Arial,sans-serif;margin:0;background:var(--surface);color:var(--ink)}
header{padding:20px 28px;background:#fff;border-bottom:1px solid var(--line)}
header h2{font-size:20px;margin:0 0 5px;font-weight:650}
.small{font-size:12px;color:var(--muted);line-height:1.5}
.layout{display:grid;grid-template-columns:275px minmax(0,1fr);min-height:calc(100vh - 100px)}
aside{padding:20px 17px;background:#fff;border-right:1px solid var(--line);overflow:auto}
main{padding:20px 22px;min-width:0}
label.caption{display:block;font-size:12px;font-weight:650;margin-bottom:8px}
select{width:100%;padding:10px 8px;border:1px solid #bac7d8;border-radius:7px;background:#fff;color:var(--ink);font-size:13px}
.actions{display:flex;gap:8px;margin:14px 0 12px}
button{padding:8px 12px;border:1px solid #ccd5e2;border-radius:7px;background:white;color:var(--ink);font-size:12px;cursor:pointer;white-space:nowrap}
button:hover{background:#eef4ff;border-color:#acc6e9}
button:focus-visible,select:focus-visible,input:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.model-list{max-height:60vh;overflow:auto}
label.model{display:flex;gap:7px;align-items:center;font-size:12px;padding:5px 0;cursor:pointer;overflow-wrap:anywhere}
label.model input{margin:0;accent-color:var(--accent)}
.panel{border:1px solid var(--line);border-radius:10px;background:var(--paper);padding:16px;box-shadow:0 1px 2px #182b480b}
.zoom{display:flex;flex-wrap:wrap;gap:7px;margin-bottom:12px}
.zoom button{margin:0}
.zoom button.active{background:#e9f2ff;border-color:#76a8e2;color:#153f79}
.axis-controls{display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin:0 0 12px;padding:9px 12px;background:#f7f9fc;border:1px solid #e0e7f0;border-radius:7px}
.axis-controls label{font-size:12px;font-weight:650;white-space:nowrap}
.axis-controls select{width:auto;max-width:100%;padding:6px 8px;font-size:12px}
.axis-help{font-size:11px;color:var(--muted)}
.axis-title{font-size:13px;font-weight:600;color:var(--ink);margin:6px 0 5px 5px}
.component-toggle{display:flex;align-items:center;gap:9px;width:fit-content;margin:0 0 12px;padding:9px 11px;border-radius:7px;border:1px solid #d5dfed;background:#f8fafc;color:#33465e;font-size:12px;cursor:pointer}
.component-toggle input{accent-color:var(--accent);margin:0}
canvas{width:100%;height:min(70vh,720px);min-height:420px;display:block}
.legend{display:flex;gap:20px;flex-wrap:wrap;margin:12px 2px 2px;font-size:12px;color:var(--muted)}
.legend span{display:inline-flex;gap:7px;align-items:center}
.dot{width:12px;height:12px;border-radius:2px;display:inline-block;border:1px solid #c9d0d7}
.line{width:24px;display:inline-block;border-top:2.5px solid var(--ink)}
.line.thin{width:10px;border-top-width:1.2px}
.line.dotted{border-top-style:dotted;border-top-width:2px}
#caption{margin-top:12px;padding:11px 13px;background:#f5f8fc;border-left:3px solid #b6cce6;border-radius:4px}
#status{font-size:12px;color:var(--muted);margin-top:9px}
@media(max-width:830px){.layout{grid-template-columns:1fr}aside{border-right:0;border-bottom:1px solid var(--line)}.model-list{max-height:160px}main{padding:12px}.panel{padding:8px}header{padding:15px}.zoom{gap:5px}canvas{min-height:320px}}
</style></head>
<body>
<header><h2>__TITLE__</h2><div class="small">Numerical integration direction: UV → IR. Saved benchmark trajectories; no matching or RGE recalculation.</div></header>
<div class="layout">
<aside>
  <label class="caption" for="quantity">Displayed quantity</label>
  <select id="quantity">__OPTIONS__</select>
  <div id="modelControls">
    <div class="actions"><button id="all">Show all models</button><button id="none">Hide all models</button></div>
    <label class="caption">Models</label><div id="models" class="model-list"></div>
  </div>
  <p id="modelNote" class="small" hidden>The LLSS kernel is independent of the model and of Yukawa/scalar benchmark values.</p>
</aside>
<main>
  <div class="panel">
    <div class="zoom">
      <button id="resetZoom">Default view</button>
      <button id="zoomAll">All scales</button>
      <button id="zoomUV">Full T3</button>
      <button id="zoomInt">Intermediate EFT</button>
      <button id="zoomFinal">SM + Weinberg</button>
    </div>
    <div class="axis-controls">
      <label for="yScaleMode">Y-axis scale</label>
      <select id="yScaleMode">
        <option value="auto">Auto — visible models / current zoom</option>
        <option value="fixed">Fixed / absolute — all models / all scales</option>
      </select>
      <span class="axis-help" id="axisHelp">Auto resizes when models are shown or hidden.</span>
    </div>
    <label class="component-toggle" id="componentControl" hidden>
      <input id="showComponents" type="checkbox">
      Show hard and direct components (dotted)
    </label>
    <div class="axis-title" id="axisTitle"></div>
    <canvas id="plot" aria-label="Interactive scale dependence of T3 running observables"></canvas>
    <div class="legend" id="regionsLegend">
      <span><i class="dot" style="background:#eaf2ff"></i>Full T3</span>
      <span><i class="dot" style="background:#fff1cf"></i>Intermediate scalar EFT</span>
      <span><i class="dot" style="background:#e5f5e9"></i>SM + Weinberg</span>
    </div>
    <div class="legend" id="weinbergLegend" hidden>
      <span><i class="line"></i>Solid: hard + direct continuation, then final C₅</span>
      <span id="componentLegend" hidden><i class="line dotted"></i>Dotted: fixed hard at Mₛ and direct term</span>
    </div>
    <div id="caption" class="small"></div>
    <div id="status"></div>
  </div>
</main>
</div>
<script>
'use strict';
const DATA=__DATA__;
const META=__META__;
const COLORS=['#2466ae','#c46c25','#278864','#bf434b','#7756a5','#93612d','#b64990','#4c8193','#8a8c24','#356d43','#a55a35','#68429c','#bc6072','#4a7eab','#963f3f','#238fa1'];
const enabled={};
DATA.models.forEach(m=>enabled[m]=true);
const select=document.getElementById('quantity');
const canvas=document.getElementById('plot');
const caption=document.getElementById('caption');
const status=document.getElementById('status');
const modelsDiv=document.getElementById('models');
const modelControls=document.getElementById('modelControls');
const modelNote=document.getElementById('modelNote');
const yScaleMode=document.getElementById('yScaleMode');
const showComponents=document.getElementById('showComponents');
const componentControl=document.getElementById('componentControl');
const componentLegend=document.getElementById('componentLegend');
showComponents.onchange=draw;
yScaleMode.onchange=draw;
const regionColours={uv:'#eaf2ff',int:'#fff1cf',sm:'#e5f5e9'};
DATA.models.forEach((m,i)=>{
 const label=document.createElement('label');label.className='model';
 const check=document.createElement('input');check.type='checkbox';check.checked=true;
 check.onchange=()=>{enabled[m]=check.checked;draw()};
 const text=document.createElement('span');text.textContent=m;text.style.color=COLORS[i%COLORS.length];
 label.append(check,text);modelsDiv.appendChild(label);
});
document.getElementById('all').onclick=()=>setAll(true);
document.getElementById('none').onclick=()=>setAll(false);
function setAll(value){DATA.models.forEach(m=>enabled[m]=value);modelsDiv.querySelectorAll('input').forEach(c=>c.checked=value);draw()}
const view={hi:null,lo:null,mode:'default'};
function isWeinberg(q){return !!(META[q]&&META[q].weinberg)}
function isLLSS(q){return q==='llss_kernel'}
function getScales(){return DATA.scales||{}}
function defaultLimits(q){
 const s=getScales();
 if(isLLSS(q)&&s.mu_fermion_threshold>s.mu_scalar_threshold)return [s.mu_fermion_threshold,s.mu_scalar_threshold];
 if(isWeinberg(q)&&s.mu_fermion_threshold>s.mu_low)return [s.mu_fermion_threshold,s.mu_low];
 return [s.mu_uv||1e7,s.mu_low||100];
}
function changeView(mode,hi,lo){
 if(!(hi>lo&&lo>0))return;
 view.hi=hi;view.lo=lo;view.mode=mode;draw();
}
function resetView(){const a=defaultLimits(select.value);changeView('default',a[0],a[1])}
select.onchange=resetView;
document.getElementById('resetZoom').onclick=resetView;
document.getElementById('zoomAll').onclick=()=>{const s=getScales();changeView('all',s.mu_uv||1e7,s.mu_low||100)};
document.getElementById('zoomUV').onclick=()=>{const s=getScales();changeView('uv',s.mu_uv||1e7,s.mu_fermion_threshold||1e5)};
document.getElementById('zoomInt').onclick=()=>{const s=getScales();changeView('int',s.mu_fermion_threshold||1e5,s.mu_scalar_threshold||1e3)};
document.getElementById('zoomFinal').onclick=()=>{const s=getScales();changeView('sm',s.mu_scalar_threshold||1e3,s.mu_low||100)};
function kernelCurve(){
 const s=getScales(),hi=s.mu_fermion_threshold,lo=s.mu_scalar_threshold;
 if(!(hi>lo&&lo>0))return null;
 const x=[],y=[],hbar=1/(16*Math.PI*Math.PI);
 for(let i=0;i<=120;i++){
  const mu=hi*Math.pow(lo/hi,i/120);
  x.push(mu);y.push(Math.log(mu/hi)*hbar);
 }
 return {x,y};
}
function draw(){
 const q=select.value,weinberg=isWeinberg(q),llss=isLLSS(q);
 const traces=DATA.quantities[q]||{},matchedCurves=DATA.matched[q]||{};
 const components=DATA.components[q]||{};
 const showParts=weinberg&&showComponents.checked;
 const kernel=llss?kernelCurve():null;
 modelControls.hidden=llss;modelNote.hidden=!llss;
 document.getElementById('weinbergLegend').hidden=!weinberg;
 componentControl.hidden=!weinberg;
 componentLegend.hidden=!showParts;
 const dpr=window.devicePixelRatio||1,rect=canvas.getBoundingClientRect();
 const W=Math.max(520,rect.width),H=Math.max(360,rect.height);
 canvas.width=Math.round(W*dpr);canvas.height=Math.round(H*dpr);
 const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);
 ctx.clearRect(0,0,W,H);
 const L=89,R=24,T=43,B=65,right=W-R,bottom=H-B;
 const s=getScales();
 if(!(view.hi>view.lo)){const v=defaultLimits(q);view.hi=v[0];view.lo=v[1]}
 const hi=Math.log10(view.hi),lo=Math.log10(view.lo);
 const X=x=>L+(hi-Math.log10(x))/(hi-lo)*(right-L);
 // Fixed mode uses every model and every stored scale for this quantity.
 // Auto mode uses only enabled models inside the current horizontal viewport.
 // Weinberg range includes the verified complex hard+direct continuation.
 // Optional hard/direct component overlays affect auto range only when visible.
 const fixedY=yScaleMode.value==='fixed';
 const allTraces=[];
 if(llss){if(kernel)allTraces.push(kernel)}
 else DATA.models.forEach(m=>{
  if(!fixedY&&!enabled[m])return;
  if(traces[m])allTraces.push(traces[m]);
  if(weinberg&&matchedCurves[m])allTraces.push(matchedCurves[m]);
  if(showParts&&components[m]){
   allTraces.push(components[m].hard,components[m].direct);
  }
 });
 const visible=[];
 allTraces.forEach(tr=>{for(let i=0;i<tr.x.length;i++){
  const x=tr.x[i],y=tr.y[i];
  if(!(x>0&&Number.isFinite(y)))continue;
  if(fixedY||(Math.log10(x)>=lo-1e-8&&Math.log10(x)<=hi+1e-8))visible.push(y);
 }});
 let ymin=visible.length?Math.min(...visible):0,ymax=visible.length?Math.max(...visible):1;
 if(weinberg){ymin=Math.min(0,ymin);ymax=ymax||1e-10}
 if(ymax===ymin){const span=Math.abs(ymax)*.1||1;ymin-=span;ymax+=span}
 const pad=(ymax-ymin)*.075;
 if(!weinberg)ymin-=pad;
 ymax+=pad;
 // Expose current limits for reproducible comparisons and browser tests.
 canvas.dataset.yMin=String(ymin);
 canvas.dataset.yMax=String(ymax);
 const Y=y=>T+(ymax-y)/(ymax-ymin)*(bottom-T);
 document.getElementById('axisHelp').textContent=fixedY?
  'Locked to the full stored range for this quantity, regardless of visible models or zoom.':
  'Rescales to visible models in the current horizontal window.';
 function rectRegion(h,l,color,name){
  if(!(h>l&&l>0))return;
  const a=Math.min(h,view.hi),b=Math.max(l,view.lo);
  if(a<=b)return;
  const x1=X(a),x2=X(b);
  ctx.fillStyle=color;ctx.fillRect(x1,T,x2-x1,bottom-T);
  ctx.font='11px system-ui';ctx.fillStyle='#596579';
  const labelWidth=ctx.measureText(name).width;
  if(x2-x1>labelWidth+18)ctx.fillText(name,x1+9,T+16);
 }
 rectRegion(s.mu_uv,s.mu_fermion_threshold,regionColours.uv,'Full T3');
 rectRegion(s.mu_fermion_threshold,s.mu_scalar_threshold,regionColours.int,'Intermediate scalar EFT');
 rectRegion(s.mu_scalar_threshold,s.mu_low,regionColours.sm,'SM + Weinberg EFT');
 ctx.strokeStyle='#d8e0eb';ctx.lineWidth=1;ctx.font='11px system-ui';
 for(let k=0;k<=5;k++){
  const value=ymin+k*(ymax-ymin)/5,py=Y(value);
  ctx.beginPath();ctx.moveTo(L,py);ctx.lineTo(right,py);ctx.stroke();
  ctx.fillStyle='#516074';ctx.fillText(llss?value.toFixed(4):value.toExponential(2),6,py+4);
 }
 for(let p=Math.floor(hi);p>=Math.ceil(lo);p--){
  const xp=X(Math.pow(10,p));ctx.strokeStyle='#e1e7ef';ctx.beginPath();ctx.moveTo(xp,T);ctx.lineTo(xp,bottom);ctx.stroke();
  ctx.fillStyle='#516074';ctx.font='11px system-ui';ctx.fillText('10^'+p,xp-12,bottom+20);
 }
 ctx.strokeStyle='#607084';ctx.beginPath();ctx.moveTo(L,T);ctx.lineTo(L,bottom);ctx.lineTo(right,bottom);ctx.stroke();
 const thresholds=[['F',s.mu_fermion_threshold],['S',s.mu_scalar_threshold]];
 thresholds.forEach(([name,mu])=>{
  if(!(mu>0&&mu<=view.hi*1.00001&&mu>=view.lo*.99999))return;
  const px=X(mu);ctx.save();ctx.setLineDash([4,4]);ctx.strokeStyle='#8792a3';
  ctx.beginPath();ctx.moveTo(px,T);ctx.lineTo(px,bottom);ctx.stroke();ctx.restore();
  ctx.font='11px system-ui';ctx.fillStyle='#536277';
  ctx.fillText(name+' threshold',Math.min(px+5,right-72),T-12);
 });
 function drawLine(tr,col,width,dash=[],opacity=1){
  if(!tr||tr.x.length<2)return;
  ctx.save();ctx.strokeStyle=col;ctx.lineWidth=width;ctx.lineJoin='round';ctx.lineCap='round';
  ctx.setLineDash(dash);ctx.globalAlpha=opacity;
  ctx.beginPath();let started=false;
  for(let i=0;i<tr.x.length;i++){
   const x=tr.x[i],y=tr.y[i];
   if(!(x>0&&Number.isFinite(y)))continue;
   const px=X(x),py=Y(y);
   if(!started){ctx.moveTo(px,py);started=true}else ctx.lineTo(px,py);
  }
  if(started)ctx.stroke();ctx.restore();
 }
 let joined=0;
 ctx.save();ctx.beginPath();ctx.rect(L,T,right-L,bottom-T);ctx.clip();
 if(llss){if(kernel)drawLine(kernel,'#245c98',2.9)}
 else DATA.models.forEach((m,i)=>{
  if(!enabled[m])return;
  const col=COLORS[i%COLORS.length];
  if(weinberg){
   const continuation=matchedCurves[m],final=traces[m],parts=components[m];
   // The optional hard component is a fixed MS matching anchor, NOT an
   // intermediate-EFT hard Wilson coefficient with its own beta function.
   if(showParts&&parts){
    drawLine(parts.hard,col,1.5,[2,5],.45);
    drawLine(parts.direct,col,1.7,[2,4],.75);
   }
   // The two solid segments meet at MS because Python validated complex
   // C5_hard + C5_direct = C5_final before taking magnitudes.
   drawLine(continuation,col,2.6);
   drawLine(final,col,2.6);
   if(continuation&&final)joined++;
  }else drawLine(traces[m],col,2.3);
 });
 ctx.restore();
 ctx.fillStyle='#26364b';ctx.font='12px system-ui';
 ctx.fillText('Renormalisation scale μ [GeV]    UV → IR',Math.max(L,(W-280)/2),H-13);
 const ylabel=llss?'LLSS running kernel η (dimensionless)':META[q].label+(META[q].unit?' ['+META[q].unit+']':'');
 document.getElementById('axisTitle').textContent=ylabel;
 if(llss){
  caption.textContent='LLSS intermediate EFT: η(μ)=ln(μ/MF)/(16π²). Fixed-order C_LLSS(μ)=C_LLSS(MF)+η(μ)β^(1). This blue line is the universal dimensionless kernel, NOT a numerically evaluated LLSS Wilson coefficient.';
 }else if(weinberg){
  const count=DATA.models.filter(m=>enabled[m]&&traces[m]).length;
  const missing=DATA.models.filter(m=>enabled[m]&&traces[m]&&!matchedCurves[m]);
  const methods=new Set(DATA.models.filter(m=>enabled[m]&&matchedCurves[m]).map(m=>DATA.continuation_methods[m]));
  caption.textContent=count+' model(s) shown. Yellow: threshold-anchored |C₅hard(Mₛ) + ΔC₅direct(μ)| (NOT an intermediate-EFT Wilson coefficient). Green: final evolved |C₅(μ)|. '+joined+' continuity checks passed. '+
   (methods.has('validated_real_log_affine_reconstruction')?'Legacy real-only trajectories reconstructed and checked against every saved intermediate magnitude. ':'')+
   (showParts?'Dotted: individual |hard(Mₛ)| and |direct(μ)| magnitudes; they must not be summed as magnitudes. ':'')+
   (missing.length?missing.length+' model(s) lack a verified complex reconstruction: '+missing.map(m=>m+' ('+(DATA.continuation_issues[m]||'missing')+')').join('; '):'');
 }else{
  caption.textContent=DATA.models.filter(m=>enabled[m]&&traces[m]).length+' model(s) shown. Curves are taken from saved diagnostics.';
 }
 status.textContent=(!DATA.models.some(m=>enabled[m])&&!llss)?
  'No models selected. Select a model to display its curve.':
  (visible.length?'':'No stored data in this selected scale window.');
 document.querySelectorAll('.zoom button').forEach(b=>b.classList.remove('active'));
 const activeId={default:'resetZoom',all:'zoomAll',uv:'zoomUV',int:'zoomInt',sm:'zoomFinal'}[view.mode];
 if(activeId)document.getElementById(activeId).classList.add('active');
}
canvas.addEventListener('wheel',event=>{
 event.preventDefault();
 const bounds=canvas.getBoundingClientRect(),L=89,R=24;
 const f=(event.offsetX-L)/(bounds.width-L-R);
 if(!(f>=0&&f<=1))return;
 const hi=Math.log10(view.hi),lo=Math.log10(view.lo),center=hi-f*(hi-lo);
 const all=getScales(),max=Math.log10(all.mu_uv||1e7),min=Math.log10(all.mu_low||100);
 const span=Math.min(max-min,Math.max(.025,(hi-lo)*(event.deltaY>0?1.25:.8)));
 let nextHi=center+span*f,nextLo=center-span*(1-f);
 if(nextHi>max){nextLo-=nextHi-max;nextHi=max}
 if(nextLo<min){nextHi+=min-nextLo;nextLo=min}
 changeView('custom',Math.pow(10,nextHi),Math.pow(10,nextLo));
},{passive:false});
canvas.addEventListener('dblclick',resetView);
window.addEventListener('resize',draw);
resetView();
</script></body></html>"""

    page = (
        template.replace("__TITLE__", title_html)
        .replace("__OPTIONS__", options)
        .replace("__DATA__", payload_json)
        .replace("__META__", metadata_json)
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(page, encoding="utf-8")
    return output_path


def refresh_full_comparison(
    input_root: Path = PROJECT_ROOT / "output" / "full" / "comparison",
    report_root: Path = PROJECT_ROOT / "Reports" / "output" / "full" / "comparison",
    display_root: Path | None = PROJECT_ROOT / "Reports" / "output" / "display" / "comparison",
) -> list[Path]:
    """Regenerate dashboards from saved diagnostics; optionally update display copies."""
    generated: list[Path] = []
    for scenario in COMPARISON_SCENARIOS:
        source = Path(input_root) / scenario
        output = Path(report_root) / scenario / "interactive_comparison.html"
        path = write_dashboard(source, output, f"T3 normalized comparison: {scenario}")
        generated.append(path)
        print(f"Updated comparison dashboard: {path}")
        if display_root is not None:
            display = Path(display_root) / scenario / path.name
            display.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, display)
            print(f"Updated display copy: {display}")
    return generated


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study_dir", nargs="?", type=Path)
    parser.add_argument("output", nargs="?", type=Path)
    parser.add_argument("--title", default="T3 model comparison")
    parser.add_argument("--refresh-full-comparison", action="store_true",
                        help="Refresh all four benchmark HTML dashboards from saved diagnostics")
    parser.add_argument("--input-root", type=Path,
                        default=PROJECT_ROOT / "output" / "full" / "comparison")
    parser.add_argument("--report-root", type=Path,
                        default=PROJECT_ROOT / "Reports" / "output" / "full" / "comparison")
    parser.add_argument("--display-root", type=Path,
                        default=PROJECT_ROOT / "Reports" / "output" / "display" / "comparison")
    parser.add_argument("--no-display-copy", action="store_true")
    args = parser.parse_args()
    if args.refresh_full_comparison:
        if args.study_dir is not None or args.output is not None:
            parser.error("Do not supply positional paths with --refresh-full-comparison")
        refresh_full_comparison(
            input_root=args.input_root,
            report_root=args.report_root,
            display_root=None if args.no_display_copy else args.display_root,
        )
        return 0
    if args.study_dir is None or args.output is None:
        parser.error("Supply study_dir and output, or --refresh-full-comparison")
    print(write_dashboard(args.study_dir, args.output, args.title))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
