"""Build a standalone interactive HTML comparison of T3 numerical diagnostics.

This is presentation-only: it reads existing running_diagnostics.json files and
never recomputes matching, RGEs, thresholds, or neutrino observables.
"""
from __future__ import annotations

import argparse
import html
import json
import math
from pathlib import Path
from typing import Any


QUANTITIES = {
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
        QUANTITIES[f"c5_{i+1}{j+1}"] = (
            f"|C5^{{{i+1}{j+1}}}|", "final_running", "c5_abs", "GeV⁻¹", (i, j)
        )


def _extract(payload: dict[str, Any], spec: tuple[Any, ...]) -> tuple[list[float], list[float]]:
    _, block_name, key, _, selector = spec
    block = payload[block_name]
    x = [float(v) for v in block["mu_gev"]]
    raw = block[key]
    if isinstance(selector, tuple):
        i, j = selector
        y = [float(row[i][j]) for row in raw]
    elif isinstance(selector, int) and not isinstance(selector, bool):
        y = [float(row[selector]) for row in raw]
    else:
        y = [float(v) for v in raw]
        if selector is True:
            y = [abs(v) for v in y]
    return x, y


def _model_name(path: Path, payload: dict[str, Any]) -> str:
    # Parent of data/running_diagnostics.json is the model output directory.
    model_dir = path.parent.parent
    rep = payload.get("benchmark_search", {}).get("model")
    if isinstance(rep, dict):
        return str(rep.get("model_key", model_dir.name))
    return model_dir.name


def collect(study_dir: Path) -> dict[str, Any]:
    result: dict[str, Any] = {"quantities": {}, "models": []}
    files = sorted(Path(study_dir).rglob("running_diagnostics.json"))
    for path in files:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            if payload.get("status") != "Success":
                continue
            name = _model_name(path, payload)
            result["models"].append(name)
            for qkey, spec in QUANTITIES.items():
                try:
                    x, y = _extract(payload, spec)
                except (KeyError, TypeError, ValueError, IndexError):
                    continue
                result["quantities"].setdefault(qkey, {})[name] = {"x": x, "y": y}
        except (OSError, json.JSONDecodeError):
            continue
    result["models"] = sorted(set(result["models"]))
    return result


def write_dashboard(study_dir: Path, output_path: Path, title: str) -> Path:
    data = collect(study_dir)
    if not data["models"]:
        raise FileNotFoundError(f"No successful running_diagnostics.json files under {study_dir}")

    metadata = {
        key: {"label": spec[0], "unit": spec[3]}
        for key, spec in QUANTITIES.items()
    }
    options = "\n".join(
        f'<option value="{html.escape(key)}">{html.escape(str(spec[0]))}</option>'
        for key, spec in QUANTITIES.items()
    )
    payload_json = json.dumps(data, separators=(",", ":"))
    metadata_json = json.dumps(metadata, separators=(",", ":"))
    title_html = html.escape(title)

    page = f'''<!doctype html>
<html><head><meta charset="utf-8"><title>{title_html}</title>
<style>
body{{font-family:system-ui,Segoe UI,Arial,sans-serif;margin:0;background:#fff;color:#111}}
header{{padding:18px 24px;border-bottom:1px solid #ccc}} .wrap{{display:grid;grid-template-columns:280px 1fr;min-height:85vh}}
aside{{padding:16px;border-right:1px solid #ddd;overflow:auto}} main{{padding:16px}}
select,button{{width:100%;padding:8px;margin:4px 0 10px}} label.model{{display:block;padding:3px 0;font-size:13px}}
canvas{{width:100%;height:70vh;border:1px solid #bbb}} .small{{font-size:12px;color:#555}}
</style></head><body>
<header><h2 style="margin:0">{title_html}</h2><div class="small">Interactive comparison from stored numerical diagnostics; no physics is recomputed.</div></header>
<div class="wrap"><aside>
<label>Quantity</label><select id="quantity">{options}</select>
<button id="all">Show all models</button><button id="none">Hide all models</button>
<div id="models"></div></aside><main><canvas id="plot"></canvas><div id="caption" class="small"></div></main></div>
<script>
const DATA={payload_json}; const META={metadata_json};
const colors=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf','#003f5c','#58508d','#bc5090','#ff6361','#ffa600','#2f4b7c','#665191','#a05195'];
const enabled={{}}; DATA.models.forEach(m=>enabled[m]=true);
const modelsDiv=document.getElementById('models');
DATA.models.forEach((m,i)=>{{let l=document.createElement('label');l.className='model';let c=document.createElement('input');c.type='checkbox';c.checked=true;c.onchange=()=>{{enabled[m]=c.checked;draw()}};l.appendChild(c);let s=document.createElement('span');s.textContent=' '+m;s.style.color=colors[i%colors.length];l.appendChild(s);modelsDiv.appendChild(l)}});
document.getElementById('all').onclick=()=>setAll(true); document.getElementById('none').onclick=()=>setAll(false);
document.getElementById('quantity').onchange=draw;
function setAll(v){{DATA.models.forEach(m=>enabled[m]=v);modelsDiv.querySelectorAll('input').forEach(c=>c.checked=v);draw()}}
function draw(){{
 const q=document.getElementById('quantity').value, traces=DATA.quantities[q]||{{}}, canvas=document.getElementById('plot');
 const dpr=window.devicePixelRatio||1, r=canvas.getBoundingClientRect();canvas.width=Math.max(600,r.width*dpr);canvas.height=Math.max(400,r.height*dpr);const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);const W=r.width,H=r.height;
 ctx.clearRect(0,0,W,H);const L=82,R=20,T=28,B=62;let xs=[],ys=[];
 DATA.models.forEach(m=>{{if(enabled[m]&&traces[m]){{xs.push(...traces[m].x.filter(v=>v>0));ys.push(...traces[m].y.filter(Number.isFinite))}}}});
 if(!xs.length||!ys.length){{ctx.fillText('No visible data',L,T+20);return}}
 const lx=xs.map(Math.log10), xmin=Math.min(...lx),xmax=Math.max(...lx); let ymin=Math.min(...ys),ymax=Math.max(...ys);if(ymin===ymax){{ymin-=1;ymax+=1}};const pad=.05*(ymax-ymin);ymin-=pad;ymax+=pad;
 const X=x=>L+(Math.log10(x)-xmin)/(xmax-xmin||1)*(W-L-R), Y=y=>T+(ymax-y)/(ymax-ymin)*(H-T-B);
 ctx.strokeStyle='#222';ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(L,T);ctx.lineTo(L,H-B);ctx.lineTo(W-R,H-B);ctx.stroke();ctx.font='12px system-ui';ctx.fillStyle='#222';
 for(let k=0;k<=5;k++){{let y=ymin+k*(ymax-ymin)/5,py=Y(y);ctx.strokeStyle='#eee';ctx.beginPath();ctx.moveTo(L,py);ctx.lineTo(W-R,py);ctx.stroke();ctx.fillStyle='#222';ctx.fillText(y.toExponential(3),4,py+4)}}
 for(let p=Math.ceil(xmin);p<=Math.floor(xmax);p++){{let px=X(10**p);ctx.strokeStyle='#eee';ctx.beginPath();ctx.moveTo(px,T);ctx.lineTo(px,H-B);ctx.stroke();ctx.fillStyle='#222';ctx.fillText('10^'+p,px-15,H-B+20)}}
 DATA.models.forEach((m,i)=>{{if(!enabled[m]||!traces[m])return;let tr=traces[m];ctx.strokeStyle=colors[i%colors.length];ctx.lineWidth=1.7;ctx.beginPath();let started=false;for(let k=0;k<tr.x.length;k++){{let x=tr.x[k],y=tr.y[k];if(!(x>0&&Number.isFinite(y)))continue;let px=X(x),py=Y(y);if(!started){{ctx.moveTo(px,py);started=true}}else ctx.lineTo(px,py)}}ctx.stroke()}});
 ctx.fillStyle='#111';ctx.font='14px system-ui';ctx.fillText('Renormalisation scale μ [GeV] (log scale)',Math.max(L,(W-300)/2),H-15);ctx.save();ctx.translate(18,H/2);ctx.rotate(-Math.PI/2);ctx.fillText(META[q].label+(META[q].unit?' ['+META[q].unit+']':''),0,0);ctx.restore();
 document.getElementById('caption').textContent=DATA.models.filter(m=>enabled[m]&&traces[m]).length+' visible model(s)';
}}
window.onresize=draw;draw();
</script></body></html>'''
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(page, encoding="utf-8")
    return output_path


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study_dir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--title", default="T3 model comparison")
    args=parser.parse_args()
    path=write_dashboard(args.study_dir,args.output,args.title)
    print(path)
    return 0

if __name__ == "__main__": raise SystemExit(main())
