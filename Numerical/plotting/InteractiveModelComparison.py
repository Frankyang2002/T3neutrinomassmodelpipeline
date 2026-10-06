"""Build a standalone interactive HTML comparison of T3 numerical diagnostics.

This is presentation-only: it reads existing running_diagnostics.json files and
never recomputes matching, RGEs, thresholds, or neutrino observables.

The dashboard displays the renormalisation scale in the numerical integration
direction (UV -> IR) and marks the EFT threshold structure.
"""
from __future__ import annotations

import argparse
import html
import json
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
    model_dir = path.parent.parent
    rep = payload.get("benchmark_search", {}).get("model")
    if isinstance(rep, dict):
        return str(rep.get("model_key", model_dir.name))
    return model_dir.name


def _is_canonical_optimal_diagnostic(study_dir: Path, path: Path) -> bool:
    """Reject legacy short-path diagnostics when reading an optimal study.

    Current optimal outputs live below representation-specific directories named
    T3_dS1_*_dS2_*_dF_*_alpha_*.  Older runs may leave short aliases such as
    optimal/T3_B_alpha_m1/data/running_diagnostics.json beside them.  Mixing the
    two generations corrupts the common scale range in the aggregate dashboard.
    """
    study_dir = Path(study_dir)
    if study_dir.name != "optimal":
        return True

    try:
        relative = path.relative_to(study_dir)
    except ValueError:
        return False

    if not relative.parts:
        return False
    return relative.parts[0].startswith("T3_dS1_")


def collect(study_dir: Path) -> dict[str, Any]:
    result: dict[str, Any] = {"quantities": {}, "models": [], "scales": {}}
    study_dir = Path(study_dir)
    files = [
        path
        for path in sorted(study_dir.rglob("running_diagnostics.json"))
        if _is_canonical_optimal_diagnostic(study_dir, path)
    ]
    for path in files:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            if payload.get("status") != "Success":
                continue
            name = _model_name(path, payload)
            result["models"].append(name)

            scales = payload.get("scales_gev", {})
            if isinstance(scales, dict) and not result["scales"]:
                result["scales"] = {
                    str(k): float(v)
                    for k, v in scales.items()
                    if isinstance(v, (int, float))
                }

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
        raise FileNotFoundError(
            f"No successful running_diagnostics.json files under {study_dir}"
        )

    metadata = {
        key: {"label": spec[0], "unit": spec[3], "block": spec[1]}
        for key, spec in QUANTITIES.items()
    }
    options = "\n".join(
        f'<option value="{html.escape(key)}">{html.escape(str(spec[0]))}</option>'
        for key, spec in QUANTITIES.items()
    )
    payload_json = json.dumps(data, separators=(",", ":"))
    metadata_json = json.dumps(metadata, separators=(",", ":"))
    title_html = html.escape(title)

    template = r"""<!doctype html>
<html><head><meta charset="utf-8"><title>__TITLE__</title>
<style>
body{font-family:system-ui,Segoe UI,Arial,sans-serif;margin:0;background:#fff;color:#111}
header{padding:18px 24px;border-bottom:1px solid #ccc}
.wrap{display:grid;grid-template-columns:280px 1fr;min-height:85vh}
aside{padding:16px;border-right:1px solid #ddd;overflow:auto}
main{padding:16px}
select,button{width:100%;padding:8px;margin:4px 0 10px}
label.model{display:block;padding:3px 0;font-size:13px}
canvas{width:100%;height:70vh;border:1px solid #bbb}
.small{font-size:12px;color:#555}
.eftlegend{display:flex;gap:18px;flex-wrap:wrap;margin-top:8px;font-size:12px;color:#444}
.swatch{display:inline-block;width:16px;height:10px;margin-right:5px;border:1px solid #aaa;vertical-align:-1px}
</style></head><body>
<header>
<h2 style="margin:0">__TITLE__</h2>
<div class="small">RG flow is displayed in the numerical integration direction: UV → IR. Stored diagnostics only; no physics is recomputed.</div>
</header>
<div class="wrap"><aside>
<label>Quantity</label><select id="quantity">__OPTIONS__</select>
<button id="all">Show all models</button><button id="none">Hide all models</button>
<div id="models"></div></aside>
<main>
<div style="display:flex;gap:8px;margin-bottom:8px">
<button id="resetZoom" style="width:auto">Reset zoom</button>
<button id="zoomUV" style="width:auto">Full T3</button>
<button id="zoomInt" style="width:auto">Intermediate EFT</button>
<button id="zoomFinal" style="width:auto">SM + Weinberg</button>
</div>
<canvas id="plot"></canvas>
<div class="eftlegend">
<span><i class="swatch" style="background:rgba(214,232,255,.65)"></i>Full T3</span>
<span><i class="swatch" style="background:rgba(255,235,190,.65)"></i>Intermediate scalar EFT</span>
<span><i class="swatch" style="background:rgba(216,242,216,.65)"></i>SM + Weinberg EFT</span>
</div>
<div id="caption" class="small"></div>
</main></div>
<script>
const DATA=__DATA__;
const META=__META__;
const colors=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf','#003f5c','#58508d','#bc5090','#ff6361','#ffa600','#2f4b7c'];
const enabled={};
DATA.models.forEach(m=>enabled[m]=true);
const modelsDiv=document.getElementById('models');
DATA.models.forEach((m,i)=>{
 let l=document.createElement('label');l.className='model';
 let c=document.createElement('input');c.type='checkbox';c.checked=true;
 c.onchange=()=>{enabled[m]=c.checked;draw()};
 l.appendChild(c);
 let s=document.createElement('span');s.textContent=' '+m;s.style.color=colors[i%colors.length];
 l.appendChild(s);modelsDiv.appendChild(l)
});
document.getElementById('all').onclick=()=>setAll(true);
document.getElementById('none').onclick=()=>setAll(false);
document.getElementById('quantity').onchange=()=>{window.viewHi=undefined;window.viewLo=undefined;draw()};
function setAll(v){
 DATA.models.forEach(m=>enabled[m]=v);
 modelsDiv.querySelectorAll('input').forEach(c=>c.checked=v);
 draw()
}
function draw(){
 const q=document.getElementById('quantity').value;
 const traces=DATA.quantities[q]||{};
 const canvas=document.getElementById('plot');
 const dpr=window.devicePixelRatio||1,r=canvas.getBoundingClientRect();
 canvas.width=Math.max(600,r.width*dpr);canvas.height=Math.max(400,r.height*dpr);
 const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);
 const W=r.width,H=r.height;
 ctx.clearRect(0,0,W,H);
 const L=82,R=20,T=52,B=68;
 let xs=[],ys=[];
 DATA.models.forEach(m=>{
   if(enabled[m]&&traces[m]){
     xs.push(...traces[m].x.filter(v=>v>0));
     ys.push(...traces[m].y.filter(Number.isFinite))
   }
 });
 if(!xs.length||!ys.length){ctx.fillText('No visible data',L,T+20);return}
 const scaleData=DATA.scales||{};
 const fullHi=Math.log10(scaleData.mu_uv||1e12);
 const fullLo=Math.log10(scaleData.mu_low||1e3);
 if(window.viewHi===undefined||window.viewLo===undefined){
   window.viewHi=fullHi; window.viewLo=fullLo;
 }
 const xmax=window.viewHi, xmin=window.viewLo;

 let visibleY=[];
 DATA.models.forEach(m=>{
   if(!enabled[m]||!traces[m])return;
   for(let k=0;k<traces[m].x.length;k++){
     const lx=Math.log10(traces[m].x[k]), y=traces[m].y[k];
     if(lx>=xmin&&lx<=xmax&&Number.isFinite(y))visibleY.push(y);
   }
 });
 if(visibleY.length)ys=visibleY;
 let ymin=Math.min(...ys),ymax=Math.max(...ys);
 if(ymin===ymax){const d=Math.abs(ymin)*.05||1;ymin-=d;ymax+=d}
 const pad=.05*(ymax-ymin);ymin-=pad;ymax+=pad;

 const X=x=>L+(xmax-Math.log10(x))/(xmax-xmin||1)*(W-L-R);
 const Y=y=>T+(ymax-y)/(ymax-ymin)*(H-T-B);

 function shadeRegion(high,low,color,label){
   if(!(high>0&&low>0))return;
   const hi=Math.min(high,10**xmax),lo=Math.max(low,10**xmin);
   if(!(hi>lo))return;
   const x1=X(hi),x2=X(lo);
   ctx.fillStyle=color;ctx.fillRect(x1,T,x2-x1,H-T-B);
   ctx.fillStyle='#555';ctx.font='11px system-ui';
   ctx.fillText(label,x1+6,T+14);
 }
 const s=DATA.scales||{};
 const uv=s.mu_uv,f=s.mu_fermion_threshold,sc=s.mu_scalar_threshold,low=s.mu_low;
 shadeRegion(uv,f,'rgba(214,232,255,.45)','Full T3');
 shadeRegion(f,sc,'rgba(255,235,190,.45)','Intermediate scalar EFT');
 shadeRegion(sc,low,'rgba(216,242,216,.45)','SM + Weinberg EFT');

 ctx.strokeStyle='#222';ctx.lineWidth=1;
 ctx.beginPath();ctx.moveTo(L,T);ctx.lineTo(L,H-B);ctx.lineTo(W-R,H-B);ctx.stroke();
 ctx.font='12px system-ui';
 for(let k=0;k<=5;k++){
   let y=ymin+k*(ymax-ymin)/5,py=Y(y);
   ctx.strokeStyle='#e5e5e5';ctx.beginPath();ctx.moveTo(L,py);ctx.lineTo(W-R,py);ctx.stroke();
   ctx.fillStyle='#222';ctx.fillText(y.toExponential(3),4,py+4)
 }
 for(let p=Math.floor(xmax);p>=Math.ceil(xmin);p--){
   let px=X(10**p);
   ctx.strokeStyle='#e5e5e5';ctx.beginPath();ctx.moveTo(px,T);ctx.lineTo(px,H-B);ctx.stroke();
   ctx.fillStyle='#222';ctx.fillText('10^'+p,px-15,H-B+20)
 }

 const thresholds=[
   ['UV start',uv],
   ['F threshold',f],
   ['S threshold',sc],
   ['low scale',low]
 ];
 thresholds.forEach(([label,value])=>{
   if(!(value>0))return;
   const lv=Math.log10(value);
   if(lv<xmin-1e-10||lv>xmax+1e-10)return;
   const px=X(value);
   ctx.save();ctx.setLineDash([5,4]);ctx.strokeStyle='#666';
   ctx.beginPath();ctx.moveTo(px,T);ctx.lineTo(px,H-B);ctx.stroke();ctx.restore();
   ctx.fillStyle='#444';ctx.font='11px system-ui';
   ctx.save();ctx.translate(px+5,T+32);ctx.rotate(-Math.PI/2);
   ctx.fillText(label+'  '+Number(value).toExponential(2)+' GeV',0,0);ctx.restore()
 });

 DATA.models.forEach((m,i)=>{
   if(!enabled[m]||!traces[m])return;
   let tr=traces[m];ctx.strokeStyle=colors[i%colors.length];ctx.lineWidth=1.7;
   ctx.beginPath();let started=false;
   for(let k=0;k<tr.x.length;k++){
     let x=tr.x[k],y=tr.y[k];
     if(!(x>0&&Number.isFinite(y)))continue;
     let px=X(x),py=Y(y);
     if(!started){ctx.moveTo(px,py);started=true}else ctx.lineTo(px,py)
   }
   ctx.stroke()
 });
 ctx.fillStyle='#111';ctx.font='14px system-ui';
 ctx.fillText('Renormalisation scale μ [GeV] (log scale): UV → IR',Math.max(L,(W-390)/2),H-15);
 ctx.save();ctx.translate(18,H/2);ctx.rotate(-Math.PI/2);
 ctx.fillText(META[q].label+(META[q].unit?' ['+META[q].unit+']':''),0,0);ctx.restore();
 document.getElementById('caption').textContent=
   DATA.models.filter(m=>enabled[m]&&traces[m]).length+
   ' visible model(s). Background shading denotes the EFT domain covered by the selected diagnostic.';
}


function setView(hi,lo){window.viewHi=Math.log10(hi);window.viewLo=Math.log10(lo);draw()}
document.getElementById('resetZoom').onclick=()=>{window.viewHi=undefined;window.viewLo=undefined;draw()};
document.getElementById('zoomUV').onclick=()=>{const z=DATA.scales||{};setView(z.mu_uv||1e12,z.mu_fermion_threshold||1e10)};
document.getElementById('zoomInt').onclick=()=>{const z=DATA.scales||{};setView(z.mu_fermion_threshold||1e10,z.mu_scalar_threshold||7e9)};
document.getElementById('zoomFinal').onclick=()=>{const z=DATA.scales||{};setView(z.mu_scalar_threshold||7e9,z.mu_low||1e3)};
canvas.addEventListener('wheel',e=>{
 e.preventDefault();
 const z=DATA.scales||{},hi0=Math.log10(z.mu_uv||1e12),lo0=Math.log10(z.mu_low||1e3);
 if(window.viewHi===undefined){window.viewHi=hi0;window.viewLo=lo0}
 const r=canvas.getBoundingClientRect(),L=82,R=20,w=r.width-L-R;
 if(e.offsetX<L||e.offsetX>r.width-R)return;
 const f=(e.offsetX-L)/w,center=window.viewHi-f*(window.viewHi-window.viewLo);
 const factor=e.deltaY>0?1.35:0.75;
 let span=Math.max(.025,Math.min(hi0-lo0,(window.viewHi-window.viewLo)*factor));
 let hi=center+f*span,lo=center-(1-f)*span;
 if(hi>hi0){lo-=hi-hi0;hi=hi0}
 if(lo<lo0){hi+=lo0-lo;lo=lo0}
 window.viewHi=Math.min(hi0,hi);window.viewLo=Math.max(lo0,lo);draw();
},{passive:false});
canvas.addEventListener('dblclick',()=>{window.viewHi=undefined;window.viewLo=undefined;draw()});

window.onresize=draw;draw();
</script></body></html>"""

    page = (
        template.replace("__TITLE__", title_html)
        .replace("__OPTIONS__", options)
        .replace("__DATA__", payload_json)
        .replace("__META__", metadata_json)
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(page, encoding="utf-8")
    return output_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study_dir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--title", default="T3 model comparison")
    args = parser.parse_args()
    path = write_dashboard(args.study_dir, args.output, args.title)
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
