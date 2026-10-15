"""Generate Study A/B interactive endpoint and prediction-space comparisons.

This is a presentation-only utility. It reads saved running/fixed CSVs and
never recalculates RGEs, masses, matching, fits, or experimental likelihoods.

From the repository root:
    python -m Numerical.plotting.BuildStudyABInteractive
    python -m Numerical.plotting.BuildStudyABInteractive --study A
    python -m Numerical.plotting.BuildStudyABInteractive --study B

The oval display uses empirical descriptive covariance of sampled predictions,
NOT parameter uncertainties or model confidence regions.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OBS = ('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')
LABELS = {'delta_m21_sq_ev2':'Δm²₂₁ [eV²]', 'delta_m3l_sq_ev2':'Δm²₃ℓ [eV²]',
          'sin2_theta12':'sin²θ₁₂','sin2_theta13':'sin²θ₁₃','sin2_theta23':'sin²θ₂₃'}
GROUPS = {(1,3,2):'A',(2,2,1):'B',(2,2,3):'C',(3,1,2):'D',(3,3,2):'E'}

def short_name(model):
    m=re.fullmatch(r'T3_dS1_(\d+)_dS2_(\d+)_dF_(\d+)_alpha_([mp])(\d+)',model)
    if not m:return model
    a,b,c,sign,n=m.groups()
    return f'T3-{GROUPS.get((int(a),int(b),int(c)), "?")}-{sign}{n}'

def clean_float(value):
    try:
        x=float(value)
        return x if math.isfinite(x) else None
    except (ValueError,TypeError):return None

def read_rows(path,study):
    if not path.is_file():raise FileNotFoundError(f'Missing Study {study} CSV: {path}')
    with path.open(encoding='utf-8-sig',newline='') as stream:
        rows=list(csv.DictReader(stream))
    required=['model','status','reference_yukawa_target']+[f'{p}_{k}' for p in ('running','fixed') for k in OBS]
    if study=='B':required.append('mass_factor')
    missing=[x for x in required if not rows or x not in rows[0]]
    if missing:raise ValueError(f'{path}: missing required columns: {missing}')
    results=[]
    for row in rows:
        if row['status']!='Success':continue
        datum={'model':row['model'],'label':short_name(row['model']),
               'yref':clean_float(row['reference_yukawa_target']),
               'factor':clean_float(row.get('mass_factor','1')),
               'running':{k:clean_float(row[f'running_{k}']) for k in OBS},
               'fixed':{k:clean_float(row[f'fixed_{k}']) for k in OBS}}
        if datum['yref'] is None or datum['factor'] is None:continue
        if any(v is None for p in ('running','fixed') for v in datum[p].values()):continue
        results.append(datum)
    if not results:raise ValueError(f'No complete successful rows found in {path}')
    return results,len(rows)

STYLE='''
:root{color-scheme:light;--ink:#17243b;--muted:#556375;--line:#d7dee8;--surface:#f7f9fc;--accent:#225ba8}
*{box-sizing:border-box}body{margin:0;font:13px system-ui,-apple-system,Segoe UI,Arial,sans-serif;background:var(--surface);color:var(--ink)}
header{padding:18px 25px;background:white;border-bottom:1px solid var(--line)}h2{font-size:19px;margin:0 0 5px}.small{font-size:12px;color:var(--muted);line-height:1.6}
.layout{display:grid;grid-template-columns:285px minmax(0,1fr);min-height:calc(100vh - 95px)}aside{padding:18px;background:white;border-right:1px solid var(--line)}main{padding:20px;min-width:0}
.panel{background:white;border:1px solid var(--line);border-radius:10px;padding:15px}select{width:100%;padding:8px;margin:5px 0 12px;background:white;border:1px solid var(--line);border-radius:6px}
button{padding:7px 10px;background:white;border:1px solid var(--line);border-radius:6px;cursor:pointer}button:hover{background:#eff4fc}.controls{display:flex;flex-wrap:wrap;gap:15px;margin-bottom:10px}
.scroll{max-height:49vh;overflow-y:auto}label.choice{display:block;margin:7px 0;font-size:12px;overflow-wrap:anywhere}canvas{width:100%;height:min(70vh,720px);min-height:420px;display:block}
#note{padding:10px;background:#f6f9fc;border-left:3px solid #b6cce6;line-height:1.6}.legend{display:flex;gap:16px;flex-wrap:wrap;color:var(--muted);font-size:12px;margin:9px 0}
@media(max-width:830px){.layout{grid-template-columns:1fr}aside{border-right:0;border-bottom:1px solid var(--line)}main{padding:10px}}
'''

# Experimental centres and 1 sigma marginal widths from the scenario fit
# target are intentionally loaded from the repository NuFIT configuration by
# the original sensitivity generator. This standalone display does not infer
# covariance from marginal errors. The optional overlays below are set by an
# explicitly provided reference JSON, or omitted when unavailable.
def reference_from_json(path):
    if not path or not Path(path).is_file():return None
    data=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    # Do not guess the NuFIT schema: only accept an explicit display envelope.
    ref=data.get('study_display_gaussian_reference')
    if not isinstance(ref,dict):return None
    if not all(k in ref and isinstance(ref[k],dict) and
               clean_float(ref[k].get('centre')) is not None and
               clean_float(ref[k].get('sigma')) not in (None,0) for k in OBS):return None
    return ref

def page(study,kind,rows,ref,source,total):
    title=(f'Study {study} — RGE running vs fixed predictions' if kind=='running'
           else f'Study {study} — prediction-space scan envelopes')
    intro=('Saved final low-energy predictions versus the scanned input parameter; no intermediate-scale trajectories are reconstructed.'
           if kind=='running' else
           '2× empirical covariance ellipses of the discrete parameter scan, not statistical model confidence intervals. Points are the actual calculated predictions.')
    models=sorted({r['model'] for r in rows})
    data={'study':study,'kind':kind,'rows':rows,'models':models,
          'labels':LABELS,'obs':OBS,'reference':ref,'source':str(source), 'total':total}
    json_data=json.dumps(data,ensure_ascii=False,allow_nan=False,separators=(',',':')).replace('</','<\\/')
    if kind=='running':
        selectors='''<label for="observable">Observable</label><select id="observable"></select>
<label for="metric">Vertical quantity</label><select id="metric"><option value="value">Prediction</option><option value="delta">Running − fixed</option><option value="percent">100 × (running − fixed) / fixed [%]</option></select>'''
    else:selectors='''<label for="x">X observable</label><select id="x"></select><label for="y">Y observable</label><select id="y"></select>'''
    study_control=('<label for="yref">Reference Yukawa setting</label><select id="yref"></select>' if study=='B' else '')
    controls=('''<label><input type="checkbox" id="showRunning" checked> Running</label><label><input type="checkbox" id="showFixed" checked> Fixed</label>'''
              if kind=='running' else
              '''<label><input type="checkbox" id="showRunning" checked> Running points/envelopes</label><label><input type="checkbox" id="showFixed" checked> Fixed points/envelopes</label>
<label><input type="checkbox" id="envelopes" checked> Descriptive ovals</label><label><input type="checkbox" id="showPoints" checked> Evaluated points</label><label><input type="checkbox" id="experiment" checked> Experimental Gaussian (if provided)</label>''')
    js=JS
    page=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{title}</title><style>{STYLE}</style></head>
<body><header><h2>{title}</h2><div class="small">{intro}</div></header><div class="layout"><aside>
{selectors}{study_control}<div><button id="all">Show all models</button> <button id="none">Hide all</button></div>
<h4>Models</h4><div id="models" class="scroll"></div><p class="small">{len(models)} models; {len(rows)}/{total} complete successful rows.</p>
<p class="small">{('x-axis: reference Yukawa (log scale)' if study=='A' else 'x-axis: heavy-fermion mass factor k')}</p></aside>
<main><div class="panel"><div class="controls"><button id="reset">Reset zoom</button>{controls}</div>
<canvas id="chart"></canvas><div class="legend"><span>Solid / circles: RGE running</span><span>Dashed / crosses: fixed parameters</span><span>Colours: T3 models</span></div>
<p id="note" class="small"></p><p class="small">Source: {source}. No RGE or matching recomputation.</p></div></main></div><script>'use strict';const D={json_data};{js}</script></body></html>'''
    return page

JS=r'''
const $=id=>document.getElementById(id),c=$('chart'),ctx=c.getContext('2d');
const colors=['#2563a6','#df6c27','#19856e','#8e50b6','#c23f60','#9c8523','#357eaa','#b45789','#637b33','#aa5736','#3a65b6','#7362b3','#aa891e','#29696b','#cc5333','#63717f'];
const col=Object.fromEntries(D.models.map((m,i)=>[m,colors[i%colors.length]]));
const selected=new Set(D.models);let bounds=null;
function option(select,value,label){let o=document.createElement('option');o.value=value;o.textContent=label;select.appendChild(o)}
function unique(arr){return [...new Set(arr)].sort((a,b)=>a-b)}
if(D.kind==='running'){D.obs.forEach(k=>option($('observable'),k,D.labels[k]));}
else{D.obs.forEach(k=>{option($('x'),k,D.labels[k]);option($('y'),k,D.labels[k]);});$('y').value=D.obs[1]}
if(D.study==='B'){unique(D.rows.map(r=>r.yref)).forEach(v=>option($('yref'),v,String(v)));$('yref').value='0.01';}
D.models.forEach(m=>{let label=document.createElement('label');label.className='choice';let box=document.createElement('input');box.type='checkbox';box.checked=true;box.addEventListener('change',()=>{box.checked?selected.add(m):selected.delete(m);bounds=null;draw()});label.append(box,' '+D.rows.find(r=>r.model===m).label);label.style.color=col[m];$('models').append(label)});
$('all').onclick=()=>{D.models.forEach(m=>selected.add(m));$('models').querySelectorAll('input').forEach(e=>e.checked=true);bounds=null;draw()};
$('none').onclick=()=>{selected.clear();$('models').querySelectorAll('input').forEach(e=>e.checked=false);bounds=null;draw()};
$('reset').onclick=()=>{bounds=null;draw()};
document.querySelectorAll('select,.controls input').forEach(x=>x.addEventListener('change',()=>{bounds=null;draw()}));
function baseRows(){return D.rows.filter(r=>selected.has(r.model)&&(D.study==='A'||r.yref===Number($('yref').value)))}
function xs(r){return D.study==='A'?Math.log10(r.yref):r.factor}
function val(r,key,kind){let a=r.running[key],b=r.fixed[key];if(kind==='running')return a;if(kind==='fixed')return b;
if(kind==='delta')return a-b;if(kind==='percent')return b===0?NaN:100*(a-b)/b;return NaN}
function format(v){return Math.abs(v)>=1e4||(Math.abs(v)>0&&Math.abs(v)<0.0001)?v.toExponential(2):v.toPrecision(4)}
function configure(){const dpr=devicePixelRatio||1;let rect=c.getBoundingClientRect();c.width=Math.max(500,Math.round(rect.width*dpr));c.height=Math.max(400,Math.round(rect.height*dpr));ctx.setTransform(dpr,0,0,dpr,0,0);return [c.width/dpr,c.height/dpr]}
function draw(){let [w,h]=configure();ctx.clearRect(0,0,w,h);const L=100,R=28,T=30,B=74,pw=w-L-R,ph=h-T-B;
let arr=baseRows();let k=D.kind==='running'?$('observable').value:null;
let xobs=D.kind==='oval'?$('x').value:null,yobs=D.kind==='oval'?$('y').value:null;
const showR=$('showRunning').checked,showF=$('showFixed').checked;
let pts=[];
if(D.kind==='running'){
 let metric=$('metric').value;
 for(const row of arr)for(const typ of (metric==='value'?['running','fixed']:['delta'])){
   if(typ==='running'&&!showR||typ==='fixed'&&!showF||typ==='delta'&&!showR)continue;
   let y=val(row,k,typ==='delta'?metric:typ);if(Number.isFinite(y))pts.push({x:xs(row),y,model:row.model,typ,row});
 }
}else{
 for(const row of arr)for(const typ of ['running','fixed']){
   if(typ==='running'&&!showR||typ==='fixed'&&!showF)continue;
   let x=row[typ][xobs],y=row[typ][yobs];if(Number.isFinite(x)&&Number.isFinite(y))pts.push({x,y,model:row.model,typ,row});
 }
}
const expRef=D.reference;
let ext=pts.map(p=>[p.x,p.y]);
if(D.kind==='oval'&&$('experiment').checked&&expRef){let rx=expRef[xobs],ry=expRef[yobs];ext.push([rx.centre-3*rx.sigma,ry.centre-3*ry.sigma],[rx.centre+3*rx.sigma,ry.centre+3*ry.sigma]);}
if(!ext.length){$('note').textContent='No selected finite predictions.';return;}
let xmin=Math.min(...ext.map(p=>p[0])),xmax=Math.max(...ext.map(p=>p[0]));
let ymin=Math.min(...ext.map(p=>p[1])),ymax=Math.max(...ext.map(p=>p[1]));
let dx=(xmax-xmin)||Math.max(Math.abs(xmin)*0.1,1e-5),dy=(ymax-ymin)||Math.max(Math.abs(ymin)*0.1,1e-5);
if(!bounds){bounds=[xmin-dx*.11,xmax+dx*.11,ymin-dy*.15,ymax+dy*.15]}
[xmin,xmax,ymin,ymax]=bounds;
const X=v=>L+(v-xmin)/(xmax-xmin)*pw,Y=v=>T+ph-(v-ymin)/(ymax-ymin)*ph;
ctx.fillStyle='#fff';ctx.fillRect(L,T,pw,ph);
ctx.strokeStyle='#e7ebf1';ctx.lineWidth=1;ctx.font='11px system-ui';
for(let i=0;i<=5;i++){let x=xmin+(xmax-xmin)*i/5,y=ymin+(ymax-ymin)*i/5;
 ctx.beginPath();ctx.moveTo(X(x),T);ctx.lineTo(X(x),T+ph);ctx.moveTo(L,Y(y));ctx.lineTo(L+pw,Y(y));ctx.stroke();
 ctx.fillStyle='#526074';ctx.textAlign='center';ctx.fillText(D.kind==='running'&&D.study==='A'?'10^'+x.toFixed(1):format(x),X(x),T+ph+18);
 ctx.textAlign='right';ctx.fillText(format(y),L-8,Y(y)+3)}
ctx.save();ctx.beginPath();ctx.rect(L,T,pw,ph);ctx.clip();
if(D.kind==='oval'&&$('experiment').checked&&expRef){let rx=expRef[xobs],ry=expRef[yobs];for(const n of [3,1]){ctx.beginPath();ctx.ellipse(X(rx.centre),Y(ry.centre),Math.abs(X(rx.centre+n*rx.sigma)-X(rx.centre)),Math.abs(Y(ry.centre+n*ry.sigma)-Y(ry.centre)),0,0,2*Math.PI);ctx.strokeStyle=n===1?'#2b5eb0':'#9bbce6';ctx.lineWidth=n===1?2:1.3;ctx.stroke()}}
for(const m of D.models){let group=pts.filter(p=>p.model===m);
 for(const typ of (D.kind==='running'&&$('metric').value!=='value'?['delta']:['running','fixed'])){
  let seq=group.filter(p=>p.typ===typ).sort((a,b)=>a.row.yref-b.row.yref||a.row.factor-b.row.factor);
  if(!seq.length)continue;ctx.strokeStyle=col[m];ctx.fillStyle=col[m];ctx.lineWidth=1.7;ctx.setLineDash(typ==='fixed'?[5,3]:[]);
  if(D.kind==='running'){ctx.beginPath();seq.forEach((p,i)=>i?ctx.lineTo(X(p.x),Y(p.y)):ctx.moveTo(X(p.x),Y(p.y)));ctx.stroke()}
  ctx.setLineDash([]);
  if(D.kind==='oval'&&$('envelopes').checked&&seq.length>=3){
   const cx=seq.reduce((s,p)=>s+p.x,0)/seq.length,cy=seq.reduce((s,p)=>s+p.y,0)/seq.length;
   let xx=0,xy=0,yy=0;for(const p of seq){let a=p.x-cx,b=p.y-cy;xx+=a*a;xy+=a*b;yy+=b*b}
   xx/=seq.length-1;xy/=seq.length-1;yy/=seq.length-1;
   const scaleX=pw/(xmax-xmin),scaleY=ph/(ymax-ymin);xx*=scaleX*scaleX;xy*=-scaleX*scaleY;yy*=scaleY*scaleY;
   let tr=xx+yy,disc=Math.sqrt(Math.max(0,(xx-yy)**2+4*xy*xy));let la=(tr+disc)/2,lb=(tr-disc)/2;
   if(lb>1e-14&&Number.isFinite(la)){let angle=.5*Math.atan2(2*xy,xx-yy);ctx.beginPath();ctx.ellipse(X(cx),Y(cy),2*Math.sqrt(la),2*Math.sqrt(lb),angle,0,2*Math.PI);ctx.globalAlpha=typ==='fixed'?.45:.85;ctx.strokeStyle=col[m];ctx.lineWidth=typ==='running'?2:1.4;ctx.setLineDash(typ==='fixed'?[5,4]:[]);ctx.stroke();ctx.setLineDash([]);ctx.globalAlpha=1}
  }
  if(D.kind==='running'||$('showPoints').checked){for(const p of seq){let x=X(p.x),y=Y(p.y);if(typ==='fixed'){ctx.beginPath();ctx.moveTo(x-3,y-3);ctx.lineTo(x+3,y+3);ctx.moveTo(x-3,y+3);ctx.lineTo(x+3,y-3);ctx.strokeStyle=col[m];ctx.lineWidth=1.4;ctx.stroke()}else{ctx.beginPath();ctx.arc(x,y,2.8,0,2*Math.PI);ctx.fill()}}}
 }
}
ctx.restore();ctx.strokeStyle='#6b7b91';ctx.strokeRect(L,T,pw,ph);
ctx.fillStyle='#17243b';ctx.font='13px system-ui';ctx.textAlign='center';ctx.fillText(D.kind==='running'?(D.study==='A'?'log₁₀(reference Yukawa)':'Heavy-fermion mass factor k'):D.labels[xobs],L+pw/2,h-18);
ctx.save();ctx.translate(19,T+ph/2);ctx.rotate(-Math.PI/2);let ylabel=D.kind==='oval'?D.labels[yobs]:$('metric').selectedOptions[0].text+' — '+D.labels[k];ctx.fillText(ylabel,0,0);ctx.restore();
$('note').textContent=D.kind==='running'?'Each point is a final low-energy prediction from the saved full RGE pipeline or the frozen-UV control at identical inputs. No μ-dependent neutrino-observable trajectory is inferred from endpoints. When displaying a correction, only the Running checkbox controls its visibility.':
'Coloured ovals: 2× empirical covariance of the discrete parameter scan (not model confidence regions and not experimental uncertainties). Only drawn when the sampled points have nonzero two-dimensional rank. '+
(expRef?'Blue ellipses: axis-aligned illustrative 1σ/3σ Gaussian reference from explicitly supplied centres and marginal errors, not official NuFIT joint contours.':'No experimental Gaussian reference is drawn because the optional explicit reference was not supplied.')+
' Scan parameter is not a statistical probability distribution.';
}
let start=null;
c.addEventListener('wheel',e=>{e.preventDefault();if(!bounds)return;let z=e.deltaY>0?1.2:.83;let [a,b,d,f]=bounds,cx=(a+b)/2,cy=(d+f)/2;bounds=[cx+(a-cx)*z,cx+(b-cx)*z,cy+(d-cy)*z,cy+(f-cy)*z];draw()},{passive:false});
c.addEventListener('dblclick',()=>{bounds=null;draw()});window.addEventListener('resize',draw);draw();
'''

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',choices=['A','B','both'],default='both')
    p.add_argument('--scenario',default='smallY_smallL')
    p.add_argument('--study-a-csv',type=Path)
    p.add_argument('--study-b-csv',type=Path)
    p.add_argument('--output-root',type=Path,default=ROOT/'Reports'/'output'/'display'/'study_comparisons')
    p.add_argument('--reference-json',type=Path,help='Optional explicit study_display_gaussian_reference JSON; otherwise experimental ellipses omitted')
    args=p.parse_args(argv)
    ref=reference_from_json(args.reference_json)
    specs={'A':args.study_a_csv or ROOT/'Reports'/'output'/'full'/'common_scalar_rg_controls'/args.scenario/'running_vs_fixed.csv',
           'B':args.study_b_csv or ROOT/'Reports'/'output'/'full'/'heavy_fermion_scale_comparison'/args.scenario/'running_vs_fixed_mass_scan.csv'}
    for study in ('A','B'):
        if args.study not in (study,'both'):continue
        source=specs[study];rows,total=read_rows(source,study)
        out=args.output_root/args.scenario/f'Study{study}'
        out.mkdir(parents=True,exist_ok=True)
        for kind,name in [('running','interactive_comparison.html'),('oval','prediction_envelopes.html')]:
            target=out/name
            target.write_text(page(study,kind,rows,ref,source,total),encoding='utf-8')
            print(f'Wrote {target} ({len(rows)}/{total} complete points, {len(set(r["model"] for r in rows))} models)')
    return 0

if __name__=='__main__':raise SystemExit(main())
