"""Offline general-scalar T3 reports from saved continuous_v1 artifacts.

Presentation only; no RG, matching or sensitivity scan is performed. Models are
kept in separate matching-scheme records to avoid silently combining unlike EFT
approximations. Genuine T1/T2/T3 scan envelopes require separately saved points.

From repository root:
    python -m scripts.BuildGeneralScalarInteractiveDisplay --run-id continuous_v1
    python -m scripts.BuildGeneralScalarInteractiveDisplay --run-id continuous_v1 \
        --sensitivity-csv path/to/verified_general_scalar_gamma_scan.csv
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import html
import json
import math
import os
from pathlib import Path
import re
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
REL_BASE = Path('output/studies/general_scalar_extended/continuous_full_scalar')
MODELS = (
    'A_+0', 'A_-2', 'A_-4', 'B_+1', 'B_-1', 'B_-3',
    'C_+1', 'C_-1', 'C_-3', 'D_+0', 'D_+2', 'D_-2',
    'E_+0', 'E_+2', 'E_-2', 'E_-4',
)
SOURCES = (
    ('baseline', 'neutrino_v1/summary.csv', 'Saved benchmark C5',
     'Initial general-scalar neutrino benchmark; self-bilinear effects omitted in six cases.'),
    ('self_massdep', 'self_anomalous_massdep_v1/combined_16_summary.csv',
     'Mass-dependent self-bilinear one-scale',
     'Approximate one-scale mass-dependent matching; not sequential scalar matching.'),
    ('direct', 'direct_wilson_sequential_v1/summary.csv',
     'Direct hierarchical Wilson transport',
     'Hard boundary and direct EFT1 Weinberg source only.'),
    ('extra', 'bc_extra_wilson_sequential_v1/summary.csv',
     'B/C selected extra-quartic source',
     'Selected extra-quartic source on the running scalar background; LLSS held fixed.'),
    ('coupled', 'bc_coupled_wilson_v1/summary.csv',
     'B/C coupled L/eC Wilson candidate',
     'Coupled implemented L/eC tensor; remaining full-flavour and finite-threshold terms unimplemented.'),
)
OBS = ('mass_hi_ev', 'mass_mid_ev', 'mass_lo_ev', 'delta_m21_sq_ev2',
       'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
OBS_LABELS = {
    'mass_hi_ev':'largest neutrino mass [eV]', 'mass_mid_ev':'middle neutrino mass [eV]',
    'mass_lo_ev':'smallest neutrino mass [eV]',
    'delta_m21_sq_ev2':'Δm²₂₁ [eV²]', 'delta_m3l_sq_ev2':'Δm²₃ℓ [eV²]',
    'sin2_theta12':'sin²θ₁₂', 'sin2_theta13':'sin²θ₁₃',
    'sin2_theta23':'sin²θ₂₃',
    'gap_small_mid_sq_ev2':'(m_middle² − m_smallest²) [eV²], derived',
    'gap_small_large_sq_ev2':'(m_largest² − m_smallest²) [eV²], derived',
}
SOURCE_ORDER = {name: i for i, (name, *_rest) in enumerate(SOURCES)}
EXPECTED_FORMAT = 't3-continuous-general-scalar-candidate-v1'


def require(condition: bool, msg: str) -> None:
    if not condition:
        raise ValueError('General scalar display: ' + msg)


def number(value: Any) -> float | None:
    if value is None or value == '':
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def read_json(path: Path) -> dict:
    require(path.is_file(), f'missing input {path}')
    payload = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(payload, dict), f'invalid JSON object {path}')
    return payload


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle))


def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def slug(model: str) -> str:
    require(model in MODELS, 'invalid model ' + str(model))
    return model.replace('+', 'p').replace('-', 'm')


def complex_abs(value: Any) -> float | None:
    if not isinstance(value, dict):
        return None
    re_value, im_value = number(value.get('real')), number(value.get('imag'))
    if re_value is None or im_value is None:
        return None
    return math.hypot(re_value, im_value)


def read_continuous(root: Path, *, partial: bool) -> tuple[dict, dict[str, str]]:
    paths = {}
    data = {}
    for model in MODELS:
        path = root / slug(model) / 'continuous_trajectory.json'
        if not path.exists():
            if partial:
                continue
            require(False, f'missing 16-model scalar trajectory for {model}: {path}')
        j = read_json(path)
        require(j.get('model') == model, f'{path}: model mismatch')
        require(j.get('format') == EXPECTED_FORMAT, f'{path}: unexpected format')
        require(j.get('status') == 'CandidateContinuousUVFermionThresholdEFT1Integrated',
                f'{path}: scalar integration not successful')
        uv, eft = j.get('uv_trajectory'), j.get('eft1_trajectory')
        require(isinstance(uv, list) and isinstance(eft, list) and len(uv) >= 2 and len(eft) >= 2,
                f'{model}: absent scale-resolved UV/EFT1 samples')
        for stage, rows in (('UV', uv), ('EFT1', eft)):
            scales = [number(row.get('mu_gev')) for row in rows]
            require(all(x is not None and x > 0 for x in scales), f'{model}: invalid {stage} scales')
            require(all(x > y for x, y in zip(scales, scales[1:])),
                    f'{model}: {stage} scales are not descending')
            for row in rows:
                require(number(row.get('mS1Sq_gev2')) is not None and
                        number(row.get('mS2Sq_gev2')) is not None,
                        f'{model}: missing scalar mass-squared trajectory')
        # Serialize only measured saved scale samples; do not invent final C5 curves.
        samples = []
        for stage, rows in (('UV', uv), ('EFT1', eft)):
            for row in rows:
                extra = row.get('extra_quartics') or {}
                require(isinstance(extra, dict), f'{model}: non-object extra quartics')
                samples.append({
                    'stage':stage, 'mu':number(row['mu_gev']),
                    'mS1Sq':number(row['mS1Sq_gev2']),
                    'mS2Sq':number(row['mS2Sq_gev2']),
                    'mu2':row.get('mu2_gev2'), 'quartics':extra,
                })
        data[model] = {
            'model':model, 'quadratic_kind':j.get('quadratic_kind', 'unknown'),
            'samples':samples, 'input_hashes':j.get('input_sha256', {}),
            'settings':j.get('settings', {}),
            'scalar_min_m2':number((j.get('scalar_endpoint_spectrum') or {}).get('minimum_mass_squared_gev2')),
            'spectrum_admissible':j.get('scalar_endpoint_physically_admissible'),
        }
        paths[str(path)] = source_hash(path)
    require(bool(data), 'no saved continuous scalar trajectories found')
    return data, paths


def convert_row(model: str, source: str, row: dict) -> dict | None:
    if source == 'baseline':
        masses = [row.get(k) for k in ('largest_neutrino_mass_ev','middle_neutrino_mass_ev','smallest_neutrino_mass_ev')]
    elif source == 'self_massdep':
        if row.get('matching_category') != 'self-bilinear-one-scale-massdependent-conditional':
            return None
        masses = [row.get(k) for k in ('mass1_ev_descending','mass2_ev_descending','mass3_ev_descending')]
    elif source == 'direct':
        if row.get('direct_wilson_source_evolved') not in (True, 'True', 'true', '1'):
            return None
        masses = [row.get(k) for k in ('max_mnu_ev','middle_mnu_ev','min_mnu_ev')]
    elif source == 'extra':
        masses = [row.get(k) for k in ('updated_largest_mnu_ev','updated_middle_mnu_ev','updated_smallest_mnu_ev')]
    elif source == 'coupled':
        masses = [row.get(k) for k in ('new_largest_mnu_ev','new_middle_mnu_ev','new_smallest_mnu_ev')]
    else:
        raise ValueError('unknown source')
    vals = [number(x) for x in masses]
    if any(v is None or v < 0 for v in vals):
        return None
    require(vals[0] >= vals[1] >= vals[2], f'{model} {source}: masses not descending')
    item = {'model':model, 'source':source, 'status':str(row.get('status','')),
            'mass_hi_ev':vals[0], 'mass_mid_ev':vals[1], 'mass_lo_ev':vals[2],
            'gap_small_mid_sq_ev2':vals[1]**2-vals[2]**2,
            'gap_small_large_sq_ev2':vals[0]**2-vals[2]**2}
    for key in ('delta_m21_sq_ev2','delta_m3l_sq_ev2', 'sin2_theta12','sin2_theta13','sin2_theta23'):
        item[key] = number(row.get(key))
    # The B/C selected and coupled sources only save singular masses. Keep their
    # observables unavailable; do not copy mixing angles from a different scheme.
    item['scheme_facts'] = {}
    for key in ('anomalous_over_normal_c5','relative_C5_change_vs_frozen_lambda',
                'relative_extra_C5_vs_direct_at_scalar_threshold',
                'relative_extra_C5_vs_frozen_extra','relative_C5_change_vs_previous_selected',
                'direct_source_parity','extra_source_parity'):
        value = number(row.get(key))
        if value is not None:
            item['scheme_facts'][key] = value
    return item


def load_variants(root: Path, models: set[str]) -> tuple[list[dict],dict[str,str],list[str]]:
    all_items = []
    hashes = {}
    availability = []
    for code, rel, title, note in SOURCES:
        path = root / rel
        if not path.is_file():
            availability.append(f'{code}: not found ({rel})')
            if code == 'baseline':
                require(False, f'benchmark mass table is required: {path}')
            continue
        hashes[str(path)] = source_hash(path)
        raw = read_csv(path)
        found = set()
        for row in raw:
            model = row.get('model')
            if model not in models:
                continue
            require(model not in found, f'duplicate {model} in {path}')
            found.add(model)
            item = convert_row(model, code, row)
            if item:
                item.update(source_title=title, source_note=note, source_relative_path=rel)
                all_items.append(item)
        if code == 'baseline':
            require(models <= {r['model'] for r in all_items if r['source'] == 'baseline'},
                    f'missing baseline neutrino masses for {sorted(models - found)}')
    all_items.sort(key=lambda x:(MODELS.index(x['model']),SOURCE_ORDER[x['source']]))
    return all_items, hashes, availability


def choose_preferred(items: list[dict], model: str) -> dict:
    candidates = {x['source']:x for x in items if x['model'] == model}
    # Most specific supported approximation; NOT an additive matched total.
    priority = ('coupled','extra','self_massdep','direct','baseline') if model in ('B_-1','C_-1') else (
        ('self_massdep','direct','baseline') if model in ('A_+0','A_-2','D_+0','D_-2','E_+0','E_-2')
        else ('direct','baseline'))
    for name in priority:
        if name in candidates:
            return dict(candidates[name], selection_reason='Highest available model-specific candidate; schemes not added')
    raise ValueError('missing mass result for '+model)


def load_parameter_scans(path: Path | None, models: set[str]) -> tuple[list[dict],str|None]:
    if path is None:
        return [], 'No general-scalar T1/T2/T3 response scan was supplied; ovals are intentionally absent.'
    require(path.is_file(), f'sensitivity CSV not found: {path}')
    rows = read_csv(path)
    require(bool(rows), 'empty sensitivity CSV')
    groups = {}
    data = []
    for row in rows:
        model = row.get('model')
        if model not in models:
            continue
        direction = str(row.get('direction',''))
        require(direction in ('T1','T2','T3'), 'unknown sensitivity direction')
        gamma = number(row.get('gamma'))
        require(gamma is not None and gamma > 0, 'nonpositive gamma')
        record = {'model':model,'direction':direction,'gamma':gamma}
        for key in OBS:
            record[key] = number(row.get(key))
        require(sum(record[k] is not None for k in OBS) >= 2,
                'each sensitivity row must contain >=2 observables')
        data.append(record)
        groups.setdefault((model,direction),[]).append(gamma)
    require(bool(data), 'no matching models in sensitivity CSV')
    for key, values in groups.items():
        require(len(values) >= 3 and len(set(values)) == len(values),
                f'{key}: need >=3 distinct gamma points')
        require(any(math.isclose(g,1.0,rel_tol=1e-9,abs_tol=1e-12) for g in values),
                f'{key}: missing gamma=1 reference')
    return data, None


def export_csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer = csv.DictWriter(stream,fieldnames=columns,extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _script_json(data: Any) -> str:
    return json.dumps(data,ensure_ascii=False,allow_nan=False,separators=(',',':')).replace('<','\\u003c')


CSS = '''
:root{color-scheme:dark;--bg:#0d1420;--panel:#172235;--panel2:#1f2d43;--line:#31425a;--text:#e7eef7;--muted:#a5b5cb;--cyan:#70d2e9;--amber:#efc37a;--violet:#b6a7ee}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--text);font:14px/1.55 system-ui,-apple-system,Segoe UI,sans-serif}
main{max-width:1460px;margin:0 auto;padding:28px 28px 64px}header{border-bottom:1px solid var(--line);background:#101b2c;padding:18px max(28px,calc((100vw - 1404px)/2))}
nav{display:flex;gap:20px;flex-wrap:wrap;align-items:center}nav a{color:var(--muted);text-decoration:none;font-weight:600}nav a.active,nav a:hover{color:var(--cyan)}
h1{font-size:28px;letter-spacing:-.025em;margin:8px 0 3px}h2{font-size:19px;margin:0 0 16px}h3{font-size:15px;margin:12px 0 8px}
p{color:var(--muted)}.eyebrow{color:var(--cyan);font-size:11px;letter-spacing:.1em;text-transform:uppercase;font-weight:700}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:14px;margin:22px 0}.metric,.panel{background:var(--panel);border:1px solid var(--line);border-radius:11px;padding:20px}
.metric strong{display:block;font-size:28px;color:var(--text);margin:6px 0 0}.metric small{color:var(--muted)}
.panel{margin-top:20px}.controls{display:flex;gap:12px;flex-wrap:wrap;align-items:end;margin-bottom:14px}.control{display:flex;flex-direction:column;gap:4px;min-width:175px}
label{font-size:12px;color:var(--muted)}select,input,button{border-radius:7px;background:var(--panel2);border:1px solid var(--line);color:var(--text);font:inherit;padding:9px 11px}
button{cursor:pointer}button:hover{border-color:var(--cyan)}input[type=checkbox]{accent-color:var(--cyan);width:16px;height:16px}
.chart{background:#101b2b;border:1px solid var(--line);border-radius:8px;min-height:380px;position:relative;overflow:hidden}
.chart svg{width:100%;height:440px;display:block}.legend{display:flex;gap:12px;flex-wrap:wrap;margin:12px 0}.legend label{color:var(--text);display:flex;align-items:center;gap:6px}.legend input{margin:0}
.note{border-left:3px solid var(--amber);padding:10px 14px;background:rgba(239,195,122,.055);color:var(--muted);margin:14px 0}
.table-scroll{overflow:auto;max-height:72vh}.table-scroll table{white-space:nowrap}table{width:100%;border-collapse:collapse;font-size:12px}th{background:#213149;color:#c9d7e8;position:sticky;top:0;z-index:1;text-align:left}td,th{padding:10px 12px;border-bottom:1px solid var(--line)}tr:hover td{background:#20314a}td.num{font-variant-numeric:tabular-nums;text-align:right}a{color:var(--cyan)}.pill{display:inline-block;background:#263a54;color:#d3e5fa;padding:2px 8px;border-radius:20px;font-size:11px}
footer{margin-top:25px;color:var(--muted);font-size:12px}.empty{color:var(--amber);padding:20px}.sub{font-size:12px;color:var(--muted)}
@media(max-width:650px){main{padding:18px 13px 40px}header{padding:14px}h1{font-size:22px}.panel{padding:12px}.chart svg{height:360px}}
'''
JS_COMMON = '''
const $=id=>document.getElementById(id);
const fmt=(v,n=5)=>v===null||v===undefined||!Number.isFinite(Number(v))?'—':Number(v).toPrecision(n);
const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const palette=['#79dce7','#f0b77b','#bbabed','#8cdb9c','#eb8ba0','#e7d576','#96aff0','#87cfbb','#dc9fc8','#c5cedd'];
function pathLine(points){return points.map((p,i)=>(i?'L':'M')+p[0].toFixed(2)+','+p[1].toFixed(2)).join(' ')}
function svgDraw(id,series,opts={}){
 const host=$(id);const w=1050,h=430,L=78,R=25,T=25,B=58;
 const all=series.flatMap(s=>s.points||[]).filter(p=>Number.isFinite(p[0])&&Number.isFinite(p[1]));
 if(!all.length){host.innerHTML='<div class="empty">No saved numerical points for this selection.</div>';return;}
 let xx=all.map(p=>p[0]), yy=all.map(p=>p[1]);let xmin=Math.min(...xx),xmax=Math.max(...xx),ymin=Math.min(...yy),ymax=Math.max(...yy);
 if(xmin===xmax){xmin-=1;xmax+=1;}if(ymin===ymax){let d=Math.max(Math.abs(ymin)*.1,1e-20);ymin-=d;ymax+=d;}
 const pad=(ymax-ymin)*.07;ymin-=pad;ymax+=pad;
 const X=x=>L+(x-xmin)/(xmax-xmin)*(w-L-R),Y=y=>h-B-(y-ymin)/(ymax-ymin)*(h-T-B);
 let parts=[`<svg viewBox="0 0 ${w} ${h}" role="img" aria-label="${esc(opts.title||'Interactive scientific plot')}">`];
 for(let k=0;k<=5;k++){let y=ymin+(ymax-ymin)*k,py=Y(y);parts.push(`<line x1="${L}" x2="${w-R}" y1="${py}" y2="${py}" stroke="#2a3a4f"/><text x="${L-9}" y="${py+4}" fill="#b8c9dc" font-size="12" text-anchor="end">${esc(opts.yFmt?opts.yFmt(y):fmt(y,3))}</text>`)}
 for(let k=0;k<=6;k++){let x=xmin+(xmax-xmin)*k,px=X(x);parts.push(`<line x1="${px}" x2="${px}" y1="${T}" y2="${h-B}" stroke="#243349"/><text x="${px}" y="${h-B+23}" fill="#b8c9dc" font-size="12" text-anchor="middle">${esc(opts.xFmt?opts.xFmt(x):fmt(x,3))}</text>`)}
 if(Number.isFinite(opts.threshold)&&opts.threshold>xmin&&opts.threshold<xmax){let tx=X(opts.threshold);parts.push(`<line x1="${tx}" x2="${tx}" y1="${T}" y2="${h-B}" stroke="#f1c479" stroke-dasharray="5 5"/><text x="${tx+7}" y="${T+14}" fill="#f1c479" font-size="12">${esc(opts.thresholdLabel||'F threshold')}</text>`)}
 series.forEach((s,i)=>{let color=s.color||palette[i%palette.length],pts=(s.points||[]).filter(p=>Number.isFinite(p[0])&&Number.isFinite(p[1]));if(!pts.length)return;
 if(opts.sort)pts=pts.slice().sort((a,b)=>a[0]-b[0]);let mapped=pts.map(p=>[X(p[0]),Y(p[1])]);
 if(mapped.length>1)parts.push(`<path d="${pathLine(mapped)}" fill="none" stroke="${color}" stroke-width="2.7" opacity=".93"/>`);
 if(opts.showEllipse && mapped.length>=3){
   const n=mapped.length,mx=mapped.reduce((a,p)=>a+p[0],0)/n,my=mapped.reduce((a,p)=>a+p[1],0)/n;
   let cxx=0,cyy=0,cxy=0;for(const p of mapped){let dx=p[0]-mx,dy=p[1]-my;cxx+=dx*dx;cyy+=dy*dy;cxy+=dx*dy;}
   cxx/=n-1;cyy/=n-1;cxy/=n-1;
   const trace=(cxx+cyy)/2,disc=Math.sqrt(Math.max(0,((cxx-cyy)/2)**2+cxy*cxy));
   const rx=Math.sqrt(Math.max(0,trace+disc)),ry=Math.sqrt(Math.max(0,trace-disc));
   const a=0.5*Math.atan2(2*cxy,cxx-cyy)*180/Math.PI;
   if(Number.isFinite(rx)&&Number.isFinite(ry))parts.push(`<ellipse cx="${mx}" cy="${my}" rx="${Math.max(rx,.7)}" ry="${Math.max(ry,.7)}" transform="rotate(${a} ${mx} ${my})" stroke="${color}" fill="none" stroke-width="1.5" stroke-dasharray="7 4" opacity=".7"><title>Descriptive one-sigma sample covariance (not a confidence interval)</title></ellipse>`);
 }
 pts.forEach((p,k)=>parts.push(`<circle cx="${mapped[k][0]}" cy="${mapped[k][1]}" r="3.3" fill="${color}"><title>${esc(s.name)}: ${esc(opts.xFmt?opts.xFmt(p[0]):fmt(p[0]))}, ${esc(fmt(p[1]))}${p[2]?' — '+esc(p[2]):''}</title></circle>`));
 });
 parts.push(`<text x="${(L+w-R)/2}" y="${h-8}" text-anchor="middle" fill="#d5e0ed" font-size="13">${esc(opts.xTitle||'x')}</text>`);
 parts.push('</svg>');host.innerHTML=parts.join('');
}
function fillSelect(el,options,selected){el.innerHTML=options.map(([k,v])=>`<option value="${esc(k)}">${esc(v)}</option>`).join('');if(selected!==undefined)el.value=selected;}
'''


def html_page(title: str, active: str, content: str, data: dict | None = None, script: str = '') -> str:
    nav = [('index.html','Overview'),('interactive_rge.html','Interactive RGE'),
           ('sensitivity_report.html','Sensitivity report'),('all_results.html','All results')]
    links = ''.join(f'<a href="{file}" class="{"active" if file == active else ""}">{label}</a>' for file,label in nav)
    encoded = f'<script type="application/json" id="payload">{_script_json(data)}</script>' if data is not None else ''
    return ('<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" '
            'content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><style>'+CSS+
            '</style></head><body><header><nav>'+links+'</nav></header><main>'+content+
            '</main>'+encoded+'<script>'+JS_COMMON+script+'</script></body></html>')


def build_index(data: dict, preferred: list[dict], variants: list[dict], availability: list[str]) -> str:
    categories = {k:sum(x['source']==k for x in preferred) for k,*_ in SOURCES}
    examples = ''.join(f'<tr><td><a href="all_results.html?model={html.escape(x["model"])}">{html.escape(x["model"])}</a></td>'
                       f'<td><span class="pill">{html.escape(x["source"])}</span></td><td class="num">{x["mass_hi_ev"]:.6g}</td>'
                       f'<td class="num">{x["mass_mid_ev"]:.6g}</td><td class="num">{x["mass_lo_ev"]:.6g}</td></tr>' for x in preferred)
    note = ('These are candidate results from different matching schemes, not a single certified full-scalar EFT. '
            'B/C coupled, self-bilinear one-scale, direct sequential and legacy control results are intentionally kept distinct. '
            'No uncomputed scalar-threshold or full-flavour effects are inserted.')
    return (f'<div class="eyebrow">T3 / general-scalar campaign</div><h1>Numerical display · {html.escape(data["run_id"])}</h1>'
            f'<p>Offline interactive displays reconstructed from saved full-scalar trajectories and matching summaries. No physics rerun.</p>'
            f'<div class="grid"><div class="metric"><small>Saved scalar trajectories</small><strong>{len(data["trajectories"])} / 16</strong></div>'
            f'<div class="metric"><small>Selected benchmark rows</small><strong>{len(preferred)}</strong></div>'
            f'<div class="metric"><small>Saved matching variants</small><strong>{len(variants)}</strong></div>'
            f'<div class="metric"><small>Coupled-Wilson model results</small><strong>{categories.get("coupled",0)}</strong></div></div>'
            f'<div class="note">{html.escape(note)}</div>'
            '<div class="panel"><h2>Current 16-model comparison</h2><p class="sub">Masses in eV, descending. '
            'The displayed source is the most model-specific available candidate, not the sum of all sources.</p>'
            '<div class="table-scroll"><table><thead><tr><th>Model</th><th>Selected source</th><th>Largest mass</th><th>Middle mass</th><th>Smallest mass</th></tr></thead><tbody>'+
            examples+'</tbody></table></div><p><a href="all_results.html">Open sortable complete table</a> · '
            '<a href="all_results.csv">Download selected results CSV</a> · <a href="matching_variants.csv">Download all saved variants</a></p></div>'
            '<div class="panel"><h2>Reports</h2><p><a href="interactive_rge.html">Interactive RGE</a> — UV/EFT1 mass-squared, holomorphic quadratic and new-quartic scale trajectories.</p>'
            '<p><a href="sensitivity_report.html">Sensitivity</a> — saved matching-approximation shifts; parameter-response points and descriptive ovals only if supplied.</p></div>'
            + ('<div class="panel"><h2>Unavailable inputs</h2><p>'+html.escape('; '.join(availability))+'</p></div>' if availability else '')
            +'<footer>Source manifest: <a href="provenance.json">provenance.json</a>. All numerics are saved candidates; complete sequential matching is not certified.</footer>')


def build_rge(trajectories: dict) -> str:
    body = '''<div class="eyebrow">Scale-resolved saved states</div><h1>Interactive RGE report</h1>
<p>Actual saved UV → fermion threshold → EFT1 scalar evolution. Horizontal axis: log₁₀ μ [GeV].
Curves are reconstructed from the stored samples, not interpolated by an RGE solver.</p>
<div class="panel"><div class="controls"><div class="control"><label for="model">Model</label><select id="model"></select></div>
<div class="control"><label for="family">Quantity</label><select id="family"><option value="masses">Scalar mass-squared [GeV²]</option><option value="quartics">New quartic couplings</option><option value="mu2">Holomorphic quadratic |μ²| [GeV²]</option></select></div>
<div class="control"><label for="component">Complex component</label><select id="component"><option value="abs">Absolute magnitude</option><option value="real">Real part</option><option value="imag">Imaginary part</option></select></div>
<label><input type="checkbox" id="uv" checked> UV</label><label><input type="checkbox" id="eft" checked> EFT1</label>
</div><div class="legend" id="legend"></div><div class="chart" id="rgechart"></div>
<p class="sub" id="details"></p></div>
<div class="note">This plot contains scalar renormalisable running only. Neither a full scale-resolved C₅ nor an SM+Weinberg trajectory is present in the saved continuous_trajectory.json source. It is not reconstructed from endpoint masses.</div>
'''
    js = '''
const trajectories=JSON.parse($('payload').textContent);
fillSelect($('model'),Object.keys(trajectories).map(k=>[k,k]),'B_-1');
let enabled={};
function component(z,mode){if(!z||typeof z!=='object')return null;const re=Number(z.real),im=Number(z.imag);if(!Number.isFinite(re)||!Number.isFinite(im))return null;return mode==='real'?re:mode==='imag'?im:Math.hypot(re,im)}
function render(){const name=$('model').value,kind=$('family').value,mode=$('component').value,d=trajectories[name];
 const keepStages=new Set([$('uv').checked?'UV':null,$('eft').checked?'EFT1':null]);
 const rows=d.samples.filter(r=>keepStages.has(r.stage));let keys=[];
 if(kind==='masses')keys=['mS1Sq','mS2Sq'];else if(kind==='mu2')keys=['mu2'];else keys=[...new Set(rows.flatMap(r=>Object.keys(r.quartics||{})))].sort();
 const legend=$('legend');legend.innerHTML=keys.map((k,i)=>`<label><input type="checkbox" data-quantity="${esc(k)}" ${enabled[k]===false?'':'checked'}><span style="color:${palette[i%palette.length]}">●</span>${esc(k)}</label>`).join('');
 legend.querySelectorAll('input').forEach(x=>x.addEventListener('change',()=>{enabled[x.dataset.quantity]=x.checked;draw()}));
 function draw(){let out=keys.filter(k=>enabled[k]!==false).map((key,i)=>({name:key,color:palette[keys.indexOf(key)%palette.length],points:rows.map(r=>{
 let val=kind==='masses'?r[key]:kind==='mu2'?component(r.mu2,mode):component((r.quartics||{})[key],mode);
 return [Math.log10(r.mu),val,r.stage+' / '+fmt(r.mu,5)+' GeV'];}).filter(p=>Number.isFinite(p[1]))}));
 svgDraw('rgechart',out,{sort:true,xTitle:'log₁₀ μ [GeV]',threshold:Math.log10(Number(d.settings.fermion_threshold_gev||1e5)),thresholdLabel:'F threshold',
 xFmt:x=>x.toFixed(2)});}
 draw();$('details').textContent=name+' · '+d.quadratic_kind+' · '+d.samples.filter(x=>x.stage==='UV').length+' UV points + '+d.samples.filter(x=>x.stage==='EFT1').length+' EFT1 points';
}
['model','family','component','uv','eft'].forEach(id=>$(id).addEventListener('change',render));render();
'''
    return html_page('General scalar · interactive RGE', 'interactive_rge.html',body,trajectories,js)


def build_sensitivity(preferred: list[dict], variants: list[dict], scans: list[dict], message: str|None) -> str:
    body = '''<div class="eyebrow">Numerical response and approximation dependence</div><h1>Sensitivity report</h1>
<p>The top panel measures the shift between saved matching approximations for a fixed benchmark. This is <strong>not</strong> a UV parameter derivative, statistical uncertainty, or an additive EFT correction.</p>
<div class="panel"><h2>Matching-scheme sensitivity</h2><div class="controls"><div class="control"><label for="model">Model</label><select id="model"></select></div>
<div class="control"><label for="obs">Observable</label><select id="obs"></select></div>
</div><div class="chart" id="schemeChart"></div><p class="sub" id="schemeDetails"></p></div>
<div class="panel"><h2>T1 / T2 / T3 parameter-response display</h2><p id="scanMessage" class="note"></p>
<div id="scanControls"><div class="controls"><div class="control"><label for="scanModel">Model</label><select id="scanModel"></select></div>
<div class="control"><label for="direction">Direction</label><select id="direction"><option>T1</option><option>T2</option><option>T3</option></select></div>
<div class="control"><label for="xobs">x observable</label><select id="xobs"></select></div><div class="control"><label for="yobs">y observable</label><select id="yobs"></select></div></div>
<label><input type="checkbox" id="showOval" checked> Show descriptive sample covariance oval (1σ, not a confidence region)</label>
<div class="chart" id="scanChart"></div><p class="sub">Each point is a saved γ sample. No points are synthesized. Dashed ovals are descriptive sample covariance, not confidence contours. Genuine general-scalar sensitivity scans must use the same UV/RGE/matching convention as the benchmark being compared.</p></div></div>
<div class="note">Do not interpret the 16 different models or different matching schemes as a T1/T2/T3 covariance sample. Historical reduced-scalar sensitivity ovals are not general-scalar ovals.</div>'''
    if message:
        body = body.replace('<p id="scanMessage" class="note"></p>', '<p id="scanMessage" class="note">'+html.escape(message)+'</p>')
    data={'preferred':preferred,'variants':variants,'scans':scans,'obs_labels':OBS_LABELS}
    js='''
const dat=JSON.parse($('payload').textContent),models=[...new Set(dat.preferred.map(x=>x.model))],lab=dat.obs_labels;
fillSelect($('model'),models.map(x=>[x,x]),models.includes('B_-1')?'B_-1':models[0]);fillSelect($('obs'),Object.entries(lab),'mass_hi_ev');
function schemer(){let model=$('model').value,obs=$('obs').value,entries=dat.variants.filter(r=>r.model===model&&Number.isFinite(r[obs]));
 const baseline=entries.find(r=>r.source==='baseline');if(!baseline){$('schemeChart').innerHTML='<div class="empty">No benchmark baseline</div>';return}
 const val=baseline[obs];let values=entries.map((r,i)=>({name:r.source,points:[[i,(r[obs]-val)/(Math.abs(val)||1)*100,r.source+' / '+fmt(r[obs])]],color:palette[i%palette.length]}));
 svgDraw('schemeChart',values,{xTitle:'Scheme order (baseline → more specific)',xFmt:x=>{let i=Math.round(x);return Math.abs(x-i)<.1?(entries[i]?.source||''):''},yFmt:y=>y.toFixed(1)+'%'});
 $('schemeDetails').textContent='Signed shift = 100 × (scheme observable − saved baseline observable) / |saved baseline observable|. Baseline: '+fmt(val)+'; '+entries.length+' saved schemes. Separate candidate approximations, not measurements of UV perturbation.';
}
$('model').addEventListener('change',schemer);$('obs').addEventListener('change',schemer);schemer();
const groups=[...new Set(dat.scans.map(r=>r.model))];if(!groups.length){$('scanControls').style.display='none';}
else{ $('scanMessage').textContent='Saved general-scalar parameter-response points supplied: '+dat.scans.length+'. Displaying measured γ response curves.';
 fillSelect($('scanModel'),groups.map(k=>[k,k])); fillSelect($('xobs'),Object.entries(lab),'sin2_theta13');fillSelect($('yobs'),Object.entries(lab),'delta_m3l_sq_ev2');
 function scanRender(){let m=$('scanModel').value,d=$('direction').value,x=$('xobs').value,y=$('yobs').value;
 let points=dat.scans.filter(r=>r.model===m&&r.direction===d&&Number.isFinite(r[x])&&Number.isFinite(r[y])).sort((a,b)=>a.gamma-b.gamma);
 svgDraw('scanChart',[{name:d+' / '+m,points:points.map(r=>[r[x],r[y],'γ = '+fmt(r.gamma)])}],{xTitle:lab[x],showEllipse:$('showOval').checked});}
 ['scanModel','direction','xobs','yobs','showOval'].forEach(id=>$(id).addEventListener('change',scanRender));scanRender();}
'''
    return html_page('General scalar · sensitivity','sensitivity_report.html',body,data,js)


def build_results(preferred: list[dict], variants: list[dict]) -> str:
    body = '''<div class="eyebrow">Saved numerical outputs</div><h1>All general-scalar results</h1>
<p>Sixteen chosen candidate benchmark rows and every separately saved matching variant. Empty cells mean the corresponding quantity was not calculated in that scheme.</p>
<div class="panel"><div class="controls"><div class="control"><label for="modelFilter">Model</label><select id="modelFilter"></select></div>
<div class="control"><label for="kindFilter">View</label><select id="kindFilter"><option value="preferred">One candidate per model</option><option value="variants">All matching variants</option></select></div>
<div class="control"><label for="sourceFilter">Source</label><select id="sourceFilter"></select></div>
<div class="control"><label for="sort">Sort</label><select id="sort"><option value="model">Model</option><option value="mass_hi_ev">Largest mass</option><option value="mass_mid_ev">Middle mass</option><option value="mass_lo_ev">Smallest mass</option><option value="delta_m21_sq_ev2">Δm²₂₁</option><option value="delta_m3l_sq_ev2">Δm²₃ℓ</option></select></div></div>
<p class="sub"><a href="all_results.csv">Selected 16-row CSV</a> · <a href="matching_variants.csv">All variants CSV</a></p>
<div class="table-scroll" id="resultsTable"></div><p class="sub" id="counts"></p></div>
<div class="note">Mass values are stored in descending order. Mass ordering and observables must be interpreted according to the original matching calculation. Mixing angles are intentionally absent for the B/C selected-source and coupled stages because their saved summaries contain masses only. Values are not silently inherited from earlier stages.</div>'''
    data={'preferred':preferred,'variants':variants,'labels':OBS_LABELS}
    js='''
const data=JSON.parse($('payload').textContent),names=[...new Set(data.preferred.map(r=>r.model))];
fillSelect($('modelFilter'),[['all','All 16 models'],...names.map(k=>[k,k])]);
fillSelect($('sourceFilter'),[['all','All sources'],...([...new Set(data.variants.map(r=>r.source))].map(k=>[k,k]))]);
const incoming=new URLSearchParams(location.search).get('model');if(names.includes(incoming))$('modelFilter').value=incoming;
function drawResults(){let use=$('kindFilter').value==='variants'?data.variants:data.preferred;
 let m=$('modelFilter').value,s=$('sourceFilter').value,key=$('sort').value;
 let rows=use.filter(r=>(m==='all'||r.model===m)&&(s==='all'||r.source===s)).slice();
 if(key!=='model')rows.sort((a,b)=>((b[key]??-Infinity)-(a[key]??-Infinity)));else rows.sort((a,b)=>names.indexOf(a.model)-names.indexOf(b.model));
 const col=['model','source',...Object.keys(data.labels)];
 let head='<thead><tr>'+col.map(k=>'<th>'+esc(data.labels[k]||k)+'</th>').join('')+'</tr></thead>';
 let body='<tbody>'+rows.map(r=>'<tr>'+col.map(k=>'<td'+(k!=='model'&&k!=='source'?' class="num"':'')+'>'+(k==='source'?'<span class="pill">'+esc(r[k])+'</span>':typeof r[k]==='number'?fmt(r[k],7):esc(r[k]??'—'))+'</td>').join('')+'</tr>').join('')+'</tbody>';
 $('resultsTable').innerHTML='<table>'+head+body+'</table>';$('counts').textContent=rows.length+' records · all masses in eV; Δm² in eV².';
}
['modelFilter','sourceFilter','kindFilter','sort'].forEach(id=>$(id).addEventListener('change',drawResults));drawResults();
'''
    return html_page('General scalar · all results','all_results.html',body,data,js)


def build(input_root: Path, output: Path, *, run_id: str,
          sensitivity_csv: Path | None = None, allow_partial: bool = False) -> dict:
    require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,48}',run_id)), 'invalid run ID')
    input_root = input_root.resolve()
    output = output.resolve()
    require(input_root.is_dir(), f'missing run directory {input_root}')
    # Avoid placing the report over numerical inputs.
    require(output != input_root and input_root not in output.parents
            and output not in input_root.parents and output != ROOT.resolve(),
            'report output must not contain or be contained by the numerical input, or be the repository root')
    trajectories, src = read_continuous(input_root, partial=allow_partial)
    models = set(trajectories)
    variants, csv_src, availability = load_variants(input_root,models)
    src.update(csv_src)
    preferred = [choose_preferred(variants,m) for m in MODELS if m in models]
    scans, sensitivity_note = load_parameter_scans(sensitivity_csv,models)
    if sensitivity_csv is not None:
        src[str(sensitivity_csv.resolve())] = source_hash(sensitivity_csv)
    # Validate all required sources BEFORE writing output.
    output.mkdir(parents=True,exist_ok=True)
    columns = ['model','source','source_title','status',*OBS,
               'gap_small_mid_sq_ev2','gap_small_large_sq_ev2']
    export_csv(output/'all_results.csv', preferred,columns)
    export_csv(output/'matching_variants.csv',variants,columns)
    if scans:
        export_csv(output/'parameter_sensitivity_points.csv',scans,['model','direction','gamma',*OBS])
    content = {'run_id':run_id,'trajectories':trajectories}
    for filename, text in (
        ('index.html',build_index(content,preferred,variants,availability)),
        ('interactive_rge.html',build_rge(trajectories)),
        ('sensitivity_report.html',build_sensitivity(preferred,variants,scans,sensitivity_note)),
        ('all_results.html',build_results(preferred,variants)),
    ):
        (output/filename).write_text(text,encoding='utf-8')
    counts={source:sum(r['source']==source for r in preferred) for source,*_ in SOURCES}
    manifest = {
        'format':'general-scalar-offline-interactive-display-v1',
        'generated_utc':datetime.now(timezone.utc).isoformat(),
        'run_id':run_id, 'model_count':len(preferred),'matching_variant_count':len(variants),
        'selected_matching_sources':counts,
        'general_scalar_parameter_sensitivity_scan_loaded':bool(scans),
        'parameter_scan_record_count':len(scans),
        'sensitivity_note':sensitivity_note,
        'rge_curves_from_saved_uv_and_eft1_samples':True,
        'scale_resolved_C5_saved_for_all_models':False,
        'physical_full_general_scalar_matching_certified':False,
        'matching_schemes_added_to_one_another':False,
        'unavailable_optional_sources':availability,
        'input_sha256':{str(Path(p).relative_to(input_root)) if Path(p).is_relative_to(input_root)
                        else str(p):v for p,v in src.items()},
        'created_files':['index.html','interactive_rge.html','sensitivity_report.html',
                         'all_results.html','all_results.csv','matching_variants.csv','provenance.json']+
                        (['parameter_sensitivity_points.csv'] if scans else []),
    }
    (output/'provenance.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return manifest


def main(argv: list[str] | None = None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',default='continuous_v1')
    p.add_argument('--input-root',type=Path,help='Override saved continuous_full_scalar/<run-id> directory')
    p.add_argument('--output',type=Path,help='Output folder (default Reports/output/display_general_scalar/<run-id>)')
    p.add_argument('--sensitivity-csv',type=Path,help='Optional actual general-scalar T1/T2/T3 gamma sample CSV')
    p.add_argument('--allow-partial',action='store_true',help='Allow a pilot subset rather than 16 models')
    args=p.parse_args(argv)
    source=args.input_root or ROOT/REL_BASE/args.run_id
    output=args.output or ROOT/'Reports/output/display_general_scalar'/args.run_id
    try:
        manifest=build(source,output,run_id=args.run_id,
                       sensitivity_csv=args.sensitivity_csv,allow_partial=args.allow_partial)
    except (OSError,ValueError,KeyError,TypeError) as error:
        p.exit(2,str(error)+'\n')
    print('General scalar display:',output)
    print('Trajectories:',manifest['model_count'],'variants:',manifest['matching_variant_count'],
          'UV sensitivity scan points:',manifest['parameter_scan_record_count'])
    print('Open:',output/'index.html')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
