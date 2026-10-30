from __future__ import annotations

"""Build self-contained HTML reports for frozen general-scalar T3 fits.

This builder consumes:
- accepted frozen fit states under output/studies/general_scalar_extended/bulk_48/
- optional saved frozen trajectory JSON files
- optional saved frozen T1/T2/T3 sensitivity CSV files

It does not rerun the physics.
"""

import argparse
import csv
import hashlib
import html
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
OUT_BASE = ROOT / 'output/studies/general_scalar_extended/bulk_48'
REPORT_BASE = ROOT / 'Reports/output/display_frozen_general_scalar'
CASES = ('small_uv_frozen_5e-9', 'large_uv_frozen_0p5')
MODELS = (
    'A_+0', 'A_-2', 'A_-4', 'B_+1', 'B_-1', 'B_-3',
    'C_+1', 'C_-1', 'C_-3', 'D_+0', 'D_+2', 'D_-2',
    'E_+0', 'E_+2', 'E_-2', 'E_-4',
)
OBS = (
    'mass_hi_ev', 'mass_mid_ev', 'mass_lo_ev', 'delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
    'sin2_theta12', 'sin2_theta13', 'sin2_theta23',
)
OBS_LABELS = {
    'mass_hi_ev': 'Largest mass [eV]',
    'mass_mid_ev': 'Middle mass [eV]',
    'mass_lo_ev': 'Smallest mass [eV]',
    'delta_m21_sq_ev2': 'Δm²21 [eV²]',
    'delta_m3l_sq_ev2': 'Δm²3ℓ [eV²]',
    'sin2_theta12': 'sin²θ12',
    'sin2_theta13': 'sin²θ13',
    'sin2_theta23': 'sin²θ23',
}
RGE_Q = {
    'y1_fro': '||Y1||F',
    'y2_fro': '||Y2||F',
    'max_abs_y': 'max |Y|',
    'effective_weinberg_proxy': 'Weinberg proxy',
    'mass_hi_ev': 'Largest mass [eV]',
    'delta_m21_sq_ev2': 'Δm²21 [eV²]',
    'delta_m3l_sq_ev2': 'Δm²3ℓ [eV²]',
}
ACC = {'Accepted', 'CandidateConstrainedAccepted', 'ExistingConditionalAccepted'}


def finite(x: Any) -> float | None:
    try:
        if x in (None, ''):
            return None
        x = float(x)
        return x if math.isfinite(x) else None
    except Exception:
        return None


def safe(model: str) -> str:
    return model.replace('+', 'p').replace('-', 'm')


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_script(x: Any) -> str:
    return json.dumps(x, ensure_ascii=False, separators=(',', ':'), allow_nan=False).replace('<', '\\u003c')


def _extract_masses(pred: dict[str, Any], cfg_path: Path) -> tuple[float | None, float | None, float | None]:
    masses = [
        finite(pred.get('largest_neutrino_mass_ev') or pred.get('mass_hi_ev')),
        finite(pred.get('middle_neutrino_mass_ev') or pred.get('mass_mid_ev')),
        finite(pred.get('smallest_neutrino_mass_ev') or pred.get('mass_lo_ev')),
    ]
    if any(v is None for v in masses) and cfg_path.is_file():
        cfg = read_json(cfg_path)
        vals = [
            finite(cfg.get('largest_neutrino_mass_ev')),
            finite(cfg.get('middle_neutrino_mass_ev')),
            finite(cfg.get('smallest_neutrino_mass_ev')),
        ]
        for i, v in enumerate(vals):
            if masses[i] is None:
                masses[i] = v
    present = [v for v in masses if v is not None]
    if not present:
        return None, None, None
    present.sort(reverse=True)
    while len(present) < 3:
        present.append(None)
    return present[0], present[1], present[2]


def summary_rows(run_ids: list[str]) -> tuple[list[dict[str, Any]], dict[str, str]]:
    rows: list[dict[str, Any]] = []
    provenance: dict[str, str] = {}
    for case in CASES:
        for model in MODELS:
            for run_id in run_ids:
                matches = [OUT_BASE / run_id / case / safe(model) / 'state.json']
                matches += sorted(OUT_BASE.glob(f'{run_id}_*/{case}/{safe(model)}/state.json'))
                selected = next((p for p in matches if p.is_file() and read_json(p).get('status') in ACC), None)
                if selected is None:
                    continue
                path = selected
                doc = read_json(path)
                best = doc.get('best') or {}
                chi2 = finite(best.get('chi2'))
                if doc.get('status') not in ACC or chi2 is None:
                    continue
                rel = str(path.relative_to(ROOT))
                provenance[rel] = sha(path)
                pred = best.get('prediction') or {}
                cfg = path.parent / 'fitted_config.json'
                if cfg.is_file():
                    provenance[str(cfg.relative_to(ROOT))] = sha(cfg)
                m_hi, m_mid, m_lo = _extract_masses(pred, cfg)
                rows.append({
                    'model': model,
                    'case': case,
                    'run_id': path.parents[2].name,
                    'status': doc.get('status'),
                    'source_path': rel,
                    'chi2': chi2,
                    'mass_hi_ev': m_hi,
                    'mass_mid_ev': m_mid,
                    'mass_lo_ev': m_lo,
                    'delta_m21_sq_ev2': finite(pred.get('delta_m21_sq_ev2')),
                    'delta_m3l_sq_ev2': finite(pred.get('delta_m3l_sq_ev2')),
                    'sin2_theta12': finite(pred.get('sin2_theta12')),
                    'sin2_theta13': finite(pred.get('sin2_theta13')),
                    'sin2_theta23': finite(pred.get('sin2_theta23')),
                })
                break
    rows.sort(key=lambda r: (r['case'], r['model']))
    return rows, provenance


def load_trajectories(path_in: str | None, provenance: dict[str, str]) -> tuple[list[dict[str, Any]], str]:
    if not path_in:
        return [], 'NONE_FOUND'
    base = ROOT / path_in if not Path(path_in).is_absolute() else Path(path_in)
    if not base.is_dir():
        return [], 'MISSING_DIRECTORY'
    rows = []
    for path in sorted(base.rglob('*.json')):
        try:
            doc = read_json(path)
        except Exception:
            continue
        if doc.get('model') not in MODELS or doc.get('case') not in CASES:
            continue
        if not isinstance(doc.get('uv_trajectory'), list):
            continue
        rows.append(doc)
        provenance[str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)] = sha(path)
    return rows, ('SAVED_FROZEN_TRAJECTORIES' if rows else 'NO_VALID_TRAJECTORIES')


def load_sensitivity(path_in: str | None, provenance: dict[str, str]) -> tuple[list[dict[str, Any]], str]:
    if not path_in:
        return [], 'NONE_PROVIDED'
    path = ROOT / path_in if not Path(path_in).is_absolute() else Path(path_in)
    if not path.is_file():
        return [], 'MISSING_CSV'
    provenance[str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)] = sha(path)
    rows = []
    with path.open(newline='', encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            model, case, direction = row.get('model'), row.get('case'), row.get('direction')
            gamma = finite(row.get('gamma'))
            if model not in MODELS or case not in CASES or direction not in ('T1', 'T2', 'T3') or gamma is None:
                continue
            item = {'model': model, 'case': case, 'direction': direction, 'gamma': gamma,
                    'status': row.get('status') or ''}
            for key in ('chi2', 'y1_fro', 'y2_fro', 'max_abs_y', *OBS):
                item[key] = finite(row.get(key))
            rows.append(item)
    return rows, ('SAVED_FROZEN_GAMMA_CSV' if rows else 'EMPTY_CSV')


def build_index(report_dir: Path, rows, traj_state, sens_state) -> None:
    txt = f'''<!doctype html><html><head><meta charset="utf-8"><title>Frozen general-scalar reports</title>
<style>body{{font-family:Arial,sans-serif;max-width:1000px;margin:2rem auto;line-height:1.45}} code{{background:#f3f3f3;padding:.1rem .3rem}}</style>
</head><body><h1>Frozen general-scalar display index</h1>
<p>This report set uses accepted frozen-UV fits already saved under <code>output/studies/general_scalar_extended/bulk_48</code>.</p>
<ul>
<li><a href="interactive_rge.html">Interactive RGE report</a> — state: <b>{html.escape(traj_state)}</b></li>
<li><a href="sensitivity_report.html">T1/T2/T3 sensitivity report</a> — state: <b>{html.escape(sens_state)}</b></li>
<li><a href="all_results.html">Accepted frozen-fit table</a> — rows: <b>{len(rows)}</b></li>
</ul></body></html>'''
    (report_dir / 'index.html').write_text(txt, encoding='utf-8')


COMMON_STYLE = """
body{font-family:Arial,sans-serif;margin:0;background:#f6f7fb;color:#111}
header{padding:1rem 1.5rem;background:#1f2a44;color:#fff} main{padding:1rem 1.5rem}
.panel{background:#fff;border:1px solid #d8deea;border-radius:10px;padding:1rem;margin-bottom:1rem;box-shadow:0 1px 3px rgba(0,0,0,.04)}
.row{display:flex;gap:1rem;flex-wrap:wrap}.col{flex:1 1 300px}
label{font-weight:700;display:block;margin-bottom:.3rem} select{padding:.35rem .5rem;width:100%;max-width:28rem}
svg{width:100%;height:460px;background:#fff;border:1px solid #d8deea;border-radius:8px}
table{border-collapse:collapse;width:100%} th,td{border:1px solid #d8deea;padding:.35rem .45rem;font-size:.92rem} th{background:#eef2fa;cursor:pointer}
small.mono, code{font-family:Consolas,monospace} .legend span{display:inline-flex;align-items:center;margin-right:1rem;font-size:.92rem}
.sw{width:14px;height:3px;display:inline-block;margin-right:.4rem}.note{font-size:.92rem;color:#444}
"""

COMMON_JS = """
function fmt(x){return (x===null||x===undefined||Number.isNaN(x))?'' : Number(x).toPrecision(6)}
function clear(el){while(el.firstChild) el.removeChild(el.firstChild)}
function line(svg, x1,y1,x2,y2, stroke, width=2, dash=''){const e=document.createElementNS('http://www.w3.org/2000/svg','line');e.setAttribute('x1',x1);e.setAttribute('y1',y1);e.setAttribute('x2',x2);e.setAttribute('y2',y2);e.setAttribute('stroke',stroke);e.setAttribute('stroke-width',width);if(dash)e.setAttribute('stroke-dasharray',dash);svg.appendChild(e);return e}
function text(svg,x,y,s,size=12,anchor='middle'){const e=document.createElementNS('http://www.w3.org/2000/svg','text');e.setAttribute('x',x);e.setAttribute('y',y);e.setAttribute('font-size',size);e.setAttribute('text-anchor',anchor);e.textContent=s;svg.appendChild(e);return e}
function polyline(svg, pts, stroke, width=2){const e=document.createElementNS('http://www.w3.org/2000/svg','polyline');e.setAttribute('fill','none');e.setAttribute('stroke',stroke);e.setAttribute('stroke-width',width);e.setAttribute('points',pts.map(p=>p.join(',')).join(' '));svg.appendChild(e);return e}
function circle(svg,x,y,r,fill){const e=document.createElementNS('http://www.w3.org/2000/svg','circle');e.setAttribute('cx',x);e.setAttribute('cy',y);e.setAttribute('r',r);e.setAttribute('fill',fill);svg.appendChild(e);return e}
"""


def build_rge(report_dir: Path, trajectories: list[dict[str, Any]]) -> None:
    html_txt = f'''<!doctype html><html><head><meta charset="utf-8"><title>Frozen interactive RGE report</title><style>{COMMON_STYLE}</style></head>
<body><header><h1>Frozen general-scalar RGE report</h1><div class="note">Frozen UV input couplings are displayed as saved numerical trajectories. For frozen models the UV couplings are flat by construction.</div></header>
<main><div class="panel row"><div class="col"><label>Model / case</label><select id="modelSel"></select></div><div class="col"><label>Quantity</label><select id="qtySel"></select></div></div>
<div class="panel"><svg id="chart" viewBox="0 0 960 460"></svg><div class="legend" id="legend"></div></div>
<div class="panel"><h2>Saved source note</h2><div class="note" id="note"></div></div>
<script>const DATA={json_script(trajectories)};const QMAP={json_script(RGE_Q)};{COMMON_JS}
const colors={{uv:'#2357d4', eft:'#d45523'}};
const modelSel=document.getElementById('modelSel'), qtySel=document.getElementById('qtySel'), svg=document.getElementById('chart');
const note=document.getElementById('note'), legend=document.getElementById('legend');
const byKey=new Map(DATA.map(d=>[[d.model,d.case].join(' | '),d]));
for(const key of [...byKey.keys()]){{const o=document.createElement('option');o.value=key;o.textContent=key;modelSel.appendChild(o)}}
for(const [k,v] of Object.entries(QMAP)){{const o=document.createElement('option');o.value=k;o.textContent=v;qtySel.appendChild(o)}}
function collectPoints(tr, field){{const rows=[];for(const seg of ['uv_trajectory','eft1_trajectory']){{for(const p of (tr[seg]||[])){{if(p[field]!==undefined&&p[field]!==null&&!Number.isNaN(Number(p[field]))&&p.scale_gev) rows.push({{seg, scale:Number(p.scale_gev), value:Number(p[field])}})}}}}return rows.sort((a,b)=>a.scale-b.scale)}}
function render(){{clear(svg); legend.innerHTML=''; const tr=byKey.get(modelSel.value); if(!tr) return; const field=qtySel.value; const pts=collectPoints(tr, field); note.textContent=tr.note||''; if(!pts.length){{text(svg,480,230,'No saved points for this quantity',16); return;}} const m={{l:80,r:20,t:20,b:50,w:860,h:390}}; const xs=pts.map(p=>Math.log10(p.scale)), ys=pts.map(p=>p.value); const xmin=Math.min(...xs), xmax=Math.max(...xs), ymin=Math.min(...ys), ymax=Math.max(...ys); const dx=(xmax-xmin)||1, dy=(ymax-ymin)||1; const X=x=>m.l+(x-xmin)/dx*m.w; const Y=y=>m.t+m.h-(y-ymin)/dy*m.h; for(let i=0;i<=5;i++){{const xx=m.l+i*m.w/5; line(svg,xx,m.t,xx,m.t+m.h,'#e5e9f2',1); const val=(xmin+i*dx/5).toFixed(2); text(svg,xx,m.t+m.h+20,'1e'+val,11)}} for(let i=0;i<=5;i++){{const yy=m.t+i*m.h/5; line(svg,m.l,yy,m.l+m.w,yy,'#e5e9f2',1); const val=(ymax-i*dy/5); text(svg,m.l-8,yy+4,Number(val).toPrecision(4),11,'end')}} line(svg,m.l,m.t,m.l,m.t+m.h,'#444',1.4); line(svg,m.l,m.t+m.h,m.l+m.w,m.t+m.h,'#444',1.4); text(svg,480,450,'Scale [GeV]',13); text(svg,15,35,QMAP[field],13,'start');
 const segs=['uv_trajectory','eft1_trajectory']; for(const seg of segs){{ const arr=(tr[seg]||[]).filter(p=>p[field]!==undefined&&p[field]!==null&&!Number.isNaN(Number(p[field]))&&p.scale_gev).sort((a,b)=>a.scale_gev-b.scale_gev); if(!arr.length) continue; const col=seg==='uv_trajectory'?colors.uv:colors.eft; const points=arr.map(p=>[X(Math.log10(Number(p.scale_gev))),Y(Number(p[field]))]); polyline(svg,points,col,2.2); for(const p of points) circle(svg,p[0],p[1],2.8,col); const s=document.createElement('span'); s.innerHTML=`<span class="sw" style="background:${{col}}"></span>${{seg==='uv_trajectory'?'Frozen UV':'EFT / low-energy proxy'}}`; legend.appendChild(s); }} }}
modelSel.onchange=render; qtySel.onchange=render; render();</script></main></body></html>'''
    (report_dir / 'interactive_rge.html').write_text(html_txt, encoding='utf-8')


def build_sensitivity(report_dir: Path, rows: list[dict[str, Any]]) -> None:
    html_txt = f'''<!doctype html><html><head><meta charset="utf-8"><title>Frozen T1/T2/T3 sensitivity</title><style>{COMMON_STYLE}</style></head>
<body><header><h1>Frozen-model T1/T2/T3 sensitivity</h1><div class="note">The ovals are built from genuine saved scan points. Directions T1/T2/T3 are defined by the generator script and documented there.</div></header>
<main><div class="panel row"><div class="col"><label>Model / case</label><select id="modelSel"></select></div><div class="col"><label>X observable</label><select id="xSel"></select></div><div class="col"><label>Y observable</label><select id="ySel"></select></div></div>
<div class="panel"><svg id="oval" viewBox="0 0 960 460"></svg><div class="legend" id="legend"></div></div>
<div class="panel"><h2>Gamma scan</h2><svg id="line" viewBox="0 0 960 460"></svg></div>
<script>const ROWS={json_script(rows)}; const OBS={json_script(OBS_LABELS)}; {COMMON_JS}
const colors={{T1:'#2357d4',T2:'#d45523',T3:'#1f9d55'}}; const modelSel=document.getElementById('modelSel'), xSel=document.getElementById('xSel'), ySel=document.getElementById('ySel'); const oval=document.getElementById('oval'), lineSvg=document.getElementById('line'), legend=document.getElementById('legend');
const keys=[...new Set(ROWS.map(r=>[r.model,r.case].join(' | ')))].sort(); for(const k of keys){{const o=document.createElement('option');o.value=k;o.textContent=k;modelSel.appendChild(o)}} for(const [k,v] of Object.entries(OBS)){{const ox=document.createElement('option');ox.value=k;ox.textContent=v;xSel.appendChild(ox); const oy=document.createElement('option');oy.value=k;oy.textContent=v;ySel.appendChild(oy)}} xSel.value='sin2_theta13'; ySel.value='sin2_theta23';
function subset(){{const [m,c]=modelSel.value.split(' | '); return ROWS.filter(r=>r.model===m && r.case===c && r.status && r.status!=='');}}
function drawPlane(){{clear(oval); legend.innerHTML=''; const xs=xSel.value, ys=ySel.value; const data=subset().filter(r=>r[xs]!=null&&r[ys]!=null).sort((a,b)=>a.gamma-b.gamma); if(!data.length){{text(oval,480,230,'No saved scan points',16); return;}} const m={{l:80,r:20,t:20,b:50,w:860,h:390}}; const xvals=data.map(r=>r[xs]), yvals=data.map(r=>r[ys]); const xmin=Math.min(...xvals), xmax=Math.max(...xvals), ymin=Math.min(...yvals), ymax=Math.max(...yvals); const dx=(xmax-xmin)||1, dy=(ymax-ymin)||1; const X=x=>m.l+(x-xmin)/dx*m.w; const Y=y=>m.t+m.h-(y-ymin)/dy*m.h; for(let i=0;i<=5;i++){{const xx=m.l+i*m.w/5; line(oval,xx,m.t,xx,m.t+m.h,'#e5e9f2',1); text(oval,xx,m.t+m.h+20,Number(xmin+i*dx/5).toPrecision(4),11)}} for(let i=0;i<=5;i++){{const yy=m.t+i*m.h/5; line(oval,m.l,yy,m.l+m.w,yy,'#e5e9f2',1); text(oval,m.l-8,yy+4,Number(ymax-i*dy/5).toPrecision(4),11,'end')}} line(oval,m.l,m.t,m.l,m.t+m.h,'#444',1.4); line(oval,m.l,m.t+m.h,m.l+m.w,m.t+m.h,'#444',1.4); text(oval,480,450,OBS[xs],13); text(oval,18,35,OBS[ys],13,'start'); for(const dir of ['T1','T2','T3']){{const arr=data.filter(r=>r.direction===dir); if(!arr.length) continue; const pts=arr.map(r=>[X(r[xs]),Y(r[ys])]); polyline(oval,pts,colors[dir],2.2); for(const p of pts) circle(oval,p[0],p[1],3,colors[dir]); const s=document.createElement('span'); s.innerHTML=`<span class="sw" style="background:${{colors[dir]}}"></span>${{dir}}`; legend.appendChild(s); }} }}
function drawLine(){{clear(lineSvg); const ys='chi2'; const data=subset().filter(r=>r.gamma!=null&&r[ys]!=null).sort((a,b)=>a.gamma-b.gamma); if(!data.length){{text(lineSvg,480,230,'No saved chi² points',16); return;}} const m={{l:80,r:20,t:20,b:50,w:860,h:390}}; const xvals=data.map(r=>r.gamma), yvals=data.map(r=>r[ys]); const xmin=Math.min(...xvals), xmax=Math.max(...xvals), ymin=Math.min(...yvals), ymax=Math.max(...yvals); const dx=(xmax-xmin)||1, dy=(ymax-ymin)||1; const X=x=>m.l+(x-xmin)/dx*m.w; const Y=y=>m.t+m.h-(y-ymin)/dy*m.h; for(let i=0;i<=5;i++){{const xx=m.l+i*m.w/5; line(lineSvg,xx,m.t,xx,m.t+m.h,'#e5e9f2',1); text(lineSvg,xx,m.t+m.h+20,Number(xmin+i*dx/5).toPrecision(4),11)}} for(let i=0;i<=5;i++){{const yy=m.t+i*m.h/5; line(lineSvg,m.l,yy,m.l+m.w,yy,'#e5e9f2',1); text(lineSvg,m.l-8,yy+4,Number(ymax-i*dy/5).toPrecision(4),11,'end')}} line(lineSvg,m.l,m.t,m.l,m.t+m.h,'#444',1.4); line(lineSvg,m.l,m.t+m.h,m.l+m.w,m.t+m.h,'#444',1.4); text(lineSvg,480,450,'γ',13); text(lineSvg,18,35,'χ²',13,'start'); for(const dir of ['T1','T2','T3']){{const arr=data.filter(r=>r.direction===dir); if(!arr.length) continue; const pts=arr.map(r=>[X(r.gamma),Y(r[ys])]); polyline(lineSvg,pts,colors[dir],2.2); for(const p of pts) circle(lineSvg,p[0],p[1],3,colors[dir]); }} }}
function render(){{drawPlane(); drawLine();}} modelSel.onchange=render; xSel.onchange=render; ySel.onchange=render; render();</script></main></body></html>'''
    (report_dir / 'sensitivity_report.html').write_text(html_txt, encoding='utf-8')


def build_table(report_dir: Path, rows: list[dict[str, Any]]) -> None:
    cols = ['model', 'case', 'chi2', *OBS, 'run_id', 'status']
    head = ''.join(f'<th onclick="sortTable({i})">{html.escape(c)}</th>' for i, c in enumerate(cols))
    body = []
    for row in rows:
        body.append('<tr>' + ''.join(f'<td>{html.escape("" if row.get(c) is None else str(row.get(c)))}</td>' for c in cols) + '</tr>')
    txt = f'''<!doctype html><html><head><meta charset="utf-8"><title>Accepted frozen fits</title><style>{COMMON_STYLE}</style></head>
<body><header><h1>Accepted frozen-fit table</h1></header><main><div class="panel"><table id="results"><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table></div>
<script>
function sortTable(n){{const t=document.getElementById('results'); const b=t.tBodies[0]; const rows=[...b.rows]; const num=rows.every(r=>!Number.isNaN(Number(r.cells[n].textContent))); rows.sort((a,b)=> num ? Number(a.cells[n].textContent)-Number(b.cells[n].textContent) : a.cells[n].textContent.localeCompare(b.cells[n].textContent)); for(const r of rows) b.appendChild(r); }}
</script></main></body></html>'''
    (report_dir / 'all_results.html').write_text(txt, encoding='utf-8')
    with (report_dir / 'all_results.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['model', 'case', 'chi2', *OBS, 'run_id', 'status', 'source_path'])
        writer.writeheader(); writer.writerows(rows)


def build(run_ids: list[str], output_dir: str, trajectory_dir: str | None, sensitivity_csv: str | None) -> dict[str, Any]:
    rows, provenance = summary_rows(run_ids)
    if not rows:
        raise ValueError('No accepted frozen states were found in the selected run IDs')
    trajectories, traj_state = load_trajectories(trajectory_dir, provenance)
    sensitivity, sens_state = load_sensitivity(sensitivity_csv, provenance)
    report_dir = REPORT_BASE / output_dir
    report_dir.mkdir(parents=True, exist_ok=True)
    build_index(report_dir, rows, traj_state, sens_state)
    build_rge(report_dir, trajectories)
    build_sensitivity(report_dir, sensitivity)
    build_table(report_dir, rows)
    manifest = {
        'accepted_row_count': len(rows),
        'trajectory_model_count': len({(x['model'], x['case']) for x in trajectories}),
        'sensitivity_model_count': len({(x['model'], x['case']) for x in sensitivity}),
        'trajectory_state': traj_state,
        'sensitivity_state': sens_state,
        'run_ids': run_ids,
    }
    (report_dir / 'provenance.json').write_text(json.dumps({'files': provenance, 'manifest': manifest}, indent=2, sort_keys=True), encoding='utf-8')
    return manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-ids', default='frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1')
    ap.add_argument('--trajectory-dir', default='')
    ap.add_argument('--sensitivity-csv', default='')
    ap.add_argument('--output-dir', default='frozen_massnorm_reports_v1')
    args = ap.parse_args(argv)
    run_ids = [x.strip() for x in args.run_ids.split(',') if x.strip()]
    manifest = build(run_ids, args.output_dir, args.trajectory_dir or None, args.sensitivity_csv or None)
    print(f"Wrote {REPORT_BASE / args.output_dir} | {manifest['accepted_row_count']} accepted rows | {manifest['trajectory_model_count']} trajectory models | {manifest['sensitivity_model_count']} sensitivity models")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
