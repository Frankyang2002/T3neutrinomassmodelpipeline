"""Rebuild A(0) shared-Yukawa comparison dashboard and cross-model trajectory envelopes.

Presentation only. No numerical physics recomputation. Preserves existing outputs.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np

from Numerical.plotting.InteractiveModelComparison import write_dashboard
from Numerical.plotting.RunningResultFigures import generate_running_result_figures

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'output/studies/small_lambda_1e-7/a0_shared_benchmark'
OUTPUT = ROOT / 'Reports/output/display_a0_shared_fixed'
STRICT = 5.86
NOMINAL95 = 11.070497693516351
FIELDS = {
    'delta_m21_sq_ev2': ('Solar splitting', 'eV²'),
    'delta_m3l_sq_ev2': ('Atmospheric splitting', 'eV²'),
    'sin2_theta12': ('sin² θ12', ''),
    'sin2_theta13': ('sin² θ13', ''),
    'sin2_theta23': ('sin² θ23', ''),
}


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def find_records(source):
    model_root = source / 'models'
    if not model_root.is_dir():
        raise FileNotFoundError(f'Missing model diagnostics root: {model_root}')
    rows = []
    for folder in sorted(model_root.iterdir()):
        if not folder.is_dir() or not folder.name.startswith('T3_'):
            continue
        matches = list(folder.rglob('running_diagnostics.json'))
        if len(matches) != 1:
            raise RuntimeError(f'{folder.name}: expected one diagnostics JSON; found {len(matches)}')
        doc = read(matches[0])
        if doc.get('status') != 'Success':
            raise RuntimeError(f'{folder.name}: diagnostics not Success')
        rows.append((folder.name, matches[0], doc))
    if len(rows) != 16:
        raise RuntimeError(f'Expected 16 distinct A(0) model diagnostics, found {len(rows)}')
    return rows


def chi2_rows(source):
    result = {}
    file = source / 'comparison.csv'
    if not file.is_file():
        raise FileNotFoundError(f'Missing comparison file {file}')
    with file.open('r', encoding='utf-8-sig', newline='') as handle:
        for row in csv.DictReader(handle):
            result[row['model']] = row
    return result


def get_chi2(row):
    for field in ('chi2', 'common_chi2'):
        if row.get(field) not in (None, ''):
            return float(row[field])
    raise KeyError('comparison.csv requires chi2 or common_chi2 field')


def create_dashboard(records, destination):
    # InteractiveModelComparison._model_name takes path.parent.parent.name.
    # Source layout is models/<name>/running_diagnostics.json, giving 'models'.
    # Stage diagnostics as <name>/data/running_diagnostics.json instead.
    with tempfile.TemporaryDirectory(prefix='t3_a0_dashboard_') as temp:
        staged = Path(temp)
        for name, path, _ in records:
            target = staged / name / 'data' / 'running_diagnostics.json'
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
        result = write_dashboard(staged, destination,
            'T3 A(0) shared UV Yukawas — 16 representation-specific trajectories')
        content = result.read_text(encoding='utf-8')
        if not all(name in content for name, _, _ in records):
            raise RuntimeError('Dashboard did not contain all 16 model identifiers')
    return destination


def build_envelope(records, selected, output):
    # Spread of five successful distinct models, NOT a statistical interval.
    # Original distinct trajectory sampling grids usually coincide; interpolate
    # in log scale when needed and refuse extrapolation.
    data = {'description': 'Pointwise range and centre across five strict-pass A(0) common-Yukawa models; not a parameter uncertainty or confidence region.',
            'model_names': selected, 'quantities': {}}
    by_name = {name: doc for name, _, doc in records}
    for field, (label, units) in FIELDS.items():
        curves = []
        for name in selected:
            block = by_name[name]['final_running']
            x = np.asarray(block['mu_gev'], dtype=float)
            y = np.asarray(block[field], dtype=float)
            if x.ndim != 1 or y.shape != x.shape or len(x) < 2 or np.any(x <= 0) or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
                raise ValueError(f'Invalid {field} for {name}')
            order = np.argsort(x)
            curves.append((x[order], y[order]))
        lo = max(x[0] for x, _ in curves)
        hi = min(x[-1] for x, _ in curves)
        if not lo < hi:
            raise ValueError('No common scale interval')
        xgrid = np.geomspace(lo, hi, 120)
        values = np.vstack([np.interp(np.log(xgrid), np.log(x), y) for x, y in curves])
        data['quantities'][field] = {
            'label': label, 'unit': units, 'mu': xgrid.tolist(),
            'min': values.min(axis=0).tolist(), 'max': values.max(axis=0).tolist(),
            'mean': values.mean(axis=0).tolist(),
            'curves': {name: values[i].tolist() for i, name in enumerate(selected)}
        }
    output.write_text(_envelope_html(data), encoding='utf-8')
    (output.parent / 'envelope_data.json').write_text(json.dumps(data, indent=2), encoding='utf-8')


def _envelope_html(data):
    blob = json.dumps(data, ensure_ascii=False).replace('</', '<\\/')
    return '''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>A(0) common Yukawa: model spread</title>
<style>body{font:15px system-ui,Arial,sans-serif;max-width:1300px;margin:22px auto;padding:0 18px;color:#182638}h1{font-size:24px}p{line-height:1.5}fieldset{border:1px solid #d4dce6;border-radius:8px;margin:12px 0;padding:12px}label{margin:6px 16px 6px 0;display:inline-block}canvas{width:100%;height:540px;border:1px solid #d4dce6;border-radius:8px}select{padding:8px}small{color:#56677a}</style></head><body>
<h1>A(0) shared Yukawa — spread across compatible T3 models</h1><p>The shaded envelope is the pointwise <b>minimum–maximum across the five strict-pass models</b> at each renormalisation scale. It is <b>not</b> a 1σ/3σ uncertainty band, a NuFIT confidence region, or a scalar-mass sensitivity scan.</p>
<label>Observable <select id="quantity"></select></label><fieldset><legend>Models (toggle individual curves)</legend><div id="checks"></div></fieldset><canvas id="plot" width="1250" height="540"></canvas><p id="foot"></p><p><a href="index.html">Back to benchmark summary</a> · <a href="interactive_comparison.html">All-model interactive running dashboard</a></p>
<script>const D=__DATA__;const q=document.getElementById('quantity'),check=document.getElementById('checks'),canvas=document.getElementById('plot'),ctx=canvas.getContext('2d');let palette=['#1d5ca8','#a1342b','#228352','#8b4ca5','#9d7300'];for(const [k,v]of Object.entries(D.quantities)){let o=document.createElement('option');o.value=k;o.textContent=v.label;q.appendChild(o)}D.model_names.forEach((name,i)=>{let label=document.createElement('label'),c=document.createElement('input');c.type='checkbox';c.checked=true;c.dataset.model=name;c.onchange=draw;label.appendChild(c);label.appendChild(document.createTextNode(name));label.style.color=palette[i%palette.length];check.appendChild(label)});q.onchange=draw;function draw(){const v=D.quantities[q.value],W=canvas.width,H=canvas.height,L=95,R=35,T=28,B=80,xs=v.mu.map(Math.log10),ys=[...v.min,...v.max],ylo=Math.min(...ys),yhi=Math.max(...ys),pad=(yhi-ylo)*.07||Math.max(Math.abs(yhi),1e-9)*.07;ylo-=pad;yhi+=pad;const X=t=>L+(t-xs[0])/(xs[xs.length-1]-xs[0])*(W-L-R),Y=t=>H-B-(t-ylo)/(yhi-ylo)*(H-T-B);ctx.clearRect(0,0,W,H);ctx.font='14px system-ui';ctx.fillStyle='#182638';ctx.strokeStyle='#dce3eb';for(let i=0;i<=5;i++){let yy=ylo+(yhi-ylo)*i/5;ctx.beginPath();ctx.moveTo(L,Y(yy));ctx.lineTo(W-R,Y(yy));ctx.stroke();ctx.fillText(yy.toExponential(3),4,Y(yy)+4)}for(let i=0;i<=5;i++){let xx=xs[0]+(xs[xs.length-1]-xs[0])*i/5;ctx.beginPath();ctx.moveTo(X(xx),T);ctx.lineTo(X(xx),H-B);ctx.stroke();ctx.fillText('10^'+xx.toFixed(2),X(xx)-27,H-B+25)}ctx.fillText('Renormalisation scale μ [GeV]',W/2-115,H-13);ctx.fillText(v.label+' ['+v.unit+']',L+10,17);ctx.fillStyle='rgba(69,129,194,.18)';ctx.beginPath();v.min.forEach((y,i)=>{i?ctx.lineTo(X(xs[i]),Y(y)):ctx.moveTo(X(xs[i]),Y(y))});for(let i=v.max.length-1;i>=0;i--)ctx.lineTo(X(xs[i]),Y(v.max[i]));ctx.closePath();ctx.fill();const drawLine=(arr,color,width)=>{ctx.strokeStyle=color;ctx.lineWidth=width;ctx.beginPath();arr.forEach((val,i)=>i?ctx.lineTo(X(xs[i]),Y(val)):ctx.moveTo(X(xs[i]),Y(val)));ctx.stroke()};drawLine(v.mean,'#183b6e',3);D.model_names.forEach((name,i)=>{if(check.querySelector('input[data-model="'+name+'"]')?.checked)drawLine(v.curves[name],palette[i%palette.length],1.5)});document.getElementById('foot').textContent='Five-model range (shaded); mean (dark line); individual model curves (toggle above). Fixed λT3=10⁻⁷ for these A/B/D models.'}draw();</script></body></html>'''.replace('__DATA__', blob)


def build(source, output, refresh_figures=True):
    records = find_records(source)
    comparisons = chi2_rows(source)
    if set(comparisons) != {name for name, _, _ in records}:
        raise ValueError('comparison.csv models differ from successful diagnostics')
    selected = [name for name, _, _ in records if get_chi2(comparisons[name]) < STRICT]
    if len(selected) != 5:
        raise ValueError(f'Expected 5 strict-pass models, found {len(selected)}')
    output.mkdir(parents=True, exist_ok=True)
    create_dashboard(records, output / 'interactive_comparison.html')
    build_envelope(records, selected, output / 'model_spread_envelope.html')
    if refresh_figures:
        for i, (name, path, _) in enumerate(records, 1):
            generate_running_result_figures(path, output/'figures'/name)
            print(f'figures [{i}/16] {name}', flush=True)
    with (output/'model_summary.csv').open('w', newline='', encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['model','chi2','strict_chi2_lt_5_86','nominal95_chi2_le_11_0705'])
        for name, _, _ in records:
            chi2=get_chi2(comparisons[name]);w.writerow([name,chi2,chi2<STRICT,chi2<=NOMINAL95])
    body=''.join('<tr><td>'+html.escape(name)+'</td><td>'+format(get_chi2(comparisons[name]),'.7g')+'</td><td>'+('Yes' if name in selected else 'No')+'</td></tr>' for name,_,_ in records)
    (output/'index.html').write_text('''<!doctype html><meta charset="utf-8"><title>T3 A(0) shared benchmark display</title><style>body{font:15px system-ui,Arial;max-width:1160px;margin:30px auto;padding:0 18px;color:#17243b}table{border-collapse:collapse;width:100%}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left}a{color:#185ca8}</style><h1>T3 A(0) shared Yukawa benchmark</h1><p>16/16 trajectories calculated; 5 strict χ²&lt;5.86, 6 nominal χ²≤11.0705. Common A(0) Y₁,Y₂ applied to every model; λT3=10⁻⁷ for A–D, 5×10⁻⁸ for E.</p><p><a href="interactive_comparison.html">All 16 models: interactive running</a> | <a href="model_spread_envelope.html">Five compatible models: trajectory envelope</a> | <a href="model_summary.csv">Summary CSV</a></p><p>The envelope represents model-to-model spread, not parameter uncertainties or a statistical confidence region. Screening χ² thresholds are nominal.</p><table><thead><tr><th>Model</th><th>χ²</th><th>Strict pass</th></tr></thead><tbody>'''+body+'</tbody></table>', encoding='utf-8')
    print('Rebuilt A(0) dashboard with 16 model selectors and five-model envelope:', output)
    print('Open', output/'index.html')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,default=SOURCE)
    ap.add_argument('--output',type=Path,default=OUTPUT)
    ap.add_argument('--skip-figures',action='store_true')
    a=ap.parse_args()
    build(a.source,a.output,not a.skip_figures)

if __name__=='__main__':main()
