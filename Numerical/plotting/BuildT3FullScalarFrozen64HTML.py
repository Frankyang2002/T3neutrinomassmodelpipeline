"""Adapter to the ORIGINAL InteractiveModelComparison.py RGE HTML canvas.

Render the saved frozen-64 SM+Weinberg trajectories and neutrino observables
in the existing interactive RGE report. Does not render scalar RGEs, substitute
old scalar report output, or manufacture uncomputed intermediate LLSS Wilson flow.
"""
from __future__ import annotations

import argparse
import ast
import csv
import html
import json
import math
from pathlib import Path

from scripts import T3Frozen64Sources as policy

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'output/studies/general_scalar_extended/full_scalar_frozen_64'
OUTPUT = ROOT / 'Reports/output/display_general_scalar_frozen_64'
WEINBERG_FORMAT = 't3-frozen-64-saved-sm-weinberg-trajectory-v1'

CSS = '''
:root{font-family:system-ui,-apple-system,Segoe UI,Arial,sans-serif;color:#203048;background:#f4f7fb}
*{box-sizing:border-box}body{margin:0}.wrap{max-width:1510px;margin:auto;padding:20px 26px}
h1{font-size:23px;margin:0 0 7px}h2{font-size:17px;margin:0 0 9px}h3{font-size:14px}
a{color:#2366a5}nav{display:flex;gap:15px;flex-wrap:wrap;margin:13px 0 20px}
.muted{color:#53657b;font-size:13px}.card{background:white;border:1px solid #d7e1eb;border-radius:9px;padding:16px;margin:14px 0}
.grid{display:grid;grid-template-columns:1fr 290px;gap:18px}.pill{font-size:11px;background:#e9eff6;padding:3px 7px;border-radius:5px}
label{font-size:13px;font-weight:600}select{width:100%;font:13px system-ui;padding:8px;margin:5px 0 12px;border:1px solid #c6d4e4;border-radius:5px;background:#fff}
input{accent-color:#2f68ac}.controls{display:grid;grid-template-columns:repeat(2,minmax(170px,1fr));gap:10px}
.canvas{width:100%;height:520px;display:block}.legend{display:flex;gap:12px;flex-wrap:wrap;margin:12px 0}
.tag{display:inline-block;padding:3px 8px;border-radius:4px;background:#f0f3f7;font-size:11px}
.scroller{overflow:auto}table{border-collapse:collapse;width:100%;background:#fff;font-size:12px}
th,td{text-align:left;padding:8px 10px;border-bottom:1px solid #e4eaf1;white-space:nowrap}th{background:#eef3f8;position:sticky;top:0}
.warning{border-left:3px solid #c77d36;padding:9px 12px;background:#fff9f0}
.success{color:#1d735a}.bad{color:#a1443f;font-weight:600}
@media(max-width:950px){.grid{grid-template-columns:1fr}.wrap{padding:16px}.canvas{height:460px}}
'''


def esc(value):
    return html.escape(str(value), quote=True)


def make_page(title, body):
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{esc(title)}</title><style>{CSS}</style></head><body>'
            f'<div class="wrap"><h1>{esc(title)}</h1>{body}</div></body></html>')


def load_validated(cell, run_id):
    """Verify frozen campaign data and source digests before using any results."""
    b, m = cell['benchmark'], cell['model']
    files = {}
    if cell['status'] in ('SCALAR_COMPLETE', 'NEUTRINO_COMPLETE'):
        expected = (f'output/studies/general_scalar_extended/full_scalar_frozen_64/'
                    f'{run_id}/{b}/{policy.tag(m)}/continuous_trajectory.json')
        policy.require(cell.get('trajectory') == expected, f'{b}/{m}: wrong scalar source')
        p = policy.validate_source(ROOT, cell['trajectory'], cell['trajectory_sha256'])
        t = policy.read_json(p)
        policy.require(t.get('model') == m and t.get('boundary_mode') == 'FROZEN_UV'
                       and t.get('settings', {}).get('benchmark') == b
                       and t.get('uv_ode_evaluations') == 0
                       and t.get('status') == 'CandidateContinuousUVFermionThresholdEFT1Integrated',
                       f'{b}/{m}: incorrect frozen scalar source')
        policy.require(t.get('input_sha256') == cell.get('source_sha256'),
                       f'{b}/{m}: source digest mapping changed')
        for file, sha in t['input_sha256'].items():
            policy.validate_source(ROOT, file, sha)
        files['trajectory'] = t
    if cell['status'] == 'NEUTRINO_COMPLETE':
        expected = (f'output/studies/general_scalar_extended/full_scalar_frozen_64/'
                    f'{run_id}/{b}/neutrino_v1/{policy.tag(m)}_neutrino.json')
        policy.require(cell.get('neutrino') == expected, f'{b}/{m}: wrong neutrino source')
        p = policy.validate_source(ROOT, cell['neutrino'], cell['neutrino_sha256'])
        n = policy.read_json(p)
        policy.require(n.get('model') == m and
                       n.get('source_continuous_sha256') == cell['trajectory_sha256']
                       and n.get('full_general_scalar_physical_c5_certified') is False,
                       f'{b}/{m}: wrong neutrino provenance')
        files['neutrino'] = n
    return files


def load_weinberg(cell, base, output_id):
    """No re-running, interpolation or substitutes in HTML; missing => unavailable."""
    if cell['status'] != 'NEUTRINO_COMPLETE':
        return None
    p = base / output_id / cell['benchmark'] / (policy.tag(cell['model']) + '.json')
    if not p.is_file():
        return None
    x = policy.read_json(p)
    src = x.get('source', {})
    generator = ROOT / 'scripts/GenerateT3Frozen64WeinbergRunning.py'
    policy.require(x.get('status') == 'VerifiedSMWeinbergEndPointParity'
                   and x.get('model') == cell['model']
                   and x.get('benchmark') == cell['benchmark']
                   and x.get('eft1_operator_c12_status') == 'NOT_SAVED_OR_INTEGRATED_FOR_THIS_FROZEN_RUN'
                   and src.get('trajectory_sha256') == cell['trajectory_sha256']
                   and src.get('neutrino_sha256') == cell['neutrino_sha256']
                   and src.get('generator_sha256') == policy.digest(generator),
                   f'{cell["benchmark"]}/{cell["model"]}: incompatible Weinberg data')
    pts = x.get('points')
    policy.require(isinstance(pts, list) and len(pts) >= 3, 'short Weinberg trajectory')
    scales = [float(r['mu_gev']) for r in pts]
    policy.require(math.isclose(scales[0], 1e3) and math.isclose(scales[-1], 100.)
                   and all(a > b > 0 for a, b in zip(scales, scales[1:])),
                   'unexpected Weinberg scales')
    for r in pts:
        nums = [float(r['c5_norm'])] + [float(v) for v in r['c5_components'].values()]
        nums += [float(v) for v in r['masses_ev_descending']]
        nums += [float(r['observables'][f]) for f in policy.FIELDS]
        policy.require(len(r['c5_components']) == 6 and len(r['masses_ev_descending']) == 3
                       and all(math.isfinite(v) for v in nums), 'bad Weinberg point')
    return x


def endpoint_screen(neutrino, target_dm3=2.511e-3):
    """A conspicuous scale diagnostic, NOT a hypothesis test or fit chi-square."""
    if not neutrino:
        return 'NOT_CALCULATED'
    masses = neutrino.get('takagi_masses_ev_descending') or []
    gap = (neutrino.get('oscillation_observables_normal_ordering') or {}).get('delta_m3l_sq_ev2')
    if len(masses) != 3 or gap is None:
        return 'INCOMPLETE'
    if float(masses[0]) > 1. or abs(float(gap)) > 10 * target_dm3:
        return 'MASS_SCALE_OUTLIER'
    return 'CANDIDATE_ONLY'


def finite_text(v, digits=5):
    try:
        value = float(v)
        return f'{value:.{digits}g}' if math.isfinite(value) else '—'
    except (ValueError, TypeError):
        return '—'


def table_rows(cells):
    fields = ['Model','Status','C12 EFT1','C5 SM RGE','Largest mν [eV]',
              'Δm²21 [eV²]','Δm²3l [eV²]','sin²θ12','sin²θ13','sin²θ23','Scale check','Matching']
    out = []
    for c, f, w in cells:
        n = f.get('neutrino')
        o = n.get('oscillation_observables_normal_ordering') if n else {}
        o = o or {}
        screen = endpoint_screen(n)
        out.append([c['model'], c['status'], 'Not computed',
                    'Verified' if w is not None else 'Not available',
                    finite_text((n or {}).get('takagi_masses_ev_descending', [None])[0]),
                    *(finite_text(o.get(k)) for k in policy.FIELDS),
                    screen, (n or {}).get('status') or '—'])
    th = '<tr>' + ''.join(f'<th>{esc(k)}</th>' for k in fields) + '</tr>'
    trs = ''.join('<tr>' + ''.join(f'<td class="{"bad" if v=="MASS_SCALE_OUTLIER" else ""}">{esc(v)}</td>'
                                  for v in row) + '</tr>' for row in out)
    return th, trs, fields, out


def original_interactive_rge_template() -> str:
    """Extract the ACTUAL original HTML canvas from InteractiveModelComparison.py.

    Do not fork the HTML/JS renderer. Fail if its source structure changes.
    """
    src = ROOT / 'Numerical/plotting/InteractiveModelComparison.py'
    tree = ast.parse(src.read_text(encoding='utf-8-sig'))
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == 'write_dashboard':
            for statement in node.body:
                if (isinstance(statement, ast.Assign)
                    and any(isinstance(t, ast.Name) and t.id == 'template'
                            for t in statement.targets)
                    and isinstance(statement.value, ast.Constant)
                    and isinstance(statement.value.value, str)):
                    template = statement.value.value
                    policy.require(all(t in template for t in
                        ('__TITLE__','__OPTIONS__','__DATA__','__META__',
                         'id="quantity"','id="plot"','id="zoomFinal"',
                         'id="yScaleMode"','id="resetZoom"')),
                         'original interactive RGE template contract changed')
                    return template
    raise RuntimeError('Missing original write_dashboard() template; no duplicate HTML fallback')


def curve_from_points(points, field):
    x, y = [], []
    for p in points:
        mu = float(p['mu_gev'])
        if field == 'c5_norm':
            v = float(p['c5_norm'])
        elif field.startswith('c5_'):
            v = float(p['c5_components'][field[3:]])
        elif field in ('m1','m2','m3'):
            # Original InteractiveModelComparison identifiers are physical NO
            # m1 <= m2 <= m3; saved Takagi masses are descending.
            v = float(p['masses_ev_descending'][3-int(field[-1])])
        else:
            source = {'dm21':'delta_m21_sq_ev2',
                      'dm3l':'delta_m3l_sq_ev2',
                      's12':'sin2_theta12','s13':'sin2_theta13',
                      's23':'sin2_theta23'}.get(field,field)
            v = float(p['observables'][source])
        policy.require(math.isfinite(mu) and mu > 0 and math.isfinite(v),
                       f'nonfinite plotted observable: {field}')
        x.append(mu)
        y.append(v)
    policy.require(len(x) >= 3 and all(a > b for a, b in zip(x,x[1:])),
                   'RGE data must be saved on a decreasing energy grid')
    return {'x':x,'y':y}


# Original-model selector keys retained where defined in InteractiveModelComparison.
RGE_OPTIONS = (
    ('Weinberg operator C5 (final SM EFT)', (
        ('c5_norm','||C₅||F','GeV⁻¹'),
        ('c5_11','|C₅₁₁|','GeV⁻¹'),('c5_12','|C₅₁₂|','GeV⁻¹'),
        ('c5_13','|C₅₁₃|','GeV⁻¹'),('c5_22','|C₅₂₂|','GeV⁻¹'),
        ('c5_23','|C₅₂₃|','GeV⁻¹'),('c5_33','|C₅₃₃|','GeV⁻¹'),
    )),
    ('Neutrino masses', (
        ('m1','m₁','eV'),('m2','m₂','eV'),('m3','m₃','eV'),
    )),
    ('Neutrino oscillation observables', (
        ('dm21','Δm²₂₁','eV²'),('dm3l','Δm²₃ℓ','eV²'),
        ('s12','sin²θ₁₂',''),('s13','sin²θ₁₃',''),
        ('s23','sin²θ₂₃',''),
    )),
)


def rge_payload(records):
    """Convert verified operator trajectories to original interactive RGE schema."""
    data = {'models':list(policy.MODELS),
            'scales':{'mu_uv':1e7,'mu_fermion_threshold':1e5,
                      'mu_scalar_threshold':1e3,'mu_low':100.},
            'quantities':{},'overlays':{},'matched':{},'components':{},
            'continuation_methods':{},'continuation_issues':{}}
    meta = {}
    options = []
    for group, items in RGE_OPTIONS:
        choices = []
        for key, label, unit in items:
            meta[key] = {'label':label,'unit':unit, 'block':'final_running',
                         # Do not trigger legacy hard+direct continuation display:
                         # only a genuine 1 TeV to 100 GeV C5 is available.
                         'weinberg':False}
            curves = {}
            for cell, _files, trajectory in records:
                if trajectory is None:
                    continue
                curves[cell['model']] = curve_from_points(trajectory['points'],key)
            if curves:
                data['quantities'][key] = curves
            choices.append(f'<option value="{esc(key)}">{esc(label)}</option>')
        options.append(f'<optgroup label="{esc(group)}">'+''.join(choices)+'</optgroup>')
    return data, meta, ''.join(options)


def render_original_interactive(template, benchmark, records):
    """Use the source HTML/JS canvas; adapt DATA and textual physics labels only."""
    policy.require(all(marker in template for marker in
                ('__TITLE__','__OPTIONS__','__DATA__','__META__')),
                'original interactive canvas placeholders missing')
    data, meta, options = rge_payload(records)
    title = 'T3 interactive RGE — C5 and neutrino observables — '+benchmark
    page = (template.replace('__TITLE__',esc(title))
            .replace('__OPTIONS__',options)
            .replace('__DATA__',json.dumps(data,allow_nan=False,separators=(',',':')).replace('<','\\u003c'))
            .replace('__META__',json.dumps(meta,allow_nan=False,separators=(',',':')).replace('<','\\u003c')))
    # The ORIGINAL report defaults to a 1e7 -> 100 GeV view, appropriate for
    # its multi-EFT dataset. Here every available trace begins at 1 TeV.
    original = 'return [s.mu_uv||1e7,s.mu_low||100];'
    policy.require(original in page, 'original RGE zoom logic changed')
    page = page.replace(original, 'return [s.mu_scalar_threshold||1e3,s.mu_low||100];')
    old_header = ('Numerical integration direction: UV → IR. Saved benchmark '
                  'trajectories; no matching or RGE recalculation.')
    policy.require(old_header in page, 'original RGE header text changed')
    page = page.replace(old_header,
        'Original T3 interactive RGE renderer. Verified saved numerical SM+Weinberg '
        'running at μ = 1 TeV → 100 GeV; no physics recomputation in HTML.')
    # Original axis controls, checkbox model overlays, zoom buttons, scroll/
    # double-click zoom, and y-axis scaling all remain owned by the old canvas.
    count = sum(w is not None for _,_,w in records)
    outliers = [c['model'] for c,f,w in records
                if endpoint_screen(f.get('neutrino')) == 'MASS_SCALE_OUTLIER']
    banner = ('<div style="margin:12px 20px;padding:12px;border:1px solid #d6dee8;'
              'font-size:13px;background:#f8fafc">'
              '<strong>Numerical scope:</strong> final Weinberg coefficient C₅(μ) '
              'and derived neutrino masses, gaps and mixing angles, 1 TeV → 100 GeV. '
              'The intermediate EFT C₁₂ (LLSS) Wilson trajectory was <strong>not computed</strong> '
              'for these frozen-64 cases; no η-kernel substitute is displayed. '
              'UV is frozen. Model-specific general-scalar matching remains conditional; '
              'six self-bilinear cases omit their mixing contribution. '
              f'Verified SM+Weinberg trajectories: {count}/16. '
              f'Mass-scale outliers: {len(outliers)}/16 '
              f'({esc(", ".join(outliers)) if outliers else "none"}). '
              '<a href="all_results.html">All numerical results</a> · '
              '<a href="sensitivity_ovals.html">Sensitivity ovals</a> · '
              '<a href="../index.html">All benchmarks</a>.</div>')
    anchor = '<div class="layout"><aside>'
    policy.require(anchor in page, 'original interactive layout marker changed')
    page = page.replace(anchor,banner+anchor,1)
    # Avoid accidentally leaving placeholders unresolved.
    policy.require(all(k not in page for k in
                  ('__TITLE__','__OPTIONS__','__DATA__','__META__')),
                  'unresolved original interactive HTML placeholders')
    return page


def operator_html(benchmark, records):
    return render_original_interactive(original_interactive_rge_template(),benchmark,records)


# Kept from the previous production report: the original empirical gamma ovals
# are sourced from exact per-benchmark files, not from historical reduced fits.
def scans_for(cell, files, root, scan_id):
    if 'neutrino' not in files:
        return None
    path=root/scan_id/cell['benchmark']/(policy.tag(cell['model'])+'.csv')
    meta=path.with_suffix('.json')
    if not path.is_file() and not meta.is_file():
        return None
    policy.require(path.is_file() and meta.is_file(),'partial gamma scan: '+str(path))
    doc=policy.read_json(meta)
    policy.require(doc['csv_sha256']==policy.digest(path),'gamma scan hash mismatch')
    source=doc.get('source',{})
    policy.require(source.get('model')==cell['model'] and
                   source.get('benchmark')==cell['benchmark'] and
                   source.get('trajectory_sha256')==cell['trajectory_sha256'] and
                   source.get('neutrino_sha256')==cell['neutrino_sha256'] and
                   source.get('generator_sha256') == policy.digest(
                       ROOT / 'scripts/GenerateT3FrozenFullScalarGamma.py'),
                   'gamma scan is from an incompatible numerical source')
    with path.open(newline='',encoding='utf-8-sig') as f:
        lines=list(csv.DictReader(f))
    required={(d,float(g)) for d in policy.DIRECTIONS for g in source['gamma']}
    policy.require(len(lines)==len(required),'incomplete gamma scan')
    indexed={}
    for r in lines:
        k=r['direction'],float(r['gamma'])
        policy.require(k in required and k not in indexed,'duplicate/unknown gamma point')
        policy.require(r['model']==cell['model'] and r['benchmark']==cell['benchmark'] and
                       r['status']=='NeutrinoEvaluated','invalid gamma result')
        p={'gamma':k[1]}
        for field in policy.FIELDS:
            v=float(r[field]);policy.require(math.isfinite(v),'nonfinite gamma response')
            p[field]=v
        indexed[k]=p
    policy.require(set(indexed)==required,'missing gamma point')
    bydir={d:sorted((indexed[d,g] for dd,g in required if dd==d),key=lambda x:x['gamma'])
           for d in policy.DIRECTIONS}
    centre=[next(p for p in bydir[d] if abs(p['gamma']-1)<1e-12) for d in policy.DIRECTIONS]
    for a in centre[1:]:
        for field in policy.FIELDS:
            policy.require(math.isclose(a[field],centre[0][field],rel_tol=1e-7,abs_tol=1e-10),
                           'gamma=1 centres disagree')
    saved=files['neutrino']['oscillation_observables_normal_ordering']
    for field in policy.FIELDS:
        policy.require(math.isclose(centre[0][field],float(saved[field]),rel_tol=1e-6,abs_tol=1e-10),
                       'gamma centre does not reproduce central neutrino')
    return {'model':cell['model'],'label':cell['model'],'status':files['neutrino']['status'],
            'fit_stage':'historical_seed_on_frozen_general_scalar',
            'centre':centre[0], 'directions':bydir}


def sensitivity_html(benchmark,rows,base,scan_id):
    plotted=[]
    missing=[]
    for cell,files in rows:
        item=scans_for(cell,files,base,scan_id)
        if item:
            plotted.append(item)
        else:
            missing.append({'model':cell['model'],'status':'No verified complete full-scalar gamma scan'})
    from Numerical.plotting import BuildFittedModelEnvelopes as original
    ref=policy.read_json(ROOT/'docs/FROZEN_LAMBDA_TARGET_REFERENCE.json')
    target=ref['target']
    for m in plotted:
        # Pure display chi-square: historical diagonal Gaussian NuFIT covariance.
        c=m['centre']; keys=target['observable_names']
        m['chi2']=sum((c[k]-target['central_values'][i])**2/target['covariance'][i][i]
                      for i,k in enumerate(keys))
    payload={'models':plotted,'skipped':missing,
             'target':target,'labels':{k:original.LABELS[k] for k in policy.FIELDS},
             'observables':list(policy.FIELDS),
             'source':'verified frozen full-scalar gamma CSVs: '+scan_id+'/'+benchmark,
             'future_sigma':ref['future_sigma'],'future_source':ref['future_source'],
             'envelope_definition':'2x empirical sample covariance of actual gamma responses; descriptive, not confidence regions'}
    text=original.HTML.replace('__DATA__',json.dumps(payload,allow_nan=False).replace('<','\\u003c'))
    policy.require(text!=original.HTML,'missing sensitivity renderer placeholder')
    text=text.replace('Individually fitted T3 model sensitivity','Frozen general-scalar sensitivity — '+benchmark)
    text=text.replace('Each model uses its own fitted λT3, Yukawa matrices and other fitted parameters. Data from saved scans only — no RGE recomputation.',
                'Genuine T1/T2/T3 perturbations at 100 TeV. Original historical Yukawa seeds are NOT reoptimised for the general scalar potential.')
    text=text.replace('fitted models plotted','frozen-general-scalar candidate models plotted')
    text=text.replace('about each fitted centre','about each saved gamma=1 centre')
    explanation=('<div style="margin:12px;padding:12px;background:#fff7e6;border:1px solid #ddc999">'
                 '<b>Physics status:</b> General-scalar candidate RGE/matching only. '
                 'Six A/D/E self-bilinear C5 calculations are restricted unmixed controls; '
                 'the B/C one-scale mixing result is conditional. No complete sequential EFT matching. '
                 'The blue NuFIT and orange JUNO ovals are historical Gaussian approximations, '
                 'not official likelihood contours. Missing gamma data are not replaced with old scans. '
                 f'Full-gamma models: {len(plotted)}/16.</div>')
    return text.replace('<div class="layout">',explanation+'<div class="layout">',1),len(plotted)



def write_report(run_id: str, scan_id: str, output_id: str, destination: Path) -> dict:
    base = DATA / run_id
    manifest_file = base / 'manifest.json'
    manifest = policy.validate_manifest(policy.read_json(manifest_file))
    # Fully validate all central outputs once, and no wildcard discovery.
    data = {}
    for cell in manifest['cells']:
        f = load_validated(cell, run_id)
        w = load_weinberg(cell, base, output_id)
        data[(cell['benchmark'], cell['model'])] = (cell, f, w)
    destination.mkdir(parents=True, exist_ok=True)
    coverage = []
    nav = []
    for bench in policy.BENCHMARKS:
        records = [data[(bench, model)] for model in policy.MODELS]
        folder = destination / bench
        folder.mkdir(parents=True, exist_ok=True)
        (folder / 'interactive_rge.html').write_text(
            operator_html(bench, records), encoding='utf-8')
        # Keep original gamma scan / oval presentation as a separate view.
        sensitive, count = sensitivity_html(
            bench, [(cell, files) for cell, files, _ in records], base, scan_id)
        (folder / 'sensitivity_ovals.html').write_text(sensitive, encoding='utf-8')
        th, trs, fields, rows = table_rows(records)
        (folder / 'all_results.html').write_text(make_page(
            bench + ': 100 GeV candidate neutrino results',
            f'<nav><a href="../index.html">All benchmarks</a>'
            f'<a href="interactive_rge.html">C5 running and observables</a></nav>'
            f'<div class="card scroller"><table><thead>{th}</thead><tbody>{trs}</tbody></table></div>'),
            encoding='utf-8')
        with (folder / 'all_results.csv').open('w', newline='', encoding='utf-8') as stream:
            wr = csv.writer(stream)
            wr.writerow(fields)
            wr.writerows(rows)
        nav.append(f'<div class="card"><h2>{esc(bench)}</h2><nav>'
                   f'<a href="{bench}/interactive_rge.html">C5 + observables RGE</a>'
                   f'<a href="{bench}/all_results.html">Neutrino results</a>'
                   f'<a href="{bench}/sensitivity_ovals.html">Sensitivity ovals</a>'
                   f'</nav><p class="muted">{sum(w is not None for _, _, w in records)}/16 '
                   f'verified Weinberg trajectories; {count}/16 genuine gamma scans.</p></div>')
        for c, f, w in records:
            n = f.get('neutrino')
            coverage.append({'benchmark': bench, 'model': c['model'],
                             'central_neutrino': bool(n),
                             'weinberg_running_verified': w is not None,
                             'c12_eft1_running': False,
                             'sensitivity_scan': bool(scans_for(c, f, base, scan_id)),
                             'mass_scale_check': endpoint_screen(n),
                             'matching_status': (n or {}).get('status', '')})
        print(f'{bench}: verified final-Weinberg/observable running '
              f'{sum(w is not None for _,_,w in records)}/16; '
              f'sensitivity scans {count}/16')
    with (destination / 'coverage.csv').open('w', newline='', encoding='utf-8') as stream:
        wr = csv.DictWriter(stream, fieldnames=list(coverage[0]))
        wr.writeheader()
        wr.writerows(coverage)
    proof = {'run_id': run_id, 'scan_id': scan_id, 'weinberg_output_id': output_id,
             'source_manifest': policy.relative(ROOT, manifest_file),
             'manifest_sha256': policy.digest(manifest_file),
             'weinberg_generator_sha256': policy.digest(
                 ROOT / 'scripts/GenerateT3Frozen64WeinbergRunning.py'),
             'full_general_scalar_sequential_matching_certified': False,
             'c12_trajectory_saved': False,
             'scope': 'Numerical SM+Weinberg running and neutrino observables from saved C5(1TeV)',
             'coverage': {'slots': 64,
                          'central_neutrino': sum(c['central_neutrino'] for c in coverage),
                          'verified_weinberg_rge': sum(c['weinberg_running_verified'] for c in coverage),
                          'c12_rge': 0,
                          'sensitivity': sum(c['sensitivity_scan'] for c in coverage),
                          'mass_scale_outliers': sum(c['mass_scale_check'] == 'MASS_SCALE_OUTLIER'
                                                     for c in coverage)}}
    (destination / 'provenance.json').write_text(
        json.dumps(proof, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    intro = ('<p class="muted">16 T3 models × 4 frozen-UV Yukawa/scalar choices. '
             'No scalar coupling or scalar mass RGE plots are produced by this generator.</p>'
             '<div class="card warning"><strong>Interpretation:</strong> 64 numerical results '
             'do not imply 64 physical neutrino fits. The low-coupling general-scalar '
             'candidates can produce masses far above 0.05 eV. The C12 intermediate EFT '
             'running was not saved, and complete general-scalar sequential matching is '
             'uncertified. Unavailable quantities are stated explicitly, never inferred.</div>')
    body = (intro + ''.join(nav) + '<p><a href="coverage.csv">Coverage audit</a> · '
            '<a href="provenance.json">Source provenance</a></p>')
    (destination / 'index.html').write_text(make_page(
        'T3 frozen general-scalar: Weinberg operators and neutrino observables', body),
        encoding='utf-8')
    return proof


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-id', default='frozen_general64_v1')
    parser.add_argument('--scan-id', default='gamma_scans_v1')
    parser.add_argument('--weinberg-output-id', default='weinberg_running_v1')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    for name in (args.run_id, args.scan_id, args.weinberg_output_id):
        policy.require(bool(policy.RUN_ID.fullmatch(name)), 'unsafe run/scan/output ID')
    dest = args.output or OUTPUT / args.run_id
    proof = write_report(args.run_id, args.scan_id, args.weinberg_output_id, dest)
    print('Report:', dest / 'index.html')
    print('Coverage:', proof['coverage'])


if __name__ == '__main__':
    main()
