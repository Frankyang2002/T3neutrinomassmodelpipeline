"""Build continuous_v1 general-scalar reports using the *existing Lambda HTML designs*.

Presentation only: reads saved JSON/CSV. Does not solve RGEs, match, fit,
perturb parameters or use historical reduced-potential gamma scans.

From the repository root:
    python -m Numerical.plotting.BuildGeneralScalarLambdaReports

Uses the literal HTML/canvas templates in InteractiveModelComparison.py and
BuildFittedModelEnvelopes.py, without modifying either original generator.
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import html
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
BASE = Path('output/studies/general_scalar_extended/continuous_full_scalar')
OUTPUT = Path('Reports/output/display_general_scalar')
MODELS = ('A_+0', 'A_-2', 'A_-4', 'B_+1', 'B_-1', 'B_-3',
          'C_+1', 'C_-1', 'C_-3', 'D_+0', 'D_+2', 'D_-2',
          'E_+0', 'E_+2', 'E_-2', 'E_-4')
FORMAT = 't3-continuous-general-scalar-candidate-v1'
STATUS = 'CandidateContinuousUVFermionThresholdEFT1Integrated'
SCHEMES = (
    ('baseline', 'neutrino_v1/summary.csv', 'Baseline benchmark',
     'Original general-scalar C5 benchmark; six self-bilinear terms omitted.'),
    ('self_massdep', 'self_anomalous_massdep_v1/combined_16_summary.csv',
     'Mass-dependent self-bilinear', 'One-scale conditional approximation; six corrected models only.'),
    ('direct', 'direct_wilson_sequential_v1/summary.csv', 'Direct Wilson transport',
     'Hierarchical hard plus time-resolved direct Weinberg source.'),
    ('extra', 'bc_extra_wilson_sequential_v1/summary.csv', 'B/C extra-quartic source',
     'Selected extra-quartic source; remaining terms not included.'),
    ('coupled', 'bc_coupled_wilson_v1/summary.csv', 'B/C coupled L/eC Wilson',
     'Implemented coupled candidate, not complete finite-order matching.'),
)
OBS = ('mass_hi_ev', 'mass_mid_ev', 'mass_lo_ev', 'delta_m21_sq_ev2',
       'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
LABELS = {'mass_hi_ev':'Largest mass [eV]', 'mass_mid_ev':'Middle mass [eV]',
          'mass_lo_ev':'Smallest mass [eV]', 'delta_m21_sq_ev2':'Δm²₂₁ [eV²]',
          'delta_m3l_sq_ev2':'Δm²₃ℓ [eV²]', 'sin2_theta12':'sin²θ₁₂',
          'sin2_theta13':'sin²θ₁₃', 'sin2_theta23':'sin²θ₂₃',
          'gap_small_mid_sq_ev2':'(middle² − smallest²) [eV²]; derived',
          'gap_small_large_sq_ev2':'(largest² − smallest²) [eV²]; derived'}
EXTRA = ('anomalous_over_normal_c5', 'relative_C5_change_vs_frozen_lambda',
         'relative_extra_C5_vs_direct_at_scalar_threshold',
         'relative_extra_C5_vs_frozen_extra',
         'relative_C5_change_vs_previous_selected',
         'direct_source_parity', 'extra_source_parity',
         'frozen_hierarchical_kernel_parity_error',
         'projection_relative_residual', 'weinberg_projection_residual',
         'baseline_largest_mnu_ev', 'old_largest_mnu_ev')
SCHEME_FLAGS = ('matching_includes_nonzero_bilinear', 'self_bilinear_c5_omitted',
                'matching_category', 'direct_wilson_source_evolved')
COLUMNS = ('model','source','source_title','status',*OBS,
           'gap_small_mid_sq_ev2','gap_small_large_sq_ev2',*EXTRA,*SCHEME_FLAGS,
           'quadratic_kind','scalar_endpoint_admissible','scalar_min_m2_gev2',
           'selection_reason','source_note','source_relative_path')


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def finite(value: Any) -> float | None:
    try:
        if value is None or value == '':
            return None
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError, OverflowError):
        return None


def complex_part(value: Any, component: str) -> float | None:
    if not isinstance(value, dict):
        return None
    real, imag = finite(value.get('real')), finite(value.get('imag'))
    if real is None or imag is None:
        return None
    return {'real':real, 'imag':imag, 'abs':math.hypot(real,imag)}[component]


def read_json(path: Path) -> dict:
    require(path.is_file(), f'Missing required JSON: {path}')
    result = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(result, dict), f'Expected JSON object in {path}')
    return result


def read_csv(path: Path) -> list[dict[str,str]]:
    with path.open(newline='',encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_script(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(',',':'),
                      allow_nan=False).replace('<','\\u003c')


def slug(model: str) -> str:
    return model.replace('+','p').replace('-','m')


def read_trajectories(run: Path, hashes: dict[str,str]) -> tuple[dict,dict]:
    models: dict[str,dict] = {}
    settings = None
    for name in MODELS:
        path = run / slug(name) / 'continuous_trajectory.json'
        data = read_json(path)
        require(data.get('model') == name, f'Model mismatch: {path}')
        require(data.get('format') == FORMAT and data.get('status') == STATUS,
                f'Not a validated continuous_v1 scalar trajectory: {path}')
        segs = {}
        for stage, key in (('UV','uv_trajectory'),('EFT1','eft1_trajectory')):
            rows = data.get(key)
            require(isinstance(rows,list) and len(rows)>=2,
                    f'Missing saved {stage} running for {name}')
            scales = [finite(v.get('mu_gev')) if isinstance(v,dict) else None for v in rows]
            require(all(v is not None and v>0 for v in scales),
                    f'Invalid {stage} scales: {name}')
            require(all(a>b for a,b in zip(scales,scales[1:])),
                    f'Non-descending {stage} scales: {name}')
            for row in rows:
                require(finite(row.get('mS1Sq_gev2')) is not None and
                        finite(row.get('mS2Sq_gev2')) is not None,
                        f'Missing scalar mass-squared running: {name}')
                require(isinstance(row.get('extra_quartics',{}),dict),
                        f'Invalid extra quartics: {name}')
            segs[stage] = rows
        model_settings = data.get('settings') or {}
        if settings is None:
            settings = {k:finite(model_settings.get(k)) for k in
                        ('uv_start_gev','fermion_threshold_gev','scalar_endpoint_gev')}
            require(all(v and v>0 for v in settings.values()),'Invalid reference scales')
            require(settings['uv_start_gev']>settings['fermion_threshold_gev']>
                    settings['scalar_endpoint_gev'], 'Invalid threshold ordering')
        else:
            for k,v in settings.items():
                require(math.isclose(finite(model_settings.get(k)) or -1,v,rel_tol=1e-9),
                        f'Inconsistent {k}: {name}')
        models[name] = {'stages':segs,'quadratic_kind':data.get('quadratic_kind'),
                        'admissible':data.get('scalar_endpoint_physically_admissible'),
                        'scalar_min_m2_gev2':finite((data.get('scalar_endpoint_spectrum') or {}).get('minimum_mass_squared_gev2'))}
        hashes[str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)] = digest(path)
    return models,settings


def convert_variant(model: str, scheme: str, row: dict) -> dict | None:
    if scheme=='baseline':
        keys=('largest_neutrino_mass_ev','middle_neutrino_mass_ev','smallest_neutrino_mass_ev')
    elif scheme=='self_massdep':
        if row.get('matching_category') != 'self-bilinear-one-scale-massdependent-conditional':
            return None  # Exclude ten pass-through controls.
        keys=('mass1_ev_descending','mass2_ev_descending','mass3_ev_descending')
    elif scheme=='direct':
        if str(row.get('direct_wilson_source_evolved')).lower() not in ('true','1'):
            return None
        keys=('max_mnu_ev','middle_mnu_ev','min_mnu_ev')
    elif scheme=='extra':
        keys=('updated_largest_mnu_ev','updated_middle_mnu_ev','updated_smallest_mnu_ev')
    elif scheme=='coupled':
        keys=('new_largest_mnu_ev','new_middle_mnu_ev','new_smallest_mnu_ev')
    else:
        raise ValueError('Unknown scheme: '+scheme)
    masses=[finite(row.get(k)) for k in keys]
    if any(x is None or x<0 for x in masses):
        return None
    require(masses[0]>=masses[1]>=masses[2], f'Non-descending {model}/{scheme} masses')
    result={'model':model,'source':scheme,'status':str(row.get('status','')),
            'mass_hi_ev':masses[0],'mass_mid_ev':masses[1], 'mass_lo_ev':masses[2],
            'gap_small_mid_sq_ev2':masses[1]**2-masses[2]**2,
            'gap_small_large_sq_ev2':masses[0]**2-masses[2]**2}
    for key in OBS[3:]:
        result[key] = finite(row.get(key))
    for key in EXTRA:
        result[key] = finite(row.get(key))
    for key in SCHEME_FLAGS:
        result[key] = row.get(key) or None
    return result


def read_variants(run: Path, hashes: dict[str,str]) -> tuple[list[dict],list[dict],list[str]]:
    rows: list[dict] = []
    absent=[]
    for scheme,rel,title,note in SCHEMES:
        path=run/rel
        if not path.is_file():
            if scheme=='baseline':
                raise FileNotFoundError('Missing required baseline: '+str(path))
            absent.append(f'{scheme}: {rel}')
            continue
        hashes[str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)] = digest(path)
        seen=set()
        for raw in read_csv(path):
            name=raw.get('model')
            if name not in MODELS:
                continue
            require(name not in seen,f'Duplicate model {name}: {path}')
            seen.add(name)
            item=convert_variant(name,scheme,raw)
            if item is not None:
                item.update(source_title=title,source_note=note,source_relative_path=rel)
                rows.append(item)
        if scheme=='baseline':
            require(set(MODELS)==seen, 'Baseline summary missing models: '+str(set(MODELS)-seen))
            require(len([v for v in rows if v['source']=='baseline'])==16,
                    'Baseline missing usable mass predictions')
    rows.sort(key=lambda d:(MODELS.index(d['model']),
                            next(i for i,s in enumerate(SCHEMES) if s[0]==d['source'])))
    preferred=[]
    for model in MODELS:
        variants={v['source']:v for v in rows if v['model']==model}
        if model in ('B_-1','C_-1'):
            priority=('coupled','extra','direct','baseline')
        elif model in ('A_+0','A_-2','D_+0','D_-2','E_+0','E_-2'):
            priority=('self_massdep','direct','baseline')
        else:
            priority=('direct','baseline')
        chosen=next((variants[x] for x in priority if x in variants),None)
        require(chosen is not None,'No saved mass results for '+model)
        preferred.append(dict(chosen,selection_reason='Presentation priority only; not certified best EFT'))
    return rows,preferred,absent


def extract_lambda_template(path: Path, variable: str, function: str|None=None) -> str:
    """Get the actual literal HTML of the existing Lambda generators (read-only)."""
    tree=ast.parse(path.read_text(encoding='utf-8-sig'))
    nodes=tree.body
    if function:
        nodes=next((n.body for n in tree.body
                    if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.name==function),[])
    for node in nodes:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id==variable
                                              for t in node.targets):
            if isinstance(node.value,ast.Constant) and isinstance(node.value.value,str):
                return node.value.value
    raise RuntimeError(f'Could not find {variable} HTML in {path}; Lambda generator changed')


def build_rge(models: dict, settings: dict, legacy_template: str) -> str:
    quartics=sorted({key for model in models.values()
                     for stage in model['stages'].values()
                     for row in stage for key in row.get('extra_quartics',{})})
    specifications=[('m1sq','mS1Sq_gev2',None,'mS1²','GeV²'),
                    ('m2sq','mS2Sq_gev2',None,'mS2²','GeV²')]
    for comp,label in [('real','Re'),('imag','Im'),('abs','| |')]:
        specifications.append((f'mu2_{comp}','mu2_gev2',comp,
                               f'{label} μ²','GeV²'))
    for k,key in enumerate(quartics):
        for comp,label in [('real','Re'),('imag','Im'),('abs','| |')]:
            specifications.append((f'q{k}_{comp}',key,comp,f'{label} {key}',''))
    payload={'models':list(MODELS),'scales':{
        'mu_uv':settings['uv_start_gev'],
        'mu_fermion_threshold':settings['fermion_threshold_gev'],
        'mu_scalar_threshold':settings['scalar_endpoint_gev'],
        'mu_low':settings['scalar_endpoint_gev']},
        'quantities':{},'matched':{},'components':{},'overlays':{},
        'continuation_methods':{},'continuation_issues':{}}
    meta={}
    for name,field,component,label,unit in specifications:
        meta[name]={'label':label,'unit':unit,'block':'saved_scalar','weinberg':False}
        curves={}
        for model,m in models.items():
            xs,ys=[],[]
            for stage in ('UV','EFT1'):
                for row in m['stages'][stage]:
                    val=(row.get('extra_quartics') or {}).get(field) if name.startswith('q') else row.get(field)
                    y=complex_part(val,component) if component else finite(val)
                    if y is not None:
                        xs.append(float(row['mu_gev']))
                        ys.append(y)
            if len(xs)>=2:
                curves[model]={'x':xs,'y':ys}
        if curves:
            payload['quantities'][name]=curves
    sections=(('Scalar mass squared',('m1sq','m2sq')),
              ('Quadratic bilinear',tuple(x for x in payload['quantities'] if x.startswith('mu2_'))),
              ('General-scalar quartics',tuple(x for x in payload['quantities'] if x.startswith('q'))))
    options=''.join('<optgroup label="'+html.escape(group,quote=True)+'">'+
                    ''.join('<option value="'+k+'">'+html.escape(meta[k]['label'])+'</option>'
                            for k in keys if k in payload['quantities'])+'</optgroup>'
                    for group,keys in sections)
    require('m1sq' in payload['quantities'] and 'm2sq' in payload['quantities'],
            'Could not construct saved scalar curves')
    replacements={
        '__TITLE__':'General-scalar RGE running — continuous_v1',
        '__OPTIONS__':options,'__DATA__':json_script(payload),'__META__':json_script(meta),
        'Numerical integration direction: UV → IR. Saved benchmark trajectories; no matching or RGE recalculation.':
            '16 saved general-scalar trajectories: UV → EFT1 → scalar endpoint. No physics recomputation; no C5(μ) or mixing-angle curves.',
        '<button id="zoomUV">Full T3</button>':'<button id="zoomUV">UV scalar</button>',
        '<button id="zoomInt">Intermediate EFT</button>':'<button id="zoomInt">Scalar EFT1</button>',
        '<span><i class="dot" style="background:#eaf2ff"></i>Full T3</span>':
            '<span><i class="dot" style="background:#eaf2ff"></i>UV scalar running</span>',
        '<span><i class="dot" style="background:#fff1cf"></i>Intermediate scalar EFT</span>':
            '<span><i class="dot" style="background:#fff1cf"></i>Scalar EFT1 running</span>',
        'rectRegion(s.mu_uv,s.mu_fermion_threshold,regionColours.uv,\'Full T3\');':
            'rectRegion(s.mu_uv,s.mu_fermion_threshold,regionColours.uv,\'UV scalar\');',
        'rectRegion(s.mu_fermion_threshold,s.mu_scalar_threshold,regionColours.int,\'Intermediate scalar EFT\');':
            'rectRegion(s.mu_fermion_threshold,s.mu_scalar_threshold,regionColours.int,\'Scalar EFT1\');',
    }
    result=legacy_template
    for old,new in replacements.items():
        require(old in result,f'Lambda RGE HTML changed: missing marker {old[:60]}')
        result=result.replace(old,new)
    result=result.replace('<button id="zoomFinal">SM + Weinberg</button>',
                          '<button id="zoomFinal" hidden>Post-scalar (not saved)</button>')
    result=result.replace('<span><i class="dot" style="background:#e5f5e9"></i>SM + Weinberg</span>','')
    result=result.replace('<header><h2>General-scalar RGE running',
                          '<header><a href="index.html" style="font-size:12px">Report index</a>'
                          '<h2>General-scalar RGE running',1)
    # Nothing to plot below the saved scalar endpoint; no synthetic continuation.
    require('const DATA='+json_script(payload)+';' in result, 'RGE data was not embedded')
    return result


def read_sensitivity(run: Path, csv_path: Path | None,
                     hashes: dict[str,str]) -> tuple[list[dict],list[str],str]:
    if csv_path is None:
        return ([], [], 'NO_GENERAL_SCALAR_GAMMA_SCAN')
    require(csv_path.is_file(),f'Sensitivity CSV not found: {csv_path}')
    # Reject historical scan inputs, even if the caller happens to supply them.
    require(csv_path.resolve().is_relative_to(run.resolve()),
            'Sensitivity CSV must be within the continuous_v1 run, not a historical reduced-potential study')
    rows=read_csv(csv_path)
    require(bool(rows),'Sensitivity CSV is empty')
    grouped: dict[tuple[str,str],list[dict]]={}
    for row in rows:
        model=row.get('model')
        if model not in MODELS:
            continue
        direction=row.get('direction')
        require(direction in ('T1','T2','T3'),f'Invalid gamma direction {direction}')
        gamma=finite(row.get('gamma'))
        require(gamma is not None and gamma>0, 'Gamma must be finite and positive')
        p={'gamma':gamma}
        for k in OBS:
            p[k]=finite(row.get(k))
        require(sum(p[k] is not None for k in OBS)>=2,
                'Gamma point requires at least two stored observables')
        grouped.setdefault((model,direction),[]).append(p)
    require(grouped,'No general-scalar models in gamma scan')
    present=[]
    skipped=[]
    for model in MODELS:
        bydir={}
        for d in ('T1','T2','T3'):
            pts=sorted(grouped.get((model,d),[]),key=lambda p:p['gamma'])
            if len(pts)<3 or any(math.isclose(a['gamma'],b['gamma'],rel_tol=1e-12)
                                for a,b in zip(pts,pts[1:])):
                break
            centre=[p for p in pts if math.isclose(p['gamma'],1,rel_tol=1e-9,abs_tol=1e-12)]
            if len(centre)!=1:
                break
            bydir[d]=pts
        if len(bydir)!=3:
            skipped.append(f'{model}: incomplete genuine T1/T2/T3 scan')
            continue
        centres=[next(p for p in bydir[d] if math.isclose(p['gamma'],1,rel_tol=1e-9,abs_tol=1e-12))
                 for d in ('T1','T2','T3')]
        obs=[k for k in OBS if all(p[k] is not None
                                  for d in ('T1','T2','T3') for p in bydir[d])]
        if len(obs)<2:
            skipped.append(f'{model}: fewer than two fully sampled observables')
            continue
        if any(not math.isclose(centres[0][k], c[k],rel_tol=1e-7,abs_tol=1e-10)
               for c in centres[1:] for k in obs):
            skipped.append(f'{model}: inconsistent gamma=1 centres')
            continue
        present.append({'model':model,'centre':centres[0],'directions':bydir,'obs':set(obs)})
    hashes[str(csv_path.relative_to(ROOT)) if csv_path.is_relative_to(ROOT) else str(csv_path)] = digest(csv_path)
    return present,skipped,'SAVED_GENERAL_SCALAR_GAMMA_CSV'


def build_sensitivity(scans: list[dict], skipped: list[str], state: str,
                      original: str, csv_path: Path | None) -> str:
    # The original HTML is the fitted-Lambda canvas plot, including gamma range,
    # T1/T2/T3 direction filters, sample paths and empirical covariance ellipses.
    # Only observables present in EVERY plotted scan may be selectable.
    common=set.intersection(*(s['obs'] for s in scans)) if scans else set()
    ordered=[k for k in OBS if k in common]
    if scans and len(ordered)<2:
        skipped=skipped+['Scans have no shared pair of fully saved observables.']
        scans=[]
    ordered=ordered if len(ordered)>=2 else list(OBS[3:])
    data={'models':[],'skipped':[{'model':s.split(':')[0], 'status':s.split(':',1)[-1].strip()}
                                 for s in skipped],
          'target':{},'labels':{k:LABELS[k] for k in ordered},'observables':ordered,
          'source':str(csv_path) if csv_path else 'continuous_v1 — no T1/T2/T3 scans',
          'future_sigma':{},'future_source':'Not available',
          'envelope_definition':'2x empirical sample covariance of saved T1/T2/T3 gamma points; descriptive, not confidence regions'}
    for s in scans:
        data['models'].append({'model':s['model'],'label':s['model'],
            'chi2':None,'fit_stage':'general_scalar_saved_scan','status':'SavedGeneralScalarScan',
            'centre':{k:s['centre'][k] for k in ordered},
            'directions':{d:[{'gamma':p['gamma'],**{k:p[k] for k in ordered}}
                             for p in s['directions'][d]] for d in ('T1','T2','T3')}})
    page=original.replace('__DATA__',json_script(data))
    require(page!=original,'Lambda sensitivity HTML data marker missing')
    page=page.replace('Individually fitted T3 model sensitivity',
                      'General-scalar T1/T2/T3 sensitivity — continuous_v1')
    page=page.replace('Each model uses its own fitted λT3, Yukawa matrices and other fitted parameters. Data from saved scans only — no RGE recomputation.',
                      'Saved general-scalar gamma scans only. No historical reduced-potential data; no RGE recomputation.')
    page=page.replace("let x='sin2_theta12',y='delta_m21_sq_ev2'",
                      "let x='"+ordered[0]+"',y='"+ordered[1]+"'")
    page=page.replace("+' (χ²='+Number(m.chi2).toPrecision(3)+')'", "+' (saved gamma scan)'")
    page=page.replace('NuFIT Gaussian 1σ/3σ', 'NuFIT contours (not supplied)')
    page=page.replace('JUNO 6-year projection (solar only)','JUNO projection (not supplied)')
    page=page.replace('<input type="checkbox" id="experimental" checked>',
                      '<input type="checkbox" id="experimental" disabled>')
    page=page.replace('<input type="checkbox" id="future" checked>',
                      '<input type="checkbox" id="future" disabled>')
    page=page.replace("D.models.length+' fitted models plotted; '",
                      "D.models.length+' saved general-scalar models plotted; '")
    explanation=('Verified input: '+html.escape(str(csv_path)) if scans else
                 'No genuine continuous_v1 T1/T2/T3 parameter scans are available to this build. '
                 'No sensitivity ovals have been manufactured. The matching-scheme differences '
                 'in the results table are NOT gamma sensitivity.')
    page=page.replace('<div class="layout"><aside>',
        '<div class="panel" style="margin:14px 20px 0;font-size:13px"><strong>Data status:</strong> '
        +explanation+'</div><div class="layout"><aside>')
    # Original note text identifies individually fitted scans; correct its physics framing.
    page=page.replace('about each fitted centre','about each saved γ=1 centre')
    page=page.replace("' fitted models plotted; '","' saved general-scalar models plotted; '")
    page=page.replace('Fitted centres','Saved γ=1 centres')
    page=page.replace('Fit trajectories','Frame trajectories')
    page=page.replace('<header><h2>General-scalar T1/T2/T3 sensitivity',
                      '<header><a href="index.html" style="font-size:12px">Report index</a>'
                      '<h2>General-scalar T1/T2/T3 sensitivity',1)
    return page


def write_csv(path: Path, items: list[dict]) -> None:
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=COLUMNS,extrasaction='ignore')
        writer.writeheader()
        writer.writerows(items)


CSS='''<style>:root{color-scheme:light;--ink:#17243b;--muted:#556375;--line:#d7dee8;--surface:#f7f9fc;--accent:#225ba8}
*{box-sizing:border-box}body{margin:0;background:var(--surface);color:var(--ink);font:13px system-ui,-apple-system,Segoe UI,Arial,sans-serif}
header{background:#fff;padding:22px 26px;border-bottom:1px solid var(--line)}h1{font-size:22px;margin:0 0 6px}h2{font-size:17px}
p{line-height:1.55}.small{color:var(--muted);font-size:12px}main{padding:22px 26px;max-width:1600px;margin:auto}
nav{display:flex;flex-wrap:wrap;gap:10px;margin:14px 0}a{color:#225ba8}nav a,.tile{background:#fff;border:1px solid var(--line);border-radius:8px;padding:12px;text-decoration:none}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px}.tile{display:block}.tile strong{display:block;font-size:15px;margin-bottom:8px}
.panel{background:#fff;border:1px solid var(--line);border-radius:10px;padding:15px;margin:14px 0}
input,select,button{border:1px solid #bac7d8;border-radius:6px;background:white;color:var(--ink);padding:8px;font:inherit}
input{min-width:210px}button{cursor:pointer}table{width:100%;border-collapse:collapse;font-size:12px;white-space:nowrap}
th,td{padding:8px 10px;border-bottom:1px solid #e1e7ef;text-align:left}th{position:sticky;top:0;background:#f1f5fb;cursor:pointer;user-select:none}
tr:hover{background:#f6f9fe}.scroll{max-width:100%;overflow:auto;max-height:72vh}.muted{color:var(--muted)}
@media(max-width:700px){main{padding:10px}header{padding:15px}}
</style>'''


def page(title: str, body: str) -> str:
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            '<title>'+html.escape(title)+'</title>'+CSS+'</head><body>'+body+'</body></html>')


def build_index(preferred: list[dict], variants: list[dict], scans: list[dict],
                absent: list[str]) -> str:
    cards=(('interactive_rge.html','Interactive RGE running',
            'All 16 saved scalar UV/EFT1 running curves; exact sampled scales, threshold and quantity selectors.'),
           ('sensitivity_report.html','T1/T2/T3 sensitivity ovals',
            f'{len(scans)} models with actual saved gamma scans. No historical sensitivity substituted.'),
           ('all_results.html','All numerical results',
            f'16 default selections and {len(variants)} distinct saved scheme records, fully sortable.'))
    header='<header><h1>T3 general-scalar · continuous_v1</h1><p class="small">Lambda-style offline HTML reports, assembled from saved numerical outputs only.</p></header>'
    nav='<nav>'+''.join(f'<a href="{p}">{html.escape(t)}</a>' for p,t,_ in cards)+'</nav>'
    tiles='<div class="grid">'+''.join(f'<a class="tile" href="{p}"><strong>{html.escape(t)}</strong><span class="small">{html.escape(desc)}</span></a>' for p,t,desc in cards)+'</div>'
    text='<p>Default selections are presentation priorities, not certified best predictions. Distinct matching approximations are never added together. Time-resolved neutrino mixing and physical C5(μ) were not stored in the scalar trajectories.</p>'
    if absent:
        text+='<p class="small">Missing optional saved sources: '+html.escape('; '.join(absent))+'</p>'
    return page('General-scalar report index',header+'<main>'+nav+tiles+
                '<div class="panel"><h2>Interpretation</h2>'+text+
                '<p><a href="all_results.csv">16-model CSV</a> · '
                '<a href="matching_variants.csv">All matching variants CSV</a> · '
                '<a href="provenance.json">Input provenance / hashes</a></p></div></main>')


def build_table(preferred: list[dict], variants: list[dict], absent: list[str]) -> str:
    opts=[('all','All matching schemes')]+[(code,title) for code,_,title,_ in SCHEMES]
    options=''.join('<option value="'+code+'">'+html.escape(title)+'</option>' for code,title in opts)
    cols=['model','source','status','quadratic_kind','scalar_endpoint_admissible','scalar_min_m2_gev2',
          *OBS,'gap_small_mid_sq_ev2','gap_small_large_sq_ev2',*EXTRA,*SCHEME_FLAGS]
    colheads={'model':'Model','source':'Matching source','status':'Status',
              'quadratic_kind':'Quadratic class',
              'scalar_endpoint_admissible':'Scalar endpoint admissible',
              'scalar_min_m2_gev2':'Min scalar mass² [GeV²]',**LABELS,
              **{k:k.replace('_',' ') for k in EXTRA},
              **{k:k.replace('_',' ') for k in SCHEME_FLAGS}}
    script='''<script>'use strict';const RECORDS=__RECORDS__;const PREFERRED=__PREFERRED__;const COLS=__COLS__;
let order='model',ascending=true;const body=document.getElementById('rows'),search=document.getElementById('search'),
    scheme=document.getElementById('scheme'),mode=document.getElementById('mode');
const fmt=x=>x===null||x===undefined||x===''?'—':typeof x==='number'?(Number.isFinite(x)?x.toPrecision(7):'—'):String(x);
function render(){let items=(mode.value==='default'?PREFERRED:RECORDS).filter(r=>
(!search.value||[r.model,r.source,r.status].join(' ').toLowerCase().includes(search.value.toLowerCase()))&&
(scheme.value==='all'||r.source===scheme.value));
items.sort((a,b)=>{let x=a[order],y=b[order];if(x===null||x===undefined)return y===null||y===undefined?0:1;
if(y===null||y===undefined)return -1;let k=typeof x==='number'&&typeof y==='number'?x-y:String(x).localeCompare(String(y));return ascending?k:-k});
body.replaceChildren();for(const row of items){const tr=document.createElement('tr');for(const k of COLS){const td=document.createElement('td');td.textContent=fmt(row[k]);if(row[k]===null||row[k]===undefined)td.className='muted';
if(k==='source')td.title=row.source_note||'';tr.append(td)}body.append(tr)}document.getElementById('count').textContent=items.length+' records';}
document.querySelectorAll('th[data-key]').forEach(th=>th.onclick=()=>{const k=th.dataset.key;
if(order===k)ascending=!ascending;else{order=k;ascending=true}render()});[search,scheme,mode].forEach(e=>e.oninput=render);render();</script>'''
    script=script.replace('__RECORDS__',json_script(variants)).replace('__PREFERRED__',json_script(preferred)).replace('__COLS__',json_script(cols))
    thead='<tr>'+''.join('<th data-key="'+k+'" title="Click to sort">'+html.escape(colheads[k])+'</th>' for k in cols)+'</tr>'
    heading='<header><h1>General-scalar numerical results</h1><p class="small">All 16 models; central predictions and independent matching-scheme variants. A dash means this scheme did not save the observable.</p></header>'
    controls='<div class="panel"><select id="mode"><option value="default">One preferred row per model (16)</option><option value="variants">All matching-scheme records</option></select> <select id="scheme">'+options+'</select> <input type="search" id="search" placeholder="Filter model / scheme / status"> <span class="small" id="count"></span></div>'
    expl='<p class="small">Masses are sorted descending, not automatically labelled physical m₁/m₂/m₃. The two derived squared-mass gaps are algebraic mass differences, not replacements for saved Δm²₂₁ and Δm²₃ℓ. Mixing angles from one scheme are never copied into another. Rows from different sources are not added.</p>'
    if absent:
        expl+='<p class="small">Optional files not found: '+html.escape('; '.join(absent))+'</p>'
    return page('All general-scalar model results',heading+'<main><nav><a href="index.html">Report index</a><a href="interactive_rge.html">RGE</a><a href="sensitivity_report.html">Sensitivity</a></nav>'+expl+controls+'<div class="panel scroll"><table><thead>'+thead+'</thead><tbody id="rows"></tbody></table></div><p><a href="all_results.csv">Download 16-model CSV</a> · <a href="matching_variants.csv">Download variant CSV</a></p></main>'+script)


def build(run: Path, out: Path, sensitivity_csv: Path|None=None) -> dict:
    hashes={}
    models,scales=read_trajectories(run,hashes)
    variants,preferred,absent=read_variants(run,hashes)
    for item in [*variants,*preferred]:
        m=models[item['model']]
        item.update(quadratic_kind=m['quadratic_kind'],
                    scalar_endpoint_admissible=m['admissible'],
                    scalar_min_m2_gev2=m['scalar_min_m2_gev2'])
    scans,skipped,scan_state=read_sensitivity(run,sensitivity_csv,hashes)
    rge_source=ROOT/'Numerical/plotting/InteractiveModelComparison.py'
    sensitivity_source=ROOT/'Numerical/plotting/BuildFittedModelEnvelopes.py'
    rge_html=extract_lambda_template(rge_source,'template','write_dashboard')
    sens_html=extract_lambda_template(sensitivity_source,'HTML')
    hashes[str(rge_source.relative_to(ROOT))]=digest(rge_source)
    hashes[str(sensitivity_source.relative_to(ROOT))]=digest(sensitivity_source)
    outputs={
        'index.html':build_index(preferred,variants,scans,absent),
        'interactive_rge.html':build_rge(models,scales,rge_html),
        'sensitivity_report.html':build_sensitivity(scans,skipped,scan_state,sens_html,sensitivity_csv),
        'all_results.html':build_table(preferred,variants,absent),
    }
    manifest={'run_root':str(run),'input_sha256':hashes,
              'model_count':len(models),'preferred_model_count':len(preferred),
              'variant_count':len(variants),'variant_counts_by_scheme':
                  {k:sum(v['source']==k for v in variants) for k,_,_,_ in SCHEMES},
              'sensitivity_status':scan_state,'sensitivity_model_count':len(scans),
              'sensitivity_skipped':skipped,'optional_sources_not_found':absent,
              'rge_saved_segments':['UV','EFT1'],
              'limitations':['no saved physical time-dependent C5 or PMNS from scalar trajectories',
                             'matching candidates have nonidentical approximations',
                             'no historical reduced-potential sensitivity inputs permitted'],
              'generator':'Numerical/plotting/BuildGeneralScalarLambdaReports.py'}
    # Build all content in memory before modifying the destination.
    out.mkdir(parents=True,exist_ok=True)
    for name,content in outputs.items():
        (out/name).write_text(content,encoding='utf-8')
    write_csv(out/'all_results.csv',preferred)
    write_csv(out/'matching_variants.csv',variants)
    (out/'provenance.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    sensitivity_copy=out/'parameter_sensitivity_points.csv'
    if sensitivity_csv is not None:
        sensitivity_copy.write_bytes(sensitivity_csv.read_bytes())
    elif sensitivity_copy.exists():
        sensitivity_copy.unlink()  # Remove stale old display artifact, not physics data.
    return manifest


def main(argv: list[str]|None=None) -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-id',default='continuous_v1')
    ap.add_argument('--input',type=Path,help='Override saved run directory')
    ap.add_argument('--output',type=Path,help='Override output report directory')
    ap.add_argument('--sensitivity-csv',type=Path,help='Real general-scalar scan CSV INSIDE the saved run')
    args=ap.parse_args(argv)
    run=(args.input or ROOT/BASE/args.run_id).resolve()
    out=(args.output or ROOT/OUTPUT/args.run_id).resolve()
    require(run.is_dir(),f'Saved numerical directory does not exist: {run}')
    result=build(run,out,args.sensitivity_csv.resolve() if args.sensitivity_csv else None)
    print(f"Wrote {out} | {result['model_count']} models | {result['variant_count']} matching rows | "
          f"{result['sensitivity_model_count']} general-scalar gamma scans")
    return 0

if __name__=='__main__':
    raise SystemExit(main())
