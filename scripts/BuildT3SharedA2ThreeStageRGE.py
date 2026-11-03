"""Add actual SMEFT C5/neutrino-mass running to saved general-T3 RGE reports.

This additive, local-only renderer reads the shared-A2 gamma=1 centres and the
original saved UV/EFT1 numerical RGE traces. It integrates ONLY the final
SM+Weinberg segment using the repository's unmodified SMWeinbergEvolution beta.

This is NOT an EFT1 Wilson RGE: the complete physical LLSS anomalous-dimension
basis and UV->EFT1->SMEFT matching are not certified in the saved general study.
Never draw an EFT1 C5 curve from these one-scale endpoint coefficients.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
STUDY_ROOT = ROOT / 'output' / 'studies' / 'general_scalar_extended' / 'shared_a2_general16_reports'
FORMAT = 't3-general-shared-a2-sm-weinberg-sampled-v1'
DEFAULT_RUN = 'shared_A2_general16_gamma08_12_v1'
SMEFT_KEYS = {
    'c5_fnorm': '||C5||F [GeV^-1]',
    'c5_ee_abs': '|C5_ee| [GeV^-1]',
    'c5_emu_abs': '|C5_eμ| [GeV^-1]',
    'c5_mumu_abs': '|C5_μμ| [GeV^-1]',
    'c5_tautau_abs': '|C5_ττ| [GeV^-1]',
    'm1_ev': 'm1 [eV]', 'm2_ev': 'm2 [eV]', 'm3_ev': 'm3 [eV]',
    'delta_m21_sq_ev2': 'Δm21² [eV²]',
    'delta_m3l_sq_ev2': 'Δm31² [eV²], NO',
    'sin2_theta12': 'sin²θ12', 'sin2_theta13': 'sin²θ13',
    'sin2_theta23': 'sin²θ23',
    'gY': 'gY', 'g2': 'g2', 'g3': 'g3',
    'lambdaH': 'λH',
}


def load(path: Path) -> dict:
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(data, dict):
        raise ValueError('Expected JSON object: ' + str(path))
    return data


def save(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n',
                    encoding='utf-8')


def matrix(raw: Any) -> np.ndarray:
    """Decode compact {real,imag} or saved per-entry complex 3x3 matrices."""
    if isinstance(raw, dict) and 'real' in raw and 'imag' in raw:
        value = np.asarray(raw['real'], dtype=float) + 1j * np.asarray(raw['imag'], dtype=float)
    elif isinstance(raw, list):
        value = np.asarray([[complex(float(x['real']), float(x['imag']))
                             if isinstance(x, dict) else complex(x) for x in row]
                            for row in raw], dtype=complex)
    else:
        raise ValueError('Unexpected matrix representation')
    if value.shape != (3, 3) or not np.isfinite(value).all():
        raise ValueError('Expected finite complex 3x3 matrix')
    return value


def real(raw: Any) -> float:
    if isinstance(raw, dict):
        z = complex(float(raw['real']), float(raw['imag']))
    else:
        z = complex(raw)
    if not math.isfinite(z.real) or not math.isfinite(z.imag) or abs(z.imag) > 1e-11 * max(abs(z.real), 1.):
        raise ValueError('Expected a finite real-valued renormalisable parameter')
    return float(z.real)


def sample_smeft(centre: dict, n: int, *, rtol: float = 1e-8, atol: float = 1e-11) -> dict:
    """Run the already-installed full-flavour one-loop SMEFT ODE at actual t_eval points."""
    from scipy.integrate import solve_ivp
    from Numerical.running.SMWeinbergEvolution import (
        SMWeinbergInitialConditions, _pack, _unpack, _beta,
    )
    from physics.NeutrinoTrajectory import (
        neutrino_mass_from_c5, charged_lepton_mass_basis_matrix,
    )
    from physics.NeutrinoObservables import calculate_neutrino_observables

    if n < 3 or n > 501:
        raise ValueError('--points must be 3..501')
    if centre.get('status') != 'Success':
        raise ValueError('Centre not marked successful')
    detail = centre.get('numerical_details')
    if not isinstance(detail, dict):
        raise ValueError('No saved gamma=1 detailed boundary')
    if detail.get('matching_scope') != 'full-general Matchete quadratic Yukawa channels, one-scale':
        raise ValueError('Only nonmixed general one-scale C5 centres supported')

    eft = detail['eft1_at_scalar_threshold']
    initial = SMWeinbergInitialConditions(
        gY=real(eft['gY']), g2=real(eft['g2']), g3=real(eft['g3']),
        lambdaH=real(eft['lambdaH']),
        ye=matrix(eft['ye']), yu=matrix(eft['yu']), yd=matrix(eft['yd']),
        K=matrix(detail['c5_at_scalar_threshold']),
    ).validated()
    mu_s, mu_low = 1e3, 1e2  # fixed in verified shared-A2 general fitter
    t = np.linspace(math.log(mu_s), math.log(mu_low), n)
    solution = solve_ivp(_beta, (float(t[0]), float(t[-1])), _pack(initial),
                         method='DOP853', t_eval=t, rtol=rtol, atol=atol)
    if not solution.success or solution.y.shape != (76, n):
        raise RuntimeError('SMEFT ODE failed: ' + str(solution.message))

    all_rows = []
    for k in range(n):
        gY, g2, g3, lam, ye, yu, yd, K = _unpack(solution.y[:, k])
        K = 0.5*(K + K.T)
        if not np.isfinite(K).all():
            raise ValueError('Nonfinite SMEFT C5 sample')
        mass = neutrino_mass_from_c5(K, vev_gev=246.22)
        flavour = charged_lepton_mass_basis_matrix(mass, ye)
        o = calculate_neutrino_observables(flavour, ordering='NO')
        if not np.isfinite(o.masses_ev).all() or not math.isfinite(o.takagi_residual):
            raise ValueError('Invalid Takagi output')
        if o.takagi_residual > 1e-7:
            raise ValueError('Takagi residual too large at SMEFT sample')
        pmns = o.pmns_abs
        den = 1. - float(pmns[0, 2]**2)
        if den <= 0.:
            raise ValueError('Invalid PMNS denominator')
        row = {
            'mu_gev': float(math.exp(t[k])),
            'c5_fnorm': float(np.linalg.norm(K)),
            'c5_ee_abs': float(abs(K[0, 0])), 'c5_emu_abs': float(abs(K[0, 1])),
            'c5_mumu_abs': float(abs(K[1, 1])), 'c5_tautau_abs': float(abs(K[2, 2])),
            'm1_ev': float(o.masses_ev[0]), 'm2_ev': float(o.masses_ev[1]),
            'm3_ev': float(o.masses_ev[2]),
            'delta_m21_sq_ev2': float(o.delta_m21_sq_ev2),
            'delta_m3l_sq_ev2': float(o.delta_m31_sq_ev2),
            'sin2_theta12': float(pmns[0, 1]**2/den),
            'sin2_theta13': float(pmns[0, 2]**2),
            'sin2_theta23': float(pmns[1, 2]**2/den),
            'gY': float(gY), 'g2': float(g2), 'g3': float(g3), 'lambdaH': float(lam),
        }
        if any(not math.isfinite(value) for value in row.values()):
            raise ValueError('Nonfinite computed SMEFT observable')
        all_rows.append(row)

    saved_low = matrix(detail['c5_at_low_scale'])
    K_out = _unpack(solution.y[:, -1])[-1]
    mismatch = float(np.linalg.norm(K_out - saved_low) / max(np.linalg.norm(saved_low), 1e-99))
    if mismatch > 3e-6:
        raise ValueError(f'SMEFT C5 endpoint fails saved centre parity: {mismatch:.3e}')
    pred = centre.get('prediction') or {}
    for key in ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23'):
        if key in pred and not math.isclose(all_rows[-1][key], float(pred[key]),
                                             rel_tol=3e-5, abs_tol=2e-10):
            raise ValueError(f'SMEFT observable parity mismatch for {key}')

    return {
        'format': FORMAT, 'model': centre['model'], 'status': 'Success',
        'scheme': 'full general one-scale C5 boundary + existing full-flavour SMEFT RGEs',
        'mu_gev': [r.pop('mu_gev') for r in all_rows],
        'series': {key: [r[key] for r in all_rows] for key in SMEFT_KEYS},
        'nfev': int(solution.nfev), 'relative_endpoint_c5_residual': mismatch,
        'normal_ordering': True, 'vev_gev_fixed': 246.22,
        'notes': [
            'All plotted samples are numerical ODE samples, not endpoint interpolation.',
            'mν(μ) = -v² C5(μ)/2 with v held fixed, as in the existing fitted pipeline.',
            'Saved C5 boundary originates from conditional full-general one-scale Matchete matching.',
            'This stage does not certify missing physical EFT1 Wilson mixing or sequential matching.',
            'The Higgs quadratic mass-squared mH² is not part of the current SMWeinbergEvolution state.',
        ]
    }


def combine(uv_eft1: dict, smeft: dict) -> dict:
    if uv_eft1.get('status') != 'Success' or smeft.get('status') != 'Success':
        raise ValueError('Only successful actual trajectories may be combined')
    return {'uv': uv_eft1['uv'], 'eft1': uv_eft1['eft1'], 'smeft': smeft}


_HTML = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>T3 general UV / EFT1 / SMEFT RGE</title>
<style>
:root{color-scheme:light}body{font:15px system-ui,-apple-system,Segoe UI,sans-serif;max-width:1100px;margin:0 auto;padding:26px;color:#1d2634;background:#fff}
h1{font-size:26px;margin-bottom:7px}p{line-height:1.5}.muted{color:#556473}.controls{display:flex;flex-wrap:wrap;gap:14px;margin:22px 0}
label{display:flex;flex-direction:column;gap:5px;font-weight:600}select{font:inherit;border:1px solid #aeb8c5;border-radius:6px;background:#fff;padding:8px;min-width:170px}
svg{display:block;width:100%;height:auto;background:#fafcff;border:1px solid #dae1e8;border-radius:9px}
.note{background:#f3f5f7;border-left:3px solid #8899aa;padding:12px 16px;margin-top:20px}
nav a{margin-right:15px}.footer{font-size:12px;margin-top:15px;color:#4b5563}
</style></head><body>
<h1>General T3 — running across three EFT stages</h1>
<p class="muted">Saved shared A(−2) Yukawa benchmark, γ=1. Curves are computed numerical solver samples. Full-general <em>one-scale</em> Matchete C5 matching remains conditional.</p>
<nav><a href="interactive_RGE.html">Original UV/EFT1 dashboard</a><a href="sensitivity_oval.html">Sensitivity ovals</a><a href="comparison.csv">Model comparisons</a></nav>
<div class="controls"><label>Model<select id="model"></select></label><label>EFT stage<select id="stage"><option value="uv">UV (10⁷ → 10⁵ GeV)</option><option value="eft1">EFT1 (10⁵ → 10³ GeV)</option><option value="smeft" selected>SMEFT (10³ → 10² GeV)</option></select></label><label>Parameter<select id="series"></select></label><label>Vertical scale<select id="yscale"><option value="linear">Linear</option><option value="log">Log absolute value</option></select></label></div>
<svg id="plot" viewBox="0 0 1010 510" role="img" aria-label="Scale-resolved renormalization group trajectory"></svg>
<p class="muted" id="info"></p>
<div class="note"><strong>EFT1 C5 is not plotted.</strong> Its physically complete LLSS operator mixing, including zero-boundary components and any antilinear C* contribution, is not certified for this general benchmark. No inferred C5 curve is added between fermion and scalar thresholds. In SMEFT, m₁,m₂,m₃ and Δm² are obtained from the running full-flavour C5 using a fixed v=246.22 GeV; these are scale-dependent running mass parameters, not pole-mass predictions. Higgs mH² is not part of the current SMEFT solver state.</div>
<p class="footer">Report generated locally from saved centre.json and rge_trajectories.json. B(−1)/C(−1) are not included because the existing dense UV/EFT1 general solver does not export equivalent samples. No fits, symbolic matching, or source data were modified.</p>
<script>
const DATA = __DATA__;
const LABELS = __LABELS__;
const model = document.getElementById('model');
const stage = document.getElementById('stage');
const series = document.getElementById('series');
const yscale = document.getElementById('yscale');
const svg = document.getElementById('plot');
const info = document.getElementById('info');
const ns = 'http://www.w3.org/2000/svg';
function node(tag, props={}, str=''){let e=document.createElementNS(ns,tag);for(let [k,v] of Object.entries(props))e.setAttribute(k,String(v));if(str)e.textContent=str;svg.appendChild(e);return e;}
function fmt(n){return Number(n).toExponential(3);}
function populateModels(){Object.keys(DATA).sort().forEach(m=>{let o=document.createElement('option');o.value=m;o.textContent=m.replace('_','(')+')';model.appendChild(o);});if(DATA['A_m2'])model.value='A_m2';}
function populateSeries(){let p=series.value;series.replaceChildren();let names=Object.keys(DATA[model.value][stage.value].series).sort();names.forEach(k=>{let opt=document.createElement('option');opt.value=k;opt.textContent=LABELS[k]||k;series.appendChild(opt);});if(names.includes(p))series.value=p;else if(stage.value==='smeft'&&names.includes('m3_ev'))series.value='m3_ev';else if(names.includes('lambdaH'))series.value='lambdaH';}
function redraw(){svg.replaceChildren();let d=DATA[model.value][stage.value];if(!d)return;let xs=d.mu_gev, ys=d.series[series.value];let logY=yscale.value==='log';let clean=ys.map(v=>logY?Math.log10(Math.abs(v)):v);let valid=clean.filter(Number.isFinite);if(valid.length<2){info.textContent='Unavailable on selected vertical scale';return;}
let ymin=Math.min(...valid),ymax=Math.max(...valid);if(ymin===ymax){ymin-=Math.abs(ymin)*.1+1e-12;ymax+=Math.abs(ymax)*.1+1e-12;}let pad=(ymax-ymin)*.1; ymin-=pad; ymax+=pad;
let left=115,right=970,top=50,bottom=436;let xmin=Math.min(...xs.map(Math.log10)),xmax=Math.max(...xs.map(Math.log10));
const px=x=>left+(Math.log10(x)-xmin)/(xmax-xmin)*(right-left),py=y=>bottom-(y-ymin)/(ymax-ymin)*(bottom-top);
for(let i=0;i<=5;i++){let y=ymin+(ymax-ymin)*i/5,Y=py(y);node('line',{x1:left,y1:Y,x2:right,y2:Y,stroke:'#dce3eb'});node('text',{x:left-10,y:Y+5,'text-anchor':'end','font-size':13,fill:'#425466'},logY?('10^'+y.toFixed(2)):fmt(y));}
for(let i=0;i<=4;i++){let x=xmin+(xmax-xmin)*i/4,X=left+(right-left)*i/4;node('line',{x1:X,y1:top,x2:X,y2:bottom,stroke:'#e4e9f0'});node('text',{x:X,y:bottom+23,'text-anchor':'middle','font-size':13,fill:'#425466'},'10^'+x.toFixed(2));}
node('line',{x1:left,y1:bottom,x2:right,y2:bottom,stroke:'#45566b','stroke-width':1.6});node('line',{x1:left,y1:top,x2:left,y2:bottom,stroke:'#45566b','stroke-width':1.6});
let pts=xs.map((x,i)=>Number.isFinite(clean[i])?`${px(x)},${py(clean[i])}`:null);let runs=[],run=[];for(let v of pts){if(v)run.push(v);else{if(run.length>1)runs.push(run);run=[];}}if(run.length>1)runs.push(run);runs.forEach(r=>node('polyline',{points:r.join(' '),fill:'none',stroke:'#1e679a','stroke-width':3,'stroke-linejoin':'round'}));
node('text',{x:545,y:30,'text-anchor':'middle','font-size':20,'font-weight':600,fill:'#172b44'},LABELS[series.value]||series.value);
node('text',{x:545,y:488,'text-anchor':'middle','font-size':15,fill:'#34495e'},'Renormalisation scale μ [GeV] (log₁₀ axis)');
info.textContent=`${model.value} · ${stage.value.toUpperCase()} · ${ys.length} genuine saved ODE samples · μ=${fmt(xs[0])} to ${fmt(xs[xs.length-1])} GeV · ${LABELS[series.value]||series.value}: ${fmt(ys[0])} → ${fmt(ys[ys.length-1])}`;
}
model.addEventListener('change',()=>{populateSeries();redraw();});stage.addEventListener('change',()=>{populateSeries();redraw();});series.addEventListener('change',redraw);yscale.addEventListener('change',redraw);populateModels();populateSeries();redraw();
</script></body></html>'''


def render_html(combined: dict) -> str:
    # The JSON is embedded in script data, so guard against literal script termination.
    encoded = json.dumps(combined, ensure_ascii=False, allow_nan=False).replace('</', '<\\/')
    labels = json.dumps(SMEFT_KEYS, ensure_ascii=False).replace('</', '<\\/')
    return _HTML.replace('__DATA__', encoded).replace('__LABELS__', labels)


def run(run_id: str, n: int, *, output_dir: Path | None = None) -> dict:
    folder = STUDY_ROOT / run_id
    if output_dir is None:
        output_dir = folder
    if not folder.is_dir():
        raise FileNotFoundError('Missing saved shared-A2 run: ' + str(folder))
    if not (folder / 'manifest.json').is_file():
        raise FileNotFoundError('Missing saved numerical report manifest')
    manifest = load(folder / 'manifest.json')
    if manifest.get('interactive_rge_generated') is not True:
        raise ValueError('Saved general run does not declare genuine UV/EFT1 RGE curves')
    rows = []
    combined = {}
    for centre_path in sorted(folder.glob('*/centre.json')):
        dest = centre_path.parent
        name = dest.name
        path = dest / 'rge_trajectories.json'
        if not path.is_file():
            rows.append({'model': name, 'status': 'Unavailable',
                         'reason': 'No genuine dense UV/EFT1 samples; no interpolation used'})
            continue
        try:
            rge = load(path)
            centre = load(centre_path)
            smeft = sample_smeft(centre, n)
            combined[name] = combine(rge, smeft)
            save(output_dir/name/'smeft_trajectories.json', smeft)
            with (output_dir/name/'smeft_observables.csv').open('w', newline='', encoding='utf-8') as f:
                w = csv.writer(f)
                columns = ['mu_gev', *SMEFT_KEYS]
                w.writerow(columns)
                for i, mu in enumerate(smeft['mu_gev']):
                    w.writerow([mu, *[smeft['series'][k][i] for k in SMEFT_KEYS]])
            rows.append({'model': name, 'status': 'Success', 'smeft_samples': n,
                         'endpoint_parity': smeft['relative_endpoint_c5_residual']})
        except Exception as exc:
            reason = f'{type(exc).__name__}: {exc}'
            save(output_dir/name/'smeft_error.json', {'model': name, 'status': 'Failed', 'reason': reason})
            rows.append({'model': name, 'status': 'Failed', 'reason': reason})
    if combined:
        (output_dir/'interactive_RGE_with_SMEFT.html').write_text(render_html(combined), encoding='utf-8')
    report = {
        'format': FORMAT, 'run_id': run_id, 'status': 'Success' if combined else 'NoSMEFTPlots',
        'available_models': list(combined), 'model_status': rows,
        'eft1_c5_running_calculated': False,
        'smeft_c5_and_neutrino_masses_numerically_sampled': bool(combined),
        'scalar_threshold_c5_matching_physical_sequential_certified': False,
        'original_reports_preserved': True,
        'new_html': str(output_dir/'interactive_RGE_with_SMEFT.html') if combined else None,
    }
    save(output_dir/'smeft_status.json', report)
    return report


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id', default=DEFAULT_RUN)
    p.add_argument('--points', type=int, default=41)
    p.add_argument('--output-dir', type=Path, help='Optional separate destination (source run unchanged)')
    a = p.parse_args(argv)
    if not a.run_id or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in a.run_id):
        p.error('--run-id must contain only letters, digits, underscore or hyphen')
    report = run(a.run_id, a.points, output_dir=a.output_dir)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if report['available_models'] and all(r['status'] != 'Failed' for r in report['model_status']) else 2

if __name__ == '__main__':
    raise SystemExit(main())
