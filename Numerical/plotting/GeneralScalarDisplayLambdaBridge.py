"""Adapt saved general-scalar results to *existing* DisplayLambda HTML generators.

Uses InteractiveModelComparison.write_dashboard and, when a fully sampled
five-observable T1/T2/T3 scan exists, BuildFittedModelEnvelopes.build.
Never computes RGEs, matches a threshold, fits parameters or invents scan data.
"""
from __future__ import annotations

from collections import defaultdict
import json
import math
from pathlib import Path
import re
from tempfile import TemporaryDirectory
from typing import Any, Mapping


CLASS_DIMS = {
    'A': (1, 3, 2), 'B': (2, 2, 1), 'C': (2, 2, 3),
    'D': (3, 1, 2), 'E': (3, 3, 2),
}
GAMMA_OBSERVABLES = (
    'delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
    'sin2_theta12', 'sin2_theta13', 'sin2_theta23',
)


def _identifier(model: str) -> str:
    m = re.fullmatch(r'([ABCDE])_([+-])(\d+)', model)
    if not m:
        raise ValueError(f'Unsupported general-scalar model {model!r}')
    cls, sign, alpha = m.groups()
    d1, d2, df = CLASS_DIMS[cls]
    return f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{"p" if sign == "+" else "m"}{alpha}'


def _finite_complex(z: Any) -> tuple[float, float] | None:
    if not isinstance(z, Mapping):
        return None
    try:
        x, y = float(z['real']), float(z['imag'])
    except (KeyError, TypeError, ValueError):
        return None
    return (x, y) if math.isfinite(x) and math.isfinite(y) else None


def scalar_curves(trajectories: Mapping[str, dict]) -> tuple[dict[str, dict], dict[str, tuple]]:
    """Only use actual scalar samples; separate UV and EFT1 explicitly."""
    if not trajectories:
        raise ValueError('No saved scalar trajectories')
    metadata: dict[str, tuple] = {}
    payloads: dict[str, dict] = {}
    for model, obj in sorted(trajectories.items()):
        by_stage = {'uv': [], 'eft1': []}
        for row in obj['samples']:
            stage = str(row['stage']).lower()
            if stage not in by_stage:
                raise ValueError(f'Unknown trajectory stage: {stage}')
            mu = float(row['mu'])
            if not math.isfinite(mu) or mu <= 0:
                raise ValueError('Nonpositive scalar trajectory scale')
            by_stage[stage].append(row)
        if any(len(seq) < 2 for seq in by_stage.values()):
            raise ValueError(f'{model}: missing UV or EFT1 samples')
        uv, eft = by_stage['uv'], by_stage['eft1']
        if not math.isclose(uv[-1]['mu'], eft[0]['mu'], rel_tol=1e-8):
            raise ValueError(f'{model}: UV/EFT1 threshold mismatch')
        scalars = {'uv_running': {}, 'final_running': {}}
        for stage, block in (('uv', 'uv_running'), ('eft1', 'final_running')):
            rows = by_stage[stage]
            scalars[block]['mu_gev'] = [float(row['mu']) for row in rows]
            label = 'UV' if stage == 'uv' else 'EFT1'
            for field, title in (('mS1Sq', 'mS1² [GeV²]'), ('mS2Sq', 'mS2² [GeV²]')):
                key = f'{stage}_{field}'
                vals = [float(row[field]) for row in rows]
                if not all(math.isfinite(v) for v in vals):
                    raise ValueError(f'{model}: nonfinite {key}')
                scalars[block][key] = vals
                metadata[key] = (f'{label}: {title}', block, key, 'GeV²', False)
            pairs = {'mu2': [(_finite_complex(row.get('mu2'))) for row in rows]}
            quartic_names = sorted(set().union(*(set(row.get('quartics', {})) for row in rows)))
            pairs.update({name: [_finite_complex(row.get('quartics', {}).get(name)) for row in rows]
                          for name in quartic_names})
            for name, zs in pairs.items():
                if any(z is None for z in zs):
                    # A missing per-model coupling is not a numerical zero.
                    continue
                units = 'GeV²' if name == 'mu2' else ''
                for comp in ('abs', 'real', 'imag'):
                    key = f'{stage}_{name}_{comp}'
                    scalars[block][key] = [math.hypot(*z) if comp == 'abs'
                                           else z[0] if comp == 'real' else z[1] for z in zs]
                    metadata[key] = (f'{label}: {name} ({comp})', block, key, units, False)
        payloads[model] = {
            'status': 'Success',
            'benchmark_search': {'model': {'model_key': _identifier(model)}},
            'scales_gev': {
                'mu_uv': float(uv[0]['mu']),
                'mu_fermion_threshold': float(uv[-1]['mu']),
                'mu_scalar_threshold': float(eft[-1]['mu']),
                # No saved SM+Weinberg curve: intentionally do not invent one.
                'mu_low': float(eft[-1]['mu']),
            },
            **scalars,
        }
    return payloads, metadata


def write_existing_rge_dashboard(trajectories: Mapping[str, dict], output_path: Path) -> dict:
    """Run the original DisplayLambda renderer on a temporary scalar adapter."""
    from Numerical.plotting import InteractiveModelComparison as upstream

    payloads, quantities = scalar_curves(trajectories)
    # The existing renderer creates Weinberg and LLSS menus unconditionally.
    # Retain metadata for them to satisfy its API, then remove those menu
    # groups from the generated page because these data are not available.
    placeholders = {'c5_norm': ('unavailable', 'final_running', 'c5_abs', 'GeV⁻¹', 'frobenius')}
    for i in range(1, 4):
        for j in range(i, 4):
            placeholders[f'c5_{i}{j}'] = ('unavailable', 'final_running', 'c5_abs', 'GeV⁻¹', (i-1, j-1))
    effective = {**placeholders, **dict(sorted(quantities.items()))}
    original = upstream.QUANTITIES
    try:
        upstream.QUANTITIES = effective
        with TemporaryDirectory(prefix='t3-general-scalar-display-') as tmp:
            root = Path(tmp)
            for model, payload in payloads.items():
                path = root / _identifier(model) / 'running_diagnostics.json'
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps(payload, allow_nan=False), encoding='utf-8')
            upstream.write_dashboard(root, output_path,
                                     'T3 general scalar: saved UV / EFT1 running')
    finally:
        upstream.QUANTITIES = original
    text = output_path.read_text(encoding='utf-8')
    for group in ('Weinberg operator — complete matching history',
                  'Weinberg-like LLSS — analytic diagnostic'):
        pattern = r'<optgroup label="' + re.escape(group) + r'">.*?</optgroup>'
        text, count = re.subn(pattern, '', text, flags=re.S)
        if count != 1:
            raise RuntimeError(f'DisplayLambda HTML selector changed: {group}')
    # Remove unavailable controls rather than offering empty or misleading data.
    text = re.sub(r'<button id="zoomFinal">.*?</button>', '', text)
    text = text.replace('<span><i class="dot" style="background:#e5f5e9"></i>SM + Weinberg</span>', '')
    text = text.replace('Saved benchmark trajectories; no matching or RGE recalculation.',
                        'Saved general-scalar UV/EFT1 scalar-parameter trajectories only; '
                        'no Weinberg or neutrino trajectory was recorded. No physics recomputation.')
    text = text.replace('Other running observables', 'Saved scalar running observables')
    text = text.replace('Other running observables', 'Saved scalar running observables')
    # A durable visible disclosure for users viewing HTML without the index.
    text = text.replace('</header>', '</header><p style="padding:10px 24px;font-size:12px;color:#526074">'
                        'Only saved scalar masses, holomorphic quadratic terms and general-scalar quartics '
                        'are shown. The UV/EFT1 threshold is real; no C5 or neutrino-mass curves are synthesized.'
                        '</p>', 1)
    output_path.write_text(text, encoding='utf-8')
    return {'renderer': 'Numerical.plotting.InteractiveModelComparison.write_dashboard',
            'models': len(payloads), 'available_scalar_quantities': len(quantities),
            'scale_resolved_c5': False}


def complete_gamma_scan_models(scans: list[dict]) -> tuple[dict, list[str]]:
    """Require all old-renderer inputs rather than invent missing observables."""
    groups = defaultdict(lambda: defaultdict(list))
    for row in scans:
        groups[row['model']][row['direction']].append(row)
    accepted, missing = {}, []
    for model, directions in sorted(groups.items()):
        if set(directions) != {'T1', 'T2', 'T3'}:
            missing.append(f'{model}: missing one or more T1/T2/T3 directions')
            continue
        complete = True
        centres = []
        for direction, rows in directions.items():
            rows.sort(key=lambda r: r['gamma'])
            if len(rows) < 3 or any(any(r.get(k) is None for k in GAMMA_OBSERVABLES)
                                    for r in rows):
                complete = False
                break
            match = [r for r in rows if abs(float(r['gamma'])-1.) <= 1e-9]
            if len(match) != 1:
                complete = False
                break
            centres.append(match[0])
        if complete:
            for k in GAMMA_OBSERVABLES:
                if any(abs(float(r[k]) - float(centres[0][k])) > 1e-7 *
                       max(abs(float(r[k])), abs(float(centres[0][k])), 1e-9)
                       for r in centres[1:]):
                    complete = False
                    break
        if not complete:
            missing.append(f'{model}: incomplete five-observable scan or inconsistent gamma centres')
            continue
        accepted[model] = directions
    return accepted, missing


def write_existing_sensitivity_ovals(scans: list[dict], output_path: Path) -> dict:
    """Re-use the historical oval renderer only for *real*, complete scans."""
    from Numerical.plotting import BuildFittedModelEnvelopes as upstream

    models, missing = complete_gamma_scan_models(scans)
    if not models:
        return {'renderer': None, 'models': 0, 'skipped': missing,
                'reason': 'No complete five-observable T1/T2/T3 saved general-scalar response scans'}
    with TemporaryDirectory(prefix='t3-general-scalar-gamma-') as tmp:
        root = Path(tmp)
        summary = []
        for model, directions in models.items():
            summary.append({'model': model, 'status': 'Success', 'verified_chi2': None})
            payload = {
                'status': 'Success', 'target': {},
                'directions': {
                    direction: {'points': [
                        {'gamma': float(r['gamma']),
                         'prediction': {k: float(r[k]) for k in GAMMA_OBSERVABLES}}
                        for r in rows]}
                    for direction, rows in directions.items()
                },
            }
            path = root / model / 'sensitivity_diagnostics.json'
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(payload, allow_nan=False), encoding='utf-8')
        (root / 'comparison_summary.json').write_text(json.dumps(summary), encoding='utf-8')
        upstream.build(root, output_path)
    page = output_path.read_text(encoding='utf-8')
    page = page.replace('Individually fitted T3 model sensitivity',
                        'General-scalar T3: saved T1/T2/T3 response ovals')
    page = page.replace('Each model uses its own fitted λT3, Yukawa matrices and other fitted parameters.',
                        'Each point comes from a supplied general-scalar gamma scan; no fitted-parameter assumptions are made.')
    page = page.replace('2× empirical standard-deviation covariance of saved UV sensitivity points',
                        '2× empirical standard-deviation covariance of supplied general-scalar gamma points')
    page = page.replace('No recomputation. Source:', 'No recomputation. Source:')
    # No verified experimental target covariance is supplied by the response CSV.
    page = re.sub(r'<label><input type="checkbox" id="experimental" checked>.*?</label>', '', page)
    page = re.sub(r'<label><input type="checkbox" id="future" checked>.*?</label>', '', page)
    page = page.replace("for(const k of ['experimental','future','points','envelopes','trajectories','centres'])$(k).onchange=draw;",
                        "for(const k of ['points','envelopes','trajectories','centres'])$(k).onchange=draw;")
    page = page.replace("$('future').checked", "false")
    page = page.replace("$('experimental').checked", "false")
    page = page.replace("' (χ²='+Number(m.chi2).toPrecision(3)+')'", "' (saved gamma scan)'")
    page = page.replace('Blue ovals: Gaussian NuFIT 1σ (Δχ²=2.30) and 3σ (Δχ²=11.83) approximations, NOT official NuFIT joint contours.',
                        'No experimental covariance overlay is supplied.')
    page = page.replace('Orange ovals: JUNO projected six-year solar precision (NuFIT centre, diagonal errors); not a global future fit.',
                        'No JUNO forecast overlay is supplied.')
    page = page.replace("D.models.length+' fitted models plotted", "D.models.length+' general-scalar scanned models plotted")
    page = page.replace("const t=D.target,n=t.observable_names||[]", "const t=D.target,n=t.observable_names||[]")
    output_path.write_text(page, encoding='utf-8')
    return {'renderer': 'Numerical.plotting.BuildFittedModelEnvelopes.build',
            'models': len(models), 'skipped': missing}
