"""Shared fitted A(-2) Yukawa study through existing general T3 evaluators.

For each selected model keep its own representation, SM/heavy masses, scalar
couplings, general RGEs and matching source. Copy only the A(-2) fitted Y1,Y2
UV matrices, and seed symmetry-allowed conjugate Yukawas from the A(-2)
reference matrices (matching A(-2)'s original full-16 fit convention).

Generate actual T1/T2/T3 numerical scans and feed the original
Numerical.plotting.BuildFittedModelEnvelopes interactive oval renderer.
No restricted C5 fallback, no invented RGE trajectories, no refitting.
Mixed B_-1/C_-1 are explicitly labelled conditional mixed-one-scale matching.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import csv
from hashlib import sha256
import html
import json
import math
from pathlib import Path
import re
import sys
import traceback

import numpy as np
from scripts import RunT3GeneralAll16FromReference as campaign
from scripts import RunT3Aminus2YukawaBenchmark as a2benchmark

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'output/studies/general_scalar_extended/shared_a2_general16_reports'
FORMAT = 't3-shared-a2-general16-t123-conditional-v1'
DIRECTIONS = ('T1', 'T2', 'T3')
DEFINITIONS = {
    'T1': 'lambdaT3 -> gamma lambdaT3; y1,y2 -> gamma^(-1/2) y1,y2',
    'T2': 'lambdaT3 -> gamma lambdaT3; y2 -> gamma^(-1) y2',
    'T3': 'y1 -> gamma y1; y2 -> gamma^(-1) y2',
}
OBS = campaign.OBS


def require(test: bool, message: str) -> None:
    if not test:
        raise ValueError('Shared A(-2) general16: ' + message)


def valid_id(value: str) -> str:
    require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', value)), 'unsafe run-id')
    return value


def gammas(lo: float, hi: float, n: int) -> list[float]:
    require(n >= 3 and n % 2 == 1 and n <= 101,
            'gamma points must be odd and in [3,101]')
    require(math.isfinite(lo) and math.isfinite(hi) and 0 < lo < 1 < hi,
            'require 0 < gamma-min < 1 < gamma-max')
    half = n // 2
    g = np.r_[np.geomspace(lo, 1.0, half + 1),
              np.geomspace(1.0, hi, half + 1)[1:]]
    g[half] = 1.0
    return [float(v) for v in g]


def perturb(ordinary: dict, direction: str, gamma: float) -> tuple[dict, dict]:
    """Same three UV directions as Numerical.fitting.BenchmarkSensitivity.

    Return new (18-real Yukawa parameter dictionary, lambdaT3 scalar spec).
    """
    require(direction in DIRECTIONS and math.isfinite(gamma) and gamma > 0,
            'invalid sensitivity deformation')
    result = dict(ordinary)
    lam_factor = gamma if direction in ('T1', 'T2') else 1.0
    if direction == 'T1':
        factors = {'y1': gamma**(-0.5), 'y2': gamma**(-0.5)}
    elif direction == 'T2':
        factors = {'y1': 1.0, 'y2': gamma**(-1)}
    else:
        factors = {'y1': gamma, 'y2': gamma**(-1)}
    for name in ('y1', 'y2'):
        for i in range(3):
            for j in range(3):
                k = f'{name}_{i}{j}'
                result[k] = float(ordinary[k]) * factors[name]
    return result, {'real': lam_factor, 'imag': 0.0}


def reference_with_lambda(reference: dict, scale: float) -> dict:
    cfg = deepcopy(reference)
    old = cfg['base_state']['ordinary']['lambdaT3']
    require(abs(float(old.get('imag', 0.0))) < 1e-13,
            'only real lambdaT3 used by this exploratory scan')
    old['real'] = float(old['real']) * scale
    return cfg


def fixed_for_model(model: str, reference: dict, uv: dict, shared_reference: dict,
                    settings: dict, *, lambda_scale: float, conj_source: str) -> dict:
    varied = reference_with_lambda(reference, lambda_scale)
    fixed = campaign.fixed_overrides(
        varied, uv, quartic_fraction=float(settings['quartic_fraction']),
        mu2_fraction=float(settings['mu2_fraction']),
        extra_yukawa_fraction=float(settings['extra_yukawa_fraction']))
    if conj_source == 'a2-reference':
        for suffix, original in (('y1Conj', 'y1'), ('y2Conj', 'y2')):
            if suffix in fixed:
                y = campaign.matrix(shared_reference['base_state']['ordinary'][original])
                fixed[suffix] = campaign.entry_matrix(
                    y * float(settings['extra_yukawa_fraction']))
    else:
        require(conj_source == 'model-reference', 'unknown extra Yukawa source')
    return fixed


def target_payload(target) -> dict:
    return {'ordering': str(target.ordering),
            'observable_names': list(target.observable_names),
            'central_values': np.asarray(target.central_values, dtype=float).tolist(),
            'covariance': np.asarray(target.covariance, dtype=float).tolist(),
            'source': str(target.source),
            'note': 'Gaussian approximation; overlay contours are not official confidence regions.'}


def donor(run_id: str) -> tuple[dict, dict, dict, dict]:
    """Verify fitted A(-2) donor and return (parameters, settings, provenance, reference)."""
    reference = campaign.read(campaign.ref_path('A_-2'))
    params, settings, provenance = a2benchmark.load_fit(run_id, reference)
    require(set(params) == set(campaign.NAMES), 'donor Yukawa parameter names differ')
    require(all(math.isfinite(x) for x in params.values()), 'nonfinite donor Yukawas')
    return params, settings, provenance, reference


def load_evaluator(model: str, settings: dict, shared_reference: dict,
                   target, conj_source: str):
    """Return model-specific evaluator with a function to build each gamma point.

    The mixed models use the exact conditional B/C evaluator from the 16-run;
    14 others use the full-general Matchete channel evaluator.
    """
    reference = campaign.read(campaign.ref_path(model))
    if model not in campaign.MIXED:
        uv, eft, provenance, kernel, source, kind = campaign.nonmixed_inputs(
            model, str(settings['canonical_run_id']), generate_general_c5=False)
        status = campaign.read(campaign.general_coefficient_status(model))
        require(status.get('status') == 'Matched' and status.get('model') == model,
                model + ': invalid general Matchete matching status')
        if model in campaign.SELF:
            require({'T31PlusMass', 'T31MinusMass'} <= kernel.symbols or
                    {'T32PlusMass', 'T32MinusMass'} <= kernel.symbols,
                    model + ': self-bilinear mass splitting missing from C5')
        require(source.is_file(), model + ': missing general C5')

        def evaluate(params, lam_scale):
            fixed = fixed_for_model(model, reference, uv, shared_reference, settings,
                                    lambda_scale=lam_scale, conj_source=conj_source)
            engine = campaign.NonMixedEvaluator(model, uv, eft, kernel, fixed, target)
            return engine.evaluate(params, detail=True)

        matching = 'full-general Matchete one-scale Yukawa channels' + (
            ' / split self-bilinear masses' if model in campaign.SELF else '')
        source_file = source
    else:
        from scripts import RunT3BCGeneralFromReference as bc
        from scripts import RunT3GeneralScalarBCExploratoryFit as bcfit
        from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE
        source = bc.build_reference_source(model)
        bc.ensure_betas(model, source, export=False, wolframscript='wolframscript')
        coefficient = DEFAULT_SOURCE / (
            {'B_-1': 'B_m1_final', 'C_-1': 'C_m1_final'}[model]
            + '_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
        compiled = bcfit.compile_source(coefficient)
        bundle = bcfit.beta_bundle(model, source, source['registry'])

        def evaluate(params, lam_scale):
            point_source = deepcopy(source)
            entry = point_source['cfg']['base_state']['ordinary']['lambdaT3']
            entry['real'] = float(entry['real']) * lam_scale
            entry['imag'] = float(entry['imag']) * lam_scale
            point_source['lambdaT3'] = complex(point_source['lambdaT3']) * lam_scale
            engine = bcfit.CandidateEvaluator(
                'reference_scalar_general', model, point_source,
                point_source['registry'], bundle, compiled, target,
                quartic_fraction=float(settings['quartic_fraction']),
                mu2_fraction=float(settings['mu2_fraction']))
            return engine.evaluate(params, detail=True)

        provenance = {'scope': 'conditional mixed scalar basis, not full-general matching'}
        matching = 'conditional B/C mixed-scalar one-scale (not full-general Matchete)'
        source_file = coefficient
    return evaluate, {'model': model, 'reference': str(campaign.ref_path(model)),
                      'matching_scheme': matching, 'matching_source': str(source_file),
                      'matching_sha256': sha256(source_file.read_bytes()).hexdigest(),
                      'symbolic_provenance': provenance}


def encode_prediction(fit) -> dict:
    prediction = dict(zip(fit.observable_names, [float(x) for x in fit.prediction]))
    require(all(k in prediction and math.isfinite(prediction[k]) for k in OBS),
            'incomplete or nonfinite neutrino observables')
    score = float(fit.chi2)
    require(math.isfinite(score) and score >= 0, 'invalid chi2')
    return {'status': 'Success', 'chi2': score, 'prediction': prediction}


def run_one(model: str, params: dict, settings: dict, shared_reference: dict,
            target, grid: list[float], conj_source: str, dest: Path) -> dict:
    evaluate, physics_source = load_evaluator(
        model, settings, shared_reference, target, conj_source)
    directions = {}
    centres = {}
    for direction in DIRECTIONS:
        points = []
        for gamma in grid:
            transformed, scale = perturb(params, direction, gamma)
            point = {'gamma': gamma}
            try:
                fit, detail = evaluate(transformed, scale['real'])
                point.update(encode_prediction(fit))
                if gamma == 1.0:
                    centres[direction] = (point, detail)
            except Exception as exc:
                point.update(status='Failed', message=f'{type(exc).__name__}: {exc}')
            points.append(point)
        directions[direction] = {'definition': DEFINITIONS[direction], 'points': points}
    good = all(len([p for p in directions[d]['points'] if p['status']=='Success']) >= 3
               and d in centres for d in DIRECTIONS)
    if good:
        zero = centres['T1'][0]['prediction']
        for d in ('T2', 'T3'):
            for k in OBS:
                if not math.isclose(zero[k], centres[d][0]['prediction'][k],
                                    rel_tol=2e-6, abs_tol=1e-11):
                    raise ValueError(model + ': inconsistent gamma=1 centres')
    status = 'Success' if good else 'IncompleteScan'
    scan = {'status': status, 'model': model,
            'gamma_min': grid[0], 'gamma_max': grid[-1], 'n_points': len(grid),
            'directions': directions, 'target': target_payload(target),
            'n_failed_points': sum(p['status']=='Failed' for d in DIRECTIONS
                                   for p in directions[d]['points']),
            'scheme': physics_source['matching_scheme'],
            'note': 'All models share A(-2) fitted Y1/Y2. These are model-specific conditional calculations; not individually refitted.'}
    dest.mkdir(parents=True, exist_ok=True)
    campaign.save(dest/'sensitivity_diagnostics.json', scan)
    centre_data = {'model': model, 'status': status, 'physics_source': physics_source,
                   'input_yukawas': {'y1': campaign.matrix_json(
                        np.array([[params[f'y1_{i}{j}'] for j in range(3)] for i in range(3)])),
                        'y2': campaign.matrix_json(
                        np.array([[params[f'y2_{i}{j}'] for j in range(3)] for i in range(3)]))},
                   'chi2': centres.get('T1', ({},))[0].get('chi2') if 'T1' in centres else None,
                   'prediction': centres.get('T1', ({},))[0].get('prediction') if 'T1' in centres else None,
                   'numerical_details': centres['T1'][1] if 'T1' in centres else None}
    campaign.save(dest/'centre.json', a2benchmark.yaml_free_json(centre_data))
    return {'model': model, 'status': status,
            'verified_chi2': centre_data['chi2'],
            'fit_stage': 'shared_A_m2_fitted_Yukawa_forward_evaluation',
            'matching_scheme': physics_source['matching_scheme'],
            'prediction': centre_data['prediction'],
            'failed_points': scan['n_failed_points']}


def report_html(rows: list[dict], target_dir: Path, *, donor_run: str, gamma_grid: list[float], oval: bool, rge: bool = False) -> None:
    """Presentation-only overview. Ovals produced by existing renderer separately."""
    out = target_dir/'index.html'
    head = ('<!doctype html><html lang="en"><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            '<title>Shared A(-2) general T3 comparison</title>'
            '<style>body{font:15px system-ui;max-width:1150px;margin:36px auto;'
            'padding:0 20px;color:#17243b}table{border-collapse:collapse;width:100%}'
            'th,td{padding:8px;border-bottom:1px solid #dfe4ed;text-align:left}'
            'th{background:#f4f7fc}a{color:#215aa4}tr:hover{background:#f9fbff}'
            'p{line-height:1.6}</style>')
    intro = (f'<h1>Shared A(-2) Yukawas: 16 T3 models</h1><p>Donor: '
             f'<code>{html.escape(donor_run)}</code>. All models use the same fitted '
             f'UV Y1/Y2, but keep their own mass parameters, scalar couplings, RGEs, '
             f'and matching. No re-optimisation. Gamma range: {gamma_grid[0]:g}–{gamma_grid[-1]:g}.</p>'
             '<p>' + ('<a href="sensitivity_oval.html">Interactive sensitivity ovals</a> · ' if oval else '') +
             ('<a href="interactive_RGE.html">Interactive UV/EFT1 RGEs</a> · ' if rge else '') +
             '<a href="comparison.csv">Observable comparison CSV</a> · '
             '<a href="gamma_points.csv">Gamma scan points</a> · '
             '<a href="manifest.json">Physics provenance</a></p>'
             '<p>B(-1) and C(-1) use conditional mixed-scalar one-scale matching. '
             'The other models use full-general one-scale Matchete matching. '
             'Sequential EFT matching and the physical coupling map remain uncertified. '
             'The optional interactive RGE report shows genuine newly saved UV/EFT1 solver points at gamma=1 for non-mixed models only; it does not invent full SMEFT C5 running.</p>')
    lines = [head, intro, '<table><thead><tr><th>Model</th><th>Status</th>'
             '<th>χ² at γ=1</th><th>Δm²₂₁ [eV²]</th><th>Δm²₃ℓ [eV²]</th>'
             '<th>sin²θ₁₂</th><th>sin²θ₁₃</th><th>sin²θ₂₃</th><th>Source</th></tr></thead><tbody>']
    for r in rows:
        p = r.get('prediction') or {}
        def cell(v):
            return html.escape(f'{v:.7g}' if isinstance(v, (int,float)) and math.isfinite(v) else '—')
        lines.append('<tr>' + ''.join(f'<td>{html.escape(str(x))}</td>' for x in
             (r['model'], r['status'])) +
             '<td>'+cell(r.get('verified_chi2'))+'</td>' +
             ''.join('<td>'+cell(p.get(k))+'</td>' for k in OBS) +
             '<td>'+html.escape(r.get('matching_scheme', ''))+'</td></tr>')
    lines.append('</tbody></table></html>')
    out.write_text(''.join(lines), encoding='utf-8')


def saved_rows(path: Path) -> list[dict]:
    """Read list-valued comparison summary (campaign.read accepts objects only)."""
    if not path.is_file():
        raise FileNotFoundError(str(path))
    rows = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(rows, list) and all(isinstance(row, dict) for row in rows),
            'invalid comparison summary: ' + str(path))
    return rows


def build_oval(source: Path) -> bool:
    """Adapt saved tag-based model folders to the existing oval renderer.

    Numerical calculations are not repeated. The renderer expects folders named
    by the literal model keys (A_-2), whereas numerical results use A_m2.
    Stage verified copies; never relabel a physical result or drop failures.
    """
    from Numerical.plotting.BuildFittedModelEnvelopes import build

    rows = saved_rows(source/'comparison_summary.json')
    require(isinstance(rows, list), 'invalid comparison summary for oval rendering')
    display = source/'_sensitivity_display_inputs'
    display.mkdir(parents=True, exist_ok=True)
    for row in rows:
        if row.get('status') != 'Success':
            continue
        model = str(row['model'])
        require(model in campaign.MODELS, 'unknown model in oval summary: ' + model)
        original = source/campaign.tag(model)/'sensitivity_diagnostics.json'
        require(original.is_file(), model + ': missing saved sensitivity scan')
        scan = campaign.read(original)
        require(scan.get('status') == 'Success' and scan.get('model') == model,
                model + ': invalid or mismatched saved sensitivity scan')
        folder = display/model
        folder.mkdir(parents=True, exist_ok=True)
        (folder/'sensitivity_diagnostics.json').write_bytes(original.read_bytes())
    campaign.save(display/'comparison_summary.json', rows)

    try:
        data = build(display, source/'sensitivity_oval.html')
    except RuntimeError as exc:
        if 'No verified successful fitted sensitivity results' not in str(exc):
            raise
        return False
    require(bool(data['models']), 'renderer produced no models')
    page = source/'sensitivity_oval.html'
    content = page.read_text(encoding='utf-8')
    content = content.replace('Individually fitted T3 model sensitivity',
                              'Shared A(-2) UV Yukawa T3 sensitivity')
    content = content.replace(
        'Each model uses its own fitted λT3, Yukawa matrices and other fitted parameters.',
        'All models use the same A(-2) fitted UV Y1/Y2; other model inputs are model-specific. Conditional one-scale matching, not a refit.')
    content = content.replace('Fitted centres', 'Shared-Yukawa centres')
    page.write_text(content, encoding='utf-8')
    return True


def render_saved_only(run_id: str) -> int:
    """Recreate existing HTML from saved results without any Matchete/RGE calls."""
    out = OUTPUT/valid_id(run_id)
    settings = campaign.read(out/'settings.json')
    rows = saved_rows(out/'comparison_summary.json')
    require(settings.get('format') == FORMAT and isinstance(rows, list),
            'invalid existing shared-A2 campaign')
    require(settings.get('run_id', run_id) == run_id, 'run-id mismatch')
    oval = build_oval(out)
    from scripts import RunT3SharedA2General16RGE as rge
    rge_ok = rge.write_html(out, rows)
    report_html(rows, out, donor_run=settings['donor_run_id'],
                gamma_grid=settings['gammas'], oval=oval, rge=rge_ok)
    manifest_path = out/'manifest.json'
    manifest = campaign.read(manifest_path) if manifest_path.is_file() else {}
    manifest['interactive_oval_generated'] = oval
    manifest['interactive_rge_generated'] = rge_ok
    manifest['successful_scan_count'] = sum(r.get('status') == 'Success' for r in rows)
    campaign.save(manifest_path, manifest)
    print('HTML:', out/'index.html', flush=True)
    print('Sensitivity ovals:', out/'sensitivity_oval.html' if oval else 'None available', flush=True)
    print('RGE HTML:', out/'interactive_RGE.html' if rge_ok else 'No sampled RGE curves saved', flush=True)
    return 0 if oval else 2


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--donor-run-id', default='thesis16_general_c5_full_v1')
    parser.add_argument('--run-id', required=True)
    parser.add_argument('--models', help='Comma-separated models; default all 16')
    parser.add_argument('--gamma-min', type=float, default=0.8)
    parser.add_argument('--gamma-max', type=float, default=1.2)
    parser.add_argument('--gamma-points', type=int, default=21)
    parser.add_argument('--conj-source', choices=('a2-reference','model-reference'),
                        default='a2-reference')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--render-only', action='store_true',
                        help='Rebuild HTML from existing saved scans; never run RGEs or matching')
    parser.add_argument('--no-render', action='store_true')
    parser.add_argument('--rge-points', type=int, default=41,
                        help='Save actual UV/EFT1 solver samples per stage at gamma=1')
    parser.add_argument('--no-rge', action='store_true',
                        help='Skip generating new dense RGE trajectories')
    parser.add_argument('--rge-only', action='store_true',
                        help='Add/recreate genuine RGE HTML for saved scan, no gamma scans')
    args = parser.parse_args(argv)
    if args.render_only:
        try:
            return render_saved_only(args.run_id)
        except (ValueError, FileNotFoundError, KeyError) as exc:
            print(f'ERROR: {exc}', file=sys.stderr)
            return 2
    try:
        if args.rge_only:
            out = OUTPUT/valid_id(args.run_id)
            saved = campaign.read(out/'settings.json')
            require(saved.get('format') == FORMAT, 'not a shared-A2 report campaign')
            require(args.donor_run_id == saved['donor_run_id'],
                    '--donor-run-id must match the saved campaign')
            params, settings, provenance, shared_ref = donor(saved['donor_run_id'])
            require(provenance['fit_result_sha256'] == saved['donor_fit_result_sha256']
                    and provenance['fit_config_sha256'] == saved['donor_fit_config_sha256'],
                    'fitted Yukawa source changed since sensitivity campaign')
            rows = saved_rows(out/'comparison_summary.json')
            from scripts import RunT3SharedA2General16RGE as rge
            generated = rge.generate(out, rows, params, settings, shared_ref,
                                     saved['conj_source'], n_points=args.rge_points)
            render_saved_only(args.run_id)
            return 0 if generated else 2
        valid_id(args.run_id)
        valid_id(args.donor_run_id)
        grid = gammas(args.gamma_min, args.gamma_max, args.gamma_points)
        selected = [x.strip() for x in args.models.split(',') if x.strip()] if args.models else list(campaign.MODELS)
        require(bool(selected) and len(selected) == len(set(selected))
                and not (set(selected) - set(campaign.MODELS)), 'invalid model list')
        selected = [m for m in campaign.MODELS if m in selected]
        params, settings, provenance, shared_ref = donor(args.donor_run_id)
        from Numerical.fitting.OscillationFit import OscillationFitTarget
        from Numerical.orchestration.RunSmallLambdaRefit import TARGET
        target = OscillationFitTarget.from_json(TARGET)
        out = OUTPUT/args.run_id
        signature = {'format': FORMAT, 'donor_run_id': args.donor_run_id,
                     'donor_fit_result_sha256': provenance['fit_result_sha256'],
                     'donor_fit_config_sha256': provenance['fit_config_sha256'],
                     'models': selected, 'gammas': grid, 'conj_source': args.conj_source,
                     'quartic_fraction': settings['quartic_fraction'],
                     'mu2_fraction': settings['mu2_fraction'],
                     'extra_yukawa_fraction': settings['extra_yukawa_fraction']}
        if (out/'settings.json').is_file():
            require(args.resume and campaign.read(out/'settings.json') == signature,
                    'existing output; use --resume with unchanged inputs or a fresh --run-id')
        else:
            require(not out.exists() or not any(out.iterdir()),
                    'output directory exists and is nonempty')
            out.mkdir(parents=True, exist_ok=True)
            campaign.save(out/'settings.json', signature)
        statuses = []
        for idx, model in enumerate(selected, 1):
            model_dir = out/campaign.tag(model)
            result = model_dir/'centre.json'
            scan = model_dir/'sensitivity_diagnostics.json'
            if args.resume and result.is_file() and scan.is_file():
                old = campaign.read(result)
                if old.get('status') == 'Success' and old.get('model') == model:
                    statuses.append({'model': model, 'status': 'Success',
                        'verified_chi2': old['chi2'],
                        'fit_stage': 'shared_A_m2_fitted_Yukawa_forward_evaluation',
                        'matching_scheme': old['physics_source']['matching_scheme'],
                        'prediction': old['prediction']})
                    print(f'[{idx}/{len(selected)}] {model}: reuse chi2={old["chi2"]:.7g}', flush=True)
                    continue
            print(f'[{idx}/{len(selected)}] {model}: shared A(-2) sensitivity', flush=True)
            try:
                row = run_one(model, params, settings, shared_ref, target,
                              grid, args.conj_source, model_dir)
                if model == 'A_-2' and args.conj_source == 'a2-reference' and row['status'] == 'Success':
                    expected = float(provenance['fit_chi2_original'])
                    actual = float(row['verified_chi2'])
                    require(math.isclose(actual, expected, rel_tol=5e-5, abs_tol=5e-4),
                        'A(-2) gamma=1 chi2 does not reproduce donor full16 fit: '
                        f'actual {actual:g}, saved {expected:g}')
                statuses.append(row)
                print(f'  {row["status"]}: chi2={row.get("verified_chi2")}', flush=True)
            except Exception as exc:
                statuses.append({'model': model, 'status': 'Blocked',
                                 'message': f'{type(exc).__name__}: {exc}'})
                model_dir.mkdir(parents=True, exist_ok=True)
                campaign.save(model_dir/'error.json', {'model': model,
                    'error': f'{type(exc).__name__}: {exc}',
                    'traceback': traceback.format_exc(limit=9)})
                print(f'  BLOCKED: {type(exc).__name__}: {exc}', flush=True)
            campaign.save(out/'comparison_summary.json', statuses)
        campaign.save(out/'comparison_summary.json', statuses)
        with (out/'comparison.csv').open('w', newline='', encoding='utf-8') as stream:
            columns = ['model','status','verified_chi2',*OBS,'matching_scheme','message']
            writer = csv.DictWriter(stream, fieldnames=columns)
            writer.writeheader()
            for row in statuses:
                writer.writerow({**{k: row.get(k, '') for k in columns},
                                 **{k: (row.get('prediction') or {}).get(k, '') for k in OBS}})
        with (out/'gamma_points.csv').open('w', newline='', encoding='utf-8') as stream:
            columns = ['model','direction','gamma','status','chi2',*OBS,'message']
            writer = csv.DictWriter(stream, fieldnames=columns)
            writer.writeheader()
            for model in selected:
                source = out/campaign.tag(model)/'sensitivity_diagnostics.json'
                if not source.is_file():
                    continue
                scan = campaign.read(source)
                for direction in DIRECTIONS:
                    for point in (scan.get('directions',{}).get(direction,{})
                                  .get('points',[])):
                        writer.writerow({'model':model,'direction':direction,
                            **{k:point.get(k,'') for k in ('gamma','status','chi2','message')},
                            **{k:(point.get('prediction') or {}).get(k,'') for k in OBS}})
        oval = build_oval(out) if not args.no_render else False
        from scripts import RunT3SharedA2General16RGE as rge
        if not args.no_rge:
            rge_ok = rge.generate(out, statuses, params, settings, shared_ref,
                                  args.conj_source, n_points=args.rge_points)
        else:
            rge_ok = rge.write_html(out, statuses)
        report_html(statuses, out, donor_run=args.donor_run_id, gamma_grid=grid,
                    oval=oval, rge=rge_ok)
        campaign.save(out/'manifest.json', {
            'format': FORMAT, 'settings': signature, 'donor': provenance,
            'total_models': len(selected),
            'successful_scan_count': sum(r['status'] == 'Success' for r in statuses),
            'interactive_oval_generated': oval,
            'interactive_rge_generated': rge_ok,
            'physics_certified_sequential_matching': False,
            'c5_matching': '14 model-specific general Matchete one-scale; B/C conditional mixed one-scale',
            'shared_boundary': 'Only Y1,Y2 fitted A(-2) common across all models',
            'extra_yukawas': 'Allowed conjugate Yukawas at fraction of '+args.conj_source+' parent',
            'scans': DEFINITIONS,
            'limitations': [
                'No refitting; chi2 values are forward predictions for common Yukawas.',
                'Fitted Y1/Y2 are real, additional conjugate Yukawas fixed and non-optimised.',
                'Original 16 run has conditional one-scale matching, not certified sequential matching.',
                'B/C use a different conditional mixed-scalar matching scheme.',
                'Gamma sensitivity points save UV/EFT endpoints; gamma=1 RGE HTML saves separate actual solver samples.',
                'B/C conditional solver does not supply comparable dense curves; no invented interpolation.',
                'No full-flavour scale-dependent C5 below the scalar threshold is plotted.',
                'Gamma ovals are descriptive covariances, not confidence regions.',
                'T1/T2 also rescale additional quartics initially defined as fractions of lambdaT3.',
            ]})
        print('Saved:', out/'comparison.csv', flush=True)
        print('HTML:', out/'index.html', flush=True)
        if oval:
            print('Sensitivity ovals:', out/'sensitivity_oval.html', flush=True)
        if rge_ok:
            print('Interactive RGEs:', out/'interactive_RGE.html', flush=True)
        return 0 if all(r['status'] == 'Success' for r in statuses) else 2
    except (ValueError, FileNotFoundError, KeyError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
