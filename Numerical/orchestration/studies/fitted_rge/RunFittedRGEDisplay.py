"""Re-evaluate the standard scale-resolved RGE display at each model's best fit.

Reuses existing fitted_candidate_config.json and comparison symbolic artifacts.
Does NOT regenerate matching/RGBeta, refit parameters, or overwrite the original
comparison diagnostics. The existing InteractiveModelComparison dashboard is used.

Usage from project root:
 python -m Numerical.orchestration.studies.fitted_rge.RunFittedRGEDisplay
 python -m Numerical.orchestration.studies.fitted_rge.RunFittedRGEDisplay --run
 python -m Numerical.orchestration.studies.fitted_rge.RunFittedRGEDisplay --run --limit 1
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import csv
import html
import tempfile
from pathlib import Path

from common.RunRecords import RunRecord
from Numerical.orchestration.support.ResumePipelineNumericalResults import (
    PROJECT_ROOT, OUTPUT_ROOT, REQUIRED, discover_cases, parse_model_key,
)
from Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison import best_candidate
from Numerical.orchestration.PipelineNumericalResults import (
    PIPELINE_NUMERICAL_CONFIG_KIND, _run_and_build_payload,
)
from Numerical.plotting.InteractiveModelComparison import write_dashboard
from Numerical.plotting.RunningResultFigures import generate_running_result_figures

SCENARIO = 'smallY_smallL'
CLASS_BY_DIMS = {(1, 3, 2): 'A', (2, 2, 1): 'B', (2, 2, 3): 'C', (3, 1, 2): 'D', (3, 3, 2): 'E'}
DEST_ROOT = PROJECT_ROOT / 'Reports' / 'output' / 'full' / 'fitted_rge'
DISPLAY_ROOT = PROJECT_ROOT / 'Reports' / 'output' / 'display'


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding='utf-8-sig'))


def signature(payload: dict, artifacts: list[Path]) -> str:
    h = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':'),
                                ensure_ascii=False).encode('utf-8'))
    for file in artifacts:
        stat = file.stat()
        h.update(str(file.resolve()).encode('utf-8'))
        h.update(str(stat.st_size).encode('ascii'))
        h.update(str(stat.st_mtime_ns).encode('ascii'))
    return h.hexdigest()


def fixed_payload(config: Path) -> dict:
    payload = load_json(config)
    if payload.get('kind') != PIPELINE_NUMERICAL_CONFIG_KIND:
        raise ValueError(f'Unexpected numerical config kind in {config}')
    # Prevent any implicit re-fitting, Sobol search or sensitivity scan.
    payload['benchmark_search'] = {**payload.get('benchmark_search', {}), 'enabled': False}
    payload['sensitivity'] = {**payload.get('sensitivity', {}), 'enabled': False}
    return payload


def cases_for_scenario(scenario: str):
    cases = {}
    for sc, model, case in discover_cases(OUTPUT_ROOT):
        if sc != scenario:
            continue
        if model in cases:
            raise RuntimeError('Ambiguous comparison artifact directories for ' + model)
        cases[model] = case
    return cases



def short_label(model: str) -> str:
    d1, d2, df, alpha = parse_model_key(model)
    cls = CLASS_BY_DIMS.get((d1, d2, df), f'{d1}-{d2}-{df}')
    return f'T3-{cls}-{"m" if alpha < 0 else "p"}{abs(alpha)}'


def ordinary_parameter_metadata(payload: dict) -> dict:
    """Preserve the fitted UV input literally. Do not interpret running masses as pole masses."""
    ordinary = payload['base_state']['ordinary']
    names = ('lambdaT3', 'MF', 'mS1Sq', 'mS2Sq', 'y1', 'y2')
    return {name: ordinary.get(name) for name in names}


def complex_scalar_string(value) -> str:
    """Format an actual scalar, never silently coerce a matrix to float."""
    if value is None:
        return '—'
    if isinstance(value, dict):
        real, imag = value.get('real'), value.get('imag', 0)
        if isinstance(real, (list, tuple)) or isinstance(imag, (list, tuple)):
            return '[matrix — see details]'
        if real is None:
            return json.dumps(value, ensure_ascii=False)
        return f'{float(real):.8g}{float(imag):+.8g}i'
    if isinstance(value, (list, tuple)):
        return '[matrix — see details]'
    return f'{float(value):.8g}' if isinstance(value, (int, float)) else str(value)


def csv_parameter_string(value) -> str:
    """Preserve matrix entries in CSV in JSON form, rather than losing them."""
    if isinstance(value, dict) and isinstance(value.get('real'), list):
        return json.dumps(value, separators=(',', ':'), ensure_ascii=False)
    if isinstance(value, (dict, list)) and (isinstance(value, list)):
        return json.dumps(value, separators=(',', ':'), ensure_ascii=False)
    return complex_scalar_string(value)


def matrix_html(value) -> str:
    if not isinstance(value, dict) or 'real' not in value:
        return '<em>Not present in fitted config</em>'
    real, imag = value['real'], value.get('imag')
    if not isinstance(real, list):
        return '<em>Invalid matrix format</em>'
    cells = []
    for i, row in enumerate(real):
        entries = []
        for j, r in enumerate(row):
            v = 0 if imag is None else imag[i][j]
            entries.append('<td>' + html.escape(f'{float(r):.7g}{float(v):+.7g}i') + '</td>')
        cells.append('<tr>' + ''.join(entries) + '</tr>')
    return '<table class="matrix">' + ''.join(cells) + '</table>'


def write_parameter_report(records: list[dict], output_dir: Path, display_dir: Path,
                           scenario: str, dashboard_name: str) -> None:
    """Write full fitted UV inputs plus a compact human-readable summary."""
    display_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    params = []
    for row in records:
        if row['status'] not in ('Success', 'AlreadyComplete'):
            continue
        config = row.get('fit_config')
        if not config:
            continue
        payload = fixed_payload(Path(config))
        params.append({
            'model': row['model'], 'label': short_label(row['model']),
            'fit_chi2': row['fit_chi2'], 'fit_stage': row['fit_stage'],
            'config': config,
            'uv_ordinary': ordinary_parameter_metadata(payload),
        })
    (output_dir / 'optimised_benchmark_parameters.json').write_text(
        json.dumps({'scenario': scenario, 'models': params}, indent=2) + '\n',
        encoding='utf-8')
    with (output_dir / 'optimised_benchmark_parameters.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['model', 'label', 'fit_chi2', 'fit_stage',
                                               'lambdaT3', 'MF', 'mS1Sq', 'mS2Sq', 'config'])
        writer.writeheader()
        for item in params:
            o = item['uv_ordinary']
            writer.writerow({k: (csv_parameter_string(o.get(k)) if k in ('lambdaT3', 'MF', 'mS1Sq', 'mS2Sq')
                                 else item.get(k, '')) for k in writer.fieldnames})
    rows = []
    for item in params:
        o = item['uv_ordinary']
        scalars = ''.join('<td>' + html.escape(complex_scalar_string(o.get(k))) + '</td>'
                          for k in ('lambdaT3', 'MF', 'mS1Sq', 'mS2Sq'))
        detail = (f'<details><summary>Show MF and complex Yukawa matrices</summary>'
                  f'<b>MF (UV, GeV)</b>{matrix_html(o.get("MF"))}'
                  f'<b>y1 (UV)</b>{matrix_html(o.get("y1"))}'
                  f'<b>y2 (UV)</b>{matrix_html(o.get("y2"))}'
                  f'<p>Saved fitted config: <code>{html.escape(str(item["config"]))}</code></p>'
                  f'</details>')
        rows.append(f'<tr><th>{html.escape(item["label"])}</th>'
                    f'<td>{item["fit_chi2"]:.5g}</td>{scalars}<td>{detail}</td></tr>')
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>Optimised T3 RGE — benchmark inputs</title><style>
body{{font:14px system-ui,Arial,sans-serif;color:#213047;max-width:1500px;margin:35px auto;padding:0 24px}}
a{{color:#2362a8}}table{{border-collapse:collapse;width:100%;font-size:12px}}
th,td{{border:1px solid #d6dde6;padding:9px;vertical-align:top;text-align:left}}
thead{{background:#edf2f8}}.matrix{{width:auto;margin:8px 0 18px;font-family:monospace}}
.matrix td{{white-space:nowrap}}details{{min-width:210px}}summary{{cursor:pointer}}
code{{white-space:normal;overflow-wrap:anywhere}}.note{{background:#f4f6f9;padding:12px}}
</style></head><body><h1>Optimised T3 — fitted UV benchmark parameters</h1>
<p><a href="{html.escape(dashboard_name, quote=True)}">Open interactive optimised RGE display →</a></p>
<p class="note">Values are the <strong>UV boundary conditions from the saved fitted candidate configs</strong>,
not their running values at another scale. mS1Sq and mS2Sq are mass-squared parameters (GeV² if the
fitted configuration uses GeV units); MF is the fermion mass matrix (GeV; click details to see entries). Their values are copied
without rescaling. The full matrices are complex and shown under each model. These are not pole masses.
{len(params)} completed model benchmark(s). Scenario: {html.escape(scenario)}.</p>
<table><thead><tr><th>Model</th><th>Fitted χ²</th><th>λT3</th><th>MF</th><th>mS1Sq</th><th>mS2Sq</th><th>Yukawas</th></tr></thead>
<tbody>{''.join(rows)}</tbody></table></body></html>"""
    (display_dir / f'optimised_rge_parameters_{scenario}.html').write_text(page, encoding='utf-8')


def add_parameters_link(dashboard: Path, scenario: str) -> None:
    """Inject a plain navigation link without changing the plotting machinery."""
    page = dashboard.read_text(encoding='utf-8')
    name = f'optimised_rge_parameters_{scenario}.html'
    link = (f'<div style="padding:9px 24px;background:#eef3fa;border-bottom:1px solid #ccd8e7">'
            f'<a href="{html.escape(name, quote=True)}">Optimised benchmark inputs: '
            f'λT3, scalar masses and full fitted Yukawa matrices</a></div>')
    if '<body' in page:
        import re
        page = re.sub(r'(<body[^>]*>)', lambda m: m.group(1) + link, page, count=1, flags=re.I)
    else:
        raise ValueError('Generated dashboard has no HTML body tag')
    dashboard.write_text(page, encoding='utf-8')

def build_display_from_completed(records: list[dict], output: Path, dashboard: Path, title: str) -> None:
    """Build the original comparison dashboard with one distinct entry per model.

    The upstream plotter expects <study>/<model>/data/running_diagnostics.json.
    Our fitted-RGE cache is <study>/<model>/running_diagnostics.json.  Use a
    temporary, presentation-only directory to satisfy that contract; never
    relocate, mutate or recalculate saved numerical diagnostics.
    """
    with tempfile.TemporaryDirectory(prefix='optimised_t3_rge_display_') as temp:
        study = Path(temp)
        included = 0
        for row in records:
            if row['status'] not in ('Success', 'AlreadyComplete'):
                continue
            original = Path(row['diagnostic_path'])
            if not original.is_file():
                raise FileNotFoundError(f'Missing completed model diagnostics: {original}')
            payload = load_json(original)
            if payload.get('status') != 'Success':
                raise RuntimeError(f'Completed model does not have successful diagnostics: {original}')
            # Display metadata only. Do not change the original saved result.
            benchmark = dict(payload.get('benchmark_search') or {})
            representation = dict(benchmark.get('model') or {})
            representation['model_key'] = short_label(row['model'])
            benchmark['model'] = representation
            payload['benchmark_search'] = benchmark
            display_file = study / row['model'] / 'data' / 'running_diagnostics.json'
            display_file.parent.mkdir(parents=True, exist_ok=True)
            display_file.write_text(json.dumps(payload, ensure_ascii=False), encoding='utf-8')
            included += 1
        if included == 0:
            raise RuntimeError('No successfully completed fitted models available for RGE display')
        write_dashboard(study, dashboard, title)
        print(f'Built display from {included} individually fitted model trajectories')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scenario', default=SCENARIO)
    p.add_argument('--run', action='store_true', help='Run missing baseline scale trajectories')
    p.add_argument('--limit', type=int, default=0, help='Maximum model cases to inspect, zero means all')
    p.add_argument('--max-fit-chi2', type=float, default=0.01)
    p.add_argument('--force', action='store_true', help='Recalculate existing fitted RGE outputs')
    p.add_argument('--no-figures', action='store_true')
    args = p.parse_args(argv)
    if args.limit < 0 or args.max_fit_chi2 < 0:
        p.error('Require nonnegative --limit and --max-fit-chi2')

    output = DEST_ROOT / args.scenario
    records = []
    case_map = cases_for_scenario(args.scenario)
    for model in sorted(case_map)[:args.limit or None]:
        case = case_map[model]
        row = dict(model=model, status='Ready', fit_stage='', fit_chi2=None,
                   diagnostic_path='', message='')
        selected = best_candidate(args.scenario, model)
        if selected is None:
            row.update(status='MissingFit', message='No saved fitted candidate')
            records.append(row); print(row['status'], model, flush=True); continue
        chi2, stage, config, summary = selected
        row.update(fit_stage=stage, fit_chi2=chi2, fit_config=str(config))
        if not math.isfinite(chi2) or chi2 > args.max_fit_chi2:
            row.update(status='FitAboveThreshold', message=f'chi2={chi2:.6g}')
            records.append(row); print(row['status'], model, flush=True); continue
        artifacts = [case / 'data' / name for name in REQUIRED]
        missing = [str(a) for a in artifacts if not a.is_file()]
        if missing:
            row.update(status='MissingArtifacts', message='; '.join(missing))
            records.append(row); print(row['status'], model, flush=True); continue
        destination = output / model / 'running_diagnostics.json'
        row['diagnostic_path'] = str(destination)
        try:
            payload = fixed_payload(config)
            key = signature(payload, artifacts)
            if destination.is_file() and not args.force:
                previous = load_json(destination)
                if previous.get('status') == 'Success' and previous.get('fitted_rge', {}).get('input_signature') == key:
                    row['status'] = 'AlreadyComplete'
                    records.append(row); print(row['status'], model, flush=True); continue
            if not args.run:
                row['status'] = 'NeedsRGE' if destination.is_file() else 'Ready'
            else:
                d1, d2, df, alpha = parse_model_key(model)
                record = RunRecord(name=case.name.replace('_', '-', 1), alpha=alpha,
                                   d_s1=d1, d_s2=d2, d_f=df, return_code=0,
                                   summary={}, output_dir=case)
                result = _run_and_build_payload(record=record, payload=payload, config_path=config)
                result['fitted_rge'] = {
                    'source_fitted_config': str(config),
                    'source_fit_summary': str(summary),
                    'fit_stage': stage,
                    'reported_fit_chi2': chi2,
                    'input_signature': key,
                    'note': 'Individually fitted UV input, standard numerical RGE diagnostics',
                }
                destination.parent.mkdir(parents=True, exist_ok=True)
                temp = destination.with_name(destination.name + '.tmp')
                temp.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
                temp.replace(destination)
                if not args.no_figures:
                    generate_running_result_figures(destination, output / model / 'figures')
                row['status'] = 'Success'
        except Exception as exc:
            row.update(status='Failed', message=f'{type(exc).__name__}: {exc}')
        records.append(row)
        print(f"{row['status']:18} {model} fit_chi2={chi2:.6g} {row['message']}", flush=True)

    if args.run:
        output.mkdir(parents=True, exist_ok=True)
        (output / 'fitted_rge_summary.json').write_text(json.dumps(records, indent=2)+'\n', encoding='utf-8')
        complete = sum(r['status'] in ('Success', 'AlreadyComplete') for r in records)
        if complete:
            dashboard = DISPLAY_ROOT / f'optimised_rge_{args.scenario}.html'
            build_display_from_completed(records, output, dashboard, f'Optimised T3 RGE: {args.scenario}')
            write_parameter_report(records, output, DISPLAY_ROOT, args.scenario, dashboard.name)
            add_parameters_link(dashboard, args.scenario)
            print('Interactive optimised RGE display:', dashboard)
            print('Optimised benchmark inputs:', DISPLAY_ROOT / f'optimised_rge_parameters_{args.scenario}.html')
        else:
            print('No completed fitted RGE diagnostics to display.')
    else:
        print('Dry run: add --run to calculate missing scale-resolved trajectories.')
    return 1 if any(r['status'] == 'Failed' for r in records) else 0


if __name__ == '__main__':
    raise SystemExit(main())
