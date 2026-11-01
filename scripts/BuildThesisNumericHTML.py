"""Build thesis HTML reports from saved numerical diagnostics, using EXISTING plotters.

No RGE integration, fitting, or observables calculation is performed here.

    python -X utf8 -m scripts.BuildThesisNumericHTML --study thesis_A_m2
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
       'sin2_theta12', 'sin2_theta13', 'sin2_theta23')


def model_key(payload: dict) -> str:
    cfg = payload.get('numerical_config')
    if cfg and Path(cfg).is_file():
        src = json.loads(Path(cfg).read_text(encoding='utf-8-sig'))
        rep = src.get('representation', {})
        try:
            a = int(rep['alpha']); d1=int(rep['d_s1']); d2=int(rep['d_s2']); df=int(rep['d_f'])
            return (f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_'
                    f'{"m" if a<0 else "p"}{abs(a)}')
        except (KeyError, ValueError, TypeError):
            pass
    search = payload.get('benchmark_search', {})
    key = search.get('model', {}).get('model_key') if isinstance(search.get('model'), dict) else None
    if key:
        return str(key)
    raise ValueError('Saved diagnostics lack a usable model identity/config')


def observable_row(payload: dict, key: str) -> dict:
    if payload.get('status') != 'Success':
        raise ValueError(f'{key}: numerical output status is not Success')
    low = payload.get('low_energy', {})
    series = payload.get('final_running', {})
    row = {'model': key, 'ordering': payload.get('ordering', 'NO'),
           'numerical_config': payload.get('numerical_config', '')}
    for n in OBS:
        v = low.get(n)
        if v is None:
            vals = series.get(n, [])
            if not vals:
                raise ValueError(f'{key}: missing observable {n}')
            v = vals[-1]
        v = float(v)
        if not math.isfinite(v):
            raise ValueError(f'{key}: nonfinite {n}')
        row[n] = v
    masses = low.get('masses_ev', [])
    for i in range(3):
        row[f'm{i+1}_ev'] = float(masses[i]) if len(masses)>i else ''
    row['takagi_residual'] = low.get('takagi_residual', '')
    sensitivity = payload.get('sensitivity', {})
    if sensitivity.get('status') == 'Success':
        pts = sensitivity.get('directions', {}).get('T1', {}).get('points', [])
        centre = [p for p in pts if abs(float(p.get('gamma', -1))-1.) < 1e-9
                  and p.get('status') == 'Success']
        if len(centre) == 1:
            row['chi2_at_gamma1'] = centre[0].get('chi2', '')
    return row


def write_outputs(input_dir: Path, report_dir: Path) -> tuple[list[dict], Path | None]:
    paths = sorted(input_dir.rglob('running_diagnostics.json'))
    if not paths:
        raise FileNotFoundError(f'No numerical running_diagnostics.json found under {input_dir}')
    rows = []
    scans = []
    for path in paths:
        payload = json.loads(path.read_text(encoding='utf-8-sig'))
        key = model_key(payload)
        rows.append(observable_row(payload, key))
        scan = payload.get('sensitivity', {})
        if isinstance(scan, dict) and scan.get('status') == 'Success':
            scans.append((key, scan))
    report_dir.mkdir(parents=True, exist_ok=True)
    with (report_dir/'observables.csv').open('w', newline='', encoding='utf-8') as fp:
        names = list(dict.fromkeys(k for row in rows for k in row))
        writer = csv.DictWriter(fp, fieldnames=names)
        writer.writeheader(); writer.writerows(rows)
    (report_dir/'observables.json').write_text(json.dumps(rows, indent=2)+'\n', encoding='utf-8')
    # Existing renderer: no physics duplicated.
    from Numerical.plotting.InteractiveModelComparison import write_dashboard
    rge_path = report_dir/'interactive_RGE.html'
    write_dashboard(input_dir, rge_path, 'T3 thesis numerical RGE and observables')
    sensitivity_path = None
    if scans:
        source = report_dir / '_sensitivity_display_inputs'
        for key, scan in scans:
            folder = source / key
            folder.mkdir(parents=True, exist_ok=True)
            (folder/'sensitivity_diagnostics.json').write_text(
                json.dumps(scan, indent=2)+'\n', encoding='utf-8')
        statuses = [{'model': key, 'status': 'Success', 'fit_stage': 'saved_reference',
                     'verified_chi2': next((r.get('chi2_at_gamma1') for r in rows if r['model']==key), None)}
                    for key, scan in scans]
        (source/'comparison_summary.json').write_text(json.dumps(statuses, indent=2)+'\n', encoding='utf-8')
        from Numerical.plotting.BuildFittedModelEnvelopes import build
        sensitivity_path = report_dir/'sensitivity_oval.html'
        build(source, sensitivity_path)
    index = ['# T3 thesis numerical pilot', '',
             'Uses the repository original numerical evaluator and HTML renderers.',
             'Benchmark observables used in the fit are calibration targets, not independent predictions.',
             '', '- [Interactive RGEs](interactive_RGE.html)',
             '- [Neutrino observables](observables.csv)']
    if sensitivity_path is not None:
        index.append('- [Interactive sensitivity ovals](sensitivity_oval.html)')
    else:
        index.append('- Sensitivity ovals unavailable: no saved successful sensitivity scan')
    index.append('')
    (report_dir/'README.md').write_text('\n'.join(index), encoding='utf-8')
    return rows, sensitivity_path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', default='thesis_A_m2')
    parser.add_argument('--input-dir', type=Path)
    parser.add_argument('--report-dir', type=Path)
    args = parser.parse_args(argv)
    inp = args.input_dir or ROOT/'output'/args.study
    rep = args.report_dir or ROOT/'Reports'/'output'/args.study/'Numerical'
    rows, oval = write_outputs(inp, rep)
    print(f'Numerical results: {len(rows)} model(s)')
    print('RGE/observables HTML:', rep/'interactive_RGE.html')
    print('Observables CSV:', rep/'observables.csv')
    if oval:
        print('Sensitivity oval HTML:', oval)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
