"""Build lambda-style HTML reports for frozen general-scalar T3 fits.

Presentation only. Reads saved state.json/fit.json files and optional saved
RGE trajectory / T1-T3 sensitivity CSV artifacts if they already exist.

No physics is rerun here.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
OUT_BASE = Path('output/studies/general_scalar_extended/bulk_48')
REPORT_BASE = Path('Reports/output/display_frozen_general_scalar')
CASES = ('small_uv_frozen_5e-9', 'large_uv_frozen_0p5')
MODELS = (
    'A_+0', 'A_-2', 'A_-4', 'B_+1', 'B_-1', 'B_-3',
    'C_+1', 'C_-1', 'C_-3', 'D_+0', 'D_+2', 'D_-2',
    'E_+0', 'E_+2', 'E_-2', 'E_-4',
)
LABELS = {
    'delta_m21_sq_ev2': 'Δm²₂₁ [eV²]',
    'delta_m3l_sq_ev2': 'Δm²₃ℓ [eV²]',
    'sin2_theta12': 'sin²θ₁₂',
    'sin2_theta13': 'sin²θ₁₃',
    'sin2_theta23': 'sin²θ₂₃',
    'mass_hi_ev': 'Largest mass [eV]',
    'mass_mid_ev': 'Middle mass [eV]',
    'mass_lo_ev': 'Smallest mass [eV]',
}
OBS = (
    'mass_hi_ev', 'mass_mid_ev', 'mass_lo_ev', 'delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
    'sin2_theta12', 'sin2_theta13', 'sin2_theta23',
)
ACCEPTED = {'Accepted', 'CandidateConstrainedAccepted', 'ExistingConditionalAccepted'}
RGE_TITLE = 'Frozen general-scalar running'
SENS_TITLE = 'Frozen-model T1/T2/T3 sensitivity'


def require(ok: bool, msg: str) -> None:
    if not ok:
        raise ValueError(msg)


def safe(model: str) -> str:
    return model.replace('+', 'p').replace('-', 'm')


def finite(x: Any) -> float | None:
    try:
        if x in (None, ''):
            return None
        x = float(x)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError, OverflowError):
        return None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_script(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        separators=(',', ':'),
        allow_nan=False,
    ).replace('<', '\\u003c')


def load_template_text(module_path: Path, *, markers: tuple[str, ...]) -> str:
    text = module_path.read_text(encoding='utf-8')
    for marker in markers:
        pos = text.find(marker)
        if pos >= 0:
            return text[pos:]
    return text


def find_state(base: Path, run_ids: list[str], model: str, case: str) -> tuple[dict | None, Path | None]:
    for run_id in run_ids:
        path = base / run_id / case / safe(model) / 'state.json'
        if path.is_file():
            doc = read_json(path)
            best = doc.get('best') or {}
            chi2 = finite(best.get('chi2'))
            if doc.get('status') in ACCEPTED and chi2 is not None:
                return doc, path
    return None, None


def _extract_masses(pred: dict[str, Any], cfg_path: Path) -> tuple[float | None, float | None, float | None]:
    masses = [
        finite(pred.get('largest_neutrino_mass_ev') or pred.get('mass_hi_ev')),
        finite(pred.get('middle_neutrino_mass_ev') or pred.get('mass_mid_ev')),
        finite(pred.get('smallest_neutrino_mass_ev') or pred.get('mass_lo_ev')),
    ]
    if any(v is None for v in masses) and cfg_path.is_file():
        cfg = read_json(cfg_path)
        masses = [
            finite(cfg.get('largest_neutrino_mass_ev')),
            finite(cfg.get('middle_neutrino_mass_ev')),
            finite(cfg.get('smallest_neutrino_mass_ev')),
        ]
    ordered = sorted(masses, reverse=True)
    return ordered[0], ordered[1], ordered[2]


def summary_rows(base: Path, run_ids: list[str]) -> tuple[list[dict[str, Any]], dict[str, str]]:
    rows: list[dict[str, Any]] = []
    provenance: dict[str, str] = {}
    for case in CASES:
        for model in MODELS:
            state, path = find_state(base, run_ids, model, case)
            if state is None or path is None:
                continue
            rel = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
            provenance[rel] = sha(path)
            best = state['best']
            pred = best.get('prediction') or {}
            cfg = path.parent / 'fitted_config.json'
            if cfg.is_file():
                rel_cfg = str(cfg.relative_to(ROOT)) if cfg.is_relative_to(ROOT) else str(cfg)
                provenance[rel_cfg] = sha(cfg)
            m_hi, m_mid, m_lo = _extract_masses(pred, cfg)
            rows.append({
                'model': model,
                'case': case,
                'run_id': path.parents[2].name,
                'chi2': finite(best.get('chi2')),
                'status': state.get('status'),
                'source_path': rel,
                'mass_hi_ev': m_hi,
                'mass_mid_ev': m_mid,
                'mass_lo_ev': m_lo,
                'delta_m21_sq_ev2': finite(pred.get('delta_m21_sq_ev2')),
                'delta_m3l_sq_ev2': finite(pred.get('delta_m3l_sq_ev2')),
                'sin2_theta12': finite(pred.get('sin2_theta12')),
                'sin2_theta13': finite(pred.get('sin2_theta13')),
                'sin2_theta23': finite(pred.get('sin2_theta23')),
            })
    rows.sort(key=lambda r: (r['case'], r['model']))
    return rows, provenance


def load_trajectories(trajectory_dir: Path | None, provenance: dict[str, str]) -> tuple[list[dict[str, Any]], str]:
    if trajectory_dir is None:
        return [], 'NONE_FOUND'
    path0 = (ROOT / trajectory_dir) if not trajectory_dir.is_absolute() else trajectory_dir
    if not path0.is_dir():
        return [], 'MISSING_DIRECTORY'
    rows: list[dict[str, Any]] = []
    for path in sorted(path0.rglob('*.json')):
        try:
            doc = read_json(path)
        except Exception:
            continue
        if not isinstance(doc, dict):
            continue
        model = doc.get('model')
        case = doc.get('case')
        uv = doc.get('uv_trajectory') or []
        eft = doc.get('eft1_trajectory') or []
        if model not in MODELS or case not in CASES:
            continue
        if not (isinstance(uv, list) and isinstance(eft, list) and uv and eft):
            continue
        rel = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
        provenance[rel] = sha(path)
        settings = doc.get('settings') or {}
        rows.append({
            'model': model,
            'case': case,
            'uv_trajectory': uv,
            'eft1_trajectory': eft,
            'scalar_endpoint_gev': finite(settings.get('scalar_endpoint_gev')),
            'fermion_threshold_gev': finite(settings.get('fermion_threshold_gev')),
            'uv_start_gev': finite(settings.get('uv_start_gev')),
        })
    return rows, ('SAVED_FROZEN_TRAJECTORIES' if rows else 'NO_VALID_TRAJECTORIES')


def load_sensitivity(csv_path: Path | None, provenance: dict[str, str]) -> tuple[list[dict[str, Any]], list[str], str]:
    if csv_path is None:
        return [], [], 'NONE_PROVIDED'
    path0 = (ROOT / csv_path) if not csv_path.is_absolute() else csv_path
    require(path0.is_file(), f'Missing sensitivity CSV: {path0}')
    rel = str(path0.relative_to(ROOT)) if path0.is_relative_to(ROOT) else str(path0)
    provenance[rel] = sha(path0)
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    skipped: list[str] = []
    with path0.open(newline='', encoding='utf-8-sig') as stream:
        for row in csv.DictReader(stream):
            model = row.get('model')
            case = row.get('case')
            direction = row.get('direction')
            gamma = finite(row.get('gamma'))
            if model not in MODELS or case not in CASES or direction not in ('T1', 'T2', 'T3'):
                continue
            if gamma is None or gamma <= 0:
                continue
            grouped.setdefault((model, case), []).append({
                'direction': direction,
                'gamma': gamma,
                'mass_hi_ev': finite(row.get('mass_hi_ev')),
                'mass_mid_ev': finite(row.get('mass_mid_ev')),
                'mass_lo_ev': finite(row.get('mass_lo_ev')),
                'delta_m21_sq_ev2': finite(row.get('delta_m21_sq_ev2')),
                'delta_m3l_sq_ev2': finite(row.get('delta_m3l_sq_ev2')),
                'sin2_theta12': finite(row.get('sin2_theta12')),
                'sin2_theta13': finite(row.get('sin2_theta13')),
                'sin2_theta23': finite(row.get('sin2_theta23')),
            })
    models: list[dict[str, Any]] = []
    for (model, case), pts in sorted(grouped.items()):
        directions = {p['direction'] for p in pts}
        if directions != {'T1', 'T2', 'T3'}:
            skipped.append(f'{model} {case}: missing one or more T1/T2/T3 directions')
            continue
        gammas = {d: sorted({p['gamma'] for p in pts if p['direction'] == d}) for d in directions}
        if not all(len(g) >= 3 and any(abs(x - 1.0) < 1e-12 for x in g) for g in gammas.values()):
            skipped.append(f'{model} {case}: insufficient gamma sampling')
            continue
        models.append({'model': model, 'case': case, 'points': pts})
    return models, skipped, 'SAVED_FROZEN_GAMMA_CSV'


def table_html(rows: list[dict[str, Any]]) -> str:
    headers = ['Model', 'Case', 'χ²', *LABELS.values(), 'Run']
    body = []
    for row in rows:
        vals = [row['model'], row['case'], row['chi2'], *(row[k] for k in OBS), row['run_id']]
        cells = ''.join(f'<td>{html.escape("" if v is None else str(v))}</td>' for v in vals)
        body.append(f'<tr>{cells}</tr>')
    return (
        '<table id="results"><thead><tr>'
        + ''.join(f'<th>{h}</th>' for h in headers)
        + '</tr></thead><tbody>'
        + ''.join(body)
        + '</tbody></table>'
    )


def build_index(rows: list[dict[str, Any]], trajectory_state: str, sensitivity_state: str) -> str:
    return f'''<!doctype html><html><head><meta charset="utf-8"><title>Frozen general-scalar reports</title>
<style>body{{font-family:Arial,sans-serif;max-width:1000px;margin:2rem auto;line-height:1.4}} code{{background:#f3f3f3;padding:.1rem .3rem}} li{{margin:.4rem 0}}</style></head>
<body><h1>Frozen general-scalar display index</h1>
<p>This report set is presentation-only. It uses accepted frozen-UV fits already saved under <code>output/studies/general_scalar_extended/bulk_48</code>.</p>
<ul>
<li><a href="interactive_rge.html">Interactive RGE report</a> — state: <b>{html.escape(trajectory_state)}</b></li>
<li><a href="sensitivity_report.html">T1/T2/T3 sensitivity report</a> — state: <b>{html.escape(sensitivity_state)}</b></li>
<li><a href="all_results.html">Accepted frozen-fit table</a> — rows: <b>{len(rows)}</b></li>
</ul></body></html>'''


def build_rge(rows: list[dict[str, Any]], trajectories: list[dict[str, Any]], template_hint: str) -> str:
    template = template_hint.strip() or '<html><body><h1>__TITLE__</h1><script>const DATA=__DATA__;</script>__NOTE__</body></html>'
    note = ('Saved frozen benchmark trajectories only. No matching or fit reruns are performed.' if trajectories
            else 'No saved frozen trajectory JSONs were supplied, so no RGE curves can be drawn yet.')
    payload = {'rows': trajectories, 'accepted_rows': rows}
    out = template.replace('__TITLE__', RGE_TITLE).replace('__DATA__', json_script(payload))
    return out.replace('__NOTE__', f'<p>{html.escape(note)}</p>')


def build_sensitivity(rows: list[dict[str, Any]], skipped: list[str], state: str, template_hint: str) -> str:
    template = template_hint.strip() or '<html><body><h1>__TITLE__</h1><script>const DATA=__DATA__;</script>__NOTE__</body></html>'
    note = ('Saved frozen T1/T2/T3 scans only — no historical reduced-potential envelopes are loaded.' if rows
            else 'No real frozen T1/T2/T3 gamma scan CSV was supplied, so no ovals have been rendered.')
    payload = {'models': rows, 'state': state}
    extra = ''.join(f'<li>{html.escape(item)}</li>' for item in skipped)
    out = template.replace('__TITLE__', SENS_TITLE).replace('__DATA__', json_script(payload))
    return out.replace('__NOTE__', f'<p>{html.escape(note)}</p>' + (f'<ul>{extra}</ul>' if extra else ''))


def build(
    run_ids: list[str],
    dest: Path,
    *,
    sensitivity_csv: Path | None = None,
    trajectory_dir: Path | None = None,
) -> dict[str, Any]:
    base = ROOT / OUT_BASE
    rows, provenance = summary_rows(base, run_ids)
    require(rows, 'No accepted frozen states were found in the selected run IDs')

    plotting = ROOT / 'Numerical' / 'plotting'
    rge_template = ''
    sens_template = ''
    rge_source = plotting / 'InteractiveModelComparison.py'
    sens_source = plotting / 'BuildFittedModelEnvelopes.py'
    if rge_source.is_file():
        rge_template = load_template_text(rge_source, markers=('const DATA=__DATA__', '<html'))
        provenance[str(rge_source.relative_to(ROOT))] = sha(rge_source)
    if sens_source.is_file():
        sens_template = load_template_text(sens_source, markers=('const D=__DATA__', '<html'))
        provenance[str(sens_source.relative_to(ROOT))] = sha(sens_source)

    trajectories, traj_state = load_trajectories(trajectory_dir, provenance)
    sens_models, skipped, sens_state = load_sensitivity(sensitivity_csv, provenance)

    dest.mkdir(parents=True, exist_ok=True)
    (dest / 'index.html').write_text(build_index(rows, traj_state, sens_state), encoding='utf-8')
    (dest / 'interactive_rge.html').write_text(build_rge(rows, trajectories, rge_template), encoding='utf-8')
    (dest / 'sensitivity_report.html').write_text(build_sensitivity(sens_models, skipped, sens_state, sens_template), encoding='utf-8')
    (dest / 'all_results.html').write_text(
        '<!doctype html><html><head><meta charset="utf-8"><title>Accepted frozen fits</title></head>'
        '<body><h1>Accepted frozen fits</h1>' + table_html(rows) + '</body></html>',
        encoding='utf-8',
    )
    with (dest / 'all_results.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=['model', 'case', 'chi2', *OBS, 'run_id', 'status', 'source_path'])
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        'accepted_row_count': len(rows),
        'trajectory_model_count': len(trajectories),
        'sensitivity_model_count': len(sens_models),
        'trajectory_state': traj_state,
        'sensitivity_state': sens_state,
        'run_ids': run_ids,
    }
    (dest / 'provenance.json').write_text(json.dumps({'files': provenance, 'manifest': manifest}, indent=2, sort_keys=True), encoding='utf-8')
    return manifest


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-ids', default='frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1')
    ap.add_argument('--sensitivity-csv', default='')
    ap.add_argument('--trajectory-dir', default='')
    ap.add_argument('--output-dir', default='')
    args = ap.parse_args(argv)

    run_ids = [x.strip() for x in args.run_ids.split(',') if x.strip()]
    if not run_ids:
        ap.error('No run IDs selected')
    out = ROOT / REPORT_BASE / (args.output_dir or '_'.join(run_ids))
    manifest = build(
        run_ids,
        out,
        sensitivity_csv=(Path(args.sensitivity_csv) if args.sensitivity_csv else None),
        trajectory_dir=(Path(args.trajectory_dir) if args.trajectory_dir else None),
    )
    print(
        f'Wrote {out} | {manifest["accepted_row_count"]} accepted rows | '
        f'{manifest["trajectory_model_count"]} trajectory models | '
        f'{manifest["sensitivity_model_count"]} sensitivity models'
    )
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
