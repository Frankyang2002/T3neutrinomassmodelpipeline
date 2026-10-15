"""Fitted T3 heavy-fermion mass/threshold sensitivity (fermion-first EFT).

Rescale the UV heavy mass matrix MF and fermion decoupling scale together,
keeping the UV starting scale, scalar threshold, other UV parameters and
matching/RGE machinery unchanged. This is a physical mass-hierarchy scan,
not an RG-on/off comparison or a matching-scale uncertainty test.

Run from project root:
 python -m Numerical.orchestration.RunFittedFermionScaleScan
 python -m Numerical.orchestration.RunFittedFermionScaleScan --run --limit 2
 python -m Numerical.orchestration.RunFittedFermionScaleScan --run
"""
from __future__ import annotations
import argparse
import csv
from copy import deepcopy
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from Numerical.orchestration.ResumePipelineNumericalResults import (
    PROJECT_ROOT, REQUIRED, parse_model_key,
)
from Numerical.orchestration.RunFittedSensitivityComparison import best_candidate
from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
from Numerical.running.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.fitting.OscillationFit import OscillationFitTarget
from Numerical.fitting.ParameterScan import make_t3_point_evaluator
from Numerical.fitting.ScanCLI import build_uv_state_from_config

OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
       'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
DEFAULT_FACTORS = (0.5, 0.75, 1.0, 1.5, 2.0)
SCENARIO = 'smallY_smallL'


ARTIFACT_SETS = ('hypercharge', 'dimensions', 'smoke')
REPRESENTATIONS = {
    'A': (1, 3, 2), 'B': (2, 2, 1), 'C': (2, 2, 3),
    'D': (3, 1, 2), 'E': (3, 3, 2),
}


def artifact_metadata(path):
    payload = read_json(path)
    metadata = payload.get('metadata')
    if not isinstance(metadata, dict):
        raise ValueError(f'No metadata object in {path}')
    keys = ('dS1', 'dS2', 'dF', 'alpha')
    if any(key not in metadata for key in keys):
        raise ValueError(f'Missing representation metadata in {path}')
    return tuple(int(metadata[key]) for key in keys)


def locate_artifacts(model, representation):
    """Find validated artifact trio without assuming the comparison layout."""
    d1, d2, df, alpha = representation
    matches = []
    examined = []
    for class_name, dimensions in REPRESENTATIONS.items():
        if dimensions != (d1, d2, df):
            continue
        code = ('m' + str(-alpha)) if alpha < 0 else ('p' + str(alpha))
        short = f'T3_{class_name}_alpha_{code}'
        for group in ARTIFACT_SETS:
            directory = PROJECT_ROOT / 'output' / 'full' / group / short / 'data'
            paths = tuple(directory / name for name in REQUIRED)
            if not all(path.is_file() for path in paths):
                examined.append(f'{group}/{short}: missing artifacts')
                continue
            try:
                actual_uv = artifact_metadata(paths[0])
                actual_eft = artifact_metadata(paths[1])
                if actual_uv != representation or actual_eft != representation:
                    examined.append(f'{group}/{short}: UV/EFT metadata {actual_uv}/{actual_eft}')
                    continue
                matches.append((group, paths))
            except Exception as exc:
                examined.append(f'{group}/{short}: {type(exc).__name__}: {exc}')
    if not matches:
        return None, 'No verified matching artifacts; ' + '; '.join(examined)
    # Prefer hypercharge when duplicate validated files occur across independent sets.
    return matches[0], ''


def discover_fitted_models(scenario):
    """Use fitted result directories as authority for model discovery."""
    from Numerical.orchestration.RunFittedSensitivityComparison import STAGES
    keys = set()
    for stage in STAGES:
        folder = PROJECT_ROOT / 'Reports' / 'output' / 'full' / stage / scenario
        if folder.is_dir():
            for model_dir in folder.iterdir():
                if model_dir.is_dir():
                    try:
                        parse_model_key(model_dir.name)
                    except ValueError:
                        continue
                    keys.add(model_dir.name)
    return sorted(keys)


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def rescale_mass_payload(source, factor):
    """Return independent config with both physical MF and F threshold rescaled."""
    payload = deepcopy(source)
    ordinary = payload['base_state']['ordinary']
    mass_spec = ordinary['MF']
    real = np.asarray(mass_spec['real'], dtype=float)
    imag = np.asarray(mass_spec.get('imag', np.zeros_like(real)), dtype=float)
    if real.shape != (3, 3) or imag.shape != (3, 3):
        raise ValueError('Expected 3x3 complex fermion mass matrix MF')
    mass_spec['real'] = (factor * real).tolist()
    mass_spec['imag'] = (factor * imag).tolist()
    payload['scales']['mu_fermion_threshold_gev'] *= factor
    return payload


def evaluate(payload, artifact_paths, target, loaded):
    uv_rge, eft_rge, c5_evaluator = loaded
    scales = payload['scales']
    numerical = payload.get('numerical', {})
    # Fitted candidate configs need not contain scan; the UV builder expects it.
    config_payload = deepcopy(payload)
    config_payload.setdefault('scan', {'bindings': {}})
    view = SimpleNamespace(raw=config_payload)
    evaluator = make_t3_point_evaluator(
        uv_state_builder=lambda _: build_uv_state_from_config(view, {}),
        c5_builder=c5_evaluator,
        uv_rgbeta_payload=uv_rge,
        eft1_rgbeta_payload=eft_rge,
        target=target,
        mu_fermion_threshold_gev=float(scales['mu_fermion_threshold_gev']),
        mu_scalar_threshold_gev=float(scales['mu_scalar_threshold_gev']),
        mu_low_gev=float(scales['mu_low_gev']),
        vev_gev=float(numerical.get('vev_gev', 246.22)),
        rtol=float(numerical.get('rtol', 1.0e-8)),
        atol=float(numerical.get('atol', 1.0e-11)),
    )
    fit = evaluator({})
    prediction = dict(zip(fit.observable_names, map(float, fit.prediction)))
    if not math.isfinite(float(fit.chi2)) or not all(math.isfinite(v) for v in prediction.values()):
        raise ValueError('Nonfinite fit output')
    return {'chi2': float(fit.chi2), 'prediction': prediction}


def plot_rows(rows, output):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib unavailable; CSV/JSON retained')
        return
    labels = {
        'delta_m21_sq_ev2': r'$\Delta m^2_{21}$',
        'delta_m3l_sq_ev2': r'$\Delta m^2_{3\ell}$',
        'sin2_theta12': r'$\sin^2\theta_{12}$',
        'sin2_theta13': r'$\sin^2\theta_{13}$',
        'sin2_theta23': r'$\sin^2\theta_{23}$',
    }
    for observable in OBS:
        fig, ax = plt.subplots(figsize=(9, 5.5))
        for model in sorted({r['model'] for r in rows}):
            subset = sorted((r for r in rows if r['model'] == model and 'factor' in r), key=lambda r:r['factor'])
            # Do not join across failed physical points.
            segment = []
            segments = []
            for r in subset:
                if r['status'] == 'Success' and math.isfinite(r.get('relative_percent', {}).get(observable, float('nan'))):
                    segment.append(r)
                elif segment:
                    segments.append(segment)
                    segment=[]
            if segment: segments.append(segment)
            for index, segment in enumerate(segments):
                ax.plot([r['mass_ratio'] for r in segment],
                        [r['relative_percent'][observable] for r in segment],
                        marker='o', ms=3, lw=1.3, label=model if index == 0 else None)
        ax.axhline(0, color='black', lw=.7)
        ax.set_xscale('log')
        ax.set_xlabel(r'$\mu_F/\mu_S$ (MF and fermion threshold varied together)')
        ax.set_ylabel(f'Change in {labels[observable]} from each model baseline (%)')
        ax.grid(alpha=.2)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=5, ncol=2, loc='best')
        fig.tight_layout()
        fig.savefig(output / f'fermion_scale_{observable}.png', dpi=180)
        plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scenario', default=SCENARIO)
    parser.add_argument('--factors', nargs='+', type=float, default=list(DEFAULT_FACTORS),
                        help='Multipliers applied to both UV MF and fermion matching scale')
    parser.add_argument('--run', action='store_true', help='Evaluate scan (default: inventory only)')
    parser.add_argument('--limit', type=int, default=0)
    args = parser.parse_args(argv)
    if args.limit < 0: parser.error('--limit must be nonnegative')
    factors = sorted(set(args.factors + [1.0]))
    if any(not math.isfinite(f) or f <= 0 for f in factors):
        parser.error('factors must be finite and >0')
    found = discover_fitted_models(args.scenario)
    if not found:
        print('No fitted model directories found for scenario:', args.scenario)
        print('Expected under Reports/output/full/{fit_feasibility,fit_feasibility_retry,fit_feasibility_focused}/')
        return 2
    output = PROJECT_ROOT / 'Reports' / 'output' / 'full' / 'fermion_scale' / args.scenario
    target_path = PROJECT_ROOT / 'configs' / 'nufit_6_1_ic24_no.json'
    rows = []
    for model in found[:args.limit or None]:
        selection = best_candidate(args.scenario, model)
        if selection is None:
            print('Missing usable fitted candidate:', model, flush=True)
            rows.append({'model': model, 'status': 'MissingFit', 'message': 'No fitted configuration/summary'})
            continue
        _, stage, config_path, _ = selection
        representation = parse_model_key(model)
        match, issue = locate_artifacts(model, representation)
        if match is None:
            print(model, issue, flush=True)
            rows.append({'model': model, 'status': 'MissingArtifacts', 'message': issue, 'fit_stage': stage})
            continue
        artifact_set, paths = match
        print(f'{model}: verified artifact set={artifact_set}', flush=True)
        if not target_path.is_file():
            print('Missing target:', target_path, flush=True)
            rows.append({'model': model, 'status': 'MissingTarget', 'message': str(target_path), 'fit_stage': stage})
            continue
        base = read_json(config_path)
        base_rep = base.get('representation', {})
        if tuple(int(base_rep[k]) for k in ('d_s1','d_s2','d_f','alpha')) != representation:
            print('Fitted representation mismatch:', model, flush=True)
            rows.append({'model': model, 'status': 'InvalidFit', 'message': 'Fitted representation mismatch', 'fit_stage': stage})
            continue
        uv_mu = float(base['base_state']['mu_gev'])
        scales = base['scales']
        mu_f = float(scales['mu_fermion_threshold_gev'])
        mu_s = float(scales['mu_scalar_threshold_gev'])
        mu_low = float(scales['mu_low_gev'])
        if not uv_mu > mu_f > mu_s > mu_low > 0:
            print('Invalid baseline ordering:', model, flush=True)
            rows.append({'model': model, 'status': 'InvalidOrdering', 'message': 'Invalid baseline EFT ordering', 'fit_stage': stage})
            continue
        valid = [f for f in factors if uv_mu > mu_f*f > mu_s]
        if not args.run:
            print(f'{model}: fit={stage}, UV={uv_mu:g}, F={mu_f:g}, S={mu_s:g}, valid factors={valid}', flush=True)
            continue
        try:
            loaded = (load_rgbeta_payload(paths[0]), load_rgbeta_payload(paths[1]), FinalC5TrajectoryEvaluator(paths[2]))
            target = OscillationFitTarget.from_json(target_path)
        except Exception as exc:
            print(model, 'load failed:', exc, flush=True)
            rows.append({'model': model, 'status': 'LoadFailed', 'message': str(exc), 'fit_stage': stage})
            continue
        model_rows = []
        for factor in factors:
            row = {'model':model,'representation':base.get('representation', {}), 'fit_stage': stage,
                   'factor':factor,'mu_f_gev':mu_f*factor,'mu_s_gev':mu_s,
                   'mass_ratio':mu_f*factor/mu_s, 'artifact_set': artifact_set, 'status':'', 'message':''}
            if not uv_mu > mu_f*factor > mu_s:
                row.update(status='OutsideEFTOrdering', message='Require mu_UV > mu_F > mu_S')
            else:
                try:
                    varied = rescale_mass_payload(base, factor)
                    result = evaluate(varied, paths, target, loaded)
                    row.update(status='Success', **result)
                except Exception as exc:
                    row.update(status='Failed', message=f'{type(exc).__name__}: {exc}')
            model_rows.append(row)
            print(model, f'factor={factor:g}', row['status'], row['message'][:100], flush=True)
        baseline = next((r for r in model_rows if r['factor'] == 1.0 and r['status']=='Success'), None)
        for row in model_rows:
            if baseline and row['status']=='Success':
                row['relative_percent'] = {o: 100*(row['prediction'][o]/baseline['prediction'][o]-1)
                                           if baseline['prediction'][o] != 0 else float('nan') for o in OBS}
            else:
                row['relative_percent'] = {}
        rows.extend(model_rows)
    if not args.run:
        print('Dry-run only. Add --run to evaluate.')
        return 0
    output.mkdir(parents=True, exist_ok=True)
    (output / 'scan.json').write_text(json.dumps(rows, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    columns = ['model','representation','fit_stage','factor','mu_f_gev','mu_s_gev',
               'mass_ratio','artifact_set','status','message','chi2'] + list(OBS) + [o+'_relative_percent' for o in OBS]
    with (output / 'scan.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for r in rows:
            flat = {k:r.get(k,'') for k in columns}
            flat['representation'] = json.dumps(r.get('representation', {}), separators=(',',':')) if 'representation' in r else ''
            for o in OBS:
                flat[o] = r.get('prediction', {}).get(o, '')
                flat[o+'_relative_percent'] = r.get('relative_percent', {}).get(o, '')
            writer.writerow(flat)
    plot_rows(rows, output)
    print('Wrote:', output)
    return 0 if rows and all(r['status']=='Success' for r in rows) else 1

if __name__ == '__main__':
    raise SystemExit(main())
