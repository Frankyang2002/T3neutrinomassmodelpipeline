"""Numerical-only heavy-fermion mass scan around saved common-input T3 benchmarks.

Uses the existing comparison RGBeta and final-C5 JSON artifacts without running
Matchete, RGBeta generation, fitting, or symbolic calculations.

From project root (PowerShell):
    python -m Numerical.orchestration.studies.common_inputs.RunCommonInputFermionScale
    python -m Numerical.orchestration.studies.common_inputs.RunCommonInputFermionScale --run --limit 2
    python -m Numerical.orchestration.studies.common_inputs.RunCommonInputFermionScale --run

Both the physical UV mass matrix MF and fermion threshold mu_F are scaled by
factor r. All remaining UV couplings, scalar parameters, mu_UV, mu_S and
mu_low are kept fixed. The frozen control preserves one-loop C5 matching,
but freezes parameters rather than solving RGEs. It is not full fixed-order EFT.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from Numerical.orchestration.support.ResumePipelineNumericalResults import (
    PROJECT_ROOT, OUTPUT_ROOT, REQUIRED, discover_cases, parse_model_key,
)
from Numerical.orchestration.studies.fermion_scales.RunHeavyFermionScaleComparison import (
    scale_heavy_mass, check_scales,
)
from Numerical.orchestration.studies.scalar_rge.RunCommonScalarComparison import evaluate
from Numerical.orchestration.studies.scalar_rge.RunCommonScalarRGControls import fixed_point
from Numerical.fitting.OscillationFit import OscillationFitTarget
from Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison import OBS

SCENARIO = 'smallY_smallL'
DEFAULT_FACTORS = (0.5, 0.75, 1.0, 1.5, 2.0)
SELECTED = (
    'T3_dS1_1_dS2_3_dF_2_alpha_m2',
    'T3_dS1_1_dS2_3_dF_2_alpha_p0',
    'T3_dS1_2_dS2_2_dF_1_alpha_m3',
    'T3_dS1_2_dS2_2_dF_1_alpha_p1',
    'T3_dS1_3_dS2_1_dF_2_alpha_m2',
    'T3_dS1_3_dS2_1_dF_2_alpha_p0',
)


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding='utf-8-sig'))


def signature(config_path: Path, artifact_paths: tuple[Path, ...], factors: list[float]) -> str:
    h = hashlib.sha256(config_path.read_bytes())
    h.update(json.dumps(factors, separators=(',', ':')).encode())
    for path in artifact_paths:
        stat = path.stat()
        h.update(str(path.resolve()).encode())
        h.update(str(stat.st_size).encode())
        h.update(str(stat.st_mtime_ns).encode())
    return h.hexdigest()


def save_plots(rows: list[dict], folder: Path) -> None:
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('Matplotlib not available: JSON and CSV saved; skipping PNG plots.')
        return
    folder.mkdir(parents=True, exist_ok=True)
    labels = {
        'delta_m21_sq_ev2': r'$\Delta m^2_{21}$ (eV$^2$)',
        'delta_m3l_sq_ev2': r'$\Delta m^2_{3\ell}$ (eV$^2$)',
        'sin2_theta12': r'$\sin^2\theta_{12}$',
        'sin2_theta13': r'$\sin^2\theta_{13}$',
        'sin2_theta23': r'$\sin^2\theta_{23}$',
    }
    palette = plt.get_cmap('tab10')
    models = list(dict.fromkeys(r['model'] for r in rows))
    for ob in OBS:
        fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 5))
        for i, model in enumerate(models):
            entries = sorted((r for r in rows if r['model'] == model), key=lambda r: r['factor'])
            # Never join through failed points.
            segments, segment = [], []
            for row in entries:
                if row['status'] == 'Success':
                    segment.append(row)
                elif segment:
                    segments.append(segment)
                    segment = []
            if segment:
                segments.append(segment)
            for j, pts in enumerate(segments):
                x = [r['mu_f_gev'] for r in pts]
                color = palette(i % 10)
                a0.plot(x, [r['running'][ob] for r in pts], '-o', lw=1.3, ms=3,
                        color=color, label=model if j == 0 else None)
                a0.plot(x, [r['fixed'][ob] for r in pts], '--', lw=1.0, color=color)
                a1.plot(x, [r['relative_percent'][ob] for r in pts], '-o', lw=1.3, ms=3,
                        color=color, label=model if j == 0 else None)
        a0.set(title='Solid: RGE; dashed: frozen parameters', ylabel=labels[ob])
        a1.set(title='RGE effect relative to frozen control',
               ylabel=r'$100(X_{\rm RG}-X_{\rm frozen})/X_{\rm frozen}$ (%)')
        for ax in (a0, a1):
            ax.set(xscale='log', xlabel=r'$\mu_F$ [GeV] ($M_F$ rescaled together)')
            ax.grid(alpha=0.25)
        a0.legend(fontsize=6, ncol=2)
        fig.suptitle('Common-input T3 heavy-fermion mass scan')
        fig.tight_layout()
        fig.savefig(folder / (ob + '_running_vs_frozen.png'), dpi=170)
        plt.close(fig)
    print('Saved plots:', folder)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scenario', default=SCENARIO)
    parser.add_argument('--run', action='store_true', help='Evaluate the numerical scan (otherwise inventory)')
    parser.add_argument('--limit', type=int, default=0, help='Number of selected models, 0=all six')
    parser.add_argument('--all-models', action='store_true', help='Scan all available common-input models')
    parser.add_argument('--factors', nargs='+', type=float, default=list(DEFAULT_FACTORS),
                        help='Positive factors scaling both MF and mu_F; 1.0 is automatically included')
    parser.add_argument('--force', action='store_true', help='Recompute even if matching cached results exist')
    args = parser.parse_args(argv)
    if args.limit < 0:
        parser.error('--limit must be nonnegative')
    if not args.factors or any(not math.isfinite(x) or x <= 0 for x in args.factors):
        parser.error('--factors must be positive and finite')
    factors = sorted(set(args.factors) | {1.0})
    common_root = PROJECT_ROOT / 'Reports' / 'output' / 'full' / 'common_input_rge' / args.scenario
    dest = PROJECT_ROOT / 'Reports' / 'output' / 'full' / 'common_input_fermion_scale' / args.scenario
    target_path = PROJECT_ROOT / 'configs' / 'nufit_6_1_ic24_no.json'
    cases = {model: path for scenario, model, path in discover_cases(OUTPUT_ROOT)
             if scenario == args.scenario}
    if not common_root.is_dir():
        raise FileNotFoundError('Common-input output not found: ' + str(common_root))
    models = sorted(p.name for p in common_root.iterdir() if p.is_dir() and
                    (p / 'common_input_config.json').is_file())
    if not args.all_models:
        models = [m for m in SELECTED if m in models]
    models = models[:args.limit or None]
    if not target_path.is_file():
        raise FileNotFoundError('NuFIT target missing: ' + str(target_path))
    target = OscillationFitTarget.from_json(target_path) if args.run else None
    all_rows = []
    failed = False
    print('Models:', len(models), 'factors:', factors)
    for model in models:
        config_path = common_root / model / 'common_input_config.json'
        if model not in cases:
            print('Missing comparison case:', model)
            failed = True
            continue
        case = cases[model]
        paths = tuple(case / 'data' / name for name in REQUIRED)
        missing = [p for p in paths if not p.is_file()]
        if missing:
            print('Missing existing numerical artifacts:', model, ', '.join(map(str, missing)))
            failed = True
            continue
        base = load(config_path)
        dims = parse_model_key(model)
        rep = base['representation']
        if tuple(int(rep[k]) for k in ('d_s1', 'd_s2', 'd_f', 'alpha')) != dims:
            raise ValueError('Config representation mismatch for ' + model)
        if base.get('benchmark_search', {}).get('enabled', False):
            raise ValueError('Refusing enabled fitting/search in ' + model)
        u, f, ss, low = check_scales(base)
        invalid = [x for x in factors if not u > f*x > ss]
        print(f'{model}: UV={u:g}, F={f:g}, S={ss:g}, valid factors=' +
              str([x for x in factors if x not in invalid]), flush=True)
        if not args.run:
            continue
        key = signature(config_path, paths, factors)
        saved = dest / model / 'mass_scan.json'
        if saved.is_file() and not args.force:
            prior = load(saved)
            if prior.get('signature') != key:
                raise RuntimeError('Different scan settings for ' + model + '. Use --force.')
            rows = prior['points']
            print('Reused cached:', model)
        else:
            rows = []
            for factor in factors:
                row = {'model': model, 'factor': factor, 'mu_f_gev': f*factor,
                       'mu_s_gev': ss, 'status': 'Failed', 'error': ''}
                if not u > f*factor > ss:
                    row.update(status='OutsideEFTOrdering', error='Require mu_UV > mu_F > mu_S')
                else:
                    try:
                        point = scale_heavy_mass(base, factor)
                        check_scales(point)
                        running = evaluate(point, list(paths), target)['prediction']
                        frozen = fixed_point(point, paths[2], target)
                        row['running'] = {o: float(running[o]) for o in OBS}
                        row['fixed'] = {o: float(frozen[o]) for o in OBS}
                        row['delta'] = {o: row['running'][o]-row['fixed'][o] for o in OBS}
                        row['relative_percent'] = {
                            o: (100*row['delta'][o]/row['fixed'][o]
                                if row['fixed'][o] != 0 else None) for o in OBS
                        }
                        if not all(math.isfinite(v) for d in (row['running'], row['fixed'], row['delta'])
                                   for v in d.values()):
                            raise ValueError('Non-finite observable')
                        row['status'] = 'Success'
                    except Exception as exc:
                        row['error'] = type(exc).__name__ + ': ' + str(exc)
                rows.append(row)
                print(f'  factor={factor:g}: {row["status"]} {row["error"]}', flush=True)
            saved.parent.mkdir(parents=True, exist_ok=True)
            saved.write_text(json.dumps({'signature': key, 'config_path': str(config_path),
                                         'artifact_paths': list(map(str, paths)),
                                         'factors': factors, 'points': rows}, indent=2) + '\n',
                             encoding='utf-8')
        all_rows.extend(rows)
        failed |= any(r['status'] != 'Success' for r in rows)
    if not args.run:
        print('Dry run only. Add --run to calculate.')
        return 1 if failed else 0
    if not all_rows:
        print('No numerical points produced.')
        return 1
    dest.mkdir(parents=True, exist_ok=True)
    (dest / 'scan.json').write_text(json.dumps(all_rows, indent=2) + '\n', encoding='utf-8')
    columns = ['model', 'factor', 'mu_f_gev', 'mu_s_gev', 'status', 'error']
    for ob in OBS:
        columns += [prefix + '_' + ob for prefix in ('running', 'fixed', 'delta', 'relative_percent')]
    with (dest / 'scan.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in all_rows:
            flat = {k: row.get(k, '') for k in columns[:6]}
            for ob in OBS:
                for prefix in ('running', 'fixed', 'delta', 'relative_percent'):
                    flat[prefix + '_' + ob] = row.get(prefix, {}).get(ob, '')
            writer.writerow(flat)
    save_plots(all_rows, dest / 'plots')
    print('Output:', dest)
    print('Successful:', sum(x['status'] == 'Success' for x in all_rows), '/', len(all_rows))
    return int(failed)


if __name__ == '__main__':
    raise SystemExit(main())
