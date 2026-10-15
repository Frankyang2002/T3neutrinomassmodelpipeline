"""Compare full SM+Weinberg running with a frozen-C5 reference.

Both predictions use the same UV/EFT1 RG trajectory, matching at mu_F and mu_S,
SM evolution, and low-scale charged-lepton basis. Only the final Weinberg K
matrix is frozen in the reference. Thus this isolates final-EFT K running,
NOT all quantum corrections or running above the scalar threshold.

Prerequisite: the artifact-aware RunFittedFermionScaleScan.py from the prior step.

PowerShell (project root):
 python -m Numerical.orchestration.RunFittedWeinbergRunningComparison
 python -m Numerical.orchestration.RunFittedWeinbergRunningComparison --run --limit 2
 python -m Numerical.orchestration.RunFittedWeinbergRunningComparison --run
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

from Numerical.orchestration.ResumePipelineNumericalResults import PROJECT_ROOT, parse_model_key
from Numerical.orchestration.RunFittedSensitivityComparison import best_candidate
from Numerical.orchestration.RunFittedFermionScaleScan import (
    OBS, SCENARIO, DEFAULT_FACTORS, discover_fitted_models,
    locate_artifacts, rescale_mass_payload, read_json,
)
from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
from Numerical.running.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.running.T3Trajectory import (
    run_renormalisable_t3_trajectory, continue_with_weinberg_running,
)
from Numerical.fitting.OscillationFit import (
    OscillationFitTarget, evaluate_mass_matrix_fit,
)
from Numerical.fitting.ParameterScan import (
    neutrino_mass_from_c5, charged_lepton_mass_basis_matrix,
)
from Numerical.fitting.ScanCLI import build_uv_state_from_config


def predictions(c5, ye, target, vev):
    mnu = neutrino_mass_from_c5(np.asarray(c5, dtype=complex), vev_gev=vev)
    ml = charged_lepton_mass_basis_matrix(mnu, ye)
    fit = evaluate_mass_matrix_fit(ml, target)
    p = {name: float(val) for name, val in zip(fit.observable_names, fit.prediction)}
    if not math.isfinite(float(fit.chi2)) or any(not math.isfinite(v) for v in p.values()):
        raise ValueError('Nonfinite oscillation prediction')
    return {'chi2': float(fit.chi2), 'prediction': p}


def evaluate(payload, uv_rge, eft_rge, c5_evaluator, target):
    scales = payload['scales']
    n = payload.get('numerical', {})
    vev = float(n.get('vev_gev', 246.22))
    rtol, atol = float(n.get('rtol', 1e-8)), float(n.get('atol', 1e-11))
    config = deepcopy(payload)
    config.setdefault('scan', {'bindings': {}})
    uv_state = build_uv_state_from_config(SimpleNamespace(raw=config), {})
    trajectory = run_renormalisable_t3_trajectory(
        uv_state, uv_rge, eft_rge,
        mu_fermion_threshold_gev=float(scales['mu_fermion_threshold_gev']),
        mu_scalar_threshold_gev=float(scales['mu_scalar_threshold_gev']),
        rtol=rtol, atol=atol,
    )
    # The fitted full calculation and reference share EXACTLY the same C5 matching.
    c5_at_s = np.asarray(c5_evaluator({}, trajectory), dtype=complex)
    full = continue_with_weinberg_running(
        trajectory, c5_at_s, float(scales['mu_low_gev']),
        rtol=rtol, atol=atol,
    )
    # Freeze K at scalar matching scale, but use the fully evolved low-scale ye
    # for BOTH predictions, avoiding a charged-lepton-basis mismatch.
    low_ye = full.weinberg.ye
    result_full = predictions(full.weinberg.K, low_ye, target, vev)
    result_frozen = predictions(c5_at_s, low_ye, target, vev)
    return result_full, result_frozen


def write_plots(rows, output):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not found; CSV/JSON are available.')
        return
    labels = {
        'delta_m21_sq_ev2': r'$\Delta m^2_{21}$',
        'delta_m3l_sq_ev2': r'$\Delta m^2_{3\ell}$',
        'sin2_theta12': r'$\sin^2\theta_{12}$',
        'sin2_theta13': r'$\sin^2\theta_{13}$',
        'sin2_theta23': r'$\sin^2\theta_{23}$',
    }
    for observable in OBS:
        fig, ax = plt.subplots(figsize=(9, 5))
        for model in sorted({r['model'] for r in rows}):
            subset = sorted((r for r in rows if r['model'] == model and r.get('factor') is not None),
                            key=lambda r: r['factor'])
            sections, section = [], []
            for r in subset:
                value = r.get('weinberg_running_percent', {}).get(observable)
                if r['status'] == 'Success' and value is not None and math.isfinite(value):
                    section.append(r)
                elif section:
                    sections.append(section); section = []
            if section: sections.append(section)
            for i, part in enumerate(sections):
                ax.plot([r['mass_ratio'] for r in part],
                        [r['weinberg_running_percent'][observable] for r in part],
                        marker='o', markersize=3, linewidth=1.2, label=model if i == 0 else None)
        ax.axhline(0, color='black', linewidth=.7)
        ax.set_xscale('log')
        ax.set_xlabel(r'$\mu_F/\mu_S$ (both $M_F$ and $\mu_F$ varied)')
        ax.set_ylabel(f'{labels[observable]}: evolved vs frozen $C_5$ (%)')
        ax.set_title('Final SM+Weinberg EFT: effect of running $C_5$ only')
        ax.grid(alpha=.2)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=5, ncol=2)
        fig.tight_layout()
        fig.savefig(output / f'weinberg_running_{observable}.png', dpi=180)
        plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scenario', default=SCENARIO)
    parser.add_argument('--factors', type=float, nargs='+', default=list(DEFAULT_FACTORS))
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--limit', type=int, default=0)
    args = parser.parse_args(argv)
    if args.limit < 0: parser.error('--limit must be >= 0')
    factors = sorted(set([float(f) for f in args.factors] + [1.0]))
    if any(not math.isfinite(f) or f <= 0 for f in factors):
        parser.error('Factors must be finite and >0')
    models = discover_fitted_models(args.scenario)
    if not models:
        print('No fitted models found. Expected Reports/output/full/fit_feasibility*/' + args.scenario)
        return 2
    target_path = PROJECT_ROOT / 'configs' / 'nufit_6_1_ic24_no.json'
    output = PROJECT_ROOT / 'Reports' / 'output' / 'full' / 'weinberg_running_comparison' / args.scenario
    rows = []
    for model in models[:args.limit or None]:
        selected = best_candidate(args.scenario, model)
        if selected is None:
            print(model, 'MissingFit', flush=True)
            rows.append({'model':model, 'status':'MissingFit', 'message':'Missing usable fitted config'})
            continue
        _, stage, config_path, _ = selected
        representation = parse_model_key(model)
        match, issue = locate_artifacts(model, representation)
        if match is None:
            print(model, 'MissingArtifacts', issue, flush=True)
            rows.append({'model':model, 'fit_stage':stage, 'status':'MissingArtifacts', 'message':issue})
            continue
        group, paths = match
        base = read_json(config_path)
        r = base['representation']
        if tuple(int(r[k]) for k in ('d_s1','d_s2','d_f','alpha')) != representation:
            print(model, 'InvalidFit: representation mismatch')
            rows.append({'model':model, 'fit_stage':stage, 'status':'InvalidFit', 'message':'representation mismatch'})
            continue
        scales = base['scales']
        mu_f, mu_s, mu_low = [float(scales[k]) for k in
                              ('mu_fermion_threshold_gev','mu_scalar_threshold_gev','mu_low_gev')]
        uv_mu = float(base['base_state']['mu_gev'])
        if not uv_mu > mu_f > mu_s > mu_low > 0:
            rows.append({'model':model, 'fit_stage':stage, 'status':'InvalidOrdering', 'message':'baseline scale ordering'})
            print(model,'InvalidOrdering')
            continue
        valid = [f for f in factors if uv_mu > mu_f*f > mu_s]
        print(f'{model}: artifact={group} valid factors={valid}', flush=True)
        if not args.run: continue
        if not target_path.is_file():
            rows.append({'model':model, 'fit_stage':stage, 'status':'MissingTarget', 'message':str(target_path)})
            continue
        try:
            uv_rge, eft_rge = load_rgbeta_payload(paths[0]), load_rgbeta_payload(paths[1])
            c5_evaluator = FinalC5TrajectoryEvaluator(paths[2])
            target = OscillationFitTarget.from_json(target_path)
        except Exception as exc:
            rows.append({'model':model, 'fit_stage':stage, 'status':'LoadFailed', 'message':f'{type(exc).__name__}: {exc}'})
            print(model, 'LoadFailed', exc, flush=True)
            continue
        for factor in factors:
            row = {'model': model, 'representation': list(representation), 'fit_stage': stage,
                   'artifact_set':group, 'factor':factor, 'mu_f_gev':mu_f*factor,
                   'mu_s_gev':mu_s, 'mass_ratio':mu_f*factor/mu_s,
                   'status':'', 'message':''}
            if not uv_mu > mu_f*factor > mu_s:
                row.update(status='OutsideEFTOrdering', message='Require mu_UV > mu_F > mu_S')
            else:
                try:
                    varied = rescale_mass_payload(base, factor)
                    full, frozen = evaluate(varied, uv_rge, eft_rge, c5_evaluator, target)
                    delta = {o: (100 * (full['prediction'][o] / frozen['prediction'][o] - 1)
                                 if frozen['prediction'][o] != 0 else None) for o in OBS}
                    row.update(status='Success', full=full, frozen_c5=frozen,
                               weinberg_running_percent=delta)
                except Exception as exc:
                    row.update(status='Failed', message=f'{type(exc).__name__}: {exc}')
            rows.append(row)
            print(model, f'factor={factor:g}',row['status'], row['message'][:120],flush=True)
    if not args.run:
        print('Dry-run only; use --run to compute.')
        return 0
    output.mkdir(parents=True,exist_ok=True)
    (output / 'scan.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    fields = ['model','representation','fit_stage','artifact_set','factor','mu_f_gev',
              'mu_s_gev','mass_ratio','status','message','chi2_full','chi2_frozen']
    fields += ['full_'+o for o in OBS]
    fields += ['frozen_'+o for o in OBS]
    fields += [o+'_weinberg_running_percent' for o in OBS]
    with (output / 'scan.csv').open('w',newline='',encoding='utf-8') as stream:
        writer = csv.DictWriter(stream,fieldnames=fields)
        writer.writeheader()
        for r in rows:
            record = {k:r.get(k,'') for k in fields}
            if 'representation' in r: record['representation'] = json.dumps(r['representation'])
            record['chi2_full'] = r.get('full',{}).get('chi2','')
            record['chi2_frozen'] = r.get('frozen_c5',{}).get('chi2','')
            for o in OBS:
                record['full_'+o] = r.get('full',{}).get('prediction',{}).get(o,'')
                record['frozen_'+o] = r.get('frozen_c5',{}).get('prediction',{}).get(o,'')
                value = r.get('weinberg_running_percent',{}).get(o)
                record[o+'_weinberg_running_percent'] = '' if value is None else value
            writer.writerow(record)
    write_plots(rows, output)
    status_count = {}
    for r in rows: status_count[r['status']] = status_count.get(r['status'],0)+1
    print('Status counts:', status_count)
    print('Wrote:',output)
    return 0 if status_count.get('Success',0) else 1


if __name__ == '__main__':
    raise SystemExit(main())
