"""Cumulative ablation of the three RGE segments in the fermion-first T3 EFT.

Four predictions use the SAME existing fermion/scalar boundary builders and
FinalC5TrajectoryEvaluator. An 'off' segment transports its complete state
unchanged to the next matching scale. All predictions are diagonalised in the
same fully evolved low-energy charged-lepton basis for a controlled comparison.

Stage order is 000 (all frozen) -> 100 (UV on) -> 110 (UV+EFT1 on)
-> 111 (UV+EFT1+Weinberg on). Contributions depend on this chosen ordering.
These are RGE-switch-off diagnostics, NOT a separate matched perturbative order.

Requires the local artifact-aware RunFittedFermionScaleScan.py and the existing
T3 numerical pipeline. Does not modify existing physics files.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
from copy import deepcopy
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from Numerical.orchestration.support.ResumePipelineNumericalResults import PROJECT_ROOT, parse_model_key
from Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison import best_candidate
from Numerical.orchestration.studies.fermion_scales.RunFittedFermionScaleScan import (
    OBS, SCENARIO, DEFAULT_FACTORS, discover_fitted_models,
    locate_artifacts, rescale_mass_payload, read_json,
)
from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
from Numerical.core.StateVector import pack_uv_state
from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state
from Numerical.running.UVRunner import UVRunningResult, run_uv_segment
from Numerical.running.IntermediateScalarRunner import (
    IntermediateScalarRunningResult, run_intermediate_scalar_segment,
)
from Numerical.running.FermionThresholdBoundary import (
    fermion_threshold_masses, build_intermediate_scalar_boundary,
)
from Numerical.running.ScalarThresholdBoundary import (
    scalar_threshold_masses, build_final_sm_boundary,
)
from Numerical.running.T3Trajectory import (
    T3RenormalisableTrajectory, run_renormalisable_t3_trajectory,
    continue_with_weinberg_running,
)
from Numerical.running.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.fitting.OscillationFit import OscillationFitTarget
from Numerical.fitting.ScanCLI import build_uv_state_from_config
from Numerical.orchestration.studies.fitted_rge.RunFittedWeinbergRunningComparison import predictions

STAGES = ('000', '100', '110', '111')
STAGE_LABELS = {
    '000': 'No RGE segments',
    '100': 'UV RGE only',
    '110': 'UV + EFT1 RGEs',
    '111': 'UV + EFT1 + Weinberg RGEs',
}
INCREMENTS = (('UV', '000', '100'), ('EFT1', '100', '110'),
              ('Weinberg', '110', '111'))


def frozen_result(state, final_mu, *, uv):
    """Construct a zero-beta endpoint using the runner's validated layout.

    No beta functions or matching routines are patched or altered.
    """
    state = state.validated()
    mu0 = float(state.mu_gev)
    mu1 = float(final_mu)
    if not (mu0 > mu1 > 0):
        raise ValueError('Frozen segment requires descending positive scales')
    v, layout = pack_uv_state(state) if uv else pack_intermediate_scalar_state(state)
    kw = dict(mu_initial_gev=mu0, mu_final_gev=mu1, layout=layout,
              mu_gev=np.array([mu0, mu1], dtype=float),
              t=np.log(np.array([mu0, mu1], dtype=float)),
              y=np.column_stack((v, v)), solver_success=True,
              solver_message='Frozen state: beta set to zero by construction',
              nfev=0, njev=0, nlu=0)
    return UVRunningResult(**kw) if uv else IntermediateScalarRunningResult(**kw)


def with_uv_eft_states(uv, eft, mu_uv, mu_f, mu_s):
    uv_at_f = uv.final_state
    eft_initial = build_intermediate_scalar_boundary(uv_at_f, matching_scale_gev=mu_f)
    # 'eft' is an already integrated/frozen result beginning at this same state.
    if not np.allclose(eft.initial_vector, pack_intermediate_scalar_state(eft_initial)[0],
                       rtol=1e-11, atol=1e-11):
        raise ValueError('Intermediate-state boundary mismatch')
    eft_at_s = eft.final_state
    return T3RenormalisableTrajectory(
        mu_uv_initial_gev=mu_uv,
        mu_fermion_threshold_gev=mu_f,
        mu_scalar_threshold_gev=mu_s,
        uv=uv,
        fermion_threshold=fermion_threshold_masses(uv_at_f),
        eft1_initial_state=eft_initial,
        eft1=eft,
        scalar_threshold=scalar_threshold_masses(eft_at_s),
        final_sm_boundary=build_final_sm_boundary(eft_at_s, matching_scale_gev=mu_s),
    )


def evaluate_staged(payload, uv_beta, eft_beta, c5_evaluator, target):
    scales = payload['scales']
    mu_f, mu_s, mu_low = (float(scales[k]) for k in
        ('mu_fermion_threshold_gev', 'mu_scalar_threshold_gev', 'mu_low_gev'))
    numerics = payload.get('numerical', {})
    vev = float(numerics.get('vev_gev', 246.22))
    rtol, atol = float(numerics.get('rtol', 1e-8)), float(numerics.get('atol', 1e-11))
    config = deepcopy(payload)
    config.setdefault('scan', {'bindings': {}})
    uv_initial = build_uv_state_from_config(SimpleNamespace(raw=config), {})
    mu_uv = float(uv_initial.mu_gev)
    if not mu_uv > mu_f > mu_s > mu_low > 0:
        raise ValueError('Require mu_UV > mu_F > mu_S > mu_low > 0')

    # Real full renormalisable trajectory computed exactly as in benchmark scan.
    full_renorm = run_renormalisable_t3_trajectory(
        uv_initial, uv_beta, eft_beta,
        mu_fermion_threshold_gev=mu_f, mu_scalar_threshold_gev=mu_s,
        rtol=rtol, atol=atol,
    )
    matched_full = np.asarray(c5_evaluator({}, full_renorm), dtype=complex)
    evolved = continue_with_weinberg_running(
        full_renorm, matched_full, mu_low, rtol=rtol, atol=atol)
    # FIX the charged-lepton basis across all four variants, including 000.
    common_ye = evolved.weinberg.ye
    results = {'111': predictions(evolved.weinberg.K, common_ye, target, vev),
               '110': predictions(matched_full, common_ye, target, vev)}

    # UV on, EFT1 frozen; same F and S boundary functions + C5 evaluator.
    uv_real = full_renorm.uv
    eft_start = build_intermediate_scalar_boundary(uv_real.final_state,
                                                   matching_scale_gev=mu_f)
    eft_frozen = frozen_result(eft_start, mu_s, uv=False)
    uv_only = with_uv_eft_states(uv_real, eft_frozen, mu_uv, mu_f, mu_s)
    c5_uv_only = np.asarray(c5_evaluator({}, uv_only), dtype=complex)
    results['100'] = predictions(c5_uv_only, common_ye, target, vev)

    # Both renormalisable EFT segments frozen; threshold functions unchanged.
    uv_frozen = frozen_result(uv_initial, mu_f, uv=True)
    eft_no_uv_start = build_intermediate_scalar_boundary(uv_frozen.final_state,
                                                        matching_scale_gev=mu_f)
    eft_both_frozen = frozen_result(eft_no_uv_start, mu_s, uv=False)
    no_running = with_uv_eft_states(uv_frozen, eft_both_frozen,
                                   mu_uv, mu_f, mu_s)
    c5_no_running = np.asarray(c5_evaluator({}, no_running), dtype=complex)
    results['000'] = predictions(c5_no_running, common_ye, target, vev)
    return results


def fractional_percent(new, previous, observable):
    x = float(previous['prediction'][observable])
    y = float(new['prediction'][observable])
    if x == 0.0:
        return None
    return 100.0 * (y / x - 1.0)


def summarise(point):
    stages = point['stages']
    increments = {}
    for label, before, after in INCREMENTS:
        increments[label] = {o: fractional_percent(stages[after], stages[before], o)
                             for o in OBS}
    total = {o: fractional_percent(stages['111'], stages['000'], o) for o in OBS}
    # EXACT additive decomposition of log ratios for positive observables.
    logs = {}
    for o in OBS:
        values = [float(stages[k]['prediction'][o]) for k in STAGES]
        if min(values) <= 0:
            logs[o] = None
        else:
            parts = [100.0 * math.log(values[i+1]/values[i]) for i in range(3)]
            logs[o] = dict(zip((x[0] for x in INCREMENTS), parts))
            logs[o]['total'] = 100.0 * math.log(values[3]/values[0])
            logs[o]['sum_residual'] = sum(parts) - logs[o]['total']
    point['increment_percent'] = increments
    point['total_percent'] = total
    point['log_contributions'] = logs
    return point


def write_outputs(rows, output):
    output.mkdir(parents=True, exist_ok=True)
    (output / 'scan.json').write_text(json.dumps(rows, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    fields = ['model','representation','fit_stage','artifact_set','factor',
              'mu_f_gev','mu_s_gev','mass_ratio','status','message']
    fields += [f'{s}_chi2' for s in STAGES]
    fields += [f'{s}_{o}' for s in STAGES for o in OBS]
    fields += [f'{k}_{o}_percent' for k in ('UV','EFT1','Weinberg','total') for o in OBS]
    fields += [f'{k}_{o}_log_percent' for k in ('UV','EFT1','Weinberg','total') for o in OBS]
    with (output/'scan.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for p in rows:
            r = {k:p.get(k,'') for k in fields}
            if 'representation' in p:
                r['representation'] = json.dumps(p['representation'], separators=(',',':'))
            if p['status']=='Success':
                for s in STAGES:
                    r[f'{s}_chi2']=p['stages'][s]['chi2']
                    for o in OBS: r[f'{s}_{o}']=p['stages'][s]['prediction'][o]
                for k in ('UV','EFT1','Weinberg','total'):
                    block = p['total_percent'] if k=='total' else p['increment_percent'][k]
                    for o in OBS:
                        val=block[o]; r[f'{k}_{o}_percent'] = '' if val is None else val
                        logval=p['log_contributions'][o]
                        r[f'{k}_{o}_log_percent'] = '' if logval is None else logval[k]
            w.writerow(r)


def write_plots(rows, output):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib unavailable: skipped plots')
        return
    usable = [p for p in rows if p['status']=='Success']
    if not usable:
        return
    for o in OBS:
        fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
        for ax, label in zip(axes, ('UV','EFT1','Weinberg')):
            for model in sorted(set(p['model'] for p in usable)):
                selection = sorted((p for p in usable if p['model']==model),
                                   key=lambda x:x['factor'])
                ax.plot([p['mass_ratio'] for p in selection],
                        [p['increment_percent'][label][o] for p in selection],
                        marker='o', ms=3, lw=1.2, label=model)
            ax.axhline(0, color='black', lw=.6)
            ax.set_ylabel(f'{label} incremental change (%)')
            ax.grid(alpha=.25)
        axes[-1].set_xlabel(r'$\mu_F/\mu_S$')
        axes[-1].set_xscale('log')
        axes[0].set_title(f'Cumulative RGE ablation: {o}')
        axes[0].legend(fontsize=5, ncol=2, loc='best')
        fig.tight_layout()
        fig.savefig(output/f'staged_rg_{o}.png', dpi=170)
        plt.close(fig)


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--scenario',default=SCENARIO)
    ap.add_argument('--factors',nargs='+',type=float,default=list(DEFAULT_FACTORS))
    ap.add_argument('--limit',type=int,default=0)
    ap.add_argument('--run',action='store_true')
    args=ap.parse_args(argv)
    if args.limit < 0: ap.error('--limit must be >= 0')
    factors=sorted(set(args.factors+[1.0]))
    if any(not math.isfinite(f) or f<=0 for f in factors):
        ap.error('All mass factors must be finite and > 0')
    models=discover_fitted_models(args.scenario)
    if not models:
        print('No fitted models: Reports/output/full/fit_feasibility*/'+args.scenario)
        return 2
    target_path=PROJECT_ROOT/'configs'/'nufit_6_1_ic24_no.json'
    output=PROJECT_ROOT/'Reports'/'output'/'full'/'staged_rg_comparison'/args.scenario
    rows=[]
    for model in models[:args.limit or None]:
        rowbase={'model':model}
        selected=best_candidate(args.scenario,model)
        if selected is None:
            print(model,'MissingFit',flush=True)
            rows.append(dict(rowbase,status='MissingFit',message='No fitted config'))
            continue
        _,stage,config_path,_=selected
        rowbase['fit_stage']=stage
        rep=parse_model_key(model)
        match,issue=locate_artifacts(model,rep)
        if match is None:
            print(model,'MissingArtifacts:',issue,flush=True)
            rows.append(dict(rowbase,status='MissingArtifacts',message=issue))
            continue
        group,paths=match
        rowbase['artifact_set']=group
        rowbase['representation']=list(rep)
        base=read_json(config_path)
        parsed=tuple(int(base['representation'][k]) for k in
                     ('d_s1','d_s2','d_f','alpha'))
        if parsed!=rep:
            print(model,'InvalidFit: representation mismatch',flush=True)
            rows.append(dict(rowbase,status='InvalidFit',message='representation mismatch'))
            continue
        mu_f=float(base['scales']['mu_fermion_threshold_gev'])
        mu_s=float(base['scales']['mu_scalar_threshold_gev'])
        mu_low=float(base['scales']['mu_low_gev'])
        mu_uv=float(base['base_state']['mu_gev'])
        if not mu_uv>mu_f>mu_s>mu_low>0:
            rows.append(dict(rowbase,status='InvalidOrdering',message='baseline ordering'))
            continue
        valid=[f for f in factors if mu_uv>mu_f*f>mu_s]
        print(f'{model}: verified {group}, valid factors={valid}',flush=True)
        if not args.run: continue
        if not target_path.is_file():
            rows.append(dict(rowbase,status='MissingTarget',message=str(target_path)))
            continue
        try:
            uv_beta,eft_beta=load_rgbeta_payload(paths[0]),load_rgbeta_payload(paths[1])
            c5_evaluator=FinalC5TrajectoryEvaluator(paths[2])
            target=OscillationFitTarget.from_json(target_path)
        except Exception as e:
            rows.append(dict(rowbase,status='LoadFailed',message=f'{type(e).__name__}: {e}'))
            print(model,'LoadFailed:',e,flush=True)
            continue
        for f in factors:
            p=dict(rowbase,factor=f,mu_f_gev=mu_f*f,mu_s_gev=mu_s,
                   mass_ratio=mu_f*f/mu_s,status='',message='')
            if not mu_uv>mu_f*f>mu_s:
                p.update(status='OutsideEFTOrdering',message='Require mu_UV > mu_F > mu_S')
            else:
                try:
                    varied=rescale_mass_payload(base,f)
                    p['stages']=evaluate_staged(varied,uv_beta,eft_beta,c5_evaluator,target)
                    p['status']='Success'
                    summarise(p)
                except Exception as e:
                    p.update(status='Failed',message=f'{type(e).__name__}: {e}')
            rows.append(p)
            print(model,f'factor={f:g}',p['status'],p['message'][:150],flush=True)
    if not args.run:
        print('Dry run only. Use --run to evaluate.')
        return 0
    write_outputs(rows,output)
    write_plots(rows,output)
    counts={s:sum(p['status']==s for p in rows) for s in sorted(set(p['status'] for p in rows))}
    print('Status counts:',counts)
    print('Wrote:',output)
    return 0 if counts.get('Success',0)>0 else 1

if __name__=='__main__':
    raise SystemExit(main())
