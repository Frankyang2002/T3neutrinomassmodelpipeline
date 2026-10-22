"""Separate fermion-threshold-start T3 benchmarks: 5e-9 and 5e-10.

Fit 18 real Yukawa entries with Sobol+staged refinement, then cross-evaluate
one individually fitted A(-2) donor across all 16 representations. UV RGEs are
FROZEN, while fermion matching, EFT1 running, scalar matching and final SM+
Weinberg running are kept. No symbolic artifacts are regenerated.

The rescaled 5e-10 seeds come from individually fitted 5e-9 configurations.
This is a new boundary-value study: these inputs are interpreted at mu_F.
"""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from Numerical.orchestration import RunSmallLambdaRefit as base
from Numerical.orchestration import RunLambda5e9Refit as old
from Numerical.orchestration.studies.fitted_rge.RunFittedStagedRGComparison import (
    frozen_result, with_uv_eft_states,
)
from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
from Numerical.fitting.ScanCLI import build_uv_state_from_config
from Numerical.fitting.SobolBenchmarkSearch import (
    BenchmarkParameter, generate_sobol_parameter_points, refine_sobol_candidates,
)
from Numerical.fitting.OscillationFit import OscillationFitTarget, evaluate_mass_matrix_fit
from Numerical.fitting.ParameterScan import neutrino_mass_from_c5, charged_lepton_mass_basis_matrix
from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
from Numerical.running.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.running.T3Trajectory import continue_with_weinberg_running
from Numerical.orchestration.ModelNames import shorten

ROOT = base.ROOT / 'output' / 'studies' / 'fermion_start_lambda_comparison'
DONOR = 'T3_dS1_1_dS2_3_dF_2_alpha_m2'
FIELDS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
CASES = {'5e-9': (5e-9, 2.5e-9), '5e-10': (5e-10, 2.5e-10)}


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temp.replace(path)


def path_for(case, folder, model):
    return ROOT / case / folder / (model + '.json')


def cfg_lambda(case, row):
    return CASES[case][base.cls(row) == 'E']


def artifacts_for(model, mapping):
    return base.locate_artifacts(model, mapping)


def prepare(case):
    """Rescale 5e-9 fitted configurations for the 5e-10 seed only."""
    summary = []
    for row in base.original_rows():
        name = row['model']
        source = old.OUT / 'fitted' / (name + '.json')
        cfg = base.load(source)
        if cfg['representation'] != row['representation']:
            raise ValueError(f'Representation mismatch: {name}')
        uv = cfg['base_state']['ordinary']
        previous = float(uv['lambdaT3']['real'])
        target = cfg_lambda(case, row)
        if previous <= 0 or target <= 0:
            raise ValueError('Positive lambda required')
        factor = np.sqrt(previous / target)
        for key in ('y1', 'y2'):
            z = base.matr(uv[key])
            if np.max(np.abs(z.imag)) > 1e-12:
                raise ValueError(f'{name}: expected real fitted Yukawa matrices')
            base.setmat(uv[key], z * factor)
        uv['lambdaT3'] = {'real': target, 'imag': 0.0}
        cfg['benchmark_search'] = {'enabled': False}
        cfg['sensitivity'] = {'enabled': False}
        cfg['fermion_start_study'] = {
            'boundary_definition': 'UV initial parameter values imposed directly at mu_F; UV RGEs disabled',
            'case': case, 'lambdaT3': target, 'seed_factor': float(factor),
            'source': str(source.relative_to(base.ROOT)),
            'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        }
        save(path_for(case, 'prepared', name), cfg)
        summary.append({'model': name, 'lambdaT3': target, 'yukawa_factor': float(factor),
                        'max_yukawa_seed': float(max(np.max(abs(base.matr(uv[k]))) for k in ('y1', 'y2')))})
    save(ROOT / case / 'prepared_manifest.json', summary)
    print(f'Prepared {len(summary)} fermion-start configurations: {case}', flush=True)


class FermionStartEvaluator:
    """One 011 calculation, with matching and both EFT running regions intact."""
    def __init__(self, cfg, artifacts):
        self.cfg = copy.deepcopy(cfg)
        self.eft_beta = load_rgbeta_payload(artifacts[1])
        self.c5 = FinalC5TrajectoryEvaluator(artifacts[2])
        self.target = OscillationFitTarget.from_json(base.TARGET)
        self.names = [f'{k}_{i}{j}' for k in ('y1','y2') for i in range(3) for j in range(3)]
        self.vev = float(cfg.get('numerical', {}).get('vev_gev', 246.22))
        n = cfg.get('numerical', {})
        self.rtol = float(n.get('rtol', 1e-8))
        self.atol = float(n.get('atol', 1e-11))
        self.mu_f = float(cfg['scales']['mu_fermion_threshold_gev'])
        self.mu_s = float(cfg['scales']['mu_scalar_threshold_gev'])
        self.mu_low = float(cfg['scales']['mu_low_gev'])
        self.mu_uv = float(cfg['base_state']['mu_gev'])
        if not self.mu_uv > self.mu_f > self.mu_s > self.mu_low > 0:
            raise ValueError('Expected mu_UV > mu_F > mu_S > mu_low > 0')

    def evaluate(self, values, *, return_trajectory=False):
        config = copy.deepcopy(self.cfg)
        for key, offset in (('y1',0), ('y2',9)):
            matrix = np.array([float(values[self.names[j]]) for j in range(offset, offset+9)])
            base.setmat(config['base_state']['ordinary'][key], matrix.reshape(3,3))
        config.setdefault('scan', {'bindings': {}})
        initial = build_uv_state_from_config(SimpleNamespace(raw=config), {})
        uv = frozen_result(initial, self.mu_f, uv=True)
        eft_initial = build_intermediate_scalar_boundary(uv.final_state, matching_scale_gev=self.mu_f)
        eft = run_intermediate_scalar_segment(eft_initial, self.eft_beta, self.mu_s,
                                              rtol=self.rtol, atol=self.atol)
        trajectory = with_uv_eft_states(uv, eft, self.mu_uv, self.mu_f, self.mu_s)
        at_s = np.asarray(self.c5({}, trajectory), dtype=complex)
        low = continue_with_weinberg_running(trajectory, at_s, self.mu_low,
                                             rtol=self.rtol, atol=self.atol)
        matrix = neutrino_mass_from_c5(np.asarray(low.weinberg.K, dtype=complex), vev_gev=self.vev)
        flavour_matrix = charged_lepton_mass_basis_matrix(matrix, low.weinberg.ye)
        result = evaluate_mass_matrix_fit(flavour_matrix, self.target)
        if not np.isfinite(result.chi2):
            raise ValueError('Non-finite chi-square')
        if return_trajectory:
            return result, trajectory, low
        return result

    def __call__(self, point):
        return self.evaluate(point)


def vector_values(cfg):
    v = old.get_vector(cfg)
    return dict(zip([f'{k}_{i}{j}' for k in ('y1','y2') for i in range(3) for j in range(3)], map(float,v)))


def make_evaluator(case, model, mapping):
    cfg = base.load(path_for(case, 'prepared', model))
    return cfg, FermionStartEvaluator(cfg, artifacts_for(model, mapping))


def result_record(result, source, point):
    return {'source': source, 'chi2': float(result.chi2),
            'parameters': {k: float(v) for k,v in point.items()},
            'ordering': str(result.ordering),
            'prediction': dict(zip(result.observable_names, map(float,result.prediction))),
            'whitened_residual': np.asarray(result.whitened_residual,dtype=float).tolist()}


def fit(case, args):
    mapping = base.load(args.artifact_map) if args.artifact_map else None
    for number,row in enumerate(base.original_rows(),1):
        model = row['model']
        if args.model and args.model not in (model, shorten(model)):
            continue
        fit_path = path_for(case, 'fits', model)
        fitted_path = path_for(case, 'fitted', model)
        if args.resume and fit_path.exists() and fitted_path.exists():
            prior = base.load(fit_path)
            if prior.get('status') == 'Success' and float(prior.get('chi2_fit',np.inf)) <= args.accept_chi2:
                print(f'[{number}/16] accepted fit reused: {shorten(model)}', flush=True)
                continue
        cfg, evaluator = make_evaluator(case, model, mapping)
        seed = vector_values(cfg)
        names = list(seed)
        if any(abs(v) >= args.bound for v in seed.values()):
            raise ValueError(f'{shorten(model)}: seed exceeds bound={args.bound}; increase --bound')
        definitions = [BenchmarkParameter(k, f'ordinary.{k.split("_")[0]}.real.{k[-2]}.{k[-1]}',
                                          -args.bound, args.bound) for k in names]
        local = []
        for p in definitions:
            middle = seed[p.name]
            half = max(args.sobol_fraction*abs(middle), args.sobol_floor)
            local.append(BenchmarkParameter(p.name,p.binding,
                        max(-args.bound,middle-half),min(args.bound,middle+half)))
        checkpoint = path_for(case, 'sobol', model)
        stamp = {'case':case,'prepared_sha256':hashlib.sha256(path_for(case,'prepared',model).read_bytes()).hexdigest(),
                 'artifact_sha256':[hashlib.sha256(x.read_bytes()).hexdigest() for x in artifacts_for(model,mapping)],
                 'target_sha256':hashlib.sha256(base.TARGET.read_bytes()).hexdigest(),
                 'points':args.sobol_points,'seed':args.sobol_seed+number,
                 'bound':args.bound,'fraction':args.sobol_fraction,'floor':args.sobol_floor}
        prior = base.load(checkpoint) if checkpoint.exists() and args.resume else None
        if prior and prior.get('signature') != stamp:
            raise ValueError(f'{checkpoint}: incompatible resume; move old checkpoint first')
        candidates = list(prior.get('candidates',[])) if prior else []
        evaluated = int(prior.get('evaluated',0)) if prior else 0
        def retain(rec):
            candidates.append(rec)
            candidates.sort(key=lambda r:r['chi2'])
            del candidates[max(8,args.seed_count):]
        def checkpoint_now():
            save(checkpoint,{'signature':stamp,'evaluated':evaluated,'candidates':candidates})
        if not prior:
            try:
                seed_result = evaluator(seed)
                retain(result_record(seed_result,'exact_fermion_start_seed',seed))
                print(f'{shorten(model)} exact seed chi2={seed_result.chi2:.7g}',flush=True)
            except Exception as exc:
                print(f'{shorten(model)} exact seed failed: {exc}',flush=True)
            checkpoint_now()
        points = generate_sobol_parameter_points(local,args.sobol_points,seed=args.sobol_seed+number)
        for j in range(evaluated,len(points)):
            try: retain(result_record(evaluator(points[j]),f'sobol_{j+1}',points[j]))
            except Exception as exc:
                if j < 3: print('Sobol failure:',str(exc)[:200],flush=True)
            evaluated=j+1
            if evaluated%64==0 or evaluated==len(points):
                checkpoint_now()
                print(f'{shorten(model)} Sobol {evaluated}/{len(points)} best={candidates[0]["chi2"]:.6g}' if candidates
                      else f'{shorten(model)} Sobol {evaluated}/{len(points)} no valid point',flush=True)
        if not candidates: raise RuntimeError(f'No successful candidates for {model}')
        best = candidates[0]
        if best['chi2'] > args.accept_chi2:
            best, history = refine_sobol_candidates(
                evaluator=evaluator,observable_names=FIELDS,definitions=definitions,
                candidates=candidates,target_chi2=args.accept_chi2,
                seed_count=args.seed_count,active_counts=(5,8,12),max_nfev=(30,40,60),
                sensitivity_step=.02)
        else:
            history={'skipped':True,'reason':'Sobol seed accepted'}
        final = evaluator(best['parameters'])
        new_cfg = copy.deepcopy(cfg)
        for key,offset in (('y1',0),('y2',9)):
            base.setmat(new_cfg['base_state']['ordinary'][key],np.array([
                best['parameters'][names[j]] for j in range(offset,offset+9)]).reshape(3,3))
        new_cfg['fermion_start_study']['fit_chi2']=float(final.chi2)
        new_cfg['fermion_start_study']['fit_method']='Sobol + ranked 5/8/12 active parameter refinement'
        record={'model':model,'status':'Success' if final.chi2<=args.accept_chi2 else 'PoorFit',
                'chi2_fit':float(final.chi2),'lambdaT3':cfg_lambda(case,row),
                'prediction':dict(zip(final.observable_names,map(float,final.prediction))),
                'sobol_evaluated':evaluated,'refinement':history}
        save(fitted_path,new_cfg)
        save(fit_path,record)
        print(f'{shorten(model)} FINAL chi2={final.chi2:.7g} ({record["status"]})',flush=True)


def compare(case,args):
    mapping=base.load(args.artifact_map) if args.artifact_map else None
    donor_path=path_for(case,'fitted',DONOR)
    if not donor_path.exists(): raise FileNotFoundError(f'Fit A(-2) first: {donor_path}')
    donor=base.load(donor_path)
    donor_y={k:base.matr(donor['base_state']['ordinary'][k]) for k in ('y1','y2')}
    if args.mode=='shared' and base.load(path_for(case,'fits',DONOR)).get('status')!='Success':
        raise ValueError('A(-2) has not passed the fit threshold; refusing to label a poor point benchmark')
    mode='shared_a_m2' if args.mode=='shared' else 'optimised'
    out=ROOT/case/'comparison'/mode
    rows=[]
    for row in base.original_rows():
        model=row['model']
        src=path_for(case,'fitted',model)
        if not src.exists():
            rows.append({'model':model,'status':'MissingFit'});continue
        cfg=base.load(src)
        if args.mode=='shared':
            for key in ('y1','y2'):base.setmat(cfg['base_state']['ordinary'][key],donor_y[key])
        ev=FermionStartEvaluator(cfg,artifacts_for(model,mapping))
        try:
            fit=ev(vector_values(cfg))
            r={'model':model,'status':'Success','chi2':float(fit.chi2),
               'lambdaT3':cfg_lambda(case,row),
               'prediction':dict(zip(fit.observable_names,map(float,fit.prediction)))}
        except Exception as exc:
            r={'model':model,'status':'Failed','message':f'{type(exc).__name__}: {exc}'}
        rows.append(r)
        save(out/'models'/model/'evaluation.json',r)
        save(out/'models'/model/'fermion_start_config.json',cfg)
        print(f'{case} {mode} {shorten(model)}: {r.get("chi2",r["status"])}',flush=True)
    save(out/'summary.json',{'mode':mode,'case':case,'start_scale':'mu_F','uv_running':False,
                             'donor':DONOR if args.mode=='shared' else None,'models':rows})
    with (out/'comparison.csv').open('w',newline='',encoding='utf-8') as fh:
        columns=['model','status','lambdaT3','chi2',*FIELDS]
        writer=csv.DictWriter(fh,fieldnames=columns);writer.writeheader()
        for r in rows:writer.writerow({**{k:r.get(k,'') for k in columns},
                                       **{k:r.get('prediction',{}).get(k,'') for k in FIELDS}})
    print('Wrote',out,flush=True)


def display(case):
    """Saved-data endpoint comparison HTML (no fitting or ODE evaluation)."""
    from html import escape
    outputs=[]
    for mode in ('optimised','shared_a_m2'):
        root=ROOT/case/'comparison'/mode
        data=base.load(root/'summary.json')
        outputs.append((mode,data))
    dest=base.ROOT/'Reports'/'output'/'display_fermion_start'/case
    dest.mkdir(parents=True,exist_ok=True)
    names=[r['model'] for r in outputs[0][1]['models']]
    header='''<!doctype html><html lang="en"><meta charset="utf-8"><title>Fermion-threshold-start T3 comparison</title><style>body{font:15px system-ui;max-width:1200px;margin:30px auto;padding:0 18px}td,th{padding:8px;border-bottom:1px solid #ccc;text-align:left}table{border-collapse:collapse;width:100%}a{color:#185ca8}</style>'''
    for mode,data in outputs:
        cols=[r for r in data['models'] if r['status']=='Success']
        rows=''.join('<tr><td>'+escape(shorten(r['model']))+'</td><td>'+f'{r["chi2"]:.6g}'+'</td>'+''.join(f'<td>{r["prediction"][key]:.7g}</td>' for key in FIELDS)+'</tr>' for r in cols)
        text=header+f'<h1>Fermion start (UV off): {escape(case)} — {escape(mode)}</h1><p>Input boundary at μF; fermion matching, EFT1 running, scalar matching and Weinberg running are enabled. This is NOT full UV running.</p><p><a href="index.html">Overview</a></p><table><tr><th>Model</th><th>χ²</th>'+''.join(f'<th>{escape(x)}</th>' for x in FIELDS)+'</tr>'+rows+'</table></html>'
        (dest/(mode+'.html')).write_text(text,encoding='utf-8')
    both={mode:{r['model']:r for r in data['models']} for mode,data in outputs}
    rows=''
    for name in names:
        a=both['optimised'].get(name,{});b=both['shared_a_m2'].get(name,{})
        def cell(r):return f'{r["chi2"]:.6g}' if r.get('status')=='Success' else escape(r.get('status','Missing'))
        rows+=f'<tr><td>{escape(shorten(name))}</td><td>{cell(a)}</td><td>{cell(b)}</td></tr>'
    text=header+f'<h1>T3 {escape(case)} — fermion-threshold-start benchmarks</h1><p><strong>UV running disabled; EFT1 and Weinberg running enabled.</strong> Class E uses half the quoted quartic. Saved-data report; no physics rerun.</p><p><a href="optimised.html">Individual optimisation</a> · <a href="shared_a_m2.html">Shared A(-2)</a></p><table><tr><th>Model</th><th>Optimised χ²</th><th>Shared A(-2) χ²</th></tr>{rows}</table></html>'
    (dest/'index.html').write_text(text,encoding='utf-8')
    print('HTML:',dest/'index.html',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=('prepare','fit','compare','display'))
    parser.add_argument('--case',choices=tuple(CASES),required=True)
    parser.add_argument('--mode',choices=('optimised','shared'),default='optimised')
    parser.add_argument('--model',help='For fit stage, exact long key or short label')
    parser.add_argument('--sobol-points',type=int,default=64)
    parser.add_argument('--sobol-seed',type=int,default=20261022)
    parser.add_argument('--sobol-fraction',type=float,default=.5)
    parser.add_argument('--sobol-floor',type=float,default=.1)
    parser.add_argument('--bound',type=float,default=4.0)
    parser.add_argument('--seed-count',type=int,default=5)
    parser.add_argument('--accept-chi2',type=float,default=5.86)
    parser.add_argument('--artifact-map',type=Path)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if args.sobol_points<2 or args.sobol_points&(args.sobol_points-1):
        parser.error('--sobol-points must be a power of two >=2')
    if args.stage=='prepare':prepare(args.case)
    elif args.stage=='fit':fit(args.case,args)
    elif args.stage=='compare':compare(args.case,args)
    else:display(args.case)

if __name__=='__main__':main()
