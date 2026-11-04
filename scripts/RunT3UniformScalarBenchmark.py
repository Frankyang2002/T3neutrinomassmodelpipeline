"""Independent qualitative T3 benchmark: uniform BSM scalar quartics at UV.

At mu=1e7 GeV, all permitted BSM scalar quartics OTHER than lambdaT3
are +0.01 (real), allowed conjugate Yukawas are ZERO, and lambdaT3 is
scanned over 0.0005, 0.005, 0.5. The SM Higgs quartic, gauge couplings,
SM Yukawas, masses and B/C mixed scalar bilinear boundary are unchanged.

Use accepted merged 16-model fits ONLY as validated starting Yukawa seeds.
The existing physics modules own the RGE, matching and oscillation fit.
B/C initial conditions are built with ABSOLUTE quartics, avoiding the
existing extra_quartic_fraction<=1 interface without changing the physics.

C5 matching is conditional one-scale, NOT certified sequential EFT matching.
UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT has
NOT been independently established by these refits.
"""
from __future__ import annotations

import argparse
import csv
from copy import deepcopy
import json
import math
from pathlib import Path
import re
import traceback
from typing import Any

import numpy as np

from scripts import RunT3QualitativeScalarRefits as base

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/studies/thesis_qualitative/uniform_scalar_001'
QUARTIC = 0.01
SCHEMA = 't3-thesis-uniform-001-scalar-benchmark-v1'
ACCEPTED = frozenset(('AcceptedApproxFit', 'AcceptedRescaledSeed'))


def real_scalar_001(template: Any) -> Any:
    """Preserve the storage type used by the canonical UV state."""
    if isinstance(template, dict) and 'real' in template and 'imag' in template:
        if isinstance(template['real'], (float, int)) and isinstance(template['imag'], (float, int)):
            return {'real': QUARTIC, 'imag': 0.0}
    if isinstance(template, (float, int)) and not isinstance(template, bool):
        return QUARTIC
    raise ValueError(f'Unexpected scalar-quartic encoding: {type(template).__name__}')


def zero_matrix_like(value: Any) -> Any:
    """Zero a complex matrix without changing its canonical JSON encoding."""
    if isinstance(value, list):
        return [zero_matrix_like(v) for v in value]
    if isinstance(value, dict):
        if set(value) == {'real', 'imag'}:
            return {k: zero_matrix_like(v) for k, v in value.items()}
        raise ValueError('Unexpected conjugate Yukawa matrix dictionary')
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return 0.0
    raise ValueError('Unexpected conjugate Yukawa matrix encoding')


def uniform_config(config: dict, lam: float) -> tuple[dict, list[str]]:
    """Uniform BSM scalar sector. SM lambdaH and other boundaries untouched."""
    if not math.isfinite(lam) or lam <= 0:
        raise ValueError('lambdaT3 must be positive and finite')
    changed = deepcopy(config)
    ordinary = changed['base_state']['ordinary']
    names = []
    for name, value in list(ordinary.items()):
        if name.startswith('lambda') and name not in ('lambdaT3', 'lambdaH'):
            ordinary[name] = real_scalar_001(value)
            names.append(name)
        elif name in ('y1Conj', 'y2Conj'):
            ordinary[name] = zero_matrix_like(value)
    ordinary['lambdaT3'] = {'real': float(lam), 'imag': 0.0}
    return changed, sorted(names)


def uniform_overrides(overrides: dict) -> tuple[dict, list[str], list[str]]:
    """Canonical full-general UV layout, including newly allowed quartics."""
    out = deepcopy(overrides)
    quartics, conjugates = [], []
    for name, value in list(out.items()):
        if name.startswith('lambda') and name not in ('lambdaT3', 'lambdaH'):
            out[name] = real_scalar_001(value)
            quartics.append(name)
        if name in ('y1Conj', 'y2Conj'):
            out[name] = zero_matrix_like(value)
            conjugates.append(name)
    return out, sorted(quartics), sorted(conjugates)


def fixed_bc_extra(record: dict) -> dict[str, complex]:
    """Seven independent absolute B/C quartics, without a fraction cap."""
    names = tuple(record['new']) + tuple(record['reused'])
    if len(names) != 7 or len(set(names)) != 7:
        raise ValueError('Expected seven independent B/C quartics from the registry')
    return {name: complex(QUARTIC, 0.0) for name in names}


def bc_prepare_initial_absolute(record: dict, registry: dict, mu2_fraction: float):
    """Same native B/C state layout as the existing pilot, absolute quartics."""
    from scripts import CoupledUVRunning as s31
    from Numerical.core.StateVector import unpack_uv_state
    from Numerical.experimental.GeneralT3ScalarRGState import (
        GeneralT3ScalarRGState, pack_general_scalar_state)

    if not (math.isfinite(mu2_fraction) and 0 <= mu2_fraction <= 0.2):
        raise ValueError('Invalid original B/C mixed scalar bilinear fraction')
    new, reused = record['new'], record['reused']
    extra = fixed_bc_extra(record)
    cfg = record['cfg']
    ordinary = cfg['base_state']['ordinary']
    m1, m2 = float(ordinary['mS1Sq']), float(ordinary['mS2Sq'])
    if min(m1, m2) <= 0 or not np.isfinite([m1, m2]).all():
        raise ValueError('Nonpositive B/C UV scalar mass-squared')
    mu2 = complex(mu2_fraction * min(m1, m2), 0.0)
    vec0, base_layout = s31.prepare_state(cfg, record['entry'], reused, extra)
    state = unpack_uv_state(vec0, base_layout, mu_gev=1e7)
    general = GeneralT3ScalarRGState(state, {name:extra[name] for name in new}).validated(registry)
    vec, layout = pack_general_scalar_state(general, registry)
    return np.concatenate((vec, [mu2.real, mu2.imag])), layout, extra, mu2


def bc_integrate_point_absolute(model: str, record: dict, registry: dict,
                                bundle: dict, *, mu2_fraction: float):
    """Existing B/C RG equations and threshold projection; new UV boundary only.

    Mirrors RunT3GeneralScalarBCExploratoryFit.integrate_point except that
    b9.prepare_initial is replaced by bc_prepare_initial_absolute.
    """
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc
    from Numerical.experimental.GeneralT3ScalarRGState import (
        unpack_general_scalar_state, project_general_after_fermion,
        pack_general_scalar_state)

    y0, ly0, extra, mu2 = bc_prepare_initial_absolute(record, registry, mu2_fraction)
    new, reused = record['new'], record['reused']
    if record['uv_on']:
        uv = bc.b7.integrate_segment(
            y0, lambda t,y: bc.b7.rhs_at(t,y,layout=ly0,registry=registry,
                                         betas=bundle['UV'],stage='UV'),
            1e7, 1e5, n_points=3)
        yf, uv_nfev = uv.y[:,-1], int(uv.nfev)
    else:
        yf, uv_nfev = y0.copy(), 0
    at_f = unpack_general_scalar_state(yf[:ly0.size],ly0,registry,uv=True,mu_gev=1e5)
    heavy_state = at_f.baseline
    intermediate = project_general_after_fermion(at_f,registry)
    ye, lye = pack_general_scalar_state(intermediate,registry)
    ye = np.concatenate((ye, yf[ly0.size:]))
    betas = bc.b7.check_payload(
        bundle['payloads']['EFT1'],model=model,stage='EFT1',
        rep=lye.representation,base_names=tuple(lye.names[:-len(new)]),
        new_names=new,reused_names=reused)
    eft = bc.b7.integrate_segment(
        ye, lambda t,y: bc.b7.rhs_at(t,y,layout=lye,registry=registry,
                                     betas=betas,stage='EFT1'),
        1e5,1e3,n_points=3)
    end = bc.b7.snapshot(eft.y[:,-1],lye,registry,stage='EFT1',mu_gev=1e3,
                         new=new,reused=reused)
    eft_state = unpack_general_scalar_state(eft.y[:-2,-1],lye,registry,
                                            uv=False,mu_gev=1e3).baseline
    trajectory = {'eft1_end':end,'uv_nfev':uv_nfev,'eft1_nfev':int(eft.nfev),
                  'boundary_extra_quartics':{
                      k:{'real':float(z.real),'imag':float(z.imag)} for k,z in extra.items()},
                  'boundary_mu2_gev2':{'real':mu2.real,'imag':mu2.imag},
                  'boundary_policy':'absolute-0.01-no-fraction-compensation'}
    return trajectory, heavy_state, eft_state


class BCFixedQuarticEvaluator:
    """Reuse existing physical integrator, C5 matching and low-energy fit."""
    def __init__(self, model: str, source: dict, registry: dict, bundle: dict,
                 coefficient: Any, target: Any, mu2_fraction: float):
        self.model, self.source, self.registry = model, source, registry
        self.bundle, self.coefficient, self.target = bundle, coefficient, target
        self.mu2_fraction = mu2_fraction

    def evaluate(self, params: dict, *, detail: bool = False):
        from scripts import RunT3GeneralScalarBCExploratoryFit as bc
        values = np.array([bc.real_finite(params[k],k) for k in base.NAMES],dtype=float)
        trial = bc.trial_record(self.source, values)
        trajectory, uv_state, eft_state = bc_integrate_point_absolute(
            self.model,trial,self.registry,self.bundle,mu2_fraction=self.mu2_fraction)
        heavy = bc.heavy_flavour_at_muF(uv_state)
        c5 = bc.c5_at_endpoint(self.model,trajectory['eft1_end'],heavy,self.coefficient)
        fit, low = bc.low_energy_fit(c5,eft_state,self.target,
            vev=float(self.source['cfg'].get('numerical',{}).get('vev_gev',246.22)))
        if detail:
            return fit, {'trajectory':trajectory,
                         'MF_muF_gev':heavy[0].tolist(),
                         'heavy_mass_basis_residual':heavy[3],
                         'C5_muS_gev_inv':bc.b10.matrix_json(c5),
                         'C5_low_gev_inv':bc.b10.matrix_json(low.K)}
        return fit

    def __call__(self, params: dict):
        return self.evaluate(params)


class ModelFactory:
    """Compile matching and load beta payloads ONCE per model, not per point."""
    def __init__(self, model: str, config: dict, fit: dict, paths: dict,
                 target: Any, canonical_run: str):
        self.model, self.config, self.fit = model, config, fit
        self.paths, self.target = paths, target
        if model in base.MIXED:
            from scripts import RunT3BCGeneralFromReference as br
            from scripts import RunT3GeneralScalarBCExploratoryFit as bc
            from scripts import RunT3GeneralZeeSensitivity as zee
            self.source = br.build_reference_source(model)
            self.mu2_fraction = zee._bc_input_fractions(fit)[1]
            initial_cfg,_ = uniform_config(config,base.real_lambda(config))
            self.source['cfg'] = initial_cfg
            self.beta = bc.beta_bundle(model,self.source,self.source['registry'])
            self.coefficient = bc.compile_source(paths['matching'])
            self.scheme = 'conditional-mixed-scalar-one-scale-signed-quartic-map'
        else:
            from scripts import RunT3GeneralAll16FromReference as general
            self.uv,self.eft,_,self.kernel,source,self.scheme = general.nonmixed_inputs(
                model,canonical_run,generate_general_c5=False)
            if source.resolve() != paths['matching'].resolve():
                raise ValueError('Canonical general C5 differs from accepted fit source')
            frozen = base.json_read(paths['overrides'])
            self.frozen,self.all_quartics,self.all_conjugates = uniform_overrides(frozen)

    def build(self, lam: float):
        changed, ordinary_quartics = uniform_config(self.config,lam)
        if self.model in base.MIXED:
            source = dict(self.source)
            source['cfg'] = changed
            source['lambdaT3'] = complex(lam,0.0)
            evaluator = BCFixedQuarticEvaluator(
                self.model,source,source['registry'],self.beta,
                self.coefficient,self.target,self.mu2_fraction)
            extra = fixed_bc_extra(source)
            boundary = {
                'ordinary_quartics':ordinary_quartics,
                'additional_quartics':sorted(extra),
                'conjugate_yukawas_zero':[],
                'extra_quartic_value':QUARTIC,
                'extra_mixed_bilinear_fraction':self.mu2_fraction,
                'extra_mixed_bilinear_gev2':self.mu2_fraction*min(
                    changed['base_state']['ordinary']['mS1Sq'],
                    changed['base_state']['ordinary']['mS2Sq']),
                'matching':self.scheme,
                'signed_direct_quartic_map_physically_certified':False}
            fixed = None
        else:
            from scripts import RunT3GeneralAll16FromReference as general
            fixed = deepcopy(self.frozen)
            fixed['lambdaT3'] = {'real':float(lam),'imag':0.0}
            evaluator = general.NonMixedEvaluator(
                self.model,self.uv,self.eft,self.kernel,fixed,self.target)
            boundary = {'ordinary_quartics':ordinary_quartics,
                        'all_quartics_in_uv_layout':self.all_quartics,
                        'conjugate_yukawas_zero':self.all_conjugates,
                        'extra_quartic_value':QUARTIC,
                        'matching':self.scheme}
        return evaluator,changed,fixed,boundary


def run_point(factory: ModelFactory, lam: float, args, original: np.ndarray) -> tuple[dict,dict,dict|None]:
    result = {'model':factory.model,'lambdaT3':float(lam),'status':'Failed',
              'source_chi2':float(factory.fit['chi2_final']),
              'benchmark':SCHEMA}
    changed = None
    fixed = None
    try:
        evaluator,changed,fixed,boundary = factory.build(lam)
        result['scalar_boundary'] = boundary
        old = base.real_lambda(factory.config)
        seed = original * math.sqrt(old/lam)
        if np.max(np.abs(seed)) >= args.bound:
            raise base.UnsupportedScalarBoundary(
                'Rescaled Yukawa seed exceeds --bound; increase --bound explicitly')
        seed_fit = None
        if args.accept_seed:
            try:
                seed_fit = evaluator(base.parameter_dict(seed))
            except Exception:
                # An invalid heuristic seed does not imply the nearby fit fails.
                seed_fit = None
        if seed_fit is not None and math.isfinite(seed_fit.chi2) and seed_fit.chi2 <= args.chi2_cutoff:
            optimum = seed
            history = {'skipped':True,'reason':'Rescaled seed already passes chi2 cutoff',
                       'physics_evaluations':1,
                       'is_independent_numerical_optimisation':False}
            status = 'AcceptedRescaledSeed'
        else:
            _,optimum,history = base.optimize(evaluator,original,old,lam,
                bound=args.bound,max_nfev=args.max_nfev,
                regularization=args.regularization)
            status = None
        final,details = evaluator.evaluate(base.parameter_dict(optimum),detail=True)
        if status is None:
            status = 'AcceptedApproxFit' if float(final.chi2) <= args.chi2_cutoff else 'PoorFit'
        prediction = base.fit_prediction(final)
        mass = base.low_masses(details,changed,str(final.ordering))
        result.update(status=status,chi2=float(final.chi2),
                      prediction=prediction,neutrino_masses_ev=mass,
                      ordering=str(final.ordering),
                      whitened_residual=np.asarray(final.whitened_residual,dtype=float).tolist(),
                      fitted_yukawas=base.yukawa_matrices(optimum),
                      yukawa_norms=base.norms(optimum),
                      norm_ratio_to_reference=base.norms(optimum)['combined_frobenius']/
                                          base.norms(original)['combined_frobenius'],
                      optimizer=history,numerical_details=details)
    except Exception as exc:
        result.update(status='Failed',reason=f'{type(exc).__name__}: {exc}'[:1500],
                      traceback=traceback.format_exc(limit=5)[-2500:])
    return result,changed,fixed


def model_summary(rows: list[dict]) -> dict:
    accepted = sorted((r for r in rows if r['status'] in ACCEPTED),
                      key=lambda r:r['lambdaT3'])
    slope = None
    if len(accepted) >= 2:
        lo,hi=accepted[0],accepted[-1]
        a,b=lo['yukawa_norms']['combined_frobenius'],hi['yukawa_norms']['combined_frobenius']
        if min(a,b)>0:
            slope=math.log(b/a)/math.log(hi['lambdaT3']/lo['lambdaT3'])
    return {'model':rows[0]['model'],
            'status':'Complete' if len(accepted)==len(base.BENCHMARKS) else 'Incomplete',
            'accepted_count':len(accepted),
            'benchmark_status':{format(r['lambdaT3'],'.12g'):r['status'] for r in rows},
            'accepted_rescaled_seeds':sum(r['status']=='AcceptedRescaledSeed' for r in rows),
            'log_log_yukawa_norm_vs_lambda_slope_from_accepted_endpoints':slope,
            'comparison_slope_reference':-0.5,
            'slope_is_only_diagnostic':True}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--models',default='pilot',
                   help='pilot = A_-4,B_-1, all = 16, or comma-separated model names')
    p.add_argument('--canonical-run-id',default='full16_general_v1')
    p.add_argument('--run-id',default='uniform001_pilot_v1')
    p.add_argument('--max-nfev',type=int,default=5)
    p.add_argument('--bound',type=float,default=0.1)
    p.add_argument('--regularization',type=float,default=0.001)
    p.add_argument('--chi2-cutoff',type=float,default=5.86)
    p.add_argument('--accept-seed',action='store_true',
                   help='Skip expensive optimisation if the heuristic seed already passes chi2; record it as unoptimised')
    p.add_argument('--preflight',action='store_true',
                   help='Read and check merged sources without any physical evaluations or output')
    args=p.parse_args(argv)
    wanted = (('A_-4','B_-1') if args.models=='pilot' else
              base.MODELS if args.models=='all' else
              tuple(s.strip() for s in args.models.split(',') if s.strip()))
    if not wanted or len(set(wanted)) != len(wanted) or set(wanted)-set(base.MODELS):
        p.error('Unknown or duplicate model name')
    if not (re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id) and
            re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.canonical_run_id)):
        p.error('Invalid run id')
    if not (np.isfinite([args.max_nfev,args.bound,args.regularization,args.chi2_cutoff]).all()
            and args.max_nfev>=1 and args.bound>0 and 0<=args.regularization<=0.1
            and args.chi2_cutoff>0):
        p.error('Invalid numeric options')

    from scripts import RunT3GeneralZeeSensitivity as zee
    selection = zee.fit_selections(list(wanted),'merged')
    destination = OUT / args.run_id
    if not args.preflight and destination.exists():
        p.error(f'Refusing to overwrite existing run: {destination}')
    if not args.preflight:
        destination.mkdir(parents=True,exist_ok=False)
    target = None if args.preflight else zee.target_object()
    model_summaries,all_rows=[],[]
    for model in wanted:
        print(f'{model}: accepted merged source',flush=True)
        try:
            source_args = argparse.Namespace(fit_run_id='merged',
                            canonical_run_id=args.canonical_run_id)
            config,fit,paths = zee.sources(model,source_args,selection[model])
            original = base.fitted_vector(config)
            _,quartics = uniform_config(config,base.BENCHMARKS[0])
            if model in base.MIXED:
                from scripts import RunT3BCGeneralFromReference as br
                source = br.build_reference_source(model)
                extras=sorted(fixed_bc_extra(source))
                conjugates=[]
            else:
                _,all_quartics,conjugates=uniform_overrides(base.json_read(paths['overrides']))
                extras=sorted(set(all_quartics)-set(quartics))
            if args.preflight:
                print(f'  source={selection[model]["source"]}; '+
                      f'ordinary_quartics={quartics}; extras={extras}; '+
                      f'zero_conjugate_yukawas={conjugates}',flush=True)
                continue
            factory=ModelFactory(model,config,fit,paths,target,args.canonical_run_id)
            rows=[]
            for lam in base.BENCHMARKS:
                print(f'  lambdaT3={lam:g}; independent extra quartics={QUARTIC:g}',flush=True)
                row,changed,fixed=run_point(factory,lam,args,original)
                row['source']=selection[model]
                row['input_sha256']={k:base.sha(v) for k,v in paths.items()}
                row['input_paths']={k:v.relative_to(ROOT).as_posix() for k,v in paths.items()}
                row['original_lambdaT3']=base.real_lambda(config)
                row['original_yukawa_norms']=base.norms(original)
                rows.append(row)
                all_rows.append(row)
                folder=destination/base.tag(model)/base.lambda_tag(lam)
                base.json_save(folder/'result.json',row)
                if row['status'] in ACCEPTED or row['status']=='PoorFit':
                    # Preserve the full new UV boundary in the per-point config.
                    changed['base_state']['ordinary'].update(row['fitted_yukawas'])
                    changed['thesis_uniform_scalar_benchmark']={
                        'schema':SCHEMA,'quartic':QUARTIC,
                        'lambdaT3':lam,'conjugate_yukawas_zero':True,
                        'matching_scope':'conditional-one-scale'}
                    base.json_save(folder/'fitted_config.json',changed)
                    if fixed is not None:
                        from scripts import RunT3GeneralAll16FromReference as general
                        overrides=general.override_yukawas(fixed,base.parameter_dict(
                            base.fitted_vector(changed)))
                        base.json_save(folder/'general_uv_overrides.json',overrides)
                print('    '+row['status']+
                      (f', chi2={row["chi2"]:.6g}' if 'chi2' in row else ', '+row.get('reason','')[:150]),
                      flush=True)
            model_summaries.append(model_summary(rows))
        except Exception as exc:
            error={'model':model,'status':'Blocked','source':selection[model],
                   'reason':f'{type(exc).__name__}: {exc}'[:1500]}
            model_summaries.append(error)
            print('  BLOCKED: '+error['reason'],flush=True)
        if not args.preflight:
            base.json_save(destination/'summary.json',{
                'schema':SCHEMA,'quartic':QUARTIC,'benchmarks':list(base.BENCHMARKS),
                'models':model_summaries,'source_manifest_sha256':base.sha(zee.MERGED_FILE),
                'conditional_one_scale_not_sequential_fermion_scalar_matching':True,
                'conjugate_yukawas_zero':True,
                'accept_seed_when_passed':bool(args.accept_seed)})
    if args.preflight:
        print('Preflight only; no physical integration or files written.',flush=True)
        return 2 if model_summaries else 0
    fields=('model','lambdaT3','status','chi2',*base.OBS,'m1_ev','m2_ev','m3_ev',
            'y1_frobenius','y2_frobenius','combined_frobenius',
            'norm_ratio_to_reference','reason')
    with (destination/'comparison.csv').open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields)
        writer.writeheader()
        for row in all_rows:
            flat={k:row.get(k,'') for k in ('model','lambdaT3','status','chi2','reason','norm_ratio_to_reference')}
            flat.update(row.get('prediction',{}))
            flat.update(row.get('yukawa_norms',{}))
            for i,m in enumerate(row.get('neutrino_masses_ev',[]),1):
                flat[f'm{i}_ev']=m
            writer.writerow({k:flat.get(k,'') for k in fields})
    n_accepted=sum(r['status'] in ACCEPTED for r in all_rows)
    blocked=sum(s.get('status')=='Blocked' for s in model_summaries)
    print(f'Accepted {n_accepted}/{len(all_rows)}; blocked models {blocked}; '
          f'outputs in {destination}',flush=True)
    return 0 if not blocked and len(all_rows)==3*len(wanted) and n_accepted==len(all_rows) else 2


if __name__=='__main__':
    raise SystemExit(main())
