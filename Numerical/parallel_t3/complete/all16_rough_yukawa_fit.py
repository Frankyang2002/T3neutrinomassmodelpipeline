"""One independent 18-real-Yukawa rough fit for every T3 model.

This is a fresh-input, source-backed CONDITIONAL pipeline, not physical MSbar
one-loop matching. Every *unique* trial goes through:
    UV RGBeta -> F ALL OUT -> EFT1(S1,S2) RGBeta -> full tensor candidate
    Wilson RG -> BOTH scalars OUT at a common MS -> SMEFT RG -> PMNS.

The all-16 Workstreams B / rough Stream1 / Stream2 / Stream3 are the existing
repository implementations. The UV Yukawa override is injected into B's
per-run initializer; no archived F/Wilson tensor is reused across trials.

EXPLICIT APPROXIMATIONS retained from those workstreams:
 - conditional GammaCC -> Weyl signed factor (default +1), not certified;
 - one-loop finite fermion-hard C5 explicitly set zero for rough baseline;
 - scalar potential may be a reduced projection of full running potential;
 - full Wilson candidate ADM, mass/basis maps, matching not fully certified;
 - 18 real entries of Y1,Y2 only (not general complex+conjugate Yukawas).

A source/mass/Weinberg-purity failure is BLOCKED, not an artificial low chi2.
Scratch paths are uniquely owned; user data and existing outputs untouched.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Any
from unittest.mock import patch
from uuid import uuid4
import json
import math
import shutil
import traceback

import numpy as np
from scipy.optimize import least_squares

FORMAT = 't3-all16-sobol-sensitivity-rough-fit-target-chi2-v2'
YUKAWA_SCALE = 1e-4
SEED_MODELS = ('A_-4', 'D_-2')
BENCHMARK = Path('configs/thesis_numerics/t3_thesis_A4_Dm2_pilot.json')


class PipelineBlocked(RuntimeError):
    """Source check or numerical calculation cannot be performed."""


class TargetReached(Exception):
    """Normal early exit from least_squares after goal is evaluated."""


class TrialBudgetReached(Exception):
    """Normal early exit when the unique-trial budget is reached."""


def all_models() -> tuple[str, ...]:
    from Numerical.parallel_t3.threshold.models import ALL_MODELS
    result = tuple(ALL_MODELS)
    if len(result) != 16 or len(set(result)) != 16:
        raise PipelineBlocked('T3 all-model inventory is not the expected 16 unique models')
    return result


def normalize_models(models: list[str] | None) -> tuple[str, ...]:
    from Numerical.parallel_t3.threshold.models import normalize_model
    result = tuple(normalize_model(m) for m in models) if models else all_models()
    if not result or len(set(result)) != len(result):
        raise PipelineBlocked('Empty or duplicate model selection')
    return result


def read_json(path: Path) -> dict:
    obj = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(obj, dict):
        raise PipelineBlocked('Expected JSON object: ' + str(path))
    return obj


def save_json(path: Path, obj: dict) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.json.partial')
    temporary.write_text(json.dumps(obj, indent=2, sort_keys=True,
                                    allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def seed_pair(cfg: dict, seed: str) -> np.ndarray:
    """Dimensionless optimiser coordinates, x=Y/(1e-4), row-major 9+9."""
    if seed not in SEED_MODELS:
        raise PipelineBlocked('Not one of the two existing benchmark seed labels')
    if not math.isclose(float(cfg['fitting']['input_Yukawa_scale']), YUKAWA_SCALE):
        raise PipelineBlocked('Original benchmark Yukawa scale differs from 1e-4')
    record = cfg['models'][seed]
    y1 = np.asarray(record['y1_UV_times_1e4'], dtype=float)
    y2 = np.asarray(record['y2_UV_times_1e4'], dtype=float)
    if y1.shape != (3, 3) or y2.shape != (3, 3) or not (
            np.isfinite(y1).all() and np.isfinite(y2).all()):
        raise PipelineBlocked('Malformed real Yukawa 3x3 benchmark')
    return np.r_[y1.reshape(-1), y2.reshape(-1)]


def physical_yukawas(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    v = np.asarray(x, dtype=float)
    if v.shape != (18,) or not np.isfinite(v).all():
        raise PipelineBlocked('Trial must contain 18 finite real Yukawa entries')
    return v[:9].reshape(3, 3)*YUKAWA_SCALE, v[9:].reshape(3, 3)*YUKAWA_SCALE


def preflight(root: Path, model: str, *, batch1: str) -> dict:
    """Metadata check only; cannot predict Ward/purity or actual solver success."""
    from Numerical.parallel_t3.background.contracts import tag
    import os
    cfg = root/BENCHMARK
    b = root/batch1
    slug = tag(model)
    paths = [cfg, b/'models'/f'{slug}.json',
             b/'rgbeta'/slug/'uv.json', b/'rgbeta'/slug/'eft1.json']
    found = [str(p.relative_to(root)) for p in paths if p.is_file()]
    missing = [str(p.relative_to(root)) for p in paths if not p.is_file()]
    return {'model':model, 'input_ready':not missing, 'found':found, 'missing':missing,
            'evaluation_stage':'Actual Matchete F boundary and Ward/Weinberg checks run per trial'}


def _uv_overrides(cfg: dict, model: str, names: list[str], opts: Any,
                  x: np.ndarray, *, mass_bilinears: str) -> dict:
    """Shared thesis UV benchmark, with explicitly selected special mass policy.

    Use the already-tested Workstream B SM Yukawa defaults. Other scalar,
    gauge and heavy-mass parameters follow the SAME fixed thesis benchmark.
    Model-specific symmetry-allowed extra couplings are zero initially.
    """
    from Numerical.parallel_t3.background import runner as b_runner
    from Numerical.parallel_t3.background.contracts import encode_complex
    # _original_overrides is bound on the calling thread before patching.
    base = b_runner._original_overrides_for_fit16(model, names, opts)
    y1, y2 = physical_yukawas(x)
    result = dict(base)
    result.update({
        'MF':encode_complex(np.diag(cfg['heavy_fermion_masses_UV_GeV']).astype(complex)),
        'y1':encode_complex(y1.astype(complex)),
        'y2':encode_complex(y2.astype(complex)),
        'lambdaH':float(cfg['lambdaH_UV']),
        'lambdaT3':float(cfg['lambdaT3_UV']),
    })
    result.update({k:float(v) for k,v in cfg['gauge_UV'].items()})
    result.update({k:float(v) for k,v in cfg['scalar_masses_squared_UV_GeV2'].items()})
    for name in names:
        if name.startswith('lambda') and name not in ('lambdaH','lambdaT3'):
            result[name] = float(cfg['lambda_other_UV'])
        elif name.startswith('muT3General'):
            if mass_bilinears=='zero':
                result[name] = {'real':0.,'imag':0.}
            elif mass_bilinears=='nonzero-exploratory':
                result[name] = base[name]  # Workstream B's explicit model-specific input
            else:
                raise PipelineBlocked('Unknown scalar bilinear policy')
        elif name in ('y1Conj','y2Conj'):
            result[name]=encode_complex(np.zeros((3,3),complex))
    if not set(result) <= set(names):
        raise PipelineBlocked('Trial initial state has non-source beta symbols')
    return result


@dataclass
class Trial:
    observables: dict
    residuals: list[float]
    chi2: float
    meta: dict


class FreshAll16Evaluator:
    """Execute actual Workstreams B -> Stream1 -> Stream2 -> Stream3 per trial."""

    def __init__(self, root: Path, model: str, cfg: dict, *,
                 batch1: str, gamma_factor: float, mu_low: float,
                 vev: float, reduced_potential: bool, mass_bilinears: str,
                 wilson_solver: str = 'RK2', wilson_steps: int = 2,
                 keep_trial_files: bool = False):
        self.root = Path(root).resolve()
        self.model = model
        self.cfg = cfg
        self.batch1 = batch1
        self.gamma_factor = float(gamma_factor)
        self.mu_low = float(mu_low)
        self.vev = float(vev)
        self.reduced_potential = bool(reduced_potential)
        self.mass_bilinears = mass_bilinears
        self.wilson_solver = wilson_solver
        self.wilson_steps = int(wilson_steps)
        self.keep_trial_files = bool(keep_trial_files)
        if model not in all_models():
            raise PipelineBlocked('Unknown model')
        if wilson_solver not in ('RK2','DOP853') or self.wilson_steps<2:
            raise PipelineBlocked('Invalid full-tensor Wilson solver controls')
        if not 0 < mu_low <= cfg['scales_GeV']['S']:
            raise PipelineBlocked('SMEFT low scale must be at/below MS')
        if not np.isfinite(self.gamma_factor) or self.gamma_factor == 0:
            raise PipelineBlocked('Explicit finite nonzero GammaCC factor required')
        self.target = read_json(self.root/cfg['fitting']['target'])
        self.trial_no = 0

    def __call__(self, x: np.ndarray) -> Trial:
        from Numerical.parallel_t3.background import runner as b_runner
        from Numerical.parallel_t3.background.contracts import tag
        from Numerical.parallel_t3.rough_inputs import builder
        from Numerical.parallel_t3.wilson import rough_stream2
        from Numerical.parallel_t3.wilson.workflow import run_wilson
        from Numerical.parallel_t3.threshold.rough_adapter import prepare_pair
        from Numerical.parallel_t3.threshold.preflight import preflight_pair
        from Numerical.parallel_t3.threshold.rough_numerics import run_rough_case
        from Numerical.parallel_t3.complete.final_one_batch_fit import residual5

        self.trial_no += 1
        token = 'fit16_' + uuid4().hex[:10]
        slug = tag(self.model)
        # All five paths are owned by this exact trial token, never existing data.
        folders = [self.root/'output/parallel_t3'/part/token for part in
                   ('background','rough_inputs','rough_wilson','rough_results','fit16_adapted')]
        for p in folders:
            if p.exists():
                raise PipelineBlocked('Unexpected trial-token collision: '+str(p))
        stage = 'UV+EFT1'
        try:
            opts = b_runner.Options(root=self.root,run_id=token,batch1=self.batch1,
                  mu_uv=float(self.cfg['scales_GeV']['UV']),
                  mu_f=float(self.cfg['scales_GeV']['F']),
                  mu_s=float(self.cfg['scales_GeV']['S']),
                  uv_samples=5,eft1_samples=5,rtol=1e-8,atol=1e-11)
            original = b_runner._bench_overrides
            # Only a single-process, sequential optimisation. Never patch in parallel.
            def overrides(model, names, options):
                return _uv_overrides(self.cfg,model,names,options,x,
                                     mass_bilinears=self.mass_bilinears)
            with patch.object(b_runner,'_original_overrides_for_fit16',original,create=True):
                with patch.object(b_runner,'_bench_overrides',overrides):
                    upstream=b_runner.run_one(opts,self.model)
            if not str(upstream['fermion_boundary']).startswith('Formal'):
                raise PipelineBlocked('F matching unavailable: '+str(upstream['fermion_boundary']))
            stage = 'source-mapped conditional F tree & scalar tensor'
            aopts=builder.Options(root=self.root,B_run_id=token,out_run_id=token,
                      gamma_factor=self.gamma_factor,assume_formal_identity=True,
                      assume_sm_weyl_yukawa=True,reduced_potential=self.reduced_potential,
                      hard_c5_zero_benchmark=True,overwrite=False)
            prepared=builder.build_model(self.model,aopts)
            if prepared['status']!='written_provisional_A_C_inputs':
                raise PipelineBlocked('Stream1 failed: '+str(prepared))
            stage = 'EFT1 full-tensor Wilson running'
            case=rough_stream2.load_canonical(self.root,
                self.root/'output/parallel_t3/rough_inputs'/token/'a'/slug,
                self.model,c5_boundary='require',assume_formal_identity=True)
            evolved=run_wilson(case.basis,case.initial,case.background,
                case.mf,case.ms,adm_variant='fonseca_eq22_six_term',
                solver=self.wilson_solver,steps=self.wilson_steps,
                rtol=1e-6,atol=1e-21)
            wfolder=self.root/'output/parallel_t3/rough_wilson'/token/slug
            rough_stream2.write_c_endpoint(case,evolved,wfolder,
                adm_variant='fonseca_eq22_six_term',steps=self.wilson_steps,
                rtol=1e-6,atol=1e-21)
            # Existing writer's solver label is hard-coded DOP853. Correct it
            # when using its supported run_wilson RK2 mode (not a physics change).
            manifest=wfolder/'wilson.json'
            doc=read_json(manifest)
            doc['solver']['method']=self.wilson_solver
            save_json(manifest,doc)
            stage = 'source/basis and Weinberg purity preflight'
            bfolder=self.root/'output/parallel_t3/rough_inputs'/token/'c'/slug
            pair=prepare_pair(self.model,manifest,bfolder/'background.json',
                              self.root/'output/parallel_t3/fit16_adapted'/token/slug,
                              allow_synthetic=False)
            pf=preflight_pair(pair.wilson,pair.background,mode='rge_improved',
                              expected_model=self.model,mu_low=self.mu_low,
                              permit_synthetic=False,allow_projection_residual=False)
            if pf['status']!='ready_conditional':
                raise PipelineBlocked('Preflight blocked: '+' | '.join(pf['blockers'][:4]))
            stage = 'simultaneous S1/S2 matching + SMEFT RG + Takagi/PMNS'
            cfolder=self.root/'output/parallel_t3/rough_results'/token/slug
            details=run_rough_case(self.model,pair,cfolder,mu_low=self.mu_low,
                                  vev=self.vev,allow_synthetic=False,
                                  allow_projection_residual=False)
            if details['status']!='RoughNumericalBenchmark':
                raise PipelineBlocked('Non-genuine rough result: '+details['status'])
            obs=read_json(cfolder/'observables.json')
            residual=np.asarray(residual5(obs,self.target),float)
            if residual.shape!=(5,) or not np.isfinite(residual).all():
                raise PipelineBlocked('Invalid five-neutrino residual vector')
            chi2=float(np.dot(residual,residual))
            return Trial(obs,residual.tolist(),chi2,{
                 'model':self.model,'trial_number':self.trial_no,
                 'wilson_solver':self.wilson_solver,
                 'wilson_rhs_evaluations':evolved['metadata']['rhs_evaluations'],
                 'uv_and_eft1_recomputed':True,
                 'formal_F_boundary_recomputed':True,
                 'full_tensor_Wilson_recomputed':True,
                 'joint_scalar_matching_recomputed':True,
                 'SMEFT_running_recomputed':True,
                 'source_model_status':upstream['fermion_boundary'],
                 'scalar_mass_kind':details['scalar_mass_kind'],
                 'GammaCC_factor_assumed':self.gamma_factor,
                 'fermion_hard_finite_C5_assumed_zero':True,
                 'reduced_potential':self.reduced_potential,
                 'mass_bilinears_policy':self.mass_bilinears,
                 'conditional_not_physical':True,
                 'scratch_token':token if self.keep_trial_files else None})
        except Exception as exc:
            raise PipelineBlocked(f'{self.model} trial {self.trial_no}, stage {stage}: '
                                  f'{type(exc).__name__}: {exc}') from exc
        finally:
            if not self.keep_trial_files:
                for p in folders:
                    # No path outside an EXACT newly generated, private run-id.
                    if p.is_dir() and p.name==token and token.startswith('fit16_'):
                        shutil.rmtree(p)


def _search_box(start: np.ndarray, *, bound: float,
                perturbation: float) -> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    """Thesis-style bounded 18D neighbourhood in x = Y / 1e-4 coordinates."""
    start=np.asarray(start,float)
    if start.shape!=(18,) or not np.isfinite(start).all():
        raise PipelineBlocked('Invalid initial 18-component Yukawa vector')
    cap=float(bound)/YUKAWA_SCALE
    if not (np.isfinite(cap) and cap>np.max(abs(start)) and
            np.isfinite(perturbation) and 0<perturbation<=20):
        raise PipelineBlocked('Invalid perturbation or global Yukawa bound')
    half=float(perturbation)*np.maximum(np.abs(start),0.2)
    lo=np.maximum(start-half,-cap)
    hi=np.minimum(start+half,cap)
    if not np.all(hi>lo):
        raise PipelineBlocked('Degenerate bounded Sobol search region')
    centre=(lo+hi)/2
    widths=(hi-lo)/2
    return lo,hi,centre,widths


def optimize_one(evaluate: Callable[[np.ndarray], Trial], cfg: dict,
                 *, model: str, target_chi2: float, max_nfev: int,
                 max_trials: int, yukawa_bound: float,
                 sobol_points: int = 16, screen_points: int = 18,
                 active_schedule: tuple[int,...] = (5,9,18),
                 starts: int = 2, perturbation: float = 1.,
                 rng_seed: int = 17,
                 callback: Callable[[dict],None] | None = None) -> dict:
    """One model fit: two seeds -> Sobol -> sensitivity -> multistart TRF.

    Every *new* x is evaluated with the supplied actual EFT pipeline.
    The unique-trial cap includes both initial benchmarks, all Sobol samples,
    sensitivity probes, objective and Jacobian evaluations; the final physical
    pipeline recheck is independently repeated and reported separately.
    The candidate remains approximate, even if chi2 reaches the target.
    """
    from scipy.stats import qmc

    if not (np.isfinite(target_chi2) and 0<target_chi2<100 and
            max_nfev>=1 and max_trials>=2 and
            YUKAWA_SCALE<yukawa_bound<1 and
            sobol_points>=0 and 0<=screen_points<=18 and
            starts>=1 and isinstance(rng_seed,int) and
            active_schedule and all(isinstance(n,int) and 1<=n<=18 for n in active_schedule)):
        raise PipelineBlocked('Invalid Sobol / sensitivity / fit limits')

    records: list[dict]=[]
    cache: dict[bytes,np.ndarray]={}
    valid_points: dict[bytes,tuple[np.ndarray,np.ndarray,float]]={}
    best: dict | None=None
    seed_results: list[dict]=[]
    chosen_seed: str | None=None
    local_runs: list[dict]=[]
    sensitivity: list[dict]=[]
    scores=np.zeros(18,float)
    active_history: list[list[int]]=[]
    message='Not started'
    optimizer_success=False

    def residual(x: np.ndarray,phase: str) -> np.ndarray:
        nonlocal best
        x=np.asarray(x,dtype=float)
        if x.shape!=(18,) or not np.isfinite(x).all():
            raise PipelineBlocked('Optimiser supplied malformed Yukawa vector')
        key=x.tobytes()
        if key in cache:
            # No second expensive physics evaluation for the same exact point.
            return cache[key].copy()
        if len(records)>=max_trials:
            raise TrialBudgetReached()
        try:
            prediction=evaluate(x)
            vec=np.asarray(prediction.residuals,dtype=float)
            if vec.shape!=(5,) or not np.isfinite(vec).all():
                raise PipelineBlocked('Five observable residuals are invalid')
            chi2=float(vec@vec)
            if not np.isfinite(chi2):
                raise PipelineBlocked('Chi2 is not finite')
            row={'eval':len(records)+1,'phase':phase,'status':'VALID','chi2':chi2}
            valid_points[key]=(x.copy(),vec.copy(),chi2)
        except Exception as exc:
            # An invalid trial is never saved as an optimum.
            vec=np.full(5,1e6,dtype=float)
            chi2=float(vec@vec)
            row={'eval':len(records)+1,'phase':phase,'status':'BLOCKED',
                 'reason':f'{type(exc).__name__}: {exc}'[:750]}
        records.append(row)
        cache[key]=vec.copy()
        if row['status']=='VALID' and (best is None or chi2<best['chi2']):
            y1,y2=physical_yukawas(x)
            best={'chi2':chi2,'UV_y1':y1.tolist(),'UV_y2':y2.tolist(),
                  'x_scaled':x.tolist(),'residuals':vec.tolist(),
                  'observables':prediction.observables,
                  'trial_metadata':prediction.meta,'eval':row['eval'],
                  'search_phase':phase}
            if callback:
                callback({'status':'IN_PROGRESS','model':model,'best':best,
                          'unique_physics_evaluations':len(records),
                          'target_chi2_strict':target_chi2})
        if row['status']=='VALID' and chi2<target_chi2:
            raise TargetReached()
        return vec.copy()

    def ranked() -> list[tuple[np.ndarray,np.ndarray,float]]:
        return sorted(valid_points.values(),key=lambda v:v[2])

    def stage_fit(active: np.ndarray, centre: np.ndarray,
                  widths: np.ndarray, stage_number: int) -> None:
        nonlocal optimizer_success, message
        ordered=ranked()
        seen=set()
        selected=[]
        for x,vec,chi2 in ordered:
            # Do not launch repeated local minimisations from identical points.
            key=x.tobytes()
            if key not in seen:
                selected.append(x.copy());seen.add(key)
            if len(selected)>=starts:
                break
        for ordinal,frozen in enumerate(selected,1):
            z0=np.clip((frozen[active]-centre[active])/widths[active],
                       -1.+1e-10,1.-1e-10)
            def make_x(z):
                x=frozen.copy()
                x[active]=centre[active]+widths[active]*np.asarray(z,float)
                return x
            def f(z):
                return residual(make_x(z),f'least_squares_s{stage_number}_start{ordinal}')
            def jac(z):
                z=np.asarray(z,float)
                base=f(z)
                matrix=np.empty((5,len(active)),float)
                for j in range(len(active)):
                    # Direct, scale-normalised finite difference, as in the
                    # existing thesis Sobol fitter (not SciPy's tiny default).
                    step=.005 if z[j]<=.995 else -.005
                    probe=z.copy()
                    probe[j]=np.clip(probe[j]+step,-1.,1.)
                    delta=probe[j]-z[j]
                    matrix[:,j]=(f(probe)-base)/delta if abs(delta)>1e-12 else 0.
                return matrix
            info={'stage':stage_number,'start':ordinal,
                  'initial_chi2':next((v[2] for v in ordered if np.array_equal(v[0],frozen)),None),
                  'active_parameter_indices':active.tolist(),
                  'solver':'bounded_scipy_least_squares_TRF_explicit_Jacobian'}
            try:
                sol=least_squares(f,z0,jac=jac,method='trf',
                     bounds=(-np.ones(len(active)),np.ones(len(active))),
                     x_scale='jac',max_nfev=max_nfev,
                     ftol=1e-7,xtol=1e-7,gtol=1e-7)
                info.update(status='SOLVER_RETURNED',success=bool(sol.success),
                            nfev=int(sol.nfev),njev=int(sol.njev),
                            message=str(sol.message))
                optimizer_success=optimizer_success or bool(sol.success)
            except (TargetReached,TrialBudgetReached):
                # Record the local start before propagating the intentional stop.
                info['status']='STOPPED_ON_GOAL_OR_TRIAL_BUDGET'
                local_runs.append(info)
                raise
            except (ValueError,FloatingPointError,RuntimeError) as exc:
                # Do not discard other independent starts because one failed.
                info.update(status='SOLVER_BLOCKED',
                            reason=f'{type(exc).__name__}: {exc}'[:700])
            local_runs.append(info)
            message='Completed local stage '+str(stage_number)

    try:
        # Start from existing thesis benchmarks, not new made-up Yukawa inputs.
        for seed in SEED_MODELS:
            x=seed_pair(cfg,seed)
            try:
                residual(x,'initial_seed_'+seed)
            except TargetReached:
                seed_results.append({'seed':seed,'status':'VALID',
                                     'chi2':best['chi2']})
                chosen_seed=seed
                raise
            r=records[-1]
            seed_results.append({'seed':seed,'status':r['status'],
                                 'chi2':r.get('chi2') if r['status']=='VALID' else None,
                                 'reason':r.get('reason')})
        accepted=[r for r in seed_results if r['status']=='VALID']
        if not accepted:
            message='Both benchmark seeds blocked in the live EFT pipeline'
        else:
            chosen_seed=min(accepted,key=lambda row:row['chi2'])['seed']
            lo,hi,centre,widths=_search_box(seed_pair(cfg,chosen_seed),
                                bound=yukawa_bound,perturbation=perturbation)
            if sobol_points:
                n_draw=2**int(math.ceil(math.log2(sobol_points)))
                sobol=qmc.Sobol(d=18,scramble=True,
                        seed=int(rng_seed)+all_models().index(model))
                for unit in sobol.random_base2(m=int(math.log2(n_draw)))[:sobol_points]:
                    residual(lo+(hi-lo)*unit,'sobol')
            # Thesis-style dimensionless coordinate screening. Use Sobol
            # response differences for all coordinates; replace selected
            # proxy scores with direct finite differences around current best.
            anchor=np.asarray(best['x_scaled'],float)
            r0=np.asarray(best['residuals'],float)
            for x,rvec,_ in ranked():
                scores+=np.abs((x-anchor)/widths)*np.linalg.norm(rvec-r0)
            scores/=max(1,len(valid_points))
            for index in np.argsort(-scores,kind='stable')[:screen_points]:
                x=anchor.copy()
                step=.05*widths[index]
                if x[index]+step>hi[index]:
                    step=-step
                x[index]=np.clip(x[index]+step,lo[index],hi[index])
                actual=x[index]-anchor[index]
                if abs(actual)<=1e-14:
                    continue
                probe=residual(x,'sensitivity')
                if x.tobytes() in valid_points:
                    scores[index]=float(np.linalg.norm(probe-r0)*widths[index]/abs(actual))
                    sensitivity.append({'index':int(index),
                           'scaled_sensitivity_score':float(scores[index])})
            ordering=np.argsort(-scores,kind='stable')
            for stage_number,k in enumerate(dict.fromkeys(active_schedule),1):
                active=np.asarray(ordering[:k],dtype=int)
                active_history.append(active.tolist())
                stage_fit(active,centre,widths,stage_number)
    except TargetReached:
        message='Target reached strictly below chi2 threshold; stopped immediately'
    except TrialBudgetReached:
        message='Unique fresh-pipeline evaluation budget exhausted'

    if best is None:
        return {'model':model,'status':'BLOCKED_NO_VALID_SEED',
                'seed_comparisons':seed_results,'target_met':False,
                'optimizer_message':message,'trials':records,
                'unique_physics_evaluations':len(records),
                'physical_C5_certified':False}
    if chosen_seed is None:
        chosen_seed=min((r for r in seed_results if r['status']=='VALID'),
                        key=lambda r:r['chi2'])['seed']

    # Always rerun the complete numerical pipeline once at the BEST Yukawas.
    # This is NOT obtained from cached residual values.
    try:
        check=evaluate(np.asarray(best['x_scaled'],float))
        delta=abs(float(check.chi2)-float(best['chi2']))
        tol=max(1e-8,1e-2*max(abs(float(check.chi2)),abs(float(best['chi2']))))
        verified=bool(np.isfinite(check.chi2) and delta<=tol and
                      check.chi2>=0 and np.asarray(check.residuals).shape==(5,))
        verification={'passed':verified,'rechecked_chi2':float(check.chi2),
                      'initial_chi2':best['chi2'],'difference':delta,'tolerance':tol}
    except Exception as exc:
        verified=False
        verification={'passed':False,'error':f'{type(exc).__name__}: {exc}'[:1000]}
    hit=bool(verified and best['chi2']<target_chi2 and
             verification['rechecked_chi2']<target_chi2)
    return {'model':model,'status':('TARGET_REACHED_CONDITIONAL' if hit else
                   'CONDITIONAL_OPTIMIZATION_INCOMPLETE' if verified else
                   'BLOCKED_NONREPRODUCIBLE_BEST'),
            'initial_seed':chosen_seed,'seed_comparisons':seed_results,
            'best':best,'best_chi2':best['chi2'],'target_met':hit,
            'target_chi2_strict':target_chi2,'verification':verification,
            'optimizer_success':optimizer_success,
            'optimizer_message':message,
            'unique_physics_evaluations':len(records),
            'valid_physics_evaluations':sum(r['status']=='VALID' for r in records),
            'invalid_physics_evaluations':sum(r['status']=='BLOCKED' for r in records),
            'sobol_points_requested':sobol_points,
            'sensitivity_probes':sensitivity,
            'coordinate_sensitivity_scores':scores.tolist(),
            'active_parameter_stages':active_history,
            'local_least_squares_runs':local_runs,
            'max_nfev_per_local_start':max_nfev,
            'max_unique_trial_evaluations':max_trials,
            'trials':records,
            'physical_C5_certified':False,
            'complete_one_loop_matching_certified':False,
            'parameterization':'18 real UV Yukawa entries (2 x 3 x 3), row-major',
            'pipeline':'fresh UV RGE -> all F at MF -> EFT1 -> both S jointly at MS -> SMEFT',
            'search_method':'scrambled Sobol + scaled finite-difference sensitivity + multistart bounded TRF',
            'classification':'ROUGH_CONDITIONAL_NOT_PHYSICAL'}
