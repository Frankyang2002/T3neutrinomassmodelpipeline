"""Read-only final reporting for the reduced T3 scalar benchmark.

The primary results are *only* the two frozen-UV, scenario-specific, newly
fitted A(-2) comparisons. The 16 fixed-boundary running-UV results are shown
as exploratory mass/hierarchy diagnostics. Optional UV trajectory replay runs
only the existing 10^7 -> 10^5 GeV renormalisable RGEs, with the historical
running-UV Yukawa inputs and exact zero additional quartic boundaries.

No optimiser, C5 rematching, physical scalar threshold retuning or GitHub write.

python -u -m scripts.BuildT3ReducedScalarFinalResults --run-id final_reduced_v1
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path

FROZEN_CASES = ('small_uv_frozen_5e-9', 'large_uv_frozen_0p5')
RUNNING_CASE = 'small_uv_running_5e-9'
FROZEN_FORMAT = 't3-reduced-zero-extra-sobol-shared-a2-v1'
FIT_CUTOFF = 5.86
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12',
       'sin2_theta13', 'sin2_theta23')
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'output/studies/reduced_scalar_final_reports'
INPUT_BASE = ROOT/'output/studies/reduced_scalar_sobol'
TARGET = ROOT/'configs/nufit_6_1_ic24_no.json'


def require(test, message):
    if not test:
        raise ValueError('T3 reduced-scalar reporting: '+str(message))


def read_json(path):
    path=Path(path)
    require(path.is_file(), 'missing file: '+str(path))
    return json.loads(path.read_text(encoding='utf-8-sig'))


def checksum(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def model_safe(model):
    return model.replace('+','p').replace('-','m')


def nonnegative_chi2(value, label):
    try:
        result=float(value)
    except (ValueError, TypeError):
        raise ValueError('T3 reduced-scalar reporting: missing chi2 '+str(label)) from None
    require(math.isfinite(result) and result>=0, 'invalid chi2 '+str(label))
    return result


def finite_or_none(value):
    if value is None or value=='':
        return None
    v=float(value)
    return v if math.isfinite(v) else None


def validate_frozen_payload(data, model_names):
    """Reject historical-A2, nonzero-extra, incomplete or duplicate campaigns."""
    settings=data.get('settings') or {}
    require(data.get('format')==FROZEN_FORMAT, 'unknown frozen comparison schema')
    require(tuple(data.get('cases',[]))==FROZEN_CASES,
            'must have exactly the two frozen-UV cases in declared order')
    require(settings.get('a2_source')=='new-fit', 'shared A(-2) is not newly fitted')
    require(settings.get('a2_as_start') is True,
            'A(-2) transfer was not available as an individual seed')
    require(settings.get('extra_quartic_boundary_real_imag')==[0.,0.] and
            settings.get('extra_quadratic_boundary')==0.,
            'additional couplings were not zero at input boundary')
    rows=data.get('rows')
    require(isinstance(rows,list) and len(rows)==32, 'expected exactly 32 frozen results')
    result={}
    for r in rows:
        key=(r.get('case'),r.get('model'))
        require(key not in result and key[0] in FROZEN_CASES and key[1] in model_names,
                'unexpected or duplicated model/case '+str(key))
        require(r.get('status')=='Accepted', 'individual case not accepted '+str(key))
        c=nonnegative_chi2(r.get('chi2_optimised'), key)
        require(c<=FIT_CUTOFF,'individual fit exceeds cutoff '+str(key))
        s=nonnegative_chi2(r.get('chi2_shared_a2'), 'shared '+str(key))
        require(r.get('shared_status')=='NeutrinoEvaluated',
                'shared reference not physically evaluated '+str(key))
        require(float(r.get('boundary_scale_gev',float('nan')))==1e5,
                'not a frozen fermion-scale input '+str(key))
        expected_lam=(5e-9 if key[0]==FROZEN_CASES[0] else .5)*(
            .5 if key[1].startswith('E_') else 1.)
        require(math.isclose(float(r.get('boundary_lambdaT3',float('nan'))),
                expected_lam,rel_tol=1e-12,abs_tol=0.),
                'scenario lambdaT3 mismatch '+str(key))
        for k in OBS:
            require(finite_or_none(r.get('fit_'+k)) is not None and
                    finite_or_none(r.get('shared_'+k)) is not None,
                    'missing oscillation observable '+k+' '+str(key))
        result[key]=r
    require(set(result)=={(case,model) for case in FROZEN_CASES for model in model_names},
            'model/case set not complete')
    for case in FROZEN_CASES:
        a=result[(case,'A_-2')]
        require(math.isclose(float(a['chi2_optimised']),float(a['chi2_shared_a2']),
                             rel_tol=1e-7,abs_tol=1e-12),
                'A(-2) shared reference does not reproduce its fit')
    return result


def validate_frozen_fit_files(source, rows):
    """Check that the campaign aggregate agrees with immutable per-model records."""
    verified={}
    for (case,model),r in rows.items():
        p=source/case/model_safe(model)/'fit.json'
        d=read_json(p)
        require(d.get('format')==FROZEN_FORMAT and d.get('model')==model and
                d.get('case')==case and d.get('status')=='Accepted',
                'fit identity/status mismatch: '+str(p))
        require(d.get('extra_quadratic_boundary_exact_zero') is True and
                d.get('quartic_one_loop_running_included') is True and
                d.get('general_scalar_one_loop_matching_included') is False,
                'reported approximation changed: '+str(p))
        for key,got,expected in (
                ('optimised', (d.get('best_mass_safe') or {}).get('chi2'),r['chi2_optimised']),
                ('shared',(d.get('shared_a2') or {}).get('chi2'),r['chi2_shared_a2'])):
            require(got is not None and math.isclose(float(got),float(expected),rel_tol=1e-9,abs_tol=1e-12),
                    key+' chi2 disagrees with aggregate: '+str(p))
        verified[(case,model)]=p
    return verified


def extract_running_endpoint(record, model):
    """The historical seed is a boundary diagnostic, not a fitted prediction."""
    require(record.get('format')==FROZEN_FORMAT and
            record.get('model')==model and record.get('case')==RUNNING_CASE,
            'invalid running trial identity '+model)
    trials=record.get('sobol_trials') or []
    first=[t for t in trials if t.get('label')=='historical_seed']
    require(len(first)==1, 'missing unique historical Yukawa seed for '+model)
    trial=first[0]
    masses=trial.get('scalar_masses_muS_tev') or []
    rho=finite_or_none(trial.get('rhoF'))
    require(len(masses)==2 and rho is not None and rho>=0 and
            all(finite_or_none(x) is not None and float(x)>0 for x in masses),
            'missing historical scalar endpoint for '+model)
    require(record.get('status')=='NoMassSafeNeutrinoFit',
            'running trial classification changed for '+model)
    require(trial.get('chi2') is None and trial.get('status')!='NeutrinoEvaluated',
            'running historical point unexpectedly has a neutrino fit '+model)
    return {'model':model,'status':record['status'],
            'initial_trial_label':'historical_seed',
            'mS1_formal_muS_tev':float(masses[0]),
            'mS2_formal_muS_tev':float(masses[1]),
            'rhoF_at_muF':rho,
            'screen_pass':str(bool(trial.get('screen',False))).lower(),
            'max_quartic_abs':finite_or_none(trial.get('max_scalar_quartic_abs')),
            'number_of_initial_trials':len(trials),
            'warning':'formal nominal-EFT1 endpoint; 1 TeV scalar matching not valid for heavy scalars'}


def load_running(source,model_names):
    rows={};paths=[]
    for model in model_names:
        p=source/RUNNING_CASE/model_safe(model)/'fit.json'
        d=read_json(p)
        rows[model]=extract_running_endpoint(d,model)
        paths.append(p)
    return rows,paths


def make_observable_rows(rows,central,errors):
    records=[]
    for case in FROZEN_CASES:
        for model in sorted({m for _,m in rows}):
            r=rows[(case,model)]
            for typ,prefix,c in (('individual','fit_',float(r['chi2_optimised'])),
                                 ('shared_A_-2','shared_',float(r['chi2_shared_a2']))):
                for obs in OBS:
                    val=float(r[prefix+obs]);sigma=float(errors[obs]);target=float(central[obs])
                    records.append({'case':case,'model':model,'yukawa_choice':typ,
                        'observable':obs,'prediction':val,'target':target,'sigma':sigma,
                        'pull':(val-target)/sigma,'chi2_five_observables':c})
    return records


def write_csv(path,rows,columns=None):
    require(bool(rows),'cannot write an empty result table')
    columns=columns or list(rows[0])
    with Path(path).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore')
        w.writeheader();w.writerows(rows)


def plot_frozen(rows,model_names,dest):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    f,axs=plt.subplots(2,1,figsize=(14,9),sharex=True,layout='constrained')
    x=np.arange(len(model_names))
    for ax,case in zip(axs,FROZEN_CASES):
        indiv=np.array([float(rows[(case,m)]['chi2_optimised']) for m in model_names])
        shared=np.array([float(rows[(case,m)]['chi2_shared_a2']) for m in model_names])
        ax.scatter(x-.15,np.maximum(indiv,1e-22),marker='o',s=36,label='Individual Yukawa fit')
        ax.scatter(x+.15,np.maximum(shared,1e-22),marker='x',s=48,label='Fixed A(−2) Yukawas')
        ax.axhline(FIT_CUTOFF,ls='--',lw=1,label=r'Acceptance $\chi^2=5.86$')
        ax.set_yscale('log');ax.set_ylabel(r'$\chi^2$')
        ax.set_title('Frozen UV — '+('small topology coupling' if case==FROZEN_CASES[0] else 'large topology coupling'))
        ax.grid(alpha=.2,axis='y');ax.legend(loc='lower left',fontsize=8,ncol=3)
    axs[-1].set_xticks(x,model_names,rotation=40,ha='right')
    axs[-1].set_xlabel('T3 representation')
    f.suptitle('Neutrino oscillation fit: independent versus common A(−2) Yukawas',fontsize=14)
    f.savefig(dest/'frozen_chi2_comparison.png',dpi=180)
    plt.close(f)


def plot_running_endpoint(running,model_names,dest):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    x=np.arange(len(model_names))
    m1=[running[m]['mS1_formal_muS_tev'] for m in model_names]
    m2=[running[m]['mS2_formal_muS_tev'] for m in model_names]
    f,ax=plt.subplots(figsize=(13.5,5.5),layout='constrained')
    ax.scatter(x-.14,m1,marker='o',s=39,label=r'$S_1$ nominal EFT1 endpoint')
    ax.scatter(x+.14,m2,marker='^',s=39,label=r'$S_2$ nominal EFT1 endpoint')
    ax.axhspan(.5,5,alpha=.12,label='Original 0.5–5 TeV window')
    ax.set_xticks(x,model_names,rotation=40,ha='right')
    ax.set_ylabel('Running scalar mass magnitude at nominal 1 TeV scale [TeV]')
    ax.set_title('Exploratory running-UV boundaries: fixed nominal EFT1 extrapolation')
    ax.grid(alpha=.2,axis='y');ax.legend(fontsize=8)
    f.savefig(dest/'running_uv_formal_mass_endpoints.png',dpi=180)
    plt.close(f)



def select_uv_betas(model, entry, new, reused, layout, *, bulk, s31):
    """Validate a UV export against the *complete* model-specific quartic inventory.

    The historical bulk fitting helper intentionally passed an empty reused
    inventory to Stage31 because its B(-1)/C(-1) fits were handled separately.
    A reduced-sector UV replay, however, evolves those reused coordinates in
    the ordinary UV state (from a zero UV boundary). The Stage31 export lists
    all six reused B/C symbols, so its validator must receive that inventory.
    Do not modify the original bulk campaign or weaken source validation.
    """
    if not reused:
        # Keeps the existing special Stage30 A(+0) handling and Stage31
        # handling for models without preexisting additional quartic slots.
        return bulk.uv_betas(model, entry, new, layout)
    require(model in ('B_-1', 'C_-1'),
            'unexpected reused UV quartics in model '+model)
    beta_path=bulk.BETA_DIR/f'stage31_{bulk.safe(model)}_full_uv_betas.json'
    payload=bulk.read(beta_path)
    betas=s31.validate_export(payload, entry, new, reused, layout.names)
    return betas,beta_path


def sample_uv_trajectory(model, run_source):
    """Actually integrate the candidate one-loop UV system with 0 extra couplings.

    No scalar-only EFT1 continuation is presented as physically consistent when
    the mass hierarchy invalidates the nominal scalar threshold.
    """
    import numpy as np
    from scipy.integrate import solve_ivp
    from scripts import RunT3GeneralScalar48Campaign as bulk
    from scripts import Stage31CoupledUVRunning as s31
    from scripts import Stage30CoupledA0UVRunning as s30
    from Numerical.core.StateVector import unpack_uv_state
    registry=bulk.read(bulk.REGISTRY)
    entry,new,reused=bulk.inventory(registry,model)
    hist=bulk.historical(RUNNING_CASE,model,entry)
    require(hist['uv_on'], 'not a UV-running historical input')
    cfg=bulk.update_yukawas(hist,hist['seed'])
    q={name:0j for name in reused}
    vec,layout=s31.prepare_state(cfg,entry,reused,q)
    betas,beta_path=select_uv_betas(model,entry,new,reused,layout,bulk=bulk,s31=s31)
    ini=np.r_[vec,np.zeros(2*len(new),dtype=float)]
    tt=np.linspace(math.log(1e7),math.log(1e5),80)
    sol=solve_ivp(lambda t,y:s30.combined_rhs(t,y,layout=layout,names=new,betas=betas),
                  (tt[0],tt[-1]),ini,t_eval=tt,rtol=1e-8,atol=1e-11,method='DOP853')
    require(sol.success and np.isfinite(sol.y).all(), 'UV trajectory failed: '+model)
    points=[]
    for i,t in enumerate(sol.t):
        mu=math.exp(float(t))
        state=unpack_uv_state(sol.y[:layout.size,i],layout,mu_gev=mu)
        m1sq=float(state.mS1Sq)/1e6
        m2sq=float(state.mS2Sq)/1e6
        points.append({'model':model,'mu_gev':mu,'mS1Sq_tev2':m1sq,
                       'mS2Sq_tev2':m2sq,
                       'mS1_mass_tev_if_positive':math.sqrt(m1sq) if m1sq>=0 else '',
                       'mS2_mass_tev_if_positive':math.sqrt(m2sq) if m2sq>=0 else ''})
    return points,[hist['cfg_path'],hist['fit_path'],beta_path],int(sol.nfev)


def plot_uv_trajectories(points,model_names,dest):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    models={m:[] for m in model_names}
    for row in points:models[row['model']].append(row)
    f,axes=plt.subplots(4,4,figsize=(15,11),layout='constrained')
    for ax,model in zip(axes.flat,model_names):
        rr=sorted(models[model],key=lambda r:-r['mu_gev'])
        require(len(rr)==80,'incomplete UV trajectory '+model)
        x=[r['mu_gev'] for r in rr]
        ax.plot(x,[float(r['mS1_mass_tev_if_positive']) if r['mS1_mass_tev_if_positive']!='' else float('nan') for r in rr],label=r'$m_{S_1}$')
        ax.plot(x,[float(r['mS2_mass_tev_if_positive']) if r['mS2_mass_tev_if_positive']!='' else float('nan') for r in rr],label=r'$m_{S_2}$')
        ax.axhspan(.5,5,alpha=.12)
        ax.set_xscale('log');ax.set_xlim(1e7,1e5)
        ax.grid(alpha=.15);ax.set_title(model,fontsize=10)
        if ax in axes[:,0]:ax.set_ylabel('mass [TeV]')
        if ax in axes[-1,:]:ax.set_xlabel(r'$\mu$ [GeV]')
    axes[0,0].legend(fontsize=8)
    f.suptitle('Candidate one-loop UV running: original scalar-mass inputs, historical Yukawas\n'
               'Only 10⁷ → 10⁵ GeV; not a physical scalar-threshold matching trajectory',fontsize=13)
    f.savefig(dest/'running_uv_scalar_trajectories.png',dpi=170)
    plt.close(f)


def report_text(frozen,running,model_names,center,all_provenance,*,trajectory_status):
    small=FROZEN_CASES[0];large=FROZEN_CASES[1]
    nacc=lambda case:sum(float(frozen[(case,m)]['chi2_shared_a2'])<=FIT_CUTOFF for m in model_names)
    lines=['# Reduced T3 scalar-sector final numerical summary','',
           '## Scope and accuracy','',
           'Primary data: **two frozen-UV benchmarks**, 16 models each, exact zero additional scalar boundary,\n'
           '18 real fitted Yukawa parameters per model, candidate one-loop scalar EFT1 running,\n'
           'model-specific existing leading fixed-order one-loop C5 matching. Frozen boundary is imposed at 10^5 GeV.',
           '',
           'Reference: **newly optimised scenario-specific A(-2) Yukawa matrices**, transferred unchanged\n'
           'to each model. Other model-specific inputs and group-theory factors are retained.',
           '',
           'UV-running information is **exploratory**: existing high-scale inputs are evolved, not refitted,\n'
           'and mass/hierarchy failures are not interpreted as neutrino-data exclusions.',
           '',
           '## Primary frozen comparison','',
           '| Scenario | Accepted individual fits | Acceptable shared A(−2) evaluations | Total |',
           '|---|---:|---:|---:|',
           f'| Small topology coupling | 16 | {nacc(small)} | 16 |',
           f'| Large topology coupling | 16 | {nacc(large)} | 16 |','',
           'Acceptance is defined by the five-observable Gaussianised NuFIT target χ² ≤ 5.86.\n'
           'The 18 real Yukawa fitting degrees of freedom exceed the five fitted observables: excellent χ²\n'
           'is fit flexibility, not a prediction or a guarantee of uniqueness.',
           '',
           '| Model | small: individual | small: A(−2) shared | large: individual | large: A(−2) shared |',
           '|---|---:|---:|---:|---:|']
    for m in model_names:
        z=[float(frozen[(c,m)][k]) for c in FROZEN_CASES
           for k in ('chi2_optimised','chi2_shared_a2')]
        lines.append('| '+m+' | '+' | '.join(f'{v:.7g}' for v in z)+' |')
    lines.extend(['','## Running-UV diagnostic at fixed inputs','',
        'Historical running-UV Yukawas and scalar input boundaries at 10^7 GeV are left unchanged.\n'
        'The nominal 1 TeV EFT1 endpoint masses are a **formal diagnostic**, not physically matched\n'
        'threshold masses when the scalars become heavy or reverse the decoupling hierarchy.\n'
        'For all 16 historical seeds the original 0.5–5 TeV scalar window fails.',
        '', '| Model | formal mS1 at 1 TeV [TeV] | formal mS2 at 1 TeV [TeV] | rhoF at 100 TeV |',
        '|---|---:|---:|---:|'])
    for m in model_names:
        r=running[m]
        lines.append(f"| {m} | {r['mS1_formal_muS_tev']:.4f} | {r['mS2_formal_muS_tev']:.4f} | {r['rhoF_at_muF']:.4f} |")
    lines.extend(['', f'UV trajectory replay: **{trajectory_status}**.', '',
       '## Interpretation and limitations','',
       '1. Model-specific Yukawa refitting reproduces oscillation data in both frozen scenarios.\n'
       '   The shared A(−2) reference is often not transferable; the comparison does not exclude models.',
       '2. The selected high-scale boundary yields heavy scalar running masses and violates the\n'
       '   *chosen benchmark window*, not a general theoretical requirement.',
       '3. One-loop RG-improved couplings are used with leading fixed-order C5; this is not a uniform\n'
       '   two-loop calculation. Radiatively generated extra B/C quartics are not inserted into C5.',
       '4. Eight additional holomorphic quadratics have zero input boundary and are not evolved;\n'
       '   the zero-boundary relevance argument is conditional, not an exact symmetry.',
       '5. Four extreme-hypercharge models use nonuniform legacy one-scale C5 artifacts.\n'
       '   Full Matchete/RGBeta component normalization and general sequential Wilson matching\n'
       '   have not been certified. No vacuum stability, confidence intervals or predictivity claim.',
       '', '## Provenance','',
       'This report was generated from preexisting local result files; no fit values were modified.'])
    for key,digest in sorted(all_provenance.items()):
        lines.append(f'- `{key}`: `{digest}`')
    return '\n'.join(lines)+'\n'


def run(*,run_id,frozen_run,uv_run,replay_uv=True):
    require(re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]{0,63}',run_id),'invalid run-id')
    require(re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]{0,63}',frozen_run),'invalid frozen-run ID')
    require(re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]{0,63}',uv_run),'invalid uv-run ID')
    dest=OUT/run_id
    require(not dest.exists(),'run-id already exists; reports are not overwritten: '+str(dest))
    frozen_src=INPUT_BASE/frozen_run
    uv_src=INPUT_BASE/uv_run
    payload_path=frozen_src/'comparison.json'
    payload=read_json(payload_path)
    from scripts import RunT3GeneralScalar48Campaign as bulk
    reg=read_json(bulk.REGISTRY)
    model_names=tuple(reg['models'])
    require(len(model_names)==16,'registry no longer has 16 models')
    frozen=validate_frozen_payload(payload,model_names)
    files=validate_frozen_fit_files(frozen_src,frozen)
    running,running_paths=load_running(uv_src,model_names)
    target=read_json(TARGET)
    require(target.get('ordering')=='NO' and all(float(target['one_sigma_errors'][k])>0
            for k in OBS),'unexpected oscillation target')
    source_paths=[payload_path,frozen_src/'comparison.csv',TARGET,bulk.REGISTRY,
                  Path(__file__),*files.values(),*running_paths]
    require(all(p.is_file() for p in source_paths),'missing provenance input')
    dest.mkdir(parents=True,exist_ok=False)
    frozen_rows=[]
    for case in FROZEN_CASES:
        for model in model_names:
            r=frozen[(case,model)]
            frozen_rows.append({'model':model,'case':case,'chi2_individual':r['chi2_optimised'],
                'chi2_shared_A_m2':r['chi2_shared_a2'],
                'individual_accepted':True,
                'shared_A_m2_accepted':float(r['chi2_shared_a2'])<=FIT_CUTOFF,
                'c5_source_kind':r['c5_source_kind'],
                'lambdaT3_input':r['boundary_lambdaT3']})
    write_csv(dest/'frozen_comparison.csv',frozen_rows)
    obs=make_observable_rows(frozen,target['central_values'],target['one_sigma_errors'])
    write_csv(dest/'oscillation_observable_pulls.csv',obs)
    write_csv(dest/'running_uv_formal_endpoints.csv',[running[m] for m in model_names])
    plot_frozen(frozen,model_names,dest)
    plot_running_endpoint(running,model_names,dest)
    tr_status='not requested';tr_nfev={}
    if replay_uv:
        tr_status='candidate renormalisable UV 10^7→10^5 GeV replays completed'
        all_points=[]
        for model in model_names:
            print('REPLAY UV',model,flush=True)
            points,paths,nfev=sample_uv_trajectory(model,uv_src)
            all_points+=points
            tr_nfev[model]=nfev
            source_paths+=paths
        write_csv(dest/'running_uv_trajectories.csv',all_points)
        plot_uv_trajectories(all_points,model_names,dest)
    provenance={str(p.relative_to(ROOT)).replace('\\','/') if p.is_relative_to(ROOT) else str(p):checksum(p)
                for p in source_paths}
    (dest/'summary.md').write_text(report_text(frozen,running,model_names,
        target['central_values'],provenance,trajectory_status=tr_status),encoding='utf-8')
    obj={'format':'t3-reduced-final-report-v1',
         'primary_frozen_run':frozen_run,'diagnostic_running_run':uv_run,
         'number_of_individual_accepted':sum(x['individual_accepted'] for x in frozen_rows),
         'number_of_shared_accepted':sum(x['shared_A_m2_accepted'] for x in frozen_rows),
         'number_of_running_historical_scalar_screen_pass':sum(
             x['screen_pass']=='true' for x in running.values()),
         'candidate_uv_replay_nfev':tr_nfev,'provenance_sha256':provenance,
         'limitations':['not complete general scalar potential',
             'one-loop RG improvement but leading one-loop C5 only',
             '4 legacy one-scale C5 sources are nonuniform',
             'running UV 1 TeV EFT1 endpoints not validated scalar-threshold matching']}
    (dest/'manifest.json').write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('FINAL REPORT',dest,flush=True)
    print('FROZEN',obj['number_of_individual_accepted'],'individual,',
          obj['number_of_shared_accepted'],'shared accepted',flush=True)
    print('RUNNING',obj['number_of_running_historical_scalar_screen_pass'],
          'historical seeds pass original window',flush=True)
    return obj


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-id',default='final_reduced_v1')
    ap.add_argument('--frozen-run',default='reduced_frozen_overnight_v2')
    ap.add_argument('--uv-run',default='reduced_48_v1')
    ap.add_argument('--skip-uv-replay',action='store_true',
                    help='Create tables/plots from saved results only; skip new UV ODE samples')
    a=ap.parse_args()
    run(run_id=a.run_id,frozen_run=a.frozen_run,uv_run=a.uv_run,
        replay_uv=not a.skip_uv_replay)

if __name__=='__main__':main()
