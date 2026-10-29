"""One-pass exploratory T3 general-scalar 16x3 campaign, not sequential EFT.

Bulk native RGBeta general-quartic EFT1 export; use model-specific Stage30/31
UV equations and existing Stage36 hierarchical/legacy C5. Preserve the 8
previously accepted B(-1), C(-1), B(+1) benchmarks by reference. Each other
model/case runs independently; an error never stops the campaign.

    python -m scripts.RunT3GeneralScalar48Campaign --run-id bulk_v1
    python -m scripts.RunT3GeneralScalar48Campaign --run-id bulk_v1 --resume
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.stats import qmc

from scripts import FitT3BPlus1GeneralThreeStudies as b21
from scripts import CoupledA0UVRunning as s30
from scripts import CoupledUVRunning as s31
from scripts import IntermediateEFTRunning as s34

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
OUT = ROOT / 'output/studies/general_scalar_extended/bulk_48'
BETA_DIR = ROOT / 'output/general_t3_stage22'
WL = ROOT / 'scripts/ExportT3GeneralScalarEFT1Bulk.wl'
STUDIES = b21.CASES
LOOP = 16.0 * math.pi**2
EXISTING = ('B_-1', 'C_-1', 'B_+1')
REUSED = ('B_-1', 'C_-1')
NAMES = b21.NAMES
FORMAT = 't3-48-model-general-quartic-conditional-campaign-v1'


def fail(msg):
    raise ValueError('T3 48-study campaign: ' + str(msg))


def read(path):
    p = Path(path)
    if not p.is_file(): fail('Missing '+str(p))
    obj = json.loads(p.read_text(encoding='utf-8-sig'))
    if not isinstance(obj, dict): fail('Invalid JSON object '+str(p))
    return obj


def sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def safe(model):
    return model.replace('+','p').replace('-','m')


def stem(entry):
    a = int(entry['alpha']); d1,d2,df=entry['dimensions']
    return f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{"p" if a>=0 else "m"}{abs(a)}.json'


def inventory(registry, model):
    entry=registry['models'][model]
    new=tuple(x['name'] for x in entry['new_rgbeta_symbols'])
    reused=tuple(x['existing_name'] for x in entry['preexisting_rgbeta_symbols'])
    if entry['total_additional_quartics'] != len(new)+len(reused):
        fail(model+' registry inventory inconsistent')
    return entry,new,reused


def historical(case, model, entry):
    sub,uv_on,nominal_lam=STUDIES[case]
    lam=scenario_lambda_t3(entry,nominal_lam)
    fn=stem(entry)
    cp=ROOT/sub/'fitted'/fn; fp=ROOT/sub/'fits'/fn
    cfg=read(cp);fit=read(fp)
    d1,d2,df=entry['dimensions']
    rep={'d_s1':d1,'d_s2':d2,'d_f':df,'alpha':entry['alpha'],'shared_scalar':False}
    if cfg.get('representation') != rep: fail('Historical representation differs: '+model)
    if cfg.get('base_state',{}).get('mu_gev') != 1e7: fail('Historical seed scale changed')
    if cfg.get('scales') != {'mu_fermion_threshold_gev':1e5,
         'mu_scalar_threshold_gev':1e3,'mu_low_gev':1e2}: fail('Historical threshold scales changed')
    ordinary=cfg['base_state']['ordinary']
    saved=b21.scalar_complex(ordinary['lambdaT3'])
    if abs(saved-imaginary_zero(lam)) > max(1e-20,abs(lam)*1e-12):
        fail('Historical lambdaT3 not fixed at original benchmark')
    if fit.get('model') != fn[:-5] or not math.isclose(float(fit['lambdaT3']),lam,rel_tol=1e-12):
        fail('Historical fitted record differs from model/lambda')
    if uv_on and 'fermion_start_study' in cfg:
        fail('UV-running input actually has frozen-UV flag')
    if not uv_on:
        meta=cfg.get('fermion_start_study',{})
        if ('UV RGEs disabled' not in meta.get('boundary_definition','') or
            meta.get('case') != ('0p5' if nominal_lam==.5 else '5e-9')):
            fail('Frozen-UV historical boundary definition changed')
    from Numerical.orchestration import RunSmallLambdaRefit as r
    y1,y2=[r.matr(ordinary[k]) for k in ('y1','y2')]
    if y1.shape != (3,3) or y2.shape!=(3,3) or max(np.max(abs(y1.imag)),np.max(abs(y2.imag)))>1e-12:
        fail('Expected 18 real historical Yukawas')
    seed=np.r_[y1.real.ravel(),y2.real.ravel()]
    if not np.all(np.isfinite(seed)): fail('Nonfinite seed')
    return {'cfg':cfg,'fit':fit,'cfg_path':cp,'fit_path':fp,'seed':seed,
            'lambda_fixed':lam,'uv_on':uv_on}


def scenario_lambda_t3(entry, nominal_lambda):
    """Exact original scenario boundary, including the E-class half-convention.

    Historical E-class fits use 2.5e-9 (small) and 0.25 (large), not
    5e-9 / 0.5. These are *not* revised inputs; they are the saved originals.
    """
    dimensions = tuple(entry['dimensions'])
    if dimensions not in {(1,3,2),(2,2,1),(2,2,3),(3,1,2),(3,3,2)}:
        fail('Unexpected T3 representation when selecting lambdaT3')
    if nominal_lambda not in (5e-9,0.5):
        fail('Unexpected historical scenario lambdaT3 label')
    return float(nominal_lambda) * (0.5 if dimensions==(3,3,2) else 1.0)


def imaginary_zero(x):
    return complex(float(x),0.)


def update_yukawas(hist,x):
    from Numerical.orchestration import RunSmallLambdaRefit as r
    cfg=deepcopy(hist['cfg']);x=np.asarray(x,dtype=float)
    if x.shape!=(18,) or not np.all(np.isfinite(x)):fail('Invalid fit vector')
    for i,k in enumerate(('y1','y2')):
        r.setmat(cfg['base_state']['ordinary'][k],x[9*i:9*i+9].reshape(3,3).astype(complex))
    if b21.scalar_complex(cfg['base_state']['ordinary']['lambdaT3']) != complex(hist['lambda_fixed']):
        fail('Original lambdaT3 boundary changed')
    return cfg


def validate_eft_betas(payload, state, model, new, reused):
    from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state
    vec,layout=pack_intermediate_scalar_state(state)
    rep=state.representation
    expected={'SharedScalar':False,'dS1':rep.d_s1,'dS2':rep.d_s2,
              'dF':rep.d_f,'alpha':rep.alpha,'IntegratedField':'F',
              'ActiveBSMFields':['S1','S2'],
              'ReportBetaConvention':s30.BETA_CONVENTION}
    if (payload.get('format') != 't3-bulk-general-scalar-eft1-rgbeta-v1' or
        payload.get('status') != 'Success' or payload.get('model') != model or
        payload.get('physical_coupling_mapping_certified') is not False or
        payload.get('wilson_matching_performed') is not False or
        payload.get('general_quadratic_terms_evolved') is not False or
        not isinstance(payload.get('metadata'),dict) or
        any(payload['metadata'].get(k)!=v for k,v in expected.items())):
        fail(model+': native EFT1 beta metadata mismatch')
    if (set(payload.get('base_beta_symbols',[]))!=set(layout.names) or
        len(payload.get('base_beta_symbols',[]))!=len(layout.names) or
        tuple(payload.get('new_quartic_symbols',[]))!=tuple(new) or
        set(payload.get('reused_quartic_symbols',[]))!=set(reused)):
        fail(model+': native EFT1 beta coordinates differ from numerical state')
    betas=payload.get('report_betas')
    if (not isinstance(betas,dict) or set(betas)!=set(layout.names)|set(new) or
        any(not isinstance(v,str) or not v.strip() or any(t in v for t in
            ('BetaTerm[','Missing[','$Failed','$Aborted','Indeterminate','ComplexInfinity'))
            for v in betas.values()) or not np.isfinite(vec).all()):
        fail(model+': native EFT1 equations incomplete')
    for n in new:
        if n not in betas[n]: fail(model+': new quartic has no self beta: '+n)
        if not any(n in betas[k] for k in layout.names):
            fail(model+': new quartic has no backreaction on native scalar/SM beta: '+n)
    return betas,layout


def export_betas(model, entry, new, reused, output_dir, wolframscript, timeout, resume):
    path=output_dir/'betas'/f'{safe(model)}.json'
    log=output_dir/'logs'/f'{safe(model)}_rgbeta.log'
    path.parent.mkdir(parents=True,exist_ok=True)
    log.parent.mkdir(parents=True,exist_ok=True)
    if not (resume and path.is_file()):
        # Reuse the locally verified Batch20 B(+1) 18+1 native export.
        # This avoids re-running a Wolfram calculation that already passed.
        proven=ROOT/'output/studies/general_scalar_extended/b_p1_eft1_running/b_p1_eft1_general_v1/native_eft1_betas.json'
        if model=='B_+1' and new==('lambdaT3GeneralQ022000I1',) and not reused and proven.is_file():
            from scripts import RunT3BPlus1GeneralEFT1 as b20
            baseline=read(proven)
            if (baseline.get('format') != b20.FORMAT or
                baseline.get('status')!='Success' or
                baseline.get('model')!='B_+1' or
                baseline.get('physical_coupling_mapping_certified') is not False):
                fail('Previously successful Batch20 B(+1) native betas are stale')
            wrapped={'format':'t3-bulk-general-scalar-eft1-rgbeta-v1',
                'status':'Success','model':model,'metadata':baseline['metadata'],
                'base_beta_symbols':baseline['base_beta_symbols'],
                'new_quartic_symbols':baseline['new_quartic_symbols'],
                'reused_quartic_symbols':[],
                'physical_coupling_mapping_certified':False,
                'general_quadratic_terms_evolved':False,
                'wilson_matching_performed':False,
                'source_batch20_sha256':sha(proven),
                'report_betas':baseline['report_betas']}
            path.write_text(json.dumps(wrapped,indent=2)+'\n',encoding='utf-8')
            print('REUSE B(+1): native Batch20 18+1 RGBeta equations',flush=True)
        elif not new and not reused:
            # These four models genuinely have no additional general quartics.
            # Stage34 already contains their complete renormalisable scalar
            # quartic system; avoid launching Wolfram merely to re-export it.
            source=BETA_DIR/f'stage34_{safe(model)}_eft1_full_betas.json'
            baseline=read(source)
            if (baseline.get('format') != s34.BETA_FORMAT or
                baseline.get('status') != 'Success' or
                baseline.get('new_stage29_quartics_integrated') is not False or
                baseline.get('physical_coupling_mapping_certified') is not False):
                fail(model+': saved zero-extra Stage34 export invalid')
            wrapped={'format':'t3-bulk-general-scalar-eft1-rgbeta-v1',
                'status':'Success','model':model,'metadata':baseline['metadata'],
                'base_beta_symbols':baseline['beta_symbols'],
                'new_quartic_symbols':[],'reused_quartic_symbols':[],
                'physical_coupling_mapping_certified':False,
                'general_quadratic_terms_evolved':False,
                'wilson_matching_performed':False,
                'source_stage34_sha256':sha(source),
                'report_betas':baseline['report_betas']}
            path.write_text(json.dumps(wrapped,indent=2)+'\n',encoding='utf-8')
            print(f'REUSE {model}: 0 new quartics, saved native Stage34 betas',flush=True)
        else:
            print(f'EXPORT {model}: {len(new)} new + {len(reused)} reused quartics',flush=True)
            result=subprocess.run([wolframscript,'-file',str(WL),model,str(path)],
                cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                timeout=timeout,check=False)
            log.write_text(result.stdout,encoding='utf-8')
            if result.returncode or not path.is_file():
                fail(f'{model} native RGBeta export failed (exit {result.returncode}); see {log}')
    obj=read(path)
    if obj.get('model')!=model or obj.get('format')!='t3-bulk-general-scalar-eft1-rgbeta-v1':
        fail(model+': stale native beta file; delete it and rerun')
    return obj,path


def uv_base(cfg,entry):
    from Numerical.core.StateVector import unpack_uv_state
    vec,layout=s31.prepare_state(cfg,entry,(),{})
    return vec,layout,unpack_uv_state


def uv_betas(model,entry,new,layout):
    if model=='A_+0':
        p=BETA_DIR/'stage30_A_p0_full_uv_betas.json'
        obj=read(p); betas=s30.validate_export(obj,new,layout.names)
    else:
        p=BETA_DIR/f'stage31_{safe(model)}_full_uv_betas.json'
        obj=read(p); betas=s31.validate_export(obj,entry,new,(),layout.names)
    return betas,p


def run_uv(vec,layout,betas,names,initial_extra):
    from Numerical.core.StateVector import unpack_uv_state
    q=np.asarray([x for z in initial_extra for x in (z.real,z.imag)],dtype=float)
    initial=np.r_[vec,q]
    rhs=lambda t,y:s30.combined_rhs(t,y,layout=layout,names=names,betas=betas)
    sol=solve_ivp(rhs,(math.log(1e7),math.log(1e5)),initial,
        method='DOP853',rtol=1e-8,atol=1e-11)
    if not sol.success or not np.isfinite(sol.y).all():fail('Full coupled UV integration failed: '+sol.message)
    end=sol.y[:,-1]
    return unpack_uv_state(end[:layout.size],layout,mu_gev=1e5),tuple(
        complex(end[layout.size+2*i],end[layout.size+2*i+1]) for i in range(len(names))),int(sol.nfev)


def eft_rhs(t,y,layout,betas,names):
    from Numerical.running.IntermediateScalarStateVector import unpack_intermediate_scalar_state
    from Numerical.running.IntermediateScalarRGBetaEvaluator import intermediate_scalar_state_environment
    from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
    from Numerical.core.BetaVector import beta_values_to_derivative
    state=unpack_intermediate_scalar_state(y[:layout.size],layout,mu_gev=math.exp(t))
    env=intermediate_scalar_state_environment(state)
    env.update({n:complex(y[layout.size+2*i],y[layout.size+2*i+1])
                for i,n in enumerate(names)})
    values={n:evaluate_inputform_expression(expr,env) for n,expr in betas.items()}
    base=beta_values_to_derivative({n:values[n] for n in layout.names},layout,
                                   divide_by_loop_factor=True)
    extra=[]
    for n in names:
        if isinstance(values[n],np.ndarray):fail('Unexpected matrix beta for '+n)
        z=complex(values[n]); extra.extend((z.real/LOOP,z.imag/LOOP))
    derivative=np.r_[base,extra]
    if not np.isfinite(derivative).all():fail('Nonfinite coupled EFT1 beta')
    return derivative


def run_eft(uv,extra,betas,model,new,reused):
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from Numerical.running.IntermediateScalarStateVector import (
        pack_intermediate_scalar_state,unpack_intermediate_scalar_state)
    from Numerical.running.ScalarThresholdBoundary import scalar_threshold_masses
    initial_state=build_intermediate_scalar_boundary(uv,matching_scale_gev=1e5)
    beta,layout=validate_eft_betas(betas,initial_state,model,new,reused)
    vec,_=pack_intermediate_scalar_state(initial_state)
    y0=np.r_[vec,[x for z in extra for x in (z.real,z.imag)]]
    sol=solve_ivp(lambda t,y:eft_rhs(t,y,layout,beta,new),
       (math.log(1e5),math.log(1e3)),y0,method='DOP853',rtol=1e-8,atol=1e-11)
    if not sol.success or not np.isfinite(sol.y).all():fail('Coupled EFT1 integration failed: '+sol.message)
    last=sol.y[:,-1]
    state=unpack_intermediate_scalar_state(last[:layout.size],layout,mu_gev=1e3)
    masses=np.asarray(scalar_threshold_masses(state).masses_gev,dtype=float)
    values=tuple(complex(last[layout.size+2*i],last[layout.size+2*i+1])
                 for i in range(len(new)))
    return state,masses,values,int(sol.nfev)


class LegacyKernel(b21.MatchedKernel):
    """Explicit older one-scale approximation for four extreme hypercharges."""
    def __init__(self,path):
        import sympy as sp
        from RGE.matching.MatchedWeinbergRGE import parse_matchete_c5
        from RGE.matching.WeinbergFlavorMatching import extract_one_generation_loop_kernel
        expr=extract_one_generation_loop_kernel(parse_matchete_c5(Path(path).read_text(encoding='utf-8')))
        self.names=('MF','MS1','MS2','MS','lambdaT3','hbar')
        syms=[sp.Symbol(n) for n in self.names]
        if not expr.free_symbols<=set(syms):fail('Legacy C5 has unresolved non-standard symbols')
        self.fun=sp.lambdify(syms,expr,modules='numpy',cse=True)
        self.path=Path(path)

    def parity_check(self,uv,eft):
        from Numerical.running.FinalC5TrajectoryAdapter import extract_final_c5_inputs
        from Numerical.running.SMWeinbergStage import evaluate_symbolic_c5
        view=SimpleNamespace(uv_threshold_state=uv,eft1_threshold_state=eft,
                              mu_scalar_threshold_gev=1e3)
        exact=evaluate_symbolic_c5(self.path,extract_final_c5_inputs(view).as_config())
        rapid=self.evaluate(uv,eft)
        err=float(np.max(abs(exact-rapid))/max(float(np.max(abs(exact))),1e-30))
        if err>3e-8:fail('Legacy fast C5 numerical parity failed')
        return err


def matching(model):
    group,alpha=model.split('_')
    token=('m' if alpha.startswith('-') else 'p')+alpha[1:]
    p=ROOT/f'output/full/hypercharge/T3_{group}_alpha_{token}/data/final_weinberg_coefficient.json'
    if p.is_file():return b21.MatchedKernel(p),p,'hierarchical-hard-plus-direct'
    legacy=ROOT/f'output/hypercharge/T3_{group}_alpha_{token}/data/c5_coefficient.txt'
    if legacy.is_file():return LegacyKernel(legacy),legacy,'legacy-common-scale'
    fail(model+': neither own hierarchical nor legacy matching artifact found')


def testscreen(uv,eft,masses,qf,qs):
    from Numerical.running.FermionThresholdBoundary import fermion_threshold_masses
    from Numerical.running.IntermediateScalarRGBetaEvaluator import intermediate_scalar_state_environment
    f=np.asarray(fermion_threshold_masses(uv).masses_gev,dtype=float)
    sq=np.asarray([uv.mS1Sq,uv.mS2Sq],dtype=float)
    if np.min(sq)<=0 or np.min(f)<=0 or np.min(masses)<=0:
        return {'screen':False,'reason':'Nonpositive running scalar or fermion mass squared'}
    rho=float(np.sqrt(sq.max())/f.min())
    y=max(np.linalg.svd(np.asarray(z,dtype=complex),compute_uv=False)[0] for z in (uv.y1,uv.y2))
    env=intermediate_scalar_state_environment(eft)
    allq=[abs(complex(v)) for k,v in env.items() if k.startswith('lambda')]
    allq += [abs(z) for z in (*qf,*qs)]
    qm=float(max(allq))
    screen=bool(rho<=.5 and np.min(masses)>=500. and np.max(masses)<=5000.
                and y<=math.sqrt(4*math.pi) and qm<=4*math.pi)
    return {'screen':screen,'rhoF':rho,'scalar_masses_muS_tev':np.sort(masses/1e3).tolist(),
      'max_yukawa_singular_muF':float(y),'max_scalar_quartic_abs':qm,
      'reason':('AcceptedScreen' if screen else 'MassHierarchyWindowOrPerturbativity')}


class Evaluator:
    def __init__(self,case,model,entry,new,history,payload,kernel,target,extra_real):
        self.case=case;self.model=model;self.entry=entry;self.new=new
        self.hist=history;self.payload=payload;self.kernel=kernel;self.target=target
        self.q0=tuple(complex(extra_real) for _ in new)
        self.count=0;self.mass_safe=0;self.failures=0;self.best=None;self.parity=None
        self.problems=[]
        self.uvspec=None
        if history['uv_on']:
            vec,layout,_=uv_base(history['cfg'],entry)
            betas,path=uv_betas(model,entry,new,layout)
            self.uvspec=(betas,path)
        else:
            vec,layout,unpack=uv_base(history['cfg'],entry)
            initial_uv=unpack(vec,layout,mu_gev=1e5)
            self.frozen=run_eft(initial_uv,self.q0,payload,model,new,())

    def evaluate(self,x,label):
        self.count+=1
        row={'label':label,'chi2':None,'status':'Rejected','reason':'','prediction':None}
        try:
            cfg=update_yukawas(self.hist,x)
            vec,layout,unpack=uv_base(cfg,self.entry)
            if self.hist['uv_on']:
                uv,qf,nuv=run_uv(vec,layout,self.uvspec[0],self.new,self.q0)
                eft,masses,qs,neft=run_eft(uv,qf,self.payload,self.model,self.new,())
            else:
                uv=unpack(vec,layout,mu_gev=1e5)
                eft,masses,qs,neft=self.frozen
                qf=self.q0;nuv=0
            chk=testscreen(uv,eft,masses,qf,qs)
            row.update(chk)
            if chk['screen']:
                if self.parity is None:self.parity=self.kernel.parity_check(uv,eft)
                c5=self.kernel.evaluate(uv,eft)
                from scripts.RunT3GeneralScalarBCExploratoryFit import low_energy_fit
                fit,_=low_energy_fit(c5,eft,self.target,
                   vev=float(cfg.get('numerical',{}).get('vev_gev',246.22)))
                res=b21.residual_from_fit(fit)
                row.update({'status':'NeutrinoEvaluated','chi2':float(fit.chi2),
                  'prediction':dict(zip(fit.observable_names,map(float,fit.prediction))),
                  'whitened_residual':res.tolist(),
                  'parameters':dict(zip(NAMES,map(float,x))),
                  'uv_nfev':nuv,'eft1_nfev':neft})
                self.mass_safe+=1
                if self.best is None or row['chi2']<self.best['chi2']:self.best=deepcopy(row)
            if len(self.problems)<4 and row['chi2'] is None:self.problems.append({'label':label,'reason':row['reason']})
        except (ValueError,TypeError,RuntimeError,ArithmeticError,OverflowError,
                FloatingPointError,np.linalg.LinAlgError) as exc:
            self.failures+=1
            row['status']='NumericalFailure';row['reason']=str(exc)[:400]
            if len(self.problems)<4:self.problems.append({'label':label,'reason':row['reason']})
        return row


def run_case(model,case,registry,payload,matching_kernel,match_path,match_kind,
             beta_path,dest,*,sobol_points,max_nfev,random_seed,extra_real):
    entry,new,reused=inventory(registry,model)
    if reused:fail('Reused quartic special models must use preexisting verified fits')
    history=historical(case,model,entry)
    from Numerical.orchestration import RunSmallLambdaRefit as refit
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    target=OscillationFitTarget.from_json(refit.TARGET)
    if target.ordering!='NO' or len(target.observable_names)!=5:fail('Five-observable NO target changed')
    evaluator=Evaluator(case,model,entry,new,history,payload,matching_kernel,target,extra_real)
    lo,hi=b21.search_box(history['seed'],local_fraction=1.5,
                  floor_fraction=.05,abs_yukawa_max=math.sqrt(4*math.pi))
    trials=[(history['seed'],'historical_seed')]
    draws=qmc.Sobol(d=18,scramble=True,seed=random_seed).random_base2(int(math.log2(sobol_points)))
    trials += [(lo+(hi-lo)*z,f'sobol_{i+1}') for i,z in enumerate(draws)]
    accepted=[]
    for x,label in trials:
        row=evaluator.evaluate(x,label)
        if row['chi2'] is not None:accepted.append(row)
    if accepted and min(z['chi2'] for z in accepted)>5.86 and max_nfev:
        best=min(accepted,key=lambda v:v['chi2'])
        x0=np.array([best['parameters'][k] for k in NAMES]);z0=np.clip((x0-lo)/(hi-lo),1e-7,1-1e-7)
        def residual(z):
            row=evaluator.evaluate(lo+(hi-lo)*z,'local')
            return np.asarray(row.get('whitened_residual',[1e4]*5),dtype=float)
        least_squares(residual,z0,bounds=(np.zeros(18),np.ones(18)),
           max_nfev=max_nfev,xtol=1e-6,ftol=1e-6,gtol=1e-6)
    best=evaluator.best
    status=('CandidateConstrainedAccepted' if best and best['chi2']<=5.86
            else 'CandidateConstrainedPoorFit' if best else 'NoMassSafeNeutrinoFit')
    report={'format':FORMAT,'model':model,'case':case,'status':status,
      'historical_lambdaT3_fixed':history['lambda_fixed'],
      'extra_quartics_rgbeta_boundary_real':extra_real,
      'new_quartic_symbols':list(new),
      'general_quadratics_included':False,
      'c5_source_kind':match_kind,'c5_source':str(match_path.relative_to(ROOT)),
      'physical_matchete_rgbeta_map_certified':False,
      'full_sequential_matching_performed':False,
      'counts':{'evaluations':evaluator.count,'mass_safe':evaluator.mass_safe,
                'numerical_failures':evaluator.failures},
      'matching_parity_error':evaluator.parity,
      'best_mass_safe':best,'first_problems':evaluator.problems,
      'inputs_sha256':{'historical':sha(history['cfg_path']),
        'historical_fit':sha(history['fit_path']),
        'matched_c5':sha(match_path),
        'native_eft1_betas':sha(beta_path)}}
    if evaluator.uvspec:report['inputs_sha256']['uv_betas']=sha(evaluator.uvspec[1])
    dest.mkdir(parents=True,exist_ok=False)
    (dest/'fit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if best:
        cfg=update_yukawas(history,np.array([best['parameters'][k] for k in NAMES]))
        cfg['conditional_general_scalar_fit']={'model':model,'case':case,
            'chi2':best['chi2'],'sequential_matching_certified':False}
        (dest/'fitted_config.json').write_text(json.dumps(cfg,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return report


def validated_bc_fit(case, model):
    """Recover the six independently replay-validated B(-1)/C(-1) fits.

    B/C frozen fits were written as MODEL_exploratory_fit.json, not fit.json.
    The Batch18 validation report records exact paths and file digests.
    Use those rather than guessing layout or rerunning expensive fits.
    """
    if model not in REUSED:
        return None
    validation=ROOT/'output/studies/general_scalar_extended/validation/bc_six_validation_v1/report.json'
    if not validation.is_file():
        return None
    proof=read(validation)
    if (proof.get('format')!='t3-bc-six-fit-conditional-validation-v1' or
        proof.get('status')!='ReplayValidated' or
        proof.get('candidate_count')!=6 or
        proof.get('perturbative_coordinate_screens_all_pass') is not True or
        proof.get('rg_replay_performed') is not True or
        proof.get('matching_physically_certified') is not False or
        proof.get('full_sequential_matching_performed') is not False):
        fail('Saved six-fit replay validation record inconsistent')
    rows=[z for z in proof.get('results',()) if
        z.get('model')==model and z.get('case')==case]
    if len(rows)!=1:
        fail('Saved six-fit validation missing or duplicates '+str((model,case)))
    row=rows[0]
    if (row.get('status')!='SavedFitValidated' or
        row.get('replay_status')!='CandidateReplayAgrees' or
        not all(row.get(flag) is True for flag in (
          'quartic_coordinate_screen_pass','yukawa_singular_screen_pass','positivity_sampled'))):
        fail('Saved B/C fit has not passed independent replay/screening')
    chi2=float(row['chi2'])
    if not (math.isfinite(chi2) and 0<=chi2<=5.86 and
            math.isclose(chi2,float(row['replayed_chi2']),rel_tol=1e-10,abs_tol=1e-12)):
        fail('Saved B/C fit chi2 differs from independent replay')
    def verified_path(raw,digest,label):
        if not isinstance(raw,str) or not raw or not isinstance(digest,str) or len(digest)!=64:
            fail('Invalid B/C '+label+' path/digest')
        # Validation JSON was written on Windows; canonicalise before Path.
        path=(ROOT / raw.replace('\\','/')).resolve()
        if not path.is_relative_to(ROOT.resolve()) or not path.is_file() or sha(path)!=digest:
            fail('B/C '+label+' missing or differs from validated SHA-256: '+str(path))
        return path
    fit_path=verified_path(row['fit_path'],row['fit_sha256'],'fit')
    config_path=verified_path(row['cfg_path'],row['config_sha256'],'config')
    prediction=row.get('prediction')
    keys=('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')
    if not isinstance(prediction,dict) or set(keys)-set(prediction) or not all(
            math.isfinite(float(prediction[k])) for k in keys):
        fail('Saved validated B/C observables incomplete')
    m1=float(row['scalar_light_muS_tev']);m2=float(row['scalar_heavy_muS_tev'])
    rho=float(row['rhoF']);lambda_uv=float(row['lambdaT3_uv'])
    if not (0.5<=m1<=m2<=5.0 and 0<=rho<=0.5 and math.isfinite(lambda_uv)):
        fail('Saved validated B/C mass or UV lambda boundary inconsistent')
    boundary_relaxed=(case=='small_uv_running_5e-9')
    expected=5e-9 if case=='small_uv_frozen_5e-9' else 0.5 if case=='large_uv_frozen_0p5' else None
    if expected is not None and not math.isclose(lambda_uv,expected,rel_tol=1e-10,abs_tol=1e-15):
        fail('Frozen B/C fit does not preserve original lambdaT3 boundary')
    if boundary_relaxed and math.isclose(lambda_uv,5e-9,rel_tol=1e-10,abs_tol=1e-15):
        fail('Full UV B/C historical boundary relaxation label changed')
    return {'model':model,'case':case,'status':'ExistingConditionalAccepted',
        'best_mass_safe':{'chi2':chi2,'rhoF':rho,
           'scalar_masses_muS_tev':[m1,m2],
           'max_yukawa_singular_muF':float(row['max_yukawa_singular']),
           'max_scalar_quartic_abs':float(row['max_abs_quartic']),
           'prediction':{k:float(prediction[k]) for k in keys}},
        'source_fit':str(fit_path.relative_to(ROOT)),
        'source_fit_sha256':sha(fit_path),
        'source_fitted_config':str(config_path.relative_to(ROOT)),
        'source_fitted_config_sha256':sha(config_path),
        'source_validation':str(validation.relative_to(ROOT)),
        'source_validation_sha256':sha(validation),
        'historical_lambdaT3_fixed':lambda_uv,
        'original_boundary_relaxed':boundary_relaxed,
        'c5_source_kind':'conditional-one-scale-mixed-scalar',
        'reason':'Independently validated Batch18; not a new bulk refit'}


def existing_fit(case,model):
    """Reference accepted saved fits without rerunning them."""
    if model in REUSED:
        return validated_bc_fit(case,model)
    base=ROOT/'output/studies/general_scalar_extended'/case
    if not base.is_dir():return None
    accepted=[]
    for path in base.glob('*/fit.json'):
        try:
            doc=read(path);best=doc.get('best_mass_safe')
            if (doc.get('model')==model and doc.get('case')==case and
                doc.get('status')=='CandidateConstrainedAccepted' and
                isinstance(best,dict) and isinstance(best.get('chi2'),(int,float))
                and math.isfinite(best['chi2']) and best['chi2']<=5.86):
                accepted.append((float(best['chi2']),path,doc))
        except (OSError,ValueError,KeyError,TypeError):continue
    if not accepted:return None
    _,p,doc=min(accepted,key=lambda v:v[0])
    return {'model':model,'case':case,'status':'ExistingConditionalAccepted',
            'best_mass_safe':doc['best_mass_safe'],
            'source_fit':str(p.relative_to(ROOT)), 'source_sha256':sha(p),
            'historical_lambdaT3_fixed':doc.get('historical_lambdaT3_fixed'),
            'c5_source_kind': 'hierarchical-hard-plus-direct',
            'original_boundary_relaxed':False}


def row_of(report):
    best=report.get('best_mass_safe') or {}
    prediction=best.get('prediction') or {}
    return {'model':report.get('model'),'study':report.get('case'),
       'status':report.get('status'),'chi2':best.get('chi2'),
       'rhoF':best.get('rhoF'),'scalar_masses_muS_tev':best.get('scalar_masses_muS_tev'),
       'delta_m21_sq_ev2':prediction.get('delta_m21_sq_ev2'),
       'delta_m3l_sq_ev2':prediction.get('delta_m3l_sq_ev2'),
       'sin2_theta12':prediction.get('sin2_theta12'),
       'sin2_theta13':prediction.get('sin2_theta13'),
       'sin2_theta23':prediction.get('sin2_theta23'),
       'mass_safe_evaluations':report.get('counts',{}).get('mass_safe'),
       'lambdaT3_input':report.get('historical_lambdaT3_fixed'),
       'c5_source_kind':report.get('c5_source_kind'),
       'boundary_relaxed':bool(report.get('original_boundary_relaxed',False)),
       'reason':report.get('reason')}


def run(args):
    if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id):fail('Invalid run ID')
    if args.sobol_points<2 or args.sobol_points & (args.sobol_points-1):
        fail('Sobol points must be a power of two >=2')
    if not (0<=args.max_nfev<=1000 and 0<abs(args.extra_real)<=.5):
        fail('Invalid optimisation budget or extra coupling')
    registry=read(REGISTRY)
    if len(registry.get('models',{}))!=16:fail('Expected 16 registered models')
    models=tuple(registry['models'])
    if args.models:
        wanted=args.models.split(',')
        if any(m not in models for m in wanted):fail('Unknown selected model(s)')
        models=tuple(m for m in models if m in wanted)
    base=OUT/args.run_id
    if base.exists() and not args.resume:fail('Run folder already exists; specify --resume or a new --run-id')
    base.mkdir(parents=True,exist_ok=True)
    # errors.log is a permanent audit trail.  This file is only the errors
    # from the current invocation, avoiding obsolete diagnostics in tail output.
    fresh_errors=base/'errors_latest.log'
    fresh_errors.write_text('',encoding='utf-8')
    def record_error(phase, identity, exc):
        line=f'{phase} {identity}: {exc}\n'
        for target in (base/'errors.log', fresh_errors):
            with target.open('a',encoding='utf-8') as stream:
                stream.write(line)
    records={}
    # Pre-export each model once (never run Wolfram inside optimiser).
    for model in models:
        entry,new,reused=inventory(registry,model)
        if model in REUSED:continue
        if model=='B_+1' and args.phase!='export':
            # Export is still necessary for the one blocked B(+1) full-UV case.
            pass
        if args.phase in ('all','export'):
            try:
                export_betas(model,entry,new,reused,base,args.wolframscript,args.timeout_sec,args.resume)
            except Exception as exc:
                record_error('EXPORT',model,exc)
                print('FAILED EXPORT',model,str(exc),flush=True)
    if args.phase=='export':return
    for model in models:
        entry,new,reused=inventory(registry,model)
        if model in REUSED:
            for case in STUDIES:
                record=existing_fit(case,model)
                records[(model,case)]=record or {'model':model,'case':case,
                    'status':'ExistingAcceptedFitNotFound','reason':'Locate previous accepted fit.json'}
            continue
        kernel=None; match_path=None;kind=None;payload=None;beta_path=None
        for case in STUDIES:
            if model=='B_+1' and case!='small_uv_running_5e-9':
                record=existing_fit(case,model)
                if record:
                    records[(model,case)]=record;continue
            dest=base/case/safe(model)
            if dest.joinpath('fit.json').is_file() and args.resume:
                records[(model,case)]=read(dest/'fit.json');continue
            try:
                if reused:fail(model+' uses special B/C mixed-quadratic implementation; no bulk refit')
                if payload is None:
                    beta_path=base/'betas'/f'{safe(model)}.json'
                    payload=read(beta_path)
                if kernel is None:
                    kernel,match_path,kind=matching(model)
                rec=run_case(model,case,registry,payload,kernel,match_path,kind,beta_path,dest,
                    sobol_points=args.sobol_points,max_nfev=args.max_nfev,
                    random_seed=args.seed,extra_real=args.extra_real)
                records[(model,case)]=rec
                print(f'FIT {model} {case}: {rec["status"]} '
                      f'chi2={(rec.get("best_mass_safe") or {}).get("chi2")}',flush=True)
            except Exception as exc:
                rec={'model':model,'case':case,'status':'FailedOrUnavailable',
                     'reason':str(exc)[:500]}
                records[(model,case)]=rec
                record_error('FIT',f'{model} {case}',exc)
                print('FAILED',model,case,str(exc),flush=True)
    # Always write a 48-slot table in full-mode, even if individual models failed.
    for model in registry['models']:
        for case in STUDIES:
            if (model,case) not in records:
                prior=base/case/safe(model)/'fit.json'
                if prior.is_file():records[(model,case)]=read(prior)
                elif (old:=existing_fit(case,model)) is not None:
                    records[(model,case)]=old
                else:
                    records[(model,case)]={'model':model,'case':case,
                         'status':'NotRunInThisCampaign','reason':'Model not selected or input unavailable'}
    rows=[row_of(records[(m,c)]) for m in registry['models'] for c in STUDIES]
    (base/'campaign.json').write_text(json.dumps({'format':FORMAT,'rows':rows,
       'model_count':16,'study_count':3,'conditional_only':True,
       'general_quadratic_extensions_not_in_bulk_runner':True,
       'no_sequential_wilson_matching':True},indent=2,allow_nan=False)+'\n',encoding='utf-8')
    columns=list(rows[0]);with_file=base/'campaign.csv'
    with with_file.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)
    accepted=sum(r['status'] in ('CandidateConstrainedAccepted','ExistingConditionalAccepted') for r in rows)
    relaxed=sum(r['status'] in ('CandidateConstrainedAccepted','ExistingConditionalAccepted')
                and r['boundary_relaxed'] for r in rows)
    lines=['# General-scalar exploratory 16 × 3 campaign','',
     f'Rows: 48; accepted conditional: {accepted}; original-boundary accepted: {accepted-relaxed}; '
     f'modified-boundary accepted: {relaxed}.',
     'E-class frozen benchmarks retain the historical half-normalisation: 2.5e-9 / 0.25.',
     'Model-specific native extra-quartic renormalisable UV/EFT1 RGEs used where exported.',
     'Not full general scalar potential: holomorphic quadratic extensions and Wilson EFT matching remain excluded.',
     'Existing B(-1)/C(-1) special fits referenced without rerun; legacy C5 only for four extreme hypercharges.',
     'errors_latest.log contains only errors from the most recent invocation; errors.log is historical.',
     '', '| Model | Frozen small | Frozen 0.5 | Full UV small |','|---|---|---|---|']
    for m in registry['models']:
        cols=[]
        for c in STUDIES:
            r=row_of(records[(m,c)])
            v=r['chi2'];tag=r['status'];cols.append(f'{v:.5g} ({tag})' if v is not None else tag)
        lines.append('| '+m+' | '+' | '.join(cols)+' |')
    (base/'summary.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(f'CAMPAIGN SUMMARY: accepted {accepted}/48; {base/"campaign.csv"}',flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-id',required=True)
    ap.add_argument('--resume',action='store_true')
    ap.add_argument('--phase',choices=('all','export','fit'),default='all')
    ap.add_argument('--models',help='Optional comma-separated model subset, all by default')
    ap.add_argument('--sobol-points',type=int,default=8)
    ap.add_argument('--max-nfev',type=int,default=25)
    ap.add_argument('--seed',type=int,default=20261028)
    ap.add_argument('--extra-real',type=float,default=-.01)
    ap.add_argument('--timeout-sec',type=int,default=3600)
    ap.add_argument('--wolframscript',default='wolframscript')
    args=ap.parse_args()
    run(args)

if __name__=='__main__':main()
