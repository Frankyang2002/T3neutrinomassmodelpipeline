"""Conditional 16-model Yukawa refit of SAVED frozen-UV full-scalar EFT1 backgrounds.

Does not vary extra quartics, diagonal masses, holomorphic quadratic masses, MF,
SM couplings or lambdaT3. No UV 10 PeV->100 TeV integration. This is a FIXED-BACKGROUND
fit: Yukawa-dependent EFT1 scalar RGE evolution is not recalculated per trial. Uses the same
conditional one-scale C5 and SM+Weinberg code as RunT3FullScalarNeutrino.
Self-bilinear C5 is omitted for six models; the B/C(-1) Matchete/RGBeta
mapping and full sequential matching are not certified.

Run from repository root:
 python -m scripts.RefitT3FrozenFullScalarLow --plan
 python -u -m scripts.RefitT3FrozenFullScalarLow --run-id frozen_full_scalar_refit_low_v1
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from scripts import RunT3FullScalarNeutrino as neut
from scripts import RunT3FullScalarContinuous as full
from scripts import T3Frozen64Sources as policy

ROOT = full.ROOT
INPUT = ROOT / 'output/studies/general_scalar_extended/full_scalar_frozen_64/frozen_general64_v1'
OUTPUT = ROOT / 'output/studies/general_scalar_extended/full_scalar_frozen_refits'
NAMES = ('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')
PARAMS = tuple(f'{field}_{i}{j}' for field in ('y1','y2') for i in range(3) for j in range(3))
FORMAT = 't3-frozen-full-scalar-conditional-yukawa-refit-v1'


def require(ok, why):
    if not ok:
        raise ValueError('Frozen full-scalar refit: '+why)


def hashfile(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def seed_from_uv(uv):
    y1, y2 = np.asarray(uv.y1),np.asarray(uv.y2)
    require(y1.shape == (3,3) and y2.shape == (3,3),'unexpected Yukawa shape')
    require(max(np.max(abs(y1.imag)), np.max(abs(y2.imag))) < 1e-11,
            'historical seed has complex Yukawas; this fitter uses 18 real entries')
    return np.concatenate((y1.real.ravel(),y2.real.ravel())).astype(float)


def uv_with_yukawas(uv, x):
    x = np.asarray(x,dtype=float)
    require(x.shape == (18,) and np.isfinite(x).all(),'invalid Yukawa fit vector')
    return replace(uv,y1=x[:9].reshape(3,3).astype(complex),
                   y2=x[9:].reshape(3,3).astype(complex)).validated()


def whiten(pred,target):
    v = np.array([float(pred[n]) for n in target.observable_names])
    require(np.isfinite(v).all(),'nonfinite oscillation observables')
    return np.linalg.solve(np.linalg.cholesky(target.covariance),v-target.central_values)


def load_case(source_root, model):
    tag=full.model_tag(model)
    folder=source_root/tag
    traj=folder/'continuous_trajectory.json'
    data=neut.verify_continuous(model,traj)
    require(data.get('boundary_mode')=='FROZEN_UV','expected frozen UV benchmark')
    require(data.get('yukawa_mode')=='optimised','expected original optimised Yukawa seed')
    uv=neut.decode_uv(data['uv_fermion_threshold_state'])
    eft=neut.decode_eft(data['eft1_scalar_endpoint_state'])
    endpoint=folder/data['saved_candidate_eft1_endpoint']
    require(endpoint.is_file(),'missing scalar-mixing endpoint sidecar')
    hashes={str(traj.relative_to(ROOT)):hashfile(traj),
            str(endpoint.relative_to(ROOT)):hashfile(endpoint)}
    return uv,eft,folder,hashes


class Evaluator:
    def __init__(self,model,uv,eft,folder,target):
        self.model,self.uv,self.eft,self.folder,self.target=model,uv,eft,folder,target
        self.calls=0
        self.failures=0
        self.last_error=''
        self.best=None
        self.matching_source=None
        self.matching_kind=None

    def predict(self,x):
        trial=uv_with_yukawas(self.uv,x)
        if self.model in neut.MIXED:
            c5,source,kind,_=neut.matched_bc_c5(self.model,trial,self.folder)
        else:
            c5,source,kind=neut.matched_nonmixed_c5(self.model,trial,self.eft)
        self.matching_source=str(source.relative_to(ROOT))
        self.matching_kind=kind
        p=neut.low_energy(self.model,c5,self.eft)['oscillation_observables_normal_ordering']
        return p

    def evaluate(self,x):
        self.calls+=1
        try:
            p=self.predict(x)
            resid=whiten(p,self.target)
            chi=float(resid@resid)
            if math.isfinite(chi) and (self.best is None or chi<self.best['chi2']):
                self.best={'chi2':chi,'x':list(map(float,x)),'prediction':p,
                           'whitened_residual':resid.tolist()}
            return resid
        except (ValueError,TypeError,RuntimeError,ArithmeticError,OverflowError,
                FloatingPointError,np.linalg.LinAlgError,OSError) as exc:
            self.failures+=1
            self.last_error=(type(exc).__name__+': '+str(exc))[:400]
            return np.full(len(self.target.observable_names),1e6,dtype=float)


def trial_seeds(seed,original_pred,experimental, bound, count, rng):
    """Neutrino mass scales quadratically with Y, hence dm^2 scales as Y^4."""
    factor=1.0
    if original_pred and original_pred.get('delta_m3l_sq_ev2',0)>0:
        factor=(float(experimental['delta_m3l_sq_ev2'])/
                float(original_pred['delta_m3l_sq_ev2']))**0.25
    factor=float(np.clip(factor,1e-4,1e4))
    candidates=[seed, seed*factor]
    for _ in range(max(0,count-2)):
        candidates.append(seed*factor*(1+rng.normal(0,.20,18)))
    return [np.clip(p,-bound*.999,bound*.999) for p in candidates]


def run_one(model,source_root,target,args):
    uv,eft,folder,hashes=load_case(source_root,model)
    seed=seed_from_uv(uv)
    require(np.max(abs(seed))<args.bound,'historical seed exceeds --bound')
    evaluator=Evaluator(model,uv,eft,folder,target)
    old_source=INPUT/'optimised_low'/'neutrino_v1'/(full.model_tag(model)+'_neutrino.json')
    old=neut.load(old_source) if old_source.is_file() else {}
    old_pred=old.get('oscillation_observables_normal_ordering')
    target_pred=dict(zip(target.observable_names,map(float,target.central_values)))
    rng=np.random.default_rng(args.seed+policy.MODELS.index(model))
    seeds=trial_seeds(seed,old_pred,target_pred,args.bound,args.starts,rng)
    for x in seeds: evaluator.evaluate(x)
    require(evaluator.best is not None,'all initial seeds failed: '+evaluator.last_error)
    starting_best=dict(evaluator.best)
    if starting_best['chi2']>args.target and args.max_nfev:
        for x in seeds:
            if evaluator.best and evaluator.best['chi2']<=args.target:break
            try:
                least_squares(evaluator.evaluate,x,bounds=(-args.bound,args.bound),
                              max_nfev=args.max_nfev,xtol=1e-6,ftol=1e-6,gtol=1e-6,
                              diff_step=args.diff_step)
            except (ValueError,RuntimeError,OverflowError,FloatingPointError) as exc:
                evaluator.failures+=1
                evaluator.last_error=str(exc)[:400]
    best=evaluator.best
    require(best is not None,'no physical prediction')
    return {'format':FORMAT,'model':model,'benchmark':'optimised_low',
            'status':'ConditionalFitAccepted' if best['chi2']<=args.target else 'ConditionalFitPoor',
            'chi2':best['chi2'],'starting_seed_chi2':starting_best['chi2'],
            'prediction':best['prediction'], 'whitened_residual':best['whitened_residual'],
            'parameters':dict(zip(PARAMS,best['x'])),
            'source_sha256':hashes,'target_sha256':hashfile(args.target_file),
            'matching_source':evaluator.matching_source,'matching_kind':evaluator.matching_kind,
            'evaluations':evaluator.calls,'failed_evaluations':evaluator.failures,
            'last_error':evaluator.last_error,
            'frozen_uv':True,'scalar_background_held_fixed_during_fit':True,
            'yukawa_dependent_scalar_backreaction_recomputed':False,
            'scalar_boundary_refitted':False,'yukawas_refitted':True,
            'complex_yukawa_phases_refitted':False,
            'physical_matchete_rgbeta_mapping_certified':False,
            'full_sequential_matching_performed':False,
            'additional_self_bilinear_c5_omitted':model in neut.SELF,
            'mass_eigenstate_one_scale_matching_only':model in neut.MIXED}


def save(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(path)


def main(argv=None):
    from Numerical.orchestration import RunSmallLambdaRefit as historical
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-root',type=Path,default=INPUT / 'optimised_low')
    p.add_argument('--target-file',type=Path,default=historical.TARGET)
    p.add_argument('--run-id',default='frozen_full_scalar_refit_low_v1')
    p.add_argument('--models',default='all')
    p.add_argument('--max-nfev',type=int,default=30)
    p.add_argument('--starts',type=int,default=3)
    p.add_argument('--bound',type=float,default=math.sqrt(4*math.pi))
    p.add_argument('--diff-step',type=float,default=1e-4)
    p.add_argument('--target',type=float,default=5.86)
    p.add_argument('--seed',type=int,default=20261030)
    p.add_argument('--plan',action='store_true')
    p.add_argument('--resume',action='store_true')
    args=p.parse_args(argv)
    require(args.max_nfev>=0 and 1<=args.starts<=16 and args.bound>0 and args.target>0,
            'invalid fitting parameters')
    require(args.diff_step>0,'diff step must be positive')
    require(policy.RUN_ID.fullmatch(args.run_id) is not None,'invalid run id')
    models=policy.MODELS if args.models=='all' else tuple(x.strip() for x in args.models.split(','))
    require(models and len(set(models))==len(models) and all(x in policy.MODELS for x in models),
            'invalid model selection')
    root=OUTPUT/args.run_id
    require(not root.exists() or args.plan or args.resume,
            'run directory exists; choose a new --run-id or --resume')
    source_root=(args.source_root if args.source_root.is_absolute() else ROOT/args.source_root).resolve()
    require(source_root.is_relative_to(ROOT.resolve()),'source must stay under local project root')
    require(source_root.is_dir(), 'missing source benchmark directory: '+str(source_root))
    require(args.target_file.is_file(),'missing target file')
    target=OscillationFitTarget.from_json(args.target_file)
    require(target.ordering=='NO' and len(target.observable_names)==5,'expected NO five-observable target')
    rows=[]
    for model in models:
        dest=root/(full.model_tag(model)+'_fit.json')
        if args.resume and dest.is_file():
            doc=neut.load(dest)
            require(doc.get('format')==FORMAT and doc.get('model')==model and
                    doc.get('target_sha256')==hashfile(args.target_file),
                    'incompatible saved fit; use fresh run id')
            # Strict provenance: reject results if cached physics inputs have changed.
            _,_,_,hashes=load_case(source_root,model)
            require(doc['source_sha256']==hashes,'source physics changed; fresh run id required')
        elif args.plan:
            load_case(source_root,model)
            print(model,'READY (frozen EFT1 + saved C5 inputs)',flush=True)
            continue
        else:
            try:
                doc=run_one(model,source_root,target,args)
            except Exception as exc:
                doc={'format':FORMAT,'model':model,'status':'Blocked',
                     'error':type(exc).__name__+': '+str(exc)}
            save(dest,doc)
            print(model,doc['status'],doc.get('chi2'),('ERROR: '+doc['error']) if doc.get('error') else '',flush=True)
        rows.append(doc)
        save(root/'summary.json',{'format':FORMAT,'run_id':args.run_id,
             'rows':rows,'source_root':str(source_root.relative_to(ROOT))})
        columns=['model','status','chi2','starting_seed_chi2','evaluations','failed_evaluations',
                 'additional_self_bilinear_c5_omitted','matching_kind','error']
        with (root/'summary.csv').open('w',newline='',encoding='utf-8') as fh:
            w=csv.DictWriter(fh,fieldnames=columns);w.writeheader()
            for row in rows:w.writerow({k:row.get(k) for k in columns})
    if args.plan:return 0
    return 0 if all(x.get('status')!='Blocked' for x in rows) else 2

if __name__=='__main__':
    raise SystemExit(main())
