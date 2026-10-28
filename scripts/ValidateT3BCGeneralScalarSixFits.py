"""Batch 18: read-only, source-checked six-fit B/C exploratory validation.

Default mode replays candidate RG + conditional saved one-scale C5 + NuFIT for
all six accepted fits, samples UV/EFT1 masses and couplings, and checks the
saved fit against the fresh result. NOT sequential matching certification.

Run from the project root:
  python -m scripts.ValidateT3BCGeneralScalarSixFits --run-id bc_six_validation_v1
  python -m scripts.ValidateT3BCGeneralScalarSixFits --run-id bc_six_saved_v1 --saved-only
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT_BASE = Path('output/studies/general_scalar_extended/validation')
CASE_UV = 'small_uv_running_5e-9'
CASES = (
    ('small_uv_running_5e-9','B_-1','b_uv_constrained_pilot_v1','uv'),
    ('small_uv_running_5e-9','C_-1','c_uv_constrained_pilot_v1','uv'),
    ('small_uv_frozen_5e-9','B_-1','b_small_fit_v2','frozen'),
    ('small_uv_frozen_5e-9','C_-1','c_small_fit_v1','frozen'),
    ('large_uv_frozen_0p5','B_-1','b_large_fit_v1','frozen'),
    ('large_uv_frozen_0p5','C_-1','c_large_fit_v2','frozen'),
)
FORMAT = 't3-bc-six-fit-conditional-validation-v1'
CUTOFF = 5.86
Q_BOUND = 4*math.pi
Y_BOUND = math.sqrt(4*math.pi)
OBS = ('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')
TRACE_COLS = ('case','model','run_id','phase','mu_gev','scalar_light_tev',
              'scalar_heavy_tev','min_fermion_tev','rhoF','max_abs_quartic',
              'max_yukawa_singular','quadratic_mu2_abs_gev2')


def require(ok:bool,msg:str)->None:
    if not ok: raise ValueError('Six-fit validation: '+msg)


def load_json(path:Path)->dict:
    require(path.is_file(),'missing file: '+str(path))
    value=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(value,dict),'JSON object required: '+str(path))
    return value


def valid_float(x:Any,label:str,*,positive:bool=False)->float:
    y=float(x)
    require(math.isfinite(y) and (y>0 if positive else True),'invalid '+label)
    return y


def expected_paths(root:Path,case:str,model:str,run_id:str,kind:str)->tuple[Path,Path]:
    base=root/'output/studies/general_scalar_extended'/case/run_id
    if kind=='uv':return base/'fit.json',base/'fitted_config.json'
    tag='B_m1' if model=='B_-1' else 'C_m1'
    return base/(tag+'_exploratory_fit.json'),base/(tag+'_exploratory_fitted_config.json')


def sha_file(path:Path)->str:
    return sha256(path.read_bytes()).hexdigest()


def verify_hashes(root:Path,manifest:dict)->int:
    """Reject source drift, including read-only historical inputs and beta exports."""
    require(isinstance(manifest,dict) and len(manifest)>=4,'missing source hashes')
    for raw,digest in manifest.items():
        require(isinstance(raw,str) and isinstance(digest,str) and len(digest)==64,
                'invalid provenance item '+str(raw))
        path=Path(raw)
        if not path.is_absolute():path=root/path
        require(path.is_file(),'source missing: '+str(path))
        require(sha_file(path)==digest,'source changed since fitting: '+str(path))
    return len(manifest)


def frozen_masses(endpoint:dict)->tuple[float,float]:
    m1=valid_float(endpoint['mS1Sq_gev2'],'mS1Sq',positive=True)
    m2=valid_float(endpoint['mS2Sq_gev2'],'mS2Sq',positive=True)
    mu=endpoint['quadratic_mu2_gev2']
    z=complex(valid_float(mu['real'],'mu2.real'),valid_float(mu['imag'],'mu2.imag'))
    eigen=np.linalg.eigvalsh(np.array([[m1,z],[z.conjugate(),m2]],dtype=complex))
    require(np.all(np.isfinite(eigen)) and np.all(eigen>0),'nonpositive scalar eigenvalue')
    return float(np.sqrt(eigen[0])/1000),float(np.sqrt(eigen[1])/1000)


def check_c5_matrix(raw:dict,label:str)->float:
    """Check saved complex symmetric Weinberg coefficient; return Frobenius norm."""
    require(isinstance(raw,dict) and set(raw)>={'real','imag'},'missing '+label)
    a=np.asarray(raw['real'],dtype=float)+1j*np.asarray(raw['imag'],dtype=float)
    require(a.shape==(3,3) and np.isfinite(a.real).all() and
            np.isfinite(a.imag).all() and
            np.allclose(a,a.T,rtol=1e-10,atol=1e-20),
            'invalid complex symmetric '+label)
    norm=float(np.linalg.norm(a,'fro'))
    require(norm>0,'zero '+label)
    return norm


def load_target(root:Path):
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    return OscillationFitTarget.from_json(root/'configs/nufit_6_1_ic24_no.json')


def checked_chi2(prediction:dict,chi2:float,target:Any)->float:
    require(set(OBS)<=set(prediction),'missing oscillation prediction')
    names=tuple(target.observable_names)
    require(set(names)==set(OBS),'NuFIT target changed its five-observable inventory')
    p=np.asarray([valid_float(prediction[k],k) for k in names])
    resid=np.linalg.solve(np.linalg.cholesky(target.covariance),p-target.central_values)
    actual=float(resid@resid)
    require(math.isclose(actual,chi2,rel_tol=2e-5,abs_tol=2e-5),
            'saved chi2 does not agree with current NuFIT target')
    return actual


def config_lam(cfg:dict)->complex:
    raw=cfg['base_state']['ordinary']['lambdaT3']
    return complex(float(raw['real']),float(raw['imag']))


def config_yukawa(cfg:dict)->dict[str,float]:
    from Numerical.orchestration import RunSmallLambdaRefit as base
    ordinary=cfg['base_state']['ordinary']
    values={}
    for name in ('y1','y2'):
        z=base.matr(ordinary[name]); require(z.shape==(3,3),'invalid Yukawa matrix shape')
        require(np.max(abs(z.imag))<=1e-12,'current fit expects real Yukawas')
        for i in range(3):
            for j in range(3):values[f'{name}_{i}{j}']=float(z[i,j].real)
    return values


def validate_saved(root:Path,case:str,model:str,run_id:str,kind:str,target:Any)->dict:
    path,cfg_path=expected_paths(root,case,model,run_id,kind)
    r=load_json(path); cfg=load_json(cfg_path)
    require(r.get('case')==case and r.get('model')==model and r.get('run_id')==run_id,
            'fit identity mismatch: '+str(path))
    require(cfg.get('representation')=={'d_s1':2,'d_s2':2,
        'd_f':1 if model=='B_-1' else 3,'alpha':-1,'shared_scalar':False},
        'configuration representation mismatch')
    require(r.get('physical_matchete_rgbeta_map_certified') is False and
            r.get('historical_outputs_overwritten') is False,
            'matching/provenance certification changed')
    if kind=='uv':
        require(r.get('format')=='t3-bc-uv-hierarchy-constrained-sobol-refit-v1',
                'unknown UV fit format')
        require(r.get('status')=='CandidateConstrainedAccepted' and
                r.get('general_scalar_sequential_matching_performed') is False,
                'UV point not a conditional accepted fit')
        best=r.get('best_mass_safe')
        require(isinstance(best,dict) and best.get('screen_pass') is True,
                'accepted UV point was not mass screened')
        chi2=valid_float(best['chi2'],'UV chi2')
        require(chi2<=CUTOFF,'UV chi2 above cut')
        require(float(best['rhoF'])<=float(r['mass_screen']['rhoF_max'])<=0.5,
                'UV rhoF screen violated')
        lo=valid_float(best['scalar_light_muS_tev'],'scalar min',positive=True)
        hi=valid_float(best['scalar_heavy_muS_tev'],'scalar max',positive=True)
        require(0.5<=lo<=hi<=5.0,'UV scalar window violated')
        params=best['parameters']
        lam=valid_float(best['lambdaT3_uv'],'UV topology lambda',positive=True)
        require(math.isclose(lam,10**float(params['log10_lambdaT3']),rel_tol=1e-10),
                'UV lambda and log parameter disagree')
        require(r.get('quartic_fraction')==0.1 and r.get('mu2_fraction')==0,
                'new UV scalar boundary prescription')
        manifest=r.get('source_sha256')
        mass={'rhoF':float(best['rhoF']),'light_tev':lo,'heavy_tev':hi}
        prediction=best['prediction']
        c5_norm=check_c5_matrix(best['candidate_C5_muS_gev_inv'],'UV C5 at muS')
        check_c5_matrix(best['candidate_C5_muLow_gev_inv'],'UV C5 at muLow')
    else:
        require(r.get('format') in {
            't3-bc-general-scalar-exploratory-yukawa-refit-v1',
            't3-bc-general-scalar-exploratory-yukawa-continuation-v1'},
            'unknown frozen fit format')
        require(r.get('status')=='ExploratoryAcceptedFit' and
                r.get('passes_approx_chi2_5p86') is True and
                r.get('UV_running_enabled') is False and
                r.get('full_sequential_general_scalar_matching_performed') is False,
                'frozen fit not accepted or wrong RG mode')
        chi2=valid_float(r['chi2_final'],'frozen chi2')
        require(chi2<=CUTOFF,'frozen chi2 above cut')
        endpoint=r['final_diagnostics']['trajectory']['eft1_end']
        lo,hi=frozen_masses(endpoint)
        heavy=[valid_float(x,'fermion mass',positive=True) for x in r['final_diagnostics']['MF_muF_gev']]
        # The scalar mass for the *100 TeV* hierarchy check must also be
        # evaluated at muF, not reused from the later 1 TeV endpoint.
        o=cfg['base_state']['ordinary']
        opts=r['source_manifest'].get('fit_options',r['source_manifest'])
        mixfrac=valid_float(opts['mu2_fraction'],'frozen mu2_fraction')
        require(mixfrac==0,'current accepted frozen inventory expects zero input mu2')
        initial_m1=valid_float(o['mS1Sq'],'frozen scalar input mS1Sq',positive=True)
        initial_m2=valid_float(o['mS2Sq'],'frozen scalar input mS2Sq',positive=True)
        _,heavy_at_f=frozen_masses({'mS1Sq_gev2':initial_m1,
                         'mS2Sq_gev2':initial_m2,
                         'quadratic_mu2_gev2':{'real':0,'imag':0}})
        ratio=heavy_at_f*1000/min(heavy)
        require(ratio<1,'frozen benchmark does not have fermion-first ordering')
        mass={'rhoF':ratio,'light_tev':lo,'heavy_tev':hi}
        params=r['best_parameters'];lam=complex(config_lam(cfg)).real
        source=r['source_manifest']
        manifest=(source.get('prior_input_file_sha256') if 'prior_input_file_sha256' in source
                  else source.get('files_sha256'))
        prediction=r['prediction']
        c5_norm=check_c5_matrix(r['final_diagnostics']['C5_muS_gev_inv'],
                                 'frozen C5 at muS')
        check_c5_matrix(r['final_diagnostics']['C5_low_gev_inv'],
                        'frozen C5 at muLow')
    require(set(params)>=set(config_yukawa(cfg)), 'missing fitted Yukawa coordinates')
    saved_y=config_yukawa(cfg)
    require(all(math.isclose(float(params[k]),v,rel_tol=0,abs_tol=5e-11)
                for k,v in saved_y.items()),'fitted config and fit coordinates differ')
    require(math.isclose(config_lam(cfg).real,lam,rel_tol=1e-10,abs_tol=0) and
            abs(config_lam(cfg).imag)<1e-14,'fitted configuration lambda mismatch')
    checked_chi2(prediction,chi2,target)
    checked=verify_hashes(root,manifest)
    return {'case':case,'model':model,'run_id':run_id,'kind':kind,
         'status':'SavedFitValidated','chi2':chi2,'prediction':prediction,
         'rhoF':mass['rhoF'],'scalar_light_muS_tev':mass['light_tev'],
         'scalar_heavy_muS_tev':mass['heavy_tev'],'lambdaT3_uv':lam,
         'yukawa_max_input':max(abs(x) for x in saved_y.values()),
         'sources_verified':checked,'c5_muS_norm_gev_inv':c5_norm,
         'fit_path':path,'cfg_path':cfg_path,
         'fit_sha256':sha_file(path),'config_sha256':sha_file(cfg_path),
         '_saved_fit':r}


def couplings_of_vector(vector:np.ndarray,layout:Any)->dict:
    """RGBeta-coordinate magnitudes, NOT a potential boundedness test."""
    q=0.;y=0.
    for block in layout.blocks:
        v=np.asarray(vector[block.start:block.stop],dtype=float)
        if block.name.startswith('lambda'):
            if block.kind=='real_scalar':z=abs(float(v[0]))
            elif block.kind=='complex_scalar':z=abs(complex(float(v[0]),float(v[1])))
            else:z=float(np.max(abs(v)))
            q=max(q,float(z))
        elif block.name in ('y1','y2'):
            n=int(np.prod(block.shape))
            z=v[:n].reshape(block.shape)+1j*v[n:].reshape(block.shape)
            y=max(y,float(np.linalg.svd(z,compute_uv=False)[0]))
    return {'max_abs_quartic':q,'max_yukawa_singular':y}


def replay_point(root:Path,row:dict,*,grid:int,target:Any)->tuple[dict,list[dict]]:
    """Full candidate RG replay, including sampled UV/EFT1 coupling diagnostics."""
    from scripts import RunT3GeneralScalarBCExploratoryFit as b11
    from scripts import RunT3GeneralScalarBCPilot as b9
    from scripts import RunT3ZeroExtraBoundaryBenchmark as b7
    from scripts import AnalyzeT3BCMatchingConventions as b10
    from Numerical.experimental.GeneralT3ScalarRGState import (
       load_registry,unpack_general_scalar_state,pack_general_scalar_state,
       project_general_after_fermion)
    from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE
    from scripts import DiagnoseT3BCUVScalarHierarchy as b14
    registry=load_registry(root/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    model=row['model'];kind=row['kind']
    source=b9.study_source(row['case'],model,registry,root=root)
    p=row['_saved_fit']['best_mass_safe']['parameters'] if kind=='uv' else row['_saved_fit']['best_parameters']
    if kind=='uv':
        from scripts import FitT3BCUVHierarchyConstrained as b17
        x=np.array([float(p[k]) for k in b17.NAMES],dtype=float)
        record,_=b17.candidate_record(source,x,lambda_low=1e-14,lambda_high=1.)
        q=float(row['_saved_fit']['quartic_fraction']);mu=float(row['_saved_fit']['mu2_fraction'])
    else:
        x=np.array([float(p[k]) for k in b11.NAMES],dtype=float)
        record=b11.trial_record(source,x)
        meta=row['_saved_fit']['source_manifest']
        o=meta.get('fit_options',meta)
        q=float(o['quartic_fraction']);mu=float(o['mu2_fraction'])
    bundle=b11.beta_bundle(model,record,registry)
    y0,ly0,_,_=b9.prepare_initial(record,registry,
             'general_quartics_and_mixing' if mu>0 else 'general_quartics',q,mu)
    trace=[]
    new,reused=record['new'],record['reused']
    uv_times=[]
    if kind=='uv':
        uv=b7.integrate_segment(y0,lambda t,y:b7.rhs_at(t,y,layout=ly0,
             registry=registry,betas=bundle['UV'],stage='UV'),1e7,1e5,n_points=grid)
        yf=np.asarray(uv.y[:,-1],float)
        uv_times=list(zip(np.geomspace(1e7,1e5,grid),uv.y.T))
    else:
        yf=y0.copy()
        uv_times=[(1e5,y0)]
    for scale,v in uv_times:
        s=b14.sample_vector(model=model,boundary='general_quartics',stage='UV',
              mu_gev=float(scale),vector=np.asarray(v),layout=ly0,
              registry=registry,new=new,reused=reused)
        trace.append(trace_row(row,'UV' if kind=='uv' else 'UV_frozen_boundary',s,
                               couplings_of_vector(v,ly0)))
    at_f=unpack_general_scalar_state(yf[:ly0.size],ly0,registry,uv=True,mu_gev=1e5)
    heavy=b11.heavy_flavour_at_muF(at_f.baseline)
    after=project_general_after_fermion(at_f,registry)
    ye,el=pack_general_scalar_state(after,registry)
    ye=np.concatenate((ye,yf[ly0.size:]))
    betas=b7.check_payload(bundle['payloads']['EFT1'],model=model,stage='EFT1',
        rep=el.representation,base_names=tuple(el.names[:-len(new)]),
        new_names=new,reused_names=reused)
    eft=b7.integrate_segment(ye,lambda t,y:b7.rhs_at(t,y,layout=el,registry=registry,
        betas=betas,stage='EFT1'),1e5,1e3,n_points=grid)
    for scale,v in zip(np.geomspace(1e5,1e3,grid),eft.y.T):
        s=b14.sample_vector(model=model,boundary='general_quartics',stage='EFT1',
              mu_gev=float(scale),vector=np.asarray(v),layout=el,
              registry=registry,new=new,reused=reused,frozen_F_masses=tuple(heavy[0]))
        trace.append(trace_row(row,'EFT1',s,couplings_of_vector(v,el)))
    end=b7.snapshot(eft.y[:,-1],el,registry,stage='EFT1',mu_gev=1e3,
                    new=new,reused=reused)
    state=unpack_general_scalar_state(eft.y[:-2,-1],el,registry,
                                      uv=False,mu_gev=1e3).baseline
    coeff=Path(DEFAULT_SOURCE)/({'B_-1':'B_m1_final','C_-1':'C_m1_final'}[model]
                     +'_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    C5=b11.c5_at_endpoint(model,end,heavy,b11.compile_source(coeff))
    fit,_=b11.low_energy_fit(C5,state,target,
          vev=float(record['cfg'].get('numerical',{}).get('vev_gev',246.22)))
    pred=dict(zip(fit.observable_names,map(float,fit.prediction)))
    replay_chi2=float(fit.chi2)
    require(math.isclose(replay_chi2,row['chi2'],rel_tol=1e-3,abs_tol=5e-5),
            f"{model}/{row['case']} chi2 replay mismatch: {replay_chi2} vs {row['chi2']}")
    for k in OBS:
        require(math.isclose(float(pred[k]),float(row['prediction'][k]),
                            rel_tol=1e-4,abs_tol=2e-7),
                f'{model}/{row["case"]} {k} replay mismatch')
    rho=trace[-grid-1]['rhoF'] if kind=='uv' else trace[0]['rhoF']
    end_t=trace[-1]
    require(math.isclose(rho,row['rhoF'],rel_tol=2e-4,abs_tol=1e-7),'rhoF replay mismatch')
    for field,compare in [('scalar_light_tev','scalar_light_muS_tev'),
                          ('scalar_heavy_tev','scalar_heavy_muS_tev')]:
        require(math.isclose(end_t[field],row[compare],rel_tol=2e-4,abs_tol=1e-6),
                'scalar mass replay mismatch: '+field)
    maxq=max(x['max_abs_quartic'] for x in trace)
    maxy=max(x['max_yukawa_singular'] for x in trace)
    return {'replay_status':'CandidateReplayAgrees','replayed_chi2':replay_chi2,
            'trajectory_samples':len(trace),'max_abs_quartic':maxq,
            'max_yukawa_singular':maxy,
            'quartic_coordinate_screen_pass':maxq<=Q_BOUND,
            'yukawa_singular_screen_pass':maxy<=Y_BOUND,
            'positivity_sampled':True,
            'rg_sampled_scales_gev':[float(r['mu_gev']) for r in trace]},trace


def trace_row(row:dict,phase:str,s:dict,couplings:dict)->dict:
    return {'case':row['case'],'model':row['model'],'run_id':row['run_id'],
        'phase':phase,'mu_gev':s['mu_gev'],
        'scalar_light_tev':s['m_scalar_light_gev']/1000,
        'scalar_heavy_tev':s['m_scalar_heavy_gev']/1000,
        'min_fermion_tev':s['min_fermion_mass_gev']/1000,
        'rhoF':s['max_scalar_over_min_fermion'],
        'quadratic_mu2_abs_gev2':s['mu2_abs_gev2'],**couplings}


def summary_md(results:list[dict],*,replay:bool)->str:
    lines=['# Six accepted B/C general-scalar candidate benchmarks','',
        'Mode: '+('live RG / conditional one-scale C5 replay' if replay else
                  'saved-file consistency only; no RG trajectory verification'),'',
        '| Study | Model | Saved chi2 | rhoF at 100 TeV | Scalars at 1 TeV (TeV) | UV lambdaT3 | Replay |',
        '|---|---|---:|---:|---|---:|---|']
    for d in results:
        lines.append(f"| {d['case']} | {d['model']} | {d['chi2']:.8g} | "
           f"{d['rhoF']:.5g} | {d['scalar_light_muS_tev']:.4g}–"
           f"{d['scalar_heavy_muS_tev']:.4g} | {d['lambdaT3_uv']:.6g} | "
           f"{d.get('replay_status','NotRun')} |")
    lines.extend(['','## Validation','',
        f"- All six saved chi2 values <= {CUTOFF}; predictions, five-variable NuFIT target, Yukawa configurations and source SHA-256 fingerprints checked.",
        '- Scalar eigenmasses are running 2×2 quadratic eigenmasses, not pole masses.',
        '- UV-running cases have the sampled hierarchy ratio from their individual constrained fits; frozen cases are checked against saved fermion masses.',
        '- The UV-running topology couplings were fitted freely and are NOT equal to the historical fixed 5e-9 boundary.',
    ])
    if replay:
        lines.extend(['- All six candidates were independently replayed with the existing native candidate RGBeta + one-scale C5 + SM/Weinberg evaluator.',
          f'- Diagnostic coordinate screens: max |lambda| <= 4pi ({Q_BOUND:.5g}) and max Yukawa singular value <= sqrt(4pi) ({Y_BOUND:.5g}).',
          '- These coordinate screens are not a general-scalar vacuum-stability proof or a physical unitarity calculation.','',
          '| Study | Model | max |lambda| | max Yukawa singular | Quartic/Yukawa screens |',
          '|---|---|---:|---:|---|'])
        for d in results:
            lines.append(f"| {d['case']} | {d['model']} | {d['max_abs_quartic']:.5g} | "
              f"{d['max_yukawa_singular']:.5g} | "
              f"{'PASS' if d['quartic_coordinate_screen_pass'] and d['yukawa_singular_screen_pass'] else 'FLAG'} |")
    lines.extend(['','## Scope boundary','',
        '- B(-1)/C(-1) only: this is 6/6 available candidate fits, not 48/48 general-scalar fits.',
        '- The Matchete–RGBeta coupling map is conditional; the saved C5 is one-scale, not fully matched sequential EFT.',
        '- Positive sampled scalar quadratic eigenvalues do not prove boundedness of the full potential.',
        '- Grouped matching at 100 TeV / 1 TeV is approximate even when masses satisfy the screen.',
        '- Historical and fitted inputs are never overwritten.',''])
    return '\n'.join(lines)


def execute(root:Path,run_id:str,*,grid:int=10,saved_only:bool=False)->dict:
    from scripts.RunT3GeneralScalarBCPilot import parse_run_id
    parse_run_id(run_id)
    require(4<=grid<=50,'grid must contain 4..50 points per running segment')
    dst=root/OUT_BASE/run_id
    require(not dst.exists(),'output already exists: '+str(dst))
    target=load_target(root)
    rows=[];trace=[]
    for case,model,fit_id,kind in CASES:
        row=validate_saved(root,case,model,fit_id,kind,target)
        if not saved_only:
            outcome,samples=replay_point(root,row,grid=grid,target=target)
            row.update(outcome);trace.extend(samples)
        print(f"{case} {model}: saved chi2={row['chi2']:.8g}; "
              f"rhoF={row['rhoF']:.5g}; replay={row.get('replay_status','NotRun')}",flush=True)
        row.pop('_saved_fit',None)
        for k in ('fit_path','cfg_path'):row[k]=str(row[k].relative_to(root))
        rows.append(row)
    perturbative_ok=bool(saved_only or all(
        x['quartic_coordinate_screen_pass'] and x['yukawa_singular_screen_pass']
        for x in rows))
    state=('SavedOnlyValidated' if saved_only else
           'ReplayValidated' if perturbative_ok else 'ReplayValidatedWithPerturbativityFlags')
    report={'format':FORMAT,'status':state,
        'perturbative_coordinate_screens_all_pass':(perturbative_ok if not saved_only else None),
        'candidate_count':len(rows),'chi2_cut':CUTOFF,'matching_physically_certified':False,
        'full_sequential_matching_performed':False,'rg_replay_performed':not saved_only,
        'quartic_coordinate_bound':Q_BOUND,'yukawa_singular_bound':Y_BOUND,
        'results':rows,'historical_outputs_overwritten':False}
    dst.mkdir(parents=True,exist_ok=False)
    (dst/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (dst/'summary.md').write_text(summary_md(rows,replay=not saved_only),encoding='utf-8')
    if not saved_only:
        with (dst/'trajectory_samples.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=TRACE_COLS);w.writeheader();w.writerows(trace)
    print('Saved:',dst/'summary.md',flush=True)
    return report


def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',required=True)
    p.add_argument('--grid',type=int,default=10,
                   help='Points per running interval; default 10')
    p.add_argument('--saved-only',action='store_true',
                   help='Validate stored records only; do not certify RG replay')
    a=p.parse_args()
    execute(ROOT,a.run_id,grid=a.grid,saved_only=a.saved_only)


if __name__=='__main__':main()
