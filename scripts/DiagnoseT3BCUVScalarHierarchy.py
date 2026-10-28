"""Read-only B(-1)/C(-1) UV-running scalar mass-hierarchy diagnostic.

Replay exactly the saved exploratory fit's candidate RGBeta UV/EFT1 state,
with the fitted 18 real Yukawa parameters, on a dense logarithmic grid.
Replay a separate zero-extra boundary on the SAME fitted Yukawas to isolate
whether additional-scalar boundary choices cause the mass growth.

Threshold crossings are approximate running-parameter diagnostics, not pole
masses, physical decoupling scales or a new matching prescription.

    python -m scripts.DiagnoseT3BCUVScalarHierarchy --run-id mass_hierarchy_v1
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path, PurePosixPath
from typing import Any

import numpy as np

from scripts import RunT3ZeroExtraBoundaryBenchmark as b7
from scripts import RunT3GeneralScalarBCPilot as b9
from scripts import RunT3GeneralScalarBCExploratoryFit as b11

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/studies/general_scalar_extended/diagnostics'
CASE = 'small_uv_running_5e-9'
MODELS = ('B_-1', 'C_-1')
DEFAULT_RUNS = {'B_-1': 'general_fresh_v1_B_m1_cont',
                'C_-1': 'general_fresh_v1_C_m1_cont'}
FORMAT = 't3-bc-uv-running-mass-hierarchy-diagnostic-v1'
COLUMNS = ('model', 'boundary', 'phase', 'mu_gev', 'mS1Sq_gev2', 'mS2Sq_gev2',
           'mu2_abs_gev2', 'm_scalar_light_gev', 'm_scalar_heavy_gev',
           'min_fermion_mass_gev', 'max_fermion_mass_gev',
           'max_scalar_over_mu', 'max_scalar_over_min_fermion',
           'x_min_scalar2_over_max_F2', 'x_max_scalar2_over_min_F2')


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('T3 B/C UV mass-hierarchy audit: ' + message)


def safe_input(root: Path, path: str) -> Path:
    p = PurePosixPath(path)
    require(bool(path) and not p.is_absolute() and '..' not in p.parts,
            'unsafe or absolute source path')
    dest = root.joinpath(*p.parts)
    require(dest.is_file(), 'missing prior evidence ' + path)
    return dest


def verify_input_hashes(root: Path, manifest: dict) -> int:
    items = manifest.get('prior_input_file_sha256')
    require(isinstance(items, dict) and items, 'missing prior input fingerprints')
    for path, expected in sorted(items.items()):
        require(isinstance(expected, str) and len(expected) == 64,
                'invalid fingerprint ' + path)
        actual = sha256(safe_input(root, path).read_bytes()).hexdigest()
        require(actual == expected, 'input SHA-256 changed: ' + path)
    return len(items)


def fit_source(root: Path, model: str, run_id: str, registry: dict) -> tuple[dict,dict,dict]:
    require(model in MODELS, 'unsupported model')
    b9.parse_run_id(run_id)
    folder = root/'output/studies/general_scalar_extended'/CASE/run_id
    tag = b9.safe_tag(model)
    fit_path = folder/f'{tag}_exploratory_fit.json'
    cfg_path = folder/f'{tag}_exploratory_fitted_config.json'
    result = b7.read_json(fit_path)
    cfg = b7.read_json(cfg_path)
    require(result.get('model') == model and result.get('case') == CASE and
            result.get('run_id') == run_id and
            result.get('UV_running_enabled') is True and
            result.get('one_scale_matchete_expression_reused') is True and
            result.get('physical_matchete_rgbeta_map_certified') is False,
            'not a compatible saved UV-running exploratory fit')
    require(cfg.get('representation') == {
        'd_s1':2,'d_s2':2,'d_f':1 if model=='B_-1' else 3,
        'alpha':-1,'shared_scalar':False},'saved fitted config representation mismatch')
    values = result.get('best_parameters')
    require(isinstance(values,dict) and set(values)==set(b11.NAMES),
            'saved fit must contain all 18 real Yukawa entries')
    require(all(math.isfinite(float(v)) for v in values.values()),
            'non-finite saved Yukawa parameter')
    require(math.isfinite(float(result['chi2_final'])) and
            float(result['chi2_final'])>=0,'invalid chi-square')
    count=verify_input_hashes(root,result.get('source_manifest',{}))
    q=float(result['source_manifest']['quartic_fraction'])
    mu=float(result['source_manifest']['mu2_fraction'])
    b9.finite_fraction(q,'prior quartic fraction',maximum=1.0)
    b9.finite_fraction(mu,'prior mixing fraction',maximum=0.20)
    historic=b9.study_source(CASE,model,registry,root=root)
    # The saved new fitted configuration is a copy of the historical initial
    # state with modified Yukawas, not an independent source of scalar inputs.
    recorded=b11.seed_vector({'cfg':cfg})
    saved=np.array([float(values[k]) for k in b11.NAMES])
    require(np.allclose(recorded,saved,rtol=0,atol=1e-12),
            'saved fitted config and best parameter vector disagree')
    require(count>=4, 'insufficient signed source fingerprints')
    return result, historic, {'quartic_fraction':q,'mu2_fraction':mu,
                               'saved_fit_path':fit_path,'saved_config_path':cfg_path,
                               'verified_source_files':count}


def scalar_masses(m1: float, m2: float, mu: complex) -> tuple[float,float]:
    require(all(math.isfinite(v) for v in (m1,m2,mu.real,mu.imag)),
            'non-finite scalar mass parameters')
    eigen=np.linalg.eigvalsh(np.array([[m1,mu],[mu.conjugate(),m2]],dtype=complex))
    require(bool(np.all(eigen>0)), 'non-positive running scalar eigenvalue')
    return tuple(math.sqrt(float(x)) for x in eigen)


def sample_vector(*, model:str, boundary:str, stage:str, mu_gev:float,
                  vector:np.ndarray, layout:Any, registry:dict,
                  new:tuple[str,...],reused:tuple[str,...],
                  frozen_F_masses:tuple[float,...]|None=None)->dict:
    from Numerical.experimental.GeneralT3ScalarRGState import unpack_general_scalar_state
    state=unpack_general_scalar_state(vector[:layout.size],layout,registry,
                                      uv=(stage=='UV'),mu_gev=mu_gev)
    if stage=='UV':
        f=tuple(sorted(float(x) for x in np.linalg.svd(state.baseline.MF,compute_uv=False)))
    else:
        require(bool(frozen_F_masses), 'missing fermion-threshold reference masses')
        f=tuple(frozen_F_masses)
    require(all(math.isfinite(x) and x>0 for x in f),'invalid heavy fermion masses')
    m1=float(state.baseline.mS1Sq)
    m2=float(state.baseline.mS2Sq)
    mu=complex(float(vector[layout.size]),float(vector[layout.size+1]))
    light,heavy=scalar_masses(m1,m2,mu)
    return {'model':model,'boundary':boundary,'phase':stage,'mu_gev':mu_gev,
            'mS1Sq_gev2':m1,'mS2Sq_gev2':m2,'mu2_abs_gev2':abs(mu),
            'm_scalar_light_gev':light,'m_scalar_heavy_gev':heavy,
            'min_fermion_mass_gev':min(f),'max_fermion_mass_gev':max(f),
            'max_scalar_over_mu':heavy/mu_gev,
            'max_scalar_over_min_fermion':heavy/min(f),
            'x_min_scalar2_over_max_F2':(light/max(f))**2,
            'x_max_scalar2_over_min_F2':(heavy/min(f))**2}


def crossing(samples:list[dict],*,threshold:str='mu_gev',
             mass:str='m_scalar_heavy_gev')->dict|None:
    """Log-interpolated equality of running positive quantities on a sampled path.

    Returns the adjacent points as an honest bracket. This is NOT a physical
    threshold equation and no self-consistent pole mass is calculated.
    """
    require(len(samples)>=2, 'at least two samples required for crossing')
    for a,b in zip(samples,samples[1:]):
        v0=math.log(float(a[mass])/float(a[threshold]))
        v1=math.log(float(b[mass])/float(b[threshold]))
        if v0==0 or v0*v1<=0:
            z=0. if v0==0 else v0/(v0-v1)
            l0=math.log(float(a['mu_gev']));l1=math.log(float(b['mu_gev']))
            return {'mu_geV_interpolated':float(math.exp(l0+z*(l1-l0))),
                    'bracket_mu_gev':[float(a['mu_gev']),float(b['mu_gev'])],
                    'numerical_crossing_only':True}
    return None


def crossings(samples:list[dict],*,threshold:str,mass:str)->list[dict]:
    out=[]
    for a,b in zip(samples,samples[1:]):
        c=crossing([a,b],threshold=threshold,mass=mass)
        if c:out.append(c)
    return out


def run_branch(model:str,source:dict,registry:dict,bundle:dict,
               *,boundary:str,quartic_fraction:float,mu2_fraction:float,
               grid:int)->dict:
    from Numerical.experimental.GeneralT3ScalarRGState import (
        unpack_general_scalar_state, project_general_after_fermion,
        pack_general_scalar_state)
    is_general=boundary=='general_quartics'
    require(boundary in ('general_quartics','restricted_zero'),'invalid boundary')
    scenario=('general_quartics_and_mixing' if is_general and mu2_fraction>0
              else 'general_quartics' if is_general else 'restricted_zero')
    y0,ly0,_,_=b9.prepare_initial(source,registry,scenario,
                  quartic_fraction if is_general else 0.,
                  mu2_fraction if is_general else 0.)
    new,reused=source['new'],source['reused']
    uv_rhs=lambda t,y:b7.rhs_at(t,y,layout=ly0,registry=registry,
                                betas=bundle['UV'],stage='UV')
    uv=b7.integrate_segment(y0,uv_rhs,1e7,1e5,n_points=grid)
    uv_samples=[sample_vector(model=model,boundary=boundary,stage='UV',
             mu_gev=float(math.exp(t)),vector=uv.y[:,i],layout=ly0,
             registry=registry,new=new,reused=reused)
             for i,t in enumerate(uv.t)]
    yf=uv.y[:,-1]
    at_f=unpack_general_scalar_state(yf[:ly0.size],ly0,registry,uv=True,mu_gev=1e5)
    f=tuple(sorted(float(x) for x in np.linalg.svd(at_f.baseline.MF,compute_uv=False)))
    e_state=project_general_after_fermion(at_f,registry)
    ye,lye=pack_general_scalar_state(e_state,registry)
    ye=np.concatenate((ye,yf[ly0.size:]))
    eft_beta=b7.check_payload(bundle['payloads']['EFT1'],model=model,stage='EFT1',
               rep=lye.representation,base_names=tuple(lye.names[:-len(new)]),
               new_names=new,reused_names=reused)
    eft_rhs=lambda t,y:b7.rhs_at(t,y,layout=lye,registry=registry,
                                  betas=eft_beta,stage='EFT1')
    eft=b7.integrate_segment(ye,eft_rhs,1e5,1e3,n_points=grid)
    eft_samples=[sample_vector(model=model,boundary=boundary,stage='EFT1',
              mu_gev=float(math.exp(t)),vector=eft.y[:,i],layout=lye,
              registry=registry,new=new,reused=reused,frozen_F_masses=f)
              for i,t in enumerate(eft.t)]
    # Record the sign and magnitude of the actual one-loop ODE mass derivatives.
    slopes={}
    for stage,layout,sol,rhs in [('UV',ly0,uv,uv_rhs),('EFT1',lye,eft,eft_rhs)]:
        slopes[stage]={}
        for label,idx in [('first',0),('last',-1)]:
            dy=rhs(float(sol.t[idx]),sol.y[:,idx])
            slopes[stage][label]={name:float(dy[layout.block(name).start])
                                   for name in ('mS1Sq','mS2Sq')}
    full=uv_samples+eft_samples[1:]
    return {'boundary':boundary,'UV':uv_samples,'EFT1':eft_samples,
            'points':full,'uv_nfev':int(uv.nfev),'eft1_nfev':int(eft.nfev),
            'mass_squared_derivatives_dlnmu':slopes,
            'heavy_mass_at_muF_gev':list(f),
            'crossings_heavy_scalar_equals_mu':crossings(full,threshold='mu_gev',
                                                          mass='m_scalar_heavy_gev'),
            'crossings_light_scalar_equals_mu':crossings(full,threshold='mu_gev',
                                                          mass='m_scalar_light_gev'),
            'crossings_heavy_scalar_equals_minF':crossings(full,
                       threshold='min_fermion_mass_gev',mass='m_scalar_heavy_gev')}


def check_saved_endpoint(saved:dict, branch:dict)->dict:
    diagnostic=saved.get('final_diagnostics',{})
    expected=diagnostic.get('trajectory',{}).get('eft1_end',{})
    got=branch['EFT1'][-1]
    checks={}
    for key,report in [('mS1Sq_gev2','mS1Sq_gev2'),('mS2Sq_gev2','mS2Sq_gev2')]:
        original=float(expected[key]);actual=float(got[report])
        require(np.isclose(original,actual,rtol=5e-6,atol=1e-3),
                'rerun endpoint not reproduced for '+key)
        checks[key]={'saved':original,'replayed':actual,'relative_error':float(
                   abs(original-actual)/max(abs(original),1.0))}
    from scripts.RunT3ZeroExtraBoundaryBenchmark import complex_of
    mu_saved=abs(complex_of(expected['quadratic_mu2_gev2']))
    require(np.isclose(mu_saved,got['mu2_abs_gev2'],rtol=1e-5,atol=1e-4),
            'mixed mass endpoint not reproduced')
    fm=sorted(map(float,diagnostic['MF_muF_gev']))
    gn=sorted(map(float,branch['heavy_mass_at_muF_gev']))
    require(np.allclose(fm,gn,rtol=1e-5,atol=1e-4),'heavy masses at muF not reproduced')
    return {'scalar_mass_squared':checks,'mu2_abs_gev2_saved':mu_saved,
            'mu2_abs_gev2_replayed':got['mu2_abs_gev2'],
            'MF_muF_gev_saved':fm,'MF_muF_gev_replayed':gn,
            'saved_candidate_RGBeta_endpoint_reproduced':True}


def diagnostics(general:dict,restricted:dict,saved:dict,model:str)->dict:
    fg=general['UV'][-1];sg=general['EFT1'][-1]
    fr=restricted['UV'][-1];sr=restricted['EFT1'][-1]
    m_f=float(fg['min_fermion_mass_gev']);m_s=float(sg['m_scalar_heavy_gev'])
    over_f=fg['m_scalar_heavy_gev']/m_f
    return {'model':model,'fitted_chi2':float(saved['chi2_final']),
            'passes_chi2_5p86':bool(float(saved['chi2_final'])<=5.86),
            'general_at_muF':fg,'general_at_muS':sg,
            'restricted_at_muF':fr,'restricted_at_muS':sr,
            'general_at_muF_heavy_scalar_over_minF':float(over_f),
            'general_at_muS_heavy_scalar_over_muF':float(m_s/1e5),
            'general_at_muS_heavy_scalar_over_muS':float(m_s/1e3),
            'general_masses_at_muS_gev':[sg['m_scalar_light_gev'],sg['m_scalar_heavy_gev']],
            'restricted_masses_at_muS_gev':[sr['m_scalar_light_gev'],sr['m_scalar_heavy_gev']],
            'delta_mass_sq_at_muS_over_general':{
                k:float((sg[k]-sr[k])/sg[k]) for k in ('mS1Sq_gev2','mS2Sq_gev2')},
            'fermion_first_ordering_at_muF':bool(over_f<1),
            'scalar_below_muF_at_muS':bool(m_s<1e5),
            'scalar_decoupling_scale_close_to_muS':bool(m_s<10e3),
            'conclusion':'Running-parameter ordering only; no physical/pole threshold certification'}


def make_summary(output:dict)->str:
    lines=['# B/C UV-running scalar mass hierarchy diagnostic','',
           'Candidate native RGBeta, saved general-scalar fit Yukawas, no refitting.',
           'Both a general-quartic and a zero-extra boundary are replayed on IDENTICAL Yukawas.',
           'All masses are running mass parameters, not physical pole masses.',
           'The nominal thresholds are muF=100 TeV and muS=1 TeV.', '',
           '| Model | Exploratory chi2 | General max scalar at muF (TeV) | Min F at muF (TeV) | General scalar range at muS (TeV) | F-first ordering at muF | Relative extra-quartic mS1²/mS2² shifts at muS |',
           '|---|---:|---:|---:|---|---|---|']
    for model,report in output['models'].items():
        x=report['diagnostic']; fg=x['general_at_muF'];sg=x['general_at_muS']
        ds=x['delta_mass_sq_at_muS_over_general']
        lines.append(f"| {model} | {x['fitted_chi2']:.6g} | "
                     f"{fg['m_scalar_heavy_gev']/1e3:.5g} | "
                     f"{fg['min_fermion_mass_gev']/1e3:.5g} | "
                     f"{sg['m_scalar_light_gev']/1e3:.5g}–{sg['m_scalar_heavy_gev']/1e3:.5g} | "
                     f"{'YES' if x['fermion_first_ordering_at_muF'] else 'NO'} | "
                     f"{ds['mS1Sq_gev2']:.3g} / {ds['mS2Sq_gev2']:.3g} |")
    lines.extend(['', '## Indicative scale crossings','',
       'Log-interpolated equality m_scalar(mu)=mu (not an on-shell threshold).',''])
    for model,r in output['models'].items():
        cross=r['general']['crossings_heavy_scalar_equals_mu']
        shown=', '.join(f"{v['mu_geV_interpolated']:.6g} GeV "
                        f"(bracket {v['bracket_mu_gev'][0]:.6g}–{v['bracket_mu_gev'][1]:.6g})"
                        for v in cross) or 'No bracketed crossing on the sampled interval'
        lines.append(f'- {model}: {shown}')
    lines.extend(['','## Interpretation','',
       '- Reproducing the saved RG endpoint confirms the replay of the implemented candidate equations, NOT their correctness as physical threshold RGEs.',
       '- The additional-quartic comparison isolates boundary dependence, not a certified Matchete coupling normalization.',
       '- If scalar masses approach or exceed heavy-fermion masses at muF, the F-first hierarchy cannot be treated as separated.',
       '- A running mass crossing with mu is a diagnostic; it is not a pole mass, decoupling prescription, or corrected matching calculation.',
       '- Existing fits, Matchete coefficients, campaign and pilot outputs are not changed.',
       '- Chi2 is the approximate five-observable NuFIT fitting statistic, not a predicted physical goodness-of-fit with 13 degrees of freedom.', ''])
    return '\n'.join(lines)


def execute(*,models:tuple[str,...],run_id:str,grid:int=61,root:Path=ROOT,
            run_ids:dict[str,str]|None=None)->dict:
    require(3<=grid<=100, 'grid must have 3–100 samples per segment')
    b9.parse_run_id(run_id)
    require(models and len(set(models))==len(models) and set(models)<=set(MODELS),
            'models must be unique supported B/C cases')
    dest=root/'output/studies/general_scalar_extended/diagnostics'/run_id
    require(not dest.exists(), 'output already exists; select a NEW --run-id')
    registry=b7.read_json(root/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    report={'format':FORMAT,'status':'CandidateRGERunningMassHierarchyNotPhysicalThreshold',
            'input_case':CASE,'nominal_fermion_threshold_gev':1e5,
            'nominal_scalar_threshold_gev':1e3,'grid_points_per_segment':grid,
            'models':{},'candidate_RGE_only':True,
            'physical_pole_masses_calculated':False,
            'physical_matching_certified':False,
            'historical_outputs_overwritten':False}
    sample_rows=[]
    for model in models:
        rid=(run_ids or DEFAULT_RUNS)[model]
        saved,source,settings=fit_source(root,model,rid,registry)
        trial=b11.trial_record(source,np.array([float(saved['best_parameters'][k])
                                            for k in b11.NAMES]))
        bundle=b11.beta_bundle(model,trial,registry)
        general=run_branch(model,trial,registry,bundle,boundary='general_quartics',
                           quartic_fraction=settings['quartic_fraction'],
                           mu2_fraction=settings['mu2_fraction'],grid=grid)
        restricted=run_branch(model,trial,registry,bundle,boundary='restricted_zero',
                              quartic_fraction=0.,mu2_fraction=0.,grid=grid)
        validation=check_saved_endpoint(saved,general)
        info=diagnostics(general,restricted,saved,model)
        for branch in (general,restricted):
            sample_rows.extend(branch['UV']+branch['EFT1'])
        report['models'][model]={'source_run_id':rid,
                                'saved_fit_sha256':sha256(settings['saved_fit_path'].read_bytes()).hexdigest(),
                                'saved_config_sha256':sha256(settings['saved_config_path'].read_bytes()).hexdigest(),
                                'verified_source_file_count':settings['verified_source_files'],
                                'quartic_fraction':settings['quartic_fraction'],
                                'mu2_fraction':settings['mu2_fraction'],
                                'reproduction':validation,'diagnostic':info,
                                'general':{k:v for k,v in general.items() if k not in ('points','UV','EFT1')},
                                'restricted':{k:v for k,v in restricted.items() if k not in ('points','UV','EFT1')}}
        print(f'{model}: muF maxScalar/minF={info["general_at_muF_heavy_scalar_over_minF"]:.5g}; '
              f'muS scalar={info["general_masses_at_muS_gev"]}; '
              f'fit chi2={info["fitted_chi2"]:.6g}',flush=True)
    # No output files are touched before every model and replay check succeeded.
    dest.mkdir(parents=True,exist_ok=False)
    with (dest/'mass_trajectory_samples.csv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=COLUMNS)
        writer.writeheader()
        for row in sample_rows:
            writer.writerow({k:row[k] for k in COLUMNS})
    (dest/'diagnostic.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (dest/'summary.md').write_text(make_summary(report),encoding='utf-8')
    print('Saved:',dest/'summary.md',flush=True)
    return report


def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',choices=MODELS,help='Default: both models')
    p.add_argument('--run-id',required=True,help='NEW isolated diagnostic run ID')
    p.add_argument('--grid',type=int,default=61,help='samples per RG segment, 3–100')
    p.add_argument('--b-run-id',default=DEFAULT_RUNS['B_-1'])
    p.add_argument('--c-run-id',default=DEFAULT_RUNS['C_-1'])
    a=p.parse_args()
    models=(a.model,) if a.model else MODELS
    execute(models=models,run_id=a.run_id,grid=a.grid,
            run_ids={'B_-1':a.b_run_id,'C_-1':a.c_run_id})


if __name__=='__main__':
    main()
