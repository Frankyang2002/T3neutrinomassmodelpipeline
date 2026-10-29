"""Paired B(-1)/C(-1) general-scalar pilot on historical fitted inputs.

Historical input/result folders are READ ONLY. Results are NEW, run-ID segregated.
Uses existing Batch7 candidate native RGBeta UV/EFT1 equations, with correct
UV-running versus fermion-start boundaries. Compares zero extra couplings with
seven nonzero extra quartics; optionally adds a mixed scalar mass coefficient.

Optional one-scale Matchete C5 insertion holds the historical Yukawas, MF,
scalar masses and lambdaT3 FIXED between the two quartic choices.
The preview deliberately projects the scalar mass matrix to mu2=0 even when
the candidate RGE generates a nonzero endpoint mu2. Generated mixing is always
recorded, never interpreted as part of the projected C5 result. The
Matchete/RGBeta nonzero coupling normalization is NOT certified. No physical
nonzero-mu2 C5 matching is attempted.

  python -m scripts.RunT3GeneralScalarBCPilot --all
  python -m scripts.RunT3GeneralScalarBCPilot --all --c5-preview
  python -m scripts.RunT3GeneralScalarBCPilot --case large_uv_frozen_0p5 --model B_-1
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Any

import numpy as np

from scripts import RunT3ZeroExtraBoundaryBenchmark as b7

ROOT = Path(__file__).resolve().parents[1]
BASE_OUT = ROOT / 'output/studies/general_scalar_extended'
MODELS = ('B_-1', 'C_-1')
CASES = {
    'small_uv_running_5e-9': ('output/studies/lambda_5e-9_refit', True, 5e-9),
    'small_uv_frozen_5e-9': ('output/studies/fermion_start_lambda_comparison/5e-9', False, 5e-9),
    'large_uv_frozen_0p5': ('output/studies/fermion_start_lambda_comparison/0p5', False, 0.5),
}
SCENARIOS = ('restricted_zero', 'general_quartics', 'general_quartics_and_mixing')
FIT_CUTOFF = 5.86
QUARTIC_FRACTION_DEFAULT = 0.1
MIXING_FRACTION_DEFAULT = 0.01
REPORT_FORMAT = 't3-historical-bc-general-scalar-paired-pilot-v1'
DIRECT = b7.DIRECT


def require(ok: bool, msg: str) -> None:
    if not ok:
        raise ValueError('General scalar B/C pilot: ' + msg)


def safe_tag(model: str) -> str:
    require(model in MODELS, 'unsupported B/C model')
    return model.replace('-', 'm')


def parse_run_id(value: str) -> str:
    require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', value or '')),
            'run ID must be 1–64 alphanumeric/_/- characters')
    return value


def finite_fraction(value: float, name: str, *, maximum: float) -> float:
    x = float(value)
    require(math.isfinite(x) and 0 <= x <= maximum, f'{name} must be finite in [0,{maximum}]')
    return x


def study_source(case: str, model: str, registry: dict, *, root: Path = ROOT):
    """Read only the saved historical fitted config and fit, never fit again."""
    from scripts import CoupledUVRunning as s31
    require(case in CASES and model in MODELS, 'unsupported study/model')
    entry, new, reused = s31.inventories(registry, model)
    source_root, uv_on, expected_lam = CASES[case]
    stem = s31.config_path(entry).name
    conf_path = Path(root) / source_root / 'fitted' / stem
    fit_path = Path(root) / source_root / 'fits' / stem
    cfg = b7.read_json(conf_path)
    fit = b7.read_json(fit_path)
    require(cfg.get('representation') == {
        'd_s1': entry['dimensions'][0], 'd_s2': entry['dimensions'][1],
        'd_f': entry['dimensions'][2], 'alpha': entry['alpha'],
        'shared_scalar': False}, 'historical representation mismatch')
    ordinary = cfg.get('base_state', {}).get('ordinary')
    require(isinstance(ordinary, dict), 'missing historical ordinary couplings')
    lam = b7.complex_of(ordinary.get('lambdaT3'))
    require(lam.imag == 0 and math.isclose(lam.real, expected_lam, rel_tol=1e-12, abs_tol=0),
            'historical topology quartic unexpected')
    require(cfg.get('base_state', {}).get('mu_gev') == 1e7, 'input metadata UV scale changed')
    require(cfg.get('scales', {}) == {
        'mu_fermion_threshold_gev': 1e5, 'mu_scalar_threshold_gev': 1e3,
        'mu_low_gev': 1e2}, 'historical scale settings changed')
    require(fit.get('model') == stem.removesuffix('.json') and
            fit.get('status') in ('Success','PoorFit') and
            math.isclose(float(fit.get('lambdaT3')), expected_lam, rel_tol=1e-12),
            'fit record mismatch')
    chi2 = float(fit['chi2_fit'])
    require(math.isfinite(chi2) and chi2 >= 0, 'invalid fitted chi2')
    if uv_on:
        require(cfg.get('fermion_start_study') is None and
                isinstance(cfg.get('lambda_5e9_study'), dict),
                'UV-running provenance invalid')
    else:
        meta=cfg.get('fermion_start_study', {})
        expected_case='0p5' if expected_lam == 0.5 else '5e-9'
        require(meta.get('case') == expected_case and
                'UV RGEs disabled' in meta.get('boundary_definition', ''),
                'UV-frozen provenance invalid')
    require(len(new) == 1 and len(reused) == 6, 'B/C scalar quartic inventory changed')
    require(all(b7.complex_of(ordinary.get(k)) == 0j for k in reused),
            'historical reused extra quartic is not zero')
    require(not (set(new) & set(ordinary)), 'historical new quartics unexpectedly represented')
    require(not any(k.startswith('muT3General') for k in ordinary),
            'historical quadratic is no longer unrepresented')
    return {'cfg':cfg,'fit':fit,'entry':entry,'new':new,'reused':reused,
            'cfg_path':conf_path,'fit_path':fit_path,'lambdaT3':lam,
            'uv_on':uv_on,'fit_accepted':fit['status']=='Success' and chi2 <= FIT_CUTOFF,
            'chi2':chi2}


def added_boundary(new: tuple[str,...], reused: tuple[str,...],
                   lambda_t3: complex, m1: float, m2: float, *,
                   scenario: str, quartic_fraction: float,
                   mixing_fraction: float) -> tuple[dict[str,complex],complex]:
    """Fixed, deterministic candidate scenario, not a new fit or physical prior."""
    require(scenario in SCENARIOS, 'invalid scenario')
    require(len(new) == 1 and len(reused) == 6 and len(set(new + reused)) == 7,
            'unexpected B/C quartic inventory')
    quartic_fraction=finite_fraction(quartic_fraction, 'quartic fraction', maximum=1.0)
    mixing_fraction=finite_fraction(mixing_fraction,'mixing fraction',maximum=0.20)
    require(math.isfinite(m1) and math.isfinite(m2) and min(m1,m2)>0,
            'non-positive initial scalar masses squared')
    require(math.isfinite(abs(lambda_t3)) and abs(lambda_t3)>0,
            'invalid topology quartic')
    # Give each quartic the same positive real coefficient, including the new
    # sidecar. No claim that this is an SU(2)-normalised physical basis choice.
    q = complex(quartic_fraction*abs(lambda_t3)) if scenario != 'restricted_zero' else 0j
    mu2 = (complex(mixing_fraction*min(m1,m2)) if
           scenario == 'general_quartics_and_mixing' else 0j)
    return {name:q for name in new + reused}, mu2


def prepare_initial(record: dict, registry: dict, scenario: str,
                    quartic_fraction: float, mixing_fraction: float):
    from scripts import CoupledUVRunning as s31
    from Numerical.core.StateVector import unpack_uv_state
    from Numerical.experimental.GeneralT3ScalarRGState import (
        GeneralT3ScalarRGState, pack_general_scalar_state)
    entry,new,reused=record['entry'],record['new'],record['reused']
    cfg=record['cfg']
    m1=float(cfg['base_state']['ordinary']['mS1Sq'])
    m2=float(cfg['base_state']['ordinary']['mS2Sq'])
    extra,mu2=added_boundary(new,reused,record['lambdaT3'],m1,m2,
                             scenario=scenario,quartic_fraction=quartic_fraction,
                             mixing_fraction=mixing_fraction)
    vec0,base_layout=s31.prepare_state(cfg,entry,reused,extra)
    base=unpack_uv_state(vec0,base_layout,mu_gev=1e7)
    common=GeneralT3ScalarRGState(base,{n:extra[n] for n in new}).validated(registry)
    vec,layout=pack_general_scalar_state(common,registry)
    y0=np.concatenate((vec,[mu2.real,mu2.imag]))
    return y0,layout,extra,mu2


def solver_one_segment(y0:np.ndarray,layout:Any, registry:dict,
                       betas:dict[str,str],stage:str):
    start,stop=b7.SCALES[stage]
    rhs=lambda t,y:b7.rhs_at(t,y,layout=layout,registry=registry,
                              betas=betas,stage=stage)
    return b7.integrate_segment(y0,rhs,start,stop,n_points=9)


def compute_trajectory(model:str, record:dict, registry:dict, *,
                       scenario:str, quartic_fraction:float,
                       mixing_fraction:float):
    """UV frozen means no UV solver invocation, including for a nonzero µ2."""
    from Numerical.experimental.GeneralT3ScalarRGState import (
        unpack_general_scalar_state, project_general_after_fermion,
        pack_general_scalar_state)
    new,reused=record['new'],record['reused']
    y0,layout,extra,mu2=prepare_initial(record,registry,scenario,
                                         quartic_fraction,mixing_fraction)
    uv_beta_file=b7.beta_path(model,'UV')
    eft_beta_file=b7.beta_path(model,'EFT1')
    uv_raw=b7.read_json(uv_beta_file)
    beta_uv=b7.check_payload(uv_raw,model=model,stage='UV',
                     rep=layout.representation,
                     base_names=tuple(layout.names[:-len(new)]),
                     new_names=new,reused_names=reused)
    initial_mu=1e7 if record['uv_on'] else 1e5
    starting=b7.snapshot(y0,layout,registry,stage='UV',mu_gev=initial_mu,
                         new=new,reused=reused)
    if record['uv_on']:
        uv=solver_one_segment(y0,layout,registry,beta_uv,'UV')
        y_f=uv.y[:,-1]
        nfev_uv=int(uv.nfev)
        uv_end=b7.snapshot(y_f,layout,registry,stage='UV',mu_gev=1e5,
                            new=new,reused=reused)
    else:
        y_f=y0.copy()
        nfev_uv=0
        uv_end=None
    state_f=unpack_general_scalar_state(y_f[:layout.size],layout,registry,
                                        uv=True,mu_gev=1e5)
    projected=project_general_after_fermion(state_f,registry)
    y1,lay_eft=pack_general_scalar_state(projected,registry)
    y1=np.concatenate((y1,y_f[layout.size:]))
    require(np.array_equal(y1[lay_eft.size:],y_f[layout.size:]),
            'quadratic coupling altered by tree-level fermion projection')
    eft_beta=b7.check_payload(b7.read_json(eft_beta_file),model=model,stage='EFT1',
                     rep=lay_eft.representation,
                     base_names=tuple(lay_eft.names[:-len(new)]),
                     new_names=new,reused_names=reused)
    eft_initial=b7.snapshot(y1,lay_eft,registry,stage='EFT1',mu_gev=1e5,
                            new=new,reused=reused)
    sol=solver_one_segment(y1,lay_eft,registry,eft_beta,'EFT1')
    last=b7.snapshot(sol.y[:,-1],lay_eft,registry,stage='EFT1',mu_gev=1e3,
                     new=new,reused=reused)
    require(min(last['mS1Sq_gev2'],last['mS2Sq_gev2'])>0,'tachyonic diagonal m²')
    keys=list(new+reused)
    for k in keys:
        group='new_quartics' if k in new else 'reused_quartics'
        require(abs(b7.complex_of(eft_initial[group][k])-b7.complex_of(
                (uv_end or starting)[group][k]))<1e-8,
                'non-continuous quartic fermion projection')
    require(abs(b7.complex_of(eft_initial['quadratic_mu2_gev2'])-
                b7.complex_of((uv_end or starting)['quadratic_mu2_gev2']))<1e-8,
            'non-continuous mixed scalar mass fermion projection')
    return {'scenario':scenario,'boundary_mu_gev':initial_mu,
            'uv_integrated':bool(record['uv_on']),
            'nfev_uv':nfev_uv,'nfev_eft1':int(sol.nfev),
            'boundary':starting,'uv_end':uv_end,'eft1_start':eft_initial,
            'eft1_end':last,
            'input_quartics':{k:b7.complex_json(v) for k,v in extra.items()},
            'input_mu2_gev2':b7.complex_json(mu2),
            'beta_system':'candidate RGBeta UV/EFT1, not certified Matchete map'}


def mass_spectrum(model:str, endpoint:dict) -> dict:
    require(model in MODELS, 'unsupported B/C scalar mass model')
    m1=float(endpoint['mS1Sq_gev2']);m2=float(endpoint['mS2Sq_gev2'])
    mu=b7.complex_of(endpoint['quadratic_mu2_gev2'])
    require(math.isfinite(m1) and math.isfinite(m2) and min(m1,m2)>0,
            'negative/nonfinite scalar mass-squared input')
    # Same A=[[m1,mu],[conj(mu),m2]] convention as the existing
    # GeneralT3MixedScalarMassBasis; eigenvalues are gauge-unbroken only.
    A=np.array([[m1,mu],[mu.conjugate(),m2]],dtype=complex)
    eig=np.linalg.eigvalsh(A)
    require(np.isfinite(eig).all() and min(eig)>0,
            'negative scalar mass-squared eigenvalue')
    norm=math.hypot(m1-m2,2*abs(mu))
    sin2=(2*abs(mu)/norm) if norm>0 else 0.0
    return {'eigenvalues_m2_gev2':[float(v) for v in eig],
            'eigenmasses_gev':[float(math.sqrt(v)) for v in eig],
            'mixing_sin_2theta':float(sin2),
            'mu2_over_min_diagonal':float(abs(mu)/min(m1,m2)),
            'candidate_su2_mass_basis_only':True,
            'physical_matchete_coupling_map_certified':False}


def compare_rg(reference:dict,candidate:dict,model:str)->dict:
    a=reference['eft1_end'];b=candidate['eft1_end']
    require(a['mu_gev']==b['mu_gev']==1e3,'mismatched RG endpoint scales')
    m1=[float(a[k]) for k in ('mS1Sq_gev2','mS2Sq_gev2')]
    m2=[float(b[k]) for k in ('mS1Sq_gev2','mS2Sq_gev2')]
    values={}
    for group in ('new_quartics','reused_quartics'):
        require(set(a[group])==set(b[group]),'quartic inventory changed')
        for key in sorted(a[group]):
            za=b7.complex_of(a[group][key]);zb=b7.complex_of(b[group][key])
            values[key]={'restricted':b7.complex_json(za),
                         'general':b7.complex_json(zb),
                         'delta_abs':float(abs(zb-za))}
    return {'quartic_changes':values,
            'scalar_diagonal_m2_relative_shifts':[
                float((m2[i]-m1[i])/m1[i]) for i in range(2)],
            'lambdaT3_relative_change':float(abs(b7.complex_of(b['lambdaT3'])-
                                          b7.complex_of(a['lambdaT3']))/
                                          abs(b7.complex_of(a['lambdaT3']))),
            'reference_spectrum':mass_spectrum(model,a),
            'candidate_spectrum':mass_spectrum(model,b),
            'isolation':'same historical fitted inputs, different added scalar boundary values'}


def _three_by_three_complex(raw:dict, name:str) -> np.ndarray:
    require(isinstance(raw,dict) and set(raw)=={'real','imag'},
            name+': missing complex 3x3 matrix')
    r=np.asarray(raw['real'],dtype=float);i=np.asarray(raw['imag'],dtype=float)
    require(r.shape==(3,3) and i.shape==(3,3) and
            np.isfinite(r).all() and np.isfinite(i).all(),
            name+': bad complex matrix')
    return r+1j*i


def unmixed_c5_projection_info(zero:dict, quartic:dict) -> dict:
    """Declare, rather than conceal, omitted RGE-induced mass mixing.

    This is a hypothetical common diagonal-mass projection, NOT the matched
    coefficient for either actual mixed-mass EFT1 endpoint. The input mixing
    is zero in both scenarios, but one-loop RGE evolution may regenerate it.
    """
    endpoints={}
    for name, scene in (('restricted_zero', zero),
                        ('general_quartics', quartic)):
        end=scene['eft1_end']
        mu=b7.complex_of(end['quadratic_mu2_gev2'])
        m1=float(end['mS1Sq_gev2'])
        m2=float(end['mS2Sq_gev2'])
        require(math.isfinite(m1) and math.isfinite(m2) and min(m1,m2)>0,
                'invalid endpoint diagonal mass² for unmixed C5 projection')
        endpoints[name]={
            'actual_endpoint_mu2_gev2':b7.complex_json(mu),
            'actual_endpoint_mu2_over_min_diagonal_m2':float(abs(mu)/min(m1,m2)),
            'endpoint_mixing_discarded_in_c5_projection':bool(mu!=0j),
        }
    return {
        'projection':'Common hypothetical unmixed scalar basis with mu2=0, '
                     'using restricted_zero endpoint diagonal masses; '
                     'not the actual mixed-mass EFT1 endpoint',
        'projected_mu2_gev2':b7.complex_json(0j),
        'endpoint_mixing':endpoints,
        'actual_mixing_contribution_to_c5_evaluated':False,
        'any_generated_mixing_omitted':any(
            v['endpoint_mixing_discarded_in_c5_projection']
            for v in endpoints.values()),
    }


def one_scale_c5_preview(model:str, record:dict, zero:dict, quartic:dict) -> dict:
    """Isolated quartic insertion in a stated mu2=0 projection only.

    Even quartics-only zero-mu2 *inputs* can generate mu2 along the RGE. The
    actual endpoint mixing remains reported and is NOT propagated into C5.
    """
    projection=unmixed_c5_projection_info(zero,quartic)
    import mpmath as mp
    from scripts.RunGeneralT3NumericalC5 import compile_coefficient, DEFAULT_SOURCE
    from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
    # Source from actual Matchete one-scale ordered LLHH matching; not guessed loop kernel.
    fname={'B_-1':'B_m1_final', 'C_-1':'C_m1_final'}[model]
    coef=DEFAULT_SOURCE/(fname+'_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    require(coef.is_file(),'missing one-scale Matchete C5 expression '+str(coef))
    expression=coef.read_text(encoding='utf-8')
    compiled=compile_coefficient(expression)
    historical=record['cfg']['base_state']['ordinary']
    masses=np.real(np.diag(_three_by_three_complex(historical['MF'],'MF')))
    M=_three_by_three_complex(historical['MF'],'MF')
    require(np.max(abs(M-np.diag(masses)))<1e-8 and min(masses)>0,
            'historical MF is not positive real diagonal; preview blocked')
    y1=_three_by_three_complex(historical['y1'],'y1')
    y2=_three_by_three_complex(historical['y2'],'y2')
    base=zero['eft1_end']
    mixing=MixedDoubletMassBasis.build(model,base['mS1Sq_gev2'],
                                       base['mS2Sq_gev2'],0j)
    require(min(mixing.eigenvalues_sq)>0,'non-positive scalar preview masses')
    mp.mp.dps=55
    U=mixing.U
    bg={f'T3MEU{i+1}{j+1}':mp.mpc(str(U[i,j].real),str(U[i,j].imag))
        for i in range(2) for j in range(2)}
    bg.update({'MSX':mp.mpf(str(math.sqrt(mixing.eigenvalues_sq[0]))),
               'MSY':mp.mpf(str(math.sqrt(mixing.eigenvalues_sq[1]))),
               'lambdaT3':mp.mpc(*[str(v) for v in
                   (b7.complex_of(base['lambdaT3']).real,
                    b7.complex_of(base['lambdaT3']).imag)]),
               'hbar':1/(16*mp.pi**2),'conj':mp.conj,'sqrt':mp.sqrt,'log':mp.log})
    def calc(extra:dict) -> np.ndarray:
        env=bg.copy()
        for rgbeta,matchete in DIRECT.items():
            v=b7.complex_of(extra[rgbeta])
            env[matchete]=mp.mpc(str(v.real),str(v.imag))
        a=[[mp.mpc(0) for _ in range(3)] for _ in range(3)]
        for i in range(3):
            for j in range(3):
                for k,m in enumerate(masses):
                    local=dict(env,MF=mp.mpf(str(m)),
                        y1a=mp.mpc(str(y1[i,k].real),str(y1[i,k].imag)),
                        y1b=mp.mpc(str(y1[j,k].real),str(y1[j,k].imag)),
                        y2a=mp.mpc(str(y2[i,k].real),str(y2[i,k].imag)),
                        y2b=mp.mpc(str(y2[j,k].real),str(y2[j,k].imag)))
                    val=eval(compiled,{'__builtins__':{}},local)
                    require(bool(mp.isfinite(val.real) and mp.isfinite(val.imag)),
                            'nonfinite one-scale matched C5')
                    a[i][j]+=val
        result=np.array([[complex(a[i][j]+a[j][i]) for j in range(3)]
                          for i in range(3)],dtype=complex)
        require(np.isfinite(result.real).all() and np.isfinite(result.imag).all(),
                'invalid one-scale C5 matrix')
        return result
    zeros={source:b7.complex_json(0j) for source in DIRECT}
    varied={source:quartic['eft1_end']['reused_quartics'][source] for source in DIRECT}
    old=calc(zeros);new=calc(varied)
    relative=b7.relative_c5(old,new)
    return {'status':'UncertifiedOneScaleDirectQuarticSensitivity',
            'epsilon_C5_frobenius_uncertified':relative,
            'C5_baseline_norm_gev_inv':float(np.linalg.norm(old,'fro')),
            'C5_direct_delta_norm_gev_inv':float(np.linalg.norm(new-old,'fro')),
            'source_ordered_coefficient':b7.relpath(coef),
            'source_ordered_coefficient_sha256':b7.digest(coef),
            'Yukawas_and_MF':'historical fitted boundary values FROZEN; no RG running',
            'scalar_masses_and_lambdaT3':'restricted_zero EFT1 endpoint held fixed',
            'identification_assumed':'lambda_Matchete = lambda_RGBeta (UNVERIFIED)',
            'quartic_insertion_only_no_mixed_scalar_mass':True,
            'c5_comparison_keeps_mu2_zero_by_projection':True,
            'actual_endpoint_mixing_included_in_c5':False,
            'unmixed_projection_info':projection,
            'physical_C5_or_uncertainty_bound':False,
            'sequential_EFT1_C5_matching_performed':False}


def check_destination(case:str,run_id:str,*,root:Path=ROOT)->Path:
    require(case in CASES, 'invalid case')
    run_id=parse_run_id(run_id)
    dest=(Path(root)/'output/studies/general_scalar_extended'/case/run_id)
    require(dest.resolve().is_relative_to((Path(root)/'output/studies/general_scalar_extended').resolve()),
            'output directory traversal')
    require(not dest.exists(), 'run directory already exists: '+str(dest)+
            ' (choose --run-id for a NEW run; historical/new runs are never overwritten)')
    return dest


def execute_case(case:str,models:tuple[str,...],*,run_id:str,
                 quartic_fraction:float,mixing_fraction:float,
                 include_mixing:bool,c5_preview:bool,root:Path=ROOT)->dict:
    from Numerical.experimental.GeneralT3ScalarRGState import load_registry
    require(case in CASES and models and set(models).issubset(MODELS),
            'invalid model/case selection')
    dest=check_destination(case,run_id,root=root)
    registry=load_registry(b7.REGISTRY)
    rows={}
    # Calculate before publishing anything. Source failures leave no partial
    # results that might be mistaken for a completed study.
    for model in models:
        source=study_source(case,model,registry,root=root)
        scenarios={}
        for name in SCENARIOS[:3 if include_mixing else 2]:
            result=compute_trajectory(model,source,registry,scenario=name,
                quartic_fraction=quartic_fraction,mixing_fraction=mixing_fraction)
            result['spectrum']=mass_spectrum(model,result['eft1_end'])
            scenarios[name]=result
        result={'format':REPORT_FORMAT,'case':case,'model':model,
                'lambdaT3_input':b7.complex_json(source['lambdaT3']),
                'uv_on':source['uv_on'],
                'fit_status':source['fit']['status'],'chi2_fit':source['chi2'],
                'passes_original_fit_screen':source['fit_accepted'],
                'historical_fitted_config_sha256':b7.digest(source['cfg_path']),
                'historical_fit_record_sha256':b7.digest(source['fit_path']),
                'historical_config':source['cfg_path'].relative_to(root).as_posix(),
                'historical_fit':source['fit_path'].relative_to(root).as_posix(),
                'quartic_input_fraction_of_abs_lambdaT3':quartic_fraction,
                'mass_mixing_input_fraction_of_min_m2':mixing_fraction if include_mixing else None,
                'scenarios':scenarios,
                'quartic_vs_zero':compare_rg(scenarios['restricted_zero'],
                                          scenarios['general_quartics'],model),
                'quadratic_vs_zero':(compare_rg(scenarios['restricted_zero'],
                      scenarios['general_quartics_and_mixing'],model)
                      if include_mixing else None),
                'candidate_basis':'Native RGBeta, physical Matchete map uncertified',
                'nonzero_quadratic_c5_computed':False,
                'matching_with_full_general_scalar_potential_completed':False,
                'fitted_yukawas_refitted':False,
                'historical_outputs_overwritten':False,
                'beta_source_sha256':{b7.relpath(b7.beta_path(model,s)):b7.digest(b7.beta_path(model,s))
                                      for s in ('UV','EFT1')},
                'registry_sha256':b7.digest(b7.REGISTRY)}
        if c5_preview:
            result['C5_preview']=one_scale_c5_preview(model,source,
                         scenarios['restricted_zero'],scenarios['general_quartics'])
        rows[model]=result
    # Publish once all 2-model computations have succeeded.
    dest.mkdir(parents=True,exist_ok=False)
    for model,report in rows.items():
        b7.save_json(dest/(safe_tag(model)+'_paired_pilot.json'),report)
    summary={'format':REPORT_FORMAT,'case':case,'run_id':run_id,
             'models':list(rows),'model_results':{
               model:{'chi2_fit':r['chi2_fit'],'passes_original_fit_screen':r['passes_original_fit_screen'],
                      'delta_lambdaT3_abs_fraction':r['quartic_vs_zero']['lambdaT3_relative_change'],
                      'delta_scalar_m2_fractions':r['quartic_vs_zero']['scalar_diagonal_m2_relative_shifts'],
                      'epsilon_C5_preview_uncertified':r.get('C5_preview',{}).get(
                          'epsilon_C5_frobenius_uncertified'),
                      'c5_preview_mu2_projected_to_zero':bool(r.get('C5_preview')),
                      'general_quartics_actual_mu2_over_min_m2':r['scenarios'][
                          'general_quartics']['spectrum']['mu2_over_min_diagonal'],
                      'mixed_mass_sin2theta':(r['scenarios'].get('general_quartics_and_mixing',{})
                                               .get('spectrum',{}).get('mixing_sin_2theta'))}
               for model,r in rows.items()},
             'no_historical_files_modified':True,
             'physical_matchete_rgbeta_coupling_map_certified':False,
             'full_sequential_general_scalar_C5_computed':False}
    b7.save_json(dest/'summary.json',summary)
    (dest/'summary.md').write_text(format_summary(summary, rows),encoding='utf-8')
    return summary


def format_summary(summary:dict,rows:dict)->str:
    lines=['# B/C historical-input general scalar paired pilot', '',
      f"Study: `{summary['case']}`; run ID: `{summary['run_id']}`", '',
      'Candidate native RGBeta UV/EFT1; distinct zero-extra versus nonzero-extra scalar boundaries.',
      'The nonzero extras are independent benchmark assumptions, not fits or physical coupling normalisations.',
      'Original fit quality is unchanged; no Yukawa refit. No general physical C5 or vacuum-stability result.', '',
      '| Model | Original fit chi2 | Pass <=5.86 | Shift abs(lambdaT3) (fraction) | Shift mS1² | Shift mS2² | Direct C5 preview (uncertified) |',
      '|---|---:|---|---:|---:|---:|---:|']
    for m,r in rows.items():
        comparison=r['quartic_vs_zero']
        val=r.get('C5_preview',{}).get('epsilon_C5_frobenius_uncertified')
        c5_display = f"{val:.6g}" if val is not None else '—'
        lines.append(f"| {m} | {r['chi2_fit']:.6g} | "
              f"{'yes' if r['passes_original_fit_screen'] else 'NO'} | "
              f"{comparison['lambdaT3_relative_change']:.6g} | "
              f"{comparison['scalar_diagonal_m2_relative_shifts'][0]:.6g} | "
              f"{comparison['scalar_diagonal_m2_relative_shifts'][1]:.6g} | "
              f"{c5_display} |")
    if all('restricted_zero' in r['scenarios'] and
           'general_quartics' in r['scenarios'] for r in rows.values()):
        lines+=['', '## RGE-generated mixed quadratic coefficient', '',
          'Zero quadratic *input* does not guarantee a zero quadratic *endpoint*.',
          '| Model | Restricted endpoint abs(mu2)/min(m²) | Quartics endpoint abs(mu2)/min(m²) |',
          '|---|---:|---:|']
        for m,r in rows.items():
            base=r['scenarios']['restricted_zero']['spectrum']
            add=r['scenarios']['general_quartics']['spectrum']
            lines.append(f"| {m} | {base['mu2_over_min_diagonal']:.6g} | "
                         f"{add['mu2_over_min_diagonal']:.6g} |")
    if any('general_quartics_and_mixing' in r['scenarios'] for r in rows.values()):
        lines+=['','## Nonzero mass-mixing scenario (C5 blocked)', '',
          '| Model | abs(mu2)/min(m²) at muS | sin(2 theta) |', '|---|---:|---:|']
        for m,r in rows.items():
            d=r['scenarios']['general_quartics_and_mixing']['spectrum']
            lines.append(f"| {m} | {d['mu2_over_min_diagonal']:.6g} | {d['mixing_sin_2theta']:.6g} |")
    lines+=['', '## Interpretation', '',
       '- This is a controlled change of candidate scalar boundary conditions on each **historical saved fitted input**. The output is NOT an updated or refitted fit.',
       '- The C5 preview, if requested, varies only the two direct quartics at the restricted endpoint mass/Yukawa background, with **mu2 explicitly projected to zero for both evaluations**. It assumes an **uncertified identity** RGBeta-to-Matchete map and uses the old one-scale expression.',
       '- One-loop running may induce mu2 even when its boundary value is zero. Such endpoint mixing is tabulated above and in the C5 preview JSON, but is **omitted** from the hypothetical unmixed one-scale C5 sensitivity. This is not the full RG-endpoint C5 change.',
       '- When mu2 is nonzero, mixing can modify C5 at leading one-loop order. No matching is claimed for the actual mixed-mass endpoint or for the optional explicit mixing scenario.',
       '- Historical input/result locations are never written; each run ID is new-only.', '']
    return '\n'.join(lines)


def main()->None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all',action='store_true',help='Run all three historical cases (default)')
    parser.add_argument('--case',choices=tuple(CASES))
    parser.add_argument('--model',choices=MODELS,help='Default: both B_-1 and C_-1')
    parser.add_argument('--quartic-fraction',type=float,default=QUARTIC_FRACTION_DEFAULT)
    parser.add_argument('--mu2-fraction',type=float,default=MIXING_FRACTION_DEFAULT)
    parser.add_argument('--include-mixing',action='store_true')
    parser.add_argument('--c5-preview',action='store_true')
    parser.add_argument('--run-id',default='pilot_q0p1',
                        help='New result ID, must not exist for chosen case(s)')
    args=parser.parse_args()
    require(not(args.all and args.case),'--all and --case are mutually exclusive')
    quartic=finite_fraction(args.quartic_fraction,'quartic fraction',maximum=1.0)
    mu2=finite_fraction(args.mu2_fraction,'mixing fraction',maximum=0.20)
    cases=(args.case,) if args.case else tuple(CASES)
    models=(args.model,) if args.model else MODELS
    for case in cases:
        summary=execute_case(case,models,run_id=args.run_id,
                    quartic_fraction=quartic,mixing_fraction=mu2,
                    include_mixing=args.include_mixing,c5_preview=args.c5_preview)
        dest=BASE_OUT/case/args.run_id
        print(f'{case}: '+', '.join(
            f"{m} (chi2={summary['model_results'][m]['chi2_fit']:.5g}, "
            f"epsilon_C5_preview={summary['model_results'][m]['epsilon_C5_preview_uncertified']})"
            for m in models))
        print('Saved:',dest/'summary.md')


if __name__=='__main__':
    main()
