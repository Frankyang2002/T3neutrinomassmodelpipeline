"""Batch 7: zero-additional-scalar-coupling UV -> EFT1 benchmark for B_-1/C_-1.

Uses pre-existing, live-exported RGBeta UV/EFT1 candidate beta equations.
Initial additional quartics and the allowed holomorphic quadratic are exactly
zero at mu_UV = 1e7 GeV. Reuses the existing frozen fitted baseline UV input
without modifying it. Evolve renormalisable UV to mu_F=1e5 GeV and scalar
EFT1 to mu_S=1e3 GeV; tree-project across the fermion threshold.

Optional --c5-preview inserts the resulting B/C quartics into the repository's
ONE-SCALE Matchete coefficient with a deliberately UNCERTIFIED identity map.
This is a sensitivity illustration with frozen Yukawas/scalar diagonals, NOT a
physical sequential C5 or a prediction/error bound.

    python -m scripts.RunT3ZeroExtraBoundaryBenchmark --all
    python -m scripts.RunT3ZeroExtraBoundaryBenchmark --all --c5-preview
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Callable, Mapping

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/scalar_relevance/zero_boundary_benchmark'
MODELS = ('B_-1', 'C_-1')
SCALES = {'UV': (1.0e7, 1.0e5), 'EFT1': (1.0e5, 1.0e3)}
REGISTRY = ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
BETA_DIR = ROOT / 'output/general_t3_quadratic_eft1'
MASS_SYMBOL = 'muT3GeneralS1S2'
BETA_FORMAT = 'general-t3-quadratic-rgbeta-experimental-v1'
BETA_STATUS = 'SuccessCandidateMassBeta'
CONVENTION = '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'
DIRECT = {
    'lambdaHHdagS2S2': 'lambdaT3ExtraQ020020I1',
    'lambdaHHdagS1barS1bar': 'lambdaT3ExtraQ020200I1',
}


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError('Batch 7 zero-boundary benchmark: ' + reason)


def read_json(path: Path) -> dict:
    require(path.is_file(), 'missing '+str(path))
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), 'not a JSON object: '+str(path))
    return obj


def digest(path: Path) -> str:
    require(path.is_file(), 'missing provenance file '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def relpath(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path.resolve()).replace('\\', '/')


def save_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n',encoding='utf-8')


def complex_json(z: complex) -> dict[str,float]:
    z=complex(z)
    require(np.isfinite([z.real,z.imag]).all(), 'nonfinite complex scalar')
    return {'real':float(z.real),'imag':float(z.imag)}


def complex_of(value: Mapping[str,float]) -> complex:
    require(isinstance(value, dict) and set(value)=={'real','imag'},
            'invalid complex JSON object')
    return complex(float(value['real']),float(value['imag']))


def beta_path(model: str, stage: str, folder: Path = BETA_DIR) -> Path:
    require(model in MODELS and stage in SCALES, 'unsupported model/stage')
    tag=model.replace('-','m').replace('+','p')
    return folder / f'{tag}_{stage.lower()}_quadratic_live_rgbeta.json'


def check_payload(payload: dict, *, model: str, stage: str,
                  rep: Any, base_names: tuple[str,...],
                  new_names: tuple[str,...], reused_names: tuple[str,...]) -> dict[str,str]:
    """Fail closed on candidate source, stage, representation and all beta keys."""
    require(model in MODELS and stage in SCALES, 'invalid beta target')
    meta=payload.get('metadata')
    require(isinstance(meta,dict), 'missing beta metadata')
    require(payload.get('format')==BETA_FORMAT and payload.get('status')==BETA_STATUS
            and payload.get('physical_mapping_certified') is False
            and payload.get('mass_eigenstate_threshold_matching_performed') is False,
            'beta export status or physical-certification flags changed')
    expected={
        'GeneralScalarModel':model, 'GeneralScalarStage':stage,
        'GeneralScalarBasisStatus':'CandidateNativeRGBetaUnmappedToMatchete',
        'GeneralScalarQuadraticMassMixing':'IncludedCandidate',
        'QuadraticKind':'mixed', 'QuadraticMassSymbol':MASS_SYMBOL,
        'QuadraticPrefactorConvention':'V2 = mu2 * invariant + h.c.; no 1/2',
        'ReportBetaConvention':CONVENTION, 'SharedScalar':False,
        'dS1':rep.d_s1, 'dS2':rep.d_s2, 'dF':rep.d_f,
        'alpha':rep.alpha,
        'PhysicalMatcheteRGBetaMassBasisCertified':False,
    }
    require(all(meta.get(k)==v for k,v in expected.items()),
            model+'/'+stage+': candidate beta metadata mismatch')
    require(set(meta.get('GeneralScalarNewBetaKeys',[]))==set(new_names)
            and meta.get('GeneralScalarAdditionalQuartics')==len(new_names)+len(reused_names)
            and meta.get('GeneralScalarExistingStructuralCount')==len(reused_names),
            model+'/'+stage+': quartic registry mismatch')
    wanted=set(base_names)|set(new_names)|{MASS_SYMBOL}
    betas=payload.get('report_betas')
    require(isinstance(betas,dict) and set(betas)==wanted
            and set(payload.get('beta_symbols',[]))==wanted
            and len(base_names)==len(set(base_names))
            and not(set(base_names)&set(new_names)),
            model+'/'+stage+': missing/extra beta symbols')
    require(all(isinstance(expr,str) and expr.strip()
                and not any(bad in expr for bad in ('BetaTerm[','Missing[','$Aborted','$Failed'))
                for expr in betas.values()),
            model+'/'+stage+': unresolved beta expressions')
    return betas


def rhs_at(t: float, vector: np.ndarray, *, layout: Any, registry: dict,
           betas: Mapping[str,str], stage: str) -> np.ndarray:
    """Evolve the *full* native RGBeta state, including complex mu2."""
    from Numerical.experimental.GeneralT3ScalarRGState import unpack_general_scalar_state
    from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression, state_environment
    from Numerical.running.IntermediateScalarRGBetaEvaluator import intermediate_scalar_state_environment
    from Numerical.core.BetaVector import beta_values_to_derivative

    v=np.asarray(vector,dtype=float).copy()
    require(v.shape==(layout.size+2,) and np.isfinite(v).all(),
            'malformed/nonfinite extended RG vector')
    require(stage in SCALES,'invalid stage')
    # Mirror the UV Majorana mass symmetry enforcement in the existing ODE code.
    if stage=='UV' and layout.representation.majorana_fermion:
        block=layout.block('MF')
        n=int(np.prod(block.shape))
        r=v[block.start:block.start+n].reshape(block.shape)
        im=v[block.start+n:block.stop].reshape(block.shape)
        r[:]=0.5*(r+r.T)
        im[:]=0.5*(im+im.T)
    state=unpack_general_scalar_state(v[:layout.size],layout,registry,
                                      uv=stage=='UV',mu_gev=float(np.exp(t)))
    env=(state_environment(state.baseline) if stage=='UV' else
         intermediate_scalar_state_environment(state.baseline))
    env.update(state.extra)
    env[MASS_SYMBOL]=complex(float(v[layout.size]),float(v[layout.size+1]))
    vals={name:evaluate_inputform_expression(expr,env)
          for name,expr in betas.items()}
    mass=vals.pop(MASS_SYMBOL)
    require(not isinstance(mass,np.ndarray), 'quadratic beta unexpectedly matrix valued')
    mass=complex(mass)
    require(np.isfinite([mass.real,mass.imag]).all(), 'nonfinite quadratic beta')
    deriv=beta_values_to_derivative(vals,layout,divide_by_loop_factor=True)
    return np.concatenate((deriv,[mass.real/(16*math.pi**2),
                                  mass.imag/(16*math.pi**2)]))


def integrate_segment(y0: np.ndarray, rhs: Callable[[float,np.ndarray],np.ndarray],
                      start: float, stop: float, *, n_points: int = 20):
    require(start>stop>0 and 3<=n_points<=100, 'invalid running interval or sample count')
    y0=np.asarray(y0,dtype=float)
    require(y0.ndim==1 and np.isfinite(y0).all(), 'invalid initial numerical state')
    first=rhs(math.log(start),y0)
    require(first.shape==y0.shape and np.isfinite(first).all(), 'bad starting RGE derivative')
    ts=np.linspace(math.log(start), math.log(stop),n_points)
    sol=solve_ivp(rhs,(ts[0],ts[-1]),y0,method='DOP853',
                  rtol=1e-8,atol=1e-11,t_eval=ts,max_step=0.35)
    require(sol.success and sol.y.shape==(len(y0),n_points)
            and np.isfinite(sol.y).all(), 'coupled RGBeta integration failed: '+sol.message)
    return sol


def snapshot(vector: np.ndarray, layout: Any, registry: dict,
             *, stage: str, mu_gev: float, new: tuple[str,...],
             reused: tuple[str,...]) -> dict:
    from Numerical.experimental.GeneralT3ScalarRGState import unpack_general_scalar_state
    state=unpack_general_scalar_state(np.asarray(vector[:layout.size],float),layout,
                                      registry,uv=stage=='UV',mu_gev=mu_gev)
    base=state.baseline
    return {
        'mu_gev':mu_gev, 'mS1Sq_gev2':float(base.mS1Sq),
        'mS2Sq_gev2':float(base.mS2Sq),
        'lambdaT3':complex_json(base.lambdaT3),
        'new_quartics':{k:complex_json(state.extra[k]) for k in new},
        'reused_quartics':{k:complex_json(getattr(base,k)) for k in reused},
        'quadratic_mu2_gev2':complex_json(complex(vector[layout.size],vector[layout.size+1])),
    }


def check_zero_boundary(first: dict, *, new: tuple[str,...], reused: tuple[str,...]) -> None:
    require(set(first['new_quartics'])==set(new) and
            set(first['reused_quartics'])==set(reused), 'wrong zero-boundary quartic inventory')
    require(all(complex_of(x)==0 for x in [*first['new_quartics'].values(),
                *first['reused_quartics'].values(),first['quadratic_mu2_gev2']]),
            'additional scalar UV boundary is not exactly zero')


def inspect_generation(uv_start: dict, uv_end: dict, eft_end: dict,
                       *, new: tuple[str,...], reused: tuple[str,...]) -> dict:
    """Compare only like-for-like native RGBeta coupling coordinates."""
    check_zero_boundary(uv_start,new=new,reused=reused)
    symbols=tuple(sorted((*new,*reused)))
    values={}
    for symbol in symbols:
        family='new_quartics' if symbol in new else 'reused_quartics'
        at_f=complex_of(uv_end[family][symbol])
        at_s=complex_of(eft_end[family][symbol])
        values[symbol]={
            'at_muF':complex_json(at_f), 'at_muS':complex_json(at_s),
            'magnitude_muF':float(abs(at_f)), 'magnitude_muS':float(abs(at_s)),
            'status':'candidate_RGBeta_coordinate_not_physical_Matchete_map'
        }
    lam=complex_of(eft_end['lambdaT3'])
    require(abs(lam)>1e-14,'lambdaT3 vanishes: ratios undefined')
    direct={k:{'rgbeta_at_muS':values[k]['at_muS'],
               'absolute_ratio_to_lambdaT3':float(values[k]['magnitude_muS']/abs(lam)),
               'matchete_target':v,'map_certified':False}
            for k,v in DIRECT.items()}
    mass=complex_of(eft_end['quadratic_mu2_gev2'])
    min_mass=min(eft_end['mS1Sq_gev2'],eft_end['mS2Sq_gev2'])
    require(min_mass>0,'non-positive final diagonal scalar mass')
    return {
        'all_added_couplings_zero_at_UV':True,
        'additional_quartics':values,
        'direct_C5_candidate_quartics':direct,
        'lambdaT3_muS':complex_json(lam),
        'quadratic_mu2_muS_gev2':complex_json(mass),
        'quadratic_fraction_of_min_diagonal_m2':float(abs(mass)/min_mass),
        'physical_matchete_rgbeta_map_certified':False,
        'full_sequential_C5_matched':False,
    }


def c5_matrix(obj: dict, name: str) -> np.ndarray:
    rows=obj.get('symmetrized_C5_gev_inv')
    require(isinstance(rows,list) and len(rows)==3 and all(isinstance(r,list) and len(r)==3 for r in rows),
            name+': missing 3x3 C5')
    arr=np.asarray([[complex_of(c) for c in row] for row in rows],dtype=complex)
    require(np.isfinite(arr.real).all() and np.isfinite(arr.imag).all(),name+': nonfinite C5')
    require(np.allclose(arr,arr.T,rtol=1e-9,atol=1e-20),name+': asymmetric C5')
    return arr


def relative_c5(a: np.ndarray,b: np.ndarray) -> float:
    a=np.asarray(a,dtype=complex);b=np.asarray(b,dtype=complex)
    require(a.shape==(3,3) and b.shape==(3,3), 'C5 arrays must be 3x3')
    den=float(np.linalg.norm(a,'fro'))
    require(np.isfinite(den) and den>0,'reference C5 norm vanishes')
    ratio=float(np.linalg.norm(b-a,'fro')/den)
    require(math.isfinite(ratio),'nonfinite C5 deviation')
    return ratio


def one_scale_preview(model: str, end: dict, *, folder: Path) -> dict:
    """UNCERTIFIED direct one-scale Matchete sensitivity; never sequential C5."""
    from scripts.RunGeneralT3FrozenScalarC5 import build_scenario
    from scripts.RunGeneralT3NumericalC5 import (
        DEFAULT_SOURCE, DEFAULT_PARAMETERS, numerical_model,
    )
    from scripts.RunGeneralT3ReferenceC5 import MANIFEST, REFERENCE_ID
    from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis

    require(model in MODELS,'unsupported preview model')
    # The frozen matching expression is only valid with diagonal scalar inputs.
    mu=complex_of(end['quadratic_mu2_gev2'])
    min_m2=min(end['mS1Sq_gev2'],end['mS2Sq_gev2'])
    require(abs(mu)/min_m2<1e-10,
            'radiative quadratic mixing is not negligible; frozen diagonal preview blocked')
    manifest=read_json(MANIFEST)
    slug=REFERENCE_ID[model][0]
    frozen_path=MANIFEST.parent/(slug+'.json')
    params_path=DEFAULT_PARAMETERS/(model.replace('-','m')+'_final_parameters.json')
    base,basis=build_scenario(model,manifest,frozen_path,params_path,0j,
                              mixing_explicit_on_cli=False)
    require(isinstance(basis,MixedDoubletMassBasis),'invalid frozen diagonal basis')
    common_lambda=complex_of(end['lambdaT3'])
    for cfg in (base,):
        cfg['quartics']['lambdaT3']=complex_json(common_lambda)
    baseline=deepcopy(base)
    corrected=deepcopy(base)
    for source,target in DIRECT.items():
        baseline['quartics'][target]=complex_json(0j)
        corrected['quartics'][target]=end['reused_quartics'][source]
    p0=folder/'one_scale_preview'/'baseline'
    p1=folder/'one_scale_preview'/'identity_map'
    old,path0=numerical_model(model,baseline,DEFAULT_SOURCE,DEFAULT_PARAMETERS,p0,
                              scalar_basis_override=basis)
    new,path1=numerical_model(model,corrected,DEFAULT_SOURCE,DEFAULT_PARAMETERS,p1,
                              scalar_basis_override=basis)
    for item in (old,new):
        require(item.get('physical_c5_certified') is False and
                item.get('sequential_EFT1_performed') is False and
                item.get('frozen_scalar_aligned_one_scale_scenario') is True,
                'one-scale preview falsely claims physical certification')
    eps=relative_c5(c5_matrix(old,'baseline'),c5_matrix(new,'corrected'))
    return {
        'status':'UncertifiedIdentityCouplingMapOneScaleC5SensitivityOnly',
        'epsilon_C5_frobenius_uncertified':eps,
        'baseline_C5_fro_norm_gev_inv':float(np.linalg.norm(c5_matrix(old,'baseline'),'fro')),
        'baseline_largest_takagi_mass_eV':old['largest_takagi_mass_eV'],
        'identity_map_largest_takagi_mass_eV':new['largest_takagi_mass_eV'],
        'baseline_file':relpath(path0),'identity_map_file':relpath(path1),
        'assumed_rgbeta_to_matchete_coupling_map':'lambda_Matchete = lambda_RGBeta (UNVERIFIED)',
        'frozen_yukawas_used_not_running_ones':True,
        'frozen_scalar_diagonals_used_not_RG_final_ones':True,
        'not_full_sequential_matching':True,
        'not_physical_uncertainty_bound':True,
    }


def run_model(model: str, *, outdir: Path = OUT, preview: bool = False,
              n_points: int = 20) -> dict:
    from scripts import CoupledUVRunning as s31
    from scripts import CoupledA0UVRunning as s30
    from Numerical.core.StateVector import unpack_uv_state
    from Numerical.experimental.GeneralT3ScalarRGState import (
        load_registry, GeneralT3ScalarRGState,
        pack_general_scalar_state, unpack_general_scalar_state,
        project_general_after_fermion,
    )
    from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
    from Numerical.core.RGBetaEvaluator import state_environment

    require(model in MODELS,'unsupported model')
    registry=load_registry(REGISTRY)
    entry,new,reused=s31.inventories(registry,model)
    require(len(new)==1 and len(reused)==6 and
            set(DIRECT).issubset(set(reused)),model+': unexpected B/C general quartic registry')
    config_path=s31.config_path(entry)
    config=read_json(config_path)
    # Stage31.prepare_state ensures the fitted benchmark and reused quartic
    # inventory agree, but we supply zeros instead of the Stage29 illustrative values.
    base_vector,base_layout=s31.prepare_state(config,entry,reused,
                                               {name:0j for name in reused})
    base=unpack_uv_state(base_vector,base_layout,mu_gev=s30.START)
    require(math.isclose(base.mu_gev,SCALES['UV'][0],rel_tol=1e-13),
            'benchmark UV scale does not match required 1e7 GeV')
    initial=GeneralT3ScalarRGState(base,{name:0j for name in new}).validated(registry)
    vector,layout=pack_general_scalar_state(initial,registry)
    y0=np.concatenate((vector,np.array([0.,0.])))
    uv_path=beta_path(model,'UV')
    eft_path=beta_path(model,'EFT1')
    uv_raw=read_json(uv_path)
    beta_uv=check_payload(uv_raw,model=model,stage='UV',rep=base.representation,
                          base_names=tuple(layout.names[:-len(new)]),
                          new_names=new,reused_names=reused)
    base_key=s31.DATA/f'stage31_{model.replace("-","m")}_full_uv_betas.json'
    old_raw=read_json(base_key)
    base_betas=s31.validate_export(old_raw,entry,new,reused,base_layout.names)
    # Check equality of beta systems at the zero-mu2 boundary. This checks
    # against the original independent UV export, not only the quadratic one.
    old_rhs=s30.combined_rhs(math.log(s30.START),y0[:-2],layout=base_layout,
                            names=new,betas=base_betas)
    new_rhs=rhs_at(math.log(s30.START),y0,layout=layout,registry=registry,
                   betas=beta_uv,stage='UV')[:-2]
    relative_beta_residual=float(np.max(abs(old_rhs-new_rhs)/np.maximum(1.,abs(old_rhs))))
    require(relative_beta_residual<1e-9,
            'zero-mu2 UV quadratic beta export differs from independent Stage31 UV equations')
    run_uv=lambda t,y:rhs_at(t,y,layout=layout,registry=registry,betas=beta_uv,stage='UV')
    uv=integrate_segment(y0,run_uv,*SCALES['UV'],n_points=n_points)
    uv_first=snapshot(uv.y[:,0],layout,registry,stage='UV',mu_gev=SCALES['UV'][0],
                      new=new,reused=reused)
    uv_end=snapshot(uv.y[:,-1],layout,registry,stage='UV',mu_gev=SCALES['UV'][1],
                    new=new,reused=reused)
    check_zero_boundary(uv_first,new=new,reused=reused)
    uv_obj=unpack_general_scalar_state(uv.y[:layout.size,-1],layout,registry,
                                       uv=True,mu_gev=SCALES['UV'][1])
    aft=project_general_after_fermion(uv_obj,registry)
    y1,ly1=pack_general_scalar_state(aft,registry)
    y1=np.concatenate((y1,uv.y[layout.size:,-1]))
    require(np.array_equal(y1[ly1.size:],uv.y[layout.size:,-1]),
            'quadratic boundary was modified at fermion threshold')
    eft_raw=read_json(eft_path)
    beta_eft=check_payload(eft_raw,model=model,stage='EFT1',rep=aft.baseline.representation,
                          base_names=tuple(ly1.names[:-len(new)]),
                          new_names=new,reused_names=reused)
    run_eft=lambda t,y:rhs_at(t,y,layout=ly1,registry=registry,betas=beta_eft,stage='EFT1')
    eft=integrate_segment(y1,run_eft,*SCALES['EFT1'],n_points=n_points)
    eft_first=snapshot(eft.y[:,0],ly1,registry,stage='EFT1',mu_gev=SCALES['EFT1'][0],
                       new=new,reused=reused)
    eft_end=snapshot(eft.y[:,-1],ly1,registry,stage='EFT1',mu_gev=SCALES['EFT1'][1],
                     new=new,reused=reused)
    # The continuous projection changes only the dynamical field content.
    for key in (*new,*reused):
        sec='new_quartics' if key in new else 'reused_quartics'
        require(abs(complex_of(eft_first[sec][key])-complex_of(uv_end[sec][key]))<1e-10,
                'quartic discontinuity at projected fermion threshold: '+key)
    require(abs(complex_of(eft_first['quadratic_mu2_gev2'])-
                complex_of(uv_end['quadratic_mu2_gev2']))<1e-10,
            'quadratic mass discontinuity at projected fermion threshold')
    comparison=inspect_generation(uv_first,uv_end,eft_end,new=new,reused=reused)
    folder=outdir/model.replace('-','m')
    inputs=(REGISTRY,config_path,uv_path,eft_path,base_key)
    # Treat a small nonzero mu2 induced at higher order as a scenario effect,
    # not permission to plug it into a diagonal one-scale C5 expression.
    report={
        'format':'t3-zero-extra-boundary-bc-renormalisable-running-v1',
        'status':'CandidateBCZeroBoundaryTwoSegmentRunning',
        'model':model, 'boundary_scale_gev':SCALES['UV'][0],
        'fermion_threshold_gev':SCALES['UV'][1],
        'scalar_threshold_gev':SCALES['EFT1'][1],
        'boundary_all_additional_couplings_zero':True,
        'uv_extra_quartic_slots':len(new)+len(reused),
        'uv_quadratic_slots':1,
        'uv_reused_independent_beta_residual':relative_beta_residual,
        'uv_solver_nfev':uv.nfev,'eft1_solver_nfev':eft.nfev,
        'uv_initial':uv_first,'uv_final':uv_end,
        'eft1_initial':eft_first,'eft1_final':eft_end,
        'generation':comparison,
        'source_sha256':{relpath(p):digest(p) for p in inputs},
        'fermion_threshold_convention':'tree-level projection; no finite threshold matching',
        'RG_beta_basis':'candidate native RGBeta, unmapped to Matchete',
        'scalar_quadratic_beta_included_both_segments':True,
        'physical_matchete_rgbeta_map_certified':False,
        'C5_preview_performed':bool(preview),
        'full_sequential_C5_computed':False,
        'physical_neutrino_mass_error_bound':False,
    }
    if preview:
        report['one_scale_C5_preview']=one_scale_preview(model,eft_end,folder=folder)
    save_json(folder/'result.json',report)
    return report


def summary_md(rows: list[dict]) -> str:
    lines=['# T3 B/C zero-extra-coupling UV-to-EFT1 benchmark (Batch 7)', '',
           'All additional quartics and the additional quadratic mass coefficient are exactly zero at 1e7 GeV.',
           'Two native candidate one-loop RGBeta stages are integrated with a tree-level fermion projection.',
           'Physical Matchete/RGBeta normalisation and full C5 matching are not certified.', '',
           '| Model | |lambda_extra(S2)| | |lambda_extra(S1dag)| | |mu2|/min(m²) |',
           '|---|---:|---:|---:|']
    for row in rows:
        g=row['generation']
        d=g['direct_C5_candidate_quartics']
        a=d['lambdaHHdagS2S2']['rgbeta_at_muS'];b=d['lambdaHHdagS1barS1bar']['rgbeta_at_muS']
        lines.append(f"| {row['model']} | {abs(complex_of(a)):.7g} | {abs(complex_of(b)):.7g} | {g['quadratic_fraction_of_min_diagonal_m2']:.4g} |")
    lines+=['', 'The two B/C quartics can contribute to leading one-loop C5 if independently nonzero.',
            'With zero UV values their radiative contributions begin beyond strict fixed-order one-loop C5.',
            'A small generated coupling or mixing coefficient is NOT a neutrino-mass error bar.', '']
    if any(x.get('C5_preview_performed') for x in rows):
        lines+=['## Uncertified identity-map one-scale C5 sensitivity', '',
                '| Model | epsilon_C5 (uncertified) | Baseline largest mass (eV) | Identity-map largest mass (eV) |',
                '|---|---:|---:|---:|']
        for x in rows:
            if not x.get('C5_preview_performed'):continue
            p=x['one_scale_C5_preview']
            lines.append(f"| {x['model']} | {p['epsilon_C5_frobenius_uncertified']:.7g} | {p['baseline_largest_takagi_mass_eV']:.7g} | {p['identity_map_largest_takagi_mass_eV']:.7g} |")
        lines+=['', '**Not a physical error bound**: the identity Matchete/RGBeta map is unverified;',
                'frozen Yukawas and scalar diagonals are used, not their RG-evolved values;',
                'this evaluates a simultaneous one-scale coefficient, not sequential EFT1 matching.', '']
    lines+=['## Reproducibility', '', 'Per-model source hashes and endpoint couplings are stored in `B_m1/result.json` and `C_m1/result.json`.',
            'The production benchmark files and GitHub repository have not been modified.', '']
    return '\n'.join(lines)


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    select=parser.add_mutually_exclusive_group(required=True)
    select.add_argument('--all',action='store_true')
    select.add_argument('--model',choices=MODELS)
    parser.add_argument('--c5-preview',action='store_true',
                        help='uncertified native-RGBeta to Matchete identity-map one-scale C5 sensitivity')
    parser.add_argument('--outdir',type=Path,default=OUT)
    parser.add_argument('--points',type=int,default=20)
    args=parser.parse_args()
    require(args.all or args.model, 'select --all or --model')
    rows=[]
    for model in (MODELS if args.all else (args.model,)):
        report=run_model(model,outdir=args.outdir,preview=args.c5_preview,
                         n_points=args.points)
        rows.append(report)
        direct=report['generation']['direct_C5_candidate_quartics']
        print(f"{model}: lambdaHHdagS2S2={direct['lambdaHHdagS2S2']['rgbeta_at_muS']}, "
              f"lambdaHHdagS1barS1bar={direct['lambdaHHdagS1barS1bar']['rgbeta_at_muS']}, "
              f"|mu2|/min(m2)={report['generation']['quadratic_fraction_of_min_diagonal_m2']:.3g}; "
              f"UV/EFT1 nfev={report['uv_solver_nfev']}/{report['eft1_solver_nfev']}")
        if args.c5_preview:
            print(f"  UNCERTIFIED one-scale epsilon_C5={report['one_scale_C5_preview']['epsilon_C5_frobenius_uncertified']:.6g}")
    args.outdir.mkdir(parents=True,exist_ok=True)
    (args.outdir/'summary.md').write_text(summary_md(rows),encoding='utf-8')
    save_json(args.outdir/'summary.json',{
        'format':'t3-bc-zero-extra-boundary-benchmark-summary-v1',
        'models':[x['model'] for x in rows],
        'all_zero_UV_boundaries_verified':True,
        'number_of_models':len(rows),
        'C5_preview_performed':args.c5_preview,
        'physical_C5_error_bound':False,
        'model_results':{x['model']:{'file':f"{x['model'].replace('-','m')}/result.json",
                          'sha256':digest(args.outdir/x['model'].replace('-','m')/'result.json')}
                         for x in rows},
    })
    print('Batch 7 result:',args.outdir/'summary.md')


if __name__ == '__main__':
    main()
