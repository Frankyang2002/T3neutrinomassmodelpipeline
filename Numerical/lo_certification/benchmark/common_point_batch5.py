"""Batch 5: source-controlled A_-4 UV loop versus fitted sequential conditional EFT.

All numerical physics kernels are imported unchanged from the repository:
BuildT3DirectOneLoopInputs.construct, direct_uv.independent_reference,
Numerical.direct_one_loop.core.evaluate and the Batch 4 source-first chain.

No certified GammaCC, quartic dictionary, hard F matching or equivalence of
UV full-theory and EFT fixed-order content is claimed.
"""
from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
from pathlib import Path
from typing import Any
import json
import math
import re
import sys
import shutil
import tempfile

import numpy as np

from Numerical.lo_certification.benchmark import common_input_batch4 as b4

ROOT = Path(__file__).resolve().parents[3]
OWNED = Path('output/general_observable/lo_certification_benchmark')
INPUT_UV = Path('output/general_observable/lo_direct_one_loop/provisional_inputs/A_-4.json')
AUDIT = Path('output/general_observable/lo_direct_one_loop/source_audits')
REF_UV = b4.UV_REFERENCE
POTENTIAL_KEYS = ('mS1Sq', 'mS2Sq', 'lambdaH1', 'lambdaH2', 'lambdaT3')
UNSUPPORTED = ('lambdaH1Adj', 'lambdaH2Adj', 'lambdaS1Adj', 'lambdaS2Adj',
               'lambdaHHdagS2S2','lambdaHHdagS1barS1bar',
               'lambdaHHdagS1barS2barCross')
FORMAT='t3-a-m4-common-input-source-controlled-uv-versus-eft-v1'

class Batch5Error(ValueError):
    pass

def require(ok: Any, why: str) -> None:
    if not bool(ok):
        raise Batch5Error(why)

def read(path: Path) -> dict:
    require(path.is_file(), f'Missing file: {path}')
    val=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(val,dict),f'Not JSON object: {path}')
    return val

def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()

def sha256_lf(path: Path) -> str:
    """Only line-ending normalization; never normalise JSON semantics."""
    content=path.read_bytes().replace(b'\r\n',b'\n')
    return sha256(content).hexdigest()

def z3(obj: dict, label: str) -> np.ndarray:
    return b4.z_matrix(obj, label)

def from_encoded(raw: Any, label: str) -> Any:
    if isinstance(raw,dict) and set(raw)=={'real','imag'}:
        a=np.asarray(raw['real'], dtype=float)+1j*np.asarray(raw['imag'], dtype=float)
        require(np.isfinite(a).all(),f'Non-finite {label}')
        return complex(a) if a.ndim==0 else a
    if isinstance(raw,list):
        return np.asarray([[complex(e['real'],e['imag']) for e in row]
                            for row in raw],complex)
    require(type(raw) in (float,int) and math.isfinite(raw),f'Invalid real {label}')
    return raw

def to_frozen_style(raw: Any, name: str, exemplar: Any) -> Any:
    """Convert encoded RGE state into the original frozen JSON value shape."""
    val=from_encoded(raw,name)
    if isinstance(exemplar,dict) and isinstance(exemplar.get('real'),list):
        require(isinstance(val,np.ndarray) and val.shape==(3,3),f'{name}: shape drift')
        return b4.z_json(val)
    if isinstance(exemplar,dict):
        require(isinstance(val,complex),f'{name}: complex scalar expected')
        return {'real':float(val.real),'imag':float(val.imag)}
    require(np.ndim(val)==0 and abs(complex(val).imag)<=1e-13,f'{name}: real scalar expected')
    return float(complex(val).real)

def assert_fitted_running(frozen: dict, running: dict, b4summary:dict) -> None:
    require(running.get('model')==b4.MODEL and running.get('status')=='NumericalRenormalisableRunningOnly',
            'Batch 4 running schema/model mismatch')
    require(b4summary.get('model')==b4.MODEL and b4summary.get('format')==b4.FORMAT,
            'Not the completed fitted Batch 4 report')
    require(b4summary.get('status')=='COMMON_FITTED_INPUT_CONDITIONAL_LLSS_SCALAR_RESPONSE',
            'Batch 4 not numerically complete')
    require(b4summary.get('source_frozen_reference_sha256')==b4.BENCH_SHA and
            running.get('source_benchmark_sha256')==b4.BENCH_SHA,
            'Batch 4 source benchmark changed')
    require(running.get('uv_solver',{}).get('success') is True and
            running.get('eft1_solver',{}).get('success') is True,
            'Batch 4 RGE unsuccessful')
    require(running.get('scales_gev')=={'UV':1e7,'F':1e5,'S':1e3},
            'Batch 4 matching scales changed')
    source={**frozen['base_state']['sm'],**frozen['base_state']['ordinary']}
    require(set(running['uv_initial'])==set(source),'Batch 4 initial keys changed')
    for key,raw in source.items():
        a=from_encoded(running['uv_initial'][key],key)
        b=from_encoded(raw,key)
        require(np.array_equal(np.asarray(a),np.asarray(b)),f'Batch 4 initial state drift: {key}')
    for name in ('MF','y1','y2'):
        require(name in running['uv_at_fermion_threshold'],f'Missing muF {name}')
    w=z3(b4summary['F_tree_native_source_and_flavour']['ordered_flavour_matrix_gev_inv'],
         'Batch 4 heavy flavour')
    m=from_encoded(running['uv_at_fermion_threshold']['MF'],'muF MF')
    y1=from_encoded(running['uv_at_fermion_threshold']['y1'],'muF y1')
    y2=from_encoded(running['uv_at_fermion_threshold']['y2'],'muF y2')
    require(m.shape==y1.shape==y2.shape==(3,3),'Changed F Yukawa shapes')
    from Numerical.matching.GeneralT3FlavourC5 import prepare_heavy_flavour
    from types import SimpleNamespace
    uv=SimpleNamespace(MF=m,y1=y1,y2=y2,general_couplings={},
                       representation=SimpleNamespace(majorana_fermion=False))
    mass,yuk,_=prepare_heavy_flavour(uv,{'y1','y2'})
    reconstructed=b4.invariant_flavour(yuk['y1'],yuk['y2'],mass)
    relative=float(np.linalg.norm(reconstructed-w)/max(np.linalg.norm(w),1e-100))
    require(relative<1e-9,f'MuF heavy flavour differs from Batch 4: {relative}')
    require(b4summary['physical_certification']['physical_C5_certified'] is False,
            'Unexpected physical certification promotion')

def state_benchmark(frozen:dict, running:dict, *, scalar_stage:str) -> dict:
    """F-evolved y/M, with scalar parameters from UV, muF or muS.

    This is a numerical re-evaluation; it does not claim a full-theory RG
    trajectory after decoupling F or an exact scalar EWSB dictionary.
    """
    require(scalar_stage in ('UV','F','S'),'Unknown scalar stage')
    b=deepcopy(frozen)
    atf=running['uv_at_fermion_threshold']
    for key in ('y1','y2','MF'):
        b['base_state']['ordinary'][key]=to_frozen_style(atf[key],key,
                                        frozen['base_state']['ordinary'][key])
    source={'UV':frozen['base_state']['ordinary'],
            'F':atf,
            'S':running['eft1_at_scalar_threshold']}[scalar_stage]
    for key in POTENTIAL_KEYS:
        require(key in source,f'Scalar stage {scalar_stage} missing {key}')
        b['base_state']['ordinary'][key]=to_frozen_style(source[key],key,
                                          frozen['base_state']['ordinary'][key])
    # The source direct-UV constructor cannot include adjoint Higgs portals.
    # Since the benchmark begins with zero adjoint values, ONLY the documented
    # restricted projection is made and the removed values are recorded.
    omitted={}
    if scalar_stage!='UV':
        for key in UNSUPPORTED:
            if key in source:
                value=from_encoded(source[key],key)
                z=complex(value)
                if abs(z)>0:
                    omitted[key]={'real':z.real,'imag':z.imag}
            if key in b['base_state']['ordinary']:
                b['base_state']['ordinary'][key]=0.0
    b['base_state']['mu_gev']=1e5 if scalar_stage!='S' else 1e3
    return {'benchmark':b,'omitted_portals':omitted,
            'scalar_parameter_scale_gev':{'UV':1e7,'F':1e5,'S':1e3}[scalar_stage],
            'fermion_parameter_scale_gev':1e5}

def check_audits(root:Path, original:dict) -> tuple[dict,dict,dict,dict]:
    """Native CG references must match the audit names/hashes in original input."""
    fields=(('component_support_audit.json','cg_support_sha256'),
            ('matchete_generator_analysis.json','generator_analysis_sha256'))
    recorded=original.get('provenance',{})
    require(recorded.get('benchmark_sha256')==b4.BENCH_SHA,
            'Original UV input was not derived from fitted source')
    require(recorded.get('cg_generator_covariance_checked') is True and
            recorded.get('validated_uv_component_conventions') is False,
            'Original UV component certificate status changed')
    paths={}
    for filename,key in fields:
        p=root/AUDIT/filename
        require(p.is_file(),f'Missing native CG source audit {p}')
        require(recorded.get(key) in {digest(p),sha256_lf(p)},
                f'Native CG audit drift: {p}')
        paths[filename]=p
    cp,gp=paths['component_support_audit.json'],paths['matchete_generator_analysis.json']
    cv=root/AUDIT/'cg_covariance_audit.json'
    require(cv.is_file(),f'Missing CG covariance audit: {cv}')
    cases,groups,cov=read(cp),read(gp),read(cv)
    require(len(cases.get('cases',[]))==16 and not cases.get('errors'),
            'Component audit does not cover 16 models')
    require(cov.get('status')=='cg_generator_covariance_verified',
            'Source CG Ward covariance certification missing')
    label=(1,3,2,-4)
    found=[x for x in cases['cases'] if tuple(int(x['model'][k]) for k in
                        ('dS1','dS2','dF','alpha'))==label]
    good=[x for x in cov['models'] if tuple(int(x['model'][k]) for k in
                        ('dS1','dS2','dF','alpha'))==label]
    require(len(found)==len(good)==1 and good[0].get('all_covariant') is True,
            'No verified A_-4 CG tensors')
    # The native builder imports its sibling RunT3DirectOneLoop as a
    # top-level module; make that existing source directory importable also
    # from unittest/module entry points (not only direct scripts/ execution).
    script_dir=str(ROOT/'scripts')
    if script_dir not in sys.path:
        sys.path.insert(0,script_dir)
    from scripts.BuildT3DirectOneLoopInputs import _rotations
    return found[0],_rotations(groups),{'cp':cp,'gp':gp,'cv':cv},cov

def encode_matrix(a:np.ndarray)->dict:
    return {'real':np.asarray(a).real.tolist(),'imag':np.asarray(a).imag.tolist()}

def comparison(uv_gev:np.ndarray, eft_perN_ev:np.ndarray)->dict:
    """Fit one COMPLEX scalar; no equality claim for unmatched field content."""
    u=np.asarray(uv_gev,complex)*1e9
    e=np.asarray(eft_perN_ev,complex)
    require(u.shape==e.shape==(3,3) and np.isfinite(u).all() and np.isfinite(e).all(),
            'UV/EFT matrices invalid')
    require(np.linalg.norm(u)>0 and np.linalg.norm(e)>0,'Zero comparison matrix')
    n=np.vdot(e,u)/np.vdot(e,e)
    residual=float(np.linalg.norm(u-n*e)/np.linalg.norm(u))
    norm_e=float(np.linalg.norm(e))
    norm_u=float(np.linalg.norm(u))
    out={
        'best_complex_multiplier_of_EFT_per_N_to_UV':{'real':float(n.real),'imag':float(n.imag)},
        'residual_after_one_complex_multiplier':residual,
        'UV_masses_eV_descending':np.linalg.svd(u,compute_uv=False).tolist(),
        'EFT_masses_eV_per_abs_N_descending':np.linalg.svd(e,compute_uv=False).tolist(),
        'offdiagonal_uv_norm_fraction':float(np.linalg.norm(u-np.diag(np.diag(u)))/norm_u),
        'offdiagonal_eft_norm_fraction':float(np.linalg.norm(e-np.diag(np.diag(e)))/norm_e),
        'UV_frobenius_norm_eV':norm_u,
        'EFT_frobenius_norm_eV_per_abs_N':norm_e,
        'matrix_agreement_certified':False,
    }
    return out

def original_reproduction(root:Path,frozen:dict,case:dict,us:dict,paths:dict,
                          original:dict, reference:dict)->tuple[dict,dict]:
    from scripts.BuildT3DirectOneLoopInputs import construct
    from Numerical.lo_certification.direct_uv.native_uv_loop_reference import independent_reference,decode
    from Numerical.direct_one_loop.core import evaluate
    got=construct(case,frozen,us,'A_-4',root/b4.BENCH,paths['cp'],paths['gp'],True)
    for key in ('fermion_mass_gev','scalar_A_gev2','scalar_B_gev2'):
        x=decode(got[key],key)
        y=decode(original[key],key)
        rel=np.linalg.norm(x-y)/max(np.linalg.norm(y),1e-90)
        require(rel<3e-13,f'Original CG constructor fails to reproduce {key}: {rel}')
    for r1,r2 in zip(got['vertices'],original['vertices']):
        require(r1['name']==r2['name'],'Original vertex inventory changed')
        k='field' if 'field' in r1 else 'conjugate'
        x=decode(r1[k],k);y=decode(r2[k],k)
        rel=np.linalg.norm(x-y)/max(np.linalg.norm(y),1e-90)
        require(rel<3e-12,f'Original CG constructor fails to reproduce {r1["name"]}: {rel}')
    a=independent_reference(got,check_rephase=False)
    gotm=decode(a['complex_symmetric_matrix_gev'],'reconstructed original UV')
    refm=decode(reference['complex_symmetric_matrix_gev'],'published independent UV')
    delta=float(np.linalg.norm(gotm-refm)/max(np.linalg.norm(refm),1e-90))
    require(delta<3e-6,f'Original 1e7 GeV direct UV mass no longer reproducible: {delta}')
    old=evaluate(got,retain_individual_terms=False)
    alt=decode(old['total']['complex_symmetric_matrix_gev'],'original UV production')
    err=float(np.linalg.norm(gotm-alt)/max(np.linalg.norm(gotm),1e-90))
    require(err<3e-7,f'Independent/production UV methods disagree: {err}')
    return a,{'reference_mass_relative_difference':delta,
              'independent_vs_production_reference_relative':err,
              'uv_reference_masses_eV_descending':np.linalg.svd(refm,compute_uv=False).tolist()}

def evaluate_scenario(root:Path,case:dict,us:dict,paths:dict,frozen:dict,
                      running:dict, stage:str, eft:np.ndarray)->tuple[dict,dict]:
    from scripts.BuildT3DirectOneLoopInputs import construct
    from Numerical.lo_certification.direct_uv.native_uv_loop_reference import independent_reference,decode
    from Numerical.direct_one_loop.core import evaluate
    package=state_benchmark(frozen,running,scalar_stage=stage)
    payload=construct(case,package['benchmark'],us,'A_-4',root/b4.BENCH,
                      paths['cp'],paths['gp'],True)
    payload['provenance']={**payload['provenance'],
       'derived_from_evolved_running':True,'frozen_original_benchmark_sha256':b4.BENCH_SHA,
       'fermions_evolved_to_muF':True,'scalar_stage':stage,
       'scalar_stage_scale_gev':package['scalar_parameter_scale_gev'],
       'not_a_single_running_scale_for_all_parameters':stage!='F',
       'unsupported_evolved_scalar_couplings_projected_to_zero':package['omitted_portals'],
       'native_component_mapping_physically_certified':False}
    payload['provenance']['input_assumptions']=[
       'An explicitly labelled numerical comparison only',
       'All heavy F masses/Yukawas from Batch4 muF; scalar masses/quartics from '+stage,
       'Unsupported evolved adjoint Higgs portals projected to zero when present',
       'Native CG phase, scalar-potential map and fermion hard matching unverified']
    independent=independent_reference(payload,check_rephase=False)
    production=evaluate(payload,retain_individual_terms=False)
    a=decode(independent['complex_symmetric_matrix_gev'],'independent UV')
    z=decode(production['total']['complex_symmetric_matrix_gev'],'production UV')
    d=float(np.linalg.norm(a-z)/max(np.linalg.norm(a),np.linalg.norm(z),1e-95))
    require(d<4e-7,f'UV independent/production numerical disagreement at {stage}: {d}')
    require(production['diagnostics']['uv_finite'] is True,'UV pole does not cancel')
    require(independent['diagnostics']['uv_pole_relative_to_natural_amplitude']<1e-9,
            'Independent UV pole failed')
    masses=np.asarray(independent['scalar_masses_gev'])
    heavy=np.asarray(independent['fermion_mass_singular_values_gev'])
    require(np.min(heavy)>np.max(masses),'Invalid hierarchy for UV heavy expansion')
    metrics=comparison(a,eft)
    metrics.update({'scalar_stage':stage,'fermion_scale_gev':1e5,
                    'scalar_parameters_scale_gev':package['scalar_parameter_scale_gev'],
                    'not_a_single_scale_full_UV_input':stage!='F',
                    'unsupported_evolved_scalar_couplings_dropped':package['omitted_portals'],
                    'num_unsupported_scalar_couplings_dropped':len(package['omitted_portals']),
                    'max_scalar_over_min_fermion':float(np.max(masses)/np.min(heavy)),
                    'independent_vs_production_relative_discrepancy':d,
                    'pole_relative':independent['diagnostics']['uv_pole_relative_to_natural_amplitude']})
    return metrics,{'input':payload,'independent':independent,'production':production}

def preflight(root:Path, run_id:str)->dict:
    root=root.resolve()
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,50}',run_id) is not None,'Unsafe run ID')
    folder=root/OWNED/('batch4_common_input_'+run_id)
    frozen,_,_,_=b4.check_source(root)
    run=read(folder/'fitted_running.json'); b4result=read(folder/'summary.json')
    assert_fitted_running(frozen,run,b4result)
    original=read(root/INPUT_UV); reference=read(root/REF_UV)
    _,_,paths,_=check_audits(root,original)
    require(reference.get('input_sha256') in {digest(root/INPUT_UV),sha256_lf(root/INPUT_UV)},
            'Original input does not match independent UV reference')
    eft=z3(b4result['conditional_mnu_eV_per_N'],'Batch4 EFT mass')
    require(np.max(np.abs(eft))>0 and len(b4result['conditional_neutrino_masses_eV_per_abs_N'])==3,
            'Batch4 EFT mass missing')
    stored=np.asarray(b4result['conditional_neutrino_masses_eV_per_abs_N'],float)
    recalculated=np.linalg.svd(eft,compute_uv=False)
    require(np.isfinite(stored).all() and
            np.linalg.norm(recalculated-stored)<1e-10*max(np.linalg.norm(stored),1e-100),
            'Batch4 EFT singular values do not match its complex matrix')
    boundary=read(folder/'A_m4_fitted_LLSS_provisional.json')
    require(boundary.get('model')=='A_-4' and
            boundary.get('spinor_normalization',{}).get('certified') is False and
            boundary.get('spinor_normalization',{}).get('gamma_cc_to_master_weyl_factor') is None,
            'Batch4 regenerated LLSS incorrectly promoted to physical')
    rows=boundary.get('summed_independent_symmetric_coordinates')
    require(isinstance(rows,list) and len(rows)>0 and
            sum(r['coordinate'][1]!=r['coordinate'][3] for r in rows)==
            b4result['F_tree_native_source_and_flavour']['new_offdiagonal_flavour_coordinate_count'],
            'Batch4 LLSS offdiagonal support differs from summary')
    scalar=read(folder/'conditional_scalar_matching.json')
    require(scalar.get('physical_C5_certified') is False and
            scalar.get('response',{}).get('physical_C5_certified') is False,
            'Batch4 conditional scalar match incorrectly certified')
    c5=z3(scalar['response']['C5_LLSS_total_one_loop_per_N_gev_inverse'],'Batch4 scalar C5')
    expected=c5*(246.22**2/2)*1e9
    require(np.linalg.norm(expected-eft)<2e-12*max(np.linalg.norm(eft),1e-100),
            'Batch4 summary and original joint scalar C5 are not consistent')
    return {'model':'A_-4','status':'PREFLIGHT_OK_NO_OUTPUT_WRITTEN',
            'batch4_path':str(folder),'benchmark_scale_gev':1e7,
            'fermion_threshold_gev':1e5,'joint_scalar_threshold_gev':1e3,
            'fitted_input_verified':True,
            'native_CG_source_files':[str(v) for v in paths.values()],
            'UV_comparison_not_physical_certification':True,'output':None}

def calculate(root:Path, batch4_run_id:str, run_id:str)->dict:
    root=root.resolve()
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,50}',run_id) is not None,'Unsafe run ID')
    destination=root/OWNED/('batch5_uv_eft_'+run_id)
    require(not destination.exists(),f'Batch5 output exists: {destination}')
    info=preflight(root,batch4_run_id)
    frozen,_,_,_=b4.check_source(root)
    folder=root/OWNED/('batch4_common_input_'+batch4_run_id)
    running=read(folder/'fitted_running.json')
    eftdoc=read(folder/'summary.json')
    eft=z3(eftdoc['conditional_mnu_eV_per_N'],'Batch4 EFT mass')
    original=read(root/INPUT_UV)
    reference=read(root/REF_UV)
    case,us,paths,_=check_audits(root,original)
    origcheck,repro=original_reproduction(root,frozen,case,us,paths,original,reference)
    scenarios={}; payloads={}
    for stage in ('UV','F','S'):
        metrics,records=evaluate_scenario(root,case,us,paths,frozen,running,stage,eft)
        scenarios[stage]=metrics
        payloads[stage]=records
    summary={
       'format':FORMAT,'status':'CONDITIONAL_NUMERIC_COMPARISON_NOT_SAME_PHYSICAL_MATCHING',
       'model':'A_-4','input_at_muUV':'fitted frozen source; all 20 couplings',
       'source_shas':{'fitted_benchmark':digest(root/b4.BENCH),
                      'original_direct_uv_input':digest(root/INPUT_UV),
                      'published_uv_reference':digest(root/REF_UV),
                      'batch4_running':digest(folder/'fitted_running.json'),
                      'batch4_summary':digest(folder/'summary.json'),
                      'native_CG_component_audit':digest(paths['cp']),
                      'native_generators':digest(paths['gp']),
                      'native_covariance':digest(paths['cv'])},
       'original_direct_uv_reproduction':repro,
       'eft_conditional_mass_per_N_eV':b4.z_json(eft),
       'eft_conditional_masses_eV_per_abs_N':np.linalg.svd(eft,compute_uv=False).tolist(),
       'UV_scenarios':scenarios,
       'interpretation':{
          'UV':'evolved F-sector at muF; original frozen UV scalar sector',
          'F':'evolved F and scalar sector at muF, unsupported portals projected',
          'S':'muF F-sector with muS EFT1 scalar sector, unsupported portals projected; hybrid inputs',
          'UV_mass_is_not_C5_and_not_same_matching_order':True,
          'F_full_EWSB_mass_is_not_unbroken_phase_EFT_C5':True,
          'EFT_includes_only_LLSS_induced_log_plus_one_simultaneous_scalar_bubble':True,
          'all_BSM_fermions_integrated_first':True,
          'both_scalars_integrated_together':True,
          'full_hard_one_loop_F_matching_available':False,
          'GammaCC_signed_factor':None,
          'Matchete_RGBeta_quartic_map_certified':False,
          'physical_neutrino_masses_certified':False},
       'projected_scalar_operator_terms_not_a_complete_full_theory':True,
       'no_gammaCC_factor_fitted_as_physical':True,
    }
    destination.parent.mkdir(parents=True,exist_ok=True)
    staging=Path(tempfile.mkdtemp(prefix='.batch5_staging_',dir=destination.parent))
    try:
        (staging/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        lines=['# Batch 5: A_-4 fitted-input UV/EFT numerical diagnostics',
               '', '**All UV amplitudes are provisional, not physical matching.**', '',
               'UV frozen source was reproduced before substituting any evolved input.',
               f"Reproduction relative matrix difference: {repro['reference_mass_relative_difference']:.3e}.",
               '', '| Scalar inputs | UV singular values (eV) | Best-scale residual | Extra portals omitted |',
               '|---|---|---:|---:|']
        for stage in ('UV','F','S'):
            x=scenarios[stage]
            lines.append(f"| {stage} | "+', '.join(f'{v:.8g}' for v in x['UV_masses_eV_descending'])+
                         f" | {x['residual_after_one_complex_multiplier']:.6f} | {x['num_unsupported_scalar_couplings_dropped']} |")
        lines.extend(['','EFT conditional singular masses per |N_GammaCC| (eV): '+
              ', '.join(f'{v:.9g}' for v in summary['eft_conditional_masses_eV_per_abs_N']),
              '', 'The one-complex-scale fit is descriptive only; it is NOT a physical GammaCC determination.',
              'The F/S scenarios omit unsupported evolved adjoint portals if present;',
              'their omitted values are preserved in summary.json. The UV scenario',
              'uses the original UV scalar sector, not a common scalar parameter scale.',
              'The S scenario is explicitly hybrid (muF heavy fermions, muS scalars).',
              '', 'UV -> integrate out ALL F at muF -> EFT1(S1,S2) -> integrate out BOTH scalars at muS -> SMEFT.',
              'Unresolved: signed GammaCC, full hard fermion matching, complete Wilson running,',
              'native scalar/portal coupling normalisation, full EFT/UV perturbative consistency.', ''])
        (staging/'report.md').write_text('\n'.join(lines),encoding='utf-8')
        for stage,parts in payloads.items():
            for name,value in parts.items():
                (staging/f'{stage}_direct_{name}.json').write_text(
                    json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        (staging/'original_uv_recomputed_reference.json').write_text(
                json.dumps(origcheck,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        require(not destination.exists(),'Refusing to overwrite existing Batch5 output')
        staging.rename(destination)
    except BaseException:
        shutil.rmtree(staging,ignore_errors=True)
        raise
    return {'model':'A_-4','status':summary['status'],'output':str(destination),
            'original_uv_reference_reproduced':True,
            'scenarios':{name:{'residual':x['residual_after_one_complex_multiplier'],
                        'UV_mass_singular_values_eV':x['UV_masses_eV_descending'],
                        'projected_portals':x['num_unsupported_scalar_couplings_dropped']}
                         for name,x in scenarios.items()},
            'physical_C5_certified':False}
