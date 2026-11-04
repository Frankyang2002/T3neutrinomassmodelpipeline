"""All-16 FULL-GENERAL UV/EFT1/Wilson/SMEFT running ablation, no refits.

Four independent switches, 16 predictions per model:
    U: general UV beta on/off (1e7 -> 1e5)
    E: general scalar EFT beta on/off (1e5 -> 1e3)
    W: LLSS-induced intermediate Weinberg leading log on/off
    S: final SM+Weinberg beta on/off (1e3 -> 1e2)

All fermions are integrated at 1e5; BOTH scalar multiplets at 1e3.

The full-general ONE-SCALE Matchete C5 is the reference C5 for W=1. The
independently calculated conditional LLSS leading log and joint-scalar bubble
are subtracted first to define a bookkeeping residual; each is then restored
exactly ONCE. W=0 removes the LLSS leading log without modifying the rest.
This guarantees no double-counting but is NOT an independently calculated
hard matching. Matchete/GammaCC and extra-HHSS basis assumptions are explicit.

14 nonmixed models use the actual muF UV Yukawa matrices and all source terms.
B(-1)/C(-1) use source-exact native-CG Majorana Yukawa channels and U(muS).
No restricted C5 kernels and no GitHub writes.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
import json
import math
from pathlib import Path
import re
import traceback
from types import SimpleNamespace

import numpy as np

from Numerical.running.GeneralT3SequentialLeadingLog import (
    MODELS, encode, decode, require, response, scalar_tensors,
)
from Numerical.running.GeneralT3DynamicFermionWilson import (
    nonmixed_boundary,
)
from scripts.RunT3GeneralSequentialLeadingLog import low_energy, hash_file, load, save, tag

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'output/studies/general_scalar_extended/general_sequential_rge_ablation'
DEFAULT_RELEASE = ROOT/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
DEFAULT_V4 = ROOT/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'
DEFAULT_SEEDS = ROOT/'output/full16_symbolic'
DEFAULT_REGISTRY = ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json'
DEFAULT_REF_RUNNING = ROOT/'output/general_observable/all16_v1'
SCALES=(1e7,1e5,1e3,1e2)


def frozen_state(state, mu: float):
    return replace(state, mu_gev=mu).validated()


def c5_checked(value: np.ndarray) -> np.ndarray:
    z=np.asarray(value,complex)
    require(z.shape==(3,3) and np.isfinite(z).all(),'invalid one-scale Matchete C5')
    return 0.5*(z+z.T)


def snapshot(model: str, state: dict, source_path: Path):
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import RunningSnapshot
    m1=float(state['mS1Sq']);m2=float(state['mS2Sq'])
    require(np.isfinite([m1,m2]).all() and min(m1,m2)>0,
            f'{model}: nonpositive scalar mass')
    return RunningSnapshot(model,source_path,hash_file(source_path),
                    SCALES[1],SCALES[2],m1,m2,state,
                    {'source':'current_general_fitted_RGE_endpoint'})


def _load_channels(model: str, static_source, *, v4: Path, seeds: Path):
    from Numerical.joint_scalar_matching.source_quartic_v4 import source_seed_path
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import read_v4_cached_channels
    return read_v4_cached_channels(v4/(tag(model)+'.json'),
                      source_seed_path(seeds,model),static_source)


def _final_cases(model: str, c5_one: np.ndarray, source, scalar,
                 target, *, uv_on: bool, eft_on: bool,
                 uv_nfev: int, eft_nfev: int,
                 spinor_factor: float, vev: float) -> dict:
    matched=c5_checked(c5_one)
    llss=response(source,scalar,mu_f=SCALES[1],mu_s=SCALES[2],
                  spinor_factor=complex(spinor_factor))
    log=decode(llss['c5_eft1_leading_log'])
    finite=decode(llss['c5_scalar_finite_bubble'])
    # Mathematical one-scale anchor: no independent claim for hard matching.
    residual=matched-(log+finite)
    cases={}
    for w_on in (False,True):
        c5=residual+finite+(log if w_on else 0.)
        c5=c5_checked(c5)
        for sm_on in (False,True):
            key=f'{int(uv_on)}{int(eft_on)}{int(w_on)}{int(sm_on)}'
            record=low_energy(scalar_snapshot_state(scalar),c5,on=sm_on,
                        mu_s=SCALES[2],mu_low=SCALES[3],vev=vev,target=target)
            record.update({'case':key,
                'uv_on':uv_on,'eft1_on':eft_on,'wilson_log_on':w_on,
                'smeft_on':sm_on,'uv_nfev':uv_nfev,'eft1_nfev':eft_nfev,
                'c5_matching_anchored_reference':encode(matched),
                'c5_induced_wilson_log':encode(log),
                'c5_joint_scalar_bubble':encode(finite),
                'c5_bookkeeping_residual_matching':encode(residual),
                'c5_at_scalar_threshold':encode(c5),
                'c5_log_frobenius':float(np.linalg.norm(log)),
                'c5_one_scale_frobenius':float(np.linalg.norm(matched)),
                'c5_log_over_one_scale_norm':float(np.linalg.norm(log)/
                                                 max(np.linalg.norm(matched),1e-40)),
                'per_general_interaction':llss['per_channel'],
                'quartic_provenance':llss['quartic_provenance'],
                'matching_scale_cancellation':llss['matching_scale_cancellation_residual'],
                'weinberg_projection_residual':max(
                    llss['weinberg_beta_projection_residual'],
                    llss['weinberg_bubble_projection_residual']),
                'matching_approximation':'one_scale_minus_LLSS_log_and_bubble_residual',
            })
            cases[key]=record
    return cases


def scalar_snapshot_state(bundle):
    # Set by a transient attribute? TensorBundle is frozen dataclass, so the
    # caller supplies the state via the explicit ScalarEndpoint wrapper below.
    return bundle.couplings


class ScalarEndpoint:
    """Preserve the snapshot alongside tensor arrays, no mutation of physics data."""
    def __init__(self, bundle, couplings):
        self.bundle=bundle
        self.couplings=couplings
        self.model=bundle.model
    def __getattr__(self, name):
        return getattr(self.bundle,name)


def _nonmixed(model: str, saved: Path, canonical_run_id: str,
              quartic_source: tuple, registry: dict,
              *, spinor_factor: float, wolframscript: str, vev: float):
    from scripts import RunT3GeneralAll16FromReference as g
    from Numerical.orchestration.GeneralRunningBenchmark import (
        initial_state, _state_values)
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
    from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.RunSmallLambdaRefit import TARGET

    beta_uv,beta_eft,_,kernel,matchete_path,kind=g.nonmixed_inputs(
        model,canonical_run_id,generate_general_c5=False,wolframscript=wolframscript)
    overrides=load(saved/'general_uv_overrides.json')
    _,initial=initial_state(beta_uv,mu_uv=SCALES[0],mu_f=SCALES[1],
                             mu_s=SCALES[2],overrides=overrides)
    static_source,channels=quartic_source
    target=OscillationFitTarget.from_json(TARGET)
    from Numerical.lo_certification.couplings.extra_hhss_geometry import model_extra_geometry
    geometry=model_extra_geometry(model,registry['models'][model],include_components=True)
    cases={};sources={};scalars={}
    for uv_on in (False,True):
        if uv_on:
            uv_sol=run_uv_segment(initial,beta_uv,SCALES[1],rtol=1e-7,atol=1e-10)
            uv_end,uv_nfev=uv_sol.final_state,int(uv_sol.nfev)
        else:
            uv_end,uv_nfev=frozen_state(initial,SCALES[1]),0
        tree,tree_diagnostics=nonmixed_boundary(model,uv_end,repo_root=ROOT,
                                                  mu_f=SCALES[1])
        sources[str(int(uv_on))]=tree_diagnostics
        eft_start=build_intermediate_scalar_boundary(uv_end,matching_scale_gev=SCALES[1])
        _,layout=pack_intermediate_scalar_state(eft_start)
        require(set(layout.names)==set(beta_eft.get('report_betas',{})),
                f'{model}: general EFT1 state/beta inventory differs')
        for eft_on in (False,True):
            if eft_on:
                eft_sol=run_intermediate_scalar_segment(eft_start,beta_eft,
                              SCALES[2],rtol=1e-7,atol=1e-10)
                eft_end,eft_nfev=eft_sol.final_state,int(eft_sol.nfev)
            else:
                eft_end,eft_nfev=frozen_state(eft_start,SCALES[2]),0
            state=_state_values(eft_end)
            snap=snapshot(model,state,saved/'fitted_config.json')
            tensor=scalar_tensors(model,tree,channels,snap,registry,geometry=geometry)
            scalar=ScalarEndpoint(tensor,state)
            c5=kernel.evaluate(uv_end,eft_end)
            block=_final_cases(model,c5,tree,scalar,target,
                uv_on=uv_on,eft_on=eft_on,uv_nfev=uv_nfev,eft_nfev=eft_nfev,
                spinor_factor=spinor_factor,vev=vev)
            cases.update(block)
    return cases,{'fermion_tree':sources,'matching_kind':kind,
                  'full_general_C5_matchete_path':str(matchete_path),
                  'full_general_C5_sha256':hash_file(matchete_path),
                  'fitted_uv_sha256':hash_file(saved/'general_uv_overrides.json')}


def mixed_fit_fractions(fit: dict) -> tuple[float, float]:
    """Read B/C boundary fractions from either original-fit or continuation metadata.

    Original B/C fitting records keep the fractions at the top level. A
    continuation changes fitted Yukawas but preserves the scalar boundary;
    its saved provenance calls these fractions quartic_fraction/mu2_fraction.
    Never silently substitute a new scalar boundary when metadata is absent.
    """
    direct = ('extra_quartic_fraction', 'extra_mixed_scalar_bilinear_fraction')
    continuation = ('quartic_fraction', 'mu2_fraction')
    meta = fit.get('matching_and_rge_provenance', {})
    require(isinstance(meta, dict), 'invalid mixed-fit matching provenance')
    if meta.get('modified_scalar_boundary') is True:
        raise ValueError('continuation changed scalar boundary; cannot reuse original fractions')

    values=[]
    for original, carried in zip(direct, continuation):
        has_direct = original in fit
        has_carried = carried in meta
        require(has_direct or has_carried,
                f'missing mixed-fit boundary fraction {original!r} or '
                f'matching_and_rge_provenance.{carried}')
        value = float(fit[original] if has_direct else meta[carried])
        require(math.isfinite(value) and value >= 0,
                f'invalid mixed-fit boundary fraction {original}: {value}')
        if has_direct and has_carried:
            other = float(meta[carried])
            require(math.isfinite(other) and math.isclose(value, other,
                    rel_tol=1e-12, abs_tol=1e-15),
                    f'conflicting saved mixed-fit boundary fractions: {original} vs {carried}')
        values.append(value)
    return values[0], values[1]


def _mixed(model: str, saved: Path, quartic_source: tuple,
           reference_running: Path, registry: dict,
           *, spinor_factor: float, vev: float):
    from scripts import RunT3BCGeneralFromReference as runner
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc
    from scripts import RunT3GeneralScalarBCPilot as pilot
    from scripts import RunT3ZeroExtraBoundaryBenchmark as core
    from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE
    from Numerical.experimental.GeneralT3ScalarRGState import (
        unpack_general_scalar_state, project_general_after_fermion,
        pack_general_scalar_state)
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.RunSmallLambdaRefit import TARGET
    from Numerical.orchestration.GeneralRunningBenchmark import _state_values
    from Numerical.running.GeneralT3MixedExactFermionWilson import mixed_exact_boundary
    from Numerical.running.GeneralT3ConditionalBCQuarticMap import (
        apply_conditional_bc_direct_quartic_map,
    )
    from Numerical.running.GeneralT3ConditionalBCQuarticNormalization import (
        apply_conditional_bc_channel_normalization,
    )

    fit=load(saved/'result.json')
    source=runner.build_reference_source(model)
    source['cfg']=load(saved/'fitted_config.json')
    registry_legacy,new,reused=source['registry'],source['new'],source['reused']
    q,m=mixed_fit_fractions(fit)
    beta_bundle=bc.beta_bundle(model,source,registry_legacy)
    coefficient=Path(DEFAULT_SOURCE)/({'B_-1':'B_m1_final',
                                       'C_-1':'C_m1_final'}[model]+
              '_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    compiled=bc.compile_source(coefficient)
    target=OscillationFitTarget.from_json(TARGET)
    initial,uv_layout,_,_=pilot.prepare_initial(source,registry_legacy,
                    'general_quartics_and_mixing' if m>0 else 'general_quartics',q,m)
    static_source,channels=quartic_source
    from Numerical.lo_certification.couplings.extra_hhss_geometry import model_extra_geometry
    geometry=model_extra_geometry(model,registry['models'][model],include_components=True)
    cases={};source_details={}
    for uv_on in (False,True):
        if uv_on:
            sol=core.integrate_segment(initial,
                  lambda t,x:core.rhs_at(t,x,layout=uv_layout,
                                        registry=registry_legacy,
                                        betas=beta_bundle['UV'],stage='UV'),
                  SCALES[0],SCALES[1],n_points=3)
            uv_vec,uv_nfev=sol.y[:,-1].copy(),int(sol.nfev)
        else:
            uv_vec,uv_nfev=initial.copy(),0
        heavy=unpack_general_scalar_state(uv_vec[:uv_layout.size],uv_layout,
                         registry_legacy,uv=True,mu_gev=SCALES[1]).baseline
        eft=project_general_after_fermion(
              unpack_general_scalar_state(uv_vec[:uv_layout.size],uv_layout,
                     registry_legacy,uv=True,mu_gev=SCALES[1]),registry_legacy)
        eft_vec,eft_layout=pack_general_scalar_state(eft,registry_legacy)
        eft_vec=np.concatenate((eft_vec,uv_vec[uv_layout.size:]))
        eft_beta=core.check_payload(beta_bundle['payloads']['EFT1'],
                   model=model,stage='EFT1',rep=eft_layout.representation,
                   base_names=tuple(eft_layout.names[:-len(new)]),
                   new_names=new,reused_names=reused)
        for eft_on in (False,True):
            if eft_on:
                sol=core.integrate_segment(eft_vec,
                     lambda t,x:core.rhs_at(t,x,layout=eft_layout,
                       registry=registry_legacy,betas=eft_beta,stage='EFT1'),
                       SCALES[1],SCALES[2],n_points=3)
                end_vec,eft_nfev=sol.y[:,-1].copy(),int(sol.nfev)
            else:
                end_vec,eft_nfev=eft_vec.copy(),0
            eft_endpoint=unpack_general_scalar_state(end_vec[:eft_layout.size],
                    eft_layout,registry_legacy,uv=False,mu_gev=SCALES[2]).baseline
            state=_state_values(eft_endpoint)
            # The native B/C extra scalar variables are carried separately in
            # the current stage runner, not in the baseline dataclass fields.
            endpoint=core.snapshot(end_vec,eft_layout,registry_legacy,
                     stage='EFT1',mu_gev=SCALES[2],new=new,reused=reused)
            for name,value in endpoint['new_quartics'].items():
                if name not in state: state[name]=value
            for name,value in endpoint['reused_quartics'].items():
                if name not in state: state[name]=value
            state['muT3GeneralS1S2']=endpoint['quadratic_mu2_gev2']
            # Mixed F-tree source depends on U(muS); rebuild per U/E endpoint.
            tree,tree_diag=mixed_exact_boundary(
                model,heavy,state,repo_root=ROOT,mu_f=SCALES[1])
            source_details[f'{int(uv_on)}{int(eft_on)}']=tree_diag
            snap=snapshot(model,state,saved/'fitted_config.json')
            tensor=scalar_tensors(model,tree,channels,snap,registry,geometry=geometry)
            # Match the existing B/C one-scale conditional signed direct-quartic
            # convention ONLY in HHSS contraction, never in the RG trajectory.
            tensor=apply_conditional_bc_direct_quartic_map(model,tensor)
            # Independent B/C mass-log check fixes channel magnitudes, not a
            # fitted global scale: topology sqrt(3), signed direct channels 2.
            tensor=apply_conditional_bc_channel_normalization(model,tensor)
            scalar=ScalarEndpoint(tensor,state)
            heavy_info=bc.heavy_flavour_at_muF(heavy)
            c5=bc.c5_at_endpoint(model,endpoint,heavy_info,compiled)
            block=_final_cases(model,c5,tree,scalar,target,
                uv_on=uv_on,eft_on=eft_on,uv_nfev=uv_nfev,eft_nfev=eft_nfev,
                spinor_factor=spinor_factor,vev=vev)
            cases.update(block)
    return cases,{'mixed_native_CG_exact_majorana_source':source_details,
                  'C5_source_path':str(coefficient),
                  'C5_source_sha256':hash_file(coefficient),
                  'fitted_config_sha256':hash_file(saved/'fitted_config.json')}


def preflight(model: str, saved: Path, args, canonical_run_id: str) -> dict:
    from Numerical.joint_scalar_matching.source_quartic_v4 import source_seed_path
    inputs={'fitted_result':saved/'result.json',
            'fitted_config':saved/'fitted_config.json',
            'source_release':args.release_dir/tag(model)/'boundary.json',
            'native_quartic_channels':args.v4_channel_dir/(tag(model)+'.json'),
            'native_quartic_seed':source_seed_path(args.quartic_seed_dir,model),
            'full_general_quartic_registry':args.registry}
    issues={}
    if model not in ('B_-1','C_-1'):
        from scripts import RunT3GeneralAll16FromReference as general
        from scripts.RunT3GeneralWilsonBoundary import _first_threshold
        from Numerical.running.GeneralT3DynamicFermionWilson import F_ROOTS
        inputs['general_uv_overrides']=saved/'general_uv_overrides.json'
        inputs['general_matchete_C5']=general.general_coefficient_path(model)
        try:
            status,source,cert=_first_threshold(model,tuple(ROOT/x for x in F_ROOTS))
            inputs['current_F_tree_source']=source
            inputs['current_F_tree_certificate']=cert
        except (FileNotFoundError,ValueError,KeyError) as e:
            issues['current_general_F_threshold']=str(e)
    else:
        from scripts import RunT3ZeroExtraBoundaryBenchmark as core
        from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE
        from Numerical.running.GeneralT3MixedExactFermionWilson import SOURCE_ROOT, source_tag
        native=ROOT/SOURCE_ROOT/source_tag(model)
        inputs['mixed_exact_source_flavour']=native/'source_flavour.json'
        inputs['mixed_exact_source_CG']=native/'mixed_native_tree_llss_cg_registry.json'
        inputs['mixed_exact_ordered_reference']=native/'ordered_tensor.json'
        inputs['mixed_general_UV_rgbeta']=core.beta_path(model,'UV')
        inputs['mixed_general_EFT1_rgbeta']=core.beta_path(model,'EFT1')
        inputs['mixed_general_C5_source']=(Path(DEFAULT_SOURCE)/
             ({'B_-1':'B_m1_final','C_-1':'C_m1_final'}[model]+
              '_one_loop_probe_weinberg_ordered_coefficient.inputform.txt'))
    missing={name:str(path) for name,path in inputs.items() if not path.is_file()}
    return {'model':model,'ready':not missing and not issues,
            'missing':missing,'issues':issues,
            'inputs':{name:str(path) for name,path in inputs.items()}}


def write_csv(rows: list[dict], folder: Path):
    fields=['model','status','case','UV','EFT1','Wilson','SMEFT','chi2',
            'm_light_eV','m_middle_eV','m_heavy_eV',
            'log_over_full_C5','delta_m21_sq_ev2','delta_m3l_sq_ev2',
            'sin2_theta12','sin2_theta13','sin2_theta23','error']
    with (folder/'all16_cases.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for row in rows:
            if row['status']!='Computed':
                w.writerow({'model':row['model'],'status':row['status'],
                            'error':row.get('error')});continue
            for case,d in row['cases'].items():
                pred=d['prediction'];m=d['masses_ev_ascending']
                w.writerow({'model':row['model'],'status':row['status'],
                  'case':case,'UV':case[0],'EFT1':case[1],
                  'Wilson':case[2],'SMEFT':case[3],
                  'chi2':d['chi2'],'m_light_eV':m[0],
                  'm_middle_eV':m[1],'m_heavy_eV':m[2],
                  'log_over_full_C5':d['c5_log_over_one_scale_norm'],
                  **pred})


def parse_args(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-run',required=True,
                   help='Saved general fit campaign under thesis_general_all16')
    p.add_argument('--run-id',required=True,help='New output directory name')
    p.add_argument('--models',default=None,
                   help='Comma separated; defaults ALL 16')
    p.add_argument('--release-dir',type=Path,default=DEFAULT_RELEASE)
    p.add_argument('--v4-channel-dir',type=Path,default=DEFAULT_V4)
    p.add_argument('--quartic-seed-dir',type=Path,default=DEFAULT_SEEDS)
    p.add_argument('--registry',type=Path,default=DEFAULT_REGISTRY)
    p.add_argument('--reference-running-dir',type=Path,default=DEFAULT_REF_RUNNING)
    p.add_argument('--wolframscript',default='wolframscript')
    p.add_argument('--spinor-factor',type=float,default=1.,
                   help='Explicit provisional Matchete GammaCC-to-Weyl factor')
    p.add_argument('--vev',type=float,default=246.22)
    p.add_argument('--resume',action='store_true',
                   help='Reuse existing computed cases; retry only failed/missing models')
    p.add_argument('--check-only',action='store_true')
    return p.parse_args(argv)


def main(argv=None):
    from scripts import RunT3GeneralAll16FromReference as general
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import (
        load_release_manifest,load_source)
    from Numerical.joint_scalar_matching.source_quartic_v4 import source_seed_path
    args=parse_args(argv)
    if not (re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.source_run)
         and re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id)):
        raise ValueError('Invalid source/run identifier')
    selected=tuple(x.strip() for x in args.models.split(',') if x.strip()) if args.models else MODELS
    require(bool(selected) and set(selected)<=set(MODELS) and len(set(selected))==len(selected),
            'invalid model selection')
    require(np.isfinite([args.vev,args.spinor_factor]).all() and args.vev>0,
            'invalid vev/normalization')
    source_root=general.OUTPUT/args.source_run
    require(source_root.is_dir(),f'missing general fit source {source_root}')
    settings=load(source_root/'settings.json')
    canonical_run_id=settings.get('canonical_run_id','full16_general_v1')
    report=[preflight(m,source_root/tag(m),args,canonical_run_id) for m in selected]
    if args.check_only:
        for r in report:
            print(r['model'],'Ready' if r['ready'] else 'Missing',r['missing'],r['issues'])
        return 0 if all(x['ready'] for x in report) else 2
    if not args.registry.is_file():raise FileNotFoundError(args.registry)
    registry=load(args.registry)
    release=load_release_manifest(args.release_dir)
    require(set(release)==set(MODELS),'expected all16 source release')
    dest=OUT/args.run_id
    if args.resume:
        require(dest.is_dir(),f'no existing run to resume: {dest}')
        previous=load(dest/'summary.json')
        require(previous.get('format')=='t3-full-general-16_U_E_W_S_ablation-v1',
                'existing run is not a compatible sequential ablation')
        require(previous.get('requested_models')==list(selected),
                '--resume model selection differs from saved run')
        require(math.isclose(float(previous['spinor_factor_assumed']),
                             args.spinor_factor,rel_tol=0,abs_tol=1e-15),
                '--resume spinor factor differs from saved run')
        if previous.get('source_run') is not None:
            require(previous['source_run']==args.source_run,
                    '--resume source campaign differs from saved run')
    else:
        require(not dest.exists(),f'output already exists: {dest}; pass --resume to retry blocked models')
        dest.mkdir(parents=True)
    results=[]
    for model in selected:
        saved=source_root/tag(model)
        cases_path=dest/tag(model)/'cases.json'
        if args.resume and cases_path.is_file():
            record=load(cases_path)
            require(record.get('model')==model and record.get('status')=='Computed'
                    and isinstance(record.get('cases'),dict) and len(record['cases'])==16,
                    f'{model}: stored cases cannot be reused')
            results.append(record)
            print(model,'Reused 16 cases',flush=True)
            continue
        try:
            ready=preflight(model,saved,args,canonical_run_id)
            require(ready['ready'],f'{model}: missing prerequisites {ready["missing"]}; issues={ready["issues"]}')
            ref=load_source(release[model],expected_model=model)
            channels=_load_channels(model,ref,v4=args.v4_channel_dir,
                                    seeds=args.quartic_seed_dir)
            if model in ('B_-1','C_-1'):
                cases,provenance=_mixed(model,saved,(ref,channels),
                    args.reference_running_dir,registry,
                    spinor_factor=args.spinor_factor,vev=args.vev)
            else:
                cases,provenance=_nonmixed(model,saved,canonical_run_id,
                    (ref,channels),registry,
                    spinor_factor=args.spinor_factor,
                    wolframscript=args.wolframscript,vev=args.vev)
            require(len(cases)==16,f'{model}: expected sixteen U/E/W/S configurations')
            record={'model':model,'status':'Computed','cases':cases,
                    'provenance':provenance,
                    'general_Lagrangian':True,
                    'C5_hard_matching_is_independent':False,
                    'full_sequential_general_C5_certified':False}
            save(cases_path,record)
            (dest/tag(model)/'error.json').unlink(missing_ok=True)
        except Exception as e:
            record={'model':model,'status':'Blocked',
                    'error':f'{type(e).__name__}: {e}',
                    'traceback':traceback.format_exc(limit=8)}
            save(dest/tag(model)/'error.json',record)
        print(model,record['status'],record.get('error','')[:180],flush=True)
        results.append(record)
        summary={'format':'t3-full-general-16_U_E_W_S_ablation-v1',
                 'source_run':args.source_run,
                 'requested_models':list(selected),
                 'computed':sum(x['status']=='Computed' for x in results),
                 'spinor_factor_assumed':args.spinor_factor,
                 'basis_assumed':True,
                 'one_scale_C5_used_as_W_on_reference':True,
                 'hard_matching_residual_independently_calculated':False,
                 'no_double_counting_by_construction':True,
                 'model_rows':[{k:v for k,v in x.items() if k!='cases'} for x in results]}
        save(dest/'summary.json',summary)
    # Rebuild summary and CSV from every stored case, including reused models.
    summary={'format':'t3-full-general-16_U_E_W_S_ablation-v1',
             'source_run':args.source_run,
             'requested_models':list(selected),
             'computed':sum(x['status']=='Computed' for x in results),
             'spinor_factor_assumed':args.spinor_factor,
             'basis_assumed':True,
             'one_scale_C5_used_as_W_on_reference':True,
             'hard_matching_residual_independently_calculated':False,
             'no_double_counting_by_construction':True,
             'model_rows':[{k:v for k,v in x.items() if k!='cases'} for x in results]}
    save(dest/'summary.json',summary)
    write_csv(results,dest)
    print('Results:',dest/'summary.json')
    return 0 if len(results)==len(selected) and all(x['status']=='Computed' for x in results) else 2

if __name__=='__main__':
    raise SystemExit(main())
