"""All-16 full-general LLSS one-loop mixing / both-scalar threshold study.

Numerical inputs: existing source-verified F-tree LLSS release, native
Matchete quartic channels, generated general RGBeta quartic registry and the
all16 renormalisable running snapshots. No restricted C5 kernel is used.

Compute three independent EFT switches (2^3 cases):
  EFT1 scalar couplings/masses:  muF frozen or muS evolved
  LLSS -> Weinberg mixing:      leading logarithm off or on
  final SMEFT:                 SM+Weinberg frozen or evolved

The full-general UV fitted-state 2^3 ablation is a SEPARATE controller. This
runner uses the exact F-tree source from the selected release, without mixing
it with a different fitted Yukawa source. It does not pretend that reference
source is the F-boundary of an unrelated numerical fit.

This is a full-general exploratory scalar-pair LLSS contribution, not a complete
physical C5. General HHSS reference invariant directions, GammaCC -> Weyl
factor and unusual scalar basis rotations remain explicitly conditional.
"""
from __future__ import annotations

import argparse
import csv
import json
from hashlib import sha256
import math
from pathlib import Path
import re
import traceback

import numpy as np

from Numerical.running.GeneralT3SequentialLeadingLog import (
    MODELS, decode, encode, require, scalar_tensors, response,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RELEASE = (ROOT/'output/general_observable/lo_fermion_boundary/releases/'
                   'source_021dcc1e8fe8ba0f')
DEFAULT_V4 = ROOT/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'
DEFAULT_SEEDS = ROOT/'output/full16_symbolic'
DEFAULT_RUNNING = ROOT/'output/general_observable/all16_v1'
DEFAULT_OUT = ROOT/'output/general_observable/lo_eft1_mixing/general_sequential_llss'
DEFAULT_REGISTRY = ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json'


def tag(model: str) -> str:
    return model.replace('+','p').replace('-','m')


def load(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def save(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def hash_file(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def source_paths(root: dict[str, Path], model: str) -> dict[str,Path]:
    from Numerical.joint_scalar_matching.source_quartic_v4 import source_seed_path
    return {'boundary': root['release'] / tag(model) / 'boundary.json',
            'channels': root['v4'] / (tag(model)+'.json'),
            'quartic_seed': source_seed_path(root['seeds'],model),
            'running': root['running'] / tag(model) / 'running.json',
            'registry': root['registry']}


def available(path_set: dict[str,Path]) -> dict:
    missing={k:str(v) for k,v in path_set.items() if not v.is_file()}
    return {'ready':not missing,'missing':missing,
            'inputs':{k:str(v) for k,v in path_set.items()}}


def snapshot_at(model: str, original: object, raw: dict, source: str):
    """Same source/boundaries, but one chosen renormalisable scale endpoint."""
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import RunningSnapshot
    require(isinstance(raw,dict) and 'mS1Sq' in raw and 'mS2Sq' in raw,
            f'{model}: missing scalar masses in {source}')
    m1=float(raw['mS1Sq']);m2=float(raw['mS2Sq'])
    require(np.isfinite([m1,m2]).all() and min(m1,m2)>0,
            f'{model}: invalid scalar mass at {source}')
    return RunningSnapshot(model, original.source, original.sha256,
                           original.mu_f, original.mu_s, m1, m2, raw,
                           {**original.evidence,'selected_state':source})


def running_states(model: str, path: Path):
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import load_running_snapshot
    doc=load(path)
    snapshot_s=load_running_snapshot(path,model)
    frozen=doc.get('eft1_at_fermion_threshold')
    require(isinstance(frozen,dict), f'{model}: missing muF scalar EFT1 endpoint')
    snapshot_f=snapshot_at(model,snapshot_s,frozen,'muF_frozen_to_muS')
    return snapshot_f,snapshot_s


def _complex_matrix(raw):
    if isinstance(raw,dict):
        return decode(raw)
    if isinstance(raw,list) and len(raw)==3:
        arr=np.array([[complex(float(z['real']),float(z['imag']))
                       if isinstance(z,dict) else complex(z)
                       for z in row] for row in raw],dtype=complex)
        require(arr.shape==(3,3) and np.isfinite(arr).all(),'bad SMEFT Yukawa')
        return arr
    raise ValueError('Missing full three-family SM Yukawa matrix')


def _smeft_state(raw: dict, c5: np.ndarray, *, on: bool,
                 mu_s: float, mu_low: float):
    from Numerical.running.SMWeinbergEvolution import (
        SMWeinbergInitialConditions, evolve_sm_weinberg)
    initial=SMWeinbergInitialConditions(
        gY=float(raw['gY']),g2=float(raw['g2']),g3=float(raw['g3']),
        lambdaH=float(raw['lambdaH']),
        ye=_complex_matrix(raw['ye']),yu=_complex_matrix(raw['yu']),
        yd=_complex_matrix(raw['yd']),K=c5).validated()
    if on:
        end=evolve_sm_weinberg(initial,mu_s,mu_low,rtol=1e-8,atol=1e-11)
        return end.K,end.ye,int(end.nfev)
    return initial.K,initial.ye,0


def low_energy(row: dict, c5: np.ndarray, *, on: bool,
               mu_s: float, mu_low: float, vev: float, target=None) -> dict:
    from Numerical.fitting.ParameterScan import (
        neutrino_mass_from_c5, charged_lepton_mass_basis_matrix)
    from Numerical.fitting.OscillationFit import evaluate_mass_matrix_fit
    K, ye, nfev = _smeft_state(row,c5,on=on,mu_s=mu_s,mu_low=mu_low)
    m = neutrino_mass_from_c5(K,vev_gev=vev)
    charged = charged_lepton_mass_basis_matrix(m,ye)
    masses_ev = sorted((np.linalg.svd(charged,compute_uv=False)*1e9).tolist())
    out={'c5_at_low_scale':encode(K),
         'masses_ev_ascending':masses_ev,
         'smeft_nfev':nfev}
    if target is not None:
        fit=evaluate_mass_matrix_fit(charged,target)
        out.update(chi2=float(fit.chi2),ordering=str(fit.ordering),
                   prediction=dict(zip(fit.observable_names,map(float,fit.prediction))))
    return out


def process(model: str, paths: dict, registry: dict, *,
            spinor_factor: complex, include_extra: bool,
            include_portal: bool, with_smeft: bool, vev: float,
            target=None) -> dict:
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import read_v4_cached_channels
    p=source_paths(paths,model)
    missing=available(p)['missing']
    require(not missing,f'{model}: missing inputs {missing}')
    source=load_source(p['boundary'],expected_model=model)
    ch=read_v4_cached_channels(p['channels'],p['quartic_seed'],source)
    snapshot_f,snapshot_s=running_states(model,p['running'])
    require(abs(source.mu_f_gev-snapshot_s.mu_f)<1e-6,
            f'{model}: source matching scale differs from RG boundary')
    from Numerical.lo_certification.couplings.extra_hhss_geometry import model_extra_geometry
    geometry=(model_extra_geometry(model,registry['models'][model],include_components=True)
              if include_extra else None)
    cases={}
    for eft_on,snapshot in ((False,snapshot_f),(True,snapshot_s)):
        bundle=scalar_tensors(model,source,ch,snapshot,registry,
                              include_general_extra=include_extra,
                              include_portal_completion=include_portal,geometry=geometry)
        result=response(source,bundle,mu_f=snapshot.mu_f,mu_s=snapshot.mu_s,
                        spinor_factor=spinor_factor)
        for w_on in (False,True):
            c5=decode(result['c5_wilson_running_on'] if w_on else
                      result['c5_wilson_frozen_at_F'])
            for sm_on in ((False,True) if with_smeft else (False,)):
                key=f'{int(eft_on)}{int(w_on)}{int(sm_on)}'
                item={'eft1_renormalisable_running_on':eft_on,
                      'eft1_wilson_leading_log_on':w_on,
                      'smeft_running_on':sm_on,
                      'c5_at_scalar_threshold':encode(c5),
                      'c5_norm_at_scalar_threshold':float(np.linalg.norm(c5)),
                      'c5_eft1_leading_log':result['c5_eft1_leading_log'],
                      'c5_scalar_finite_bubble':result['c5_scalar_finite_bubble'],
                      'per_channel':result['per_channel'],
                      'quartic_provenance':result['quartic_provenance'],
                      'matching_scale_cancellation_residual':result['matching_scale_cancellation_residual'],
                      'weinberg_projection_residual':max(
                           result['weinberg_beta_projection_residual'],
                           result['weinberg_bubble_projection_residual'])}
                if with_smeft:
                    item.update(low_energy(snapshot.couplings,c5,on=sm_on,
                               mu_s=snapshot.mu_s,mu_low=1e2,vev=vev,target=target))
                else:
                    item['conditional_masses_ev_ascending']=(
                       sorted((np.linalg.svd(c5,compute_uv=False)*vev**2/2*1e9).tolist()))
                cases[key]=item
    return {'model':model,'status':'ExploratoryGeneralLLSSComputed',
            'source_release':str(paths['release']),
            'source_boundary_sha256':hash_file(p['boundary']),
            'running_state_sha256':hash_file(p['running']),
            'general_extra_hhss_enabled':include_extra,
            'full_rgbeta_portal_completion_enabled':include_portal,
            'spinor_convention_factor':{'real':spinor_factor.real,'imag':spinor_factor.imag},
            'cases':cases,
            'full_general_Lagrangian_renormalisable_state_used':True,
            'source_boundary_matched_to_THIS_reference_snapshot':True,
            'additional_HHSS_basis_candidate':True,
            'complete_one_loop_neutrino_mass':False}


def comparison(rows: list[dict], out: Path) -> None:
    cols=['model','status','case','eft1_on','wilson_on','smeft_on','C5_norm_GeV_inv',
          'm_light_ev','m_mid_ev','m_heavy_ev','log_to_full_ratio','error']
    with (out/'comparison.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader()
        for item in rows:
            if item['status']!='ExploratoryGeneralLLSSComputed':
                w.writerow({'model':item['model'],'status':item['status'],
                            'error':item.get('error')});continue
            cases=item['cases']
            for key,row in cases.items():
                m=row.get('masses_ev_ascending',row.get('conditional_masses_ev_ascending'))
                on=cases[f'{key[0]}1{key[2]}']['c5_norm_at_scalar_threshold']
                off=cases[f'{key[0]}0{key[2]}']['c5_norm_at_scalar_threshold']
                w.writerow({'model':item['model'],'status':item['status'],
                       'case':key,'eft1_on':key[0],'wilson_on':key[1],
                       'smeft_on':key[2],
                       'C5_norm_GeV_inv':row['c5_norm_at_scalar_threshold'],
                       'm_light_ev':m[0],'m_mid_ev':m[1],'m_heavy_ev':m[2],
                       'log_to_full_ratio':(None if on==0 else (on-off)/on)})


def parse_args(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--models',default=None,
                   help='Comma-separated model subset; omitted = all 16')
    p.add_argument('--release-dir',type=Path,default=DEFAULT_RELEASE)
    p.add_argument('--v4-channel-dir',type=Path,default=DEFAULT_V4)
    p.add_argument('--quartic-seed-dir',type=Path,default=DEFAULT_SEEDS)
    p.add_argument('--running-dir',type=Path,default=DEFAULT_RUNNING)
    p.add_argument('--registry',type=Path,default=DEFAULT_REGISTRY)
    p.add_argument('--output-dir',type=Path,default=DEFAULT_OUT/'batch_v1')
    p.add_argument('--spinor-factor',type=float,default=1.,
                   help='Explicit provisional Matchete GammaCC -> Weyl multiplier')
    p.add_argument('--skip-general-extra-hhss',action='store_true')
    p.add_argument('--skip-portal-completion',action='store_true')
    p.add_argument('--with-smeft',action='store_true',
                   help='Also evolve the generated conditional C5 in SM+Weinberg')
    p.add_argument('--vev',type=float,default=246.22)
    p.add_argument('--check-only',action='store_true')
    return p.parse_args(argv)


def main(argv=None) -> int:
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_release_manifest
    a=parse_args(argv)
    requested=tuple(x.strip() for x in a.models.split(',') if x.strip()) if a.models else MODELS
    if not requested or len(set(requested))!=len(requested) or set(requested)-set(MODELS):
        raise ValueError('Invalid all16 model selection')
    require(np.isfinite([a.spinor_factor,a.vev]).all() and a.vev>0,
            'nonfinite scalar normalization or vev')
    paths={'release':a.release_dir,'v4':a.v4_channel_dir,
           'seeds':a.quartic_seed_dir,'running':a.running_dir,
           'registry':a.registry}
    release=load_release_manifest(a.release_dir)
    require(set(release)==set(MODELS),'expected all-16 Matchete tree release')
    registry=load(a.registry)
    require('models' in registry and set(registry['models']) >= set(requested),
            'incomplete general scalar registry')
    if a.check_only:
        all_ok=True
        for model in requested:
            info=available(source_paths(paths,model))
            print(f'{model}: '+('Ready' if info['ready'] else 'Missing inputs'))
            for k,v in info['missing'].items():print('  ',k,v)
            all_ok = all_ok and info['ready']
        return 0 if all_ok else 2
    dest=a.output_dir.resolve()
    require(not dest.exists(),f'output already exists (choose new --output-dir): {dest}')
    dest.mkdir(parents=True)
    target=None
    if a.with_smeft:
        from Numerical.fitting.OscillationFit import OscillationFitTarget
        from Numerical.orchestration.RunSmallLambdaRefit import TARGET
        target=OscillationFitTarget.from_json(TARGET)
    records=[]
    for model in requested:
        try:
            result=process(model,paths,registry,spinor_factor=complex(a.spinor_factor),
                           include_extra=not a.skip_general_extra_hhss,
                           include_portal=not a.skip_portal_completion,
                           with_smeft=a.with_smeft,vev=a.vev,target=target)
            save(dest/tag(model)/'cases.json',result)
            print(model,result['status'])
        except Exception as e:
            result={'model':model,'status':'Blocked',
                    'error':f'{type(e).__name__}: {e}',
                    'traceback':traceback.format_exc(limit=7)}
            save(dest/tag(model)/'error.json',result)
            print(model,'Blocked:',result['error'])
        records.append(result)
    comparison(records,dest)
    summary={'format':'t3-general-16-conditional-sequential-leading-log-campaign-v1',
             'model_count':len(records),
             'computed_count':sum(x['status']=='ExploratoryGeneralLLSSComputed' for x in records),
             'source_release_root':str(paths['release']),
             'requested_models':list(requested),
             'settings':{'spinor_factor':a.spinor_factor,
                         'with_smeft':a.with_smeft,
                         'general_extra_hhss':not a.skip_general_extra_hhss,
                         'full_portal_completion':not a.skip_portal_completion},
             'note':('This campaign switches EFT1 scalar parameter running and the LLSS '
                      'leading C5 log independently. The F-tree source is tied to the '
                      'selected all16 running reference and MUST NOT be interpreted as '
                      'a different fitted Yukawa source. It is a conditional LLSS sector, '
                      'not a complete hard-matched physical C5.'),
             'models':[{k:v for k,v in x.items() if k!='cases'} for x in records]}
    save(dest/'summary.json',summary)
    print('Summary:',dest/'summary.json')
    return 0 if summary['computed_count']==len(requested) else 2

if __name__=='__main__':
    raise SystemExit(main())
