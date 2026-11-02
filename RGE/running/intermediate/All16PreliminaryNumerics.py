"""All-16 source-backed *conditional* T3 lambdaT3-only numerical comparison.

Each selected source is the saved F-tree LLSS boundary at muF; the two scalar
multiplets are integrated out jointly at muS. The scalar input is either an
independent illustrative running record (8 original-basis models), or an
explicit F/muS mixed-basis handoff (B_-1/C_-1). Six split-real cases are NOT
silently identified with the original spherical scalar field basis.

This work does not establish a physical C5: N_GammaCC and the absolute
Matchete <-> RGBeta lambdaT3 map remain uncertified; other quartic channels,
hard F loop and all other EFT1 operators are not included.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from hashlib import sha256
from math import isfinite, log, pi
from pathlib import Path
from typing import Any

import numpy as np

from RGE.running.intermediate.B1PreliminaryNumerics import (
    PilotError, cjson, complex_array, digest, mass_original_from_eigenbasis,
    q_in_fermion_basis, read, real_lift_b1, require, validate_basis_record,
)

FORMAT = 't3-all16-lambdaT3-one-channel-conditional-batch-v1'
LOOP = 1.0 / (16.0 * pi * pi)
SPLIT_MODELS = frozenset(('A_-2', 'A_+0', 'D_-2', 'D_+0', 'E_-2', 'E_+0'))
MIXED_MODELS = frozenset(('B_-1', 'C_-1'))
DATA_ROOT = 'output/general_observable'


def tag(model: str) -> str:
    return model.replace('-', 'm').replace('+', 'p')


class MissingCertifiedBasis(PilotError):
    """The scalar source basis cannot be identified with the LLSS basis."""


@dataclass(frozen=True)
class ScalarPoint:
    model: str
    source_to_original_R: np.ndarray
    mass2_boundary_muS: np.ndarray
    muS: float
    lambda_candidate: complex
    inputs: dict[str, Any]
    status: str


def _real_from_number(x: Any, label: str) -> float:
    require(type(x) in (int, float) and isfinite(x), f'{label} is not a finite real input')
    return float(x)


def _lambda(x: Any) -> complex:
    if isinstance(x, dict):
        require(set(x) >= {'real', 'imag'}, 'coupling requires real/imag')
        z = complex(_real_from_number(x['real'], 'lambdaT3.real'),
                    _real_from_number(x['imag'], 'lambdaT3.imag'))
    else:
        z = complex(_real_from_number(x, 'lambdaT3'))
    require(np.isfinite(z) and abs(z)>0, 'lambdaT3 candidate zero/nonfinite')
    return z


def _ordinary_scalar_point(root: Path, model: str, source, ch) -> ScalarPoint:
    """Explicit *illustrative cross-workstream* muS point, not a certified fit."""
    path = root / DATA_ROOT / 'all16_v1' / tag(model) / 'running.json'
    j = read(path)
    require(j.get('model') == model and
            j.get('status') == 'NumericalRenormalisableRunningOnly' and
            j.get('illustrative_benchmark') is True and
            j.get('general_matching_performed') is False and
            j.get('physical_mass_parameters_available') is False and
            j.get('eft1_wilson_running_performed') is False,
            f'{model}: illustrative running provenance unexpectedly changed')
    scales=j.get('scales_gev',{})
    require(float(scales['F'])==source.mu_f_gev and
            float(scales['S'])>0 and float(scales['S'])<source.mu_f_gev,
            f'{model}: inconsistent EFT1 F/S scale')
    s=j['eft1_at_scalar_threshold']
    require(tuple(source.labels)==ch.source_basis_labels,
            f'{model}: original and LLSS scalar layouts DIFFER; cannot use identity')
    d1, d2={'A':(1,3),'B':(2,2),'C':(2,2),'D':(3,1),'E':(3,3)}[model[0]]
    masses=[_real_from_number(s['mS1Sq'], 'mS1Sq')]*(2*d1)
    masses+= [_real_from_number(s['mS2Sq'], 'mS2Sq')]*(2*d2)
    require(len(masses)==len(source.labels) and min(masses)>0,
            f'{model}: invalid muS scalar masses')
    U=np.eye(len(masses))
    return ScalarPoint(model,U,np.diag(masses),float(scales['S']),_lambda(s['lambdaT3']),
        {'muS_running_json_path':str(path.relative_to(root)).replace('\\','/'),
         'muS_running_json_sha256':digest(path),
         'muS_scalar_quadratic_provenance':'illustrative general renormalisable running only',
         'muF_matched_to_this_running_benchmark':False,
         'muS_scalar_eigenbasis':'original S1/S2, common-mass unbroken multiplets',
         'native_RGBeta_lambdaT3_equality_certified':False},
        'IllustrativeCrossWorkstreamLambdaT3Point_NotSyncedPhysicalBenchmark')


def _mixed_scalar_point(root: Path, model: str, source, ch) -> ScalarPoint:
    slug=tag(model)
    fpath=root/DATA_ROOT/'mixed_f_synchronised_v1'/slug/'mass_basis'/f'{slug}_final_mass_basis.json'
    spart=root/'output/general_t3_mass_basis'/f'{slug}_final_mass_basis.json'
    ordered=root/DATA_ROOT/'ordered_wilson_tensor_v4'/slug/'ordered_tensor.json'
    fparams=root/DATA_ROOT/'mixed_f_synchronised_v1'/slug/'matchete'/f'{slug}_final_parameters.json'
    jf,js,jo,jp=(read(p) for p in (fpath,spart,ordered,fparams))
    require(jo.get('format')=='t3-mixed-native-matchete-ordered-tensor-v1' and
            jo.get('model')==model and float(jo.get('mu_F_gev',-1))==source.mu_f_gev,
            f'{model}: ordered F boundary provenance changed')
    require(jp.get('format')=='t3-experimental-mass-eigen-matchete-parameters-v1' and
            jp.get('model')==model and
            float(jp.get('source_eft1_scale_gev',-1))==source.mu_f_gev and
            jp.get('source_sha256')==jo.get('source_mass_basis_sha256') and
            str(jp.get('source_path', '')).replace('\\','/')==str(fpath.relative_to(root)).replace('\\','/'),
            f'{model}: Matchete native F handoff does not identify the F scalar basis')
    require(jf.get('status')=='ExperimentalMassBasisTensorTransportNoPhysicalMatching' and
            js.get('status')=='ExperimentalMassBasisTensorTransportNoPhysicalMatching',
            f'{model}: unknown F/muS mass-record status')
    Uf,_=validate_basis_record(jf,model,source.mu_f_gev)
    require(np.allclose(Uf,complex_array(jp['unitary_U'],(2,2)),atol=1e-13,rtol=0),
            f'{model}: native Matchete F-unitary mismatch')
    mus=float(js['scale_gev'])
    require(0<mus<source.mu_f_gev, f'{model}: invalid muS')
    Us,eval_s=validate_basis_record(js,model,mus)
    require(tuple(source.labels)==('NewScalarX',)*4+('NewScalarY',)*4 and ch.n==8,
            f'{model}: mixed source basis changed')
    Rf=real_lift_b1(Uf)      # source original real = Rf @ F X/Y real
    mass_orig=mass_original_from_eigenbasis(Us,eval_s)
    mass_F=Rf.T @ mass_orig @ Rf
    require(np.allclose(np.sort(np.linalg.eigvalsh(mass_F)),
                        np.sort(np.repeat(eval_s,4)),rtol=2e-11,atol=1e-7),
            f'{model}: F-to-muS scalar masses inconsistent')
    return ScalarPoint(model,Rf,mass_F,mus,
        _lambda(js['quartics']['lambdaT3']['source_rgbeta_coupling']),
        {'muF_mass_basis_json_path':str(fpath.relative_to(root)).replace('\\','/'),
         'muF_mass_basis_sha256':digest(fpath),
         'muF_basis_source_sha256':jp['source_sha256'],
         'muF_native_parameters_sha256':digest(fparams),
         'muS_mass_basis_json_path':str(spart.relative_to(root)).replace('\\','/'),
         'muS_mass_basis_sha256':digest(spart),
         'muS_scalar_quadratic_provenance':'numeric mixed-doublet mass basis at muS',
         'muF_matched_to_this_running_benchmark':True,
         'native_RGBeta_lambdaT3_equality_certified':False},
        'SourceSyncedMixedDoubletCandidate_NotPhysical')


def resolve_scalar_point(root: Path, model: str, source, ch) -> ScalarPoint:
    if model in SPLIT_MODELS:
        raise MissingCertifiedBasis(f'{model}: native Matchete S1/S2 quartic is in an original complex-scalar '
            'basis while fermion LLSS is in a real self-bilinear minus/plus basis. '
            'No independently source-matched real intertwiner/phase at muF and muS; '
            'no identity or zero-splitting substitute permitted.')
    if model in MIXED_MODELS:
        return _mixed_scalar_point(root,model,source,ch)
    return _ordinary_scalar_point(root,model,source,ch)


def _mass_c5(C: np.ndarray, quartic: np.ndarray, mass: np.ndarray,
             muF: float, muS: float, *, vev: float=246.22) -> dict:
    """One channel fixed-order O(hbar): LLSS log + all-scalar bubble, no C5 tree."""
    from Numerical.joint_scalar_matching.core import (
        scalar_bubble_matching, project_weinberg, check_tensors,
    )
    n=C.shape[-1]
    C,q=check_tensors(C,quartic,n)
    require(q.shape==(4,4,n,n) and mass.shape==(n,n) and
            0<muS<muF, 'incompatible source/scalar tensor or scales')
    beta=np.einsum('fgijab,ABab->fgijAB',C,q,optimize=True)
    bp,br=project_weinberg(beta,strict=False)
    amp=scalar_bubble_matching(C,q,mass,muS,reference_scale=muS)
    sp,sr=project_weinberg(amp['amplitude'],strict=False)
    ll=bp*log(muS/muF)*LOOP
    total=ll+sp
    other_mu=muS*np.exp(0.03125)
    other=scalar_bubble_matching(C,q,mass,other_mu,reference_scale=muS)
    op,_=project_weinberg(other['amplitude'],strict=False)
    shifted=bp*log(other_mu/muF)*LOOP+op
    norm=max(float(np.linalg.norm(total)),float(np.linalg.norm(ll)),
             float(np.linalg.norm(sp)),1e-30)
    cancel=float(np.linalg.norm(shifted-total)/norm)
    require(cancel<=2e-9, f'fixed-coupling log/bubble scale cancellation fails: {cancel:g}')
    invariant=bool(br<=1e-8 and sr<=1e-8)
    require(float(np.linalg.norm(total-total.T)) <= 1e-8*norm,
            'flavour asymmetry in projected C5')
    out={'status':'ConditionalWeinbergLambdaT3OneChannel' if invariant
                  else 'NonWeinbergTensorBlockedMassInterpretation',
         'full_SU2_Weinberg_projection_pass':invariant,
         'beta_SU2_residual_relative':float(br),
         'bubble_SU2_residual_relative':float(sr),
         'matching_scale_cancellation_relative':cancel,
         'beta_C5_master_16pi2_per_N_gev_inverse':cjson(2*bp),
         'C5_master_leading_log_per_N_gev_inverse':cjson(2*ll),
         'C5_master_joint_scalar_bubble_per_N_gev_inverse':cjson(2*sp),
         'C5_master_total_per_N_gev_inverse':cjson(2*total),
         'joint_to_master_C5_coefficient_factor':2,
         'scalar_mass_eigenvalues_gev':np.sqrt(np.linalg.eigvalsh(mass)).tolist(),
         'physical_C5_certified':False,
         'neutrino_mass_singular_values_ev_per_abs_NGamma':None}
    if invariant:
        sv=np.linalg.svd(.5*vev*vev*(2*total)*1e9,compute_uv=False)
        out['neutrino_mass_singular_values_ev_per_abs_NGamma']=sorted(float(x) for x in sv)
    return out


def evaluate_model(root: Path, model: str) -> dict:
    """No writes. Read and validate the actual 16-model source release."""
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.source_quartic_v4 import (
        build_source_channels, source_seed_path, evaluate_channels,
    )
    root=Path(root)
    sr=root/DATA_ROOT/'lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'/tag(model)/'boundary.json'
    native=root/DATA_ROOT/'lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'/f'{tag(model)}.json'
    source=load_source(sr,expected_model=model)
    require(source.mu_f_gev==100000. and source.normalization_certified is False,
            'unexpected/physical-normalization source release; inspect before proceeding')
    seed=source_seed_path(root/'output/full16_symbolic',model)
    ch=build_source_channels(seed,model)
    n=2*({'A':4,'B':4,'C':4,'D':4,'E':6}[model[0]])
    require(n==len(source.labels)==ch.n, 'heavy scalar dimension disagrees')
    native_j=read(native)
    require(native_j.get('format')=='t3-source-hhss-quartic-channel-export-v4' and
            native_j.get('model')==model and
            native_j.get('seed_quartic_sha256')==ch.seed_sha256 and
            native_j.get('source_boundary_sha256')==source.hash and
            native_j.get('couplings_required_at_muS')==sorted(ch.by_name) and
            native_j.get('source_basis_matches_fermion_boundary')==(tuple(source.labels)==ch.source_basis_labels),
            'native quartic v4 export not tied to source/seed or layout')
    point=resolve_scalar_point(root,model,source,ch)
    names=tuple(sorted(ch.by_name))
    require('lambdaT3' in names, f'{model}: no native lambdaT3 channel')
    # Explicit one-channel selection; no claim that other quartics physically vanish.
    supplied={name:(cjson(point.lambda_candidate) if name=='lambdaT3' else 0.) for name in names}
    qchannels,parse_names=evaluate_channels(ch,supplied)
    require(parse_names==names and np.isfinite(qchannels).all(), 'native quartic channel assembly failure')
    R=point.source_to_original_R
    require(R.shape==(n,n) and np.allclose(R.T@R,np.eye(n),atol=2e-11,rtol=2e-11),
            f'{model}: invalid original-to-source scalar transformation')
    q=np.einsum('kABab,ai,bj->kABij',qchannels,R,R,optimize=True)[names.index('lambdaT3')]
    result=_mass_c5(source.raw,q,point.mass2_boundary_muS,source.mu_f_gev,point.muS)
    # Independent total == sum over every native PL source, including suppressed sources.
    from Numerical.joint_scalar_matching.core import project_weinberg
    cumulative=np.zeros((3,3),complex)
    per_source=[]
    for i,ident in enumerate(source.source_indices):
        beta=np.einsum('fgijab,ABab->fgijAB',source.by_source[i],q,optimize=True)
        projected,_=project_weinberg(beta,strict=False)
        cumulative+=2*projected
        per_source.append({'PL_source_index':int(ident),
                           'beta_C5_master_per_N_gev_inverse':cjson(2*projected)})
    got=np.asarray(result['beta_C5_master_16pi2_per_N_gev_inverse']['real'])+1j*np.asarray(result['beta_C5_master_16pi2_per_N_gev_inverse']['imag'])
    require(float(np.linalg.norm(got-cumulative))<=2e-9*max(float(np.linalg.norm(got)),1e-30),
            f'{model}: PL-source beta sum mismatch')
    return {'format':FORMAT,'model':model,'status':result['status'],
            'source_point_class':point.status,
            'fermion_boundary_sha256':source.hash,
            'native_quartic_seed_sha256':ch.seed_sha256,
            'native_quartic_export_sha256':digest(native),
            'muF_gev':source.mu_f_gev,'muS_gev':point.muS,
            'source_count':len(source.source_indices),
            'n_real_heavy_scalars':n,
            'all_native_quartic_channels':list(names),
            'included_native_quartic_channel':'lambdaT3',
            'other_quartic_channels_excluded_not_zero':[x for x in names if x!='lambdaT3'],
            'lambdaT3_RGBeta_as_Matchete_hypothesis':cjson(point.lambda_candidate),
            'cross_workstream_input_provenance':point.inputs,
            'per_PL_source_beta':per_source,
            'response':result,
            'unknown_GammaCC_factor':None,
            'Matchete_RGBeta_coupling_normalization_certified':False,
            'physical_C5_computed':False,
            'full_eft1_seven_term_running_included':False,
            'hard_F_one_loop_matching_included':False,
            'both_scalars_integrated_together_at_muS':True,
            'interpretation':'Conditional one-channel lambdaT3 source response only, per unit unknown GammaCC; '
              'not a full EFT1 RGE/one-loop UV matching/physical mass. Source geometry and numeric '
              'coupling may originate in different illustrative running records.'}


def summary_row(doc: dict) -> dict:
    m=doc['model']
    row={'model':m,'status':doc['status'],
         'input_class':doc.get('source_point_class'),
         'source_count':doc.get('source_count'),
         'SU2_Ward_pass':doc['response']['full_SU2_Weinberg_projection_pass'],
         'SU2_beta_residual':doc['response']['beta_SU2_residual_relative'],
         'matching_scale_cancellation':doc['response']['matching_scale_cancellation_relative'],
         'physical_mass_certified':False}
    vals=doc['response']['neutrino_mass_singular_values_ev_per_abs_NGamma']
    row['masses_ev_per_abs_NGamma']=vals
    return row


def output_summary(rows: list[dict]) -> dict:
    evaluated=sum(r.get('status')=='ConditionalWeinbergLambdaT3OneChannel' for r in rows)
    nonward=sum(r.get('status')=='NonWeinbergTensorBlockedMassInterpretation' for r in rows)
    blocked=sum(r.get('status')=='BlockedMissingCertifiedScalarBasis' for r in rows)
    failures=sum(r.get('status')=='InputError' for r in rows)
    return {'format':FORMAT,'status':'ConditionalPartialBatch_NotPhysical' if blocked or nonward or failures
                     else 'AllConditionalOneChannel_NotPhysical',
            'models_selected':len(rows),'models_with_Ward_valid_one_channel_results':evaluated,
            'models_blocked_special_scalar_basis':blocked,
            'models_non_Weinberg':nonward,
            'models_input_errors':failures,
            'physical_neutrino_masses_computed':False,
            'source_to_Weyl_factor':None,
            'input_class_warning':('8 ordinary models use illustrative EFT1 scalar trajectories not '
              'certified as synchronized to the F-tree source benchmark; B/C mixed models use '
              'explicit synchronized F scalar bases and separately recorded muS mass bases.'),
            'other_quartic_channels_included':False,
            'EFT_sequence':'UV -> ALL F @muF -> EFT1(S1,S2) -> ALL scalars @one muS -> SMEFT',
            'models':rows}
