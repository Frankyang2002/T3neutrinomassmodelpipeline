"""Batch 6: native Matchete HHSS EWSB Hessian, candidate maps, and UV/EFT diagnostics.

Read-only input: completed Batch4 & Batch5 A_-4 plus native F-only quartic seed.
Outputs only a new lo_certification_benchmark/batch6_portal_* directory.
No ΓCC, physical matching, or Matchete-RGBeta alias certification is assumed.
"""
from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any

import numpy as np

from Numerical.lo_certification.benchmark import common_input_batch4 as b4
from Numerical.lo_certification.benchmark import common_point_batch5 as b5

ROOT = Path(__file__).resolve().parents[3]
OWNED = Path('output/general_observable/lo_certification_benchmark')
FORMAT = 't3-a-m4-native-hhss-hessian-two-alias-branches-batch6-v1'
BRANCHES = {'adjoint_alias': 'lambdaH2Adj', 'ordinary_alias': 'lambdaH2'}
NATIVE = ('lambdaH1Inv1', 'lambdaH2Inv2', 'lambdaT3')
STAGES = ('UV', 'F', 'S')


class Batch6Error(ValueError):
    pass


def require(condition: Any, reason: str) -> None:
    if not bool(condition):
        raise Batch6Error(reason)


def read(path: Path) -> dict:
    require(path.is_file(), 'Missing input ' + str(path))
    d = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(d, dict), 'Not a JSON object: ' + str(path))
    return d


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def encode(z: np.ndarray) -> dict:
    v = np.asarray(z, complex)
    require(bool(np.isfinite(v).all()), 'Nonfinite complex tensor')
    return {'real': v.real.tolist(), 'imag': v.imag.tolist()}


def value(raw: Any, name: str) -> complex:
    if isinstance(raw, dict) and set(raw) == {'real', 'imag'}:
        z = complex(float(raw['real']), float(raw['imag']))
    elif type(raw) in (float, int):
        z = complex(raw)
    else:
        raise Batch6Error('Unrecognised scalar encoding: ' + name)
    require(math.isfinite(z.real) and math.isfinite(z.imag), 'Nonfinite ' + name)
    return z


def scalar_stage(frozen: dict, running: dict, stage: str) -> dict:
    require(stage in STAGES, 'Unknown scalar stage')
    data = {
        'UV': frozen['base_state']['ordinary'],
        'F': running['uv_at_fermion_threshold'],
        'S': running['eft1_at_scalar_threshold'],
    }[stage]
    names = ('mS1Sq', 'mS2Sq', 'lambdaH1', 'lambdaH2',
             'lambdaH2Adj', 'lambdaS2Adj', 'lambdaT3')
    require(set(names) <= set(data), stage + ': incomplete scalar running state')
    for k in names:
        value(data[k], k)
    require(value(data['mS1Sq'], 'mS1Sq').imag == 0 and
            value(data['mS2Sq'], 'mS2Sq').imag == 0, 'Complex scalar mass squared')
    return data


def candidate_couplings(state: dict, branch: str) -> dict:
    require(branch in BRANCHES, 'Unknown candidate mapping')
    return {'lambdaH1Inv1': state['lambdaH1'],
            'lambdaH2Inv2': state[BRANCHES[branch]],
            'lambdaT3': state['lambdaT3']}


def mass_hessian_from_quartic(q: np.ndarray, m1sq: float, m2sq: float,
                              vev: float) -> tuple[np.ndarray, np.ndarray]:
    """Native V4=lambda_{ABab} h_A h_B s_a s_b/4!, H^0=(v+h+iG)/sqrt2.

    Interleaved scalar real axes: (S1_R,S1_I,S2_1R,S2_1I,...).
    Higgs neutral real axis is H[2] (zero-based; H=(H+[0],H0[1])).
    Hessian_ij = mass^2_ij + 1/2 v^2 lambda_{H0R H0R ij}.
    """
    t = np.asarray(q)
    require(t.shape == (4, 4, 8, 8) and np.isrealobj(t) and
            bool(np.isfinite(t).all()), 'Native HHSS tensor must be real (4,4,8,8)')
    require(np.max(np.abs(t-t.swapaxes(0, 1))) <= 1e-11 * max(np.max(abs(t)), 1e-30),
            'Higgs leg symmetry violated')
    require(np.max(np.abs(t-t.swapaxes(2, 3))) <= 1e-11 * max(np.max(abs(t)), 1e-30),
            'Heavy scalar leg symmetry violated')
    require(math.isfinite(vev) and vev > 0 and math.isfinite(m1sq) and
            math.isfinite(m2sq) and min(m1sq, m2sq) > 0, 'Invalid scalar mass/VEV')
    bare = np.diag([float(m1sq)]*2 + [float(m2sq)]*6)
    delta = 0.5 * vev**2 * t[2, 2]
    h = bare + delta
    require(np.linalg.norm(h-h.T) < 1e-11 * max(np.linalg.norm(h), 1.),
            'EWSB scalar Hessian not symmetric')
    require(np.linalg.eigvalsh(h).min() > 0, 'Nonpositive EWSB scalar mass squared')
    return (h+h.T)/2, delta


def hessian_to_complex(h_interleaved: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Invert cartesian_hessian(A,B), preserving full complex A†=A and B^T=B."""
    h = np.asarray(h_interleaved, float)
    require(h.shape == (8,8) and np.isfinite(h).all() and
            np.max(abs(h-h.T)) < 1e-9 * max(np.max(abs(h)),1.),
            'Invalid interleaved real scalar Hessian')
    r = np.arange(0,8,2); i = np.arange(1,8,2)
    rr = h[np.ix_(r,r)]; ii = h[np.ix_(i,i)]
    ri = h[np.ix_(r,i)]; ir = h[np.ix_(i,r)]
    a = (rr+ii)/2 + 0.5j*(ir-ri)
    b = (rr-ii)/2 - 0.5j*(ri+ir)
    require(np.linalg.norm(a-a.conj().T) < 1e-11*max(np.linalg.norm(a),1.),
            'Complex Hermitian A reconstruction failed')
    require(np.linalg.norm(b-b.T) < 1e-11*max(np.linalg.norm(b),1.),
            'Complex symmetric B reconstruction failed')
    # Verify by the same Cartesian-Hessian identity as the existing UV source.
    back = np.block([[a.real+b.real, -a.imag-b.imag],
                     [a.imag-b.imag, a.real-b.real]])
    perm = np.ravel(np.column_stack((np.arange(4),np.arange(4,8))))
    back_interleaved = back[np.ix_(perm,perm)]
    require(np.linalg.norm(h-back_interleaved) < 2e-11*max(np.linalg.norm(h),1.),
            'Complex A/B conversion does not reproduce native real Hessian')
    return a,b


def classify_omissions(omitted: dict) -> dict:
    require(isinstance(omitted,dict), 'Missing Batch5 projection record')
    extra = set(omitted) - {'lambdaH2Adj','lambdaS2Adj'}
    require(not extra, 'Unresolved additional projected interactions: ' + str(sorted(extra)))
    return {
        'candidate_higgs_scalar_portal': omitted.get('lambdaH2Adj'),
        'pure_scalar_self_quartic': omitted.get('lambdaS2Adj'),
        'pure_scalar_zero_background_tree_hessian': True,
        'note': 'lambdaS2Adj has four heavy scalar legs and no H; zero scalar VEV implies no tree EWSB quadratic term. Its RGE effects remain in evolved mass parameters.',
    }


def validate_old(run4: dict, run5: dict) -> None:
    require(run4.get('format') == b4.FORMAT and
            run4.get('status') == 'COMMON_FITTED_INPUT_CONDITIONAL_LLSS_SCALAR_RESPONSE',
            'Unfinished fitted Batch4')
    require(run5.get('format') == b5.FORMAT and
            run5.get('status') == 'CONDITIONAL_NUMERIC_COMPARISON_NOT_SAME_PHYSICAL_MATCHING',
            'Unfinished projected Batch5')
    require(run4.get('physical_certification',{}).get('physical_C5_certified') is False and
            run5.get('interpretation',{}).get('GammaCC_signed_factor') is None,
            'Physical certification unexpectedly changed')
    require(set(run5.get('UV_scenarios',{})) == set(STAGES),
            'Batch5 missing projected stages')


def preflight(root: Path, batch4_id: str, batch5_id: str) -> dict:
    root = root.resolve()
    for name in (batch4_id, batch5_id):
        require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,50}', name) is not None,
                'Unsafe source run ID')
    info = b5.preflight(root, batch4_id)
    f4 = root/OWNED/('batch4_common_input_'+batch4_id)
    f5 = root/OWNED/('batch5_uv_eft_'+batch5_id)
    a4 = read(f4/'summary.json'); a5 = read(f5/'summary.json')
    validate_old(a4,a5)
    require(a5['source_shas']['batch4_running'] == sha(f4/'fitted_running.json'),
            'Batch5 was evaluated from different Batch4 running input')
    require(a5['source_shas']['batch4_summary'] == sha(f4/'summary.json'),
            'Batch5 was evaluated from different Batch4 summary')
    originals = {}
    for stage in STAGES:
        p=f5/(stage+'_direct_input.json')
        obj=read(p)
        require(obj.get('model')=='A_-4' and
                obj.get('provenance',{}).get('scalar_stage')==stage,
                'Batch5 UV scenario/input mismatch: '+stage)
        oldom=a5['UV_scenarios'][stage]['unsupported_evolved_scalar_couplings_dropped']
        require(obj['provenance']['unsupported_evolved_scalar_couplings_projected_to_zero']==oldom,
                'Batch5 projection record drift: '+stage)
        classify_omissions(oldom)
        originals[stage]=p
    frozen,tree,seed,old=b4.check_source(root)
    names = sorted({re.search(r'Coupling\[([^,\]]+),',t['TermInputForm']).group(1)
                    for t in seed['ScalarQuarticTerms']
                    if t['HiggsFieldCount']==2 and t['Scalar1FieldCount']+t['Scalar2FieldCount']==2})
    require(tuple(names)==NATIVE, 'Native A_-4 HHSS coupling inventory changed: '+str(names))
    for stage in STAGES:
        x=scalar_stage(frozen,read(f4/'fitted_running.json'),stage)
        require(value(x['mS1Sq'],'m1').real > 0 and value(x['mS2Sq'],'m2').real>0,
                'Nonpositive scalar masses at '+stage)
    return {'model':'A_-4','status':'PREFLIGHT_OK_NO_OUTPUT_WRITTEN',
            'native_hhss_channels':names,'portal_alias_hypotheses':BRANCHES,
            'pure_scalar_lambdaS2Adj_not_an_HHSS_portal':True,
            'batch4':str(f4),'batch5':str(f5),
            'physical_matching_certified':False,'output':None}


def calculate(root: Path, batch4_id: str, batch5_id: str, run_id: str) -> dict:
    root = root.resolve()
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,50}', run_id) is not None,
            'Unsafe new run ID')
    target = root/OWNED/('batch6_portal_'+run_id)
    require(not target.exists(), 'Output directory exists: '+str(target))
    preflight(root,batch4_id,batch5_id)
    from Numerical.joint_scalar_matching.source_quartic_v4 import (
        build_source_channels,evaluate_channels)
    from Numerical.lo_certification.direct_uv.native_uv_loop_reference import (
        independent_reference,decode,cartesian_hessian)
    from Numerical.direct_one_loop.core import evaluate

    f4 = root/OWNED/('batch4_common_input_'+batch4_id)
    f5 = root/OWNED/('batch5_uv_eft_'+batch5_id)
    running=read(f4/'fitted_running.json'); fitted=read(f4/'summary.json')
    previous=read(f5/'summary.json')
    frozen=read(root/b4.BENCH)
    seedpath=root/b4.QUARTIC_SEED
    ch=build_source_channels(seedpath,'A_-4')
    require(tuple(sorted(ch.by_name))==NATIVE and ch.n==8,
            'Source HHSS tensor changed or is incomplete')
    eft=b5.z3(fitted['conditional_mnu_eV_per_N'],'Fitted conditional EFT mass')
    require(np.linalg.norm(eft)>0,'EFT comparison is zero')
    states={stage:scalar_stage(frozen,running,stage) for stage in STAGES}
    complete={}; records={}
    vev=float(frozen['numerical']['vev_gev'])
    require(240 < vev < 255, 'Unexpected benchmark Higgs VEV')
    for stage in STAGES:
        legacy=read(f5/(stage+'_direct_input.json'))
        olduv=read(f5/(stage+'_direct_independent.json'))
        oldmass=decode(olduv['complex_symmetric_matrix_gev'],'Batch5 UV mass')
        m1=value(states[stage]['mS1Sq'],'mS1Sq');m2=value(states[stage]['mS2Sq'],'mS2Sq')
        require(m1.imag==m2.imag==0,'Complex scalar mass')
        per={}
        for branch in BRANCHES:
            couplings=candidate_couplings(states[stage],branch)
            channels,names=evaluate_channels(ch,couplings)
            require(set(names)==set(NATIVE) and channels.shape==(3,4,4,8,8),
                    'Native source HHSS numerical component inventory drift')
            summed=np.sum(channels,axis=0)
            h,delta=mass_hessian_from_quartic(summed,float(m1.real),float(m2.real),vev)
            a,b=hessian_to_complex(h)
            payload=deepcopy(legacy)
            payload['scalar_A_gev2']=encode(a)
            payload['scalar_B_gev2']=encode(b)
            payload['provenance']['batch6_native_scalar_potential_hessian']=True
            payload['provenance']['batch6_native_source_quartic_sha256']=sha(seedpath)
            payload['provenance']['batch6_hhss_coupling_alias_hypothesis']=dict(
                lambdaH1Inv1='lambdaH1',lambdaT3='lambdaT3',
                lambdaH2Inv2=BRANCHES[branch])
            payload['provenance']['batch6_prior_legacy_scalar_A_B_replaced']=True
            payload['provenance']['batch6_quartic_alias_and_spinor_physically_certified']=False
            payload['provenance']['batch6_pure_scalar_lambdaS2Adj_at_zero_scalar_vev_not_in_EWSB_Hessian']=True
            uv=independent_reference(payload,check_rephase=False)
            calc=evaluate(payload,retain_individual_terms=False)
            nu=decode(uv['complex_symmetric_matrix_gev'],'native portal UV')
            nu2=decode(calc['total']['complex_symmetric_matrix_gev'],'second UV')
            relative=float(np.linalg.norm(nu-nu2)/max(np.linalg.norm(nu),np.linalg.norm(nu2),1e-100))
            require(relative<4e-7,'Independent/production UV numerical discrepancy')
            require(calc['diagnostics']['uv_finite'] is True and
                    uv['diagnostics']['uv_pole_relative_to_natural_amplitude']<1e-9,
                    'One-loop UV pole did not cancel')
            require(np.min(uv['fermion_mass_singular_values_gev'])>
                    np.max(uv['scalar_masses_gev']), 'Heavy-fermion hierarchy broken')
            oldh=cartesian_hessian(decode(legacy['scalar_A_gev2'],'old A'),
                                   decode(legacy['scalar_B_gev2'],'old B'))
            # Grouped (R_1,..,R_4,I_1,..,I_4) -> interleaved (R1,I1,...)
            perm=np.ravel(np.column_stack((np.arange(4),np.arange(4,8))))
            oldh=oldh[np.ix_(perm,perm)]
            mass_difference=float(np.linalg.norm(h-oldh)/max(np.linalg.norm(oldh),1e-100))
            result=b5.comparison(nu,eft)
            result.update({
                'scalar_stage':stage, 'candidate_mapping':branch,
                'alias_lambdaH2Inv2':BRANCHES[branch],
                'nu_mass_singular_values_eV_descending':np.linalg.svd(nu*1e9,compute_uv=False).tolist(),
                'scalar_eigenmasses_gev':uv['scalar_masses_gev'],
                'source_native_quartic_names':list(names),
                'quartic_channel_values':{k:{'real':value(v,k).real,'imag':value(v,k).imag}
                                          for k,v in couplings.items()},
                'native_EWSB_mass_hessian_shift_frobenius_GeV2':float(np.linalg.norm(delta)),
                'holomorphic_scalar_B_frobenius_GeV2':float(np.linalg.norm(b)),
                'native_vs_Batch5_legacy_EWSB_mass_hessian_relative':mass_difference,
                'complex_UV_mass_vs_projected_Batch5_relative':float(
                    np.linalg.norm(nu-oldmass)/max(np.linalg.norm(oldmass),1e-100)),
                'recomputed_UV_pole_relative':uv['diagnostics']['uv_pole_relative_to_natural_amplitude'],
                'independent_vs_production_UV_relative':relative,
                'pure_scalar_quartic_omission_classification':classify_omissions(
                    previous['UV_scenarios'][stage]['unsupported_evolved_scalar_couplings_dropped']),
                'physical_C5_certified':False,
            })
            per[branch]=result
            records[(stage,branch)]={'input':payload,'independent':uv,'production':calc,
                'Hessian_real_interleaved_GeV2':h.tolist(),
                'delta_Hessian_from_HHSS_GeV2':delta.tolist()}
        complete[stage]=per
    summary={
        'format':FORMAT,'model':'A_-4',
        'status':'NATIVE_HHSS_HESSIAN_CANDIDATE_ALIASES_UNCERTIFIED',
        'source_benchmark_SHA256':b4.BENCH_SHA,
        'source_quartic_seed_SHA256':sha(seedpath),
        'batch4_summary_SHA256':sha(f4/'summary.json'),
        'batch5_summary_SHA256':sha(f5/'summary.json'),
        'native_potential_convention':'V4=lambda_ABab h_A h_B s_a s_b /4!, source V=-L',
        'Higgs_vev_H0R_index_zero_based':2,
        'native_HHSS_channels':list(NATIVE),
        'candidate_alias_hypotheses':BRANCHES,
        'lambdaS2Adj_classification':'pure-scalar quartic: zero heavy-scalar background gives no tree Higgs-VEV scalar quadratic term',
        'old_projected_legacy_replaced_by_native_source_HHSS_Hessian':True,
        'scenarios':complete,
        'conditional_EFT_mnu_eV_per_unknown_N':encode(eft),
        'fermion_parameters_at_muF_all_scenarios':True,
        'scalar_mass_quartics_at_muUV_or_muF_or_muS_as_labeled':True,
        'all_F_integrated_out_at_muF':True,
        'both_scalars_integrated_out_simultaneously_at_muS_in_Batch4':True,
        'native_quartic_alias_certified':False,
        'gammaCC_to_Weyl_signed_factor':None,
        'full_hard_fermion_one_loop_matching_included':False,
        'physical_C5_certified':False,
        'full_physical_one_loop_mnu_certified':False,
    }
    target.parent.mkdir(parents=True,exist_ok=True)
    staging=Path(tempfile.mkdtemp(prefix='.batch6_portal_stage_',dir=target.parent))
    try:
        (staging/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        lines=['# Batch 6 A_-4: native HHSS EWSB Hessian and alias branches','','**NOT PHYSICAL MATCHING.**',
               'The native Matchete scalar-quartic seed defines V4=-L4, fully symmetric lambda_ABab,',
               'with V4=lambda_ABab h_A h_B s_a s_b /4!. The neutral vev is h[2]=v.',
               'Thus Delta M2_ab = v^2 lambda_[2,2,a,b]/2 (real interleaved scalar basis).',
               'The full scalar EWSB mass Hessian replaces the provisional direct-UV scalar A/B.',
               'This remains conditional because the Matchete -> RGBeta quartic aliases are unverified.',
               'lambdaS2Adj is a pure-scalar quartic; with zero heavy scalar background it has no tree EWSB mass term.',
               '','| Scalar stage | lambdaH2Inv2 alias | UV masses (eV; descending) | Best complex-scale residual | Hessian delta vs Batch5 |',
               '|---|---|---|---:|---:|']
        for stage in STAGES:
            for branch in BRANCHES:
                m=complete[stage][branch]
                lines.append('| '+stage+' | '+BRANCHES[branch]+' | '+
                             ', '.join(f'{v:.8g}' for v in m['nu_mass_singular_values_eV_descending'])+
                             f" | {m['residual_after_one_complex_multiplier']:.6g} | " +
                             f"{m['native_vs_Batch5_legacy_EWSB_mass_hessian_relative']:.6g} |")
        lines.extend(['','**No physical GammaCC factor was fit.** Best complex scalar is matrix-shape diagnostic only.',
                      'UV scenario mixes muF fermions with frozen muUV scalars; S mixes muF fermions with muS EFT1 scalars.',
                      'Direct UV is broken-phase one-loop. EFT conditional C5 is unbroken-phase LLSS-induced fixed order.',
                      'Missing hard fermion one-loop matching, complete EFT1 mixing and source normalization.',
                      'EFT ordering: UV -> ALL F at muF -> EFT1(S1,S2) -> BOTH scalars at single muS -> SMEFT.',''])
        (staging/'report.md').write_text('\n'.join(lines),encoding='utf-8')
        for (stage,branch),data in records.items():
            d=staging/(stage+'_'+branch)
            d.mkdir()
            for k,v in data.items():
                (d/(k+'.json')).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        require(not target.exists(), 'Output became occupied during Batch6 calculation')
        staging.rename(target)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return {'model':'A_-4','status':summary['status'], 'output':str(target),
            'scenarios':{s:{b:{'residual_after_best_complex_scale':round(complete[s][b]['residual_after_one_complex_multiplier'],12),
                                 'UV_masses_ev':complete[s][b]['nu_mass_singular_values_eV_descending']}
                            for b in BRANCHES} for s in STAGES},
            'physical_C5_certified':False}
