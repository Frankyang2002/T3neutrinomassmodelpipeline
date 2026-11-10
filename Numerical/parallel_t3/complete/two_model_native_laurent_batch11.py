"""A(-4), D(-2): evaluate the *actual native Matchete scalar-soft Laurent*.

Source convention only. A successful hard+soft pole cancellation is NOT a
certificate of physical MSbar matching, GammaCC/Weyl normalization, a complete
EFT1 anomalous dimension, or of the full set of one-loop operators.

At strict O(hbar), use the F-scale values of dimensionless couplings in the
F and scalar-soft loop, since differences due to their one-loop RG evolution
multiply an already O(hbar) graph. The scalar spectrum is taken from muS.
The pole-derived EFT1 log is reconstructed from the *native* logarithm's
matching-scale derivative and inserted ONCE; no Batch9 scalar log is added.
"""
from __future__ import annotations

from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from typing import Any
import json
import math
import re

import numpy as np
import sympy as sp

from Numerical.parallel_t3.complete.two_model_joint_scalar_batch9 import (
    B1_SUBDIR, MODELS, load_snapshot, insist,
)
from Numerical.parallel_t3.complete.two_model_hard_general_batch10 import (
    matrix, matrix_json, symmetric, matrix_metrics,
)
from Numerical.parallel_t3.background.general_batch1_api import decode_native, read_json

FORMAT = 't3-thesis-two-model-native-hard-soft-laurent-batch11-v1'
B10 = Path('output/parallel_t3/thesis_hard_general_hhss_batch10_fix_v2')
SOFT_PATHS = {
    'A_-4': Path('output/parallel_t3/full_general_final_one_batch/'
                 'thesis_a4_fresh_f_v3/A_-4/A_scalar_soft_only_unrenormalized.txt'),
    'D_-2': Path('output/parallel_t3/full_general_final_one_batch/'
                 'thesis_physical_scalar_continuation_v1/D_-2/'
                 'D_scalar_soft_only_unrenormalized.txt'),
}
NATIVE_STATUS = {
    'A_-4': Path('output/parallel_t3/full_general_final_one_batch/'
                 'thesis_a4_fresh_f_v3/A_-4/native_two_threshold_status.json'),
    'D_-2': Path('output/parallel_t3/full_general_final_one_batch/'
                 'thesis_physical_scalar_continuation_v1/D_-2/'
                 'native_two_threshold_status.json'),
}
ONELOOP = 1.0/(16.0*math.pi**2)
EPS, MUBAR2 = sp.symbols('EPS MUBAR2', positive=True)
NAMES = {
    'T3EFT1W12': ('W12', True),
    'T3EFT1WExtra': ('WExtra', True),
    'lambdaT3': ('lam', False),
    'lambdaT3ExtraQ020101I1': ('extra', False),
    'MS1': ('m1', False),
    'MS2': ('m2', False),
    'T32MinusMass': ('mminus', False),
    'T32PlusMass': ('mplus', False),
}
ALL_SYMBOLS = {
    k:sp.Symbol(k) for k in (
        'W12', 'W12bar', 'WExtra', 'WExtrabar', 'lam', 'lambar',
        'extra', 'extrabar', 'm1', 'm2', 'mminus', 'mplus', 'hbar',
    )
}
ALL_SYMBOLS.update({'EPS':EPS, 'MUBAR2':MUBAR2})


def sha(path:Path)->str:
    insist(path.is_file(), 'Missing native source '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def _canonical_soft(text:str,model:str)->str:
    """Replace only documented native coupling objects; no guessed CG factors.

    Both native soft sources have already projected all SU2/Dirac structures
    to ordered holomorphic LLHH, and the only flavour labels are d$$1,d$$2.
    """
    insist(model in MODELS, 'Unknown model')
    s = re.sub(r'\s+', '', text)
    insist(40 < len(s) < 60000 and 'ϵ' in s and 'μbar2' in s and
           'T3EFT1W12' in s and 'hbar' in s,
           'Native soft source lacks regulator, scale, tree Wilson, or loop marker')
    allowed = set(('T3EFT1W12','lambdaT3','MS1','MS2')) if model=='A_-4' else set(
        ('T3EFT1W12','T3EFT1WExtra','lambdaT3',
         'lambdaT3ExtraQ020101I1','MS1','T32MinusMass','T32PlusMass'))
    for name,(symbol,has_flavour) in NAMES.items():
        suffix = '{Index[d$$1,Flavor],Index[d$$2,Flavor]}' if has_flavour else '{}'
        token = f'Coupling[{name},{suffix},0]'
        if name not in allowed:
            insist(token not in s, 'Unexpected coupling in '+model+': '+name)
            continue
        # Replace conjugate first to retain complex phases, without guessing
        # a second ordering/factor for the scalar or flavour indices.
        s = s.replace('Bar['+token+']', symbol+'bar')
        s = s.replace(token, symbol)
    insist('Coupling[' not in s and 'Bar[' not in s and
           'Index[' not in s and 'Field[' not in s and 'CG[' not in s,
           'Native soft source has unrecognised operator/coupling; block')
    s = (s.replace('μbar2','MUBAR2').replace('ϵ','EPS')
           .replace('Sqrt[','sqrt(').replace('Log[','log(')
           .replace('^','**').replace('[','(').replace(']',')'))
    insist(re.fullmatch(r'[0-9A-Za-z_+*/().-]+',s) is not None,
           'Unrecognised symbols/characters in native Laurent source')
    return s


@lru_cache(maxsize=4)
def parse_laurent(text:str, model:str)->tuple[sp.Expr,sp.Expr]:
    """Exact coefficient extraction: E = R + P/epsilon, no other powers."""
    source=_canonical_soft(text,model)
    from sympy.parsing.sympy_parser import parse_expr
    local={**ALL_SYMBOLS,'sqrt':sp.sqrt,'log':sp.log}
    # No file execution, arbitrary function evaluation or free variables.
    whole=parse_expr(source,local_dict=local,
        global_dict={'Integer':sp.Integer,'Rational':sp.Rational,'Float':sp.Float,
                     'Mul':sp.Mul,'Add':sp.Add,'Pow':sp.Pow},evaluate=True)
    insist(set(whole.free_symbols) <= set(ALL_SYMBOLS.values()),
           'Native soft expression uses an unknown symbolic parameter')
    # Extract Laurent coefficients as rational-function limits. Matchete
    # writes 1/(EPS*(m1^2-m2^2)*...) so expand().coeff(EPS,-1)
    # is WRONG: it expands the EPS-containing denominator into a sum.
    pole=sp.cancel(EPS*whole).subs(EPS,0)
    regular=sp.cancel(whole-pole/EPS).subs(EPS,0)
    insist(not (EPS in pole.free_symbols or EPS in regular.free_symbols),
           'Unsupported higher-dimensional-regulator Laurent order')
    # The unexpanded formula has no O(epsilon) remainder in either saved source.
    remainder=sp.cancel(whole - regular - pole/EPS)
    insist(remainder==0,
           'Native scalar soft expression has unrecognized epsilon dependence')
    insist(pole != 0 and regular != 0,
           'Native Laurent pole/regular coefficient vanished unexpectedly')
    return regular,pole


def _numeric_entry(expr:sp.Expr, symbols:dict[str,Any])->complex:
    mapping={ALL_SYMBOLS[k]:v for k,v in symbols.items() if k in ALL_SYMBOLS}
    mapping[MUBAR2]=symbols['MUBAR2']
    out=complex(sp.N(expr.subs(mapping),18))
    insist(np.isfinite([out.real,out.imag]).all(),
           'Nonfinite numeric native scalar Laurent coefficient')
    return out


def evaluate_soft(regular:sp.Expr,pole:sp.Expr, *,
                  w12:np.ndarray,wextra:np.ndarray,lam:complex,
                  extra:complex,m1:float,m2:float,
                  mminus:float,mplus:float,mu_bar:float)->tuple[np.ndarray,np.ndarray]:
    """All nine ordered flavour entries, no artificial flavour symmetrisation."""
    W=np.asarray(w12,complex);X=np.asarray(wextra,complex)
    insist(W.shape==X.shape==(3,3) and np.isfinite(W).all() and
           np.isfinite(X).all(), 'Expected two full 3x3 flavour Wilson matrices')
    masses=[m1,m2,mminus,mplus,mu_bar]
    insist(all(math.isfinite(float(x)) and x>0 for x in masses),
           'Nonpositive/nonfinite native matching scale or scalar masses')
    vals={'lam':complex(lam),'lambar':complex(lam).conjugate(),
          'extra':complex(extra),'extrabar':complex(extra).conjugate(),
          'm1':m1,'m2':m2,'mminus':mminus,'mplus':mplus,
          'MUBAR2':mu_bar**2,'hbar':ONELOOP}
    R=np.zeros((3,3),complex); P=np.zeros_like(R)
    for p in range(3):
        for q in range(3):
            args={**vals,'W12':W[p,q],'W12bar':W[p,q].conjugate(),
                  'WExtra':X[p,q],'WExtrabar':X[p,q].conjugate()}
            R[p,q]=_numeric_entry(regular,args)
            P[p,q]=_numeric_entry(pole,args)
    return R,P


def hard_from_pole_at_scale(h_regular:np.ndarray,h_pole:np.ndarray,
                            mu_f:float,new_mu_f:float)->np.ndarray:
    """Native hard log L=log(muF^2/MF^2) has derivative 2*hard_pole."""
    insist(mu_f>0 and new_mu_f>0,'F scales must be positive')
    return np.asarray(h_regular,complex)+2*math.log(new_mu_f/mu_f)*np.asarray(h_pole,complex)


def combine_native(h_regular:np.ndarray,h_pole:np.ndarray,
                   s_regular:np.ndarray,s_pole:np.ndarray,
                   mu_f:float,mu_s:float)->dict[str,np.ndarray|float]:
    """Strict-O(hbar) source algebra, no separately inserted B9 scalar log.

    From d soft-regular / d log(mu_bar)=2 soft-pole, the source-sector
    EFT1 log ensuring independence from the arbitrary split matching scales
    is -2 soft-pole log(muS/muF). This is not a certified full ADM.
    """
    insist(0<mu_s<mu_f,'Integrate F then both scalars, not reverse')
    hreg,hpol,sreg,spol=(np.asarray(z,complex) for z in
                        (h_regular,h_pole,s_regular,s_pole))
    insist(all(z.shape==(3,3) and np.isfinite(z).all() for z in
               (hreg,hpol,sreg,spol)),'Bad native 3x3 coefficient')
    pole_total=hpol+spol
    source_scale_log=-2*spol*math.log(mu_s/mu_f)
    finite=hreg+sreg+source_scale_log
    relative=float(np.linalg.norm(pole_total)/max(np.linalg.norm(hpol),
                                                  np.linalg.norm(spol),1e-100))
    return {'native_total_pole':pole_total,
            'source_pole_closure_relative':relative,
            'source_derived_eft1_scale_log':source_scale_log,
            'native_finite_source_sector':finite}


def _candidate_soft_values(root:Path,model:str):
    snap=load_snapshot(root,model)
    f=read_json(snap.b1_folder/'fermion_matching_input.json')
    a=f['native_uv_rgbeta_couplings_at_muF']
    masses=np.asarray(decode_native(a['MF']),complex)
    insist(masses.shape==(3,3),'Invalid F mass matrix')
    d=np.diag(masses)
    insist(np.min(d.real)>0 and np.max(abs(d.imag))<1e-8*np.max(d.real) and
           np.linalg.norm(masses-np.diag(d)) <= 1e-8*np.linalg.norm(masses),
           'Native F hard uses diagonal-mass approximation: out-of-scope MF')
    m=d.real
    Y1=np.asarray(decode_native(a['y1']),complex)
    Y2=np.asarray(decode_native(a['y2']),complex)
    insist(Y1.shape==Y2.shape==(3,3),'F UV Yukawa matrix must be 3x3')
    W=(Y1/m[None,:])@Y2.T
    X=np.zeros((3,3),complex)
    if model=='D_-2':
        Yextra=np.asarray(decode_native(a['y2Conj']),complex)
        insist(Yextra.shape==(3,3) and np.count_nonzero(Yextra)==0,
               'D extra heavy Yukawa nonzero: native hard module blocks unverified map')
        X=(Y1/m[None,:])@Yextra.T
    sf=snap.trajectory.at_fermion_threshold()
    ss=snap.trajectory.at_scalar_threshold()
    lam=complex(np.asarray(sf['lambdaT3']))
    ex=0j
    if model=='D_-2':
        ex=complex(np.asarray(sf['lambdaT3GeneralQ020101I1']))
    m1=float(np.sqrt(np.real(ss['mS1Sq'])))
    m2=float(np.sqrt(np.real(ss['mS2Sq'])))
    minus=plus=m2
    if model=='D_-2':
        delta=complex(np.asarray(ss['muT3GeneralS2Sq']))
        insist(abs(delta.imag)<1e-10*max(1.,abs(delta.real)),
               'Native D real mass basis does not support complex self-bilinear')
        m_minus_sq=m2**2-2*delta.real
        m_plus_sq=m2**2+2*delta.real
        insist(m_minus_sq>0 and m_plus_sq>0,'Tachyonic D split real scalar')
        minus,plus=math.sqrt(m_minus_sq),math.sqrt(m_plus_sq)
    return snap,{'w12':W,'wextra':X,'lam':lam,'extra':ex,
                 'm1':m1,'m2':m2,'mminus':minus,'mplus':plus}


def run_one(root:Path,model:str,*,batch10:Path=B10,
            soft_source:Path|None=None,native_status:Path|None=None,
            scale_probe_fraction:float=0.13,include_illustrative:bool=False)->dict:
    """Read saved native source and Batch10 F ledger. No other files changed."""
    insist(model in MODELS,'Unknown model')
    root=Path(root).resolve()
    b10=root/batch10/model/'result.json'
    old=read_json(b10)
    insist(old.get('model')==model and old.get('status')==
           'HardAndGeneralScalarSourceComponentsComputed_ConditionalNotPhysicalC5',
           'Successful Batch10 two-model hard/scalar output is required')
    native=root/(soft_source or SOFT_PATHS[model])
    status_file=root/(native_status or NATIVE_STATUS[model])
    report=read_json(status_file)
    insist(report.get('model')==model and report.get('status')==(
           'A_PhysicalScalarSoftOnly_HardInheritanceAndWilsonRGPending' if model=='A_-4'
           else 'D_PhysicalScalarSoftOnly_HardInheritanceAndWilsonRGPending') and
           report.get('joint_scalar_match_performed',model=='A_-4') is True,
           'Saved native joint scalar threshold did not succeed in the documented source')
    expected_terms=10 if model=='A_-4' else 68
    terms=report.get('A_physical_scalar_soft_Weinberg_terms',
                     report.get('physical_scalar_soft_Weinberg_terms'))
    insist(terms==expected_terms,
           model+': native scalar soft source term inventory changed')
    text=native.read_text(encoding='utf-8-sig')
    reg,pole=parse_laurent(text,model)
    snap,values=_candidate_soft_values(root,model)
    f,s=snap.running.mu_f,snap.running.mu_s
    hreg=matrix(old['F_hard_holomorphic_Weinberg_formal_symmetrized_regular'])
    hpol=matrix(old['F_hard_holomorphic_Weinberg_formal_symmetrized_pole'])
    soft_reg_ordered,soft_pole_ordered=evaluate_soft(reg,pole,**values,mu_bar=s)
    sreg=symmetric(soft_reg_ordered)
    spol=symmetric(soft_pole_ordered)
    combined=combine_native(hreg,hpol,sreg,spol,f,s)
    insist(combined['source_pole_closure_relative'] < 2e-7,
           model+': native F-hard + native scalar-soft pole cancellation FAILED')
    # Exact source identity: the coefficient of log(mu_bar^2) in soft regular
    # equals the epsilon pole residue; fixed coupling/mass values at one loop.
    shift=1.17
    shifted_reg,_=evaluate_soft(reg,pole,**values,mu_bar=s*shift)
    test_derivative=symmetric(shifted_reg)-sreg
    log_pred=2*spol*math.log(shift)
    log_diff=float(np.linalg.norm(test_derivative-log_pred)/max(
        np.linalg.norm(log_pred),np.linalg.norm(spol),1e-100))
    insist(log_diff<5e-7,model+': native scalar pole/log coefficient disagree')
    # Two independent matching-scale changes; source finite must not depend
    # on artificial placement of F/S thresholds at fixed parameters.
    new_f=f*(1+scale_probe_fraction)
    new_s=s*(1-scale_probe_fraction)
    shifted_soft,_=evaluate_soft(reg,pole,**values,mu_bar=new_s)
    shifted=combine_native(hard_from_pole_at_scale(hreg,hpol,f,new_f),
                           hpol,symmetric(shifted_soft),spol,new_f,new_s)
    invariance=float(np.linalg.norm(shifted['native_finite_source_sector']-
                                    combined['native_finite_source_sector'])/
                     max(np.linalg.norm(combined['native_finite_source_sector']),
                         np.linalg.norm(hreg),np.linalg.norm(sreg),1e-100))
    insist(invariance<5e-7,model+': native pole-derived scale partition inconsistent')
    result={
      'format':FORMAT,'model':model,
      'status':'NativeHardSoftLaurentAndScaleLogClosed_ConditionalOperatorConvention',
      'decoupling_order':[['F'],['S1','S2']],
      'scales_GeV':{'UV':1e7,'F':f,'S':s,'SMEFT_end':100},
      'native_physical_scalar_source_term_count':
        report.get('A_physical_scalar_soft_Weinberg_terms',
                   report.get('physical_scalar_soft_Weinberg_terms')),
      'source_soft_ordered_regular_GeVinv':matrix_json(soft_reg_ordered),
      'source_soft_ordered_pole_GeVinv':matrix_json(soft_pole_ordered),
      'source_soft_symmetric_regular_GeVinv':matrix_json(sreg),
      'source_soft_symmetric_pole_GeVinv':matrix_json(spol),
      'source_hard_symmetric_regular_GeVinv':matrix_json(hreg),
      'source_hard_symmetric_pole_GeVinv':matrix_json(hpol),
      'source_native_pole_sum_GeVinv':matrix_json(combined['native_total_pole']),
      'source_pole_closure_relative':combined['source_pole_closure_relative'],
      'batch10_scalar_vs_native_soft_same_convention_not_assumed':True,
      'batch10_candidate_full_7channel_status':(old.get(
          'A4_full_seven_channel_candidate_at_muF') or {}).get('status')
          if model=='A_-4' else 'not_computed_for_D2',
      'native_pole_log_identity_relative':log_diff,
      'source_EFT1_scale_partition_log_GeVinv':matrix_json(
          combined['source_derived_eft1_scale_log']),
      'source_finite_sector_candidate_GeVinv':matrix_json(
          combined['native_finite_source_sector']),
      'source_finite_sector_candidate_Frobenius':float(np.linalg.norm(
          combined['native_finite_source_sector'])),
      'source_scale_partition_invariance_relative':invariance,
      'scale_probe':{'F_prime_GeV':new_f,'S_prime_GeV':new_s},
      'F_and_soft_dimensionless_couplings_frozen_at_F_strict_one_loop':True,
      'scalar_spectrum_taken_from_same_muS_EFT1_state':True,
      'D2_extra_quartic_native_alias_candidate_not_certified': model=='D_-2',
      'D2_extra_y2Conj_zero_in_this_benchmark':model=='D_-2',
      'full_seven_term_EFT1_Wilson_tensor_RGE_included':False,
      'full_finite_F_matching_counterterms_certified':False,
      'complete_one_loop_scalar_operator_inventory_certified':False,
      'GammaCC_to_Weyl_normalization_certified':False,
      'physical_MSbar_C5_claimed':False,
      'full_18_real_yukawa_fit_performed':False,
      'native_source_hash_sha256':sha(native),
      'native_status_hash_sha256':sha(status_file),
      'batch10_result_sha256':sha(b10),
      'unmatched_F_operator_inventory_not_promoted_to_C5':True,
      'other_renormalized_operator_counterterms_not_included':True,
      'comment':'Pole/log closure is the actual source LLHH channel only. '
        'Formal finite matrix is NOT a complete physical C5 and is not '
        'used as the project final neutrino mass input.',
    }
    if include_illustrative:
        # The source normalization is NOT known; retaining the original
        # source coefficient here allows a clearly labelled sensitivity.
        from Numerical.parallel_t3.complete.full_general_completion import (
            existing_evidence,finish_smeft,
        )
        from scripts.RunT3ThesisFinalOneBatch import BENCHMARK, B1, B2, LAURENT
        evidence=existing_evidence(model,root/B1/model,root/B2/model,root/LAURENT/model)
        cfg=read_json(root/BENCHMARK)
        target_path=root/cfg['fitting']['target']
        target=read_json(target_path)
        example=finish_smeft(evidence,combined['native_finite_source_sector'],
                             vev=float(cfg['vev_GeV']),target=target)
        result['illustrative_smeft_NOT_PHYSICAL']=example
    return result
