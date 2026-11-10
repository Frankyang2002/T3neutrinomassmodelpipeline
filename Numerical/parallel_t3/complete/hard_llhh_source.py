"""Numerical FORMAL LLHH F-threshold source channel, not a physical C5.

Reads the actual complete Matchete one-loop operator ledger (six LLHH terms
for A_-4, twelve for D_-2), WITHOUT calling an analytic triangle function.
The selected terms have common raw factor 2/Sqrt[3] and organize into
constant, 1/eps and log(mu_bar^2/M_r^2).  GammaCC/chiral conventions, EFT1
soft subtraction and MSbar counterterms remain unresolved. These matrices
must never be substituted for a full matched C5(MS).

Source-supported comparison: the B1 mass matrix off-diagonals are tiny.
A diagonal-MF approximation is therefore allowed only below a hard relative
threshold 1e-8; otherwise this provisional source calculation blocks.
"""
from __future__ import annotations

import math
from pathlib import Path
import re
from typing import Any

import numpy as np

from Numerical.parallel_t3.background.general_batch1_api import decode_native,encode_native
from Numerical.parallel_t3.complete.full_general_completion import insist


class HardSourceError(ValueError):
    pass


_HARD_TERMS = {'A_-4':6, 'D_-2':12}
# Matchete holomorphic regular first-threshold terms have the following
# source string structures (verified in current committed operator ledgers).
_TYPES = ('regular_constant', 'explicit_pole', 'regular_log')


def _classify(term: dict, model: str) -> tuple[str,str,str]:
    text=term.get('raw_inputform')
    if not isinstance(text,str) or term.get('sector')!='CandidateLeptonNumberViolatingLLHH':
        raise HardSourceError(model+': invalid LLHH ledger entry')
    # Source has this exact coefficient for these operators; do not infer
    # arbitrary prefactors from unrelated graphs.
    if (not text.startswith('(2*hbar*') or text.count('hbar')!=1 or
            text.count('GammaCC')!=1 or
            text.count('Coupling[MF,')!=(2 if 'Log[' in text else 1) or text.count('Sqrt[3]')!=1):
        raise HardSourceError(model+': noncanonical hard source term; manual map needed')
    if 'Coupling[y1,' not in text or text.count('Coupling[y1,')!=1:
        raise HardSourceError(model+': changed y1 hard source')
    is_extra='Coupling[y2Conj,' in text
    ykind='y2Conj' if is_extra else 'y2'
    if text.count('Coupling['+ykind+',')!=1:
        raise HardSourceError(model+': changed y2 hard source')
    if is_extra:
        if model!='D_-2' or 'lambdaT3ExtraQ020101I1' not in text:
            raise HardSourceError('Unexpected extra-quartic hard channel')
        channel='extra_y2Conj'
    else:
        if 'Coupling[lambdaT3,' not in text:
            raise HardSourceError('Missing baseline quartic in hard channel')
        channel='baseline_y2'
    # Chiralities reported separately; never add H and Hdag terms together.
    if 'Proj[-1]' in text and 'Proj[1]' not in text:
        chirality='holomorphic_HH'
    elif 'Proj[1]' in text and 'Proj[-1]' not in text:
        chirality='anti_holomorphic_HdagHdag'
    else:
        raise HardSourceError('Ambiguous GammaCC chirality')
    has_log='Log[' in text
    has_pole=('ϵ' in text or 'ε' in text or '\\[CurlyEpsilon]' in text)
    if has_log and has_pole:
        raise HardSourceError('Mixed log/pole source not supported by this exact inventory')
    kind='regular_log' if has_log else 'explicit_pole' if has_pole else 'regular_constant'
    if has_log and 'Log[μbar2/Coupling[MF,' not in text:
        raise HardSourceError('Unknown scale-log source')
    if kind=='explicit_pole' and '/(Sqrt[3]*ϵ*Coupling[MF,' not in text:
        raise HardSourceError('Unexpected regulator/pole normalization')
    if kind!='explicit_pole' and '/(Sqrt[3]*Coupling[MF,' not in text:
        raise HardSourceError('Unexpected regular normalization')
    return channel,chirality,kind


def extract_inventory(model: str, ledger: dict) -> dict:
    if model not in _HARD_TERMS or ledger.get('status')!='UnrenormalisedOneLoopTermsEnumerated':
        raise HardSourceError('Unsupported or unsuccessful full Matchete ledger')
    entries=[t for t in ledger.get('terms',())
             if t.get('sector')=='CandidateLeptonNumberViolatingLLHH']
    if len(entries)!=_HARD_TERMS[model]:
        raise HardSourceError('LLHH term inventory changed from source')
    channels={}
    for term in entries:
        chan,chi,kind=_classify(term,model)
        key=(chan,chi)
        channel=channels.setdefault(key,[])
        channel.append(kind)
    expected=[('baseline_y2',c) for c in
              ('holomorphic_HH','anti_holomorphic_HdagHdag')]
    if model=='D_-2':
        expected += [('extra_y2Conj',c) for c in
                     ('holomorphic_HH','anti_holomorphic_HdagHdag')]
    if set(channels)!=set(expected) or any(sorted(channels[x])!=sorted(_TYPES) for x in expected):
        raise HardSourceError('Changed scalar/Yukawa/loop inventory in full source')
    return {f'{ch}:{chi}':tuple(sorted(channels[(ch,chi)])) for ch,chi in expected}


def numerical_unsubtracted_source(model: str, ledger: dict,
                                  native_at_F: dict, *, mu_f:float=1e5) -> dict:
    inventory=extract_inventory(model,ledger)
    values=native_at_F
    mf=np.asarray(decode_native(values['MF']),complex)
    y1=np.asarray(decode_native(values['y1']),complex)
    y2=np.asarray(decode_native(values['y2']),complex)
    if any(x.shape!=(3,3) or not np.isfinite(x).all() for x in (mf,y1,y2)):
        raise HardSourceError('Must have finite 3x3 native MF/y1/y2')
    diag=np.diag(mf)
    if np.min(np.abs(diag))<=0:
        raise HardSourceError('MF diagonal has zero mode')
    offdiag=float(np.linalg.norm(mf-np.diag(diag))/np.linalg.norm(mf))
    if offdiag>1e-8:
        raise HardSourceError(f'Non-diagonal MF needs documented chiral basis rotation: {offdiag:.4g}')
    if np.max(abs(diag.imag))>1e-8*np.max(abs(diag)):
        raise HardSourceError('Complex MF phases require mass basis rotation')
    mass=np.asarray(diag.real,float)
    if np.min(mass)<=0 or not math.isfinite(mu_f) or mu_f<=0:
        raise HardSourceError('Require three positive masses and matching scale')
    lam=complex(np.asarray(decode_native(values['lambdaT3'])))
    extra=np.zeros((3,3),complex)
    if model=='D_-2':
        extra=np.asarray(decode_native(values['y2Conj']),complex)
        if extra.shape!=(3,3):
            raise HardSourceError('D2 y2Conj must be 3x3')
        # Physical extra-quartic normalization is NOT yet a verified
        # Matchete<->RGBeta coordinate. Block rather than assuming equality
        # when this channel is active.
        if np.any(extra != 0):
            raise HardSourceError('D2 nonzero y2Conj needs independent Matchete/RGBeta extra-quartic normalization')
    ratio=1+np.log(mu_f**2/mass**2)
    if not np.isfinite(ratio).all():
        raise HardSourceError('Invalid heavy flavour log')
    # These are only formal RAW ordered GammaCC components; no C5 sign,
    # spinor normalization, or physical hard/soft subtraction is assigned.
    common=2/(math.sqrt(3)*16*math.pi**2)
    ordered_holo=common*lam* (y1.conj() * (ratio/mass)[None,:]) @ y2.conj().T
    ordered_anti=common*lam.conjugate()* (y1 * (ratio/mass)[None,:]) @ y2.T
    pole_holo=common*lam*(y1.conj()/mass[None,:]) @ y2.conj().T
    pole_anti=common*lam.conjugate()*(y1/mass[None,:]) @ y2.T
    return {
        'status':'NumericalRawHardLLHHSourceOnly_NotMSbarWilson',
        'model':model,
        'source_term_count':_HARD_TERMS[model],
        'source_channel_inventory':inventory,
        'matching_scale_GeV':float(mu_f),
        'MF_diagonal_coordinate_approximation_used':True,
        'MF_offdiagonal_Frobenius_relative':offdiag,
        'baseline_ordered_holomorphic_regular_GeVinv':encode_native(ordered_holo),
        'baseline_ordered_holomorphic_pole_GeVinv':encode_native(pole_holo),
        'baseline_ordered_anti_holomorphic_regular_GeVinv':encode_native(ordered_anti),
        'baseline_ordered_anti_holomorphic_pole_GeVinv':encode_native(pole_anti),
        'baseline_ordered_holomorphic_regular_Frobenius_GeVinv':float(np.linalg.norm(ordered_holo)),
        'baseline_ordered_holomorphic_pole_Frobenius_GeVinv':float(np.linalg.norm(pole_holo)),
        'pole_cancelled':False,
        'one_loop_MSbar_renormalized':False,
        'GammaCC_to_Weyl_normalization_verified':False,
        'full_EFT1_Wilson_RGE_performed':False,
        'scalar_threshold_matched':False,
        'physical_C5_computed':False,
        'formula_source':'Each verified raw GammaCC LLHH ledger channel consists of 2*hbar/Sqrt[3] times (1, 1/epsilon, Log[mu_bar^2/MF[r]^2]), times y1*y2/MF[r].',
        'warning':'Only a formal Matchete ordered field coefficient. Not a neutrino mass and not a finite EFT Wilson.'
    }
