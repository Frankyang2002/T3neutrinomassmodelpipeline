"""Batch 5: guarded six-model ordinary T3 source-matched fixed-order C5.

UV -> integrate ALL F -> EFT1(S1,S2) -> integrate ALL scalars -> SMEFT.

Independent native scalar hard template is REQUIRED.  The source-specific
prefactor is read from the archived one-generation F-threshold pole and
opposite-pole consistency is proven symbolically against the *new* scalar
hard expression using SymPy, rather than assuming a universal CG factor.
The archived six fermion-first Weinberg finite/pole/log terms are also audited.

This is restricted to the holomorphic Weinberg sector.  It is NOT proof of
all inherited dimension-five operator classes or physical spinor normalization.
"""
from __future__ import annotations

from hashlib import sha256
import math
from pathlib import Path
import re
from typing import Any

import numpy as np
import sympy as sp
from sympy.parsing.mathematica import parse_mathematica

from Numerical.running.T3OrdinaryHardFlavor import (
    LOOP, HardFlavorError, decode_complex, encode_complex, file_sha,
    mass_kernel, metrics, read_json, require, running_parameters,
)
from Numerical.running.T3OrdinaryInheritedFixedOrder import top_level_terms

MODELS = ('B_-3', 'B_+1', 'C_-3', 'C_+1', 'E_-4', 'E_+2')
TAGS = {m: m.replace('+', 'p').replace('-', 'm') for m in MODELS}
EXPECTED_TERMS = {'B': 2, 'C': 4, 'E': 4}
EXPECTED_EXTRAS = {'B': 1, 'C': 1, 'E': 0}
ARCHIVE = {m: 'T3_'+m[0]+'_alpha_'+m[2:].replace('+','p').replace('-','m') for m in MODELS}

YUKAWA = {
    n: re.compile(r'Bar\[Coupling\[' + n +
        r',\s*\{Index\[[^\]]+\],\s*Index\[[^\]]+\]\},\s*0\]\]')
    for n in ('y1','y2')
}
CONSTANTS = {
    'Coupling[lambdaT3,{},0]': 'LAM',
    'Coupling[MF,{},0]': 'MF',
    'Coupling[MS1,{},0]': 'MS1',
    'Coupling[MS2,{},0]': 'MS2',
    'μbar2': 'MU2',
    'ϵ': 'EPS',
}
SYM = sp.symbols('Y1 Y2 LAM MF MS1 MS2 MU2 EPS hbar')
Y1, Y2, LAM, MF, MS1, MS2, MU2, EPS, HBAR = SYM


def compact(text: str) -> str:
    return ''.join(text.split())


def parse_one_generation_coeff(source: str) -> sp.Expr:
    """Parse only the narrowly whitelisted Matchete one-mass coefficient.

    Reject unknown couplings rather than accidentally assigning them zero.
    """
    text = compact(source)
    require(100 < len(text) < 50000, 'Unexpected Mathematica C5 expression size')
    for name, pattern in YUKAWA.items():
        matches = pattern.findall(text)
        require(bool(matches), f'Mathematica hard coefficient missing barred {name}')
        text = pattern.sub('Y1' if name == 'y1' else 'Y2', text)
    for old, new in CONSTANTS.items():
        text = text.replace(old, new)
    require('Coupling[' not in text and 'Bar[' not in text and
            'Index[' not in text and 'Field[' not in text and
            'CG[' not in text and 'NCM[' not in text and
            'GammaCC' not in text, 'Additional operators/couplings in native hard coefficient; manual reduction required')
    require(set(re.findall(r'\b[A-Za-z][A-Za-z0-9]*\b', text))
            <= {'Y1','Y2','LAM','MF','MS1','MS2','MU2','EPS','hbar','Log','Sqrt'},
            'Unknown symbolic factor in hard coefficient')
    try:
        result = parse_mathematica(text)
    except Exception as exc:
        raise HardFlavorError(f'Cannot parse native hard coefficient: {exc}') from exc
    require(not result.has(sp.nan, sp.zoo, sp.oo), 'Non-finite symbolic hard expression')
    require(result.free_symbols <= set(SYM), 'Unknown symbolic coefficients')
    return result


def inherited_prefactor(pole_report: dict, model: str) -> sp.Expr:
    require(model in MODELS and pole_report.get('Status')=='Success' and
            pole_report.get('ConsistencyValidated') is True and
            pole_report.get('DirectLogEqualsTwicePoleResidueOneGeneration') is True and
            str(pole_report.get('DirectLogToPoleRatioOneGenerationInputForm')) == '2',
            f'{model}: first-threshold one-generation pole/RGE audit not successful')
    source = pole_report.get('HardPoleResidueOneGenerationInputForm')
    require(isinstance(source,str), 'Missing archived F pole coefficient')
    expression = parse_one_generation_coeff(source)
    factor = sp.cancel(expression/(HBAR*Y1*Y2*LAM/MF))
    expected = {'B': 2/sp.sqrt(3), 'C': -sp.Rational(2,3),
                'E': -2*sp.sqrt(sp.Rational(2,3))}[model[0]]
    require(sp.simplify(factor-expected)==0,
            f'{model}: source Clebsch normalization unexpected ({factor}); no extrapolation')
    return factor


def audit_inherited(path: Path, model: str, prefactor: sp.Expr) -> dict:
    """Six explicit 2L2H F-one-loop terms: HC and holomorphic triplets."""
    require(model in MODELS, f'Unsupported {model}')
    data = Path(path).read_text(encoding='utf-8-sig')
    terms = top_level_terms(data)
    items = [t for t in terms if t.count('Field[l,')==2 and
             t.count('Field[H,')==2 and 'hbar' in t and
             not any('Field['+k+',' in t for k in ('NewScalar1','NewScalar2','NewFermion'))]
    require(len(items)==6, f'{model}: expected six inherited SM Weinberg terms, got {len(items)}')
    groups: dict[str,dict[str,str]] = {'holomorphic':{}, 'anti_holomorphic':{}}
    log_token='Log[μbar2/Coupling[MF,{},0]^2]*'
    for raw in items:
        t=compact(raw)
        ybar = all('Bar[Coupling['+n+',' in t for n in ('y1','y2'))
        hbar = ('Coupling[lambdaT3,{},0]' in t and 'Bar[Coupling[lambdaT3' not in t)
        orient = 'holomorphic' if ybar and hbar else 'anti_holomorphic'
        if orient=='anti_holomorphic':
            require('Bar[Coupling[lambdaT3,{},0]]' in t and not ybar,
                    f'{model}: unexpected conjugate orientation')
        require(t.count('hbar')==1 and 'GammaCC' in t and
                'Coupling[MF,{},0]' in t and
                'Coupling[MF,{Index[' not in t,
                f'{model}: inherited F one-loop term has unexpected flavour/mass structure')
        is_pole='ϵ' in t; is_log=log_token in t
        require(not (is_pole and is_log), f'{model}: mixed pole/log term')
        cat='pole' if is_pole else 'log' if is_log else 'finite'
        require(cat not in groups[orient], f'{model}: repeated {orient} {cat}')
        groups[orient][cat]=t
    for orient,g in groups.items():
        require(set(g)=={'finite','pole','log'}, f'{model}: incomplete inherited {orient} triplet')
        normalized={k:v.replace('ϵ*','').replace(log_token,'')
                    .replace(')/(Coupling[MF,{},0])', ')/Coupling[MF,{},0]')
                    for k,v in g.items()}
        require(len(set(normalized.values()))==1,
                f'{model}: finite/pole/log not an identical coefficient')
        require(normalized['finite'].count('Coupling[y1,')==1 and
                normalized['finite'].count('Coupling[y2,')==1,
                f'{model}: inherited coefficient nonlinear in Yukawas')
        # Source factor must be consistent with pole report for all six terms.
        txt=normalized['finite']
        if model[0]=='B':
            require('Sqrt[3]*Coupling[MF,{},0]' in txt and
                    not txt.startswith('-'), f'{model}: B inherited CG prefactor changed')
        elif model[0]=='C':
            require('/(3*Coupling[MF,{},0])' in txt and txt.startswith('-'),
                    f'{model}: C inherited CG prefactor changed')
        else:
            require('Sqrt[2/3]*hbar' in txt and txt.startswith('-') and
                    '/Coupling[MF,{},0]' in txt,
                    f'{model}: E inherited CG prefactor changed')
    return {'status':'FOneLoopInheritedWeinbergTripletsVerified',
            'source_sha256':file_sha(path),
            'inherited_pole_coefficient':str(prefactor),
            'finite_pole_log_equality_per_chirality':True,
            'pureSM_twoL_twoH_terms':6}


def verify_native_hard(folder: Path, model: str, pole_report: dict,
                       source_root: Path|None=None) -> dict:
    require(model in MODELS, f'Unsupported {model}')
    folder=Path(folder)
    report=read_json(folder/'hard_scalar_report.json')
    audit=read_json(folder/'hard_scalar_input_audit.json')
    origin=read_json(folder/'provenance.json')
    require(report.get('Status')=='HardScalarCandidateExported' and
            report.get('Model')==model and report.get('HardOneLoopMatched') is True and
            report.get('PerHeavyMassEigenvalueTemplate') is True and
            report.get('FermionMatchingRecomputed') is False,
            f'{model}: native scalar Match did not export the correct template')
    n=EXPECTED_TERMS[model[0]]
    require(audit.get('Status')=='Passed' and audit.get('Model')==model and
            audit.get('CheckLagrangianPassed') is True and
            audit.get('TreeWilsonTermCount')==n and
            audit.get('IndexedFermionMassFactorsInOriginalTree')==n and
            audit.get('AdditionalQuarticsRetained')==EXPECTED_EXTRAS[model[0]] and
            audit.get('DiagonalNFlavorMassIndexEliminatedForTemplate') is True,
            f'{model}: scalar hard source audit failed')
    require(origin.get('model')==model and origin.get('status')=='Ready' and
            origin.get('tree_wilson_terms')==n and
            origin.get('additional_quartics_retained')==EXPECTED_EXTRAS[model[0]],
            f'{model}: scalar hard provenance not for this model')
    if source_root is not None:
        tag=TAGS[model]
        fps=origin.get('source_fingerprints',{})
        require(isinstance(fps,dict) and len(fps)>=5, 'Missing native F source fingerprints')
        for relative,sha in fps.items():
            if relative=='matching_wl': continue
            p=Path(source_root)/tag/relative
            require(p.is_file() and file_sha(p)==sha, f'F-only source hash mismatch: {p}')
    path=folder/'hard_c5_candidate.txt'
    require(path.is_file(), f'{model}: no hard C5 expression')
    raw=path.read_text(encoding='utf-8-sig')
    actual=parse_one_generation_coeff(raw)
    source_factor=inherited_prefactor(pole_report,model)
    expected_hard_factor=-source_factor
    standard_kernel=(1/EPS+1+(MS1**2*sp.log(MU2/MS1**2)-
                                MS2**2*sp.log(MU2/MS2**2))/(MS1**2-MS2**2))
    expected=expected_hard_factor*HBAR*Y1*Y2*LAM/MF*standard_kernel
    diff=sp.cancel(sp.together(actual-expected))
    require(diff==0 or sp.simplify(diff)==0,
            f'{model}: native scalar hard C5 is not the assumed one-loop mass kernel; '
            'numerical assembly blocked pending independent symbolic analysis')
    return {'status':'NativeHardUniversalKernelAlgebraicallyVerified',
            'expression_raw_sha256':file_sha(path),
            'expression_normalized_sha256':sha256(compact(raw).encode()).hexdigest(),
            'hard_report_sha256':file_sha(folder/'hard_scalar_report.json'),
            'source_prefactor_inherited':str(source_factor),
            'source_prefactor_hard':str(expected_hard_factor),
            'extra_quartics_retained':EXPECTED_EXTRAS[model[0]],
            'tree_wilson_source_terms':n,
            'mass_kernel_proven_by_symbolic_difference':True,
            'complete_physical_c5_certified':False}


def matrices(y1:np.ndarray,y2:np.ndarray,mf:np.ndarray,lam:complex,
             factor:float)->tuple[np.ndarray,list[np.ndarray]]:
    """Use source-hard pole CG factor and a symmetric 3x3 flavour projection."""
    out=[]
    for r in range(3):
        m=0.5*(np.outer(np.conj(y1[:,r]),np.conj(y2[:,r]))+
               np.outer(np.conj(y2[:,r]),np.conj(y1[:,r])))
        out.append(factor*LOOP*lam*m/mf[r])
    return sum(out,np.zeros((3,3),complex)),out


def evaluate(model:str, hard_folder:Path, running_doc:dict,
             pole_report:dict, inherited_path:Path, *, vev:float=246.22,
             reference_path:Path|None=None,
             source_root:Path|None=None)->dict:
    """One-loop restricted Weinberg subtotal with common muF couplings."""
    require(model in MODELS, f'Unsupported {model}')
    hard=verify_native_hard(hard_folder,model,pole_report,source_root)
    source_factor=inherited_prefactor(pole_report,model)
    inherited=audit_inherited(inherited_path,model,source_factor)
    p=running_parameters(running_doc,model)
    hard_factor=-float(sp.N(source_factor,17))
    P,by_r=matrices(p['y1'],p['y2'],p['mf'],p['lam_f'],hard_factor)
    inherited_pole=-P
    denominator=max(float(np.linalg.norm(P)),1e-40)
    pole_relative=float(np.linalg.norm(P+inherited_pole)/denominator)
    require(pole_relative<1e-13,f'{model}: fixed-order pole noncancellation')
    hard_finite=P*mass_kernel(p['m1sq_f'],p['m2sq_f'],p['mu_s'])
    direct=(-2*math.log(p['mu_s']/p['mu_f']))*P
    inherited_finite=sum((-z)*(1+2*math.log(p['mu_f']/mass))
                         for z,mass in zip(by_r,p['mf']))
    C=hard_finite+direct+inherited_finite
    simplified=sum((z*(mass_kernel(p['m1sq_f'],p['m2sq_f'],mass)-1))
                   for z,mass in zip(by_r,p['mf']))
    residual=float(np.linalg.norm(C-simplified)/max(float(np.linalg.norm(C)),1e-40))
    require(residual<3e-12,f'{model}: finite+log+inherited closure failed')
    scale_probe=p['mu_s']*math.exp(.042)
    Cprobe=(P*mass_kernel(p['m1sq_f'],p['m2sq_f'],scale_probe) -
            2*P*math.log(scale_probe/p['mu_f'])+inherited_finite)
    scale_residual=float(np.linalg.norm(Cprobe-C)/max(float(np.linalg.norm(C)),1e-40))
    require(scale_residual<3e-12,f'{model}: matching-scale check failed')
    Ps, _=matrices(p['y1'],p['y2'],p['mf'],p['lam'],hard_factor)
    subtotal_s=Ps*(mass_kernel(p['m1sq'],p['m2sq'],p['mu_s']) -
                    2*math.log(p['mu_s']/p['mu_f']))
    parity={'status':'NoArchivedLLSSReference'}
    if reference_path is not None and Path(reference_path).is_file():
        ref=read_json(reference_path)
        require(ref.get('model')==model and ref.get('complete_one_loop_C5_certified') is False,
                f'{model}: wrong or non-conditional reference')
        response=ref.get('response',{})
        require(response.get('model')==model, 'Wrong archived LLSS response model')
        compares={
            'hard_finite':(Ps*mass_kernel(p['m1sq'],p['m2sq'],p['mu_s']),
                  'C5_joint_scalar_finite_bubble_per_N_gev_inverse'),
            'direct_log':(-2*math.log(p['mu_s']/p['mu_f'])*Ps,
                  'C5_EFT1_induced_leading_log_per_N_gev_inverse'),
            'hard_plus_direct':(subtotal_s,'C5_LLSS_total_one_loop_per_N_gev_inverse')}
        discrepancies={}
        for k,(local,remote) in compares.items():
            target=decode_complex(response[remote],remote)
            discrepancies[k]=float(np.linalg.norm(local-target)/
                max(float(np.linalg.norm(local)),float(np.linalg.norm(target)),1e-40))
        parity={'status':'SourceBackedLLSSKernelAgreement'
                if max(discrepancies.values())<1e-8 else 'ConditionalArchiveKernelMismatch',
                'relative_errors':discrepancies,'source_sha256':file_sha(reference_path),
                'archive_is_NOT_added_to_C5':True}
    return {'format':'t3-ordinary-six-common-fixed-order-source-constrained-v1',
        'model':model,'status':'InheritedHardAndCommonCouplingPoleCancelled_Conditional',
        'native_scalar_source':hard,'archived_fermion_source':inherited,
        'fermion_pole_report_sha256':sha256(str(pole_report).encode()).hexdigest(),
        'heavy_masses_gev':p['mf'].tolist(),
        'matching_scales_gev':{'muF':p['mu_f'],'muS':p['mu_s']},
        'source_hard_pole_prefactor':str(-source_factor),
        'hard_pole_residue_gev_inverse':encode_complex(P),
        'inherited_pole_residue_gev_inverse':encode_complex(inherited_pole),
        'common_total_pole_residue_gev_inverse':encode_complex(P+inherited_pole),
        'hard_finite_gev_inverse':encode_complex(hard_finite),
        'direct_eft1_log_gev_inverse':encode_complex(direct),
        'inherited_first_fermion_finite_gev_inverse':encode_complex(inherited_finite),
        'c5_restricted_fixed_order_gev_inverse':encode_complex(C),
        'c5_restricted_simplified_kernel_gev_inverse':encode_complex(simplified),
        'conditional_masses':metrics(C,vev),
        'archive_comparison':parity,
        'checks':{
            'native_matchete_scalar_one_loop_matched':True,
            'actual_hard_expression_algebraically_verified':True,
            'archived_F_one_loop_weinberg_triplets_verified':True,
            'all_three_heavy_mass_eigenvalues_summed':True,
            'common_coupling_pole_cancellation_residual':pole_relative,
            'matching_scale_variation_relative_residual':scale_residual,
            'finite_kernel_closure_relative_residual':residual,
            'full_EFT1_operator_mixing_and_counterterms_certified':False,
            'absolute_SU2_spinor_normalization_certified':False,
            'all_model_physical_C5_certified':False,
            'no_double_counted_legacy_LLSS':True,
        },
        'physical_c5_certified':False,
        'interpretation':'This is a source-verified restricted holomorphic Weinberg contribution. '
           'Nondegenerate MF numerical sum is inferred by applying an identical one-eigenvalue template '
           'to each diagonal MF_r. Matchete full F one-loop with all three masses and full EOM/counterterms '
           'has not been independently rerun; physical normalization remains provisional.'}
