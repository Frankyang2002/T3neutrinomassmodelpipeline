"""Independent common-log coefficient audit of general T3 one-scale Matchete C5.

The one-scale coefficient is evaluated twice, with every explicitly sourced
Log[MF^2/MScalar^2] replaced by Log[2*MF^2/MScalar^2] in the second evaluation.
The difference divided by log(2) gives the one-scale *common-log spurion*
coefficient with masses and couplings held fixed (for the one-loop log-linear
expressions shipped in the current repository).

It is compared to EFT1's scalar-pair prediction C5_LLSS_log/(2 log(muF/muS)).
The factor 2 converts log(MF^2/MS^2) to log(MF/MS). All operator signs and
unknown Matchete GammaCC -> Weyl normalisations are left as measured results;
no hard-matching subtraction is used anywhere in this audit.

UV -> integrate ALL heavy fermions -> scalar EFT1 -> integrate ALL scalars -> SMEFT.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Any, Mapping

import numpy as np

LOG_RATIO = re.compile(
    r"\s*Coupling\[MF,\s*\{\},\s*0\]\^2\s*/\s*"
    r"Coupling\[(MS1|MS2|MSX|MSY|T3[12](?:Minus|Plus)Mass),\s*\{\},\s*0\]\^2\s*"
)
MODELS = tuple(f'{f}_{a:+d}' for f, alphas in {
    'A':(-4,-2,0), 'B':(-3,-1,1), 'C':(-3,-1,1),
    'D':(-2,0,2), 'E':(-4,-2,0,2)}.items() for a in alphas)
MIXED = frozenset(('B_-1','C_-1'))
FORMAT = 't3-general-matchete-independent-common-log-audit-v1'


def require(value: Any, message: str) -> None:
    if not bool(value):
        raise ValueError('T3 independent Matchete log audit: '+message)


def read_json(path: Path) -> dict:
    doc=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    require(isinstance(doc,dict),f'{path}: JSON object required')
    return doc


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def replace_common_mass_logs(source: str, *, multiplier: int = 2) -> tuple[str, dict]:
    """Replace only FULL-THEORY heavy/scalar mass-ratio logs; reject all others.

    For a one-loop expression *linear* in its logarithms, the finite difference
    (C5_shifted-C5_original)/log(multiplier) is the exact coefficient of a
    simultaneous additive shift of all explicitly present mass-ratio logs.
    This does NOT rescale MF, the scalar masses, Yukawas, or physical scales.
    """
    require(multiplier >= 2 and isinstance(multiplier,int),
            'log spurion multiplier must be an integer >= 2')
    require(isinstance(source,str) and 0 < len(source) <= 2_000_000,
            'missing or oversized Matchete coefficient')
    output=[]; i=0; n=0; counts={}
    while True:
        at=source.find('Log[',i)
        if at < 0:
            output.append(source[i:]); break
        output.append(source[i:at])
        start=at+4
        j=start; depth=1
        while j < len(source) and depth:
            if source[j] == '[': depth+=1
            elif source[j] == ']': depth-=1
            j+=1
        require(depth==0,'unbalanced Log[...] in Matchete coefficient')
        argument=source[start:j-1]
        m=LOG_RATIO.fullmatch(argument)
        require(m is not None,
                'non-heavy/scalar mass ratio Log[...] encountered: '+argument[:180])
        scalar=m.group(1)
        counts[scalar]=counts.get(scalar,0)+1
        output.append(f'Log[{multiplier}*({argument})]')
        i=j; n+=1
    require(n>0,'Matchete C5 contains no recognised heavy/scalar ratio logarithms')
    transformed=''.join(output)
    return transformed, {'log_occurrences':n,'scalar_mass_symbols':counts,
                         'common_log_multiplier':multiplier,
                         'additive_log_shift':math.log(multiplier),
                         'every_original_log_is_MF_over_scalar_mass':True}


def decode_matrix(raw: Any) -> np.ndarray:
    from Numerical.running.GeneralT3SequentialLeadingLog import decode
    arr=np.asarray(decode(raw),complex)
    require(arr.shape==(3,3) and np.isfinite(arr).all(),
            'expected finite complex 3x3 C5 matrix')
    return arr


def encode_matrix(arr: np.ndarray) -> dict:
    from Numerical.running.GeneralT3SequentialLeadingLog import encode
    arr=np.asarray(arr,complex)
    require(arr.shape==(3,3) and np.isfinite(arr).all(),
            'invalid complex 3x3 matrix')
    return encode(arr)


def compare_coefficients(matchete_derivative: np.ndarray,
                         eft_one_log_coefficient: np.ndarray) -> dict:
    """Direct factor-1 check plus best-fit complex signed normalisation.

    The fitted factor is DIAGNOSTIC, not an allowed automatic correction.
    A small best-fit residual with factor far from 1 diagnoses normalization;
    a large residual diagnoses a missing or basis-mismatched matrix component.
    """
    x=np.asarray(eft_one_log_coefficient,complex)
    y=np.asarray(matchete_derivative,complex)
    require(x.shape==y.shape==(3,3) and np.isfinite(x).all() and
            np.isfinite(y).all(),'both log coefficients must be finite 3x3')
    xx=float(np.vdot(x,x).real)
    yn=float(np.linalg.norm(y))
    xn=float(np.linalg.norm(x))
    require(xx>1e-75 and yn>1e-38,'zero or tiny log coefficient: no normalization determined')
    factor=complex(np.vdot(x,y)/xx)
    return {
        'matchete_spurion_frobenius_gev_inverse':yn,
        'eft_spurion_frobenius_gev_inverse':xn,
        'matchete_over_eft_norm':yn/xn,
        'best_fit_complex_factor':{'real':float(factor.real),'imag':float(factor.imag)},
        'best_fit_factor_phase_rad':float(np.angle(factor)),
        'relative_matrix_difference_at_assumed_factor_1':float(np.linalg.norm(y-x)/yn),
        'relative_matrix_residual_after_best_fit_factor':float(np.linalg.norm(y-factor*x)/yn),
    }


def original_one_scale_parity(original: np.ndarray, anchored: np.ndarray) -> float:
    a=np.asarray(original,complex);b=np.asarray(anchored,complex)
    require(a.shape==b.shape==(3,3) and np.isfinite(a).all() and
            np.isfinite(b).all(),'one-scale matrix parity shape/finite')
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-38))


def saved_case_input(path: Path, model: str) -> dict:
    doc=read_json(path)
    require(doc.get('model')==model and doc.get('status')=='Computed',
            f'{model}: expected completed general_seq case record')
    cases=doc.get('cases')
    require(isinstance(cases,dict) and len(cases)==16 and '1111' in cases,
            f'{model}: need all 16 saved ablation configurations')
    point=cases['1111']
    require(point.get('case')=='1111' and
            all(point.get(k) is True for k in
                ('uv_on','eft1_on','wilson_log_on','smeft_on')),
            'saved full-running reference flags differ')
    full=decode_matrix(point['c5_matching_anchored_reference'])
    eft=decode_matrix(point['c5_induced_wilson_log'])
    mu_f=1e5; mu_s=1e3
    log_range=math.log(mu_f/mu_s)
    require(log_range>0, 'invalid EFT scale ratio')
    return {'full':full,'eft_spurion':eft/(2*log_range),
            'eft_log':eft,'mu_f_gev':mu_f,'mu_s_gev':mu_s,
            'log_range':log_range,
            'saved_spinor_factor':point.get('spinor_normalization_assumption',None),
            'case_record_sha256':sha256_file(path)}


def _nonmixed_parameter_data(model: str, fit: dict, coefficient: Path):
    """Use saved fit's ACTUAL muF UV and muS EFT1 endpoints; no new RGE run."""
    from Numerical.matching.GeneralT3MultiModelNeutrinos import (
        COUPLING_PATTERN,YUK,kernel_parameters,
    )
    from Numerical.matching.GeneralT3FlavourC5 import (
        prepare_heavy_flavour,validate_real_self_bilinear,
        signed_real_mass_basis_parameters,SELF,
    )
    from Numerical.running.GeneralT3DynamicFermionWilson import uv_from_json
    details=fit.get('numerical_details',{})
    require(isinstance(details,dict) and 'uv_at_fermion_threshold' in details
            and 'eft1_at_scalar_threshold' in details,
            f'{model}: fitted result lacks actual muF/muS state snapshots')
    uv=details['uv_at_fermion_threshold']
    eft=details['eft1_at_scalar_threshold']
    source=coefficient.read_text(encoding='utf-8-sig')
    symbols=set(COUPLING_PATTERN.findall(source))
    required=symbols&YUK
    require({'y1','y2'}<=required,'source lacks topology Yukawas')
    masses,yukawas,basis_diag=prepare_heavy_flavour(uv_from_json(model,uv),required)
    kappa=None
    if model in SELF:
        species='2' if model.endswith('-2') else '1'
        raw=eft.get('muT3GeneralS'+species+'Sq')
        require(raw is not None,f'{model}: missing self-scalar bilinear')
        kappa=complex(raw['real'],raw['imag']) if isinstance(raw,dict) else complex(raw)
        kappa=validate_real_self_bilinear(model,kappa)
    points=[]; mapping=[]
    for i, mass in enumerate(masses):
        couplings,notes=kernel_parameters(symbols,model,uv,eft,float(mass),
                                         kappa_override=kappa,strict_mapping=False)
        if model in SELF:
            couplings=signed_real_mass_basis_parameters(model,couplings,kappa,
                                                        coefficient_symbols=symbols)
        points.append({'index':i,'couplings':couplings})
        mapping.extend(notes)
    require(len(points)==3,'expected three heavy fermion matching eigenstates')
    return source,points,yukawas,sorted(set(mapping)),basis_diag


def _run_wolfram_channels(coefficient: Path,points: list[dict],
                          yukawas: Mapping[str,np.ndarray],*,
                          wolframscript: str, timeout_sec: int) -> np.ndarray:
    from Numerical.matching.GeneralT3FlavourC5 import (
        EVALUATOR_WL,flavour_c5_from_channels,
    )
    payload={'points':points,'hbar_value':1/(16*math.pi**2)}
    with tempfile.TemporaryDirectory(prefix='t3_matchete_log_spurion_') as temp:
        inputs=Path(temp)/'points.json';outputs=Path(temp)/'channels.json'
        inputs.write_text(json.dumps(payload,allow_nan=False),encoding='utf-8')
        command=[wolframscript,'-file',str(EVALUATOR_WL),
                 str(coefficient),str(inputs),str(outputs)]
        proc=subprocess.run(command,cwd=EVALUATOR_WL.parents[1],text=True,
                            capture_output=True,timeout=timeout_sec,check=False)
        require(proc.returncode==0 and outputs.is_file(),
                'Wolfram Matchete channel evaluation failed ('+
                str(proc.returncode)+'): '+(proc.stdout+proc.stderr)[-2500:])
        result=read_json(outputs)
        return flavour_c5_from_channels(result['results'],yukawas,result['channels'])


def independent_nonmixed(model: str,fit: dict,coefficient: Path,*,
                         wolframscript: str,timeout_sec: int,
                         hierarchy_factors: tuple[int,...]) -> tuple[list[dict],dict]:
    """Evaluate original/shifted exact one-scale expression at each hierarchy."""
    source,points,yukawas,mappings,basis=_nonmixed_parameter_data(model,fit,coefficient)
    shifted, log_inventory=replace_common_mass_logs(source)
    results=[]
    with tempfile.TemporaryDirectory(prefix='t3_matchete_shifted_coefficient_') as temp:
        alternative=Path(temp)/'shifted_c5_coefficient.txt'
        alternative.write_text(shifted,encoding='utf-8')
        for scale in hierarchy_factors:
            scaled=[{'index':p['index'],'couplings':dict(p['couplings'],
                           MF=float(p['couplings']['MF'])*scale)} for p in points]
            a=_run_wolfram_channels(coefficient,scaled,yukawas,
                     wolframscript=wolframscript,timeout_sec=timeout_sec)
            b=_run_wolfram_channels(alternative,scaled,yukawas,
                     wolframscript=wolframscript,timeout_sec=timeout_sec)
            results.append({'factor':scale,'original':a,
                            'spurion':(b-a)/math.log(2),
                            'delta_2_multiplier':b-a})
    return results,{'one_scale_source_sha256':sha256_file(coefficient),
                   'one_scale_path':str(coefficient),
                   'log_inventory':log_inventory,
                   'yukawa_channel_basis':'same GeneralT3FlavourC5 channel assembly as fitted one-scale',
                   'candidate_parameter_mappings':mappings,
                   'heavy_basis_diagnostics':basis,
                   'heavy_fermion_masses_gev':[float(p['couplings']['MF']) for p in points]}


def _mixed_running_endpoint(model: str, saved: Path) -> tuple[dict,tuple,dict]:
    """Reconstruct ONLY U=E=1 mixed benchmark using original solver modules."""
    from scripts import RunT3BCGeneralFromReference as runner
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc
    from scripts import RunT3GeneralScalarBCPilot as pilot
    from scripts import RunT3ZeroExtraBoundaryBenchmark as core
    from scripts.RunT3GeneralSequentialFitAblation import mixed_fit_fractions
    from Numerical.experimental.GeneralT3ScalarRGState import (
        unpack_general_scalar_state,project_general_after_fermion,
        pack_general_scalar_state,
    )
    fit=read_json(saved/'result.json')
    q,m=mixed_fit_fractions(fit)
    source=runner.build_reference_source(model)
    source['cfg']=read_json(saved/'fitted_config.json')
    reg,new,reused=source['registry'],source['new'],source['reused']
    betas=bc.beta_bundle(model,source,reg)
    y0,uv_layout,_,_=pilot.prepare_initial(source,reg,
                    'general_quartics_and_mixing' if m>0 else 'general_quartics',q,m)
    uvsol=core.integrate_segment(y0,
        lambda t,y:core.rhs_at(t,y,layout=uv_layout,registry=reg,
                              betas=betas['UV'],stage='UV'),
        1e7,1e5,n_points=3)
    end_uv=uvsol.y[:,-1].copy()
    heavy_state=unpack_general_scalar_state(
             end_uv[:uv_layout.size],uv_layout,reg,uv=True,mu_gev=1e5).baseline
    projected=project_general_after_fermion(unpack_general_scalar_state(
             end_uv[:uv_layout.size],uv_layout,reg,uv=True,mu_gev=1e5),reg)
    eft0,eft_layout=pack_general_scalar_state(projected,reg)
    eft0=np.concatenate((eft0,end_uv[uv_layout.size:]))
    eft_beta=core.check_payload(betas['payloads']['EFT1'],model=model,
             stage='EFT1',rep=eft_layout.representation,
             base_names=tuple(eft_layout.names[:-len(new)]),
             new_names=new,reused_names=reused)
    eftsol=core.integrate_segment(eft0,
        lambda t,y:core.rhs_at(t,y,layout=eft_layout,registry=reg,
                              betas=eft_beta,stage='EFT1'),
        1e5,1e3,n_points=3)
    endpoint=core.snapshot(eftsol.y[:,-1],eft_layout,reg,
                           stage='EFT1',mu_gev=1e3,new=new,reused=reused)
    masses=bc.heavy_flavour_at_muF(heavy_state)
    return endpoint,masses,{'mixed_boundary_quartic_fraction':q,
                            'mixed_boundary_mu2_fraction':m,
                            'mixed_uv_nfev':int(uvsol.nfev),
                            'mixed_eft1_nfev':int(eftsol.nfev),
                            'heavy_mass_basis_residual':float(masses[3]),
                            'heavy_fermion_masses_gev':np.asarray(masses[0],float).tolist()}


def independent_mixed(model: str,saved: Path,coefficient: Path,*,
                      hierarchy_factors: tuple[int,...]) -> tuple[list[dict],dict]:
    """Use existing mixed one-scale *compiled expression* for two mass-log inputs."""
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc
    source=coefficient.read_text(encoding='utf-8-sig')
    shifted,log_inventory=replace_common_mass_logs(source)
    compiled=bc.compile_source(coefficient)
    from scripts.RunGeneralT3NumericalC5 import compile_coefficient
    changed=compile_coefficient(shifted)
    endpoint,heavy,info=_mixed_running_endpoint(model,saved)
    results=[]
    for scale in hierarchy_factors:
        masses=np.asarray(heavy[0],float)*scale
        newheavy=(masses,heavy[1],heavy[2],heavy[3])
        original=bc.c5_at_endpoint(model,endpoint,newheavy,compiled)
        alternative=bc.c5_at_endpoint(model,endpoint,newheavy,changed)
        results.append({'factor':scale,'original':original,
                        'spurion':(alternative-original)/math.log(2),
                        'delta_2_multiplier':alternative-original})
    return results,{'one_scale_source_sha256':sha256_file(coefficient),
                    'one_scale_path':str(coefficient),
                    'log_inventory':log_inventory,
                    'yukawa_channel_basis':'existing BC c5_at_endpoint ordered source evaluator',
                    'mixed_quartic_mapping':'existing conditional lambda_M=-lambda_RGBeta',
                    **info}


def analyse_one_model(model: str, saved_fit: Path, saved_case: Path,
                      coefficient: Path, *,
                      wolframscript: str='wolframscript',
                      timeout_sec: int=300,
                      hierarchy_factors: tuple[int,...]=(1,2,4)) -> dict:
    require(model in MODELS, 'unknown T3 model '+model)
    require(hierarchy_factors and hierarchy_factors[0]==1 and
            len(set(hierarchy_factors))==len(hierarchy_factors) and
            all(type(k)==int and k>=1 for k in hierarchy_factors),
            'hierarchy factors must be distinct positive integers starting with 1')
    fit=read_json(saved_fit/'result.json')
    require(fit.get('model')==model,'mismatched fitted result')
    saved=saved_case_input(saved_case,model)
    if model in MIXED:
        pieces,provenance=independent_mixed(model,saved_fit,coefficient,
                                           hierarchy_factors=hierarchy_factors)
    else:
        pieces,provenance=independent_nonmixed(model,fit,coefficient,
            wolframscript=wolframscript,timeout_sec=timeout_sec,
            hierarchy_factors=hierarchy_factors)
    require(len(pieces)==len(hierarchy_factors), 'incomplete source log extraction')
    scan=[]
    for part in pieces:
        scale=part['factor']
        predicted=saved['eft_spurion']/scale
        comparison=compare_coefficients(part['spurion'],predicted)
        scan.append({'heavy_mass_scale_factor':scale,
                     'matchete_common_log_coefficient':encode_matrix(part['spurion']),
                     'eft_scalar_pair_common_log_prediction':encode_matrix(predicted),
                     'matchete_original_one_scale_C5':encode_matrix(part['original']),
                     'comparison':comparison})
    parity=original_one_scale_parity(pieces[0]['original'],saved['full'])
    return {'format':FORMAT,'model':model,'status':'Computed',
            'saved_fit_result_sha256':sha256_file(saved_fit/'result.json'),
            'saved_fitted_config_sha256':sha256_file(saved_fit/'fitted_config.json'),
            'saved_case_record_sha256':saved['case_record_sha256'],
            'mu_F_gev':saved['mu_f_gev'],'mu_S_gev':saved['mu_s_gev'],
            'common_log_L':saved['log_range'],
            'matchete_log_argument_shift':'Log[MF^2/Mscalar^2] -> Log[2*MF^2/Mscalar^2]',
            'coefficient_conversion':'C5_LLSS_EFT1/(2*log(muF/muS))',
            'saved_one_scale_reference_C5':encode_matrix(saved['full']),
            'saved_EFT1_leading_log_C5':encode_matrix(saved['eft_log']),
            'one_scale_reproduction_relative_residual':parity,
            'one_scale_anchor_reproduced_within_1e-4':parity<=1e-4,
            'hierarchy_scan':scan,'source_provenance':provenance,
            'no_subtracted_hard_matching_used':True,
            'independent_one_scale_mass_log_source_used':True,
            'spinor_normalization_and_quartic_basis_certified':False,
            'independent_sequential_hard_matching_calculated':False,
            'interpretation':('The common logarithmic coefficient in the saved '
              'full-theory Matchete expression is extracted independently by '
              'a log-argument spurion. Its matrix is compared to the separately '
              'generated tree-LLSS x HHSS scalar-pair EFT1 mixing coefficient. '
              'Masses and couplings are fixed for each log shift. A residual '
              'at finite MF/MS can be power-suppressed. Non-unity best-fit '
              'normalization or a nonzero orthogonal matrix residual may indicate '
              'additional operators or incorrect mapping. Not a full sequential '
              'hard matching certificate.')}
