"""Strict projection of pole-cancelled Matchete LLHH onto ordered A5, and
optional frozen-benchmark numeric fixed-order evaluation.

This is the restricted leading-C5 scalar-threshold sector at a common mu=MF;
NOT EFT1 running, a full renormalised MSbar sequential calculation, or a
physical C5 prediction. Source operator and pole certificates are mandatory.

    python -m scripts.RunGeneralT3TwoThresholdNumericalC5 --all
    python -m scripts.RunGeneralT3TwoThresholdNumericalC5 --all --numeric
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any

import mpmath as mp
import numpy as np

from scripts.RunGeneralT3TwoThresholdPoleAssembly import (
    ROOT, OUTPUT_DIR as ASSEMBLY_DIR, MODELS, split_terms, render_terms,
    load_json, load_text, sha,
)
from scripts.RunGeneralT3NumericalC5 import compile_coefficient

FORMAT = 't3-restricted-two-threshold-fixed-order-weinberg-projection-v1'
STATUS = 'RestrictedFixedOrderCommonMuProjectionNotFullSequentialC5'
OUT_DIR = ROOT / 'output/general_t3_mass_eigen_matchete/two_threshold_numerical_projection'
REFERENCE_DIR = ROOT / 'output/general_t3_numerical_c5_frozen_scalar_scenario'
EPS = ('CG[Bar[eps[SU2L]], {Bar[Index[d$$1, SU2L[fund]]], '
       'Bar[Index[d$$2, SU2L[fund]]]}]'
       '*CG[Bar[eps[SU2L]], {Bar[Index[d$$3, SU2L[fund]]], '
       'Bar[Index[d$$4, SU2L[fund]]]}]')
# Two equivalent closed-index orientations occur in the C_-1 matched EFT.
# They differ by two antisymmetric epsilon swaps: net sign +1.
H_PATTERNS = (
    ('Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]',
     'Field[H, Scalar, {Index[d$$3, SU2L[fund]]}, {}]',
     'NCM[Transp[Field[l, Fermion, {Index[d$$4, SU2L[fund]], Index[d$$2, Flavor]}, {}]], '
     'DiracProduct[GammaCC, Proj[-1]], '
     'Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$1, Flavor]}, {}]]'),
    ('Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]',
     'Field[H, Scalar, {Index[d$$4, SU2L[fund]]}, {}]',
     'NCM[Transp[Field[l, Fermion, {Index[d$$3, SU2L[fund]], Index[d$$2, Flavor]}, {}]], '
     'DiracProduct[GammaCC, Proj[-1]], '
     'Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]]'),
)
PR_NCM = ('NCM[Bar[Field[l, Fermion, {Index[d$$3, SU2L[fund]], Index[d$$1, Flavor]}, {}]], '
          'DiracProduct[GammaCC, Proj[1]], '
          'Transp[Bar[Field[l, Fermion, {Index[d$$4, SU2L[fund]], Index[d$$2, Flavor]}, {}]]]]')
PR_H = ('Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]',
        'Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]')
PR_EPS = (
    ('CG[eps[SU2L], {Index[d$$1, SU2L[fund]], Index[d$$3, SU2L[fund]]}]',
     'CG[eps[SU2L], {Index[d$$2, SU2L[fund]], Index[d$$4, SU2L[fund]]}]'),
    ('CG[eps[SU2L], {Index[d$$3, SU2L[fund]], Index[d$$1, SU2L[fund]]}]',
     'CG[eps[SU2L], {Index[d$$4, SU2L[fund]], Index[d$$2, SU2L[fund]]}]'),
)


def check_conjugate_operator(value: str) -> str:
    if (value.count('Proj[1]') != 1 or value.count('GammaCC') != 1
        or value.count('Field[l, Fermion') != 2 or value.count('CG[') != 2
        or value.count('Field[H, Scalar') != 2):
        fail('Invalid PR conjugate spinor/Higgs/CG content')
    if value.count(PR_NCM) != 1 or any(value.count(x)!=1 for x in PR_H):
        fail('Unexpected PR spinor/Higgs contraction')
    choices=[i for i, eps in enumerate(PR_EPS) if all(value.count(x)==1 for x in eps)]
    if len(choices)!=1:fail('Unexpected PR SU(2) epsilon contraction')
    return 'ordinary' if choices[0] == 0 else 'double_epsilon_swap'


COUPLINGS = {'MF','MSX','MSY','T3MEU11','T3MEU12','T3MEU21','T3MEU22',
             'lambdaT3','lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1','y1','y2'}
ERROR = 'T3 restricted two-threshold C5 projection fail-closed: '


def fail(message: str) -> None:
    raise ValueError(ERROR + message)


def stem(model: str) -> str:
    if model not in MODELS:
        fail('Unsupported model ' + model)
    return model.replace('-', 'm').replace('+','p') + '_final'


def _remove_once(value: str, literal: str) -> str:
    if value.count(literal) != 1:
        fail('Expected exactly one operator factor ' + literal[:55])
    return value.replace(literal, '1', 1)


def project_holomorphic_term(value: str) -> tuple[str, str]:
    """Strip precisely the Matchete Weinberg operator factors.

    Mirrors RunMatching.wl StripWeinbergOperatorStructure: NCM, two H legs,
    two antisymmetric SU(2) tensors replaced with 1.  No guessed spinor or
    group factors and no guessing of external lepton indices is allowed.
    """
    if value.count('Proj[-1]') != 1 or value.count('GammaCC') != 1:
        fail('Expected a single PL GammaCC spin chain')
    if value.count('CG[') != 2 or value.count('Field[H, Scalar') != 2:
        fail('Unexpected SU2/Higgs tensor structure')
    if value.count('Field[l, Fermion') != 2 or value.count('NCM[') != 1:
        fail('Unexpected lepton spinor structure')
    if not EPS in value or value.count(EPS) != 1:
        fail('Unexpected epsilon contraction or ordering')
    recognized = [(i, triplet) for i, triplet in enumerate(H_PATTERNS)
                  if all(value.count(x) == 1 for x in triplet)]
    if len(recognized) != 1:
        fail('Unknown LLHH operator contraction / external-index orientation')
    orientation, operators = recognized[0]
    text = _remove_once(value, EPS)
    for x in operators:
        text = _remove_once(text, x)
    # Remaining index-bearing objects must be exactly two Yukawa spurions with
    # external flavours 1,2 and the *same* heavy-fermion flavour dummy.
    if any(q in text for q in ('GammaCC','NCM[','CG[','Field[','Proj[','SU2L','Lorentz')):
        fail('Operator factors survived stripping')
    for flavor in ('1', '2'):
        if text.count(f'Index[d$${flavor}, Flavor]') != 1:
            fail('Incorrect flavour-index multiplicity')
    if text.count('Index[d$$1, NFlavor]') != 2 or 'Index[d$$2, NFlavor]' in text:
        fail('Unexpected heavy fermion flavour contraction')
    used = set(re.findall(r'Coupling\[([A-Za-z0-9]+)', text))
    if not used <= COUPLINGS or {'MF'} - used:
        fail('Unknown or missing coefficient coupling: ' + repr(used - COUPLINGS))
    # Assert all remaining factors are coefficient-only, using the previously
    # validated AST compiler; do not evaluate untrusted Mathematica text.
    compile_coefficient(text)
    return text, ('ordinary' if orientation == 0 else 'double_epsilon_swap')


def project_expression(raw: str, model: str) -> dict[str, Any]:
    signs: list[tuple[int,str]]=[]
    orientations: dict[str,int]={}
    conjugate=0
    conjugate_orientations:dict[str,int]={}
    for sign, part in split_terms(raw):
        if (part.count('GammaCC') != 1 or part.count('hbar') != 1
            or part.count('Field[H, Scalar') != 2
            or part.count('Field[l, Fermion') != 2
            or 'Coupling[MF, {}, 0]' not in part
            or any(f'Field[{x},' in part for x in ('NewScalarX','NewScalarY','NewFermion'))
            or (part.count('Proj[-1]') + part.count('Proj[1]')) != 1):
            fail('Invalid fixed-muF LLHH operator')
        fixed_logs = (
            'Log[Coupling[MF, {}, 0]^2/Coupling[MSX, {}, 0]^2]',
            'Log[Coupling[MF, {}, 0]^2/Coupling[MSY, {}, 0]^2]',
        )
        if 'Log[' in part and (sum(part.count(x) for x in fixed_logs) != 1
                                   or part.count('Log[') != 1):
            fail('Unexpected fixed-muF matching logarithm')
        if 'ϵ' in part or 'μbar2' in part:
            fail('Remaining UV poles or common-scale logarithms in fixed-muF input')
        if 'Proj[1]' in part:
            conjugate += 1
            ori=check_conjugate_operator(part)
            conjugate_orientations[ori]=conjugate_orientations.get(ori,0)+1
            continue
        stripped, orientation = project_holomorphic_term(part)
        orientations[orientation] = orientations.get(orientation,0)+1
        signs.append((sign,stripped))
    if len(signs) != 116 or conjugate != 116:
        fail(f'Unexpected PL/conjugate counts {len(signs)}/{conjugate}')
    if model == 'B_-1' and (orientations != {'ordinary':116} or conjugate_orientations != {'ordinary':116}):
        fail('Unexpected B_-1 orientation')
    if model == 'C_-1' and (orientations != {'ordinary':34,'double_epsilon_swap':82} or conjugate_orientations != {'ordinary':34,'double_epsilon_swap':82}):
        fail('Unexpected C_-1 orientation')
    return {'ordered_coefficient':render_terms(signs),
            'PL_terms':len(signs),'conjugate_terms':conjugate,
            'orientations':orientations,'conjugate_orientations':conjugate_orientations}


def projection_report(model: str, assembly_dir: Path, outdir: Path) -> tuple[dict[str,Any],Path]:
    base = stem(model) + '_two_threshold_raw_llhh'
    audit_path=assembly_dir/(base+'_audit.json')
    audit=load_json(audit_path)
    # Verifies scope and ensures the exact 68/68 cancellation upstream, with
    # all 4 original source SHA-256 values retained as provenance.
    if (audit.get('model') != model or
        audit.get('status') != 'PoleCancellationVerifiedCommonScaleUnprojectedNotPhysicalC5' or
        audit.get('hard_and_scalar_pole_residues_opposite') is not True or
        audit.get('matched_canonical_pole_term_count') != 68 or
        audit.get('fixed_muF_nonpole_terms') != 232 or
        audit.get('excluded_scalar_self_quartic_terms') != 122 or
        audit.get('complete_scalar_potential_hermiticity_verified') is not False or
        audit.get('weinberg_flavour_su2_projection_performed') is not False or
        audit.get('physical_sequential_C5_computed') is not False or
        set(audit.get('input_provenance_sha256',{})) != set(
            [stem(model)+'_eft1_fermion_first_oneloop_raw_hard_llhh_sector.inputform.txt',
             stem(model)+'_eft1_fermion_first_oneloop_eft1_one_loop_correction.inputform.txt',
             stem(model)+'_eft1_scalar_threshold_probe_c5_relevant_raw_scalar_llhh.inputform.txt',
             stem(model)+'_eft1_scalar_threshold_probe_c5_relevant_scalar_one_loop_correction.inputform.txt'])):
        fail('Missing or inconsistent upstream pole certificate')
    expression_path=assembly_dir/(base+'_fixed_muF_unprojected.inputform.txt')
    raw=load_text(expression_path)
    result=project_expression(raw,model)
    outdir.mkdir(parents=True,exist_ok=True)
    dest=outdir/(stem(model)+'_restricted_fixed_muF_ordered_A5.inputform.txt')
    dest.write_text(result['ordered_coefficient']+'\n',encoding='utf-8')
    cert={
        'format':FORMAT,'status':STATUS,'model':model,
        'fixed_muF_operator_source_sha256':sha(expression_path),
        'source_pole_certificate_sha256':sha(audit_path),
        'original_4_matchete_source_sha256':audit['input_provenance_sha256'],
        'upstream_pole_cancellation_verified':True,
        'common_matching_scale':'mu_bar=MF',
        'PL_holomorphic_term_count':result['PL_terms'],
        'PR_hermitian_conjugate_term_count':result['conjugate_terms'],
        'PL_operator_SU2_orientations':result['orientations'],
        'PR_operator_SU2_orientations':result['conjugate_orientations'],
        'projected_ordered_A5_sha256':sha(dest),
        'projected_ordered_A5_file':dest.name,
        'normalization_convention':'RunMatching.wl StripWeinbergOperatorStructure; C5=A5+A5^T; mnu=-v^2*C5/2',
        'hbar':'1/(16*pi^2)',
        'physical_C5_computed':False,
        'full_EFT1_running_performed':False,
        'full_scalar_EFT_Hermiticity_certified':False,
        'restricted_four_scalar_quartics_excluded':122,
        'full_sequential_MSbar_threshold_normalization_certified':False,
    }
    out=outdir/(stem(model)+'_restricted_fixed_muF_projection_audit.json')
    out.write_text(json.dumps(cert,indent=2)+'\n',encoding='utf-8')
    return cert,dest


def _number(obj: Any) -> mp.mpc:
    if not isinstance(obj,dict) or set(obj) != {'real','imag'}:
        fail('Complex value missing real/imag parts')
    try: z=mp.mpc(str(obj['real']),str(obj['imag']))
    except (TypeError,ValueError) as exc:raise ValueError(ERROR+'Invalid complex input') from exc
    if not mp.isfinite(z.real) or not mp.isfinite(z.imag):fail('Nonfinite input')
    return z


def numerical_evaluation(model: str, expression: str, reference_path: Path,
                         outdir: Path, source_report: dict[str,Any]) -> dict[str,Any]:
    if not reference_path.is_file():
        fail('Missing frozen-aligned one-scale reference result: ' + str(reference_path) +
             '; first run python -m scripts.RunGeneralT3FrozenScalarC5 --all')
    old=load_json(reference_path)
    if (old.get('model') != model or old.get('status') !=
        'FrozenScalarDiagonalDirectOneScaleC5NotPhysicalEFT1' or
        old.get('frozen_scalar_aligned_one_scale_scenario') is not True or
        old.get('mass_basis_compatible_with_frozen_reference') is not True or
        old.get('physical_c5_certified') is not False or
        old.get('sequential_EFT1_performed') is not False or
        old.get('source_mass_parameters_used_as_numeric_scalar_input') is not False or
        not isinstance(old.get('fermion_masses_gev'),list) or
        len(old.get('scalar_masses_gev',[])) != 2):
        fail('Incompatible existing frozen-aligned benchmark')
    config=old['benchmark_config']
    if config.get('format') != 't3-frozen-scalar-diagonal-direct-C5-input-v1':
        fail('Wrong frozen reference input format')
    v=mp.mpf(str(config['vev_gev']))
    if not mp.isfinite(v) or v<=0:fail('Non-positive vev')
    U=old['effective_scalar_basis_U']
    re_U, im_U = U['real'],U['imag']
    if any(len(row)!=2 for row in re_U) or len(re_U)!=2 or len(im_U)!=2:
        fail('Wrong scalar mixing matrix size')
    Umat=np.array(re_U,dtype=float)+1j*np.array(im_U,dtype=float)
    if max(abs(x) for x in (Umat.conj().T@Umat-np.eye(2)).flatten())>1e-10:
        fail('Non-unitary mass-basis U')
    env={f'T3MEU{i+1}{j+1}':mp.mpc(str(re_U[i][j]),str(im_U[i][j]))
         for i in range(2) for j in range(2)}
    for k in ('lambdaT3','lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1'):
        env[k]=_number(config['quartics'][k])
    for i,name in enumerate(('MSX','MSY')):
        env[name]=mp.mpf(str(old['scalar_masses_gev'][i]))
        if env[name]<=0:fail('Invalid scalar mass')
    fermion=[mp.mpf(str(m)) for m in old['fermion_masses_gev']]
    if any(m<=0 for m in fermion):fail('Invalid fermion mass')
    y1=[[_number(v) for v in row] for row in config['y1']]
    y2=[[_number(v) for v in row] for row in config['y2']]
    if len(y1)!=3 or len(y2)!=3 or any(len(row)!=len(fermion) for row in y1+y2):
        fail('Three-flavour Yukawa dimensions must match fermion masses')
    mp.mp.dps=55
    env.update({'hbar':1/(16*mp.pi**2),'conj':mp.conj,'sqrt':mp.sqrt,'log':mp.log})
    compiled=compile_coefficient(expression)
    ordered=[[mp.mpc(0) for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k,m in enumerate(fermion):
                variables=dict(env,MF=m,y1a=y1[i][k],y1b=y1[j][k],
                               y2a=y2[i][k],y2b=y2[j][k])
                try: z=eval(compiled,{'__builtins__':{}},variables)
                except Exception as exc:
                    fail(f'Coefficient evaluation failed at ({i},{j},{k}): {exc}')
                if not mp.isfinite(z.real) or not mp.isfinite(z.imag):
                    fail('Nonfinite projected ordered A5')
                ordered[i][j]+=z
    c5=[[ordered[i][j]+ordered[j][i] for j in range(3)] for i in range(3)]
    mnu=[[-v*v*c5[i][j]/2 for j in range(3)] for i in range(3)]
    sv=np.linalg.svd(np.array([[complex(x) for x in r] for r in mnu],dtype=complex),compute_uv=False)*1e9
    old_sv=old['takagi_singular_values_eV_descending']
    def complex_json(z:mp.mpc) -> dict[str,float]:
        return {'real':float(z.real),'imag':float(z.imag)}
    result={
       'format':FORMAT,'status':'RestrictedFixedOrderCommonMuFNumericalC5Exploratory',
       'model':model,'source_projection_audit_sha256':sha(outdir/(stem(model)+'_restricted_fixed_muF_projection_audit.json')),
       'source_frozen_reference_one_scale_result_sha256':sha(reference_path),
       'source_ordered_A5_sha256':source_report['projected_ordered_A5_sha256'],
       'scope':'restricted C5-relevant scalar threshold; pole-cancelled fixed order at mu_bar=MF, no EFT1 RGE',
       'renormalized_physical_C5_certified':False,
       'full_EFT1_running_performed':False,
       'full_scalar_EFT_hermiticity_certified':False,
       'same_frozen_scalar_and_flavor_scenario_as_one_scale':True,
       'comparison_is_not_independent_neutrino_mass_prediction':True,
       'weinberg_conversion':'C5=A5+A5^T; mnu=-v^2*C5/2; hbar=1/(16*pi^2)',
       'fermion_masses_gev':[float(m) for m in fermion],
       'scalar_masses_gev':[float(env['MSX']),float(env['MSY'])],
       'ordered_A5_gev_inv':[[complex_json(x) for x in r] for r in ordered],
       'symmetric_C5_gev_inv':[[complex_json(x) for x in r] for r in c5],
       'mnu_gev':[[complex_json(x) for x in r] for r in mnu],
       'takagi_singular_values_eV_descending':sv.tolist(),
       'largest_takagi_mass_eV':float(sv[0]),
       'one_scale_largest_takagi_mass_eV':float(old_sv[0]),
       'largest_takagi_ratio_to_one_scale':float(sv[0]/old_sv[0]),
       'note':'Numerical projection of a restricted, one-loop-leading fixed-order operator, NOT a fully renormalized sequential EFT1 result.'
    }
    path=outdir/(stem(model)+'_restricted_fixed_muF_numerical_c5.json')
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    selector=parser.add_mutually_exclusive_group(required=True)
    selector.add_argument('--model',choices=MODELS)
    selector.add_argument('--all',action='store_true')
    parser.add_argument('--numeric',action='store_true',help='Evaluate frozen-aligned benchmark using existing one-scale output')
    parser.add_argument('--assembly-dir',type=Path,default=ASSEMBLY_DIR)
    parser.add_argument('--reference-dir',type=Path,default=REFERENCE_DIR)
    parser.add_argument('--output-dir',type=Path,default=OUT_DIR)
    args=parser.parse_args()
    summary=[]
    for model in MODELS if args.all else (args.model,):
        report,coef_path=projection_report(model,args.assembly_dir,args.output_dir)
        item={'model':model,'PL_terms':report['PL_holomorphic_term_count'],
              'PR_terms':report['PR_hermitian_conjugate_term_count'],
              'pole_cancellation_verified':True,'numerical_evaluation_performed':False}
        if args.numeric:
            result=numerical_evaluation(model,coef_path.read_text(encoding='utf-8'),
                          args.reference_dir/(stem(model)+'_numerical_c5.json'),args.output_dir,report)
            item.update({'numerical_evaluation_performed':True,
                         'largest_takagi_mass_eV':result['largest_takagi_mass_eV'],
                         'one_scale_largest_takagi_mass_eV':result['one_scale_largest_takagi_mass_eV'],
                         'ratio_to_one_scale':result['largest_takagi_ratio_to_one_scale']})
        summary.append(item)
        print(f'{model}: 116 PL+116 conjugate LLHH projected; '
              + (f'fixed-order leading mass={item["largest_takagi_mass_eV"]:.8g} eV, '
                 f'one-scale={item["one_scale_largest_takagi_mass_eV"]:.8g} eV; '
                 if args.numeric else '')
              + 'full physical sequential C5 NOT certified',flush=True)
    path=args.output_dir/'summary_final.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps({'status':STATUS,'models':summary,
       'physical_sequential_C5_computed':False,
       'full_EFT1_running_performed':False},indent=2)+'\n',encoding='utf-8')
    print('Summary: '+str(path))

if __name__ == '__main__':
    main()
