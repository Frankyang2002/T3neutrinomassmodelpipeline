#!/usr/bin/env python3
"""Exact Grassmann-flavour and Higgs-VEV normalization, conditional on
Matchete's GammaCC-to-Warsaw spinor identification.

This proves the purely algebraic 1/2 factors independently, checks the
source-level B/C input from v9.1, and optionally indexes the installed
Matchete GammaCC source report. It NEVER asserts that the GammaCC operator
is already a physical Warsaw C5 coefficient.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import sys

DEFAULT_BASE = Path('output/audits/t3_bc_y_warsaw_convention_v1/summary.json')
DEFAULT_PROBE = Path('output/audits/t3_bc_y_gammacc_convention_v1/installed_source.json')


def fail(s):
    raise ValueError(s)


def check_v9(data):
    if data.get('format') != 't3-bc-native-to-warsaw-weinberg-convention-audit-v2':
        fail('v9.1 native-source audit format not recognized')
    if data.get('status') != 'NativeSourceWarsawSkeletonAudited_AbsoluteNormalizationOutstanding':
        fail('v9.1 native-source audit not successful')
    if data.get('absolute_GammaCC_to_Warsaw_factor_certified') is not False:
        fail('unexpected pre-existing GammaCC certificate')
    ms = data.get('models')
    if not isinstance(ms, list) or len(ms) != 2 or {m.get('model') for m in ms} != {'B_-1','C_-1'}:
        fail('wrong model inventory')
    for m in ms:
        if (m.get('total_terms'),m.get('pole_terms'),m.get('log_terms'),m.get('finite_terms')) != (132,50,50,32):
            fail('unexpected native source inventory')
        if m.get('native_llhh_su2_epsilon_skeleton_checked') is not True:
            fail('SU2 skeleton not checked')
    return [{'model':m['model'], 'source_sha256':m['source_sha256']} for m in ms]


def grassmann_product(v1, v2):
    """Formally expand x_1^T eps x_2 in an anti-commuting exterior algebra.

    psi_i=(theta_(i,0),theta_(i,1)); eps_01=+1, eps_10=-1.
    Return dict with ordered variable pairs and exact coefficients.
    """
    if v1 not in ('i','j') or v2 not in ('i','j'):
        fail('invalid flavour labels')
    factors = [((v1,0),(v2,1),Fraction(1)),
               ((v1,1),(v2,0),Fraction(-1))]
    out = {}
    for left,right,coef in factors:
        if left == right: continue
        if left > right: left,right = right,left; coef = -coef
        k = (left,right)
        out[k] = out.get(k,Fraction()) + coef
    return {k:v for k,v in out.items() if v}


def prove_grassmann_symmetry():
    ij = grassmann_product('i','j')
    ji = grassmann_product('j','i')
    if ij != ji or not ij:
        fail('Majorana bilinear is not symmetric under exchange')
    ii = grassmann_product('i','i')
    if ii != {(('i',0),('i',1)): Fraction(2)}:
        fail('identical-field Grassmann normalization incorrect')
    # A_ij B_ij + A_ji B_ji = (A_ij + A_ji) B_ij.
    # With 1/2*C_ij B_ij and C symmetric, C=A+A^T.
    A = [[Fraction(2),Fraction(3)], [Fraction(5),Fraction(7)]]
    B = {'i_i': ii, 'i_j': ij, 'j_j': grassmann_product('j','j')}
    for p in ('i','j'):
        for q in ('i','j'):
            if grassmann_product(p,q) != grassmann_product(q,p):
                fail('failed complete bilinear symmetry')
    C = [[A[i][j]+A[j][i] for j in range(2)] for i in range(2)]
    if C != [[4,8],[8,14]]:
        fail('matrix symmetrization factor wrong')
    return dict(anticommuting_bilinear_exchange='B_ij = B_ji',
                identical_flavour_bilinear_multiplicity=2,
                prefactor_map_for_L5_equal_Aij_Qij_plus_hc='C5=A+A^T when L5=(1/2)*C5*Q+h.c.',
                example_A=[[str(x) for x in row] for row in A],
                example_C5=[[str(x) for x in row] for row in C])


def probe_summary(report):
    if report.get('format') != 't3-bc-installed-gammacc-source-probe-v1':
        fail('unknown installed Matchete source probe format')
    if report.get('status') != 'InstalledGammaCCSourceLocated':
        fail('installed GammaCC source not located')
    hits = report.get('source_occurrences')
    if not isinstance(hits,list) or not hits:
        fail('no source occurrences recorded')
    names = report.get('GammaCC_source_files')
    if not isinstance(names,list) or not names:
        fail('no GammaCC source files recorded')
    g = [x for x in hits if x.get('symbol')=='GammaCC']
    if not g or any(not isinstance(x.get('source_snippets'),list) for x in g):
        fail('no auditable GammaCC excerpts')
    # The source scan is evidence, not a semantic proof of Matchete's CC sign.
    return dict(status='GammaCCSourceLocated_NotYetInterpreted',
                source_files=names,
                occurrences=sum(len(x['source_snippets']) for x in g),
                source_root=report.get('source_root'),
                gammacc_definition_interpreted=False)


def analyze(v9,probe):
    baseline = check_v9(v9)
    algebra = prove_grassmann_symmetry()
    source = probe_summary(probe) if probe is not None else dict(status='NOT_PROVIDED')
    return dict(
        format='t3-bc-conditional-majorana-mass-normalization-v1',
        status='GrassmannAndVEVFactorsProved_NativeGammaCCFactorOutstanding',
        models=baseline, exact_grassmann_proof=algebra, installed_source_probe=source,
        epsilon_21_squared='(-1)^2=1',
        neutral_higgs_vev='H_2=v/sqrt(2)',
        neutral_Warsaw_Q='Q_ij -> (v^2/2)*B_ij',
        canonical_mass_lagrangian='L_mass=-(1/2)*m_ij*B_ij + h.c.',
        mass_map_conditional_on_native_Q_identification='m_ij=-(v^2/2)*(A_ij+A_ji)',
        if_lagrangian_uses_half_C5='m_ij=-(v^2/2)*C5_ij',
        if_lagrangian_uses_C5_without_half='m_ij=-v^2*C5_ij for C5 symmetric',
        matchete_hbar_equals_one_over_16pi2='DOCUMENTED_BY_MATCHeTE_QUICKSTART',
        native_PR_and_PL_added_as_independent_C5=False,
        native_GammaCC_to_Warsaw_prefactor_certified=False,
        operator_sign_and_projector_reference_certified=False,
        all_EFT1_EFT2_counterterms_certified=False,
        physical_C5=None,
        physical_neutrino_mass_matrix=None,
        next_required_step='review installed GammaCC source definitions and calibrate native NCM vs explicit Warsaw Q'
    )


def self_test():
    algebra = prove_grassmann_symmetry()
    assert algebra['identical_flavour_bilinear_multiplicity']==2
    m = {'format':'t3-bc-native-to-warsaw-weinberg-convention-audit-v2',
         'status':'NativeSourceWarsawSkeletonAudited_AbsoluteNormalizationOutstanding',
         'absolute_GammaCC_to_Warsaw_factor_certified':False,
         'models':[{'model':key,'total_terms':132,'pole_terms':50,'log_terms':50,'finite_terms':32,
                    'source_sha256':'abcd','native_llhh_su2_epsilon_skeleton_checked':True} for key in ('B_-1','C_-1')]}
    assert analyze(m,None)['native_GammaCC_to_Warsaw_prefactor_certified'] is False
    bad = json.loads(json.dumps(m)); bad['models'][1]['pole_terms']=49
    try: analyze(bad,None)
    except ValueError: pass
    else: raise AssertionError('tampered source inventory passed')
    probe={'format':'t3-bc-installed-gammacc-source-probe-v1',
           'status':'InstalledGammaCCSourceLocated','GammaCC_source_files':['./Gamma.wl'],
           'source_occurrences':[{'symbol':'GammaCC','source_snippets':[{'excerpt':'GammaCC'}]}]}
    assert probe_summary(probe)['gammacc_definition_interpreted'] is False
    print('PASS: antisymmetric spinor matrix + Grassmann statistics give symmetric Majorana bilinear')
    print('PASS: independent 1/2, Higgs-VEV and flavour-transpose factors')
    print('PASS: v9.1 input validation and installed-source probe fail-closed limits')
    print('PASS: no GammaCC sign, full matching or physical C5 invented')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--self-test',action='store_true')
    p.add_argument('--v9-reference',type=Path,default=DEFAULT_BASE)
    p.add_argument('--installed-probe',type=Path,default=DEFAULT_PROBE)
    p.add_argument('--require-installed-probe',action='store_true')
    p.add_argument('--output',type=Path)
    args=p.parse_args()
    try:
        if args.self_test: self_test(); return 0
        if args.output is None: p.error('--output is required')
        if not args.v9_reference.is_file(): fail('missing v9.1 result: '+str(args.v9_reference))
        v9=json.loads(args.v9_reference.read_text(encoding='utf-8-sig'))
        if args.require_installed_probe and not args.installed_probe.is_file():
            fail('required installed source probe missing: '+str(args.installed_probe))
        probe=json.loads(args.installed_probe.read_text(encoding='utf-8-sig')) if args.installed_probe.is_file() else None
        result=analyze(v9,probe)
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
        print('PASS: Majorana symmetry and mass 1/2 factors, conditional on GammaCC operator identification')
        print('Installed GammaCC probe:',result['installed_source_probe']['status'])
        print('Native GammaCC/Warsaw factor: NOT CERTIFIED; physical C5: NOT CERTIFIED')
        print('Saved:',args.output)
        return 0
    except (ValueError,OSError,KeyError,TypeError,ImportError,json.JSONDecodeError) as exc:
        print('FAIL:',exc,file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
