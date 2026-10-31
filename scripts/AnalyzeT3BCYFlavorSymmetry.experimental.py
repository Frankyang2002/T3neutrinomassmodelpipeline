#!/usr/bin/env python3
"""Native B/C Y/LLHH PR/PL conjugacy and exact symmetric flavour factorization.

Read existing native neutral-component ledgers. This is NOT a GammaCC-to-Weyl
operator-normalization calculation or a physical Weinberg C5 coefficient.
No original files are modified; no Mathematica/Matchete calls are made.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path('output/audits/t3_bc_y_renormalization_ledger_v1')
STEMS = {'B_m1': 'B_-1', 'C_m1': 'C_-1'}
EXPECTED_QUARTICS = {'lambdaT3', 'lambdaT3ExtraQ020020I1',
                     'lambdaT3ExtraQ020200I1'}
EXPECTED_U = {'T3MEU11', 'T3MEU12', 'T3MEU21', 'T3MEU22'}


class KernelError(ValueError):
    pass


def require(ok, message):
    if not ok:
        raise KernelError(message)


def coupl(name, conjugated):
    return sp.Symbol(('bar_' if conjugated else '') + name, commutative=True)


def factor_key(record, flip=False):
    """Canonical exact monomial excluding chirality, including rational weight."""
    return (
        record['quartic'], bool(record['quartic_conjugated']) ^ flip,
        record['yukawa_i'], record['yukawa_j'],
        bool(record['yukawas_conjugated']) ^ flip,
        tuple(sorted((u['name'], bool(u['conjugated']) ^ flip,
                      u['exponent']) for u in record['mixing'])),
        record['sqrt3_exponent'],
    )


def coeff(record):
    a = record['coefficient_numerator']
    b = record['coefficient_denominator']
    require(type(a) is int and type(b) is int and b > 0,
            'Non-rational monomial coefficient')
    return Fraction(a, b)


def validate_record(record, stem, chirality):
    require(record['chirality'] == chirality, 'Chirality assignment changed')
    require(record['quartic'] in EXPECTED_QUARTICS, 'Unknown source quartic')
    require(record['yukawa_i'] in ('y1','y2') and
            record['yukawa_j'] in ('y1','y2'), 'Unknown Yukawa')
    require(record['yukawa_i'] <= record['yukawa_j'],
            'Unsupported reverse-order Yukawa term; inspect native source')
    require(record['yukawas_conjugated'] == (chirality == 'PL'),
            'Yukawa bar/chirality mismatch')
    require(record['sqrt3_exponent'] == (-1 if stem == 'C_m1' else 0),
            'Changed B/C group factor')
    require(type(record['quartic_conjugated']) is bool,
            'Nonboolean quartic conjugation')
    require(len(record['mixing']) >= 1, 'Missing mixing spurions')
    total=0
    for u in record['mixing']:
        require(u['name'] in EXPECTED_U and type(u['conjugated']) is bool
                and type(u['exponent']) is int and u['exponent'] > 0,
                'Unrecognized mixing spurion')
        total += u['exponent']
    require(total == 4, 'Unexpected quartic spurion degree')
    coeff(record)


def parse_ledger(payload, stem, kind):
    expected=32 if kind=='constant' else 50
    n_per=16 if kind=='constant' else 25
    require(payload.get('model') == STEMS[stem], 'Model label mismatch')
    require(payload.get('source_additive_terms') == expected,
            'Native source-term count changed')
    require(payload.get('one_common_factor') == 'hbar/Coupling[MF,{},0]',
            'Unexpected loop/mass normalization')
    require(payload.get('neutral_SU2_slice') ==
            'H_component_2__lepton_component_1__epsilon12_plus_1',
            'Unknown SU2 neutral slice')
    records=payload.get('monomials',[])
    require(len(records)==expected, 'Unexpected number of neutral monomials')
    require(payload['chirality_monomial_counts']=={'PR':n_per,'PL':n_per},
            'PR/PL monomial counts changed')
    part={'PR':defaultdict(Fraction),'PL':defaultdict(Fraction)}
    for rec in records:
        chir=rec['chirality']
        require(chir in part,'Unexpected Majorana projector')
        validate_record(rec,stem,chir)
        part[chir][factor_key(rec)]+=coeff(rec)
    require(all(x for x in part.values()), 'Empty chirality polynomial')
    # PL must be the conjugate of PR; never sum the two Hermitian conjugates.
    hc=defaultdict(Fraction)
    for rec in records:
        if rec['chirality']=='PR':
            hc[factor_key(rec,flip=True)]+=coeff(rec)
    require(dict(hc)==dict(part['PL']), 'PR/PL complex conjugacy failed')
    # Source left-to-right Yukawa slots are y1y1, y1y2, y2y2, with NO y2y1.
    groups=defaultdict(lambda: sp.Integer(0))
    for key,value in part['PR'].items():
        quartic,qbar,yi,yj,ybar,mixing,sqrt3exp=key
        require(not ybar, 'PR Yukawas unexpectedly conjugated')
        mon=sp.Rational(value.numerator,value.denominator)
        mon *= coupl(quartic,qbar)
        mon *= sp.sqrt(3)**sqrt3exp
        for nm,bar,exp in mixing:
            mon*=coupl(nm,bar)**exp
        groups[(yi,yj)] += mon
    for key in [('y1','y1'),('y1','y2'),('y2','y2')]:
        require(sp.expand(groups[key])!=0, 'Missing Yukawa bilinear '+str(key))
    # K is ORDERED, in native GammaCC/Proj[1] flavour ordering.
    # G=(K+K^T)/2 is a pure algebraic projection, NOT a normalization map.
    a11=sp.expand(groups[('y1','y1')])
    a12=sp.expand(groups[('y1','y2')])
    a22=sp.expand(groups[('y2','y2')])
    ordered=sp.Matrix([[a11,a12],[0,a22]])
    sym_average=(ordered+ordered.T)/2
    require(sym_average == sym_average.T, 'Symmetric projection failed')
    return {'ordered':ordered,'sym_average':sym_average,
            'coefficients':(a11,a12,a22), 'pr_monomials':n_per,
            'pl_monomials':n_per,'conjugacy':True,'records':records}


def check_between_models(b,c,kind):
    """Empirical exact parity found in the actual native B/C source.

    At fixed native mixing spurion and quartic coupling symbols, the C
    coefficient equals B / sqrt(3) for y_a*y_a and -B/sqrt(3) for y1*y2.
    This is NOT a relation between their numerical benchmark couplings.
    """
    BB=b['sym_average']
    CC=c['sym_average']
    require(sp.simplify(sp.sqrt(3)*CC[0,0]-BB[0,0])==0,
            kind+': 11 B/C relation failed')
    require(sp.simplify(sp.sqrt(3)*CC[1,1]-BB[1,1])==0,
            kind+': 22 B/C relation failed')
    require(sp.simplify(sp.sqrt(3)*CC[0,1]+BB[0,1])==0,
            kind+': 12 B/C sign relation failed')
    return True


def analytical_rank_self_test():
    x=sp.symbols('x1:4')
    y=sp.symbols('y1:4')
    a,b,c=sp.symbols('a b c')
    Y=sp.Matrix.hstack(sp.Matrix(x),sp.Matrix(y))
    G=sp.Matrix([[a,b],[b,c]])
    null=Y[:,0].cross(Y[:,1])
    require(all(sp.expand(z)==0 for z in Y.T*null),
            'Null-vector test failed')
    require(all(sp.expand(z)==0 for z in (Y*G*Y.T)*null),
            'Rank-two factorization failed')
    # Distinguish transposition from Hermitian conjugation for complex Yukawas.
    yy=sp.Matrix([[1,sp.I],[2,3],[4,5]])
    mm=yy*sp.Matrix([[2,1],[1,3]])*yy.T
    require(mm==mm.T and mm!=mm.conjugate().T,
            'Flavour transpose vs adjoint convention failed')
    return True


def self_test():
    def rec(ch,quartic,yi,yj,c=1):
        return {'chirality':ch,'quartic':quartic,
                'quartic_conjugated':ch=='PR','yukawa_i':yi,
                'yukawa_j':yj,'yukawas_conjugated':ch=='PL',
                'mixing':[{'name':'T3MEU12','conjugated':ch=='PL','exponent':4}],
                'sqrt3_exponent':0,
                'coefficient_numerator':c,'coefficient_denominator':1}
    core=[rec(ch,'lambdaT3',yi,yj,k) for ch in ('PR','PL')
          for yi,yj,k in [('y1','y1',2),('y1','y2',3),('y2','y2',5)]]
    # Synthetic shape test only, not a fake benchmark: exercise H.c. map.
    part={'PR':defaultdict(Fraction),'PL':defaultdict(Fraction)}
    for r in core:
        validate_record(r,'B_m1',r['chirality'])
        part[r['chirality']][factor_key(r)]+=coeff(r)
    require(set(part['PL'])=={factor_key(r,True) for r in core
                             if r['chirality']=='PR'},'Synthetic H.c. failed')
    tampered=[dict(x) for x in core]
    tampered[-1]['coefficient_numerator']+=1
    pld={factor_key(x):coeff(x) for x in tampered if x['chirality']=='PL'}
    prd={factor_key(x,True):coeff(x) for x in tampered if x['chirality']=='PR'}
    require(pld!=prd, 'Tampered complex-conjugate test was accepted')
    analytical_rank_self_test()
    print('PASS: exact PR/PL Hermitian-conjugacy and tamper rejection')
    print('PASS: transpose-not-adjoint symmetry, rank<=2 factorization')
    print('PASS: no GammaCC normalization or physical C5 certification')


def run(ledger_dir: Path, output: Path):
    report={'format':'t3-bc-y-ordered-flavour-factorization-v1',
            'input_directory':str(ledger_dir),'models':{},
            'common_native_factor':'hbar/MF; epsilon12=+1, ordered PR neutral component',
            'normalization':'sym_average=(ordered+transpose(ordered))/2 ONLY as algebraic projection',
            'yukawa_vector_convention':'y1_i,y2_i; no complex conjugation in PR',
            'physical_C5':None,'physical_neutrino_mass_matrix':None,
            'GammaCC_to_two_component_Weyl_normalization_certified':False,
            'full_EFT1_EFT2_operator_renormalization_certified':False,
            'finite_Y_coefficient_physical_certified':False,
            'scope':'neutral scalar-bubble pole and formal finite flavour structure only'}
    rows={}
    for stem in STEMS:
        rows[stem]={}
        item={'model':STEMS[stem],'ledgers':{}}
        for kind in ('constant','pole'):
            file=ledger_dir/f'{stem}_y_llhh_{kind}_neutral_flavour.json'
            require(file.is_file(),f'Missing native ledger {file}')
            buf=file.read_bytes()
            obj=json.loads(buf)
            parsed=parse_ledger(obj,stem,kind)
            rows[stem][kind]=parsed
            a11,a12,a22=parsed['coefficients']
            item['ledgers'][kind]={
                'source_file':str(file),'sha256':sha256(buf).hexdigest(),
                'PR_terms':parsed['pr_monomials'],
                'PL_terms':parsed['pl_monomials'],
                'PR_PL_exact_conjugacy':True,
                'ordered_A11':str(a11),'ordered_A12':str(a12),
                'ordered_A21':'0','ordered_A22':str(a22),
                'symmetric_G11':str(a11),
                'symmetric_G12':str(sp.expand(a12/2)),
                'symmetric_G22':str(a22),
                'Y_tensor_factorization':'Y @ G @ transpose(Y)',
                'rank_at_most_two_for_two_yukawa_vectors':True,
                'symmetric_by_flavour_transpose':True,
                'PR_and_PL_not_summed':True,
            }
        report['models'][stem]=item
    for kind in ('constant','pole'):
        check_between_models(rows['B_m1'][kind],rows['C_m1'][kind],kind)
    report['B_vs_C_native_symbolic_relation']={
        'checked_for':['constant','pole'],
        'diagonal_y1y1_y2y2':'C_native = B_native / sqrt(3)',
        'mixed_y1y2':'C_native = -B_native / sqrt(3)',
        'verified_exact':True,
        'same_quartic_and_mixing_symbols_only':True,
        'not_numeric_benchmark_relation':True,
    }
    analytical_rank_self_test()
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    for stem in STEMS:
        print(f'{STEMS[stem]}: finite 16 PR/16 PL; pole 25 PR/25 PL; '
              'exact H.c. and rank<=2 flavour factorization')
    print('PASS: C vs B diagonal +1/sqrt(3), mixed -1/sqrt(3) (native symbols)')
    print('Physical C5: NOT CERTIFIED')
    print('Saved: '+str(output))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--self-test',action='store_true')
    p.add_argument('--ledger-dir',type=Path,default=ROOT)
    p.add_argument('--output',type=Path)
    args=p.parse_args()
    try:
        if args.self_test:
            self_test();return 0
        if not args.output:p.error('--output is required')
        run(args.ledger_dir,args.output);return 0
    except (KernelError,KeyError,ValueError,TypeError,OSError,sp.SympifyError) as exc:
        print('FAIL (no physical C5): '+str(exc),file=sys.stderr)
        return 1

if __name__=='__main__':
    sys.exit(main())
