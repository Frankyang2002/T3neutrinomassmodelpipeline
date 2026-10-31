#!/usr/bin/env python3
"""Source-anchored F-hard one-loop LLHH boundary candidate for T3 B(-1)/C(-1).

The original F-first one-loop Matchete LLHH output contains 68 pole, 68
MF-log and 68 constant additive terms per model. Parse their native SU(2)
component coefficients with the existing, strict project parser; use the
v11.2-confirmed PL=Warsaw orientation, symmetrize flavour, and KEEP the F
boundary separate from the v12 conditional Y scalar-bubble contribution.

No new renormalization calculation: the F-hard constant sector is FORMAL,
not yet a renormalized finite matching coefficient. UV vs IR pole assignment,
independent F/EFT1 counterterms, inherited running and EOM shifts are unknown.
Do not sum the F and Y boundary coefficients or call the result physical C5.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import importlib
import importlib.util
import json
from pathlib import Path
import sys

STEMS = {'B_m1':'B_-1', 'C_m1':'C_-1'}
LOG_F='Log[μbar2/Coupling[MF, {}, 0]^2]'

class BoundaryError(ValueError):
    pass


def require(ok, message):
    if not ok: raise BoundaryError(message)


def load(path, fmt=None):
    require(path.is_file(), 'Missing local input: '+str(path))
    o=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(o,dict), 'Expected JSON object: '+str(path))
    if fmt is not None: require(o.get('format')==fmt, 'Unexpected evidence format: '+str(path))
    return o


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def configure_project_imports(root: Path) -> None:
    """Resolve the native Numerical namespace from the requested local repo.

    Python starts with scripts/ on sys.path when invoked by filename. Unlike
    `python -m`, that does not guarantee the sibling Numerical/ directory is
    importable. Use only the explicitly selected local repository root; fail
    if its audited parser sources are absent instead of importing a different
    package with the same name from the environment.
    """
    required = (
        root / 'Numerical' / 'experimental' / 'T3BCPoleLogConsistency.py',
        root / 'Numerical' / 'experimental' / 'T3BCNativeYLLHHFlavorKernel.py',
    )
    missing = [str(path) for path in required if not path.is_file()]
    require(not missing, 'Required local Numerical source missing: ' + ', '.join(missing))
    root_str = str(root.resolve())
    if root_str in sys.path:
        sys.path.remove(root_str)
    sys.path.insert(0, root_str)


def check_project_imports(root: Path) -> None:
    """Exercise the imports used by the actual B/C run, not only synthetic algebra."""
    configure_project_imports(root)
    for module_name in (
        'Numerical.experimental.T3BCPoleLogConsistency',
        'Numerical.experimental.T3BCNativeYLLHHFlavorKernel',
    ):
        module = importlib.import_module(module_name)
        expected = (root / Path(*module_name.split('.'))).with_suffix('.py').resolve()
        actual = Path(module.__file__).resolve()
        require(actual == expected,
                'Unexpected Numerical source imported: ' + str(actual) +
                '; expected ' + str(expected))


def split_sectors(text):
    # Reuse the existing fail-closed project InputForm parser.
    from Numerical.experimental.T3BCPoleLogConsistency import terms
    ts=terms(text)
    require(len(ts)==204, f'F-hard LLHH inventory changed: {len(ts)} terms')
    out={'pole':[], 'log':[], 'constant':[]}
    for term in ts:
        eps=term.count('ϵ')
        lg=term.count('Log[')
        require((eps,lg) in {(1,0),(0,1),(0,0)}, 'Unexpected F-hard pole/log structure')
        if eps:
            require(term.count('ϵ*')==1, 'Unsupported F-hard epsilon denominator')
            out['pole'].append(term.replace('ϵ*','',1))
        elif lg:
            require(term.count(LOG_F)==1, 'Unexpected F-hard mass logarithm')
            out['log'].append(term.replace(LOG_F,'1',1))
        else:
            out['constant'].append(term)
    require({k:len(v) for k,v in out.items()}=={'pole':68,'log':68,'constant':68},
            'Unexpected F-hard term count by sector')
    return out


def exact_pole_log_check(sectors):
    from Numerical.experimental.T3BCPoleLogConsistency import polynomial
    table={}
    p=polynomial(tuple(sectors['pole']),table=table)
    l=polynomial(tuple(sectors['log']),table=table)
    require(p==l, 'Source F-hard pole/log coefficient mismatch')
    return len(p)


def make_ledger(terms,model,expected=68):
    # This parser already audits original native fields, epsilons, ΓCC,
    # chirality, quartic, MF and the C(-1) 1/sqrt(3) group factor.
    from Numerical.experimental.T3BCNativeYLLHHFlavorKernel import _record_from_term
    require(len(terms)==expected, 'Incorrect source sector size')
    coeffs=defaultdict(Fraction)
    records={}
    counts=Counter()
    for term in terms:
        rec,c=_record_from_term(term,model)
        counts[rec['chirality']]+=1
        key=(rec['chirality'],rec['quartic'],rec['quartic_conjugated'],
             rec['yukawa_i'],rec['yukawa_j'],rec['yukawas_conjugated'],
             tuple(rec['mixing']),rec['sqrt3_exponent'])
        coeffs[key]+=c
        records[key]=rec
    require(counts=={'PR':34,'PL':34}, 'Incorrect PR/PL F-hard inventory')
    ledger=[]
    for key in sorted(records):
        n=coeffs[key]
        if not n: continue
        rec=records[key]
        ledger.append({
            'chirality':rec['chirality'], 'quartic':rec['quartic'],
            'quartic_conjugated':rec['quartic_conjugated'],
            'yukawa_i':rec['yukawa_i'], 'yukawa_j':rec['yukawa_j'],
            'yukawas_conjugated':rec['yukawas_conjugated'],
            'mixing':[{'name':name,'conjugated':bar,'exponent':exp}
                      for name,bar,exp in rec['mixing']],
            'sqrt3_exponent':rec['sqrt3_exponent'],
            'coefficient_numerator':n.numerator,
            'coefficient_denominator':n.denominator,
        })
    require(ledger, 'All source monomials vanished')
    return {'model':model,'monomials':ledger},dict(counts)


def import_v12(root):
    path=root/'scripts'/'BuildT3BCYFiniteBubbleInterface.experimental.py'
    require(path.is_file(), 'Requires existing v12 interface at '+str(path))
    spec=importlib.util.spec_from_file_location('t3bc_v12_f_boundary',path)
    require(spec is not None and spec.loader is not None, 'Failed to import v12')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def construct_polynomials(ledger,v12):
    import sympy as sp
    pl=v12.polynomial_coefficients(ledger,'PL')
    pr=v12.polynomial_coefficients(ledger,'PR')
    require(all(sp.expand(v12.polynomial_conjugate(x)-y)==0 for x,y in zip(pr,pl)),
            'F hard PR and PL are not exact H.c. coefficient polynomials')
    a11,a12,a22=pl
    A=sp.Matrix([[a11,a12],[0,a22]])
    C5=A+A.T  # convention of installed Warsaw model, fixed by v11.2
    require(sp.simplify(C5-C5.T)==sp.zeros(2), 'Non-symmetric normalized tensor')
    return {
        'PL_ordered_Warsaw_A_2x2_without_hbar_MF':[[str(A[i,j]) for j in range(2)] for i in range(2)],
        'conventional_C5_2x2_without_hbar_MF':[[str(C5[i,j]) for j in range(2)] for i in range(2)],
        'independent_PR_PL_conjugacy':True,
        'common_factor':'1/(16*pi**2*MF)',
        'matrix_flavour_orientation':'L=[conj(y1),conj(y2)]; A=(hbar/MF)*L@A2x2@L.T; C5=A+A.T',
    }


def analyze_one(root,audits,stem,v12):
    model=STEMS[stem]
    path=(root/'output'/'general_t3_mass_eigen_matchete'/'eft1_fermion_first_oneloop'/
          f'{stem}_final_eft1_fermion_first_oneloop_raw_hard_llhh_sector.inputform.txt')
    require(path.is_file(), 'Original F-hard LLHH source not available: '+str(path))
    text=path.read_text(encoding='utf-8-sig')
    sectors=split_sectors(text)
    opaque_monomials=exact_pole_log_check(sectors)
    ledgers={}
    counts={}
    polynomials={}
    for name in ('pole','log','constant'):
        ledgers[name],counts[name]=make_ledger(sectors[name],model)
        polynomials[name]=construct_polynomials(ledgers[name],v12)
    # Pole and log source algebra was already checked. Cross-check also
    # independently in the neutral Majorana component reduction.
    require(polynomials['pole']==polynomials['log'],
            'F-hard neutral SU2 pole and MF-log coefficients disagree')
    # No equality between the finite term and pole is assumed: report it.
    same=polynomials['constant']['conventional_C5_2x2_without_hbar_MF']==\
         polynomials['pole']['conventional_C5_2x2_without_hbar_MF']
    return {
      'model':model,
      'source_path':str(path.relative_to(root)),
      'source_sha256':sha256(path),
      'raw_additive_terms':204,
      'raw_sector_counts':{'pole':68,'log':68,'constant':68},
      'PR_PL_counts_each_sector':counts,
      'pole_equals_MF_log_coefficient':True,
      'opaque_pole_polynomial_terms':opaque_monomials,
      'constant_equals_pole_coefficient_observed_not_assumed':same,
      'sectors':polynomials,
      'constant_is_FORMAL_not_independently_renormalized':True,
      'UV_IR_residues_assigned':False,
      'F_to_EFT1_independent_counterterms_calculated':False,
      'F_hard_finite_physical_matching_certified':False,
      'physical_C5':None,
    }


def self_test():
    import sympy as sp
    # Algebraic normalized-Warsaw regression, no fabricated physical CT.
    root=Path(__file__).resolve().parents[1]
    v12=import_v12(root)
    def mon(ch, yi,yj,scale):
        return {'chirality':ch,'quartic':'lambdaT3',
                'quartic_conjugated':ch=='PR',
                'yukawa_i':yi,'yukawa_j':yj,
                'yukawas_conjugated':ch=='PL',
                'mixing':[{'name':'T3MEU11','conjugated':ch=='PL','exponent':4}],
                'sqrt3_exponent':0,
                'coefficient_numerator':scale,'coefficient_denominator':1}
    l={'monomials':[mon(ch,yi,yj,n) for ch in ('PL','PR')
                    for yi,yj,n in [('y1','y1',2),('y1','y2',3),('y2','y2',5)]]}
    c=construct_polynomials(l,v12)
    require(c['conventional_C5_2x2_without_hbar_MF'][0][1] is not None,
            'Synthetic mixed Yukawa missing')
    require(sp.expand(sp.sympify(c['conventional_C5_2x2_without_hbar_MF'][0][1]) -
                         3*sp.Symbol('lambdaT3')*sp.Symbol('bar_T3MEU11')**4)==0,
            'Mixed Yukawa factor wrongly doubled')
    tampered={'monomials':[dict(x) for x in l['monomials']]}
    tampered['monomials'][-1]['coefficient_numerator']=6
    try:construct_polynomials(tampered,v12)
    except BoundaryError:pass
    else:raise BoundaryError('Tampered H.c. source accepted')
    print('PASS: independent F boundary has ordered Warsaw A; mixed Yukawa not doubled')
    # This is the same import environment used during the B/C production run.
    check_project_imports(root)
    print('PASS: complex-conjugation and tamper rejection; no physical matching claimed')
    print('PASS: local Numerical parser imports from repository root')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--root',type=Path,default=Path.cwd())
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    try:
        if args.self_test:self_test();return 0
        require(args.output is not None,'--output required')
        root=args.root.resolve();audits=root/'output'/'audits'
        check_project_imports(root)
        v12path=audits/'t3_bc_y_conditional_finite_interface_v1'/'summary.json'
        v11path=audits/'t3_bc_y_warsaw_direct_identity_v11'/'summary_live.json'
        v12summary=load(v12path,'t3bc-y-restricted-conditional-finite-warsaw-interface-v1')
        v11summary=load(v11path,'t3-bc-warsaw-spinor-orientation-v11')
        require(v12summary.get('full_physical_matching_certified') is False,
                'v12 status unexpectedly changed')
        require(v11summary.get('live_Matchete_operator_identity',{}).get('direct_equivalence_status')=='ExactPlusEquality',
                'Native Warsaw sign identity is not verified')
        v12=import_v12(root)
        models={stem:analyze_one(root,audits,stem,v12) for stem in STEMS}
        report={
          'format':'t3-bc-F-hard-formal-boundary-plus-Y-conditional-ledger-v1',
          'scope':'original fermion-first one-loop F-hard LLHH, separately from already verified conditional Y bubble',
          'source_v12_summary_sha256':sha256(v12path),
          'source_v11_summary_sha256':sha256(v11path),
          'models':models,
          'F_hard_constant_and_Y_conditional_are_SEPARATE_boundaries':True,
          'no_naive_F_plus_Y_sum_or_double_counting_assumed':True,
          'requires_F_EFT1_CT_and_UV_IR_finite_matching':True,
          'requires_running_and_inherited_operator_matching':True,
          'physical_C5':None,
          'full_physical_matching_certified':False,
        }
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        for stem,item in models.items():
            print(f"{item['model']}: F hard 68 pole + 68 log + 68 finite; PR/PL component coefficients extracted; pole=log PASS")
        print('F hard finite is FORMAL; no F/Y addition, physical C5=NOT CERTIFIED')
        print('Saved: '+str(args.output.resolve()))
        return 0
    except (BoundaryError,KeyError,ValueError,TypeError,OSError,ImportError) as exc:
        print('FAIL (no physical C5): '+str(exc),file=sys.stderr)
        return 1

if __name__=='__main__':
    sys.exit(main())
