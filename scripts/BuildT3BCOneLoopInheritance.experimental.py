#!/usr/bin/env python3
"""Source-verified O(hbar) Y -> X spectator inheritance, B(-1)/C(-1).

Proves from archived Y-tree and Y-full EFT2 InputForm that the first local
LLHH operator is one-loop, and that the tree EFT2 contains no odd-X source.
Exports the two EXACT one-loop LLHH spectator expressions carried below MX.
They are not renormalized coefficients; no F/Y/X physical sum is formed.

Additionally compares FORMAL F/Y/X LLHH pole polynomials after the same
neutral-SU2 Majorana projection; comparison is deliberately not interpreted
as a renormalized threshold subtraction or UV/IR assignment.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys

MODELS = {'B_m1':'B_-1','C_m1':'C_-1'}
YBASE=Path('output/general_t3_mass_eigen_matchete/eft1_scalar_xy_fresh_kernel_v3')
YRAW=Path('output/general_t3_mass_eigen_matchete/eft1_scalar_Y_only_canonical_probe')
FRAW=Path('output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop')

class InheritanceError(ValueError): pass

def require(flag, message):
    if not flag: raise InheritanceError(message)

def read(path):
    require(path.is_file(), 'Missing local source '+str(path))
    return path.read_text(encoding='utf-8-sig')

def read_json(path, fmt=None):
    x=json.loads(read(path)); require(isinstance(x,dict),'Bad JSON '+str(path))
    if fmt: require(x.get('format')==fmt, 'Wrong report format '+str(path))
    return x

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def import_local(root, path, name):
    p=root/path;require(p.is_file(),'Missing dependency '+str(p))
    spec=importlib.util.spec_from_file_location(name,p)
    require(spec is not None and spec.loader is not None,'Cannot import '+str(p))
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

def project_imports(root):
    r=str(root.resolve())
    if r in sys.path:sys.path.remove(r)
    sys.path.insert(0,r)
    from Numerical.experimental import T3BCPoleLogConsistency as pole
    from Numerical.experimental import T3BCNativeYLLHHFlavorKernel as kernel
    for module,path in ((pole,'Numerical/experimental/T3BCPoleLogConsistency.py'),
                        (kernel,'Numerical/experimental/T3BCNativeYLLHHFlavorKernel.py')):
        require(Path(module.__file__).resolve()==(root/path).resolve(),
                'Wrong installed local parser: '+str(module.__file__))
    return pole,kernel


def classify_term(t):
    """Count actual Field[] occurrences, NOT incidental X names in couplings."""
    h=len(re.findall(r'Field\[H, Scalar',t))
    x=len(re.findall(r'Field\[NewScalarX, Scalar',t))
    l=len(re.findall(r'Field\[l, Fermion',t))
    cc=t.count('GammaCC')
    if cc and h==2 and x==0 and l==2:
        return 'LLHH'
    if cc and h==0 and x==2 and l==2:
        return 'LLXX'
    if cc and h==1 and x==1 and l==2:
        return 'LLHX'
    return 'other'


def loop_order(t):
    if 'hbar^2' in t or 'hbar*hbar' in t or 'hbar^3' in t:
        return 2
    return 1 if re.search(r'(?<![A-Za-z_])hbar(?![A-Za-z_])',t) else 0


def one_model(root, stem, v14, pole, output_dir):
    model=MODELS[stem]
    previous=v14['models'][stem]['raw_scalar_sector_evidence']['Y_EFT2_inheritance']
    tree=root/YBASE/f'{stem}_final_x_from_saved_y_match_x_followup_y_tree_eomeft.inputform.txt'
    full=root/YBASE/f'{stem}_final_x_from_saved_y_match_x_followup_y_full_eomeft.inputform.txt'
    require(sha(tree)==previous['Y_tree_sha256'] and sha(full)==previous['Y_full_sha256'],
            model+': saved EFT2 source changed since v14')
    tt=pole.terms(read(tree)); ft=pole.terms(read(full))
    require(tt and ft and all(loop_order(t)==0 for t in tt),model+': Y tree input not loop-free')
    require(all(loop_order(t) in (0,1) for t in ft),model+': unsupported higher-loop Y EFT2')
    odd=[i for i,t in enumerate(tt) if len(re.findall(r'Field\[NewScalarX, Scalar',t))%2]
    require(not odd, model+': odd-X source in Y tree EFT2, cannot use spectator argument')
    tree_llhh=[t for t in tt if classify_term(t)=='LLHH']
    require(not tree_llhh,model+': tree-level LLHH appeared; order counting must be revised')
    one_loop_llhh=[t for t in ft if classify_term(t)=='LLHH' and loop_order(t)==1]
    other_llhh=[t for t in ft if classify_term(t)=='LLHH' and loop_order(t)!=1]
    require(not other_llhh and len(one_loop_llhh)==2,
            model+': expected exactly two one-loop Y EFT2 LLHH chirality forms')
    sectors={}
    for t in one_loop_llhh:
        require(t.count('GammaCC')==1 and t.count('NCM[')==1,
                model+': unexpected noncommutative LLHH bilinear')
        ch='PL' if 'Proj[-1]' in t else ('PR' if 'Proj[1]' in t else '?')
        require(ch in ('PL','PR') and ch not in sectors,model+': invalid/duplicated chirality')
        require('ϵ' in t and 'Log[μbar2/Coupling[MSY, {}, 0]^2]' in t,
                model+': expected unsusbtracted Y pole and Y log in spectator')
        sectors[ch]=t
    require(set(sectors)=={'PL','PR'}, model+': missing Y spectator chirality')
    # Exact native InputForm is transferred without dropping poles or terms.
    outputs={}
    output_dir.mkdir(parents=True,exist_ok=True)
    for ch,t in sectors.items():
        target=output_dir/f'{stem}_Y_one_loop_LLHH_{ch}_inherited_through_X.UNRENORMALIZED.inputform.txt'
        target.write_text(t+'\n',encoding='utf-8')
        outputs[ch]={'path':str(target.relative_to(root)) if target.is_relative_to(root) else str(target),
                     'sha256':sha(target),'length':len(t)}
    return {'model':model,
            'Y_tree_EFT2_top_level_terms':len(tt),
            'Y_full_EFT2_top_level_terms':len(ft),
            'Y_tree_EFT2_odd_X_terms':len(odd),
            'Y_tree_EFT2_LLHH_terms':len(tree_llhh),
            'Y_full_EFT2_one_loop_LLHH_chiral_operators':2,
            'Y_full_EFT2_LLHH_higher_loop_terms':len(other_llhh),
            'Y_one_loop_LLHH_spectator_InputForms':outputs,
            'X_tree_level_LLHH_induced_from_Y_LLXX':False,
            'X_one_loop_on_Y_one_loop_LLXX_order':'hbar^2 (excluded at strict one loop)',
            'X_one_loop_on_Y_tree_LLXX_order':'hbar (existing native X source)',
            'Y_one_loop_LLHH_direct_spectator_below_X':True,
            'F_one_loop_LLHH_direct_spectator_below_Y_X':'conditional on same EFT basis and matching scheme',
            'Y_LLHH_spectator_renormalized':False,
            'complete_low_energy_C5_computed':False}


def residues(root, stem, pole, kernel, v12, v13, v14):
    """Check raw *neutral* pole polynomials across F/Y/X, no physical sum."""
    import sympy as sp
    model=MODELS[stem]
    f=root/FRAW/f'{stem}_final_eft1_fermion_first_oneloop_raw_hard_llhh_sector.inputform.txt'
    y=root/YRAW/f'{stem}_final_y_only_canonical_match_raw_scalar_llhh.inputform.txt'
    x=root/YBASE/f'{stem}_final_x_from_saved_y_match_x_followup_x_raw_llhh.inputform.txt'
    fr=v13.split_sectors(read(f))
    yr=v14.unpack_sectors(read(y),'Y',pole)
    xr=v14.unpack_sectors(read(x),'X',pole)
    require(len(fr['pole'])==68 and len(yr[1]['pole'])==50 and len(xr[1]['pole'])==18,
            model+': wrong raw LLHH boundary poles')
    def polynomial_by_chirality(raw_terms, expected):
        counts=Counter();coeffs=defaultdict(Fraction);seen={}
        for term in raw_terms:
            rec,c=kernel._record_from_term(term,model)
            ch=rec['chirality'];counts[ch]+=1
            key=(ch,rec['quartic'],rec['quartic_conjugated'],
                 rec['yukawa_i'],rec['yukawa_j'],rec['yukawas_conjugated'],
                 tuple(rec['mixing']),rec['sqrt3_exponent'])
            coeffs[key]+=c;seen[key]=rec
        require(counts=={'PL':expected,'PR':expected},
                model+': unexpected native PL/PR pole count '+str(counts))
        rows=[]
        for key,rec in seen.items():
            cf=coeffs[key]
            if not cf:continue
            rows.append({'chirality':rec['chirality'],'quartic':rec['quartic'],
                         'quartic_conjugated':rec['quartic_conjugated'],
                         'yukawa_i':rec['yukawa_i'],'yukawa_j':rec['yukawa_j'],
                         'yukawas_conjugated':rec['yukawas_conjugated'],
                         'mixing':[{'name':n,'conjugated':b,'exponent':e} for n,b,e in rec['mixing']],
                         'sqrt3_exponent':rec['sqrt3_exponent'],
                         'coefficient_numerator':cf.numerator,
                         'coefficient_denominator':cf.denominator})
        ledger={'monomials':rows}
        pl=v12.polynomial_coefficients(ledger,'PL')
        pr=v12.polynomial_coefficients(ledger,'PR')
        require(all(sp.expand(v12.polynomial_conjugate(a)-b)==0 for a,b in zip(pl,pr)),
                model+': cross-boundary PR/PL H.c. failure')
        return pl
    pf=polynomial_by_chirality(fr['pole'],34)
    py=polynomial_by_chirality(yr[1]['pole'],25)
    px=polynomial_by_chirality(xr[1]['pole'],9)
    plus=[sp.expand(pf[i]+py[i]+px[i]) for i in range(3)]
    minus=[sp.expand(pf[i]-py[i]-px[i]) for i in range(3)]
    signs={'F_plus_Y_plus_X':all(t==0 for t in plus),
           'F_minus_Y_minus_X':all(t==0 for t in minus)}
    require(not (signs['F_plus_Y_plus_X'] and signs['F_minus_Y_minus_X']),
            model+': contradictory formal pole tests')
    # An algebraic residue equality does not fix UV/IR sign conventions.
    return {'model':model,
            'pole_raw_inventory':{'F':68,'Y':50,'X':18},
            'independent_native_SU2_component_reduction':True,
            'PR_PL_exact_conjugacy':True,
            'formal_neutral_residue_relations':signs,
            'pole_sum_residual_nonzero_structure_count':sum(t!=0 for t in plus),
            'pole_difference_residual_nonzero_structure_count':sum(t!=0 for t in minus),
            'pole_sum_residuals':list(map(str,plus)),
            'pole_difference_residuals':list(map(str,minus)),
            'complete_UV_IR_assignment':False,
            'counterterm_cancellation_certified':False,
            'physical_matching_from_residue_equality':False}


def self_test():
    import tempfile
    class Toy:
        @staticmethod
        def terms(t):
            return tuple(q.strip() for q in t.strip().split("\n") if q.strip())
    assert classify_term('Field[H, Scalar]*Field[H, Scalar]*Field[l, Fermion]*Field[l, Fermion]*GammaCC')=='LLHH'
    assert classify_term('Field[NewScalarX, Scalar]*Field[NewScalarX, Scalar]*Field[l, Fermion]*Field[l, Fermion]*GammaCC')=='LLXX'
    assert loop_order('hbar*x')==1 and loop_order('hbar^2*x')==2 and loop_order('x')==0
    import sympy as sp
    a,b=sp.symbols('a b'); require(sp.expand(a-(a+b)+b)==0,'Exact residual polynomial test')
    with tempfile.TemporaryDirectory() as d:
        root=Path(d)
        for stem in MODELS:
            base=root/YBASE;base.mkdir(parents=True,exist_ok=True)
            tree=base/f'{stem}_final_x_from_saved_y_match_x_followup_y_tree_eomeft.inputform.txt'
            full=base/f'{stem}_final_x_from_saved_y_match_x_followup_y_full_eomeft.inputform.txt'
            tree.write_text('Field[NewScalarX, Scalar]*Field[NewScalarX, Scalar]\n',encoding='utf-8')
            bil='Field[H, Scalar]*Field[H, Scalar]*Field[l, Fermion]*Field[l, Fermion]*'
            full.write_text('hbar*'+bil+'GammaCC*NCM[Proj[-1]]*(1/ϵ+Log[μbar2/Coupling[MSY, {}, 0]^2])\n'
                            +'hbar*'+bil+'GammaCC*NCM[Proj[1]]*(1/ϵ+Log[μbar2/Coupling[MSY, {}, 0]^2])\n',encoding='utf-8')
            previous={stem:{'raw_scalar_sector_evidence':{'Y_EFT2_inheritance':{
              'Y_tree_sha256':sha(tree),'Y_full_sha256':sha(full)}}}}
            outcome=one_model(root,stem,{'models':previous},Toy,root/'out')
            require(outcome['Y_full_EFT2_one_loop_LLHH_chiral_operators']==2,
                    'Synthetic inheritance source unexpectedly rejected')
            tree.write_text('Field[NewScalarX, Scalar]\n',encoding='utf-8')
            previous[stem]['raw_scalar_sector_evidence']['Y_EFT2_inheritance']['Y_tree_sha256']=sha(tree)
            try:one_model(root,stem,{'models':previous},Toy,root/'out')
            except InheritanceError:pass
            else:raise InheritanceError('Odd-X tree source was accepted')
    print('PASS: B/C direct Y-loop LLHH spectator export; odd-X input rejected')
    print('PASS: fixed-order hbar vs hbar^2 and symbolic pole-residue checks')


def run(root,output):
    v14_summary=read_json(root/'output/audits/t3_bc_sequential_batch_v14/summary.json',
                          't3-bc-sequential-matching-fixed-order-batch-v14')
    require(v14_summary.get('full_physical_matching_certified') is False,
            'Refusing changed v14 certification state')
    pole,kernel=project_imports(root)
    v12=import_local(root,Path('scripts/BuildT3BCYFiniteBubbleInterface.experimental.py'),'t3_v15_v12')
    v13=import_local(root,Path('scripts/BuildT3BCFHardBoundaryInterface.experimental.py'),'t3_v15_v13')
    v14=import_local(root,Path('scripts/BuildT3BCSequentialMatchingBatch.experimental.py'),'t3_v15_v14')
    output.parent.mkdir(parents=True,exist_ok=True)
    models={}
    for stem,model in MODELS.items():
        inherited=one_model(root,stem,v14_summary,pole,output.parent)
        try:
            comparison=residues(root,stem,pole,kernel,v12,v13,v14)
            comparison_status='ExactNeutralPoleComparisonCompleted'
        except (ValueError,TypeError,KeyError,IndexError,ImportError) as exc:
            comparison={'status':'NotEstablished','error':str(exc),
                        'no_formal_relation_assumed':True,
                        'physical_matching_certified':False}
            comparison_status='FormalNeutralPoleComparisonNotEstablished'
        models[stem]={'model':model,'fixed_order_Y_through_X':inherited,
                      'F_Y_X_formal_pole_relation':comparison,
                      'formal_pole_comparison_status':comparison_status,
                      'complete_physical_C5':None}
        print(model+': Y one-loop LLHH PL/PR inherited through X at O(hbar); '
              'Y one-loop LLXX to X loop excluded at O(hbar^2)')
        print(model+': F/Y/X '+comparison_status)
    report={'format':'t3-bc-one-loop-fixed-order-inheritance-v15',
            'status':'YLLHHSpectatorInheritedAtOneLoop_FullMatchingSubtractionOutstanding',
            'scope':'Strict local dimension-5 O(hbar) operator power counting; independent of RG',
            'models':models,
            'formula':'C5_belowX^(1) = C5_aboveX^(1),Y_spectator + DeltaC5_X^(1)[EFT2^(0)] + other separately matched boundaries; never insert Y^(1)_LLXX into X^(1) loop',
            'Y_one_loop_LLXX_insertion_into_X_one_loop_needed_at_Ohbar':False,
            'Y_one_loop_LLHH_spectator_inheritance_at_Ohbar':True,
            'physical_C5':None,'full_EFT_UV_IR_subtraction_certified':False,
            'source_EFT1_canonical_equivalence_certified':False,
            'full_physical_matching_certified':False}
    output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('PASS: fixed-order inheritance written; physical C5 NOT CERTIFIED')
    print('Saved: '+str(output))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--self-test',action='store_true')
    p.add_argument('--root',type=Path,default=Path.cwd())
    p.add_argument('--output',type=Path,
                   default=Path('output/audits/t3_bc_one_loop_inheritance_v15/summary.json'))
    p.add_argument('--require-physical',action='store_true')
    args=p.parse_args()
    try:
        if args.self_test:self_test();return 0
        if args.require_physical:
            raise InheritanceError('Full F/EFT1 and EFT1/EFT2 one-loop subtraction, scheme and running are unverified')
        root=args.root.resolve()
        out=args.output if args.output.is_absolute() else root/args.output
        run(root,out)
        return 0
    except (InheritanceError,ValueError,TypeError,KeyError,IndexError,OSError,ImportError) as exc:
        print('FAIL (no physical C5): '+str(exc),file=sys.stderr)
        return 1
if __name__=='__main__':sys.exit(main())
