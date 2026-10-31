#!/usr/bin/env python3
"""v22: compare independent direct all-heavy Matchete ordered LLHH A with
source-matched F/Y/X full-mass triangle A, for archived B/C(-1).

Three Yukawa flavour monomials (y1_i y1_j, y1_i y2_j, y2_i y2_j) and all
three holomorphic quartics are checked independently, at three MF values.
Only exact scalar masses and a single loop order are compared. This is NOT
an independent Feynman-rule derivation or complete running/matching certificate.

Reads source files; never writes to GitHub or production code.
"""
from __future__ import annotations
import argparse
import ast
from collections import Counter, defaultdict
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import re
import sys

MODELS={'B_m1':('B_-1',1),'C_m1':('C_-1',3)}
BASE=Path('output/general_t3_mass_eigen_matchete')
AUD=Path('output/audits')
Q_NAMES=('lambdaT3','lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1')
Y_RE=re.compile(r'Bar\[Coupling\[(y[12]),\s*\{Index\[(d\$\$[12]), Flavor\],\s*Index\[d\$\$1, NFlavor\]\},\s*0\]\]')
BAR_C_RE=re.compile(r'Bar\[Coupling\[([A-Za-z][A-Za-z0-9]*),\s*\{\},\s*0\]\]')
C_RE=re.compile(r'Coupling\[([A-Za-z][A-Za-z0-9]*),\s*\{\},\s*0\]')
LOG_RE=re.compile(r'Log\[Coupling\[MF,\s*\{\},\s*0\]\^2/Coupling\[(MSX|MSY),\s*\{\},\s*0\]\^2\]')
LEGAL_SYMBOLS=frozenset({'hbar','SQRT3','L_FX','L_FY','MF','MSX','MSY',
    *(f'T3MEU{i}{j}' for i in (1,2) for j in (1,2)),
    *(f'bar_T3MEU{i}{j}' for i in (1,2) for j in (1,2)),
    *Q_NAMES,*(f'bar_{x}' for x in Q_NAMES)})
SHAPES=('y1_i*y1_j','y1_i*y2_j','y2_i*y2_j')

class DirectComparisonError(ValueError): pass

def need(condition,message):
    if not condition:raise DirectComparisonError(message)

def file_sha(path:Path)->str:
    need(path.is_file(),'Required source is missing: '+str(path))
    return sha256(path.read_bytes()).hexdigest()

def read_json(path:Path,fmt:str|None=None):
    need(path.is_file(),'Required JSON evidence missing: '+str(path))
    x=json.loads(path.read_text(encoding='utf-8-sig'))
    need(isinstance(x,dict),'Required JSON object: '+str(path))
    if fmt:need(x.get('format')==fmt,'Unexpected evidence format: '+str(path))
    return x

def import_module(path:Path,name:str):
    need(path.is_file(),'Required prior script missing: '+str(path))
    spec=importlib.util.spec_from_file_location(name,path)
    need(spec is not None and spec.loader is not None,'Cannot import '+str(path))
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def matching_close(text:str, start:int=0, target='('):
    need(text[start]==target,'Expected opening '+target)
    stack=[];pairs={')':'(',']':'[','}':'{'}
    for i in range(start,len(text)):
        c=text[i]
        if c in '([{':stack.append(c)
        elif c in ')]}':
            need(stack and stack.pop()==pairs[c],f'Unbalanced InputForm at offset {i}')
            if not stack:return i
    raise DirectComparisonError('Unterminated grouped InputForm')

def split_sum(text:str):
    """Split native InputForm on top-level + / -; no eval, no regex-only split."""
    stack=[];parts=[];start=0;pairs={')':'(',']':'[','}':'{'}
    for i,c in enumerate(text):
        if c in '([{':stack.append(c)
        elif c in ')]}':
            need(stack and stack.pop()==pairs[c],f'Unbalanced InputForm at offset {i}')
        elif c in '+-' and not stack and i>start:
            s=text[start:i].strip()
            if s:parts.append(s)
            start=i
    need(not stack,'Unterminated InputForm sum')
    last=text[start:].strip()
    if last:parts.append(last)
    return parts

def direct_terms(source:str,model:str):
    text=source.strip()
    need(text and 'hbar' in text,'Missing direct one-loop source')
    # C triplet source has one outer group /3 from the exact isospin trace.
    global_factor=1
    if model=='C_-1':
        need(text.startswith('(') and text.endswith(')/3'),
             'C direct coefficient normalization changed: expected outer /3')
        need(matching_close(text)==len(text)-3,
             'C source has unexpected top-level structure')
        text=text[1:-3].strip();global_factor=3
    parts=split_sum(text)
    need(len(parts)==480,f'{model}: expected 480 full-mass direct coefficient monomials, got {len(parts)}')
    return parts,global_factor

def compile_direct_term(source:str):
    """Classify Yukawa monomial and compile only a numerical arithmetic AST.
    No Python code execution, no arbitrary source calls, and no implicit symbols.
    """
    matches=list(Y_RE.finditer(source))
    need(len(matches)==2,'Direct term does not contain exactly two barred Yukawas')
    slots=sorted(((m.group(2),m.group(1)) for m in matches))
    need([fl for fl,_ in slots]==['d$$1','d$$2'], 'Unsupported direct lepton flavour index order')
    shape=slots[0][1]+'_i*'+slots[1][1]+'_j'
    need(shape in SHAPES,'Unexpected direct Yukawa pair: '+shape)
    quartet=re.findall(r'(?:Bar\[)?Coupling\[(lambdaT3(?:ExtraQ[0-9]+I[0-9]+)?),\s*\{\},\s*0\]',source)
    need(len(quartet)==1 and quartet[0] in Q_NAMES,'Unexpected quartic source in direct term')
    expr=LOG_RE.sub(lambda m:'L_FX' if m.group(1)=='MSX' else 'L_FY',source)
    expr=expr.replace('Sqrt[3]','SQRT3')
    expr=Y_RE.sub('1',expr)
    expr=BAR_C_RE.sub(lambda m:'bar_'+m.group(1),expr)
    expr=C_RE.sub(lambda m:m.group(1),expr)
    need(not any(c in expr for c in '[]{}$'),
         'Unknown Matchete InputForm function after strict replacements: '+expr[:170])
    expr=expr.replace('^','**')
    try:tree=ast.parse(expr,mode='eval').body
    except SyntaxError as exc:raise DirectComparisonError('Direct expression is not a safe arithmetic AST: '+expr[:150]) from exc
    check_arithmetic_ast(tree)
    return shape,quartet[0],tree

def check_arithmetic_ast(node):
    if isinstance(node,ast.BinOp):
        need(isinstance(node.op,(ast.Add,ast.Sub,ast.Mult,ast.Div,ast.Pow)),'Unsupported AST arithmetic')
        check_arithmetic_ast(node.left);check_arithmetic_ast(node.right)
    elif isinstance(node,ast.UnaryOp):
        need(isinstance(node.op,(ast.UAdd,ast.USub)),'Unsupported AST unary operator')
        check_arithmetic_ast(node.operand)
    elif isinstance(node,ast.Name):
        need(node.id in LEGAL_SYMBOLS,'Unrecognized Matchete coupling or function '+node.id)
    elif isinstance(node,ast.Constant):
        need(type(node.value) in (int,float),'Only numeric constants allowed')
    else:raise DirectComparisonError('Forbidden expression syntax: '+node.__class__.__name__)

def arith(node,variables,mp):
    if isinstance(node,ast.Name):return variables[node.id]
    if isinstance(node,ast.Constant):return mp.mpf(str(node.value))
    if isinstance(node,ast.UnaryOp):
        z=arith(node.operand,variables,mp)
        return z if isinstance(node.op,ast.UAdd) else -z
    x=arith(node.left,variables,mp);y=arith(node.right,variables,mp)
    if isinstance(node.op,ast.Add):return x+y
    if isinstance(node.op,ast.Sub):return x-y
    if isinstance(node.op,ast.Mult):return x*y
    if isinstance(node.op,ast.Div):return x/y
    if isinstance(node.op,ast.Pow):
        need(y==int(y) and abs(int(y))<=30,'Unexpected noninteger or large exponent')
        return x**int(y)
    raise DirectComparisonError('Unhandled safe AST type')

def compile_source(source,model):
    terms,denom=direct_terms(source,model)
    compiled=[compile_direct_term(term) for term in terms]
    counts=Counter(t[0] for t in compiled)
    need(counts=={'y1_i*y1_j':144,'y1_i*y2_j':192,'y2_i*y2_j':144},
         'Direct A5 flavour monomial counts changed: '+str(counts))
    byq=Counter(t[1] for t in compiled)
    need(byq=={'lambdaT3':192,'lambdaT3ExtraQ020020I1':144,
               'lambdaT3ExtraQ020200I1':144},'Unexpected direct quartic census '+str(byq))
    return compiled,denom,dict(counts),dict(byq)

def mp_complex(x,mp):return mp.mpc(str(float(x.real)),str(float(x.imag)))

def direct_ordered(compiled,denom,MF,MX,MY,U,quartics):
    """Return direct A monomial coefficients with hbar=1/(16 pi^2), units GeV^-1."""
    import mpmath as mp
    with mp.workdps(65):
        A=mp.mpf(str(MF));X=mp.mpf(str(MX));Y=mp.mpf(str(MY))
        need(A>Y>X>0,'Require MF > MY > MX > 0 for this source comparison')
        variables={'MF':A,'MSX':X,'MSY':Y,'hbar':1/(16*mp.pi**2),
                    'SQRT3':mp.sqrt(3),'L_FX':mp.log(A*A/(X*X)),
                    'L_FY':mp.log(A*A/(Y*Y))}
        for i in (1,2):
            for j in (1,2):
                key=f'T3MEU{i}{j}'
                z=mp_complex(complex(U[i-1,j-1]),mp)
                variables[key]=z;variables['bar_'+key]=mp.conj(z)
        for name in Q_NAMES:
            z=mp_complex(complex(quartics[name]),mp)
            variables[name]=z;variables['bar_'+name]=mp.conj(z)
        sums=defaultdict(lambda:mp.mpc(0))
        # Group same term structures before converting to float: the source has
        # large cancelling rational mass terms, so binary64 is unsafe here.
        for name,quartic,tree in compiled:
            sums[name]+=arith(tree,variables,mp)
        return tuple(complex(sums[key]/denom) for key in SHAPES)

def sequential_weights(root,stem,v16,v17,v18):
    """Recompute full-mass XX/XY/YY PL source kernels; never trust a stored ratio."""
    import sympy as sp
    v14=v16.read_json(root/AUD/'t3_bc_sequential_batch_v14/summary.json',
                      't3-bc-sequential-matching-fixed-order-batch-v14')
    extraction=v16.read_json(root/AUD/'t3_bc_scalar_finite_term_extraction_v1/summary.json',
                      't3-bc-native-local-matchete-pole-log-split-v1')
    deps=v16.import_deps(root)
    _,kernel,v12,v13,_=deps
    _,pol=v16.compare_finite_and_poles(root,stem,v14,extraction,deps)
    model=MODELS[stem][0]
    path=root/BASE/'eft1_fermion_first_oneloop'/f'{stem}_final_eft1_fermion_first_oneloop_raw_hard_llhh_sector.inputform.txt'
    fsec=v13.split_sectors(path.read_text(encoding='utf-8-sig'))
    kf=v16.component_polynomials(fsec['constant'],model,kernel,v12,68)
    pf,py,px=pol['F'],pol['Y'],pol['X']
    ky=tuple(sp.expand(pol['constant'][i]-kf[i]) for i in range(3))
    return tuple(tuple(sp.expand(z) for z in val) for val in (
       tuple(-v for v in px),tuple(-py[i]+ky[i] for i in range(3)),tuple(-v for v in ky)))

def expected_ordered(weights,MF,MX,MY,v17):
    import numpy as np
    A=MF*MF;X=MX*MX;Y=MY*MY
    phis=np.asarray([float(v17.signed_triangle_master(A,X,X)),
                     float(v17.signed_triangle_master(A,Y,Y)),
                     float(v17.signed_triangle_master(A,X,Y))])
    k=sum(weights[i]*phis[i] for i in range(3))/(16*math.pi**2*MF)
    return tuple(complex(x) for x in k)

def substitution(U,quartics):
    import sympy as sp
    sub={}
    for i in range(2):
        for j in range(2):
            n=f'T3MEU{i+1}{j+1}'
            sub[sp.Symbol(n)]=U[i,j]
            sub[sp.Symbol('bar_'+n)]=U[i,j].conjugate()
    for n,z in quartics.items():
        sub[sp.Symbol(n)]=z
        sub[sp.Symbol('bar_'+n)]=z.conjugate()
    return sub

def consistency(direct,expected):
    import numpy as np
    a=np.asarray(direct,complex);b=np.asarray(expected,complex)
    scale=max(float(np.linalg.norm(b)),1e-250)
    return (float(np.linalg.norm(a-b)/scale),
            [[float(x.real),float(x.imag)] for x in a],
            [[float(x.real),float(x.imag)] for x in b])

def source_gates(root,stem):
    model=MODELS[stem][0]
    v21=read_json(root/AUD/'t3_bc_topology_closure_v21/summary.json',
                  't3-bc-one-loop-renormalizable-topology-closure-v21')
    v19=read_json(root/AUD/'t3_bc_hard_soft_subtraction_v19/summary.json',
                  't3-bc-one-loop-hard-soft-msbar-restricted-v19')
    need(v21.get('one_loop_diagram_topology_complete_for_archived_UV_under_scope') is True,
         'v21 restricted graph completeness not established')
    need(v21['models'][stem].get('one_loop_two_trivalent_graphs_excluded_for_archived_UV') is True,
         model+': two-trivalent closure missing')
    need(v19.get('full_T3_renormalized_matching_certified') is False and
         v19.get('full_T3_physical_C5') is None,
         'v19 certification unexpectedly changed')
    direct_path=root/BASE/'one_loop_probe'/f'{stem}_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt'
    report_path=root/BASE/'one_loop_probe'/f'{stem}_final_one_loop_probe.json'
    report=read_json(report_path,'t3-matchete-native-heavy-mass-eigenbasis-v1')
    need(report.get('model')==model and report.get('one_loop_match_stage_completed') is True
         and report.get('weinberg_extraction_status')=='Success'
         and report.get('weinberg_term_count')==236
         and report.get('weinberg_holomorphic_term_count')==118
         and report.get('weinberg_conjugate_term_count')==118,
         model+': one-loop all-heavy native matching source changed')
    need(report.get('ordered_coefficient_is_canonical_flavour_symmetric_c5') is False,
         'Double counting risk: direct A is ordered, not symmetric C5')
    need(report.get('production_sequential_matching_performed') is False and
         report.get('numerical_c5_computed') is False,
         model+': direct benchmark certification semantics changed')
    need(report.get('all_17_quartics_present') is True and
         report.get('both_yukawa_interactions_present') is True and
         report.get('heavy_mass_basis_detailed_check_passed') is True,
         model+': inconsistent direct UV mass basis')
    need(report.get('weinberg_ordered_coefficient_inputform_file')==direct_path.name,
         model+': direct source filename mismatch')
    return direct_path,report_path,v19,v21

def scenarios(baseline):
    return {
       'baseline_only':{Q_NAMES[0]:baseline,Q_NAMES[1]:0j,Q_NAMES[2]:0j},
       'extra_y2y2_only':{Q_NAMES[0]:0j,Q_NAMES[1]:0.021+0.007j,Q_NAMES[2]:0j},
       'extra_y1y1_only':{Q_NAMES[0]:0j,Q_NAMES[1]:0j,Q_NAMES[2]:-0.014+0.005j},
       'all_three_complex':{Q_NAMES[0]:baseline,Q_NAMES[1]:0.021+0.007j,
                            Q_NAMES[2]:-0.014+0.005j},
    }

def run(root,out):
    import numpy as np
    v16=import_module(root/'scripts/AssembleT3BCFormalFiniteScaleEnvelope.experimental.py','t3v22_v16')
    v17=import_module(root/'scripts/CompareT3BCFiniteEnvelopeToTriangle.experimental.py','t3v22_v17')
    v18=import_module(root/'scripts/EvaluateT3BCThreeFermionTriangle.experimental.py','t3v22_v18')
    v18report=read_json(root/AUD/'t3_bc_three_fermion_triangle_v18/summary.json',
            't3-bc-frozen-reference-three-fermion-triangle-v18')
    need(v18report.get('physical_C5') is None and v18report.get('running_performed') is False,
         'Unexpected v18 certification flags')
    result={'format':'t3-bc-independent-all-heavy-vs-sequential-finite-amplitude-v22',
        'status':'DirectAllHeavyVsSequentialComparisonCompleted',
        'scope':'archived B/C one-loop ordered PL Weinberg, fixed-order all-heavy vs F/Y/X full three-mass triangle',
        'source_direct_calculation':'Matchete Match[UV,5,1] all F/X/Y heavy; independent matching hierarchy from sequential Matchete sectors',
        'normalization':'A source ordered; C5=A+A^T; hbar=1/(16*pi^2)',
        'benchmark_type':'3 nondegenerate fermion masses frozen at reference input 1e7 GeV; synthetic independent complex quartic axes',
        'full_running_performed':False,
        'independent_non_Matchete_UV_Feynman_rule_calculation_performed':False,
        'complete_EFT1_EFT2_finite_field_EOM_scheme_certified':False,
        'physical_C5':None,'physical_neutrino_masses':None,
        'models':{}}
    all_pass=True
    for stem,(model,_) in MODELS.items():
        direct_path,report_path,v19,v21=source_gates(root,stem)
        data=direct_path.read_text(encoding='utf-8-sig')
        compiled,denom,counts,quartic_counts=compile_source(data,model)
        weights_sym=sequential_weights(root,stem,v16,v17,v18)
        params=read_json(root/BASE/f'{stem}_final_parameters.json')
        need(params.get('model')==model,'Scalar mass-basis model changed')
        mx,my=map(float,params['heavy_masses_gev'])
        U=v18.as_matrix(params['unitary_U'],(2,2),'unitary_U')
        need(np.linalg.norm(U.conj().T@U-np.eye(2))<1e-10,'Unitarity failure')
        ref=v18.load_frozen_reference(root,stem)
        per_case=[];max_error=0;baseline_direct={}
        for label,quartics in scenarios(ref['quartics'][Q_NAMES[0]]).items():
            sub=substitution(U,quartics)
            numeric_weights=tuple(v18.polynomial_values(z,sub) for z in weights_sym)
            for index,mf in enumerate(ref['mf']):
                direct=direct_ordered(compiled,denom,float(mf),mx,my,U,quartics)
                if label=='baseline_only':baseline_direct[index]=direct
                expected=expected_ordered(numeric_weights,float(mf),mx,my,v17)
                rel,raw,comparison=consistency(direct,expected)
                max_error=max(max_error,rel)
                per_case.append({'scenario':label,'MF_index':index,'MF_GeV':float(mf),
                   'relative_euclidean_residual':rel,
                   'direct_ordered_A_coefficients_GeV_inv':[SHAPES[i] for i in range(3)],
                   'direct':raw,'sequential_full_mass':comparison,
                   'passed':rel<2e-8})
        # Independent direct-A -> complex-symmetric 3x3 Weinberg map, compared
        # with the frozen full-mass v18 result (not a fitted physical mass).
        direct_C5=np.zeros((3,3),complex)
        for col in range(3):
            a11,a12,a22=baseline_direct[col]
            K=np.array([[2*a11,a12],[a12,2*a22]],complex)
            Y=np.column_stack((ref['y1'][:,col].conj(),ref['y2'][:,col].conj()))
            direct_C5+=Y@K@Y.T
        old=v18report['models'][stem]['strict_and_full_triangle']['full_mass']
        stored=np.array([[complex(*z) for z in row]
                         for row in old['reference_C5_3x3_GeV_inverse_NOT_PHYSICAL']],complex)
        c5_residual=float(np.linalg.norm(direct_C5-stored)/max(np.linalg.norm(stored),1e-250))
        # Direct A vs v18 mass eigenvalues is an absolute sign and 1/2 test
        # at the matrix level; singular values alone would miss an overall sign.
        direct_mnu=-0.5*ref['vev_gev']**2*direct_C5*1e9
        singular=np.sort(np.linalg.svd(direct_mnu,compute_uv=False))
        passed=all(x['passed'] for x in per_case) and c5_residual<2e-8
        all_pass&=passed
        result['models'][stem]={'model':model,
          'one_loop_direct_coefficient_sha256':file_sha(direct_path),
          'all_heavy_direct_json_sha256':file_sha(report_path),
          'native_direct_term_count':len(compiled),
          'source_flavour_inventory':counts,'source_quartic_inventory':quartic_counts,
          'total_independent_checks':len(per_case),
          'separate_quartic_axes_checked':True,
          'three_fermion_masses_checked':True,
          'direct_vs_sequential_max_relative_residual':max_error,
          'direct_full_3x3_C5_vs_v18_relative_residual':c5_residual,
          'direct_full_3x3_C5_GeV_inverse_NOT_PHYSICAL':v18.export_matrix(direct_C5),
          'direct_reference_neutrino_singular_masses_eV_NOT_PHYSICAL':[float(x) for x in singular],
          'direct_vs_sequential_full_mass_agreement':passed,
          'checks':per_case,
          'full_physical_C5_certified':False}
        print(f'{model}: direct 480 terms; {len(per_case)} independent comparisons; '
              f'max relative residual {max_error:.5g}; agreement '+('PASS' if passed else 'FAIL'))
    result['all_direct_vs_sequential_full_mass_agree']=bool(all_pass)
    if not all_pass:result['status']='DirectVsSequentialFiniteResidualNonzero'
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('Saved:',out)
    if not all_pass:raise DirectComparisonError('DIRECT all-heavy finite coefficient disagrees with F/Y/X full-mass result: inspect JSON residuals. No physical C5.')
    print('PASS: direct all-heavy vs sequential full-mass A numerical equality in tested channels; NO complete physical C5')

def self_test():
    import mpmath as mp
    specimen='- (hbar*Bar[Coupling[lambdaT3ExtraQ020020I1, {}, 0]]*Bar[Coupling[y2, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]]*Bar[Coupling[y2, {Index[d$$2, Flavor], Index[d$$1, NFlavor]}, 0]]*Coupling[MF, {}, 0]^3*Log[Coupling[MF, {}, 0]^2/Coupling[MSY, {}, 0]^2])/(Coupling[MF, {}, 0]^2-Coupling[MSX, {}, 0]^2)^2'
    shape,quartic,astnode=compile_direct_term(specimen)
    need(shape=='y2_i*y2_j' and quartic==Q_NAMES[1], 'Direct safe parser Yukawa/quartic failure')
    with mp.workdps(55):
        vars={'MF':mp.mpf(10),'MSX':mp.mpf(3),'MSY':mp.mpf(4),
             'hbar':mp.mpf(1),'L_FX':mp.log(mp.mpf(100)/9),
             'L_FY':mp.log(mp.mpf(100)/16),'SQRT3':mp.sqrt(3),
             'lambdaT3ExtraQ020020I1':mp.mpf(2),
             'bar_lambdaT3ExtraQ020020I1':mp.mpf(2)}
        found=arith(astnode,vars,mp)
        expected= -2*mp.mpf(10)**3*mp.log(mp.mpf(100)/16)/(100-9)**2
        need(abs(found-expected)<mp.mpf('1e-40'),'Independent arithmetic parser miscalculates')
    # Full 480-term synthetic census, including native C_-1 outer /3.
    shapes=[]
    for count,yi,yj,qn in (
        (144,'y1','y1','lambdaT3ExtraQ020200I1'),
        (192,'y1','y2','lambdaT3'),
        (144,'y2','y2','lambdaT3ExtraQ020020I1')):
        t=specimen.replace('lambdaT3ExtraQ020020I1',qn)
        # first flavour slot d$$1, second d$$2
        t=t.replace('Coupling[y2, {Index[d$$1',f'Coupling[{yi}, {{Index[d$$1')
        t=t.replace('Coupling[y2, {Index[d$$2',f'Coupling[{yj}, {{Index[d$$2')
        # Source terms have either a leading + or a leading - before each group.
        shapes.extend(['('+t.lstrip('- ').strip()+')']*count)
    for name,synthetic in [('B_-1',' + '.join(shapes)),
                           ('C_-1','('+' + '.join(shapes)+')/3')]:
        parsed,denom,counts,qcounts=compile_source(synthetic,name)
        need(len(parsed)==480 and denom==(3 if name=='C_-1' else 1),
             'Synthetic direct 480-term census or SU(2) prefactor failed')
        need(counts['y1_i*y2_j']==192 and qcounts['lambdaT3']==192,
             'Wrong synthetic Yukawa/quartic split')
    # forbid arbitrary code, wrong Lorentz function, wrong Yukawa conjugation
    for broken in (specimen.replace('Log[Coupling[MF, {}, 0]^2/Coupling[MSY, {}, 0]^2]','Run["x"]'),
                   specimen.replace('Bar[Coupling[y2,','Coupling[y2,'),
                   specimen.replace('Coupling[MF, {}, 0]^3','Coupling[ForbiddenSymbol, {}, 0]^3')):
        try:compile_direct_term(broken)
        except DirectComparisonError:pass
        else:raise DirectComparisonError('Tampered direct InputForm was accepted')
    fake='(('+specimen+') + ('+specimen+'))/3'
    # Test only outer C normalization; actual direct source has 480 terms.
    need(fake.startswith('(') and matching_close(fake)==len(fake)-3,'C normalization shell failed')
    import numpy as np
    # A ordered -> Weinberg symmetric (not A+A.conj().T).
    a11=1+2j;a12=3-4j;a22=5+6j
    A=np.array([[a11,a12],[0,a22]],complex)
    C=A+A.T
    need(np.allclose(C,C.T) and C[0,1]==a12 and C[0,0]==2*a11,'Weinberg symmetrization incorrect')
    print('PASS: safe Matchete direct InputForm arithmetic and C triplet prefactor')
    print('PASS: 3 ordered flavour channels, arbitrary code/tamper rejection, Majorana transpose')
    print('PASS: physical-C5 certification remains blocked')

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--output',type=Path,default=AUD/'t3_bc_direct_vs_sequential_full_mass_v22/summary.json')
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--require-physical',action='store_true')
    a=ap.parse_args()
    try:
        need(not a.require_physical,'Full physical C5 requires matched running and full renormalized UV/EFT scheme')
        if a.self_test:self_test()
        else:
            root=a.root.resolve();out=a.output if a.output.is_absolute() else root/a.output
            run(root,out)
        return 0
    except (DirectComparisonError,KeyError,TypeError,ValueError,IndexError,ImportError,OSError,
            ZeroDivisionError,OverflowError) as exc:
        print('FAIL (no physical C5): '+str(exc),file=sys.stderr)
        return 2

if __name__=='__main__':raise SystemExit(main())
