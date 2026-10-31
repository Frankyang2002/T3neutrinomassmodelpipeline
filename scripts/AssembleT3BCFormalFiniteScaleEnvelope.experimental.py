#!/usr/bin/env python3
"""B/C F+Y+X FORMAL finite common-scale envelope (NOT a physical Wilson coefficient).

Uses exactly the F/Y/X sources already audited in v12--v15. Proves in the
same neutrino/H2 component basis:
  P_F + P_Y + P_X == 0,
  K_formal = K_F + K_Y + K_X + P_F ln(M_X^2/M_F^2)
                                + P_Y ln(M_X^2/M_Y^2)
with common hbar/MF factored out. This is only a source-amplitude formal
identity; absent full high-minus-low 1PI, wavefunction, EOM and parameter
counterterms it is NOT delta C5^ren, physical C5, or m_nu.

Do not use this output as a production matching boundary.  The optional
numeric mode supplies a explicitly-labelled formal 3x3 matrix for debugging.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import sys

MODELS={'B_m1':'B_-1','C_m1':'C_-1'}
AUD=Path('output/audits')
BASE=Path('output/general_t3_mass_eigen_matchete')

class EnvelopeError(ValueError):pass

def require(ok, message):
    if not ok:raise EnvelopeError(message)

def read_json(path, fmt=None):
    require(path.is_file(), 'Missing source evidence: '+str(path))
    obj=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj,dict),'Expected JSON object: '+str(path))
    if fmt:require(obj.get('format')==fmt,'Evidence schema changed: '+str(path))
    return obj

def digest(path):return sha256(path.read_bytes()).hexdigest()

def import_script(root, rel, name):
    path=root/rel
    require(path.is_file(),'Missing prior Workstream C script: '+str(path))
    spec=importlib.util.spec_from_file_location(name,path)
    require(spec is not None and spec.loader is not None,'Cannot load '+str(path))
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def import_deps(root):
    here=str(root.resolve())
    if here in sys.path:sys.path.remove(here)
    sys.path.insert(0,here)
    from Numerical.experimental import T3BCPoleLogConsistency as pole
    from Numerical.experimental import T3BCNativeYLLHHFlavorKernel as kernel
    for mod,rel in ((pole,'Numerical/experimental/T3BCPoleLogConsistency.py'),
                    (kernel,'Numerical/experimental/T3BCNativeYLLHHFlavorKernel.py')):
        require(Path(mod.__file__).resolve()==(root/rel).resolve(),
                'Wrong repository dependency '+str(mod.__file__))
    v12=import_script(root,'scripts/BuildT3BCYFiniteBubbleInterface.experimental.py','v16_v12')
    v13=import_script(root,'scripts/BuildT3BCFHardBoundaryInterface.experimental.py','v16_v13')
    v14=import_script(root,'scripts/BuildT3BCSequentialMatchingBatch.experimental.py','v16_v14')
    return pole,kernel,v12,v13,v14

def component_polynomials(terms, model, kernel, v12, expect):
    """Native field/CG/Warsaw-PL projection; preserves signed source monomials."""
    import sympy as sp
    require(len(terms)==expect,'Unexpected source term count')
    if expect==0:return (sp.Integer(0),)*3
    coeffs=defaultdict(Fraction); records={};counts=Counter()
    for term in terms:
        record, coefficient=kernel._record_from_term(term,model)
        ch=record['chirality'];counts[ch]+=1
        key=(ch,record['quartic'],record['quartic_conjugated'],
             record['yukawa_i'],record['yukawa_j'],record['yukawas_conjugated'],
             tuple(record['mixing']),record['sqrt3_exponent'])
        coeffs[key]+=coefficient;records[key]=record
    require(counts=={'PL':expect//2,'PR':expect//2},
            'Missing Majorana conjugate partner '+str(counts))
    rows=[]
    for key in sorted(records):
        value=coeffs[key]
        if value==0:continue
        r=records[key]
        rows.append({**{n:r[n] for n in (
            'chirality','quartic','quartic_conjugated','yukawa_i','yukawa_j',
            'yukawas_conjugated','sqrt3_exponent')},
            'mixing':[{'name':n,'conjugated':b,'exponent':e} for n,b,e in r['mixing']],
            'coefficient_numerator':value.numerator,
            'coefficient_denominator':value.denominator})
    ledger={'monomials':rows}
    pl=v12.polynomial_coefficients(ledger,'PL')
    pr=v12.polynomial_coefficients(ledger,'PR')
    require(all(sp.expand(v12.polynomial_conjugate(x)-y)==0 for x,y in zip(pr,pl)),
            'PL/PR exact conjugacy failed')
    return pl

def compare_finite_and_poles(root, stem, previous, extract, mods):
    import sympy as sp
    pole,kernel,v12,v13,v14=mods
    model=MODELS[stem]
    f=root/BASE/'eft1_fermion_first_oneloop'/f'{stem}_final_eft1_fermion_first_oneloop_raw_hard_llhh_sector.inputform.txt'
    y=root/BASE/'eft1_scalar_Y_only_canonical_probe'/f'{stem}_final_y_only_canonical_match_raw_scalar_llhh.inputform.txt'
    x=root/BASE/'eft1_scalar_xy_fresh_kernel_v3'/f'{stem}_final_x_from_saved_y_match_x_followup_x_raw_llhh.inputform.txt'
    earlier=previous['models'][stem]
    require(f.is_file() and y.is_file() and x.is_file(),model+': missing raw boundary')
    require(digest(f)==earlier['F_boundary_input_sha256'], model+': F source changed')
    for stage, source in [('Y',y),('X',x)]:
        items=[z for z in extract['sector_results'] if
               z['model']==model and z['threshold']==stage and z['sector']=='LLHH']
        require(len(items)==1 and items[0]['source_sha256']==digest(source),
                model+': '+stage+' raw sector hash differs from v14 evidence')
    fsec=v13.split_sectors(f.read_text(encoding='utf-8-sig'))
    _,ysec=v14.unpack_sectors(y.read_text(encoding='utf-8-sig'),'Y',pole)
    _,xsec=v14.unpack_sectors(x.read_text(encoding='utf-8-sig'),'X',pole)
    expected={'F':(68,68,68),'Y':(50,50,32),'X':(18,18,0)}
    result={};polys={}
    for name,sec in [('F',fsec),('Y',ysec),('X',xsec)]:
        require(tuple(len(sec[k]) for k in ('pole','log','constant'))==expected[name],
                model+': '+name+' term inventory changed')
        poly={k:component_polynomials(sec[k],model,kernel,v12,len(sec[k]))
              for k in ('pole','log','constant')}
        require(all(sp.expand(a-b)==0 for a,b in zip(poly['pole'],poly['log'])),
                model+': '+name+' exact neutral pole/log mismatch')
        polys[name]=poly
        result[name]={'pole_log_exact_equal':True,
                      'terms':{k:len(sec[k]) for k in ('pole','log','constant')},
                      'PL_ordered_pole_3spurions':[str(p) for p in poly['pole']],
                      'PL_ordered_formal_constant_3spurions':[str(p) for p in poly['constant']],
                      'raw_sha256':digest({'F':f,'Y':y,'X':x}[name])}
    residue=tuple(sp.expand(sum(polys[name]['pole'][i] for name in ('F','Y','X')))
                  for i in range(3))
    require(all(t==0 for t in residue),model+': F+Y+X neutral pole sum did not vanish')
    def add3(*vectors):return tuple(sp.expand(sum(v[i] for v in vectors)) for i in range(3))
    constants=add3(*(polys[n]['constant'] for n in ('F','Y','X')))
    # The individual one-loop logarithms are L_a=ln(mu^2/M_a^2).
    # Since P_F+P_Y+P_X=0, rewriting about M_X is an exact identity:
    # K(mu)=K_F+K_Y+K_X+P_F ln(M_X^2/M_F^2)+P_Y ln(M_X^2/M_Y^2).
    ordered={'F':polys['F']['pole'],'Y':polys['Y']['pole'],
             'X':polys['X']['pole'],'constant':constants}
    result['exact_formal_sum']={
        'residue_sum_3spurions':[str(v) for v in residue],
        'zero_pole_sum_exact':True,
        'mass_log_reference':'M_X',
        'log_coefficients_without_hbar_over_MF':{
           'ln_MX2_over_MF2':[str(t) for t in ordered['F']],
           'ln_MX2_over_MY2':[str(t) for t in ordered['Y']]},
        'common_formal_constant_3spurions':[str(t) for t in constants],
        'mass_log_mu_dependence_cancelled_exact':True,
        'formula':'K_formal=K_F^const+K_Y^const+K_X^const+P_F*ln(MX^2/MF^2)+P_Y*ln(MX^2/MY^2); C5_formal=(hbar/MF)*L[[2K11,K12],[K12,2K22]]L^T',
        'is_renormalized_physical_matching':False}
    return result,ordered

def to_complex_matrix(values):
    return [[[float(z.real),float(z.imag)] for z in row] for row in values]

def numeric_evaluate(root,stem,record,poly,source_params,v12):
    """Explicitly FORMAL scale-log envelope; cannot return physical C5."""
    import numpy as np
    import sympy as sp
    model=MODELS[stem]
    require(record.get('model')==model,'Numeric input model mismatch')
    MF=record.get('MF_GeV')
    require(type(MF) in (float,int) and math.isfinite(MF) and MF>0,
            'MF_GeV must be finite positive')
    masses=source_params['heavy_masses_gev']
    require(len(masses)==2 and all(v>0 for v in masses) and masses[0]<masses[1],
            'Invalid ordered scalar masses')
    MX,MY=masses
    inputs=source_params['unitary_U']
    U=np.asarray(inputs['real'],dtype=float)+1j*np.asarray(inputs['imag'],dtype=float)
    require(U.shape==(2,2) and np.linalg.norm(U.conj().T@U-np.eye(2))<1e-9,
            'Invalid benchmark scalar mixing')
    yuk=[]
    for name in ('y1','y2'):
        vals=record.get(name)
        require(isinstance(vals,list) and len(vals)==3,name+' must have 3 entries')
        yuk.append(np.array([v12.parse_complex(v,name) for v in vals],complex))
    quartics=record.get('quartics')
    require(isinstance(quartics,dict) and set(quartics)==v12.QUARTICS,
            'Must supply exactly the three source quartics')
    lam={name:v12.parse_complex(value,name) for name,value in quartics.items()}
    substitutions={}
    for a in range(2):
        for b in range(2):
            key=f'T3MEU{a+1}{b+1}'
            substitutions[sp.Symbol(key)]=U[a,b]
            substitutions[sp.Symbol('bar_'+key)]=U[a,b].conjugate()
    for name,v in lam.items():
        substitutions[sp.Symbol(name)]=v
        substitutions[sp.Symbol('bar_'+name)]=v.conjugate()
    def evaluate(vec):
        out=[]
        for expr in vec:
            unbound=expr.free_symbols-set(substitutions)
            require(not unbound,'Unknown spurions in numeric matrix: '+str(unbound))
            out.append(complex(sp.N(expr.subs(substitutions),16)))
        return np.array(out,dtype=complex)
    logF=2*math.log(MX/MF);logY=2*math.log(MX/MY)
    k=evaluate(poly['constant'])+logF*evaluate(poly['F'])+logY*evaluate(poly['Y'])
    L=np.column_stack((yuk[0].conjugate(),yuk[1].conjugate()))
    K=np.array([[2*k[0],k[1]],[k[1],2*k[2]]],complex)
    C=(1/(16*math.pi**2*MF))*L@K@L.T
    require(np.max(np.abs(C-C.T))<=1e-12*max(1,float(np.max(np.abs(C)))),
            'Formal 3x3 matrix not complex symmetric')
    # Probe the same expression at another arbitrary common scale mu0:
    # must give identical result to relative machine precision.
    mu0=MY
    pF=evaluate(poly['F']);pY=evaluate(poly['Y']);pX=evaluate(poly['X'])
    kcheck=evaluate(poly['constant']) + 2*math.log(mu0/MF)*pF +\
           2*math.log(mu0/MY)*pY + 2*math.log(mu0/MX)*pX
    require(np.max(np.abs(kcheck-k))<=1e-10*max(1,float(np.max(np.abs(k)))),
            'Formal mass-log scale cancellation failed')
    return {'model':model,'status':'FORMAL_RAW_BOUNDARY_SUM_ONLY',
            'input_MF_GeV':MF,'input_MX_GeV':MX,'input_MY_GeV':MY,
            'formal_common_scale_C5_3x3_GeV_inverse_NOT_PHYSICAL':to_complex_matrix(C),
            'mu_scale_cancellation_numerically_checked':True,
            'finite_matching_counterterms_included':False,
            'physical_C5':None,'physical_neutrino_masses':None}

def self_test():
    import sympy as sp
    a,b,c,d,e=sp.symbols('a b c d e')
    pf=(a,b,c);py=(d,e,a);px=tuple(-x-y for x,y in zip(pf,py))
    require(all(sp.expand(sum(t))==0 for t in zip(pf,py,px)),
            'Toy zero-pole test')
    LF,LY,LX=sp.symbols('LF LY LX')
    for f,y,x in zip(pf,py,px):
        assert sp.expand(f*LF+y*LY+x*LX - f*(LF-LX)-y*(LY-LX))==0
    bad=tuple(x+1 for x in px)
    try:require(all(sp.expand(sum(t))==0 for t in zip(pf,py,bad)),
                'Correctly reject incorrect pole cancellation')
    except EnvelopeError:pass
    else:raise EnvelopeError('Tampered pole sum accepted')
    import numpy as np
    k=np.array([[1+1j,2-3j],[2-3j,4j]])
    L=np.array([[1,2],[0,1j],[3,0]],complex)
    C=L@k@L.T
    require(np.max(np.abs(C-C.T))<1e-13,'Transpose, not dagger, required')
    print('PASS: exact common-scale logarithm identity; wrong pole sum rejected')
    print('PASS: conventional symmetric Majorana 3x3 projection; no physical C5')

def run(root,out,numeric_input):
    import sympy as sp
    mods=import_deps(root)
    v14=read_json(root/AUD/'t3_bc_sequential_batch_v14/summary.json',
                 't3-bc-sequential-matching-fixed-order-batch-v14')
    v15=read_json(root/AUD/'t3_bc_one_loop_inheritance_v15/summary.json',
                 't3-bc-one-loop-fixed-order-inheritance-v15')
    extract=read_json(root/AUD/'t3_bc_scalar_finite_term_extraction_v1/summary.json',
                      't3-bc-native-local-matchete-pole-log-split-v1')
    require(v14.get('full_physical_matching_certified') is False and
            v15.get('full_physical_matching_certified') is False,
            'Prior physical certification flag changed')
    require(v15.get('Y_one_loop_LLHH_spectator_inheritance_at_Ohbar') is True,
            'Fixed-order inheritance not certified')
    results={};polydata={}
    for stem,model in MODELS.items():
        previous=v15['models'][stem]['F_Y_X_formal_pole_relation']
        require(previous['formal_neutral_residue_relations']['F_plus_Y_plus_X'] is True and
                previous['pole_sum_residual_nonzero_structure_count']==0,
                model+': v15 formal pole cancellation is not established')
        stages,p=compare_finite_and_poles(root,stem,v14,extract,mods)
        polydata[stem]=p
        results[stem]={'model':model,'stages':stages,
                       'source_EFT1_component_equivalence_status':
                         'Previously verified separately; Matchete-native equality not established',
                       'formal_mass_log_combination_is_physical_C5':False,
                       'scheme_finite_counterterms_available':False,
                       'full_1PI_high_minus_low_subtraction_done':False,
                       'physical_C5':None}
        print(model+': exact P_F+P_Y+P_X=0; formal common-scale finite envelope constructed')
    report={'format':'t3-bc-formal-finite-common-scale-envelope-v16',
            'status':'FormalScaleLogCancellationAndFiniteEnvelope_CompletePhysicalSubtractionMissing',
            'models':results,
            'mass_log_scale_cancellation_exact':True,
            'operator_normalization':'Warsaw PL positive sign; C5=A+A.T (v11.2)',
            'critical_distinction':'Raw F+Y+X amplitude sum is not the EFT high-minus-low 1PI matching difference; finite counterterms, EOM, field and parameter redefinitions are unverified',
            'required_next_calculation':'Explicit UV-minus-EFT one-loop 1PI matching subtraction with matching common renormalization scheme and finite CTs',
            'physical_matching_certified':False,'physical_C5':None,
            'physical_neutrino_mass_matrix':None}
    if numeric_input is not None:
        data=read_json(numeric_input)
        if all(k in MODELS.values() for k in data) and set(data)==set(MODELS.values()):
            records={stem:data[model] for stem,model in MODELS.items()}
        else:
            require(data.get('model') in MODELS.values(),'Unknown numerical model')
            stem=next(s for s,m in MODELS.items() if m==data['model'])
            records={stem:data}
        report['numeric_FORMAL_NONPHYSICAL']={}
        for stem,record in records.items():
            params=read_json(root/BASE/f'{stem}_final_parameters.json')
            report['numeric_FORMAL_NONPHYSICAL'][stem]=numeric_evaluate(
                root,stem,record,polydata[stem],params,mods[2])
        report['numeric_input_sha256']=digest(numeric_input)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('PASS: B/C scale-independent FORMAL envelope; renormalized physical C5 NOT CERTIFIED')
    print('Saved: '+str(out))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--root',type=Path,default=Path.cwd())
    ap.add_argument('--output',type=Path,default=Path('output/audits/t3_bc_formal_finite_envelope_v16/summary.json'))
    ap.add_argument('--numeric-input',type=Path,help='Optional v14-compatible one-model or B/C manifest JSON')
    ap.add_argument('--require-physical',action='store_true')
    args=ap.parse_args()
    try:
        if args.self_test:self_test();return 0
        if args.require_physical:
            raise EnvelopeError('Physical C5 cannot be certified: full EFT 1PI subtraction and finite renormalization missing')
        root=args.root.resolve()
        out=args.output if args.output.is_absolute() else root/args.output
        numeric=args.numeric_input
        if numeric is not None and not numeric.is_absolute():numeric=root/numeric
        run(root,out,numeric)
        return 0
    except (EnvelopeError,ValueError,TypeError,KeyError,IndexError,OSError,ImportError) as exc:
        print('FAIL (no physical C5): '+str(exc),file=sys.stderr)
        return 1

if __name__=='__main__':sys.exit(main())
