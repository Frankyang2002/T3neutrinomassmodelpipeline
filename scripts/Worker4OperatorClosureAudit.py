"""Read-only EFT1 dimension-five closure audit; does not solve or alter production IVP.

Run from repository root: python scripts/Worker4OperatorClosureAudit.py
Optionally supply --boundary certified-or-synthetic-master-boundary.json
and --allow-synthetic, to probe the actual seven-term tensor at the boundary.
"""
from __future__ import annotations
import argparse
import ast
import csv
import hashlib
import importlib
import json
from pathlib import Path
import sys
import numpy as np

SOURCE = Path('RGE/running/intermediate/GeneralEFT1C5Mixing.py')
BRIDGE = Path('Numerical/running/T3QuarticWilsonBridge.py')
JOINT = Path('Numerical/running/T3ExtendedSimultaneousEFT1.py')
FULL = Path('Numerical/running/Aminus4EFT1Couplings.py')
TERMS = ('yukawa_wavefunction','scalar_pair','mixed_yukawa_gauge','crossed_yukawa',
         'conjugate_coefficient_yukawa','scalar_anomalous_dimension','fermion_anomalous_dimension')
OTHER_CLASSES = [
 ('psi2phi3','4.84','dimension-six coefficient times dimension-one cubic/mass; requires checking EFT power counting'),
 ('phi2F2','4.86','dimension-six bosonic class with two fermion-mass insertions'),
 ('phiF2','4.87','dimension-five scalar-gauge class with fermion mass/Yukawa contraction'),
 ('psi2phiF','4.88','dimension-six fermion-scalar-field-strength class with mass insertion'),
 ('psi2F','4.89','dimension-five fermion dipole class; gauge/Yukawa selection required'),
 ('psi2psi2','4.90','dimension-six four-fermion class with scalar trilinear'),
 ('psi4','4.91','dimension-six four-fermion class with scalar trilinear'),
 ('D2phi4','4.92','dimension-six derivative-scalar class; mass insertions'),
 ('Dpsipsiphi2','4.93','dimension-six derivative-fermion class'),
]


def _read(path):
    if not path.is_file():
        raise FileNotFoundError('Required project file missing: '+str(path))
    return path.read_text(encoding='utf-8')


def source_inventory():
    src = _read(SOURCE)
    tree = ast.parse(src)
    funcs = {n.name:n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    if not {'calculate_master_beta','_collinear','load_boundary'} <= set(funcs):
        raise RuntimeError('Actual seven-term source API changed; re-audit required')
    calc = ast.get_source_segment(src,funcs['calculate_master_beta'])
    names = [x for x in TERMS if ('dst["'+x+'"]') in calc]
    if tuple(names) != TERMS:
        raise RuntimeError('Seven-term implementation changed; audit refuses to certify')
    bridge = _read(BRIDGE)
    joint = _read(JOINT)
    full = _read(FULL)
    from Numerical.running.GeneralT3OrderedTensor import scalar_basis
    from Numerical.running.Aminus4EFT1Couplings import F_LABELS, S_LABELS
    if len(F_LABELS)!=45 or len(S_LABELS)!=12:
        raise RuntimeError('Full SM or A(-4) basis cardinality changed')
    model_ids=[f'{family}_{alpha:+d}' for family,alphas in {
        'A':(-4,-2,0),'B':(-3,-1,1),'C':(-3,-1,1),
        'D':(-2,0,2),'E':(-4,-2,0,2)}.items() for alpha in alphas]
    scalars={m:{'n_real':scalar_basis(m).total_dimension,
                'blocks':dict(scalar_basis(m).dimensions),
                'real_field_names':sorted(scalar_basis(m).real_fields)}
             for m in model_ids}
    return {
        'scalar_inventory_by_model':scalars,
        'full_45_weyl_indexed_labels_0based':list(enumerate(F_LABELS)),
        'full_A_minus4_scalar_labels_0based':list(enumerate(S_LABELS)),
        'source_path':str(SOURCE).replace('\\','/'),
        'source_sha256':hashlib.sha256(src.encode()).hexdigest(),
        'source_term_names':names,
        'source_uses_complex_conjugation':'cc = c.conj()' in calc,
        'source_dimensions_externally_specified': 'n_fermions' in src and 'n_real_scalars' in src,
        'explicit_quark_higgs_gamma_hook': 'extra_scalar_anomalous_dimension' in bridge,
        'quark_gamma_added_to_scalar_anomalous_dimension': "terms['scalar_anomalous_dimension']" in bridge,
        'simultaneous_solver_unchanged': 'solve_ivp' in joint and 'running_seven_term_beta' in joint,
        'full_45_weyl_dictionary_exists': 'assert len(F_LABELS)==45' in full,
        'full_weyl_species': ['Q','uc','dc','L','ec'],
        'nine_weyl_leptons':['L[1,1]','L[1,2]','L[2,1]','L[2,2]','L[3,1]','L[3,2]',
                              'ec[1]','ec[2]','ec[3]'],
        'full_A_minus4_real_scalar_basis':'H doublet (4 real), S1 singlet (2 real), S2 triplet (6 real)',
        'generic_psi2phi2_self_mixing_exhaustive':False,
        'source_only_implements_psi2phi2_self_mixing':True,
    }


def _write_json(p,obj):
    p.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')


def _witnesses():
    # Construct a minimal algebraic example, not an SM- or T3-gauge-certified model.
    from types import SimpleNamespace
    from RGE.running.intermediate.GeneralEFT1C5Mixing import calculate_master_beta, _collinear
    nf, ns = 3, 3
    y=np.zeros((nf,nf,ns),complex)
    q=np.zeros((ns,)*4,float)
    c=np.zeros((nf,nf,ns,ns),complex)
    # Retained observable C[0,0,0,0]; omitted same lepton pair C[0,0,1,1].
    c[0,0,1,1]=1
    for inds in {(0,0,1,1),(0,0,1,1),(0,1,0,1),(1,0,0,1),(1,1,0,0),(0,1,1,0),(1,0,1,0)}:
        q[inds]=1
    ctx=SimpleNamespace(nf=nf,ns=ns,yukawa=y,quartic=q,gauge=())
    a=calculate_master_beta(ctx,c)
    direct={k:float(abs(v[0,0,0,0])) for k,v in a.items()}
    # An omitted fermion slot can enter gamma_f and the external retained coefficient
    # only when a renormalizable flavour/representation mixing is allowed.
    c2=np.zeros_like(c);c2[0,2,0,0]=c2[2,0,0,0]=1
    y[0,2,2]=y[2,0,2]=1; y[2,2,2]=1
    b=calculate_master_beta(ctx,c2)
    return {'classification':'algebraic_support_witness_only_not_a_gauge_invariant_T3_point',
            'omitted_scalar_pair_into_retained':direct,
            'omitted_fermion_index_into_retained':{k:float(abs(v[0,0,0,0])) for k,v in b.items()},
            'zero_yukawa_scalar_gamma_in_first_example':bool(np.allclose(_collinear(SimpleNamespace(ns=ns,nf=nf,yukawa=np.zeros_like(y),gauge=()))[0],0))}


def _boundary_audit(filename,allow_synthetic):
    from RGE.running.intermediate.GeneralEFT1C5Mixing import load_boundary,calculate_master_beta,TERMS as live_terms
    if tuple(live_terms)!=TERMS: raise RuntimeError('Source terms changed')
    payload=json.loads(Path(filename).read_text(encoding='utf-8'))
    ctx=load_boundary(payload,allow_synthetic=allow_synthetic)
    n=ctx.nf; ns=ctx.ns
    retained=set(i for pair in ctx.leptons for i in pair)
    h=set(ctx.higgs)
    rows=[]
    for name,c in ctx.operators:
        positions=np.argwhere(abs(c)>0)
        for ijab in positions:
            i,j,a,b=map(int,ijab)
            if i in retained and j in retained and a not in h and b not in h:
                label='retained_LLSS'
            elif i in retained and j in retained and a in h and b in h:
                label='retained_LLHH'
            else:label='other_psi2phi2_coordinate'
            rows.append({'source':name,'i':i+1,'j':j+1,'a':a+1,'b':b+1,'class':label})
    # Probe beta using the *actual* coefficients, keeping source terms separate.
    terms={name:0.0 for name in TERMS}
    nu=[pair[0] for pair in ctx.leptons]
    rh=ctx.higgs[2]
    for _,c in ctx.operators:
        for term,arr in calculate_master_beta(ctx,c).items():
            terms[term]+=float(np.linalg.norm(arr[np.ix_(nu,nu,[rh],[rh])]))
    return {'model':ctx.model,'authority':ctx.authority,'n_fermions':n,'n_real_scalars':ns,
            'lepton_indices_1based':[[i+1 for i in p] for p in ctx.leptons],
            'higgs_indices_1based':[i+1 for i in ctx.higgs],
            'coordinate_rows':rows,'per_term_LLHH_norm_at_boundary':terms,
            'closure_proved':False,
            'warning':'No omitted-dimension-five coefficient inventory is supplied by this schema; cannot infer absence or closure.'}


def run(root:Path,boundary=None,allow_synthetic=False):
    root.mkdir(parents=True,exist_ok=True)
    inv=source_inventory()
    _write_json(root/'source_inventory.json',inv)
    wit=_witnesses()
    _write_json(root/'algebraic_contraction_witnesses.json',wit)
    with (root/'mixing_blocks.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['target','source','reference_equation','order_and_condition'])
        w.writerow(['psi2phi2','psi2phi2','4.85','direct one-loop self-class mixing; seven terms implemented'])
        for source,eq,desc in OTHER_CLASSES:w.writerow(['psi2phi2',source,eq,desc])
    if boundary:
        b=_boundary_audit(boundary,allow_synthetic)
        positions=b.pop('coordinate_rows')
        _write_json(root/'boundary_audit.json',b)
        with (root/'boundary_coordinates.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=['source','i','j','a','b','class']);w.writeheader();w.writerows(positions)
    else:
        _write_json(root/'boundary_audit.json',{'status':'NO_CERTIFIED_BOUNDARY_SUPPLIED','closure_proved':False})
    _write_json(root/'verdict.json',{
        'classification':'UNRESOLVED_FULL_CLOSURE; CONDITIONALLY_SUPPORTED_FIXED_ONE_LOOP_TRUNCATION',
        'strict_closure_proved':False,
        'key_conditions':['No omitted tree-level dimension-five boundary sources feeding retained targets',
                          'No nonzero symmetry-allowed omitted-to-retained same-order anomalous-dimension block',
                          'Dimension-six sources and mass insertions consistently excluded by stated 1/Lambda power counting',
                          'All heavy fermions integrated out before scalars, with complete matching inventory',
                          'No leptoquark or baryon/lepton violating renormalizable couplings unless explicitly audited'],
        'order_statement':'C_O(MF)=0 and C_O first generated at one-loop => feedback into C_R first at two-loop leading-log; C_O(MF) nonzero => one-loop feedback possible.',
        'source_guarantee':'Seven self-mixing terms alone do not prove closure against omitted psi2F, phiF2 or other classes.'})
    return root


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',default='output/audits/worker4_operator_closure')
    p.add_argument('--boundary',type=Path)
    p.add_argument('--allow-synthetic',action='store_true')
    a=p.parse_args()
    root=run(Path(a.output),a.boundary,a.allow_synthetic)
    print('Worker4 audit:',root,'\nVerdict: unresolved strict closure; conditional fixed-one-loop truncation')

if __name__=='__main__':main()
