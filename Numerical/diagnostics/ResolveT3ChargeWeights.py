"""Resolve T3 charge-eigenstate channel coefficients from verified Matchete CGs.

Inputs: component_support_audit.json, matchete_generator_analysis.json,
        cg_covariance_audit.json, charge_channels.json in Reports/output/full/multiplicity.
Output: charge_resolved_loop_weights.json in the same directory.

Computes the SU(2) coefficient of Y1 * Y2 * conjugate(Mixing), where
Y1 = Bar[L] F^c S1, Y2 = Bar[L] F Bar[S2], Mixing = H H S1 Bar[S2].
The conjugate Mixing vertex connects scalar propagators. Fermion F^c and F
carry opposing electric charges, with equal underlying F weight index.

Does NOT fix overall sign, identical-Higgs factors, Majorana/Dirac masses,
loop-function differences, or normalization conventions for physical Yukawas.
"""
from __future__ import annotations
import argparse
import itertools
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / 'Reports/output/full/multiplicity'
FLAGS = {'T3Y1CG': (True, True, False),
         'T3Y2CG': (True, False, True),
         'T3MixCG': (False, False, False, True)}
CLASSES = {(1,3,2):'A',(2,2,1):'B',(2,2,3):'C',(3,1,2):'D',(3,3,2):'E'}

def exact(value):
    return sp.sympify(value, locals={'I':sp.I})

def clean(value):
    return sp.simplify(sp.expand(value))

def matrix(data):
    return sp.Matrix([[exact(x) for x in row] for row in data])

def rotations(analysis):
    us = {1:sp.eye(1)}
    for entry in analysis['representations']:
        if entry['status'] != 'generator_algebra_verified':
            raise ValueError('Generator verification failed')
        d = entry['dimension']; u = matrix(entry['U_columns_are_weight_eigenvectors'])
        if clean(u.H*u-sp.eye(d)) != sp.zeros(d):
            raise ValueError(f'Nonunitary transformation for dimension {d}')
        us[d] = u
    if set(us) != {1,2,3}:
        raise ValueError('Required 1,2,3 representations missing')
    return us

def transformed(vertex, flags, us):
    dims = vertex['dimensions']
    assert len(dims) == len(flags)
    # field_old[i] = sum_a U[i,a] field_new[a] (U* for conjugate field).
    transforms = [us[d].conjugate() if conj else us[d] for d,conj in zip(dims,flags)]
    result = {}
    for component in vertex['components']:
        old = [x-1 for x in component['indices']]
        value = exact(component['cg'])
        for new in itertools.product(*(range(d) for d in dims)):
            coefficient = value
            for axis in range(len(dims)):
                coefficient *= transforms[axis][old[axis],new[axis]]
            if coefficient != 0:
                result[new] = result.get(new,sp.S.Zero)+coefficient
    return {idx:clean(val) for idx,val in result.items() if clean(val) != 0}

def charge_of(d,y,k):
    return clean(sp.Rational(d-1,2)-k+y)

def evaluate(case, us, expected_charges):
    model=case['model']; ds1=model['dS1'];ds2=model['dS2'];df=model['dF'];alpha=model['alpha']
    hyper={'F':sp.Rational(alpha+1,2), 'S1':sp.Rational(alpha,2),
           'S2':sp.Rational(alpha+2,2)}
    cg={name:transformed(case['vertices'][name],FLAGS[name],us) for name in FLAGS}
    # external neutral L^0 is weight +1/2 (index 0); neutral H^0 is -1/2 (index 1).
    channels=[]
    for f in range(df):
        q=charge_of(df,hyper['F'],f)
        for s1 in range(ds1):
            if charge_of(ds1,hyper['S1'],s1)!=q: continue
            for s2 in range(ds2):
                if charge_of(ds2,hyper['S2'],s2)!=q: continue
                y1=cg['T3Y1CG'].get((0,f,s1),sp.S.Zero)
                y2=cg['T3Y2CG'].get((0,f,s2),sp.S.Zero)
                mixing=cg['T3MixCG'].get((1,1,s1,s2),sp.S.Zero)
                # The quartic's Hermitian conjugate is required for scalar propagators.
                weight=clean(y1*y2*sp.conjugate(mixing))
                channels.append({'charge':str(q),'indices_weight_basis_one_based':{'F':f+1,'S1':s1+1,'S2':s2+1},
                                 'y1':str(y1),'y2':str(y2),'mixing':str(mixing),
                                 'complex_weight':str(weight),'nonzero':weight!=0})
    got={x['charge'] for x in channels}
    if got!=set(expected_charges):
        raise ValueError(f'Charge enumeration mismatch for {model}: {got} vs {expected_charges}')
    active=[ch for ch in channels if ch['nonzero']]
    weights=[exact(ch['complex_weight']) for ch in active]
    equal=bool(weights) and all(clean(w-weights[0])==0 for w in weights)
    return {'model':model,'class':CLASSES[(ds1,ds2,df)],'charge_allowed_count':len(channels),
            'nonzero_weight_count':len(active),'channels':channels,
            'total_group_coefficient':str(clean(sum(weights,sp.S.Zero))),
            'equal_nonzero_weights':equal,
            'multiplicity_if_degenerate_and_common_loop_integral':len(active) if equal else None,
            'verified_physical_loop_multiplicity':None,'compensation_factor':None}

def run(components,analysis,covariance,charges):
    if covariance.get('status')!='cg_generator_covariance_verified' or not all(x['all_covariant'] for x in covariance['models']):
        raise ValueError('CG covariance precondition not verified')
    if components['errors'] or components['cases'] is None:
        raise ValueError('Component support input invalid')
    us=rotations(analysis)
    lookup={(m['model']['dS1'],m['model']['dS2'],m['model']['dF'],m['model']['alpha']):m for m in charges['models']}
    results=[]
    for case in components['cases']:
        md=case['model'];key=tuple(md[k] for k in ('dS1','dS2','dF','alpha'))
        results.append(evaluate(case,us,lookup[key]['internal_charge_channels']))
    return {'status':'charge_resolved_group_coefficients', 'count':len(results), 'models':results,
            'contraction':'Y1 * Y2 * conjugate(T3MixCG) using same physical internal weight indices',
            'assumptions':['Scalar propagators connect fields to their Hermitian conjugates',
                           'F^c and F are paired by the same F weight index',
                           'Mass eigenstates do not mix different electric charges'],
            'warning':'Equal CG weights alone do not establish identical loop functions, spinor factors or physical multiplicity. Do not rescale lambdaT3 from this report without that check.'}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,default=DEFAULT_OUT)
    args=parser.parse_args(); base=args.directory
    load=lambda filename:json.loads((base/filename).read_text(encoding='utf-8'))
    result=run(load('component_support_audit.json'),load('matchete_generator_analysis.json'),
               load('cg_covariance_audit.json'),load('charge_channels.json'))
    path=base/'charge_resolved_loop_weights.json'
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Models:',result['count'])
    for rec in result['models']:
        print(rec['class'],rec['model']['alpha'], 'weights:',[(c['charge'],c['complex_weight']) for c in rec['channels']], 'total:',rec['total_group_coefficient'])
    print('Report:',path)

if __name__=='__main__':main()
