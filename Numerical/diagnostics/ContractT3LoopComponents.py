"""Explore T3 component contractions without assigning physical multiplicity.

Input: Reports/output/full/multiplicity/component_support_audit.json.
An SU(2) index is not necessarily a charge eigenstate (notably real triplets).
The available CG audit omits the propagator invariant tensors; therefore this
module compares explicitly labelled candidate index maps, and never promotes
one to a verified physical contraction or compensation factor.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sympy as sp

DEFAULT_INPUT = Path('Reports/output/full/multiplicity/component_support_audit.json')
DEFAULT_OUTPUT = Path('Reports/output/full/multiplicity/contraction_candidates.json')


def expression(raw: str) -> sp.Expr:
    return sp.sympify(raw.replace('^','**'), locals={'I':sp.I,'sqrt':sp.sqrt}, evaluate=True)


def terms(vertex: dict, prefix: tuple[int, ...] = ()) -> list[tuple[tuple[int,...],sp.Expr]]:
    return [(tuple(c['indices']),expression(c['cg'])) for c in vertex['components']
            if tuple(c['indices'][:len(prefix)]) == prefix]


def pairing_map(dimension: int, kind: str) -> dict[tuple[int,int],sp.Expr]:
    """Diagnostic candidate propagator-index pairings, not physical assertions.

    identity: array-index equality; reversed: anti-diagonal unit mapping;
    epsilon: antisymmetric pseudoreal 2x2 metric. Only use epsilon for d=2.
    """
    if kind == 'identity':
        return {(i,i):sp.Integer(1) for i in range(1,dimension+1)}
    if kind == 'reversed':
        return {(i,dimension+1-i):sp.Integer(1) for i in range(1,dimension+1)}
    if kind == 'epsilon':
        if dimension != 2:
            raise ValueError('epsilon metric is only defined here for dimension 2')
        return {(1,2):sp.Integer(1),(2,1):-sp.Integer(1)}
    raise ValueError(kind)


def candidates(case: dict, *, fermion_metric: str = 'identity',
               scalar1_metric: str = 'identity', scalar2_metric: str = 'identity',
               quartic_conjugated: bool = False) -> dict:
    model=case['model']; v=case['vertices']
    d1,d2,df=(int(model[x]) for x in ('dS1','dS2','dF'))
    maps=[pairing_map(d,k) for d,k in ((df,fermion_metric),(d1,scalar1_metric),(d2,scalar2_metric))]
    # External lepton 1 = neutral lepton in conventional (nu,e) ordering.
    # External H=2 denotes neutral Higgs in (H+,H0) ordering.
    ys1=terms(v['T3Y1CG'],(1,)); ys2=terms(v['T3Y2CG'],(1,))
    mix=terms(v['T3MixCG'],(2,2))
    out=[]
    for (l1,w1) in ys1:
        for (l2,w2) in ys2:
            fm=maps[0].get((l1[1],l2[1]))
            if fm is None: continue
            for (q,wq) in mix:
                sm1=maps[1].get((l1[2],q[2]))
                sm2=maps[2].get((l2[2],q[3]))
                if sm1 is None or sm2 is None: continue
                weight=sp.simplify(w1*w2*(sp.conjugate(wq) if quartic_conjugated else wq)*fm*sm1*sm2)
                if weight==0:continue
                out.append({'y1_indices':list(l1),'y2_indices':list(l2),
                            'mix_indices':list(q),'weight':sp.sstr(weight),
                            'weight_abs_squared':sp.sstr(sp.simplify(weight*sp.conjugate(weight)))})
    total=sp.simplify(sum((expression(x['weight']) for x in out),sp.S.Zero))
    return {'candidate_only':True,'model':model,'index_metrics':{
                'fermion':fermion_metric,'scalar1':scalar1_metric,'scalar2':scalar2_metric,
                'quartic_conjugated':quartic_conjugated},
            'external_indices_assumed':{'lepton':1,'Higgs':2},
            'channel_count':len(out),'complex_sum':sp.sstr(total),
            'channels':out, 'loop_multiplicity':None,'compensation_factor':None}


def analyze(payload: dict) -> dict:
    rows=[]
    for case in payload['cases']:
        df=int(case['model']['dF'])
        possibilities=['identity','reversed']+(['epsilon'] if df==2 else [])
        variants=[]
        for fm in possibilities:
            variants.append(candidates(case,fermion_metric=fm))
        rows.append({'model':case['model'],'variants':variants,
                     'verified_physical_multiplicity':None,
                     'reason':'Propagator index metrics, conjugation, and charge-eigenstate transformations are not supplied by the CG support audit.'})
    return {'status':'candidate_contractions_only','count':len(rows),
            'models':rows,'physics_warning':
            'These results are NOT physical neutrino-mass channel counts. In particular, real triplet indices are not charge eigenstates, and charge selection cannot be inferred from nonzero CG entries. Do not rescale lambdaT3 from these counts.'}


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',type=Path,default=DEFAULT_INPUT)
    ap.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    a=ap.parse_args()
    if not a.input.exists():
        ap.error(f'Missing component support audit: {a.input}')
    result=analyze(json.loads(a.input.read_text(encoding='utf-8-sig')))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(f"Models: {result['count']}; contraction status: {result['status']}")
    for r in result['models']:
        m=r['model']; desc=' '.join(f"{v['index_metrics']['fermion']}={v['channel_count']}" for v in r['variants'])
        print(f"{m['dS1']},{m['dS2']},{m['dF']},alpha={m['alpha']}: {desc}")
    print('Physical multiplicities: not assigned. No coupling changes.')
    print(f'JSON: {a.output.resolve()}')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
