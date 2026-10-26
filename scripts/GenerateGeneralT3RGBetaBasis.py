"""Generate candidate *complete* SU(2) extra-quartic bases for 16 ordinary T3 models.

The output is a standalone RGBeta extension, NOT a certified Matchete-to-RGBeta
normalisation map. It uses RGBeta's native eps, del, delS2 and S2 triple
invariants. Existing six B/C(-1) vertices are not duplicated.

No existing production Wolfram/Python file is edited by this generator.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import permutations, product
from math import prod
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'docs/general_scalar_trial_operator_manifest.json'
OUTPUT = ROOT / 'RGE/general/generated'
SPECIES = ('H', 'Hdag', 'S1', 'S1dag', 'S2', 'S2dag')
FIELD = {'H': 'H', 'Hdag': 'Bar @ H', 'S1': 'S1', 'S1dag': 'Bar @ S1',
         'S2': 'S2', 'S2dag': 'Bar @ S2'}
REP = {'A': (1, 3, 2), 'B': (2, 2, 1), 'C': (2, 2, 3),
       'D': (3, 1, 2), 'E': (3, 3, 2)}
# The 6 ordinary B/C(-1) RGBeta couplings are already in T3RGBetaModel.wl.
# One of the two H Hdag S1dag S2dag invariants remains missing.
BC_EXISTING = {
    (0, 0, 0, 1, 1, 2): 'lambdaS1barS2S2bar2',
    (0, 0, 0, 2, 0, 2): 'lambdaS1bar2S2bar2',
    (0, 0, 1, 2, 0, 1): 'lambdaS1S1bar2S2bar',
    (0, 2, 0, 0, 2, 0): 'lambdaHHdagS2S2',
    (0, 2, 0, 2, 0, 0): 'lambdaHHdagS1barS1bar',
    (1, 1, 0, 1, 0, 1): 'lambdaHHdagS1barS2barCross',
}

Pairing = tuple[tuple[int, int], ...]

def pairings(indices: tuple[int, ...]) -> list[Pairing]:
    if not indices:
        return [()]
    first = indices[0]
    out = []
    for partner in indices[1:]:
        rest = tuple(i for i in indices if i not in (first, partner))
        for tail in pairings(rest):
            out.append(tuple(sorted(((first, partner), *tail))))
    return sorted(set(out))


def symmetry_perms(names: list[str]) -> list[tuple[int, ...]]:
    groups = defaultdict(list)
    for i, name in enumerate(names):
        groups[name].append(i)
    perm_list = [tuple(range(len(names)))]
    for positions in groups.values():
        if len(positions) == 1:
            continue
        new = []
        for perm in perm_list:
            for swap in permutations(positions):
                nxt = list(perm)
                for p, q in zip(positions, swap):
                    nxt[p] = perm[q]
                new.append(tuple(nxt))
        perm_list = new
    return sorted(set(perm_list))


def sym_pairing(pairing: Pairing, perms: list[tuple[int, ...]], names: list[str], rep: int) -> dict[Pairing, Fraction]:
    result: dict[Pairing, Fraction] = defaultdict(Fraction)
    for perm in perms:
        edges = tuple(sorted(tuple(sorted((perm[a], perm[b]))) for a,b in pairing))
        sign = (-1) ** sum(perm[a] > perm[b] and rep == 2 and (names[a].endswith("dag") == names[b].endswith("dag")) for a,b in pairing)
        result[edges] += Fraction(sign, len(perms))
    return {k: v for k, v in result.items() if v}


def signed_pair(a: int, b: int, names: list[str], rep: int) -> str:
    if rep == 3:
        return f'del[SU2L @ S2, #{a+1}, #{b+1}]'
    if rep != 2:
        raise ValueError('Pairing only works for doublets or triplets.')
    first_bar = names[a].endswith('dag')
    second_bar = names[b].endswith('dag')
    if first_bar == second_bar:
        return f'eps[SU2L @ fund, #{a+1}, #{b+1}]'
    if first_bar:
        return f'del[SU2L @ fund, #{a+1}, #{b+1}]'
    return f'del[SU2L @ fund, #{b+1}, #{a+1}]'


def numeric_pair(a: int, b: int, indices: tuple[int, ...], names: list[str], rep: int) -> int:
    ia, ib = indices[a], indices[b]
    if rep == 3 or (names[a].endswith('dag') != names[b].endswith('dag')):
        return int(ia == ib)
    return {(0,1): 1, (1,0): -1}.get((ia,ib),0)


def numeric_vector(weights: dict[Pairing, Fraction], names: list[str], rep: int) -> tuple[Fraction, ...]:
    return tuple(
        sum((coef * prod(numeric_pair(a,b,indices,names,rep) for a,b in pairs)
             for pairs,coef in weights.items()), Fraction())
        for indices in product(range(rep), repeat=len(names))
    )


def independent_weights(names: list[str], rep: int) -> list[dict[Pairing, Fraction]]:
    from sympy import Matrix
    perms = symmetry_perms(names)
    options = [sym_pairing(p,perms,names,rep) for p in pairings(tuple(range(len(names))))]
    selected, cols = [], []
    rank = 0
    for candidate in options:
        vec = numeric_vector(candidate, names, rep)
        if not any(vec):
            continue
        attempt = cols + [vec]
        new_rank = Matrix(attempt).rank()
        if new_rank > rank:
            cols,rank = attempt,new_rank
            selected.append(candidate)
    return selected


def pairing_expr(weights: dict[Pairing, Fraction], names: list[str], rep: int) -> str:
    pieces = []
    for edges, coefficient in sorted(weights.items()):
        factors = [signed_pair(a,b,names,rep) for a,b in edges]
        item = ' * '.join(factors) if factors else '1'
        if coefficient != 1:
            item = f'({coefficient.numerator}/{coefficient.denominator}) * ({item})'
        pieces.append(item)
    return '('+' + '.join(pieces)+')'


def invariants(names: list[str], dims: list[int]) -> list[str]:
    active = [i for i,d in enumerate(dims) if d != 1]
    types = [dims[i] for i in active]
    if not types:
        return ['1']
    if len(active)==2 and types == [2,2]:
        a,b=active
        if names[a]==names[b]:
            raise ValueError('Antisymmetric doublet singlet vanishes for identical bosons')
        return [signed_pair(a,b,names,2)]
    if len(active)==2 and types == [3,3]:
        a,b=active
        return [signed_pair(a,b,names,3)]
    if len(active)==3 and sorted(types)==[2,2,3]:
        t=next(i for i in active if dims[i]==3)
        a,b=(i for i in active if dims[i]==2)
        return [f'delS2[SU2L, #{t+1}, #{a+1}, #{b+1}]']
    if len(active)==4 and len(set(types))==1:
        return [pairing_expr(w,names,types[0]) for w in independent_weights(names,types[0])]
    if len(active)==4 and sorted(types)==[2,2,3,3]:
        db = [i for i in active if dims[i]==2]
        tr = [i for i in active if dims[i]==3]
        result = []
        if names[db[0]] != names[db[1]] or names[db[0]].endswith('dag')!=names[db[1]].endswith('dag'):
            result.append('(' + signed_pair(db[0],db[1],names,2) + ' * '
                          +signed_pair(tr[0],tr[1],names,3) + ')')
        if names[tr[0]] != names[tr[1]]:
            result.append('(delS2[SU2L, a, #'+str(db[0]+1)+', #'+str(db[1]+1)+'] * '
                          +'T3RGBetaS2TripleInvariant[a, #'+str(tr[0]+1)+', #'+str(tr[1]+1)+'])')
        return result
    raise ValueError(f'Unsupported SU(2) factor content {names}: {types}')


def make_model(model: str, data: dict) -> dict[str,Any]:
    group, suffix = model.split('_', 1)
    d1,d2,df = REP[group]
    alpha=int(suffix)
    emitted = []
    reused = []
    for row in data['extra_quartic_structures']:
        counts=tuple(map(int,row['counts']))
        names=[name for name,n in zip(SPECIES,counts) for _ in range(n)]
        dims=[2 if name.startswith('H') else d1 if name.startswith('S1') else d2 for name in names]
        invs = invariants(names,dims)
        expected = int(row['invariants'])
        if len(invs)!=expected:
            raise AssertionError(f'{model} {counts} invariant rank {len(invs)} != {expected}')
        existing = BC_EXISTING.get(counts) if (group in ('B','C') and alpha == -1) else None
        if existing:
            # Existing basis coefficient is not identified with a particular Matchete I#.
            # One SU(2) structure has two singlets: keep existing and add one independent.
            if expected==1:
                reused.append({'counts': list(counts), 'existing_name': existing})
                continue
            if expected!=2:
                raise AssertionError('Unexpected B/C existing multiplicity')
            # The new contraction pairs H with the conjugate S2 via delta,
            # and the two conjugate doublets Hdag/S1dag via epsilon.
            # Epsilon(H, S2dag) would not be SU(2)-invariant: these indices
            # transform in mutually dual representations.
            # The existing RGBeta mixed quartic is only *structurally* reused;
            # its exact tensor orientation still needs independent validation.
            invs = [
                '(del[SU2L @ fund, #4, #1] * eps[SU2L @ fund, #2, #3])'
            ]
            reused.append({'counts':list(counts),'existing_name':existing})
        for idx, inv in enumerate(invs, start=1):
            digit=''.join(map(str,counts))
            name=f'lambdaT3GeneralQ{digit}I{idx}'
            emitted.append({'name':name,'counts':list(counts), 'fields':[FIELD[s] for s in names],
                            'invariant':inv,'matchete_count':expected,
                            'matchete_basis_name':f'lambdaT3ExtraQ{digit}I#',
                            'mapping_status':'NotNormalised'})
    if len(emitted)+len(reused)!=int(data['extra_quartic_complex_couplings']):
        raise AssertionError(f'{model}: incomplete emitted+reused count')
    return {'model':model,'dimensions':[d1,d2,df], 'alpha':alpha,
            'total_additional_quartics':int(data['extra_quartic_complex_couplings']),
            'new_rgbeta_symbols':emitted,'preexisting_rgbeta_symbols':reused,
            'quadratics_not_implemented':data['quadratic_complex_couplings_not_enabled']}


def wl_model(r:dict[str,Any]) -> str:
    lines = [f'(* AUTO-GENERATED RGBeta extension for {r["model"]}.',
             '   Invariant coverage checked against the independent Matchete count.',
             '   Do not identify these coupling values with Matchete couplings:',
             '   tensor basis orientation and normalisation remain uncertified.',
             '   Do not activate native RGBeta gauge option: v1.2.0 on',
             '   Windows rejects it at AddQuartic option validation.',
             '   Gauge-structure checks are in the independent Python audit. *)',
             'T3GeneralQuarticAddedSymbols = {};',
             'T3GeneralQuarticExpectedCount = '+str(r['total_additional_quartics'])+';',
             'T3GeneralQuarticReusedCount = '+str(len(r['preexisting_rgbeta_symbols']))+';',
             '']
    for v in r['new_rgbeta_symbols']:
        name=v['name']
        lines += [f'AddQuartic[{name},',
                  '  {'+', '.join(v['fields'])+'},',
                  '  GroupInvariant -> ('+v['invariant']+' &),',
                  '  SelfConjugate -> False',
                  '];',
                  f'AppendTo[T3GeneralQuarticAddedSymbols, {name}];','']
    lines += ['If[Length[T3GeneralQuarticAddedSymbols] + T3GeneralQuarticReusedCount =!= T3GeneralQuarticExpectedCount,',
              '  Print["ERROR: General RGBeta quartic coverage mismatch"]; Abort[]',
              '];','']
    return '\n'.join(lines)



def audit_candidate_gauge_structures(registry:dict[str,Any]) -> dict[str,Any]:
    """Check U(1)_Y, Z2, field-typing and invariant-rank prerequisites.

    This is NOT RGBeta's internal CheckInvariance, nor a numerical SU(2)
    Ward-identity verification of the compiled native tensor. Keeping those
    distinct prevents an accepted RGBeta beta export being called a proof of
    gauge invariance or a Matchete-normalisation certificate.
    """
    checked = 0
    for model, entry in registry['models'].items():
        group, suffix = model.split('_',1)
        alpha = int(suffix)
        if (entry['dimensions'],entry['alpha']) != (list(REP[group]),alpha):
            raise ValueError(f'{model}: inconsistent dimensions/alpha')
        charge_twice = {
            'H': 1, 'Hdag': -1,
            'S1': alpha, 'S1dag': -alpha,
            'S2': alpha+2, 'S2dag': -(alpha+2),
        }
        for v in entry['new_rgbeta_symbols']:
            fields = [species for species, count in zip(SPECIES,v['counts']) for _ in range(count)]
            if len(fields) != 4 or sum(charge_twice[f] for f in fields) != 0:
                raise ValueError(f'{model} {v["name"]}: nonzero hypercharge or nonquartic')
            if sum(f.startswith(('S1','S2')) for f in fields)%2:
                raise ValueError(f'{model} {v["name"]}: Z2 odd')
            if [FIELD[f] for f in fields] != v['fields']:
                raise ValueError(f'{model} {v["name"]}: incorrect field orientation')
            dims = [2 if f.startswith('H') else REP[group][0] if f.startswith('S1') else REP[group][1] for f in fields]
            candidates = invariants(fields,dims)
            if len(candidates) != v['matchete_count']:
                raise ValueError(f'{model} {v["name"]}: rank mismatch')
            # The special mixed B/C contraction uses an additional valid
            # delta*epsilon tensor, not a mixed-index epsilon.
            if model in ('B_-1','C_-1') and v['counts']==[1,1,0,1,0,1]:
                required='(del[SU2L @ fund, #4, #1] * eps[SU2L @ fund, #2, #3])'
                if v['invariant'] != required:
                    raise ValueError(f'{model}: invalid mixed doublet index pairing')
            elif v['invariant'] not in candidates:
                raise ValueError(f'{model} {v["name"]}: candidate invariant absent from rank basis')
            checked += 1
    if checked != registry['totals']['rgbeta_new_vertices']:
        raise ValueError('Gauge candidate audit did not check all generated quartics')
    return {
        'status':'AnalyticChargeParityFieldTypingAndInvariantRankVerified',
        'candidate_quartics_checked':checked,
        'hypercharge_zero':True,'z2_even':True,
        'su2_invariant_candidate_rank_checked':True,
        'native_rgbeta_checkinvariance_used':False,
        'native_rgbeta_su2_ward_identity_verified':False,
        'matchete_rgbeta_tensor_normalisation_verified':False,
    }


def generate(manifest:dict[str, Any], destination:Path=OUTPUT) -> dict[str,Any]:
    if manifest.get('model_count')!=16 or len(manifest.get('models',{}))!=16:
        raise ValueError('Refusing incomplete or non-16-model manifest')
    destination.mkdir(parents=True,exist_ok=True)
    result = {}
    for model,body in sorted(manifest['models'].items()):
        info = make_model(model,body)
        result[model] = info
        (destination/f'{model}_extra_quartics.wl').write_text(wl_model(info),encoding='utf-8')
    totals = {'matchete_extra_quartics':sum(d['total_additional_quartics'] for d in result.values()),
              'rgbeta_preexisting_structural':sum(len(d['preexisting_rgbeta_symbols']) for d in result.values()),
              'rgbeta_new_vertices':sum(len(d['new_rgbeta_symbols']) for d in result.values()),
              'quadratic_couplings_blocked':sum(d['quadratics_not_implemented'] for d in result.values())}
    if totals != {'matchete_extra_quartics':50,'rgbeta_preexisting_structural':12,
                   'rgbeta_new_vertices':38,'quadratic_couplings_blocked':8}:
        raise AssertionError(f'Unexpected coverage totals {totals}')
    report={'status':'GeneratedCandidateBasisNotWolframVerified', 'manifest_source_commit':manifest.get('source_commit'),
            'scope':'16 ordinary T3 models; experimental RGBeta quartic interactions',
            'totals':totals,'models':result,
            'caveats':['SU(2) candidate rank passes independent checks; native RGBeta Ward identities and Matchete tensor normalisation remain unverified.',
                       'Existing B/C couplings are counted once; their Matchete coupling maps are not established globally.',
                       '37? additional structures cannot be assumed to share Matchete invariant basis coefficients.',
                       'Eight extra complex quadratic mass operators require a separate mass-matrix extension.',
                       'No complete running of the new coupling set is yet validated.']}
    # Avoid editorial placeholder in persisted machine report.
    report['caveats'][2]='Exact basis changes, factors, phases and normalisations need an independent tensor conversion.'
    (destination/'general_rgbeta_quartic_registry.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return report


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest',type=Path,default=MANIFEST)
    ap.add_argument('--output-dir',type=Path,default=OUTPUT)
    opts=ap.parse_args()
    data=json.loads(opts.manifest.read_text(encoding='utf-8'))
    result=generate(data,opts.output_dir)
    print(json.dumps(result['totals'],indent=2))
    print('Registry:',opts.output_dir/'general_rgbeta_quartic_registry.json')
    return 0

if __name__=='__main__':raise SystemExit(main())
