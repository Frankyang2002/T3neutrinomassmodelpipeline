"""Stage 22 SU(2) covariance audit for B_-1/C_-1 quartic block 110101.

Reads the native RGBeta invariant expressions exported by Wolfram, checks exact
infinitesimal SU(2) covariance under H -> U H and conjugate -> U* conjugate,
and tests a proposed replacement for one *reused* RGBeta invariant against
Matchete's native component tensors. This produces diagnostic evidence ONLY.
No running code, vertices, or physical coupling maps are changed/certified.

python -m scripts.Stage22RGBeta110101Covariance
"""
from __future__ import annotations

import argparse
from itertools import product
import hashlib
import json
from pathlib import Path

import sympy as sp

from scripts.Stage22RGBetaExplicitBasis import eval_component, parse_fields, verify_sources
from scripts.Stage22TensorAlgebra import exact, solve_basis, TensorCertificateError

MODELS = ('B_-1', 'C_-1')
COUNTS = (1, 1, 0, 1, 0, 1)
FIELDS = ('H', 'Hdag', 'S1dag', 'S2dag')
ORIGINAL_FIRST = 'lambdaT3GeneralQ110101I1'
ORIGINAL_REUSED = 'lambdaHHdagS1barS2barCross'
# The expected native expressions are checked explicitly. Do NOT silently
# accept a changed or non-matching source, or infer native source from a fix.
EXPECTED = {
    ORIGINAL_FIRST: 'del[SU2L[fund],t3i1,t3i4]*eps[SU2L[fund],t3i2,t3i3]',
    ORIGINAL_REUSED: 'eps[SU2L[fund],t3i1,t3i3]*eps[SU2L[fund],t3i2,t3i4]',
}
PROPOSED = 'del[SU2L[fund],t3i1,t3i3]*eps[SU2L[fund],t3i2,t3i4]'
_SIGMA = (
    sp.Matrix([[0, 1], [1, 0]]),
    sp.Matrix([[0, -sp.I], [sp.I, 0]]),
    sp.Matrix([[1, 0], [0, -1]]),
)
INDEX_GRID = tuple(product(range(2), repeat=4))


def components(expression: str) -> dict[tuple[int, ...], sp.Expr]:
    return {idx: sp.simplify(eval_component(expression, idx)) for idx in INDEX_GRID}


def generator_residual(tensor: dict[tuple[int,...], sp.Expr], generator: sp.Matrix):
    """Exact Lie-algebra variation of T_ijkl Phi^i (Phi*)^j (Phi*)^k (Phi*)^l.

    Each field index contracts its generator on the *right*: sum_j T_j G_ji.
    The three conjugate fields use -conjugate(G). All indices are native order.
    """
    gs = (generator, -sp.conjugate(generator), -sp.conjugate(generator), -sp.conjugate(generator))
    out = {}
    for idx in INDEX_GRID:
        res = sp.Integer(0)
        for axis, g in enumerate(gs):
            for b in range(2):
                other = list(idx)
                other[axis] = b
                res += g[b,idx[axis]] * tensor[tuple(other)]
        res = sp.simplify(res)
        if res != 0:
            out[idx] = res
    return out


def covariance(expression: str) -> dict:
    tensor = components(expression)
    witnesses = {}
    for i,sigma in enumerate(_SIGMA, start=1):
        res = generator_residual(tensor,sigma/2)
        witnesses[f'T{i}'] = {
            'nonzero_residual_components': len(res),
            'first_nonzero': [
                {'indices_1_based': [j+1 for j in idx], 'residual': str(value)}
                for idx,value in list(res.items())[:4]
            ],
        }
    return {
        'su2_covariant': all(x['nonzero_residual_components']==0 for x in witnesses.values()),
        'generator_checks':witnesses,
        'field_transformations': ['U', 'conjugate(U)', 'conjugate(U)', 'conjugate(U)'],
    }


def _native_tensor(native: dict, mat_ref: dict, expression: str):
    if tuple(parse_fields(native['native_fields_inputform'])) != FIELDS:
        raise ValueError('Native field order changed; do not assume covariance convention')
    if tuple(mat_ref['species']) != FIELDS or tuple(mat_ref['dimensions']) != (2,2,2,2):
        raise ValueError('Matchete field/component order changed; diagnostic not applicable')
    return {
        'symbol':native['symbol'], 'counts':list(COUNTS), 'dimensions':[2,2,2,2],
        'components':[str(components(expression)[idx]) for idx in INDEX_GRID],
    }


def analyze(matchete: dict, invariants: dict, declarations: dict) -> dict:
    verify_sources(matchete,invariants,declarations)
    result={
        'format':'stage22-su2-110101-native-covariance-diagnostic-v1',
        'status':'NativeExpressionSU2CovarianceAndProposedGeometryNotPhysicalCertification',
        'source':'Wolfram-exported native RGBeta invariant expressions and Matchete tensors',
        'component_convention':'Fundamental H transforms U; Bar[H], Bar[S1], Bar[S2] transform conjugate(U); epsilon[1,2]=+1',
        'native_wolfram_runtime_covariance_checked':False,
        'native_rgbeta_component_basis_verified':False,
        'potential_prefactors_verified':False,
        'physical_coupling_mapping_certified':False,
        'production_modified':False,
        'proposed_source_change':'Only lambdaHHdagS1barS2barCross invariant; requires independent RGBeta convention checks',
        'models':{},
    }
    for label in MODELS:
        mt = [t for t in matchete['models'][label]['tensors'] if tuple(t['counts'])==COUNTS]
        nt = [t for t in invariants['models'][label]['invariants']
              if t['symbol'] in EXPECTED]
        if len(mt)!=2 or len(nt)!=2 or {t['symbol'] for t in nt}!=set(EXPECTED):
            raise ValueError(f'{label}: 110101 block inventory unexpected')
        nt = sorted(nt,key=lambda t:[ORIGINAL_FIRST,ORIGINAL_REUSED].index(t['symbol']))
        for t in nt:
            current=t['native_invariant_expression_inputform']
            if ''.join(current.split()) != EXPECTED[t['symbol']]:
                raise ValueError(f'{label}/{t["symbol"]}: source invariant changed; inspect before trusting old diagnosis')
        c0,c1=(covariance(t['native_invariant_expression_inputform']) for t in nt)
        cp=covariance(PROPOSED)
        if not c0['su2_covariant'] or c1['su2_covariant'] or not cp['su2_covariant']:
            raise ValueError(f'{label}: native SU2 covariance result differs from verified case')
        original_records=[_native_tensor(t,mt[0],t['native_invariant_expression_inputform']) for t in nt]
        proposed_records=[original_records[0],_native_tensor(nt[1],mt[0],PROPOSED)]
        # Existing pair cannot be compared as an invertible basis without
        # changing a species' index convention; candidate MUST match exactly.
        old_M=sp.Matrix.hstack(*(sp.Matrix([exact(v) for v in x['components']]) for x in mt))
        old_R=sp.Matrix.hstack(*(sp.Matrix([exact(v) for v in x['components']]) for x in original_records))
        old_column_in_mat_span=[old_M.row_join(old_R[:,i]).rank()==old_M.rank() for i in range(2)]
        if old_column_in_mat_span != [True,False]:
            raise ValueError(f'{label}: old native geometry status changed')
        new_map=solve_basis(mt,proposed_records)
        expect=[['1/2','1'],['-sqrt(3)/2','0']]
        A=sp.Matrix([[exact(x) for x in row] for row in new_map['r_in_m']])
        E=sp.Matrix([[exact(x) for x in row] for row in expect])
        if any(sp.simplify(x)!=0 for x in A-E):
            raise ValueError(f'{label}: candidate changed geometry; verify independently')
        if sp.simplify(A.det()-sp.sqrt(3)/2)!=0:
            raise ValueError(f'{label}: candidate determinant changed')
        result['models'][label]={
            'block':'110101', 'matchete_symbols':[t['symbol'] for t in mt],
            'registered_rgbeta_symbols':[t['symbol'] for t in nt],
            'native_first':c0,
            'native_reused':c1,
            'proposed_replacement':{'expression_inputform':PROPOSED, **cp},
            'original_rgbeta_columns_in_matchete_span':old_column_in_mat_span,
            'proposed_exact_geometry':new_map,
            'proposed_matrix_determinant':str(sp.simplify(A.det())),
            'physical_mapping_certified':False,
        }
    result['counts']={'models_checked':2,'native_symbols_checked':4,
                      'native_pass':2,'native_fail':2,'proposed_pass':2,
                      'proposed_full_rank_geometry_blocks':2}
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,default=Path('output/general_t3_stage22'))
    parser.add_argument('--output',type=Path,default=None)
    args=parser.parse_args()
    files=['matchete_tensors.json','rgbeta_native_invariants.json','rgbeta_native_declarations.json']
    values=[json.loads((args.input_dir/f).read_text(encoding='utf-8-sig')) for f in files]
    report=analyze(*values)
    report['input_sha256']={f:hashlib.sha256((args.input_dir/f).read_bytes()).hexdigest() for f in files}
    output=args.output or args.input_dir/'stage22_110101_su2_covariance.json'
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Stage 22 native-expression SU2 covariance:', report['counts'])
    print('Proposed 110101 basis transformation:',report['models']['B_-1']['proposed_exact_geometry']['r_in_m'])
    print('Physical coupling mapping: UNCERTIFIED; production unchanged')
    print('Report:',output)


if __name__=='__main__':main()
