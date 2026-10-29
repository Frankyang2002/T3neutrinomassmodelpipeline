"""Stage 22 independent diagnostic: resolve triplet CG and selected pseudoreal
index identifications against actual Matchete and RGBeta invariant exports.

This is NOT a physical mapping or a native RGBeta component export. In particular,
per-block doublet reindexings have not been justified as one global field basis.
No production files are changed and no certification status is promoted.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
import re

import sympy as sp

from scripts.RGBetaExplicitBasis import (
    SOURCE_ORDER, eval_component, parse_fields, reconstructed_tensor, verify_sources,
)
from scripts.RGBetaTensorAlgebra import TensorCertificateError, exact, solve_basis

CALIBRATION = 'A_+0/021010'
TWIST_BLOCKS = frozenset({'000112', '000202', '001201', '020200'})
TWIST_MODELS = frozenset({'B_-1', 'C_-1'})
BAD_BLOCK = '110101'
FUNDAMENTAL_REINDEX = ('S1', 'S1dag')
CG_SCALE = sp.sqrt(2) / 3**sp.Rational(1,4)
ROOT2 = sp.sqrt(2)

# Explicit unit-normalised symmetric doublet x doublet -> triplet tensors,
# calibrated from Matchete A_+0 and then used to predict five OTHER cases.
# C[a,i,j]; no claim this is RGBeta's native component coordinate choice.
CG = (
    sp.Matrix([[1, 0], [0, -1]]) / ROOT2,
    sp.Matrix([[0, -1], [-1, 0]]) / ROOT2,
    -sp.I * sp.eye(2) / ROOT2,
)
_EPS = sp.Matrix([[0, 1], [-1, 0]])


def _plain(x: sp.Expr) -> str:
    return str(sp.simplify(x))


def _fields_and_dims(native, mat):
    fields = parse_fields(native['native_fields_inputform'])
    if fields != mat['species']:
        raise ValueError('RGBeta/Matchete species ordering mismatch')
    return fields, mat['dimensions']


def triplet_orthonormality():
    return all(sp.simplify(sum(sp.conjugate(CG[a][i,j])*CG[b][i,j]
                                for i,j in product(range(2),repeat=2)) - (1 if a==b else 0)) == 0
               for a,b in product(range(3),repeat=2))


@lru_cache(maxsize=1)
def triple_cg():
    """Native T3RGBetaS2TripleInvariant contracted from symmetric C and epsilon."""
    out = {}
    for a,b,c in product(range(3), repeat=3):
        out[a,b,c] = sp.simplify(sum(
            CG[a][i,j]*CG[b][k,l]*CG[c][m,n] *
            _EPS[i,k]*_EPS[j,m]*_EPS[l,n]
            for i,j,k,l,m,n in product(range(2), repeat=6)
        ))
    return out


def _e_native_pattern():
    """Canonicalised exact InputForm for the native E-model triple contraction.

    The $digits names arise from Mathematica Module and differ between runs.
    Refuse any change of the actual index pairings or numerical prefactors.
    """
    return (
        '((del[SU2L[fund], i$G, t3i2]*del[SU2L[fund], j$G, t3i1])/2 + '
        '(del[SU2L[fund], i$G, t3i1]*del[SU2L[fund], j$G, t3i2])/2)*'
        'delS2[SU2L, t3i3, k$G, l$G]*delS2[SU2L, t3i4, m$G, n$G]*'
        'eps[SU2L[fund], i$G, k$G]*eps[SU2L[fund], j$G, m$G]*'
        'eps[SU2L[fund], l$G, n$G]'
    )


def triplet_tensor(native, mat):
    fields, dims = _fields_and_dims(native,mat)
    expr = native['native_invariant_expression_inputform']
    simple = re.fullmatch(r'delS2\[SU2L, t3i([34]), t3i1, t3i2\]', expr)
    mixed = re.sub(r'\$\d+', '$G', expr) == _e_native_pattern()
    if not simple and not mixed:
        raise ValueError('Unknown native triplet contraction structure; refuse to infer')
    if simple:
        t = int(simple.group(1))-1
        if dims[t] != 3 or dims[0:2] != [2,2]:
            raise ValueError('Invalid SU2 symmetric triplet axes')
        vals = [CG[idx[t]][idx[0],idx[1]] for idx in product(*(range(d) for d in dims))]
        construction = 'Unit-normalised symmetric-doublet CG; A_+0 calibrated'
    else:
        if tuple(dims) != (2,2,3,3):
            raise ValueError('Native triple-triplet field dimensions changed')
        g3 = triple_cg()
        vals = [sp.simplify(sum(CG[a][idx[0],idx[1]]*g3[a,idx[2],idx[3]]
                                for a in range(3)))
                for idx in product(*(range(d) for d in dims))]
        construction = 'Predicted with calibrated CG and source T3RGBetaS2TripleInvariant epsilon contractions'
    return {
        'symbol': native['symbol'], 'species': fields,
        'counts': mat['counts'], 'dimensions': dims,
        'components': [_plain(v) for v in vals]
    }, construction


def reindexed_doublet_tensor(native, mat, rotated=FUNDAMENTAL_REINDEX):
    fields, dims = _fields_and_dims(native,mat)
    if any(d != 2 for d in dims):
        raise ValueError('Pseudoreal reindexing requires four fundamental doublets')
    if any(f in rotated and d != 2 for f,d in zip(fields,dims)):
        raise ValueError('Attempting doublet reindexing on a non-doublet field')
    expr = native['native_invariant_expression_inputform']
    components = []
    for idx in product(*(range(d) for d in dims)):
        in_idx = tuple(1-i if fields[j] in rotated else i for j,i in enumerate(idx))
        sign = (-1)**sum(i for f,i in zip(fields,idx) if f in rotated)
        components.append(_plain(sign * eval_component(expr,in_idx)))
    return {
        'symbol':native['symbol'], 'species':fields,
        'counts':mat['counts'], 'dimensions':dims,
        'components':components,
    }


def _span_witness(mats,invs):
    """Enumerate ALL 2^4 choices of I/J on the active 110101 species."""
    species = mats[0]['species']
    unique = tuple(dict.fromkeys(species))
    if len(unique)!=4 or len(invs)!=2:
        raise ValueError('Obstruction witness expects four distinct doublets/two RGBeta slots')
    M=sp.Matrix.hstack(*(sp.Matrix([exact(x) for x in t['components']]) for t in mats))
    independent_per_native={t['symbol']:[] for t in invs}
    common=[]
    for cardinal in range(5):
        for group in combinations(unique,cardinal):
            cols=[];good=[]
            for t in invs:
                candidate = reindexed_doublet_tensor(t,mats[0],group)
                v=sp.Matrix([exact(x) for x in candidate['components']])
                fits = M.row_join(v).rank()==M.rank()
                if fits:
                    independent_per_native[t['symbol']].append(list(group))
                good.append(fits)
                cols.append(candidate)
            if all(good):
                try:
                    ans=solve_basis(mats,cols)
                    common.append({'rotation':list(group),'basis_change':ans['r_in_m']})
                except TensorCertificateError:
                    pass
    return {
        'rotation_choices_exhausted':2**len(unique),
        'native_coupling_rotation_witnesses':independent_per_native,
        'common_full_rank_rotations':common,
        'consistent_I_or_J_alignment_exists':bool(common),
    }


def _by_blocks(entries, is_mat):
    out=defaultdict(list)
    for e in entries:
        if is_mat:
            counts=tuple(e['counts'])
        else:
            fs=parse_fields(e['native_fields_inputform'])
            counts=tuple(fs.count(s) for s in SOURCE_ORDER)
        out[counts].append(e)
    return out


def analyze(mats, invariants, declarations):
    verify_sources(mats,invariants,declarations)
    if not triplet_orthonormality():
        raise ValueError('Chosen CG basis is not unit-normalised')
    report={
        'format':'stage22-aligned-geometry-diagnostics-v1',
        'status':'ExactConditionalGeometryNotPhysicalCouplingCertification',
        'native_rgbeta_component_basis_verified':False,
        'physical_prefactors_verified':False,
        'physical_coupling_map_certified':False,
        'matchete_source_calibration':CALIBRATION,
        'triplet_symmetric_cg_component_basis':[ [[_plain(CG[a][i,j]) for j in range(2)] for i in range(2)] for a in range(3)],
        'triplet_calibration_scale':_plain(CG_SCALE),
        'doublet_reindexing':'Historical B/C J candidate is not a verified global basis; Stage24 replacements must match directly',
        'models':{}, 'totals':defaultdict(int),
    }
    for model in sorted(mats['models']):
        mb=_by_blocks(mats['models'][model]['tensors'],True)
        rb=_by_blocks(invariants['models'][model]['invariants'],False)
        if mb.keys()!=rb.keys():raise ValueError(f'{model}: block contents disagree')
        blocks=[]
        for ct,mat_block in sorted(mb.items()):
            code=''.join(map(str,ct)); inv_block=rb[ct]
            if len(mat_block)!=len(inv_block): raise ValueError(f'{model}/{code}: basis multiplicity mismatch')
            base={'counts':list(ct),'block':code,
                  'matchete_symbols':[e['symbol'] for e in mat_block],
                  'rgbeta_symbols':[e['symbol'] for e in inv_block],
                  'physical_mapping_certified':False}
            try:
                if model in TWIST_MODELS and code==BAD_BLOCK:
                    if [t['symbol'] for t in inv_block] != ['lambdaT3GeneralQ110101I1', 'lambdaHHdagS1barS2barCross']:
                        raise ValueError('Unexpected 110101 native symbol order')
                    compact = [''.join(t['native_invariant_expression_inputform'].split()) for t in inv_block]
                    fixed_first = 'del[SU2L[fund],t3i1,t3i4]*eps[SU2L[fund],t3i2,t3i3]'
                    old_second = 'eps[SU2L[fund],t3i1,t3i3]*eps[SU2L[fund],t3i2,t3i4]'
                    fixed_second = 'del[SU2L[fund],t3i1,t3i3]*eps[SU2L[fund],t3i2,t3i4]'
                    if compact[0] != fixed_first or compact[1] not in (old_second,fixed_second):
                        raise ValueError('110101 source changed outside reviewed contraction; refuse inference')
                    if compact[1] == fixed_second:
                        expanded=[reconstructed_tensor(t,mat_block[0]) for t in inv_block]
                        ans=solve_basis(mat_block,expanded)
                        if ans['r_in_m'] != [['1/2','1'],['-sqrt(3)/2','0']]:
                            raise ValueError('Corrected 110101 tensor geometry changed')
                        base.update({'status':'CorrectedReusedInvariantCandidate',
                                     'geometric':ans})
                    else:
                        witness=_span_witness(mat_block,inv_block)
                        if witness['consistent_I_or_J_alignment_exists']:
                            raise ValueError('Historical 110101 unexpectedly has a common alignment')
                        base.update({'status':'IncompatibleNativePair',
                                     'reason':'Historical reused epsilon-epsilon invariant has no common I/J alignment.',
                                     'exhaustive_witness':witness})
                elif model in TWIST_MODELS and code in TWIST_BLOCKS:
                    # Historical source used delS2 x delS2, necessitating a
                    # provisional per-block J rotation.  Stage24 repairs the
                    # four registered GroupInvariants.  The corrected native
                    # tensors must now match Matchete directly, with no J.
                    if len(inv_block) != 1:
                        raise ValueError('Unexpected B/C reused quartet multiplicity')
                    expr = inv_block[0]['native_invariant_expression_inputform']
                    if 'delS2[' in expr:
                        expanded=[reindexed_doublet_tensor(t,mat_block[0]) for t in inv_block]
                        ans=solve_basis(mat_block,expanded)
                        base.update({'status':'PseudorealReindexCandidate',
                                     'component_rotation_applied':list(FUNDAMENTAL_REINDEX),
                                     'geometric':ans})
                    else:
                        expanded=[reconstructed_tensor(t,mat_block[0]) for t in inv_block]
                        ans=solve_basis(mat_block,expanded)
                        if ans['r_in_m'] != [['1']]:
                            raise ValueError('Stage24 B/C corrected native tensor must match Matchete exactly')
                        base.update({'status':'CorrectedBCReusedInvariantCandidate',
                                     'geometric':ans})
                elif any('delS2[' in t['native_invariant_expression_inputform'] for t in inv_block):
                    tensors, sources=zip(*(triplet_tensor(t,mat_block[0]) for t in inv_block))
                    ans=solve_basis(mat_block,tensors)
                    base.update({'status':'CalibratedTripletCGCandidate',
                                 'construction':list(sources), 'geometric':ans})
                else:
                    expanded=[reconstructed_tensor(t,mat_block[0]) for t in inv_block]
                    ans=solve_basis(mat_block,expanded)
                    base.update({'status':'DirectGeometryCandidate','geometric':ans})
            except (TensorCertificateError,ValueError) as exc:
                if 'status' in base and base['status']=='IncompatibleNativePair':raise
                raise ValueError(f'{model}/{code}: unexpected geometric failure: {exc}') from exc
            blocks.append(base)
            report['totals'][base['status']]+=len(mat_block)
            report['totals']['slots']+=len(mat_block)
            report['totals']['blocks']+=1
        report['models'][model]=blocks
    report['totals']=dict(report['totals'])
    core={'DirectGeometryCandidate':32,'PseudorealReindexCandidate':8,
          'CalibratedTripletCGCandidate':6,'slots':50,'blocks':46}
    expected_before={**core,'IncompatibleNativePair':4}
    expected_after={**core,'CorrectedReusedInvariantCandidate':4}
    expected_stage24={'DirectGeometryCandidate':32,
        'CorrectedBCReusedInvariantCandidate':8,
        'CalibratedTripletCGCandidate':6,
        'CorrectedReusedInvariantCandidate':4,
        'slots':50,'blocks':46}
    if report['totals'] not in (expected_before,expected_after,expected_stage24):
        raise ValueError(f'Unexpected totals {report["totals"]!r}; expected a reviewed whole-model source phase')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,default=Path('output/general_t3_stage22'))
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    names=('matchete_tensors.json','rgbeta_native_invariants.json','rgbeta_native_declarations.json')
    def load(name):return json.loads((args.input_dir/name).read_text(encoding='utf-8-sig'))
    report=analyze(*(load(n) for n in names))
    report['input_sha256']={n:hashlib.sha256((args.input_dir/n).read_bytes()).hexdigest() for n in names}
    dest=args.output or args.input_dir/'stage22_aligned_geometry_diagnostics.json'
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Stage 22 conditional aligned geometry:',report['totals'])
    print('Physical coupling normalisation: UNCERTIFIED')
    print('Report:',dest)


if __name__=='__main__':main()
