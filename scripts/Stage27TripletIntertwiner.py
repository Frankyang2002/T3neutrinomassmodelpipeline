"""Stage 27: exact six-case Matchete triplet CG / SU(2) generator witness.

This is a *global triplet reference-basis candidate*, not a measured native
RGBeta Cartesian component map. Native RGBeta delS2 is an abstract invariant;
Stage 25 did not export a component-valued S2 generator. Consequently this
stage cannot certify the physical Matchete--RGBeta coupling map.

Run (after Stage 26): python -m scripts.Stage27TripletIntertwiner
"""
from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
import re
from pathlib import Path

import sympy as sp

from scripts.Stage22RGBetaExplicitBasis import parse_fields
from scripts.Stage22RGBetaAlignedGeometry import _e_native_pattern
from scripts.Stage22TensorAlgebra import exact
from scripts.Stage26QuarticActionIdentity import run as stage26_run

DATA = Path('output/general_t3_stage22')
SIMPLE = ('A_+0', 'A_-2', 'D_+0', 'D_-2')
DOUBLE = ('E_+0', 'E_-2')
ALL = SIMPLE + DOUBLE
EPS = sp.Matrix([[0, 1], [-1, 0]])
REFERENCE = (
    sp.Matrix([[1, 0], [0, 0]]),
    sp.Matrix([[0, 1], [1, 0]]) / sp.sqrt(2),
    sp.Matrix([[0, 0], [0, 1]]),
)
PAULI = (
    sp.Matrix([[0, 1], [1, 0]]) / 2,
    sp.Matrix([[0, -sp.I], [sp.I, 0]]) / 2,
    sp.Matrix([[1, 0], [0, -1]]) / 2,
)
# SU(2) generator order is (sigma1, sigma2, sigma3)/2.


def fail(message):
    raise ValueError('Stage27 fail-closed: ' + message)


def zero(expr):
    if isinstance(expr, sp.MatrixBase):
        return all(sp.simplify(z) == 0 for z in expr)
    return sp.simplify(expr) == 0


def ser(expr):
    return str(sp.simplify(expr))


def tensor_for(model, matchete, invariants):
    if model not in matchete.get('models', {}) or model not in invariants.get('models', {}):
        fail('Absent triplet model ' + model)
    mat = matchete['models'][model]['tensors']
    rgb = invariants['models'][model]['invariants']
    candidates = [x for x in rgb if 'delS2[' in x['native_invariant_expression_inputform']]
    if len(candidates) != 1:
        fail('Expected precisely one native triplet CG quartic for ' + model)
    invariant = candidates[0]
    expression = invariant['native_invariant_expression_inputform']
    if model in SIMPLE:
        if not re.fullmatch(r'delS2\[SU2L, t3i[34], t3i1, t3i2\]', expression):
            fail('Unrecognised simple native triplet expression: ' + model)
    elif re.sub(r'\$\d+', '$G', expression) != _e_native_pattern():
        fail('Unrecognised source E-series triple-triplet expression: ' + model)
    species = parse_fields(invariant['native_fields_inputform'])
    matching = [x for x in mat if x['species'] == species]
    if len(matching) != 1:
        fail('Expected precisely one Matchete triplet CG quartic for ' + model)
    record = matching[0]
    dims = record['dimensions']
    if len(dims) != 4 or dims[:2] != [2, 2] or (dims[2:] not in ([3, 1], [1, 3], [3, 3])):
        fail('Unrecognised triplet quartic dimensions for ' + model)
    expected = sp.prod(dims)
    if len(record['components']) != expected:
        fail('Malformed Matchete triplet vector ' + model)
    return invariant, record


def entry(record, index):
    offset = 0
    for k, d in zip(index, record['dimensions']):
        offset = offset * d + k
    return exact(record['components'][offset])


def simple_cg(record):
    dims = record['dimensions']
    if dims[:2] != [2, 2] or sorted(dims[2:]) != [1, 3]:
        fail('Not a one-triplet calibration tensor')
    trip_axis = 2 if dims[2] == 3 else 3
    out = []
    for a in range(3):
        c = sp.Matrix(2, 2, lambda i, j:
            entry(record, (i, j, a, 0) if trip_axis == 2 else (i, j, 0, a)))
        if not zero(c - c.T):
            fail('Native Matchete triplet CG is not symmetric')
        out.append(c)
    norms = [sp.simplify(sum(sp.conjugate(c[i, j]) * c[i, j]
                             for i, j in product(range(2), repeat=2))) for c in out]
    if any(not zero(n - norms[0]) for n in norms) or zero(norms[0]):
        fail('Native Matchete triplet CG norms disagree or vanish')
    # Norm is >0 for this source: fixed exact positive root, no arbitrary phases.
    magnitude = sp.sqrt(norms[0]); normalized = tuple(sp.simplify(c / magnitude) for c in out)
    gram = sp.Matrix([[sp.simplify(sum(sp.conjugate(normalized[a][i,j]) * normalized[b][i,j]
                                   for i,j in product(range(2),repeat=2)))
                       for b in range(3)] for a in range(3)])
    if not zero(gram - sp.eye(3)):
        fail('Native Matchete calibration CG is not orthonormal')
    return normalized, magnitude


def triple_contraction(cg):
    """RGBeta's sourced T3RGBetaS2TripleInvariant epsilon contractions.

    The result uses the *inferred symmetric-fundamental reference basis* cg,
    NOT native RGBeta generator-component measurements.
    """
    return {(a,b,c): sp.simplify(sum(
        cg[a][i,j]*cg[b][k,l]*cg[c][m,n]
        * EPS[i,k]*EPS[j,m]*EPS[l,n]
        for i,j,k,l,m,n in product(range(2), repeat=6)))
        for a,b,c in product(range(3),repeat=3)}


def predicted_components(record, cg, simple_case):
    dims = record['dimensions']; t = []
    g3 = triple_contraction(cg) if not simple_case else None
    for idx in product(*(range(d) for d in dims)):
        if simple_case:
            t.append(cg[idx[2] if dims[2]==3 else idx[3]][idx[0],idx[1]])
        else:
            if dims != [2, 2, 3, 3]:
                fail('E-series triplet quartic dimensions changed')
            t.append(sp.simplify(sum(cg[a][idx[0],idx[1]] * g3[a,idx[2],idx[3]]
                                     for a in range(3))))
    return t


def fit_source_scale(record, predicted):
    physical = [exact(x) for x in record['components']]
    if len(physical) != len(predicted):
        fail('Incorrect triplet predictor dimensions')
    indices = [i for i, p in enumerate(predicted) if not zero(p)]
    if not indices or zero(physical[indices[0]]):
        fail('Zero or missing native Matchete CG component')
    # P = k*M: k is a source component-tensor ratio, not a physical coupling map.
    scale = sp.simplify(predicted[indices[0]] / physical[indices[0]])
    if any(not zero(v-scale*m) for v,m in zip(predicted, physical)):
        fail('A source Matchete triplet tensor disagrees with one common CG')
    return scale


def generators(cg):
    """Sym^2 representation of the standard Pauli generators.

    Two barred external doublets transform with -conjugate(t). The complex
    conjugate field orientation is kept explicit: this does NOT assert native
    RGBeta S2 generator components.
    """
    U = sp.Matrix([[sp.simplify(sum(cg[a][i,j]*sp.conjugate(REFERENCE[b][i,j])
                  for i,j in product(range(2),repeat=2))) for b in range(3)] for a in range(3)])
    if not zero(U*U.H-sp.eye(3)):
        fail('Derived triplet change of orthonormal reference basis not unitary')
    if zero(sp.det(U)):
        fail('Derived triplet basis is singular')
    ref_generators=[]
    for t in PAULI:
        tbar = -sp.conjugate(t)
        gram=sp.Matrix(3,3,lambda a,b:sp.simplify(sum(
            sp.conjugate(REFERENCE[a][i,j])
            *(tbar*REFERENCE[b]+REFERENCE[b]*tbar.T)[i,j]
            for i,j in product(range(2),repeat=2))))
        ref_generators.append(gram)
    mat_generators=[sp.simplify(U*g*U.H) for g in ref_generators]
    for gens in (ref_generators, mat_generators):
        if any(not zero(g-g.H) for g in gens):
            fail('Spin-one generator fails Hermiticity')
        if not zero(sum((g*g for g in gens), sp.zeros(3))-2*sp.eye(3)):
            fail('Spin-one Casimir is not 2')
        if any(not zero(gens[i]*gens[j]-gens[j]*gens[i]-sp.I*gens[k])
               for i,j,k in ((0,1,2),(1,2,0),(2,0,1))):
            fail('SU(2) generator commutator is nonzero')
    if any(not zero(mat_generators[i]*U-U*ref_generators[i]) for i in range(3)):
        fail('Global triplet intertwiner residual nonzero')
    return U, ref_generators, mat_generators


def audit(mat, native):
    if mat.get('model_count') != 16 or mat.get('quartic_slot_count') != 50 or native.get('total_slots') != 50:
        fail('Incomplete Matchete/RGBeta source tensor inventory')
    if len(mat.get('models', {})) != 16 or set(mat['models']) != set(native.get('models', {})):
        fail('Incoherent 16-model triplet evidence')
    records = {name:tensor_for(name, mat, native) for name in ALL}
    c0, magnitude = simple_cg(records[SIMPLE[0]][1])
    expected_source_norm = sp.sqrt(2 / sp.sqrt(3))
    if not zero(magnitude-expected_source_norm):
        fail('Matchete native CG calibration normalisation changed')
    for label in SIMPLE[1:]:
        c, mag = simple_cg(records[label][1])
        if not zero(mag-magnitude) or any(not zero(c[i]-c0[i]) for i in range(3)):
            fail('No common Matchete triplet CG across the A/D branches')
    # Check every source record without selecting a new CG for E. A *single*
    # extracted Matchete CG must predict all six contractions and H.c. cases.
    witnesses=[]
    for label in ALL:
        native_record, mat_record = records[label]
        pred = predicted_components(mat_record,c0,label in SIMPLE)
        scale = fit_source_scale(mat_record,pred)
        expected = sp.root(3, 4)/sp.sqrt(2) if label in SIMPLE else 1/sp.sqrt(2)
        if not zero(scale-expected):
            fail('Native Matchete tensor calibration ratio changed: ' + label)
        witnesses.append({'model':label,'matchete_symbol':mat_record['symbol'],
                          'rgbeta_symbol':native_record['symbol'],
                          'dimensions':mat_record['dimensions'],
                          'reference_cg_to_matchete_component_tensor_ratio':ser(scale),
                          'all_components_exact':True,
                          'same_triplet_cg_no_model_specific_rephasing':True})
    U, tref, tmat = generators(c0)
    # U U^T is a basis-conjugation diagnostic, *not* the physical barred-field rule.
    K = sp.simplify(U*U.T)
    return {
       'format':'stage27-triplet-generator-intertwiner-candidate-v1',
       'status':'ExactSingleMatcheteCGSixTripletsGlobalRGBetaBasisUncertified',
       'source_native_matchete_tensor_count':6,
       'models_triplet_checked':list(ALL),
       'matchete_reference_cg_normalization':ser(magnitude),
       'one_global_triplet_intertwiner_candidate':True,
       'triplet_cg_component_basis':[[[ser(c0[a][i,j]) for j in range(2)] for i in range(2)] for a in range(3)],
       'reference_weight_to_matchete_unitary_U':[[ser(U[i,j]) for j in range(3)] for i in range(3)],
       'bar_doublet_reference_generator_basis':'-conjugate(sigma_a/2)',
       'reference_weight_generators':[[[ser(g[i,j]) for j in range(3)] for i in range(3)] for g in tref],
       'matchete_triplet_candidate_generators':[[[ser(g[i,j]) for j in range(3)] for i in range(3)] for g in tmat],
       'su2_commutators_casimir_hermiticity_intertwiner_exact':True,
       'triplet_basis_complex_conjugation_diagnostic_U_UT':[[ser(K[i,j]) for j in range(3)] for i in range(3)],
       'native_rgbeta_s2_generator_components_exported':False,
       'native_matchete_generator_components_independently_exported':False,
       'global_matchete_rgbeta_cartesian_basis_certified':False,
       'physical_coupling_mapping_certified':False,
       'production_modified':False,
       'blocking_issue':'Need independent native Matchete/RGBeta triplet generator/component coordinates (including barred-field conjugation); this is an exact reference-basis candidate, not a physical map.',
       'six_triplet_checks':witnesses,
    }


def run(folder=DATA, output=None, check_stage26=True):
    mpath=folder/'matchete_tensors.json'; rpath=folder/'rgbeta_native_invariants.json'
    for path in (mpath,rpath):
        if not path.is_file(): fail('Missing ' + str(path))
    if check_stage26:
        previous=folder/'stage26_conditional_quartic_action.json'
        if not previous.is_file(): fail('Missing Stage26 output; rerun Stage26 first')
        saved=json.loads(previous.read_text(encoding='utf-8-sig'))
        fresh=stage26_run(folder,None)
        if saved!=fresh or fresh.get('physical_coupling_mapping_certified') is not False or \
           fresh.get('slots_checked')!=50 or fresh.get('blocks_checked')!=46:
            fail('Stage26 evidence missing, modified, stale, or unexpectedly certified')
    mat=json.loads(mpath.read_text(encoding='utf-8-sig'))
    nat=json.loads(rpath.read_text(encoding='utf-8-sig'))
    result=audit(mat,nat)
    result['input_sha256']={f:hashlib.sha256(p.read_bytes()).hexdigest()
                           for f,p in [('matchete_tensors.json',mpath),('rgbeta_native_invariants.json',rpath)]}
    if check_stage26:
        result['input_sha256']['stage26_conditional_quartic_action.json'] = hashlib.sha256(previous.read_bytes()).hexdigest()
    result['stage26_fresh_evidence_checked']=bool(check_stage26)
    if output is not None:
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,default=DATA)
    parser.add_argument('--output',type=Path,default=DATA/'stage27_triplet_generator_witness.json')
    ns=parser.parse_args()
    try:
        report=run(ns.input_dir,ns.output)
    except (ValueError, OSError, KeyError, TypeError, sp.SympifyError) as exc:
        parser.exit(2,str(exc)+'\n')
    print('Stage27 one Matchete triplet CG / SU2 generator witness: 6/6 contractions')
    print('SU2 commutator, Casimir, unitary intertwiner residuals: ZERO (exact)')
    print('Native RGBeta S2 component basis and physical coupling mapping: UNCERTIFIED')
    print('Report:',ns.output)


if __name__=='__main__':
    main()
