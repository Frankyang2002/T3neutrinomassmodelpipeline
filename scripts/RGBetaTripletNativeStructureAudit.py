"""Stage 28: live RGBeta symbolic S2 structure and Matchete triplet metric.

Uses six live RGBeta representation probes and freshly recomputed Stage 27
Matchete CG witness. An abstract tGen[S2] is not an independent numeric
component generator; the physical cross-package basis therefore stays blocked.

    python -m scripts.RGBetaTripletNativeStructureAudit
"""
from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path

import sympy as sp

from scripts.TripletIntertwiner import ALL, EPS, run as stage27_run, zero

DATA = Path('output/general_t3_stage22')
EXPECTED = {
    'A_+0': ((1, 3), ('S2',)), 'A_-2': ((1, 3), ('S2',)),
    'D_+0': ((3, 1), ('S1',)), 'D_-2': ((3, 1), ('S1',)),
    'E_+0': ((3, 3), ('S1', 'S2')), 'E_-2': ((3, 3), ('S1', 'S2')),
}


def fail(msg):
    raise ValueError('Stage28 fail-closed: ' + msg)


def scalar(expr):
    """Input is source-controlled Stage 27 expression strings, not user Python."""
    if not isinstance(expr, str) or len(expr) > 120:
        fail('Malformed Stage27 exact scalar')
    try:
        return sp.sympify(expr, locals={'I': sp.I, 'sqrt': sp.sqrt}, evaluate=True)
    except (sp.SympifyError, ValueError, SyntaxError) as exc:
        fail('Cannot read Stage27 scalar: ' + str(exc))


def mat3(rows):
    if not isinstance(rows, list) or len(rows) != 3 or any(
        not isinstance(row, list) or len(row) != 3 for row in rows
    ):
        fail('Expected three-by-three exact matrix')
    return sp.Matrix([[scalar(x) for x in row] for row in rows])


def report_for(stage27: dict, live: dict) -> dict:
    if not isinstance(stage27, dict) or not isinstance(live, dict):
        fail('Missing Stage27 or live Stage28 report')
    if stage27.get('source_native_matchete_tensor_count') != 6 or \
       stage27.get('global_matchete_rgbeta_cartesian_basis_certified') is not False or \
       stage27.get('physical_coupling_mapping_certified') is not False:
        fail('Stage27 record incomplete or unexpectedly certified')
    if stage27.get('models_triplet_checked') != list(ALL):
        fail('Stage27 six triplet labels changed')
    if live.get('format') != 'stage28-live-rgbeta-s2-structure-v1' or \
       live.get('models_checked') != 6 or \
       live.get('all_group_representation_checks_passed') is not True or \
       live.get('all_native_barred_generator_rules_passed') is not True or \
       live.get('all_native_delS2_contractions_passed') is not True or \
       live.get('native_numerical_triplet_generators_available') is not False or \
       live.get('physical_cross_package_component_basis_certified') is not False or \
       live.get('physical_coupling_mapping_certified') is not False:
        fail('Live RGBeta source structure absent or falsely certified')
    rows = live.get('models')
    if not isinstance(rows, list) or len(rows) != 6:
        fail('Native six-case structure list missing')
    got = {}
    for row in rows:
        if not isinstance(row, dict):
            fail('Malformed native model entry')
        label = row.get('model')
        if label not in EXPECTED or label in got:
            fail('Unexpected or repeated native triplet model')
        dims, fields = EXPECTED[label]
        if row.get('scalar_su2_dimensions') != list(dims) or \
           row.get('triplet_fields') != list(fields) or \
           row.get('native_S2_representation') != 'SU2L[S2]':
            fail('Native triplet field assignment mismatch: ' + label)
        if any(row.get(key) != val for key, val in {
            'dimension': 3,
            'quadratic_casimir': 2,
            'trace_normalisation': 2,
            'native_generator_is_symbolic': True,
            'native_conjugate_generator_minus_transpose': True,
            'native_delS2_orthogonality': True,
            'native_complex_conjugate_metric_two': True,
            'native_numerical_generator_components_exported': False,
        }.items()):
            fail('Native S2 algebra/metric probe failed: ' + label)
        if 'tGen[' not in row.get('native_generator_inputform', '') or \
           'S2' not in row['native_generator_inputform'] or \
           't28A' not in row['native_generator_inputform']:
            fail('Native generator is not an abstract tGen: ' + label)
        got[label] = row
    if set(got) != set(ALL):
        fail('Incomplete native triplet model set')

    source_cg = stage27.get('triplet_cg_component_basis')
    if not isinstance(source_cg, list) or len(source_cg) != 3 or any(
        not isinstance(q, list) or len(q) != 2 or any(
            not isinstance(r, list) or len(r) != 2 for r in q
        ) for q in source_cg
    ):
        fail('Stage27 Matchete CG array malformed')
    cg = [sp.Matrix([[scalar(z) for z in row] for row in q]) for q in source_cg]
    metric = sp.Matrix(3, 3, lambda a, b: sp.simplify(sum(
        cg[a][i, j] * cg[b][k, l] * EPS[i, k] * EPS[j, l]
        for i, j, k, l in product(range(2), repeat=4))))
    if not zero(metric + sp.eye(3)):
        fail('Stage27 Matchete triplet epsilon metric no longer equals -I_3')
    if metric.det() == 0 or not zero(metric.T - metric):
        fail('Derived triplet invariant metric degenerate or nonsymmetric')

    u = mat3(stage27['reference_weight_to_matchete_unitary_U'])
    gm = [mat3(g) for g in stage27['matchete_triplet_candidate_generators']]
    if not zero(u * u.H - sp.eye(3)) or not zero(u * u.T - sp.diag(1, 1, -1)):
        fail('Stage27 witness not consistent with the expected conjugation diagnostic')
    for gen in gm:
        # Invariant quadratic form satisfies T^T g + g T = 0. This identifies
        # an abstract self-dual real SU(2) representation in Matchete's basis,
        # NOT a native RGBeta numerical matrix or physical Bar transformation.
        if not zero(gen.T * metric + metric * gen):
            fail('Generator does not preserve Matchete triplet epsilon metric')
        if not zero(gen.H - gen):
            fail('Triplet Hermiticity lost')
    if any(not zero(gm[i]*gm[j]-gm[j]*gm[i]-sp.I*gm[k])
           for i,j,k in ((0,1,2),(1,2,0),(2,0,1))):
        fail('SU2 commutator check lost')
    if not zero(sum((g*g for g in gm),sp.zeros(3))-2*sp.eye(3)):
        fail('Triplet Casimir check lost')

    return {
        'format': 'stage28-triplet-native-structural-audit-v1',
        'status': 'RGBetaNativeAbstractS2AndMatcheteInvariantMetricVerifiedPhysicalBasisUncertified',
        'models_checked': 6,
        'triplet_slots_from_stage27': 6,
        'native_triplet_fields_in_six_models': sum(len(EXPECTED[x][1]) for x in ALL),
        'rgbeta_s2_dimension': 3,
        'rgbeta_s2_casimir': 2,
        'rgbeta_s2_trace_normalisation': 2,
        'rgbeta_barred_generator_minus_transpose_verified': True,
        'rgbeta_delS2_orthogonality_verified': True,
        'rgbeta_complex_field_contraction_metric_two_verified': True,
        'matchete_triplet_epsilon_invariant_metric':
            [[str(metric[i,j]) for j in range(3)] for i in range(3)],
        'matchete_triplet_invariant_metric_and_generators_exact': True,
        'native_rgbeta_numerical_triplet_generators_exported': False,
        'native_matchete_numerical_triplet_generators_independently_exported': False,
        'rgbeta_abstract_generator_does_not_select_numerical_coordinates': True,
        'barred_triplet_component_intertwiner_physically_certified': False,
        'global_matchete_rgbeta_cartesian_basis_certified': False,
        'physical_coupling_mapping_certified': False,
        'production_modified': False,
        'remaining_blocker': 'RGBeta tGen[S2] is an abstract representation. Need an independently fixed numerical basis AND the physical triplet Bar/conjugation map on both sides before certifying inter-package couplings.',
    }


def run(folder: Path = DATA, output: Path | None = None, *, check_stage27: bool = True):
    source = folder / 'stage28_live_rgbeta_triplet_structure.json'
    if not source.is_file():
        fail('Missing native Stage28 JSON; run its Wolfram exporter first')
    witness = folder / 'stage27_triplet_generator_witness.json'
    if not witness.is_file():
        fail('Missing Stage27 witness; rerun Stage27')
    if check_stage27:
        fresh = stage27_run(folder, None, check_stage26=True)
        saved = json.loads(witness.read_text(encoding='utf-8-sig'))
        if fresh != saved:
            fail('Stale Stage27 report: rerun Stage27')
    else:
        fresh = json.loads(witness.read_text(encoding='utf-8-sig'))
    native = json.loads(source.read_text(encoding='utf-8-sig'))
    result = report_for(fresh, native)
    result['stage27_fresh_evidence_checked'] = check_stage27
    result['input_sha256'] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (source, witness)
    }
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=DATA)
    parser.add_argument('--output', type=Path,
                        default=DATA/'stage28_triplet_native_structure_audit.json')
    ns = parser.parse_args()
    try:
        result = run(ns.input_dir, ns.output)
    except (ValueError, KeyError, TypeError, OSError, json.JSONDecodeError, sp.SympifyError) as exc:
        parser.exit(2, str(exc)+'\n')
    print('Stage28 RGBeta native abstract S2 structure: 6/6 models')
    print('Dimension 3, C2=2, T(S2)=2, Bar generator, delS2 and scalar metric: VERIFIED')
    print('Matchete triplet epsilon metric: -I_3; generator invariance: EXACT')
    print('Native numerical RGBeta generator coordinates: NOT DEFINED by the abstract package')
    print('Physical cross-package component basis and coupling map: UNCERTIFIED')
    print('Report:', ns.output)


if __name__ == '__main__':
    main()
