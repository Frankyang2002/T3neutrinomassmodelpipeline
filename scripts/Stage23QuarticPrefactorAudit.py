"""Stage 23: read-only, exact quartic prefactor audit and CONDITIONAL coupling maps.

The native RGBeta quartic tensor is normalised as 24*(1/sqrt(2))**4=6
for four complex scalar slots (RGBeta v1.2 FieldsAndCouplings.m AddQuartic).
CouplingTensorMap averages over native UniqueArrangements.  This check DOES NOT
establish a global Matchete/RGBeta field-component identification or certify
potential signs. No running/production files are modified.

python -m scripts.Stage23QuarticPrefactorAudit
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

import sympy as sp

FOLDER = Path('output/general_t3_stage22')
INPUTS = ('matchete_tensors.json', 'rgbeta_native_declarations.json',
          'rgbeta_native_invariants.json', 'stage22_aligned_geometry_diagnostics.json')
EXPECTED_CLASSES = {
    'DirectGeometryCandidate': 32,
    'PseudorealReindexCandidate': 8,
    'CalibratedTripletCGCandidate': 6,
    'CorrectedReusedInvariantCandidate': 4,
}
POST_STAGE24_CLASSES = {
    'DirectGeometryCandidate': 32,
    'CorrectedBCReusedInvariantCandidate': 8,
    'CalibratedTripletCGCandidate': 6,
    'CorrectedReusedInvariantCandidate': 4,
}
VALID_CLASS_PROFILES = (EXPECTED_CLASSES, POST_STAGE24_CLASSES)
CORRECTED = 'del[SU2L[fund],t3i1,t3i3]*eps[SU2L[fund],t3i2,t3i4]'
CORRECTED_FUNCTION = 'del[SU2L[fund],#1,#3]*eps[SU2L[fund],#2,#4]&'
PREFIX = '6*Matrix[{}][##1]&'


def fail(msg):
    raise ValueError('Stage23 fail-closed: ' + msg)


def compact(s):
    return re.sub(r'\s+', '', s)


def parse_permutations(s):
    if not isinstance(s, str) or not re.fullmatch(r'[\s{},1-4]+', s):
        fail('Malformed native UniqueArrangements')
    arr = re.findall(r'\{\s*([1-4]\s*,\s*[1-4]\s*,\s*[1-4]\s*,\s*[1-4])\s*\}', s)
    reconstructed = '{' + ','.join('{' + compact(x) + '}' for x in arr) + '}'
    if compact(reconstructed) != compact(s):
        fail('Unexpected UniqueArrangements serialization')
    vals = [tuple(int(y) for y in x.split(',')) for x in arr]
    if not vals or len(vals) != len(set(vals)) or any(sorted(v) != [1, 2, 3, 4] for v in vals):
        fail('Non-unique or non-permutation native arrangements')
    return vals


def scalar_component_prefactor(native_factor, nperms):
    """Coefficient per term after 1/4! real-tensor contraction and averaging.

    This proves 6/24=1/4 as a real-field *component* coefficient for one
    unsymmetrised contraction. It is not a full physical tensor comparison.
    """
    if native_factor != 6 or not isinstance(nperms, int) or not 1 <= nperms <= 24:
        fail('Nonstandard RGBeta quartic coupling normalisation')
    return sp.Rational(native_factor, nperms) * nperms / sp.factorial(4)


def exact_matrix(raw):
    if not isinstance(raw, list) or not raw or any(not isinstance(r, list) or
                                                       len(r) != len(raw) for r in raw):
        fail('Geometry matrix must be square')
    if len(raw) > 2:
        fail('Unexpected dimension beyond the certified Stage22 geometric inventory')
    vals = []
    for row in raw:
        v = []
        for elem in row:
            if not isinstance(elem, str) or not re.fullmatch(r'[0-9 +\-*/()sqrtI]+', elem):
                fail('Invalid exact algebraic matrix entry')
            expr = sp.sympify(elem, locals={'sqrt': sp.sqrt, 'I': sp.I})
            if expr.free_symbols or sp.simplify(sp.im(expr)) != 0:
                fail('Unsupported or complex geometry matrix entry')
            v.append(expr)
        vals.append(v)
    mat = sp.Matrix(vals)
    if sp.simplify(mat.det()) == 0:
        fail('Singular block: cannot recover coupling coefficients')
    return mat


def plain(expr):
    return str(sp.simplify(expr))


def required_inputs(base):
    data = {}
    for name in INPUTS:
        path = base / name
        if not path.is_file():
            fail(f'Missing {path}; regenerate native Stage22 exports and aligned geometry')
        data[name] = json.loads(path.read_text(encoding='utf-8-sig'))
    return data


def calculate(data):
    mat = data[INPUTS[0]]
    dec = data[INPUTS[1]]
    inv = data[INPUTS[2]]
    geo = data[INPUTS[3]]
    if mat.get('quartic_slot_count') != 50 or inv.get('total_slots') != 50:
        fail('Incomplete Matchete/RGBeta quartic inventories')
    models = mat.get('models', {})
    if set(models) != set(dec.get('models', {})) or set(models) != set(inv.get('models', {})) or set(models) != set(geo.get('models', {})) or len(models) != 16:
        fail('Missing or inconsistent model keys')
    profile = next((p for p in VALID_CLASS_PROFILES
                    if geo.get('totals') == {'slots': 50, 'blocks': 46, **p}), None)
    if profile is None:
        fail('Aligned geometry is not a reviewed pre/post Stage24 whole-model class profile')
    if geo.get('physical_coupling_map_certified') is not False or geo.get('physical_prefactors_verified') is not False:
        fail('Geometry certification fields inconsistent with conditional pipeline')
    status_counts = Counter()
    rows = []
    coupled_symbols = []
    for model in sorted(models):
        mt = models[model]['tensors']
        dd = dec['models'][model]['declarations']
        ii = inv['models'][model]['invariants']
        gg = geo['models'][model]
        mby = {x['symbol']: x for x in mt}
        dby = {x['symbol']: x for x in dd}
        iby = {x['symbol']: x for x in ii}
        if len(mby) != len(mt) or len(dby) != len(dd) or len(iby) != len(ii):
            fail('Duplicate slot symbol in ' + model)
        for slot in dd:
            symbol = slot['symbol']
            if symbol not in iby or slot.get('self_conjugate') is not False or iby[symbol].get('self_conjugate') is not False:
                fail('Inconsistent native quartic metadata: ' + model + '/' + symbol)
            if compact(slot['coupling_factor_inputform']) != PREFIX.format(symbol) or compact(iby[symbol]['native_coupling_function_inputform']) != PREFIX.format(symbol):
                fail('Native normalisation is not exactly 6*Matrix: ' + model + '/' + symbol)
            if compact(slot['unique_arrangements_inputform']) != compact(iby[symbol]['native_unique_arrangements_inputform']):
                fail('Native UniqueArrangements disagree: ' + model + '/' + symbol)
            if compact(slot['invariant_inputform']) != compact(iby[symbol]['native_invariant_function_inputform']):
                fail('Native invariant functions disagree: ' + model + '/' + symbol)
            if model in ('B_-1','C_-1') and symbol == 'lambdaHHdagS1barS2barCross':
                if (compact(iby[symbol]['native_invariant_expression_inputform']) != CORRECTED or
                        compact(slot['invariant_inputform']) != CORRECTED_FUNCTION):
                    fail('B/C reused invariant still pre-fix or stale: ' + model)
            perms = parse_permutations(slot['unique_arrangements_inputform'])
            if scalar_component_prefactor(6, len(perms)) != sp.Rational(1, 4):
                fail('RGBeta real-component quartic prefactor inconsistency')
            coupled_symbols.append((model, symbol, len(perms)))
        if set(dby) != set(iby):
            fail('Native declaration and invariant symbol sets differ')
        used_m, used_r = set(), set()
        for block in gg:
            cls = block['status']
            if cls not in profile:
                fail('Unrecognised provisional geometry class')
            status_counts[cls] += len(block['rgbeta_symbols'])
            msyms = block['matchete_symbols']
            rsyms = block['rgbeta_symbols']
            if len(msyms) != len(rsyms) or len(set(msyms)) != len(msyms) or len(set(rsyms)) != len(rsyms):
                fail('Invalid block symbol count')
            if any(s not in mby or mby[s]['potential_prefactor'] != '1' for s in msyms):
                fail('Missing Matchete scalar potential coefficient 1')
            if any(s not in dby for s in rsyms):
                fail('Missing native RGBeta block symbols')
            if used_m.intersection(msyms) or used_r.intersection(rsyms):
                fail('Coupling symbol used in multiple geometry blocks')
            used_m.update(msyms); used_r.update(rsyms)
            geom = block['geometric']
            if (geom['matchete_symbols'] != msyms or geom['rgbeta_symbols'] != rsyms or
                geom['counts'] != block['counts'] or geom.get('maximum_exact_component_residual') != '0' or
                not geom.get('tensor_equality') or not geom.get('zero_coupling_recovery')):
                fail('Inexact or incoherent basis geometry')
            A = exact_matrix(geom['r_in_m'])
            inverse = A.inv()
            candidate = -inverse  # Matchete Lagrangian = - RGBeta quartic potential
            if any(sp.simplify(x) != 0 for x in A*candidate + sp.eye(A.rows)):
                fail('Failed signed map recovery')
            rows.append({
                'model': model, 'block': block['block'], 'geometry_class': cls,
                'matchete_couplings': msyms, 'rgbeta_couplings': rsyms,
                'r_in_m': [[plain(A[i,j]) for j in range(A.cols)] for i in range(A.rows)],
                'lambda_rgbeta_from_lambda_matchete_if_assumptions_hold':
                    [[plain(candidate[i,j]) for j in range(candidate.cols)] for i in range(candidate.rows)],
                'native_rgbeta_permutation_counts': {s:len(parse_permutations(dby[s]['unique_arrangements_inputform'])) for s in rsyms},
                'component_prefactor_one_ordered_orientation': '1/4',
                'assumed_relative_potential_sign': '-1',
                'physical_mapping_certified': False,
            })
        if used_m != set(mby) or used_r != set(dby):
            fail('Geometric blocks do not exhaust native coupling slots in ' + model)
    if status_counts != Counter(profile) or len(coupled_symbols) != 50 or len(rows) != 46:
        fail('Unexpected 16-model slot/block total')
    return {
        'format': 'stage23-exact-source-prefactors-conditional-signed-maps-v1',
        'status': 'SourceQuarticNormalisationVerifiedConditionalSignedMapsNotPhysicalCertification',
        'source_verified': {
            'rgbeta_addquartic_complex_slot_count': 4,
            'rgbeta_quartic_real_tensor_factor': '24*(1/sqrt(2))**4=6',
            'rgbeta_permutation_averaging': 'CouplingTensorMap/AvgPermutations uses 1/N unique arrangements',
            'rgbeta_real_component_coefficient_per_orientation': '6/N*N/24=1/4',
            'matchete_export_potential_prefactors_all': '1',
            'matchete_export_records': 50,
            'native_rgbeta_coupling_functions_checked': 50,
            'native_rgbeta_unique_arrangements_checked': 50,
            'upstream_rgbeta_source': 'aethomsen/RGBeta Package/FieldsAndCouplings.m lines 751-787, 931-941, 1035-1046',
            'matchete_source': 'scripts/ExportGeneralT3Stage22Matchete.wl and Lagrangian/interactions/ScalarPotential.wl',
        },
        'conditional_assumptions': [
            'The physical Matchete interaction uses +lambda_M multiplied by the exported invariant in the Lagrangian',
            'RGBeta defines quartic interactions in the Lagrangian with the conventional -V sign and V=Lam real-scalar quartic/24',
            'A single global conjugate-doublet component transformation must be justified if the historical Stage22 B/C J reindexing profile is used',
            'The symmetric-doublet triplet CG has the same physical normalisation as RGBeta native delS2 for all six cases',
            'The corrected native B/C 110101 invariant is used in live RGBeta (source exports must be fresh)',
            'Stage24 remaining four reused invariants use their reviewed exact component tensors when the Stage24 corrected profile is active',
        ],
        'source_prefactors_verified': True,
        'global_su2_component_basis_verified': False,
        'live_real_scalar_quartic_action_verified': False,
        'signed_physical_coupling_map_certified': False,
        'production_modified': False,
        'totals': {'models': 16, 'slots': 50, 'blocks': 46, **dict(status_counts)},
        'blocks': rows,
    }


def run(base=FOLDER, output=None):
    data = required_inputs(base)
    for filename, value in zip(INPUTS[:3], [data[x] for x in INPUTS[:3]]):
        expected = data[INPUTS[3]].get('input_sha256', {}).get(filename)
        if not expected or hashlib.sha256((base / filename).read_bytes()).hexdigest() != expected:
            fail('Aligned geometry input hashes stale: ' + filename)
    report = calculate(data)
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input-dir', type=Path, default=FOLDER)
    ap.add_argument('--output', type=Path,
                    default=FOLDER / 'stage23_quartic_prefactor_candidates.json')
    args = ap.parse_args()
    try:
        result = run(args.input_dir, args.output)
    except (KeyError, ValueError, TypeError, OSError) as exc:
        ap.exit(2, str(exc) + '\n')
    print('Stage 23 exact native quartic prefactors: 50/50; 46/46 blocks')
    print('Conditional signed coupling-map matrices:', len(result['blocks']))
    print('Physical Matchete/RGBeta coupling mapping: UNCERTIFIED')
    print('Report:', args.output)


if __name__ == '__main__':
    main()
