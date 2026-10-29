"""Stage 26: exact CONDITIONAL scalar-action polynomial / H.c. identity audit.

Uses fresh Stage 22--25 evidence; does not presume a global SU(2) or Cartesian
real-component intertwiner. RGBeta's published MBL Lagrangian Eq. (2.3)
uses -lambda phi^4. Eq. (2.9) phi=Delta.Phi/sqrt(2), while Eq. (2.13) uses
-Delta^4.Phi^4/4. Accordingly the native 6/4! = 1/4 cancels against
(sqrt(2))**4=4, yielding -lambda phi^4 for each ordered orientation.

We explicitly construct each 4-field monomial (including repeated-species
multiplicities) from *native* Matchete tensors and Stage22's reconstructed
RGBeta tensor basis; add the Hermitian conjugate as a separate polynomial;
then verify L_M + V_RGBeta == 0 under the Stage23 candidate signed map.

This is a source-convention and conditional action-algebra check, NEVER a
certification of the global Matchete/RGBeta component basis.

python -m scripts.QuarticActionIdentity
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import product
import json
from pathlib import Path

import sympy as sp

from scripts.RGBetaTensorAlgebra import exact
from scripts.RGBetaAlignedGeometry import triplet_tensor
from scripts.RGBetaExplicitBasis import reconstructed_tensor
from scripts.QuarticPrefactorAudit import exact_matrix, run as stage23_run
from scripts.NativeQuarticActionAudit import run as stage25_run

ROOT = Path(__file__).resolve().parents[1]
DATA = Path('output/general_t3_stage22')
# See RGBeta (Thomsen, EPJC 81 (2021) 408), Eqs. 2.3, 2.9, 2.13, 2.15.
RGBETA_REFERENCE = 'https://doi.org/10.1140/epjc/s10052-021-09142-4'


def fail(message: str):
    raise ValueError('Stage26 fail-closed: ' + message)


def scalar_action_conversion():
    """Exactly check the two independent 4-field normalisations (no floats)."""
    from_real_tensor = sp.factorial(4) * sp.Pow(sp.sqrt(2), -4) / sp.factorial(4)
    real_to_complex = sp.sqrt(2) ** 4
    if sp.simplify(from_real_tensor * real_to_complex - 1) != 0 or \
       sp.simplify(from_real_tensor - sp.Rational(1, 4)) != 0:
        fail('RGBeta real-to-complex quartic conversion is not unity')
    return {'rgbeta_native_lambda_tensor_factor': '6',
            'real_tensor_action_factor': '1/24',
            'complex_projection_product_factor': '4',
            'net_monomial_action_factor': '1',
            'rgbeta_quartic_lagrangian_sign': '-1',
            'matchete_extra_quartic_lagrangian_sign': '+1',
            'relative_signed_candidate_map': '-A^-1'}


def verify_matchete_source(root=ROOT):
    """Guard the project source assumptions; NOT a live Matchete action proof."""
    builder = root/'Lagrangian/model/LagrangianBuilder.wl'
    exporter = root/'scripts/ExportGeneralT3ScalarTensorMatchete.wl'
    if not builder.is_file() or not exporter.is_file():
        fail('Missing current LagrangianBuilder.wl or Stage22 Matchete exporter')
    b = builder.read_text(encoding='utf-8-sig')
    e = exporter.read_text(encoding='utf-8-sig')
    if 'LInt = Total[Lookup[valid, "Expression", {}]]' not in b or \
       'LUV = (LSM + LBSM)' not in b or \
       '"potential_prefactor_convention" -> "+lambda[] times contracted fields in PlusHc[contribution]"' not in e or \
       '"potential_prefactor" -> "1"' not in e:
        fail('Matchete Lagrangian construction or exported quartic coefficient convention changed')
    return True


def tensor_polynomial(tensor):
    """Exact commuting 4-field monomial dictionary, with NO extra 4! factor.

    Identical fields naturally collect: (phi_a phi_b)=(phi_b phi_a).
    This models the exporter-specified contracted monomial, not a new
    physical component-coordinate transformation.
    """
    species = tensor['species']
    dimensions = tensor['dimensions']
    values = tensor['components']
    if len(species) != 4 or len(dimensions) != 4 or any(type(x) is not int or x not in (1, 2, 3) for x in dimensions):
        fail('Non-quartic tensor or unsupported dimension')
    if len(values) != sp.prod(dimensions):
        fail('Component vector has wrong dimensions')
    out = defaultdict(lambda: sp.Integer(0))
    for indices, component in zip(product(*(range(d) for d in dimensions)), values):
        monomial = tuple(sorted(zip(species, indices)))
        out[monomial] += exact(component)
    return {key: sp.simplify(value) for key, value in out.items() if sp.simplify(value) != 0}


def hermitian_conjugate(poly):
    """Formal complex conjugation in the provisional *same* component labels."""
    conjugates = {'H':'Hdag', 'Hdag':'H', 'S1':'S1dag', 'S1dag':'S1',
                  'S2':'S2dag', 'S2dag':'S2'}
    result = defaultdict(lambda: sp.Integer(0))
    for monomial, coeff in poly.items():
        if any(field not in conjugates for field, _ in monomial):
            fail('Unknown field in polynomial conjugation')
        key = tuple(sorted((conjugates[field], component) for field, component in monomial))
        result[key] += sp.conjugate(coeff)
    return {key: sp.simplify(value) for key, value in result.items() if sp.simplify(value) != 0}


def assemble_action(tensors, coefficients, plus_hc=True):
    if len(tensors) != len(coefficients):
        fail('Missing scalar action coefficients')
    result = defaultdict(lambda: sp.Integer(0))
    for tensor, coefficient in zip(tensors, coefficients):
        monomials = tensor_polynomial(tensor)
        if not monomials:
            fail('Vanishing quartic operator')
        for monomial, factor in monomials.items():
            result[monomial] += coefficient * factor
        if plus_hc:
            for monomial, factor in hermitian_conjugate(monomials).items():
                result[monomial] += sp.conjugate(coefficient) * factor
    return {key: sp.simplify(val) for key, val in result.items() if sp.simplify(val) != 0}


def actions_equal(first, second):
    return all(sp.simplify(first.get(key, 0)-second.get(key, 0)) == 0
               for key in set(first)|set(second))


def reconstruct(tensor_native, matchete_ref, cls):
    if cls == 'CalibratedTripletCGCandidate':
        tensor, _ = triplet_tensor(tensor_native, matchete_ref)
        return tensor
    if cls in ('DirectGeometryCandidate', 'CorrectedBCReusedInvariantCandidate',
               'CorrectedReusedInvariantCandidate'):
        return reconstructed_tensor(tensor_native, matchete_ref)
    fail('Uncertified geometry class: ' + str(cls))


def audit_block(block, mat_by_symbol, rgb_by_symbol):
    mat_names = block['matchete_couplings']
    rgb_names = block['rgbeta_couplings']
    cls = block['geometry_class']
    if not isinstance(mat_names, list) or not isinstance(rgb_names, list) or not mat_names or \
       len(mat_names) != len(rgb_names):
        fail('Invalid block coupling inventory')
    if block.get('physical_mapping_certified') is not False and \
       block.get('physical_mapping_certified') is not None:
        fail('Unexpected physical certification in conditional block')
    if block.get('assumed_relative_potential_sign') != '-1' or \
       block.get('component_prefactor_one_ordered_orientation') != '1/4':
        fail('Stage23 conditional sign or component factor changed')
    try:
        mats = [mat_by_symbol[s] for s in mat_names]
        rgbs = [reconstruct(rgb_by_symbol[s], mats[0], cls) for s in rgb_names]
    except (KeyError, IndexError) as e:
        fail('Missing block tensor: ' + str(e))
    for m, r in zip(mats, rgbs):
        if m['species'] != r['species'] or m['dimensions'] != r['dimensions']:
            fail('Mismatched species or dimensions in action basis')
    geometry = exact_matrix(block['r_in_m'])
    candidate = exact_matrix(block['lambda_rgbeta_from_lambda_matchete_if_assumptions_hold'])
    if geometry.shape != (len(mats), len(mats)) or \
       any(sp.simplify(v) != 0 for v in (geometry*candidate + sp.eye(len(mats)))):
        fail('Signed block map does not invert the geometric tensor map')
    couplings = sp.Matrix([sp.Symbol(f't26m{i}', complex=True) for i in range(len(mats))])
    rgbeta_couplings = candidate * couplings
    l_mat = assemble_action(mats, list(couplings), plus_hc=True)
    v_rgbeta = assemble_action(rgbs, list(rgbeta_couplings), plus_hc=True)
    if not actions_equal(l_mat, {k:-v for k,v in v_rgbeta.items()}):
        fail('Full quartic action identity FAILED: ' + block['model'] + '/' + block['block'])
    # Negative control: sign or 1/4 cannot silently disappear.
    if actions_equal(l_mat, v_rgbeta):
        fail('Quartic sign control failed: Lagrangian agrees with potential without minus sign')
    if actions_equal(l_mat, {k:-sp.Rational(1,4)*v for k,v in v_rgbeta.items()}):
        fail('Quartic normalisation control failed: spurious factor 1/4 was accepted')
    return {
        'model': block['model'], 'block': block['block'], 'geometry_class': cls,
        'quartic_couplings': len(mat_names),
        'ordered_orientation_monomials': sum(len(tensor_polynomial(x)) for x in mats),
        'full_action_monomials_matchete': len(l_mat),
        'full_action_monomials_rgbeta': len(v_rgbeta),
        'exact_signed_action_identity_in_provisional_basis': True,
        'hermitian_conjugate_included': True,
        'global_component_basis_verified': False,
        'physical_mapping_certified': False,
    }


def audit(data, stage23, stage25):
    if stage25.get('slots_checked') != 50 or stage25.get('blocks_covered') != 46 or \
       stage25.get('physical_coupling_mapping_certified') is not False or \
       stage25.get('live_rgbeta_coupling_tensor_map_checked') is not True or \
       stage25.get('live_complex_field_metric_two_checked') is not True:
        fail('Stage25 live native assembly is missing or uncertified')
    if stage23.get('totals', {}).get('slots') != 50 or stage23.get('totals', {}).get('blocks') != 46 or \
       len(stage23.get('blocks', [])) != 46 or \
       stage23.get('signed_physical_coupling_map_certified') is not False:
        fail('Stage23 candidate maps are stale or incomplete')
    mats, rgb = data['matchete_tensors.json'], data['rgbeta_native_invariants.json']
    if mats.get('quartic_slot_count') != 50 or rgb.get('total_slots') != 50 or \
       set(mats.get('models',{})) != set(rgb.get('models',{})) or len(mats['models']) != 16:
        fail('Incomplete 16-model native tensor evidence')
    indexed_m = {m: {row['symbol']: row for row in mats['models'][m]['tensors']} for m in mats['models']}
    indexed_r = {m: {row['symbol']: row for row in rgb['models'][m]['invariants']} for m in rgb['models']}
    if any(len(indexed_m[m]) != len(mats['models'][m]['tensors']) or
           len(indexed_r[m]) != len(rgb['models'][m]['invariants']) for m in mats['models']):
        fail('Duplicate native tensor symbol')
    rows, visited_m, visited_r = [], set(), set()
    for block in stage23['blocks']:
        model = block['model']
        if model not in indexed_m or model not in indexed_r:
            fail('Unknown block model: ' + model)
        row = audit_block(block, indexed_m[model], indexed_r[model])
        rows.append(row)
        visited_m.update((model, sym) for sym in block['matchete_couplings'])
        visited_r.update((model, sym) for sym in block['rgbeta_couplings'])
    expected_m = {(model, symbol) for model, table in indexed_m.items() for symbol in table}
    expected_r = {(model, symbol) for model, table in indexed_r.items() for symbol in table}
    if visited_m != expected_m or visited_r != expected_r or \
       sum(r['quartic_couplings'] for r in rows) != 50:
        fail('Conditional quartic action check did not exhaust exactly 50 slots')
    return {
        'format': 'stage26-conditional-quartic-action-polynomials-v1',
        'status': 'SourceSignAndConditionalQuarticActionCheckedGlobalBasisStillUncertified',
        'source_convention': scalar_action_conversion(),
        'rgbeta_published_source': RGBETA_REFERENCE,
        'published_source_equations': ['2.3','2.9','2.13','2.15'],
        'matchete_project_builder_plus_interaction_checked': True,
        'models_checked': len(mats['models']),
        'slots_checked': 50, 'blocks_checked': len(rows),
        'native_stage25_permutation_rules_checked': stage25['native_permutation_rules_checked'],
        'signed_polynomial_identity_conditional_verified': True,
        'hermitian_conjugate_formal_polynomial_checked': True,
        'source_relative_action_sign_supported': True,
        'global_su2_component_basis_verified': False,
        'actual_cartesian_real_field_action_verified': False,
        'physical_coupling_mapping_certified': False,
        'production_modified': False,
        'blocking_issue': 'Independently identify native RGBeta/Matchete triplet component basis and verify the 6 calibrated CG cases; 46 conditional maps must not be promoted to physical couplings.',
        'blocks': rows,
    }


def run(folder=DATA, output=None, root=ROOT):
    verify_matchete_source(root)
    # Both regenerate and check against stored reports: stale evidence fails.
    s25 = stage25_run(folder, None)
    stored25_path = folder/'stage25_native_quartic_assembly_diagnostics.json'
    if not stored25_path.is_file() or json.loads(stored25_path.read_text(encoding='utf-8-sig')) != s25:
        fail('Stage25 stored diagnostic missing or stale')
    s23 = stage23_run(folder, None)
    input_files = ('matchete_tensors.json','rgbeta_native_invariants.json')
    data = {}
    for f in input_files:
        fp = folder/f
        if not fp.is_file():
            fail('Missing ' + str(fp))
        data[f] = json.loads(fp.read_text(encoding='utf-8-sig'))
    report = audit(data, s23, s25)
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    return report


def main():
    pa = argparse.ArgumentParser(description=__doc__)
    pa.add_argument('--input-dir', type=Path, default=DATA)
    pa.add_argument('--project-root', type=Path, default=ROOT)
    pa.add_argument('--output', type=Path, default=DATA/'stage26_conditional_quartic_action.json')
    args = pa.parse_args()
    try:
        d = run(args.input_dir, args.output, args.project_root)
    except (ValueError, OSError, KeyError, TypeError, sp.SympifyError) as e:
        pa.exit(2, str(e) + '\n')
    print(f"Stage26 signed quartic action + H.c.: {d['slots_checked']}/50 slots, {d['blocks_checked']}/46 blocks")
    print('RGBeta published relative sign and 4-field normalisation: SOURCE-SUPPORTED')
    print('Global Matchete/RGBeta Cartesian field basis: UNCERTIFIED')
    print('Physical coupling map: UNCERTIFIED')
    print('Report:', args.output)


if __name__ == '__main__':
    main()
