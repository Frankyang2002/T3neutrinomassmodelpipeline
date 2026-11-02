"""Fail-closed export of the 16 corrected fermion-first LLSS boundaries.

Project destination: output/general_observable/lo_fermion_boundary/build_boundary.py
Only existing corrected local sources are read. No Matchete, RG, scalar or
observable worker is modified or invoked. Absolute GammaCC/Weyl spinor
normalization is deliberately NOT asserted.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import sys

MODELS = tuple(f'{f}_{a:+d}' for f, aa in (
    ('A', (-4, -2, 0)), ('B', (-3, -1, 1)),
    ('C', (-3, -1, 1)), ('D', (-2, 0, 2)),
    ('E', (-4, -2, 0, 2)),
) for a in aa)
MIXED = {'B_-1': 10, 'C_-1': 12}
ORDERED_FORMAT = 't3-ordered-source-flavour-and-real-scalar-tensor-v1'
POLY_FORMAT = 't3-general-16-polynomial-wilson-bridge-v1'
AUDIT_FORMAT = 't3-16-matchete-to-master-normalization-audit-v1'
AUDIT_OK = 'PolynomialAndPLStructureVerified_SpinorNormalizationOpen'
LEPTON = re.compile(
    r'Field\[l,\s*Fermion,\s*\{\s*'
    r'Index\[(?P<su2>[^,\]]+),\s*SU2L\[fund\]\],\s*'
    r'Index\[(?P<flav>[^,\]]+),\s*Flavor\]')
PL = re.compile(r'Proj\s*\[\s*-1\s*\]')
PR = re.compile(r'Proj\s*\[\s*1\s*\]')
FORMAT = 't3-lo-fermion-boundary-source-resolved-provisional-v1'


def ensure(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    ensure(isinstance(obj, dict), f'Expected JSON object: {path}')
    return obj


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tag(model: str) -> str:
    return model.replace('-', 'm').replace('+', 'p')


def canon(key: tuple[int, ...]) -> tuple[int, ...]:
    ensure(len(key) == 6, 'Wrong coordinate rank')
    i, p, j, q, a, b = key
    ensure(i in (1, 2) and j in (1, 2) and p in (1, 2, 3) and q in (1, 2, 3),
           'Invalid SU(2)/flavour indices')
    ensure(a >= 1 and b >= 1, 'Invalid scalar indices')
    if (i, p) > (j, q):
        i, p, j, q = j, q, i, p
    return (i, p, j, q, min(a, b), max(a, b))


def orbit(key: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    i, p, j, q, a, b = canon(key)
    return tuple(sorted({(i, p, j, q, a, b), (j, q, i, p, a, b),
                         (i, p, j, q, b, a), (j, q, i, p, b, a)}))


def complex_matrix(raw: dict) -> list[list[complex]]:
    ensure(isinstance(raw, dict) and set(('real', 'imag')) <= set(raw), 'Matrix encoding missing')
    real, imag = raw['real'], raw['imag']
    ensure(len(real) == len(imag) == 3, 'Wrong flavour matrix shape')
    result = []
    for r, s in zip(real, imag):
        ensure(len(r) == len(s) == 3, 'Wrong flavour matrix width')
        row = [complex(float(a), float(b)) for a, b in zip(r, s)]
        ensure(all(math.isfinite(z.real) and math.isfinite(z.imag) for z in row),
               'Nonfinite flavour coefficient')
        result.append(row)
    return result


def rows_from_aligned(raw: dict[tuple[int, ...], complex]) -> list[dict]:
    grouped: dict[tuple[int, ...], complex] = defaultdict(complex)
    for key, value in raw.items():
        grouped[canon(key)] += complex(value)
    rows = []
    for key, total in sorted(grouped.items()):
        if total == 0j:
            continue
        value = total / len(orbit(key))
        rows.append({'coordinate': list(key), 'coefficient_real': float(value.real),
                     'coefficient_imag': float(value.imag), 'orbit_size': len(orbit(key))})
    return rows


def as_map(rows: list[dict]) -> dict[tuple[int, ...], complex]:
    result = {}
    for row in rows:
        k = tuple(row['coordinate'])
        ensure(k == canon(k), 'Noncanonical coordinate')
        ensure(row['orbit_size'] == len(orbit(k)), 'Incorrect orbit size')
        ensure(k not in result, 'Duplicate coordinate')
        z = complex(float(row['coefficient_real']), float(row['coefficient_imag']))
        ensure(math.isfinite(z.real) and math.isfinite(z.imag), 'Nonfinite Wilson entry')
        result[k] = z
    return result


def compare_rows(left: list[dict], right: list[dict], label: str) -> None:
    a, b = as_map(left), as_map(right)
    scale = max((abs(x) for x in b.values()), default=0)
    ensure(scale > 0, label + ': empty or zero Wilson boundary')
    keys = set(a) | set(b)
    error = max((abs(a.get(k, 0j) - b.get(k, 0j)) for k in keys), default=0)
    ensure(error <= 4e-12 * scale, f'{label}: numerical mismatch {error/scale:g}')


def source_file(model: str, root: Path, ordered_root: Path, ordered: dict) -> Path:
    if model in MIXED:
        return ordered_root / tag(model) / 'source_flavour.json'
    name = ordered.get('source_boundary')
    ensure(isinstance(name, str) and bool(name), f'{model}: missing original source location')
    location = Path(name)
    if not location.is_absolute():
        location = root / location
    ensure(location.is_file(), f'{model}: original source file absent: {location}')
    ensure(digest(location) == ordered.get('source_boundary_sha256'),
           f'{model}: original Matchete boundary SHA mismatch')
    return location


def chirality(raw: str) -> str:
    ensure(isinstance(raw, str) and 'GammaCC' in raw, 'Missing GammaCC LLSS bilinear')
    is_pl, is_pr = bool(PL.search(raw)), bool(PR.search(raw))
    ensure(is_pl != is_pr, 'Missing or ambiguous Matchete chiral projector')
    if is_pl:
        ensure(bool(re.search(r'Transp\s*\[\s*Field\s*\[\s*l\s*,\s*Fermion', raw))
               and not bool(re.search(r'Transp\s*\[\s*Bar\s*\[\s*Field\s*\[\s*l', raw)),
               'PL lepton transpose is not the plain bilinear')
        return 'PL'
    ensure(bool(re.search(r'Transp\s*\[\s*Bar\s*\[\s*Field\s*\[\s*l\s*,\s*Fermion', raw)),
           'PR conjugate transpose does not contain Bar')
    return 'PR'


def checked_terms(model: str, source: dict, ordered: dict) -> dict[int, dict]:
    terms = source.get('operator_terms')
    ensure(source.get('model') == model and isinstance(terms, list) and bool(terms),
           f'{model}: invalid native source term inventory')
    selected = {}
    count = Counter()
    for term in terms:
        natural = chirality(term['source_term_inputform'])
        declared = term.get('chirality', 'PL' if model not in MIXED else None)
        ensure(declared == natural, f'{model}: incorrect {declared} on Matchete {natural} source')
        count[natural] += 1
        if natural == 'PL':
            ix = term['source_index']
            ensure(isinstance(ix, int) and ix > 0 and ix not in selected,
                   f'{model}: duplicate/invalid source index')
            selected[ix] = term
    ensure(count['PL'] == ordered['tree_operator_count_PL'], f'{model}: incomplete PL inventory')
    ensure(count['PR'] == (MIXED.get(model, 0)), f'{model}: incomplete/excess PR sector')
    ensure(len(selected) == len(ordered['operators']) and
           set(selected) == {op['source_index'] for op in ordered['operators']},
           f'{model}: tensor and source term indices differ')
    return selected


def validate_provenance(model: str, root: Path, ordered_root: Path,
                        ordered: dict, source: dict) -> None:
    ensure(ordered.get('status') == 'OrderedNumericalRealComponentTensorFromMatchete',
           f'{model}: tensor not verified')
    ensure(ordered.get('model') == model and source.get('model') == model,
           f'{model}: inconsistent source model')
    ensure(ordered.get('canonical_C_ijab_normalization_certified') is False,
           f'{model}: unexpected previous physical normalisation claim')
    ensure(float(ordered['mu_F_gev']) == float(source['mu_F_gev']),
           f'{model}: threshold scales disagree')
    if model in MIXED:
        ensure(ordered['tree_operator_count_PL'] == MIXED[model],
               f'{model}: mixed chiral inventory mismatch')
        ensure(ordered.get('format') == 't3-mixed-native-matchete-ordered-tensor-v1' and
               source.get('format') == 't3-mixed-factorised-LLSS-source-v1' and
               source.get('status') == 'NumericalOrderedSourceFlavourCGExported',
               f'{model}: obsolete B/C mixed tensor converter')
        for a, b in [('source_status_sha256', 'source_status_sha256'),
                     ('source_mass_basis_sha256', 'source_mass_basis_sha256'),
                     ('source_native_fresh_llss_sha256', 'llss_sha256')]:
            ensure(ordered.get(a) == source.get(b) and isinstance(ordered.get(a), str)
                   and len(ordered[a]) == 64, f'{model}: old or mismatched native source {a}')
    else:
        ensure(ordered.get('format') == ORDERED_FORMAT and
               source.get('format') == 't3-source-exact-F-tree-flavour-operator-handoff-v1' and
               source.get('status') == 'SourceFactoredOrderedFlavourMatricesAtF_NotCanonicalWilsonTensor',
               f'{model}: unexpected old ordinary source')
    ensure(ordered.get('scalar_basis_explicit') is True and
           ordered.get('matchete_cg_contractions_retained') is True and
           ordered.get('external_three_family_flavour_order_not_symmetrized') is True,
           f'{model}: tensor metadata incomplete')


def source_alignment(operator: dict, term: dict) -> list[int]:
    axes = term.get('flavour_matrix_axes_external_indices')
    ensure(isinstance(axes, list) and len(axes) == 2 and axes[0] != axes[1],
           'Invalid source Yukawa flavour axes')
    ensure(operator.get('flavour_matrix_axes_external_indices') == axes,
           'Tensor flavour axes differ from Matchete source')
    factors = term.get('factors')
    ensure(isinstance(factors, list) and len(factors) == 2 and
           [f.get('external_index') for f in factors] == axes,
           'Source Yukawa factor indices differ from axes')
    template = term.get('operator_template_unit_flavor_inputform')
    ensure(isinstance(template, str) and 'GammaCC' in template, 'No native operator template')
    lepton_legs = [(m.group('su2'), m.group('flav')) for m in LEPTON.finditer(template)]
    ensure(len(lepton_legs) == 2 and len(set(x[0] for x in lepton_legs)) == 2,
           'Wrong or repeated lepton SU2 index')
    flavours = [f for _, f in lepton_legs]
    ensure(set(flavours) == set(axes), 'Lepton field and Yukawa flavour labels mismatch')
    return [axes.index(flavours[0]), axes.index(flavours[1])]


def operator_rows(operator: dict, term: dict, scalar_dimension: int) -> tuple[list[dict], list[int]]:
    perm = source_alignment(operator, term)
    raw = defaultdict(complex)
    comps = operator['components']
    ensure(operator['component_count'] == len(comps) and len(comps) > 0,
           'Incomplete projected real scalar components')
    for comp in comps:
        i, j = comp['lepton_su2']
        a, b = comp['scalar_real']
        ensure(4 < a <= scalar_dimension and 4 < b <= scalar_dimension,
               'F boundary contains Higgs-only or out-of-range real scalar pair')
        matrix = complex_matrix(comp['ordered_flavour_coefficient_gev_inv'])
        for p in range(3):
            for q in range(3):
                z = matrix[p][q]
                if z:
                    flavours = (p + 1, q + 1)
                    raw[(i, flavours[perm[0]], j, flavours[perm[1]], a, b)] += z
    return rows_from_aligned(raw), perm


def compile_model(model: str, root: Path, ordered_root: Path,
                  poly_root: Path, audit_root: Path) -> tuple[dict, dict]:
    directory = tag(model)
    ofile = ordered_root / directory / 'ordered_tensor.json'
    pfile = poly_root / directory / 'polynomial_bridge.json'
    afile = audit_root / directory / 'normalization_audit.json'
    ordered, polynomial, audit = read(ofile), read(pfile), read(afile)
    sfile = source_file(model, root, ordered_root, ordered)
    source = read(sfile)
    validate_provenance(model, root, ordered_root, ordered, source)
    ensure(polynomial.get('format') == POLY_FORMAT and
           polynomial.get('status') == 'PolynomialRepresentativeVerified_NotSpinorNormalisation',
           f'{model}: missing verified polynomial boundary')
    ensure(audit.get('format') == AUDIT_FORMAT and audit.get('status') == AUDIT_OK,
           f'{model}: missing corrected normalization audit')
    for name, obj in [('polynomial', polynomial), ('audit', audit)]:
        ensure(obj.get('model') == model and float(obj['mu_F_gev']) == float(ordered['mu_F_gev']),
               f'{model}: {name} model/threshold mismatch')
    ensure(polynomial.get('ordered_boundary_sha256') == digest(ofile) and
           polynomial.get('operator_source_sha256') == digest(sfile),
           f'{model}: polynomial provenance stale')
    ensure(audit.get('ordered_tensor_sha256') == digest(ofile) and
           audit.get('polynomial_bridge_sha256') == digest(pfile) and
           audit.get('matchete_operator_source_sha256') == digest(sfile),
           f'{model}: normalization audit provenance stale')
    ch = audit['native_matchete_chirality']
    ensure(ch['native_PL_matched_term_count'] == ordered['tree_operator_count_PL'] and
           ch['native_PR_conjugate_term_count'] == MIXED.get(model, 0) and
           ch['native_chirality_structure_verified'] is True,
           f'{model}: native B/C projector audit failed')
    ensure(polynomial['overall_matchete_spinor_normalization_certified'] is False and
           audit['normalization_findings']['matchete_to_master_spinor_normalization_certified'] is False,
           f'{model}: spinor status is inconsistent')
    terms = checked_terms(model, source, ordered)
    dim = polynomial['real_scalar_dimension']
    all_coordinates: dict[tuple[int, ...], complex] = defaultdict(complex)
    source_rows = []
    for op in ordered['operators']:
        term = terms[op['source_index']]
        rows, perm = operator_rows(op, term, dim)
        # A symmetry-allowed source may vanish at this benchmark. Retain it.
        for key, value in as_map(rows).items():
            all_coordinates[key] += value
        source_rows.append({
            'matchete_source_index': op['source_index'],
            'native_chirality': 'PL',
            'zero_at_this_benchmark': not bool(rows),
            'scalar_and_su2_contractions': 'native Matchete CG already evaluated in ordered component input',
            'cg_registry_used': op.get('cg_registry_used', op.get('cg_names')),
            'flavour_permutation_lepton_order': perm,
            'inputform_operator': term['source_term_inputform'],
            'independent_symmetric_coordinates': rows,
        })
    reconstructed = [dict(coordinate=list(k), coefficient_real=float(v.real),
                          coefficient_imag=float(v.imag), orbit_size=len(orbit(k)))
                     for k, v in sorted(all_coordinates.items()) if v != 0j]
    compare_rows(reconstructed, polynomial['independent_symmetric_coordinates'], model)
    # Use the fully checked original bridge as the numerical summed authority.
    # Separate per-source sums can differ at ~machine epsilon after cancellation.
    merged = polynomial['independent_symmetric_coordinates']
    ensure(len(merged) == polynomial['independent_symmetric_coordinate_count'],
           f'{model}: coordinate count discrepancy')
    payload = {
        'format': FORMAT, 'model': model,
        'status': 'NumericalLLSSBySourceVerified_SpinorNormalizationOpen',
        'mu_F_gev': float(ordered['mu_F_gev']), 'units': 'GeV^-1',
        'field_content': {'heavy_fermions_integrated_out': 'all at mu_F',
                          'heavy_scalar_multiplets': 'both retained in EFT1'},
        'index_order': ['lepton_su2_i', 'lepton_flavour_p', 'lepton_su2_j',
                        'lepton_flavour_q', 'real_scalar_a', 'real_scalar_b'],
        'index_base': 1,
        'real_scalar_dimension': dim,
        'scalar_field_offsets_one_based': polynomial['scalar_field_offsets_one_based'],
        'scalar_field_complex_dimensions': polynomial['scalar_field_complex_dimensions'],
        'symmetric_tensor_orbit_rule': 'All distinct i,p<->j,q and a<->b orbit members equal representative. No extra orbit factor.',
        'formal_polynomial': 'sum C[i,p,j,q,a,b] psi[i,p] psi[j,q] phi[a] phi[b]',
        'spinor_normalization': {
            'certified': False, 'gamma_cc_to_master_weyl_factor': None,
            'gamma_cc_to_master_weyl_sign': None,
            'do_not_use_as_physical_master_wilson_boundary': True,
        },
        'completeness': {
            'F_tree_LLSS_included': True, 'distinct_PL_sources_retained': True,
            'PR_Hermitian_conjugate_not_double_counted': True,
            'one_loop_F_hard_matching_included': False,
            'complete_physical_leading_one_loop_neutrino_mass_certified': False,
            'EFT1_running_performed': False, 'mu_S_scalar_matching_performed': False,
        },
        'sources': source_rows,
        'summed_independent_symmetric_coordinates': merged,
        'checksums_sha256': {'ordered': digest(ofile), 'polynomial': digest(pfile),
                            'normalization_audit': digest(afile), 'native_source': digest(sfile)},
    }
    row = {'model': model, 'status': payload['status'],
           'mu_F_gev': payload['mu_F_gev'],
           'PL_sources': len(source_rows), 'PR_conjugates_excluded': MIXED.get(model, 0),
           'summed_nonzero_coordinates': len(merged),
           'source_resolved_nonzero_coordinates': sum(len(x['independent_symmetric_coordinates'])
                                                       for x in source_rows)}
    return payload, row


def save_json_new(path: Path, obj: dict) -> None:
    ensure(not path.exists(), f'Refusing to overwrite existing artifact: {path}')
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as handle:
        json.dump(obj, handle, indent=2, allow_nan=False)
        handle.write('\n')


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo-root', type=Path,
                    default=Path(__file__).resolve().parents[3])
    ap.add_argument('--check-only', action='store_true',
                    help='Validate all 16 local input sets; write nothing')
    ap.add_argument('--release-id', help='New directory name under lo_fermion_boundary/releases/')
    args = ap.parse_args(argv)
    root = args.repo_root.resolve()
    base = root / 'output/general_observable'
    ordered = base / 'ordered_wilson_tensor_v4'
    polynomial = base / 'general_wilson_polynomial_batch_v2'
    audit = base / 'general_wilson_normalization_audit_v1'
    try:
        compiled = []
        for model in MODELS:
            result, row = compile_model(model, root, ordered, polynomial, audit)
            compiled.append((result, row))
            print(f'{model}: verified {row["PL_sources"]} PL sources, '
                  f'{row["summed_nonzero_coordinates"]} coordinates')
        snapshot = hashlib.sha256(''.join(
            item[0]['checksums_sha256']['ordered'] + item[0]['checksums_sha256']['polynomial']
            + item[0]['checksums_sha256']['normalization_audit'] for item in compiled
        ).encode('ascii')).hexdigest()
        summary = {'format': FORMAT, 'status': '16of16NumericalLLSSTreeVerified_PhysicalNormalizationOpen',
                   'n_models_verified': len(compiled), 'models': [x[1] for x in compiled],
                   'source_snapshot_sha256': snapshot,
                   'absolute_spinor_normalization_certified': False,
                   'complete_leading_one_loop_mass_boundary_certified': False,
                   'physical_eft1_workers_must_reject_this_as_certified_input': True}
        if args.check_only:
            print('CHECK ONLY: 16/16, no files written; absolute spinor factor NOT certified')
            return 0
        base_out = base / 'lo_fermion_boundary' / 'releases'
        identifier = args.release_id or 'source_' + snapshot[:16]
        ensure(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,60}', identifier) is not None,
               'Invalid release-id')
        destination = base_out / identifier
        ensure(not destination.exists(), f'Existing release untouched: {destination}')
        destination.mkdir(parents=True, exist_ok=False)
        for obj, row in compiled:
            relative = tag(obj['model']) + '/boundary.json'
            save_json_new(destination / relative, obj)
            row['path'] = relative
        save_json_new(destination / 'summary.json', summary)
        print(f'CREATED: {destination}')
        print('NOTE: numerical source-resolved LLSS only; GammaCC/Weyl normalization remains open.')
        return 0
    except (OSError, KeyError, TypeError, ValueError, IndexError) as exc:
        print(f'BLOCKED (no physical certification): {type(exc).__name__}: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
