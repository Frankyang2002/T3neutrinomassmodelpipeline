"""Source-resolved first threshold LLSS validation and fail-closed handoff.

Inputs are the existing immutable ``lo_fermion_boundary/releases/source_*/``
files (not reconstructed from any historical or wrong-chirality data).
No physical phase or additional hard-one-loop matching term is invented.
"""
from __future__ import annotations
from collections import defaultdict
import codecs
from hashlib import sha256
from pathlib import Path
import json
import math
import re
from typing import Any

MODELS = tuple(f'{f}_{a:+d}' for f, aa in (
    ('A', (-4, -2, 0)), ('B', (-3, -1, 1)),
    ('C', (-3, -1, 1)), ('D', (-2, 0, 2)),
    ('E', (-4, -2, 0, 2)),
) for a in aa)
# Verified source snapshot inventory from the immutable source_021dcc... summary.
# Fields: PL source terms, excluded PR HC, independent summed rows,
# total source-resolved rows. These are not calculated from assumptions.
EXPECTED = {
    'A_-4': (1, 0, 60, 60), 'A_-2': (4, 0, 40, 120),
    'A_+0': (4, 0, 40, 120), 'B_-3': (1, 0, 48, 48),
    'B_-1': (10, 10, 108, 372), 'B_+1': (1, 0, 48, 48),
    'C_-3': (2, 0, 48, 48), 'C_-1': (12, 12, 108, 372),
    'C_+1': (2, 0, 48, 48), 'D_-2': (4, 0, 30, 120),
    'D_+0': (4, 0, 30, 120), 'D_+2': (1, 0, 60, 60),
    'E_-4': (2, 0, 120, 120), 'E_-2': (8, 0, 60, 240),
    'E_+0': (8, 0, 60, 240), 'E_+2': (2, 0, 120, 120),
}
SOURCE_FORMAT = 't3-lo-fermion-boundary-source-resolved-provisional-v1'
OUTPUT_FORMAT = 't3-lo-spinor-eft1-boundary-provisional-v1'
SUMMARY_FORMAT = 't3-lo-spinor-all16-handoff-v1'
INDEX_ORDER = ['lepton_su2_i', 'lepton_flavour_p', 'lepton_su2_j',
               'lepton_flavour_q', 'real_scalar_a', 'real_scalar_b']
RELEASE_REL = 'output/general_observable/lo_fermion_boundary/releases'
PROVENANCE_REL = {
    'ordered': 'output/general_observable/ordered_wilson_tensor_v4/{tag}/ordered_tensor.json',
    'polynomial': 'output/general_observable/general_wilson_polynomial_batch_v2/{tag}/polynomial_bridge.json',
    'normalization_audit': 'output/general_observable/general_wilson_normalization_audit_v1/{tag}/normalization_audit.json',
}
_GAMMACC = re.compile(r'DiracProduct\[GammaCC,\s*Proj\[-1\]\]')
_PL = re.compile(r'Proj\[\s*-1\s*\]')
_PR = re.compile(r'Proj\[\s*1\s*\]')
_YUK = re.compile(r'(?:Bar\[)?Coupling\[(y[12]),')


class BoundaryError(ValueError):
    pass


def require(ok: bool, why: str) -> None:
    if not ok:
        raise BoundaryError(why)


def tag(model: str) -> str:
    return model.replace('-', 'm').replace('+', 'p')


def file_sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    require(path.is_file(), f'Missing source: {path}')
    d = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(d, dict), f'Expected object: {path}')
    return d


def orbit(key: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    i, p, j, q, a, b = key
    return tuple(sorted({(i, p, j, q, a, b), (j, q, i, p, a, b),
                         (i, p, j, q, b, a), (j, q, i, p, b, a)}))


def canonical(k: tuple[int, ...]) -> tuple[int, ...]:
    i, p, j, q, a, b = k
    if (j, q) < (i, p):
        i, p, j, q = j, q, i, p
    return (i, p, j, q, min(a, b), max(a, b))


def extract_rows(rows: list[dict], scalar_dim: int) -> dict[tuple[int, ...], complex]:
    require(isinstance(rows, list), 'Missing canonical coefficient rows')
    dest = {}
    for row in rows:
        k = tuple(row['coordinate'])
        require(len(k) == 6 and all(isinstance(x, int) and not isinstance(x, bool) for x in k),
                'Nonintegral coordinate')
        i, p, j, q, a, b = k
        require(i in (1, 2) and j in (1, 2) and p in (1, 2, 3) and q in (1, 2, 3),
                'Bad SU2 or flavour coordinate')
        require(5 <= a <= scalar_dim and 5 <= b <= scalar_dim,
                'Non-BSM or out-of-range scalar')
        require(k == canonical(k) and k not in dest, 'Duplicate/noncanonical index')
        require(row.get('orbit_size') == len(orbit(k)), 'Wrong 1/2/4 symmetry orbit')
        r, im = row['coefficient_real'], row['coefficient_imag']
        require(isinstance(r, (int, float)) and isinstance(im, (int, float))
                and math.isfinite(r) and math.isfinite(im), 'Bad complex coefficient')
        require(r != 0 or im != 0, 'Explicit zero row should be omitted')
        dest[k] = complex(r, im)
    return dest


def compare_maps(a: dict, b: dict, context: str, tol: float = 5e-12) -> float:
    keys = set(a) | set(b)
    scale = max((abs(x) for x in b.values()), default=0.)
    require(scale > 0, f'{context}: zero total source tensor')
    err = max((abs(a.get(k, 0j) - b.get(k, 0j)) for k in keys), default=0.)
    require(err <= scale*tol, f'{context}: mismatch {err/scale:.3g} relative')
    return err / scale


def validate_model(d: dict, model: str, strict_counts: bool = True) -> dict:
    """Numerically reconstruct every canonical Wilson entry from *all* PL sources.

    No additional orbit factors or Hermitian-conjugate PR terms are summed.
    """
    require(model in MODELS, f'Unknown T3 model {model}')
    require(d.get('model') == model and d.get('format') == SOURCE_FORMAT,
            f'{model}: unexpected historical boundary format')
    require(d.get('status') == 'NumericalLLSSBySourceVerified_SpinorNormalizationOpen',
            f'{model}: source is not the expected independently verified provisional release')
    require(d.get('index_order') == INDEX_ORDER and d.get('index_base') == 1,
            f'{model}: incompatible master coordinate order')
    require(d.get('units') == 'GeV^-1' and d.get('mu_F_gev', 0) > 0,
            f'{model}: missing common fermion scale or units')
    require(d.get('field_content', {}).get('heavy_fermions_integrated_out') == 'all at mu_F'
            and d.get('field_content', {}).get('heavy_scalar_multiplets') == 'both retained in EFT1',
            f'{model}: forbidden partial threshold ordering')
    complete = d.get('completeness', {})
    require(complete.get('F_tree_LLSS_included') is True and
            complete.get('distinct_PL_sources_retained') is True and
            complete.get('PR_Hermitian_conjugate_not_double_counted') is True,
            f'{model}: incomplete first-threshold tree source manifest')
    require(complete.get('one_loop_F_hard_matching_included') is False and
            complete.get('complete_physical_leading_one_loop_neutrino_mass_certified') is False,
            f'{model}: unsupported assertion of hard matching/completeness')
    norm = d.get('spinor_normalization', {})
    require(norm.get('certified') is False and norm.get('gamma_cc_to_master_weyl_factor') is None
            and norm.get('gamma_cc_to_master_weyl_sign') is None,
            f'{model}: external source silently fixed undefined Weyl phase')
    dims = d.get('scalar_field_complex_dimensions')
    offsets = d.get('scalar_field_offsets_one_based')
    require(isinstance(dims, dict) and isinstance(offsets, dict) and
            2 <= len(dims) <= 3 and set(dims) == set(offsets),
            f'{model}: invalid scalar basis block inventory')
    # For charged non-self-conjugate multiplets the original complex field can be
    # represented as separate T3Scalar{1,2}Plus/Minus real blocks. Those named
    # blocks each carry *d* real coordinates, not 2*d. They are NOT an
    # additional independent heavy scalar multiplet or an extra threshold.
    spans = {name: (dims[name] if ('Plus' in name or 'Minus' in name)
                    else 2*dims[name]) for name in dims}
    scalar_dim = d.get('real_scalar_dimension')
    require(scalar_dim == 4 + sum(spans.values()), f'{model}: real scalar basis dimension')
    eligible = set()
    for name, start in offsets.items():
        require(isinstance(start, int) and dims[name] > 0, 'Invalid scalar field block')
        eligible.update(range(start, start+spans[name]))
    require(eligible == set(range(5, scalar_dim+1)), 'Overlapping/missing scalar components')
    sources = d.get('sources')
    require(isinstance(sources, list) and sources, f'{model}: no native sources')
    total_rows = 0
    seen = set()
    sums: dict[tuple[int, ...], list[complex]] = defaultdict(list)
    grouped_sources = []
    for s in sources:
        index = s.get('matchete_source_index')
        require(isinstance(index, int) and index > 0 and index not in seen,
                f'{model}: duplicate native source index')
        seen.add(index)
        native = s.get('inputform_operator', '')
        require(isinstance(native, str) and _GAMMACC.search(native) is not None
                and _PL.search(native) is not None and _PR.search(native) is None
                and 'NCM[Transp[Field[l, Fermion' in native,
                f'{model} source {index}: not a holomorphic Matchete GammaCC PL bilinear')
        require(s.get('native_chirality') == 'PL', 'PR conjugate included as extra C5')
        require(s.get('flavour_permutation_lepton_order') in ([0, 1], [1, 0]),
                f'{model}: unknown native flavour permutation')
        cg = s.get('cg_registry_used')
        require(isinstance(cg, list) and bool(cg) and all(isinstance(x, str) and x for x in cg),
                f'{model}: missing native CG registry')
        source_map = extract_rows(s['independent_symmetric_coordinates'], scalar_dim)
        require((len(source_map) == 0) == bool(s['zero_at_this_benchmark']),
                f'{model}: inconsistent nonzero/zero source label {index}')
        for key, z in source_map.items():
            require(key[4] in eligible and key[5] in eligible, f'{model}: illegal scalar component')
            sums[key].append(z)
        total_rows += len(source_map)
        grouped_sources.append({
            'source_id': f'{tag(model)}:MatchetePL:{index}',
            'matchete_source_index': index, 'chirality': 'PL',
            'native_operator_inputform': native,
            'native_CG_registries': cg,
            'yukawa_channels_in_native_operator': sorted(set(_YUK.findall(native))),
            'lepton_order_flavour_axis_permutation': s['flavour_permutation_lepton_order'],
            'zero_at_benchmark': bool(s['zero_at_this_benchmark']),
            'independent_coordinates': s['independent_symmetric_coordinates'],
        })
    summed = {k: complex(math.fsum(x.real for x in vals), math.fsum(x.imag for x in vals))
              for k, vals in sums.items()}
    summed = {k: v for k, v in summed.items() if v != 0}
    target = extract_rows(d['summed_independent_symmetric_coordinates'], scalar_dim)
    residual = compare_maps(summed, target, model)
    if strict_counts:
        exp_pl, exp_pr, exp_sum, exp_src = EXPECTED[model]
        require((len(sources), len(target), total_rows) == (exp_pl, exp_sum, exp_src),
                f'{model}: native source/coefficient count does not match pinned all16 inventory')
    else:
        exp_pr = EXPECTED[model][1]
    return {
        'model': model,
        'mu_F_gev': d['mu_F_gev'],
        'verified_source_count_PL': len(sources),
        'native_PR_conjugates_excluded_by_upstream': exp_pr,
        'source_coordinate_rows': total_rows,
        'summed_coordinate_rows': len(target),
        'max_relative_source_sum_residual': residual,
        'scalar_multiplets_active_in_EFT1': list(dims),
        'source_terms': grouped_sources,
        'summed_independent_coordinates': d['summed_independent_symmetric_coordinates'],
    }


def verify_reference_hashes(root: Path, model: str, d: dict) -> None:
    checksums = d.get('checksums_sha256')
    require(isinstance(checksums, dict), f'{model}: missing reference SHA256 provenance')
    for key, template in PROVENANCE_REL.items():
        path = root / template.format(tag=tag(model))
        require(path.is_file(), f'{model}: missing checked reference {path}')
        require(file_sha(path) == checksums.get(key), f'{model}: changed source reference {key}')
    # The original source JSON is usually under general_wilson_boundary_v1,
    # while B/C mixed sources are stored in ordered_wilson_tensor_v4.
    if model in ('B_-1', 'C_-1'):
        native_path = root / f'output/general_observable/ordered_wilson_tensor_v4/{tag(model)}/source_flavour.json'
    else:
        native_path = root / f'output/general_observable/general_wilson_boundary_v1/{tag(model)}/wilson_boundary.json'
    require(native_path.is_file(), f'{model}: missing original matched source {native_path}')
    require(file_sha(native_path) == checksums.get('native_source'),
            f'{model}: native source was replaced since original release')


def weave_source_handoff(d: dict, result: dict, source_sha: str) -> dict:
    """Stable, complete numerical three-generation *provisional* interface."""
    return {
        'format': OUTPUT_FORMAT,
        'model': result['model'],
        'status': 'CONDITIONAL_TREE_LLSS_BLOCKED_PHYSICAL_NORMALIZATION',
        'source_boundary_sha256': source_sha,
        'scale_gev': result['mu_F_gev'],
        'matching_history': [
            {'transition': 'UV->EFT1', 'threshold': 'mu_F',
             'heavy_BSM_fermions': 'ALL simultaneously',
             'heavy_BSM_scalars': 'BOTH remain active', 'SM_fermions': 'active'},
            {'transition': 'EFT1->SMEFT', 'threshold': 'mu_S',
             'heavy_BSM_scalars': 'BOTH simultaneously', 'status': 'NOT_PERFORMED'},
        ],
        'operator': 'LLSS',
        'index_order': INDEX_ORDER, 'index_base': 1,
        'flavour_generations': 3,
        'scalar_real_dimension': d['real_scalar_dimension'],
        'scalar_field_offsets_one_based': d['scalar_field_offsets_one_based'],
        'scalar_field_complex_dimensions': d['scalar_field_complex_dimensions'],
        'units': 'GeV^-1',
        'relative_formal_tensor_convention': d['formal_polynomial'],
        'permutation_orbit': 'swap BOTH lepton SU2+flavour legs and/or scalar real indices; no additional orbit factor',
        'source_resolved_unscaled': result['source_terms'],
        'summed_unscaled': result['summed_independent_coordinates'],
        'verified_relative_numeric_checks': {
            'PL_source_count': result['verified_source_count_PL'],
            'PR_excluded': result['native_PR_conjugates_excluded_by_upstream'],
            'source_coordinate_rows': result['source_coordinate_rows'],
            'summed_coordinate_rows': result['summed_coordinate_rows'],
            'max_relative_source_sum_residual': result['max_relative_source_sum_residual'],
        },
        'physical_master_weyl_boundary': None,
        'normalization': {
            'matchete_GammaCC_to_master_Weyl_complex_factor': None,
            'certified': False, 'may_multiply_unscaled_tensor': False,
            'candidate_C_PL_over_epsilon': '+1',
            'candidate_Cinv_PL_over_epsilon': '-1',
            'Matchete_GammaCC_is_C_or_Cinv': None,
        },
        'inclusions': {'fermion_tree_LLSS': True, 'fermion_hard_one_loop': False,
                       'EFT1_running': False, 'scalar_matching': False,
                       'complete_leading_one_loop_mass': False},
        'downstream_usage': {
            'schema': OUTPUT_FORMAT,
            'safe_for_nonphysical_EFT1_operator_structural_tests': True,
            'safe_for_physical_EFT1_or_scalar_matching': False,
            'accept_if_physical_normalization_required': False,
            'requires_original_source_terms_retained': True,
        },
        'blocking_evidence': [
            'No independent installed-Matchete GammaCC to master-Weyl signed bilinear identification',
            'No native all-model connected four-point Dirac+Majorana+CG amplitude certification',
            'Raw hard-one-loop fermion matching not included in source release',
        ],
    }


def checked_write(path: Path, data: dict) -> None:
    """Immutable write, refuse to modify another worker's existing output."""
    rendered = (json.dumps(data, indent=2, sort_keys=True, ensure_ascii=False,
                           allow_nan=False) + '\n').encode('utf-8')
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        require(path.read_bytes() == rendered, f'Refusing to overwrite existing {path}')
        return
    with path.open('xb') as file:
        file.write(rendered)


def certify_source_release(root: Path, release: Path, output: Path | None,
                           strict_references: bool = True,
                           strict_counts: bool = True) -> dict:
    """Check all sixteen models as a single atomic *validation* batch.

    All files are validated before the first output file is written. Unlike
    a physical release, status remains non-certified on every path.
    """
    summary_path = release / 'summary.json'
    upstream = load(summary_path)
    require(upstream.get('format') == SOURCE_FORMAT and
            upstream.get('n_models_verified') == 16 and
            upstream.get('absolute_spinor_normalization_certified') is False and
            upstream.get('complete_leading_one_loop_mass_boundary_certified') is False,
            'Input release is not the sixteen-model provisional source snapshot')
    listed = upstream.get('models')
    require(isinstance(listed, list) and len(listed) == 16 and
            {m.get('model') for m in listed} == set(MODELS), 'Incorrect model inventory')
    mu = set()
    built = {}
    checks = []
    for item in listed:
        model = item['model']
        path = release / item['path']
        require(path.resolve().is_relative_to(release.resolve()), 'Unsafe input release path')
        d = load(path)
        if strict_references:
            verify_reference_hashes(root, model, d)
        result = validate_model(d, model, strict_counts=strict_counts)
        require((result['verified_source_count_PL'], result['native_PR_conjugates_excluded_by_upstream'],
                 result['summed_coordinate_rows'], result['source_coordinate_rows']) ==
                (item['PL_sources'], item['PR_conjugates_excluded'],
                 item['summed_nonzero_coordinates'], item['source_resolved_nonzero_coordinates']),
                f'{model}: upstream summary does not match source')
        mu.add(result['mu_F_gev'])
        checks.append({'model': model, 'PL': result['verified_source_count_PL'],
                       'excluded_PR': result['native_PR_conjugates_excluded_by_upstream'],
                       'summed_rows': result['summed_coordinate_rows'],
                       'source_rows': result['source_coordinate_rows'],
                       'max_relative_residual': result['max_relative_source_sum_residual']})
        built[tag(model)] = weave_source_handoff(d, result, file_sha(path))
    require(len(mu) == 1, 'Fermion thresholds are not simultaneous')
    if strict_counts:
        require(sum(c['PL'] for c in checks) == sum(v[0] for v in EXPECTED.values()),
                'Missing PL terms')
    manifest = {
        'format': SUMMARY_FORMAT,
        'status': '16of16_PROVISIONAL_SOURCE_NUMERICS_CHECKED_PHYSICAL_BLOCKED',
        'source_release_summary_sha256': file_sha(summary_path),
        'source_snapshot_sha256': upstream.get('source_snapshot_sha256'),
        'mu_F_gev': next(iter(mu)),
        'n_models_structurally_verified': len(built),
        'n_models_physically_normalized': 0,
        'n_models_with_complete_fermion_hard_one_loop_matching': 0,
        'n_models_ready_for_physical_EFT1_running': 0,
        'models': checks,
        'model_files': {m: f'{tag(m)}/boundary.json' for m in MODELS},
        'preserves_source_terms': True,
        'gamma_cc_to_master_weyl_factor': None,
        'all_heavy_fermions_one_mu_F': True,
        'both_heavy_scalars_one_later_mu_S': True,
        'SM_fermions_retained': True,
        'no_files_outside_exclusive_worker_paths_changed': True,
    }
    if output is not None:
        # Computed source-addressed path avoids overwriting a parallel worker.
        folder = output / ('source_' + file_sha(summary_path)[:16])
        for model in MODELS:
            checked_write(folder / tag(model) / 'boundary.json', built[tag(model)])
        checked_write(folder / 'summary.json', manifest)
    return manifest


def _decode_native_evidence(raw: bytes) -> tuple[str, str]:
    """Read Mathematica/PowerShell text without silently accepting broken UTF-8.

    Windows PowerShell 5.1 redirection can produce UTF-16LE even when the
    Mathematica text is ASCII. In UTF-8 with errors='replace', every second
    byte then contains NUL, causing all source-convention checks to be false.
    """
    if raw.startswith(codecs.BOM_UTF8):
        return raw.decode('utf-8-sig'), 'utf-8-sig'
    if raw.startswith((codecs.BOM_UTF16_LE, codecs.BOM_UTF16_BE)):
        return raw.decode('utf-16'), 'utf-16'
    head = raw[:1024]
    if len(head) >= 8 and len(raw) % 2 == 0:
        even = head[0::2]
        odd = head[1::2]
        odd_nuls = odd.count(0) / max(len(odd), 1)
        even_nuls = even.count(0) / max(len(even), 1)
        if odd_nuls > 0.4 and even_nuls < 0.1:
            return raw.decode('utf-16-le'), 'utf-16-le'
        if even_nuls > 0.4 and odd_nuls < 0.1:
            return raw.decode('utf-16-be'), 'utf-16-be'
    return raw.decode('utf-8'), 'utf-8'


def native_evidence_probe(path: Path | None) -> dict:
    """Report textual package evidence, NEVER a physical C/GammaCC phase."""
    if path is None or not path.is_file():
        return {'status': 'NOT_AVAILABLE', 'GammaCC_Phase_certified': False,
                'GammaCC_to_master_weyl_factor': None}
    raw_bytes = path.read_bytes()
    provenance = {'source_sha256': sha256(raw_bytes).hexdigest()}
    try:
        raw, encoding = _decode_native_evidence(raw_bytes)
    except UnicodeError:
        return {'status': 'UNDECODABLE_PACKAGE_SOURCE', **provenance,
                'native_items': {}, 'Matchete_pl_is_left_projector': False,
                'GammaCC_Phase_certified': False,
                'GammaCC_to_master_weyl_factor': None}
    found = ['GammaCC::usage =', 'PackageExport["GammaCC"]',
             'CConj[', 'Proj[-1]', 'DiracProduct[GammaCC]']
    presence = {x: x in raw for x in found}
    return {'status': 'PACKAGE_SOURCE_LOCATED_ONLY_NOT_NUMERICAL_MATRIX_PHASE',
            **provenance, 'detected_text_encoding': encoding,
            'native_items': presence,
            # Matchete source usage explicitly specifies the +/- -> R/L mapping.
            'Matchete_pl_is_left_projector':
                'Proj[+/-1]' in raw and 'right/left' in raw,
            'GammaCC_Phase_certified': False,
            'GammaCC_to_master_weyl_factor': None}
