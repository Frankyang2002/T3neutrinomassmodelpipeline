"""Read-only 3-study/16-model scalar-boundary audit; no fitting or RG integration.

    python -m scripts.AuditT3SavedStudyScalarBoundaries

The original study trees are read-only inputs. Outputs are confined to a new
output/studies/scalar_potential_boundary_audit/ directory. A distinct
output/studies/general_scalar_extended/ directory is reserved for later runs.
"""
from __future__ import annotations

import csv
from hashlib import sha256
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = Path('RGE/general/generated/general_rgbeta_quartic_registry.json')
OUTPUT = Path('output/studies/scalar_potential_boundary_audit')
# Case labels are deliberately explicit: UV-on and fermion-start are different
# numerical experiments, even when their fitted lambdaT3 values agree.
STUDIES = {
    'small_uv_running_5e-9': {
        'root': Path('output/studies/lambda_5e-9_refit'),
        'case': None, 'lambda_ad': 5e-9, 'boundary_scale_gev': 1e7,
        'uv_running': True,
    },
    'small_uv_frozen_5e-9': {
        'root': Path('output/studies/fermion_start_lambda_comparison/5e-9'),
        'case': '5e-9', 'lambda_ad': 5e-9, 'boundary_scale_gev': 1e5,
        'uv_running': False,
    },
    'large_uv_frozen_0p5': {
        'root': Path('output/studies/fermion_start_lambda_comparison/0p5'),
        'case': '0p5', 'lambda_ad': 0.5, 'boundary_scale_gev': 1e5,
        'uv_running': False,
    },
}
MODEL_NAMES = (
    'A_+0', 'A_-2', 'A_-4', 'B_+1', 'B_-1', 'B_-3',
    'C_+1', 'C_-1', 'C_-3', 'D_+0', 'D_+2', 'D_-2',
    'E_+0', 'E_+2', 'E_-2', 'E_-4',
)
REQUIRED_TOTALS = {
    'matchete_extra_quartics': 50,
    'rgbeta_preexisting_structural': 12,
    'rgbeta_new_vertices': 38,
    'quadratic_couplings_blocked': 8,
}
FIT_SCREEN = 5.86


def fail(message: str):
    raise ValueError('T3 saved-study scalar-boundary audit: ' + message)


def source(root: Path, rel: Path, hashes: dict) -> dict:
    file = root / rel
    if not file.is_file():
        fail('Missing saved source: ' + rel.as_posix())
    raw = file.read_bytes()
    hashes[rel.as_posix()] = sha256(raw).hexdigest()
    try:
        result = json.loads(raw.decode('utf-8-sig'))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        fail('Invalid JSON: ' + rel.as_posix() + ': ' + str(exc))
    if not isinstance(result, dict):
        fail('Expected JSON object: ' + rel.as_posix())
    return result


def finite_number(v, label: str) -> float:
    if type(v) not in (float, int) or not math.isfinite(float(v)):
        fail('Expected finite numerical ' + label)
    return float(v)


def complex_scalar(value, name: str) -> complex:
    if not isinstance(value, dict) or set(value) != {'real', 'imag'}:
        fail(name + ': expected real/imag scalar')
    return complex(finite_number(value['real'], name + '.real'),
                   finite_number(value['imag'], name + '.imag'))


def close_to(x, y, label: str, *, rel=1e-10, abs_tol=1e-14):
    if not math.isclose(finite_number(x, label), finite_number(y, label),
                        rel_tol=rel, abs_tol=abs_tol):
        fail(label + ': mismatch: ' + repr(x) + ' != ' + repr(y))


def model_file(entry: dict) -> str:
    dims = entry.get('dimensions')
    a = entry.get('alpha')
    if (not isinstance(dims, list) or len(dims) != 3 or
            any(type(x) is not int or x not in (1, 2, 3) for x in dims) or
            type(a) is not int):
        fail('Invalid registry representation')
    d1, d2, df = dims
    token = ('m' if a < 0 else 'p') + str(abs(a))
    return f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{token}.json'


def registry_inventory(registry: dict):
    if (registry.get('status') != 'GeneratedCandidateBasisNotWolframVerified'
            or any(registry.get('totals', {}).get(k) != v
                   for k, v in REQUIRED_TOTALS.items())
            or set(registry.get('models', {})) != set(MODEL_NAMES)):
        fail('Unexpected general scalar registry/version/coverage')
    inventory = {}
    names_seen = set()
    total_reused = total_new = total_quadratics = 0
    for name in MODEL_NAMES:
        entry = registry['models'][name]
        if not isinstance(entry, dict) or entry.get('model') != name:
            fail(name + ': registry model mismatch')
        new = tuple(x['name'] for x in entry['new_rgbeta_symbols'])
        reused = tuple(x['existing_name'] for x in entry['preexisting_rgbeta_symbols'])
        nquad = entry.get('quadratics_not_implemented')
        if (len(set(new + reused)) != len(new + reused) or
                len(new + reused) != entry.get('total_additional_quartics') or
                any(not x.startswith('lambdaT3GeneralQ') for x in new) or
                type(nquad) is not int or nquad not in (0, 1)):
            fail(name + ': malformed general scalar inventory')
        stem = model_file(entry)
        if stem in names_seen:
            fail('Duplicate registered fitted file stem: ' + stem)
        names_seen.add(stem)
        inventory[name] = (entry, stem, new, reused, nquad)
        total_new += len(new)
        total_reused += len(reused)
        total_quadratics += nquad
    if (total_new, total_reused, total_quadratics) != (38, 12, 8):
        fail('Incorrect global operator inventory')
    return inventory


def check_case(root: Path, study_name: str, spec: dict, inventory: dict, hashes: dict):
    directory = spec['root']
    comparison = None
    if spec['case'] is not None:
        comparison = source(root, directory / 'comparison/optimised/summary.json', hashes)
        if (comparison.get('mode') != 'optimised' or
                comparison.get('case') != spec['case'] or
                comparison.get('uv_running') is not False or
                comparison.get('start_scale') != 'mu_F' or
                not isinstance(comparison.get('models'), list) or
                len(comparison['models']) != 16):
            fail(study_name + ': invalid frozen comparison summary')
        by_name = {row.get('model'): row for row in comparison['models']}
        if len(by_name) != 16:
            fail(study_name + ': duplicate comparison model entries')
    else:
        by_name = {}
    records = []
    for short_name in MODEL_NAMES:
        entry, stem, new, reused, nquad = inventory[short_name]
        model = stem.removesuffix('.json')
        config_rel = directory / 'fitted' / stem
        fit_rel = directory / 'fits' / stem
        config = source(root, config_rel, hashes)
        fit = source(root, fit_rel, hashes)
        rep = config.get('representation')
        d1, d2, df = entry['dimensions']
        if rep != {'d_s1': d1, 'd_s2': d2, 'd_f': df,
                   'alpha': entry['alpha'], 'shared_scalar': False}:
            fail(study_name + '/' + model + ': fitted representation differs from registry')
        original = config.get('base_state')
        if not isinstance(original, dict) or not isinstance(original.get('ordinary'), dict):
            fail(model + ': missing ordinary fitted state')
        if finite_number(original.get('mu_gev'), model + '.mu_gev') != 1e7:
            fail(model + ': unexpected stored UV reference scale')
        scales = config.get('scales', {})
        for field, scale in [('mu_fermion_threshold_gev', 1e5),
                             ('mu_scalar_threshold_gev', 1e3), ('mu_low_gev', 1e2)]:
            close_to(scales.get(field), scale, model + '.' + field)
        values = original['ordinary']
        lam = complex_scalar(values.get('lambdaT3'), model + '.lambdaT3')
        expected = spec['lambda_ad'] * (0.5 if short_name[0] == 'E' else 1.0)
        close_to(lam.real, expected, model + '.lambdaT3', rel=1e-10, abs_tol=0)
        if lam.imag != 0.0:
            fail(model + ': unexpected complex phase of topology quartic')
        if (fit.get('model') != model or
                fit.get('status') not in {'Success', 'PoorFit'}):
            fail(study_name + '/' + model + ': fit record mismatch/unknown status')
        close_to(fit.get('lambdaT3'), expected, model + '.fit.lambdaT3', abs_tol=0)
        chi2 = finite_number(fit.get('chi2_fit'), model + '.chi2_fit')
        if chi2 < 0:
            fail(model + ': negative fit chi2')
        passed = fit['status'] == 'Success' and chi2 <= FIT_SCREEN
        if spec['case'] is None:
            if config.get('fermion_start_study') is not None or not isinstance(config.get('lambda_5e9_study'), dict):
                fail(model + ': expected UV-on fitting provenance')
        else:
            meta = config.get('fermion_start_study')
            if (not isinstance(meta, dict) or meta.get('case') != spec['case'] or
                    not isinstance(meta.get('boundary_definition'), str) or
                    'UV RGEs disabled' not in meta['boundary_definition']):
                fail(model + ': expected frozen fermion-start provenance')
            close_to(meta.get('lambdaT3'), expected, model + '.metadata.lambdaT3', abs_tol=0)
            row = by_name.get(model)
            if not isinstance(row, dict):
                fail(model + ': missing frozen optimised comparison row')
            close_to(row.get('lambdaT3'), expected, model + '.comparison.lambdaT3', abs_tol=0)
            if row.get('chi2') is not None:
                close_to(row['chi2'], chi2, model + '.comparison.chi2', rel=1e-7, abs_tol=1e-6)
        # Original EFT state contains some extra quartics (the reused 12 model
        # slots). The 38 new quartic slots are not represented in these old
        # input files: absence is NOT an explicitly stored zero boundary.
        explicit_zeros = []
        nonzero = []
        for field in reused:
            if field not in values:
                fail(study_name + '/' + model + ': missing represented extra quartic ' + field)
            v = complex_scalar(values[field], model + '.' + field)
            (explicit_zeros if v == 0 else nonzero).append(field)
        unexpected_new = sorted(set(new) & set(values))
        if unexpected_new:
            fail(model + ': new general quartic appears in legacy fitted state: ' + ','.join(unexpected_new))
        # Do not silently certify zero for an unrepresented bilinear.
        unexpected_mass = sorted(k for k in values if k.startswith('muT3General'))
        if unexpected_mass:
            fail(model + ': candidate general quadratic appears in legacy state: ' + ','.join(unexpected_mass))
        records.append({
            'study': study_name, 'model': short_name, 'fitted_model': model,
            'lambdaT3': lam.real, 'boundary_scale_gev': spec['boundary_scale_gev'],
            'uv_running': spec['uv_running'],
            'fit_status': fit['status'], 'chi2_fit': chi2,
            'passes_chi2_screen': passed,
            'extra_quartic_slots': len(new) + len(reused),
            'explicitly_zero_reused_count': len(explicit_zeros),
            'explicitly_zero_reused_names': explicit_zeros,
            'nonzero_reused_names': nonzero,
            'new_quartics_unrepresented_count': len(new),
            'new_quartics_unrepresented_names': list(new),
            'quadratic_slots_unrepresented': nquad,
            'fitted_config_path': config_rel.as_posix(),
            'fit_record_path': fit_rel.as_posix(),
        })
    return records


def build_audit(root: Path) -> dict:
    root = Path(root)
    hashes = {}
    registry = source(root, REGISTRY, hashes)
    inventory = registry_inventory(registry)
    cases = {}
    for name, spec in STUDIES.items():
        rows = check_case(root, name, spec, inventory, hashes)
        cases[name] = {
            'boundary_scale_gev': spec['boundary_scale_gev'],
            'uv_running': spec['uv_running'],
            'lambdaT3_A_to_D': spec['lambda_ad'],
            'lambdaT3_E': spec['lambda_ad'] / 2,
            'models_checked': len(rows),
            'passes_chi2_screen': sum(r['passes_chi2_screen'] for r in rows),
            'explicit_zero_reused_slots': sum(r['explicitly_zero_reused_count'] for r in rows),
            'new_slots_unrepresented': sum(r['new_quartics_unrepresented_count'] for r in rows),
            'quadratic_slots_unrepresented': sum(r['quadratic_slots_unrepresented'] for r in rows),
            'nonzero_reused_slots': sum(len(r['nonzero_reused_names']) for r in rows),
            'models': rows,
        }
    if len(cases) != 3 or any(v['models_checked'] != 16 or
                              v['explicit_zero_reused_slots'] + v['nonzero_reused_slots'] != 12 or
                              v['new_slots_unrepresented'] != 38 or
                              v['quadratic_slots_unrepresented'] != 8
                              for v in cases.values()):
        fail('Incomplete three-study operator inventory')
    return {'format': 't3-restricted-three-study-scalar-boundary-audit-v1',
            'scope': 'Saved restricted fitted inputs; no numerical rerun',
            'status': ('ExplicitReusedQuarticsZeroAcrossAllCases'
                       if all(v['nonzero_reused_slots'] == 0 for v in cases.values())
                       else 'NonzeroReusedQuarticsFound'),
            'historical_source_files_modified': False,
            'new_general_scalar_quartics_explicitly_certified_zero': False,
            'missing_general_scalar_mass_parameters_explicitly_certified_zero': False,
            'physical_matchete_rgbeta_mapping_certified': False,
            'cases': cases, 'source_sha256': hashes}


def summary_md(report: dict) -> str:
    lines = [
        '# Three historical T3 studies — scalar boundary audit', '',
        'Read-only inspection of saved fitted configurations and fit records. '
        'No fits, RGBeta integrations, Matchete matching or historical outputs were rerun.', '',
        '| Study | lambdaT3 A–D (E) | Input boundary | Models | Chi2 <= 5.86 and fit Success | Explicit reused quartics zero | New quartics absent | Quadratics absent |',
        '|---|---:|---|---:|---:|---:|---:|---:|',
    ]
    for name, c in report['cases'].items():
        interpretation = ('mu_UV=1e7 GeV (UV on)' if c['uv_running']
                          else 'mu_F=1e5 GeV (UV frozen)')
        lines.append(f"| `{name}` | {c['lambdaT3_A_to_D']:.8g} ({c['lambdaT3_E']:.8g}) | "
                     f"{interpretation} | {c['models_checked']} | {c['passes_chi2_screen']} | "
                     f"{c['explicit_zero_reused_slots']} / 12 | "
                     f"{c['new_slots_unrepresented']} / 38 | "
                     f"{c['quadratic_slots_unrepresented']} / 8 |")
    lines.extend(['',
        '## Meaning of the checks', '',
        '- **Explicitly zero:** one of the 12 model-slots of additional quartics already represented in the old fitted input state, recorded as exactly zero real and imaginary parts.',
        '- **Absent / unrepresented:** one of the 38 newer quartic slots or eight additional bilinears is not a parameter in the saved fitted input. It is **not** an independently checked zero-valued coefficient in a general-potential model.',
        '- UV-frozen fitted inputs retain `base_state.mu_gev=1e7` metadata, but their numerical values are **imposed at mu_F=1e5**; UV equations were not integrated.',
        '- Failed / poor-fitting cases remain in the output, identified by the fit record and chi2 threshold. Matching/fit success and operator completeness are different tests.',
        '', '## Scope of the scalar-potential conclusion', '',
        'The original restricted scalar potential remains a well-defined fixed-order leading-one-loop benchmark when absent additional terms are **assumed zero by model restriction**. The audit does not reconstruct a fully general scalar potential or certify quantum corrections, nonzero coupling maps or scalar vacuum stability.',
        '', '## Separate future general-scalar numerical studies', '',
        'Store new results only under `output/studies/general_scalar_extended/` with distinct run identifiers, input hashes, RG/matching scope and original-comparison references. **Never write to** `output/studies/lambda_5e-9_refit/`, `output/studies/fermion_start_lambda_comparison/` or the saved Stage 36 / Batch 7–8 output trees. See `docs/T3_GENERAL_SCALAR_STUDY_SEPARATION.md`.',
        '', '## Reproducibility', '',
        f"- Source JSON files fingerprinted: **{len(report['source_sha256'])}** (SHA-256 keys use portable `/` separators).",
        f"- Audit status: `{report['status']}`.",
        '- A separate CSV contains the 48 model-study records, including every poor-fit status and the distinction between explicit zero and absent parameters.',
        ''])
    return '\n'.join(lines)


def write_audit(root: Path, report: dict):
    root = Path(root)
    out = root / OUTPUT
    # The audit may overwrite its OWN regenerated summary, never source cases.
    if not (out.resolve().is_relative_to((root / 'output/studies').resolve()) and
            out.resolve() != (root / 'output/studies').resolve()):
        fail('Audit output is outside dedicated study namespace')
    out.mkdir(parents=True, exist_ok=True)
    (out / 'restricted_input_audit.json').write_text(
        json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    (out / 'restricted_input_audit.md').write_text(summary_md(report), encoding='utf-8')
    with (out / 'restricted_input_audit.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['study', 'model', 'lambdaT3', 'boundary_scale_gev', 'uv_running',
                         'fit_status', 'chi2_fit', 'passes_chi2_screen',
                         'extra_quartic_slots', 'explicitly_zero_reused_count',
                         'new_quartics_unrepresented_count',
                         'quadratic_slots_unrepresented', 'nonzero_reused_names',
                         'fitted_config_path', 'fit_record_path'])
        for case in report['cases'].values():
            for row in case['models']:
                writer.writerow([row[k] if k != 'nonzero_reused_names' else
                                 ';'.join(row[k]) for k in
                                 ('study', 'model', 'lambdaT3', 'boundary_scale_gev',
                                  'uv_running', 'fit_status', 'chi2_fit',
                                  'passes_chi2_screen', 'extra_quartic_slots',
                                  'explicitly_zero_reused_count',
                                  'new_quartics_unrepresented_count',
                                  'quadratic_slots_unrepresented', 'nonzero_reused_names',
                                  'fitted_config_path', 'fit_record_path')])
    return out


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT,
                        help='Local repository root (default: this script parent)')
    ns = parser.parse_args()
    try:
        report = build_audit(ns.root)
        out = write_audit(ns.root, report)
    except (OSError, ValueError) as exc:
        parser.exit(2, str(exc) + '\n')
    print('Three-study scalar audit: 48 model-case configurations; '
          + str(len(report['source_sha256'])) + ' source files SHA-256 fingerprinted')
    for name, case in report['cases'].items():
        print(name + ': ' + str(case['explicit_zero_reused_slots']) + '/12 explicit reused zeros, '
              + str(case['new_slots_unrepresented']) + '/38 new quartics absent, '
              + str(case['quadratic_slots_unrepresented']) + '/8 quadratic slots absent; '
              + str(case['passes_chi2_screen']) + '/16 pass fit screen')
    print('Output:', out / 'restricted_input_audit.md')
    if report['status'] != 'ExplicitReusedQuarticsZeroAcrossAllCases':
        parser.exit(2, 'NONZERO additional represented quartic(s) found: inspect report.\n')


if __name__ == '__main__':
    main()
