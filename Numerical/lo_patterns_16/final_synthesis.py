"""Read-only, fail-closed consolidation of the *saved* Chat A Batches 1-8.

No loop re-evaluation, no edits to previous outputs, no EFT and no assertion
of physical Matchete/RGBeta-to-Weyl certification. Only standard library.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


BASE = Path('output/general_observable/lo_pattern_study16')
DEST = BASE / 'analysis_final_synthesis_01'
AUDIT = Path('output/general_observable/lo_direct_one_loop/source_audits/charge_resolved_loop_weights.json')
REPORT = Path('Numerical/lo_patterns_16/FINAL_TECHNICAL_SYNTHESIS.md')
SOURCES = {
    'batch01': BASE / 'provisional_01/all16_summary.json',
    'batch02': BASE / 'analysis_batch02_01/batch02_analysis.json',
    'batch03': BASE / 'analysis_batch03_01/batch03_analytic_explanation.json',
    'batch04': BASE / 'analysis_batch04_01/batch04_flavour_analysis.json',
    'batch05': BASE / 'analysis_batch05_01/batch05_synthesis.json',
    'batch06': BASE / 'analysis_batch06_01/batch06_factorization.json',
    'batch07': BASE / 'analysis_batch07_01/batch07_theorem.json',
    'batch08': BASE / 'analysis_batch08_01/batch08_convergence.json',
}
FAMILIES = {
    'A': (-4, -2, 0),
    'B': (-3, -1, 1),
    'C': (-3, -1, 1),
    'D': (-2, 0, 2),
    'E': (-4, -2, 0, 2),
}
MODELS = [f'{fam}_{alpha:+d}' for fam, vals in FAMILIES.items() for alpha in vals]
REFS = ('A_-2', 'C_-1', 'C_+1', 'E_-2')
SCENARIOS = ('baseline', 'scalar_common_shift', 'heavy_generation_tilt',
             'large_scalar_mixing', 'charge_split_diagnostic')
CG_TOTAL = {
    'A': '-2*sqrt(3)/3', 'B': '-2*sqrt(3)/3',
    'C': '2/3', 'D': '-2*sqrt(3)/3', 'E': '2*sqrt(6)/3',
}
CG_WEIGHT = {
    'A': '-2*sqrt(3)/3', 'B': '-2*sqrt(3)/3',
    'C': '2/3', 'D': '-2*sqrt(3)/3', 'E': 'sqrt(6)/3',
}
CG_RATIO = {'A': 1., 'B': 1., 'C': 1. / math.sqrt(3),
            'D': 1., 'E': math.sqrt(2)}
FORMATS = {
    'batch01': 't3-conditional-pattern-study-all16-v1',
    'batch02': 't3-chatA-readonly-batch02-numerical-patterns-v1',
    'batch03': 't3-chatA-batch03-analytic-scalar-kernel-v1',
    'batch04': 't3-chatA-batch04-flavour-geometry-v1',
    'batch05': 't3-chatA-batch05-analytical-yukawa-texture-v1',
    'batch06': 't3-chatA-batch06-heavy-kernel-factorization-v1',
    'batch07': 't3-chatA-batch07-factorization-theorem-v1',
    'batch08': 't3-chatA-batch08-splitting-convergence-v1',
}


class SynthesisError(ValueError):
    """Saved-source inconsistency; do not create a report."""


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise SynthesisError(reason)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), f'Expected JSON object: {path}')
    return obj


def rows(path: Path) -> list[dict[str, str]]:
    with path.open('r', encoding='utf-8-sig', newline='') as fp:
        return list(csv.DictReader(fp))


def unique_model_rows(data: list[dict], origin: str) -> dict[str, dict]:
    require(len(data) == 16 and all(isinstance(r, dict) and 'model' in r for r in data),
            f'{origin}: missing model rows')
    found = [r['model'] for r in data]
    require(len(set(found)) == 16 and set(found) == set(MODELS),
            f'{origin}: duplicate/missing/mislabelled models')
    return {r['model']: r for r in data}


def finite(value: Any, origin: str, *, nonnegative: bool = False) -> float:
    try:
        x = float(value)
    except (TypeError, ValueError, OverflowError) as e:
        raise SynthesisError(f'{origin}: not a finite number: {value!r}') from e
    require(math.isfinite(x) and (not nonnegative or x >= 0), f'{origin}: invalid number')
    return x


def close(x: float, y: float, *, rtol: float = 1e-8, atol: float = 0) -> bool:
    return math.isclose(x, y, rel_tol=rtol, abs_tol=atol)


def verify_inputs(root: Path) -> dict[str, Any]:
    """Validate saved result coverage, counterfactual controls and hash lineage.

    This verifies integrity of *reported* conditional results, not UV Feynman rules.
    """
    root = root.resolve()
    data = {k: read(root / p) for k, p in SOURCES.items()}
    study, b2, b3, b4, b5, b6, b7, b8 = (data[k] for k in SOURCES)
    for name, obj in data.items():
        require(obj.get('format') == FORMATS[name], f'{name}: schema drift')
        require(obj.get('physical_normalisation_certified') is False,
                f'{name}: false/unresolved physical certification status')
        if 'safe_for_physical_source_certification' in obj:
            require(obj['safe_for_physical_source_certification'] is False,
                    f'{name}: false physical-certification promotion')
    require(study['counts']['conditional_evaluated'] == 16 and
            study['counts']['blocked'] == 0 and
            study['counts']['physical_source_certified'] == 0,
            'Batch 1 no longer has 16/16 provisional, 0/16 physical')
    require(b2['counts']['models'] == 16 and b2['counts']['source_scan_rows'] == 288,
            'Batch 2 scan coverage changed')
    row1 = unique_model_rows(study['models'], 'Batch 1')
    row4 = unique_model_rows(b4['models'], 'Batch 4')
    for name in ('batch02', 'batch03', 'batch05'):
        unique_model_rows(data[name]['models'], name)
    for name in ('batch06', 'batch07', 'batch08'):
        require(len(data[name]['models']) == 16 and set(data[name]['models']) == set(MODELS),
                f'{name}: incomplete model list')

    hashes = {key: sha(root / p) for key, p in SOURCES.items()}
    hashes['cg_audit'] = sha(root / AUDIT)
    hashes['final_template'] = sha(root / REPORT)
    h6 = sha(root / BASE / 'analysis_batch06_01/SOURCE_MANIFEST.json')
    h7 = sha(root / BASE / 'analysis_batch07_01/SOURCE_MANIFEST.json')
    h5 = sha(root / BASE / 'analysis_batch05_01/SOURCE_MANIFEST.json')
    h8 = sha(root / BASE / 'analysis_batch08_01/SOURCE_MANIFEST.json')
    hashes.update({'batch05_manifest': h5, 'batch06_manifest': h6,
                   'batch07_manifest': h7, 'batch08_manifest': h8})
    require(b2['input_run_sha256'] == hashes['batch01'] and
            b3['source_summary_sha256'] == hashes['batch01'],
            'Batch 2/3 study SHA drift')
    require(b4['source_sha256']['source_run_sha256'] == hashes['batch01'] and
            b4['source_sha256']['source_batch03_sha256'] == hashes['batch03'],
            'Batch 4 upstream SHA drift')
    require(b5['source_manifest_sha256'] ==
            {'study': hashes['batch01'], 'batch03': hashes['batch03'],
             'batch04': hashes['batch04']}, 'Batch 5 upstream SHA drift')
    require(b6['source_sha256'] ==
            {'study': hashes['batch01'], 'batch05': hashes['batch05'],
             'batch05_manifest': h5}, 'Batch 6 upstream SHA drift')
    require(b7['source_sha256']['batch06_factorization'] == hashes['batch06'] and
            b7['source_sha256']['batch06_manifest'] == h6 and
            b7['source_sha256']['cg_weight_audit'] == hashes['cg_audit'] and
            b7['source_sha256']['study'] == hashes['batch01'],
            'Batch 7 upstream SHA drift')
    require(b8['source_sha256'] == {
                'batch06': hashes['batch06'], 'batch06_manifest': h6,
                'batch07': hashes['batch07'], 'batch07_manifest': h7,
                'study': hashes['batch01']}, 'Batch 8 upstream SHA drift')
    for n in ('batch06', 'batch07', 'batch08'):
        h = read(root / BASE / ('analysis_' + n + '_01') / 'HANDOFF.json')
        require(h.get('physical_normalisation_certified') is False,
                n + ': handoff physical certification drift')
        require(h.get('source_manifest_sha256', h.get('manifest_sha256')) ==
                hashes[n + '_manifest'], n + ': handoff manifest SHA drift')

    audit = read(root / AUDIT)
    require(audit.get('count') == 16 and
            audit.get('status') == 'charge_resolved_group_coefficients',
            'Source CG audit missing or unexpected')
    audit_rows: dict[str, dict] = {}
    for entry in audit['models']:
        model = f"{entry['class']}_{entry['model']['alpha']:+d}"
        require(model not in audit_rows and model in MODELS,
                'Source CG audit duplicate/unexpected model')
        audit_rows[model] = entry
    require(set(audit_rows) == set(MODELS), 'Incomplete source CG audit')
    audited = []
    for model in MODELS:
        family = model[0]
        entry = audit_rows[model]
        active = [ch for ch in entry['channels'] if ch['nonzero']]
        require(len(active) == (2 if family == 'E' else 1),
                f'{model}: source CG active charge count changed')
        require(entry['total_group_coefficient'] == CG_TOTAL[family] and
                all(ch['complex_weight'] == CG_WEIGHT[family] for ch in active),
                f'{model}: source CG convention or group coefficient drift')
        raw = row1[model]
        require(raw.get('status') == 'conditional_uv_evaluated' and
                raw.get('physical_normalisation_certified') is False and
                row4[model].get('physical_normalisation_certified') is False,
                f'{model}: conditional model changed')
        masses = [finite(raw[f'm{i}_ev'], f'{model}.m{i}', nonnegative=True)
                  for i in (1, 2, 3)]
        require(0 <= masses[0] <= masses[1] <= masses[2],
                f'{model}: Takagi masses not sorted')
        require(all(close(masses[i], finite(row4[model][f'm{i+1}_ev'], model), rtol=3e-7)
                    for i in range(3)), f'{model}: Batch 4 mass mismatch')
        audit_sha = read(root / BASE / 'provisional_01' / 'models' / (model + '.json'))
        require(audit_sha['source']['audit_hashes']['charge_resolved_loop_weights.json'] == hashes['cg_audit'],
                f'{model}: source CG audit SHA changed')
        audit_prob = len(active)
        require(int(raw['candidate_topology_charge_channel_count']) == audit_prob,
                f'{model}: charge count differs from audit')
        audited.append({'model': model,
                        'charge_channels': ', '.join(x['charge'] for x in active),
                        'source_CG_total': CG_TOTAL[family],
                        'source_CG_magnitude_ratio': CG_RATIO[family],
                        'm1_ev': masses[0], 'm2_ev': masses[1], 'm3_ev': masses[2],
                        'm1_over_m3': masses[0]/masses[2] if masses[2] > 0 else None,
                        'fermion': row4[model]['fermion_type'],
                        'light_mode_A': finite(row4[model]['lightest_mode']['geometric_projection_fraction'], model),
                        'light_mode_C': finite(row4[model]['lightest_mode']['coherence_fraction'], model)
                            if row4[model]['lightest_mode']['coherence_fraction'] is not None else None})

    # Cross-check the actual saved 16x5 heavy-coefficient scenario results.
    test_rows = rows(root / BASE / 'analysis_batch06_01/batch06_coefficient_tests.csv')
    require(len(test_rows) == 80 and
            {(r['scenario'], r['model']) for r in test_rows} ==
            {(scenario, model) for scenario in SCENARIOS for model in MODELS},
            'Batch 6 coefficient scenarios incomplete or duplicated')
    by_case = {(r['scenario'], r['model']): r for r in test_rows}
    scenarios = []
    for scenario in SCENARIOS:
        worst = max((r for r in test_rows if r['scenario'] == scenario),
                    key=lambda r: finite(r['coefficient_residual'], 'scenario', nonnegative=True))
        scenarios.append({'scenario': scenario, 'worst_model': worst['model'],
                          'max_residual': finite(worst['coefficient_residual'], 'worst', nonnegative=True)})
    require(scenarios[0]['max_residual'] < 1e-8, 'Lost baseline factorization')
    require(scenarios[-1]['max_residual'] > scenarios[0]['max_residual']*1e3,
            'Charge-splitting signature no longer present')
    for r in audited:
        k = finite(by_case[('baseline', r['model'])]['absolute_scale'], 'kappa', nonnegative=True)
        r['measured_kappa_magnitude'] = k
        r['baseline_coefficient_residual'] = finite(by_case[('baseline', r['model'])]['coefficient_residual'], 'baseline', nonnegative=True)
        require(abs(k - r['source_CG_magnitude_ratio']) < 2e-5,
                r['model'] + ': source CG/loop coefficient magnitude mismatch')

    # Full-matrix comparisons remove the 'same singular ratios' ambiguity.
    pairs = rows(root / BASE / 'analysis_batch05_01/batch05_matched_pairs.csv')
    require(len(pairs) == 7 and len({(r['model_a'], r['model_b']) for r in pairs}) == 7,
            'Batch 5 matched-pair coverage changed')
    pairwise = []
    for r in pairs:
        for ref in REFS:
            key = f'common_{ref}_matrix_residual'
            value = finite(r[key], key, nonnegative=True)
            require(value < 1e-8, 'Common texture matrices no longer proportional: ' + key)
        pairwise.append({'model_a': r['model_a'], 'model_b': r['model_b'],
                         'frozen_residual': finite(r['frozen_matrix_proportional_residual'], 'frozen', nonnegative=True),
                         'worst_matched_texture_residual': max(float(r[f'common_{ref}_matrix_residual']) for ref in REFS)})

    controls = rows(root / BASE / 'analysis_batch05_01/batch05_common_texture_controls.csv')
    require(len(controls) == 16*len(REFS) and
            {(r['model'], r['texture_reference']) for r in controls} ==
            {(m, ref) for m in MODELS for ref in REFS},
            'Batch 5 texture controls incomplete or duplicated')
    for ref in ('C_-1', 'C_+1'):
        values = [finite(r['m1_over_m3'], 'matched texture ratio', nonnegative=True)
                  for r in controls if r['texture_reference'] == ref]
        require(max(values) - min(values) < 1e-5 * max(values),
                'Common-texture spectrum universality no longer holds: ' + ref)

    b7_rows = rows(root / BASE / 'analysis_batch07_01/batch07_first_corrections.csv')
    require(len(b7_rows) == 16, 'Batch 7 correction coverage incomplete')
    for r in b7_rows:
        require(finite(r['large_mix_cubic_relative_error'], 'cubic', nonnegative=True) <
                finite(r['large_mix_linear_relative_error'], 'linear', nonnegative=True),
                'Cubic correction no longer improves: ' + r['model'])

    b8_rows = rows(root / BASE / 'analysis_batch08_01/batch08_model_summary.csv')
    slope = unique_model_rows(b8_rows, 'Batch 8')
    for model in MODELS:
        s = slope[model]
        q = finite(s['exact_remainder_order'], model + '.remainder_order')
        require(1.5 < q < 2.5, model + ': lost quadratic Taylor remainder')
        t = s['transverse_signal_order']
        if t not in ('', None, 'None', 'null'):
            t = finite(t, model + '.transverse_order')
            require(.6 < t < 1.4, model + ': lost linear transverse onset')
        else:
            t = None
        s['exact_remainder_order'] = q
        s['transverse_signal_order'] = t
    deriv = rows(root / BASE / 'analysis_batch08_01/batch08_derivative_validation.csv')
    require(len(deriv) == 16, 'Batch 8 derivative coverage incomplete')
    max_deriv = max(finite(r['central_difference_slope_relative_error'], 'central difference', nonnegative=True)
                    for r in deriv)
    require(max_deriv < 5e-4, 'Exact derivative verification failed')

    return {'format': 't3-chatA-final-synthesis-evidence-v1',
            'physical_normalisation_certified': False,
            'safe_for_physical_source_certification': False,
            'model_count': 16, 'source_sha256': hashes, 'model_table': audited,
            'scenario_maxima': scenarios, 'matched_pairs': pairwise,
            'batch07_first_corrections': b7_rows,
            'batch08_convergence': [slope[m] for m in MODELS],
            'max_exact_derivative_finite_difference_error': max_deriv,
            'conditional_only': True}


def render_appendix(d: dict[str, Any]) -> str:
    out = [
        '# Local verified evidence appendix — Chat A', '',
        '**The following tables are verified directly against this checkout.** '
        'They remain source-conditional; `physical_normalisation_certified=false`.', '',
        '## All sixteen frozen spectra', '',
        '| Model | Active charges | CG magnitude ratio | m1 [eV] | m2 [eV] | m3 [eV] | m1/m3 | A1 | C1 |',
        '|---|---|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for r in d['model_table']:
        fmt = lambda k: f'{r[k]:.9g}' if r[k] is not None else 'undefined'
        out.append(f'| {r["model"]} | {r["charge_channels"]} | '
                   f'{fmt("source_CG_magnitude_ratio")} | {fmt("m1_ev")} | '
                   f'{fmt("m2_ev")} | {fmt("m3_ev")} | {fmt("m1_over_m3")} | '
                   f'{fmt("light_mode_A")} | {fmt("light_mode_C")} |')
    out += ['', '## Factorization stress-test maxima', '',
            '| Scenario | Worst model | Maximum coefficient residual |',
            '|---|---|---:|']
    for x in d['scenario_maxima']:
        out.append(f'| {x["scenario"]} | {x["worst_model"]} | {x["max_residual"]:.9g} |')
    out += ['', '## Matched texture matrix tests', '',
            '| Model A | Model B | Frozen matrix residual | Worst controlled matrix residual |',
            '|---|---|---:|---:|']
    for x in d['matched_pairs']:
        out.append(f'| {x["model_a"]} | {x["model_b"]} | '
                   f'{x["frozen_residual"]:.9g} | {x["worst_matched_texture_residual"]:.9g} |')
    out += ['', '## Convergence by model', '',
            '| Model | Exact Taylor remainder order | Transverse onset order | Endpoint nonproportionality |',
            '|---|---:|---:|---:|']
    for x in d['batch08_convergence']:
        order = x['transverse_signal_order']
        q = 'noise floor' if order is None else f'{order:.8g}'
        out.append(f'| {x["model"]} | {x["exact_remainder_order"]:.8g} | {q} | '
                   f'{float(x["endpoint_nonproportionality"]):.9g} |')
    out += ['', f'Largest analytic-vs-central-difference derivative mismatch: '
            f'`{d["max_exact_derivative_finite_difference_error"]:.9g}`.', '',
            '## Status', '',
            'This verifies saved numbers and hash provenance, **not** the pending '
            'four-spinor/CG physical normalization. No new amplitudes, fits, UV loop '
            'computations, Matchete runs or EFT results were produced.', '']
    return '\n'.join(out)


def _write_csv(path: Path, records: list[dict], fields: list[str]) -> None:
    with path.open('w', encoding='utf-8', newline='') as fp:
        w = csv.DictWriter(fp, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(records)


def generate(root: Path, *, output_dir: Path | None = None) -> dict[str, Any]:
    root = root.resolve()
    owner = (root / BASE).resolve()
    target = (output_dir or root / DEST).resolve()
    require(target.is_relative_to(owner) and target != owner,
            'Output must be inside the new Chat A analysis tree')
    require(not target.exists() or (target.is_dir() and not any(target.iterdir())),
            'Output directory already populated: ' + str(target))
    # Read, validate and render BEFORE making any new directory.
    d = verify_inputs(root)
    text = (root / REPORT).read_text(encoding='utf-8')
    require(text.startswith('# Chat A — All-16') and '0/16' in text,
            'Technical report template changed')
    source_paths = {k: str(p) for k, p in SOURCES.items()}
    manifest = {'format': 't3-chatA-final-synthesis-manifest-v1',
                'source_sha256': d['source_sha256'],
                'source_paths': source_paths,
                'physically_certified_count': 0,
                'physical_normalisation_certified': False,
                'report': 'FINAL_VERIFIED_TECHNICAL_REPORT.md'}
    handoff = {'format': 't3-chatA-final-synthesis-handoff-v1',
               'models': 16, 'conditional_model_count': 16,
               'physically_certified_count': 0, 'physical_normalisation_certified': False,
               'safe_for_conditional_interpretation': True,
               'safe_for_physical_source_certification': False,
               'theorem': 'If B_R,g = kappa_R B_0,g for all heavy groups, then Mnu_R = kappa_R Mnu_0 at matched Yukawas.',
               'unresolved': ['B1 Matchete-to-Weyl physical spinor/normalization map',
                              'B2 cross-package SU2/coupling phase normalization'],
               'source_digest_manifest': 'SOURCE_MANIFEST.json'}
    target.mkdir(parents=True, exist_ok=True)
    (target / 'FINAL_VERIFIED_TECHNICAL_REPORT.md').write_text(
        text + '\n\n---\n\n' + render_appendix(d), encoding='utf-8')
    (target / 'consolidated_evidence.json').write_text(
        json.dumps(d, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    (target / 'SOURCE_MANIFEST.json').write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    handoff['source_manifest_sha256'] = sha(target / 'SOURCE_MANIFEST.json')
    (target / 'HANDOFF.json').write_text(
        json.dumps(handoff, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    (target / 'HANDOFF.md').write_text(
        '# Chat A final analytical synthesis handoff\n\n'
        '**16/16 conditional; 0/16 physically source certified.**\n\n'
        'The complete results and theorem are in `FINAL_VERIFIED_TECHNICAL_REPORT.md`. '
        'No existing file was modified. Checks include source hash lineage from all '
        'eight batches, all-16 CG weights, 80 heavy-coefficient comparisons, '
        '28 matched-texture matrix residuals and 16 splitting convergence studies.\n\n'
        'The source-derived Weyl normalization, relative phases and cross-package '
        'SU2 dictionary are NOT certified. The deliberate charge-splitting path '
        'is not a demonstrated gauge-invariant UV interaction. No EFT work.\n',
        encoding='utf-8')
    _write_csv(target / 'all16_verified.csv', d['model_table'],
               ['model', 'charge_channels', 'source_CG_total', 'source_CG_magnitude_ratio',
                'measured_kappa_magnitude', 'baseline_coefficient_residual',
                'fermion', 'm1_ev', 'm2_ev', 'm3_ev', 'm1_over_m3',
                'light_mode_A', 'light_mode_C'])
    _write_csv(target / 'scenario_maxima.csv', d['scenario_maxima'],
               ['scenario', 'worst_model', 'max_residual'])
    _write_csv(target / 'matched_pairs.csv', d['matched_pairs'],
               ['model_a', 'model_b', 'frozen_residual', 'worst_matched_texture_residual'])
    return d
