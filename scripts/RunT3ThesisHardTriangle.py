#!/usr/bin/env python3
"""Recalculate thesis 32 fixed-Yukawa transfers with the computed F-hard triangle.

This wraps the EXISTING thesis table generator and the EXISTING source-resolved
hard-triangle injector. No existing source is overwritten, no archived output
is reused as a new result, and no `complete one-loop matching` claim is made.

Required user-supplied arguments are the same as BuildT3GeneralisedThesisTables.py:
  --reference-results <existing verified fit directory> [additional directories]
  --output output/parallel_t3/<NEW_DIRECTORY>

The injector uses the freshly generated per-trial B background and A/C scalar
input, not the archived B_all16 background of the exploratory triangle batch.
Only a full successful 32/32 run is certified as containing the triangle in
*every* thesis table row. Stage-frozen figures are NOT regenerated here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import traceback
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.partial')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n', encoding='utf-8')
    tmp.replace(path)


def main(argv: list[str] | None = None) -> int:
    from Numerical.parallel_t3.rough_inputs import builder as builder_module
    from Numerical.parallel_t3.wilson import fermion_hard_triangle_batch1 as hard_module
    from scripts import BuildT3GeneralisedThesisTables as thesis

    args = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--output', required=True)
    parser.add_argument('--reference-results', nargs='+', required=True)
    parsed, _unknown = parser.parse_known_args(args)
    output = (ROOT / parsed.output).resolve()
    allowable = (ROOT / 'output' / 'parallel_t3').resolve()
    if allowable not in output.parents:
        parser.error('Output must be a NEW directory underneath output/parallel_t3')
    if output.exists():
        parser.error('Refusing to overwrite an existing output directory: ' + str(output))

    # A separate, immutable record of exactly which trial was injected.
    # No generic hard-triangle run's Yukawas are reused or transplanted.
    entries: list[dict] = []
    count_by_model: dict[str, int] = {}
    references = ('D_-2', 'A_-4')
    original = builder_module.build_model
    expected = 32

    def build_with_hard_triangle(model, opts):
        built = original(model, opts)
        if built.get('status') != 'written_provisional_A_C_inputs':
            raise RuntimeError('Thesis boundary build failed before hard matching: ' + str(built))
        ordinal = count_by_model.get(model, 0)
        if ordinal >= 2:
            raise RuntimeError('Unexpected repeated thesis model: ' + model)
        count_by_model[model] = ordinal + 1
        ref = references[ordinal]
        trial_id = str(opts.out_run_id)
        if trial_id != str(opts.B_run_id) or not trial_id.startswith('fit16_'):
            raise RuntimeError('Hard matching requires same-run B and rough inputs: ' + trial_id)
        if model not in thesis.all_models():
            raise RuntimeError('Unknown model in hard-matching thesis batch: ' + model)
        destination = output / 'hard_fermion' / ref / model
        # Original Batch1 hard module is hard-wired to a saved B_RUN. The
        # temporary override below redirects it to THIS trial's freshly
        # computed B(MF) source, masses and Yukawas. Runs are sequential.
        with patch.object(hard_module, 'B_RUN', trial_id):
            evidence = hard_module.install_hard_boundary(
                ROOT, model,
                input_run_id=trial_id,
                output=destination,
                gamma_factor=float(opts.gamma_factor),
            )
        magnitude = float(evidence['MF_hard_C5_frobenius'])
        if not np.isfinite(magnitude) or magnitude <= 0:
            raise RuntimeError('No nonzero finite hard C5 was computed for ' + model)
        if not (destination / 'fermion_hard_triangle.npz').is_file():
            raise RuntimeError('Hard C5 tensor evidence not saved for ' + model)
        entries.append({
            'reference': ref, 'model': model,
            'trial_id': trial_id,
            'hard_C5_frobenius_gev_inverse': magnitude,
            'hard_manifest': evidence['fermion_hard_manifest'],
            'hard_npz': evidence['fermion_hard_npz'],
            'status': 'HARD_TRIANGLE_INJECTED_BEFORE_EFT1_EVOLUTION',
            'complete_one_loop_matching': False,
        })
        _write_json(output / 'hard_triangle_injection_ledger.json', {
            'status': 'IN_PROGRESS',
            'runs_injected': len(entries),
            'required_runs': expected,
            'entries': entries,
            'hard_scope': 'scalar-quartic fermion triangle only',
        })
        return built

    # The normal generator regenerates all 32 low-scale C5's, neutrino
    # observables and LaTeX tables; only its boundary-builder is wrapped.
    # The existing figures are not modified by that generator.
    runner_rc = 1
    try:
        with patch.object(builder_module, 'build_model', build_with_hard_triangle):
            runner_rc = thesis.main(args)
    except Exception:
        _write_json(output / 'hard_triangle_injection_ledger.json', {
            'status': 'FAILED', 'runs_injected': len(entries),
            'required_runs': expected, 'entries': entries,
            'traceback': traceback.format_exc(limit=20),
        })
        raise

    results_path = output / 'results.json'
    results = json.loads(results_path.read_text(encoding='utf-8')) if results_path.is_file() else {}
    completed = all(
        results.get(ref, {}).get(model, {}).get('status') == 'VALID_EXPLORATORY'
        for ref in references for model in thesis.all_models()
    )
    exact = (
        runner_rc == 0 and len(entries) == expected and completed
        and all(count_by_model.get(model) == 2 for model in thesis.all_models())
        and len({(e['reference'], e['model']) for e in entries}) == expected
    )
    # Explicitly fail closed: never describe a partial run as having hard
    # matching in the full thesis benchmark.
    _write_json(output / 'hard_triangle_injection_ledger.json', {
        'status': 'ALL_32_RECOMPUTED_WITH_TRIANGLE' if exact else 'INCOMPLETE_NOT_FOR_THESIS',
        'runs_injected': len(entries), 'required_runs': expected,
        'entries': entries, 'complete_one_loop_matching': False,
        'all_32_low_scale_results_recomputed': bool(exact),
        'thesis_figures_recomputed': False,
        'hard_scope': 'scalar-quartic fermion triangle only; conditional normalization',
    })
    if exact:
        # The unmodified evaluator reports the OLD hard-zero assumption in its
        # Trial.meta. Correct every result and per-case JSON only after a
        # full successful rerun and matching evidence for each pair.
        by_key = {(e['reference'], e['model']): e for e in entries}
        for ref in references:
            for model in thesis.all_models():
                entry = by_key[(ref, model)]
                record = results[ref][model]
                meta = record.get('meta')
                if not isinstance(meta, dict) or meta.get('fermion_hard_finite_C5_assumed_zero') is not True:
                    raise RuntimeError('Unexpected original trial metadata: ' + ref + '/' + model)
                if meta.get('scratch_token') != entry['trial_id']:
                    raise RuntimeError('Result/injection trial mismatch: ' + ref + '/' + model)
                meta['fermion_hard_finite_C5_assumed_zero'] = False
                meta['fermion_hard_triangle_included'] = True
                meta['fermion_hard_triangle_scope'] = 'scalar-quartic triangle only'
                meta['fermion_hard_triangle_manifest'] = entry['hard_manifest']
                meta['complete_one_loop_fermion_matching'] = False
                _write_json(output / ref / (model + '.json'), record)
        _write_json(results_path, results)
        provenance_path = output / 'provenance.json'
        provenance = json.loads(provenance_path.read_text(encoding='utf-8'))
        provenance['finite_F_hard_C5_zero_benchmark'] = False
        provenance['F_hard_triangle_included_all_32'] = True
        provenance['F_hard_triangle_scope'] = 'one-loop heavy-fermion HHSS triangle only'
        provenance['complete_one_loop_fermion_matching'] = False
        provenance['F_hard_triangle_evidence'] = 'hard_triangle_injection_ledger.json'
        provenance['thesis_figures_recomputed'] = False
        provenance['conditional_not_physical'] = True
        _write_json(provenance_path, provenance)
        summary_path = output / 'summary.json'
        summary = json.loads(summary_path.read_text(encoding='utf-8'))
        summary['F_hard_triangle_included_all_32'] = True
        summary['complete_one_loop_fermion_matching'] = False
        summary['figures_updated'] = False
        _write_json(summary_path, summary)
        print('PASS: all 32 thesis transfers recalculated with the finite fermion-hard triangle.')
        print('New tables: ' + str(output))
        print('NOTE: all published figures remain OLD until separately regenerated.')
        return 0
    print('INCOMPLETE: some cases failed. DO NOT claim hard matching in the full thesis benchmark.')
    print('Inspect ' + str(output / 'hard_triangle_injection_ledger.json'))
    return 1


if __name__ == '__main__':
    raise SystemExit(main())
